#include "algorithm.hpp"
#include "typed_state.hpp"
#include "catalog.hpp"
#include "scheduler.hpp"
#include <chrono>
#include <iostream>
#include <memory>

namespace zipper_verify {
namespace {
constexpr const char *ABI="zipper-independent-verifier-v1";
const char *identity_fields[]={"schema","scope","source_sha256","python_reference_source_sha256",
                              "catalog_sha256","coefficient_ring","native_build"};
bool bool_field(const J &object,const char *name) {
    const J &value=object.at(name);
    if (!value.is_boolean()) throw ReplayError(std::string(name)+" must be boolean");
    return value.get<bool>();
}
std::size_t budget(const J &request,const char *name,I fallback) {
    I value=request.contains(name)?as_i(request.at(name)):fallback;
    if (value<0) throw ReplayError("cache budgets must be nonnegative integers");
    return static_cast<std::size_t>(value);
}
J collected_output(const J &terms) {
    if (!terms.is_array()) throw ReplayError("completed output must be an array");
    std::map<std::string,Laurent> collected;
    for (const J &term:terms) {
        if (!term.is_object() || term.size()!=3 || !term.contains("coefficient") ||
            !term.contains("word") || !term.contains("components") ||
            !term.at("word").is_array() || !term.at("components").is_array())
            throw ReplayError("completed output term has an invalid shape");
        J identity=term;
        identity.erase("coefficient");
        auto &sum=collected[canonical_json(identity)];
        for (const auto &[q,c]:coefficient(term.at("coefficient"))) sum[q]+=c;
    }
    J result=J::object();
    for (const auto &[key,sum]:collected) {
        J scalar=encoded_coefficient(sum);
        if (!scalar.empty()) result[key]=std::move(scalar);
    }
    return result;
}
class Protocol {
    std::unique_ptr<Catalog> catalog_;
    std::unique_ptr<ReplaySession> session_;
    J identity_;
    J reconstruct(const J &request) {
        I n=as_i(request.at("n"));
        const J &word=request.at("input_word");
        Limits limits(request.at("limits"));
        I batch=as_i(request.at("batch_size"));
        bool trace=bool_field(request,"trace_enabled");
        return session_->replay(n,word,limits,batch,trace);
    }
public:
    J dispatch(const J &request) {
        if (!request.is_object() || !request.contains("op") || !request.at("op").is_string())
            throw ReplayError("request requires a string operation");
        std::string op=request.at("op");
        if (op=="initialize") {
            if (session_) throw ReplayError("verifier is already initialized");
            if (request.contains("abi") && request.at("abi")!=ABI) throw ReplayError("verifier ABI differs");
            const J &identity=request.at("native_identity");
            if (!identity.is_object()) throw ReplayError("native identity must be an object");
            for (auto key:identity_fields) if (!identity.contains(key)) throw ReplayError("native identity is missing a required field");
            if (identity.at("schema")!="zipper-native-run-v1" ||
                identity.at("scope")!="literal_word_native_collected_spanning_procedure" ||
                identity.at("coefficient_ring")!="Z[q,q^-1]") throw ReplayError("native certificate contract differs");
            if (!request.at("catalog_sha256").is_string() ||
                !strict_equal(request.at("catalog_sha256"),identity.at("catalog_sha256")) ||
                !strict_equal(request.at("reference_sources"),identity.at("python_reference_source_sha256")))
                throw ReplayError("initialization source or catalog identities differ");
            if (!request.at("catalog_text").is_string()) throw ReplayError("exact catalog text is required");
            std::string text=request.at("catalog_text").get<std::string>();
            if (sha256(text)!=request.at("catalog_sha256").get<std::string>() ||
                !strict_equal(J::parse(text),request.at("catalog"))) throw ReplayError("catalog bytes differ from their bound identity");
            set_reference_sources(request.at("reference_sources"));
            catalog_=std::make_unique<Catalog>(request.at("catalog"),request.at("catalog_sha256").get<std::string>());
            session_=std::make_unique<ReplaySession>(*catalog_,budget(request,"cache_entries",262144),
                                                     budget(request,"cache_bytes",128*1024*1024));
            identity_=identity;
            return {{"ready",true},{"abi",ABI}};
        }
        if (!session_) throw ReplayError("verifier initialization is required");
        if (op=="cache_info") return session_->cache_info();
        if (op=="clear_cache") { session_->clear_cache(); return {{"cleared",true}}; }
        if (op=="replay") return reconstruct(request);
        if (op=="replay_reused")
            return session_->replay_reused(as_i(request.at("n")),request.at("input_word"),Limits(request.at("limits")));
        // The typed diagnostic accepts the owned controller schema. The legacy
        // operation below preserves unrestricted historical diagnostic inputs;
        // replay and replay_reused use the typed controller exclusively.
        if (op=="typed_transition") {
            bool trace=request.contains("trace_enabled")?bool_field(request,"trace_enabled"):
                       request.contains("trace")?bool_field(request,"trace"):false;
            const J normalized=normalize_coefficients(request.at("work"));
            auto work=typed::decode_work(normalized);
            // A typed diagnostic has a closed schema at every owned record.
            // Round-trip equality also checks nested contexts and component
            // records without introducing JSON fields into the owned state.
            if (!strict_equal(normalized,typed::encode_work(work)))
                throw ReplayError("typed controller state differs from its owned schema");
            return typed::encode_transition(typed::algorithm_transition(work,*catalog_,
                                          Limits(request.at("limits")),trace));
        }
        if (op=="transition") {
            bool trace=request.contains("trace_enabled")?bool_field(request,"trace_enabled"):
                       request.contains("trace")?bool_field(request,"trace"):false;
            return algorithm_transition(normalize_coefficients(request.at("work")),*catalog_,
                                        Limits(request.at("limits")),trace);
        }
        if (op=="verify_reused") {
            const auto started=std::chrono::steady_clock::now();
            const J &certificate=request.at("certificate");
            if (!certificate.is_object()) throw ReplayError("reused result must be an object");
            if (certificate.size()!=17 || !certificate.contains("metrics") || !certificate.at("metrics").is_object())
                throw ReplayError("reused result has an invalid shape");
            for (auto key:identity_fields) {
                if (std::string(key)=="schema" || std::string(key)=="scope") continue;
                if (!certificate.contains(key) || !strict_equal(certificate.at(key),identity_.at(key)))
                    throw ReplayError(std::string("reused result identity differs: ")+key);
            }
            const J execution={{"schedule","intrinsic-completion-dag-v1"},
                {"collection","full-output-except-scalar-v1"},{"local_term_floor",4},{"local_support_floor",2}};
            if (certificate.at("schema")!="zipper-native-reuse-v1" ||
                certificate.at("scope")!="literal_word_native_reused_spanning_procedure" ||
                bool_field(certificate,"trace_enabled") || certificate.at("status")!="complete" ||
                !certificate.at("diagnostic").is_null() ||
                !strict_equal(certificate.at("pending"),J::array()) ||
                !strict_equal(certificate.at("execution"),execution))
                throw ReplayError("complete reused-result contract differs");
            // Producer telemetry is deliberately not trusted or replayed.
            J actual=collected_output(certificate.at("output"));
            J expected=session_->replay_reused(as_i(certificate.at("n")),certificate.at("input_word"),Limits(certificate.at("limits")));
            const bool complete=expected.at("status")=="complete";
            J normal=complete?collected_output(expected.at("output")):J::object();
            const bool valid=complete && strict_equal(normal,actual);
            J result={{"ok",valid},{"metrics",expected.at("metrics")},{"cache",session_->cache_info()},
                {"verification_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()},
                {"verification_mode","independent-intrinsic-completion-v1"},
                {"replay_status",expected.at("status")},{"replay_diagnostic",expected.at("diagnostic")}};
            if (valid) result["output_sha256"]=digest(normal);
            else result["error"]=complete?"reused result differs from independent C++ completion":
                "independent C++ completion did not finish: "+expected.at("status").get<std::string>();
            return result;
        }
        if (op=="verify") {
            const auto started=std::chrono::steady_clock::now();
            const J &certificate=request.at("certificate");
            if (!certificate.is_object()) throw ReplayError("native certificate must be an object");
            for (auto key:identity_fields)
                if (!certificate.contains(key) || !strict_equal(certificate.at(key),identity_.at(key)))
                    throw ReplayError(std::string("native certificate identity differs: ")+key);
            J parameters={{"n",certificate.at("n")},{"input_word",certificate.at("input_word")},
                {"limits",certificate.at("limits")},{"trace_enabled",certificate.at("trace_enabled")},
                {"batch_size",certificate.at("execution").at("batch_size")}};
            J expected=reconstruct(parameters);
            J actual=certificate;
            actual["schema"]="zipper-collected-run-v1";
            actual["scope"]="literal_word_collected_spanning_procedure";
            actual["source_sha256"]=actual.at("python_reference_source_sha256");
            actual.erase("python_reference_source_sha256"); actual.erase("native_build");
            bool valid=strict_equal(expected,actual);
            double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
            J result={{"ok",valid},{"metrics",expected.at("metrics")},{"cache",session_->cache_info()},
                      {"verification_seconds",seconds}};
            if (!valid) result["error"]="native certificate differs from independent C++ replay";
            return result;
        }
        throw ReplayError("unknown verifier operation");
    }
};
} // namespace
} // namespace zipper_verify

int main() {
    using namespace zipper_verify;
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    Protocol protocol;
    std::string line;
    while (std::getline(std::cin,line)) {
        try {
            auto response=protocol.dispatch(J::parse(line));
            std::cout<<response.dump(-1,' ',true)<<'\n'<<std::flush;
        } catch (const std::bad_alloc &) {
            // Fixed literals avoid allocating another JSON tree while reporting
            // resource exhaustion. No equality decision was completed.
            std::cout<<"{\"ok\":false,\"error\":\"verifier allocation failed\",\"replay_status\":\"incomplete\",\"replay_diagnostic\":{\"message\":\"verifier allocation failed\"}}\n"<<std::flush;
        } catch (const std::length_error &) {
            std::cout<<"{\"ok\":false,\"error\":\"verifier container size limit exceeded\",\"replay_status\":\"incomplete\",\"replay_diagnostic\":{\"message\":\"verifier container size limit exceeded\"}}\n"<<std::flush;
        } catch (const std::exception &error) {
            std::cout<<J({{"ok",false},{"error",error.what()}}).dump(-1,' ',true)<<'\n'<<std::flush;
        }
    }
    return std::cin.bad()?1:0;
}
