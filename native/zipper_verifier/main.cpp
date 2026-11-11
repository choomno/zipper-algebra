#include "algorithm.hpp"
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
        if (op=="transition") {
            bool trace=request.contains("trace_enabled")?bool_field(request,"trace_enabled"):
                       request.contains("trace")?bool_field(request,"trace"):false;
            return algorithm_transition(normalize_coefficients(request.at("work")),*catalog_,
                                        Limits(request.at("limits")),trace);
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
        } catch (const std::exception &error) {
            std::cout<<J({{"ok",false},{"error",error.what()}}).dump(-1,' ',true)<<'\n'<<std::flush;
        }
    }
    return std::cin.bad()?1:0;
}
