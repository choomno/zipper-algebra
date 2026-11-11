#include "scheduler.hpp"
#include <iostream>
#include <memory>
#include <limits>
#include <string>

namespace {
zipper::Int protocol_int(const zipper::Json &value) {
    if (!value.is_number_integer() ||
        (value.is_number_unsigned() &&
         value.get<std::uint64_t>() >
             static_cast<std::uint64_t>(std::numeric_limits<zipper::Int>::max())))
        throw zipper::EngineError("Protocol integer must fit signed 64-bit coordinates");
    return value.get<zipper::Int>();
}
zipper::Braid protocol_braid(const zipper::Json &value) {
    if (!value.is_array())
        throw zipper::EngineError("Protocol word must be an array");
    zipper::Braid result;
    for (const auto &token : value)
        result.push_back(protocol_int(token));
    return result;
}
} // namespace

// One JSON record per request/result; no Python calls occur inside a reduction.
int main() {
    using namespace zipper;
    constexpr const char *abi = "zipper-native-engine-v1";
    std::unique_ptr<Catalog> catalog;
    std::unique_ptr<EngineSession> session;
    std::string line;
    while (std::getline(std::cin, line)) {
        Json response;
        try {
            const Json request = Json::parse(line);
            const auto op = request.at("op").get<std::string>();
            if (op == "initialize") {
                if (catalog)
                    throw EngineError("Native engine is already initialized");
                if (request.at("abi") != abi)
                    throw EngineError("Native ABI differs");
                catalog = std::make_unique<Catalog>(request.at("catalog"));
                set_reference_sources(request.at("reference_sources"),
                                      request.at("catalog_sha256").get<std::string>());
                if (protocol_int(request.at("cache_entries")) < 0 ||
                    protocol_int(request.at("cache_bytes")) < 0)
                    throw EngineError("Cache budgets must be nonnegative");
                session = std::make_unique<EngineSession>(
                    *catalog, static_cast<size_t>(protocol_int(request.at("cache_entries"))),
                    static_cast<size_t>(protocol_int(request.at("cache_bytes"))));
                response = {{"abi", abi}, {"ready", true}};
            } else if (op == "reduce") {
                if (!session)
                    throw EngineError("Native engine is not initialized");
                const auto n = protocol_int(request.at("n"));
                const auto word = protocol_braid(request.at("word"));
                const auto limits = limits_from_json(request.at("limits"));
                const bool trace = request.at("trace").get<bool>();
                const auto batch = protocol_int(request.at("batch_size"));
                if (n < 1 || batch < 1)
                    throw EngineError("Positive rank and batch size required");
                validate_word(word_from_braid(word), n);
                response = session->reduce(initial_work(n, word), limits, trace,
                                           static_cast<size_t>(batch));
                response["native_runtime"] = session->last_telemetry();
                response["native_runtime"]["cache"] = session->cache_info();
            } else if (op == "cache_info" || op == "clear_cache") {
                if (!session)
                    throw EngineError("Native engine is not initialized");
                if (op == "clear_cache")
                    session->clear_cache();
                response = session->cache_info();
            } else if (op == "transition") {
                if (!session)
                    throw EngineError("Native engine is not initialized");
                // Diagnostic requests preserve complete unit evidence even when
                // their full trace is disabled; ordinary scheduling keeps it compact.
                response = transition_json(algorithm_transition(
                    work_from_json(request.at("work")), *catalog,
                    limits_from_json(request.at("limits")), request.value("trace", false)));
            } else if (op == "turnaround") {
                if (!session)
                    throw EngineError("Native engine is not initialized");
                const auto branch = branch_from_json(request.at("branch"));
                const auto n = protocol_int(request.at("n"));
                const auto limits = limits_from_json(request.at("limits"));
                if (n < 1)
                    throw EngineError("Positive rank required");
                validate_word(branch.word, n);
                check_branch(branch, limits);
                const auto result =
                    turnaround(branch, n, *catalog, limits, request.value("trace", false));
                response = {{"branch", branch_json(result.first)}, {"certificate", result.second}};
            } else if (op == "upward_step" || op == "framed_step") {
                if (!session)
                    throw EngineError("Native engine is not initialized");
                const auto branch = branch_from_json(request.at("branch"));
                const auto n = protocol_int(request.at("n"));
                const auto limits = limits_from_json(request.at("limits"));
                if (n < 1)
                    throw EngineError("Positive rank required");
                validate_word(branch.word, n);
                check_branch(branch, limits);
                const auto step =
                    op == "upward_step"
                        ? upward_step(branch, n, *catalog, limits, request.value("trace", false))
                        : framed_step(branch, n, *catalog, limits, request.value("trace", false));
                Json branches = Json::array();
                for (const auto &child : step.branches)
                    branches.push_back(branch_json(child));
                response = {
                    {"branches", branches}, {"rule_id", step.rule_id}, {"evidence", step.evidence}};
            } else {
                throw EngineError("Unknown native engine operation");
            }
        } catch (const std::exception &error) {
            response = {{"protocol_error", error.what()}};
        }
        std::cout << response.dump() << '\n' << std::flush;
        if (!std::cout)
            return 2;
    }
    return 0;
}
