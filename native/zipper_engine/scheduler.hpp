#pragma once
#include "algorithm.hpp"
#include <memory>
namespace zipper {
class EngineSession {
  public:
    // The referenced immutable catalog must outlive the session. Successful
    // unit-coefficient transitions are retained between complete input runs.
    explicit EngineSession(Catalog &catalog, size_t cache_entries = 65536,
                           size_t cache_bytes = 256 * 1024 * 1024);
    ~EngineSession();
    EngineSession(const EngineSession &) = delete;
    EngineSession &operator=(const EngineSession &) = delete;
    Json reduce(const Work &initial, const Limits &limits, bool trace = false,
                size_t batch_size = 64);
    Json cache_info() const;
    Json last_telemetry() const;
    void record_timing(const char *phase, double seconds);
    void clear_cache();
    Transition transition(const Work &work, const Limits &limits, bool trace);

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
// All mutable frontiers, Laurent sums, and transition states remain native.
// JSON is produced only for optional trace recording and the returned boundary.
Json reduce_collected(const Work &initial, Catalog &catalog, const Limits &limits,
                      bool trace = false, size_t batch_size = 64, EngineSession *session = nullptr);
} // namespace zipper
