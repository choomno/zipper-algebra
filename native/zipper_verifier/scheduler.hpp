#pragma once
#include "core.hpp"
#include <list>
#include <string>
#include <unordered_map>

namespace zipper_verify {
class Catalog;
// Independent JSON-state replay. The retained cache contains only transitions
// computed by this verifier and retains every caller context and endpoint label.
class ReplaySession {
    struct Entry { J result; std::size_t bytes; I input; };
    using CacheList = std::list<std::pair<std::string, Entry>>;
    Catalog &catalog_;
    std::size_t max_entries_, max_bytes_, bytes_ = 0;
    CacheList cache_;
    std::unordered_map<std::string, CacheList::iterator> index_;
    I input_ = 0, hits_ = 0, misses_ = 0, evictions_ = 0, oversized_ = 0;
    I cross_input_hits_ = 0, complement_hits_ = 0, peak_entries_ = 0, peak_bytes_ = 0;
    J evaluate(const J &work, const Limits &limits, bool trace);
public:
    explicit ReplaySession(Catalog &catalog, std::size_t cache_entries=262144,
                           std::size_t cache_bytes=128*1024*1024);
    J replay(I n, const Word &word, const Limits &limits=Limits{},
             I batch_size=64, bool trace=false);
    J cache_info() const;
    void clear_cache();
};
std::string continuation_key(const J &item, bool pending=true);
} // namespace zipper_verify
