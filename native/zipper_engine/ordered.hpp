#pragma once

#include <list>
#include <string>
#include <unordered_map>
#include <utility>

namespace zipper {

// Stable insertion order with O(1) lookup and erasure. Updating an existing
// key preserves its position; deleting and reinserting places it at the end.
// Those details are part of the collected executor's scheduling semantics.
template <class Value> class OrderedMap {
  public:
    using Entry = std::pair<std::string, Value>;
    using Storage = std::list<Entry>;
    using iterator = typename Storage::iterator;
    using const_iterator = typename Storage::const_iterator;

    OrderedMap() = default;
    OrderedMap(const OrderedMap &other) {
        for (const auto &entry : other.items_)
            set(entry.first, entry.second);
    }
    OrderedMap &operator=(const OrderedMap &other) {
        if (this != &other) {
            clear();
            for (const auto &entry : other.items_)
                set(entry.first, entry.second);
        }
        return *this;
    }
    OrderedMap(OrderedMap &&) noexcept = default;
    OrderedMap &operator=(OrderedMap &&) noexcept = default;

    bool empty() const {
        return items_.empty();
    }
    std::size_t size() const {
        return items_.size();
    }
    iterator begin() {
        return items_.begin();
    }
    iterator end() {
        return items_.end();
    }
    const_iterator begin() const {
        return items_.begin();
    }
    const_iterator end() const {
        return items_.end();
    }
    bool contains(const std::string &key) const {
        return index_.count(key) != 0;
    }
    Value *get(const std::string &key) {
        auto found = index_.find(key);
        return found == index_.end() ? nullptr : &found->second->second;
    }
    const Value *get(const std::string &key) const {
        auto found = index_.find(key);
        return found == index_.end() ? nullptr : &found->second->second;
    }
    // Return and promote an existing node without allocating or invalidating
    // either the index iterator or pointers to the stored value.
    Value *touch(const std::string &key) {
        auto found = index_.find(key);
        if (found == index_.end())
            return nullptr;
        items_.splice(items_.end(), items_, found->second);
        return &found->second->second;
    }
    void set(std::string key, Value value) {
        auto found = index_.find(key);
        if (found != index_.end()) {
            found->second->second = std::move(value);
            return;
        }
        items_.emplace_back(std::move(key), std::move(value));
        auto last = std::prev(items_.end());
        try {
            index_.emplace(last->first, last);
        } catch (...) {
            // Keep list/index membership identical if hashing, key copying,
            // or allocation fails after the list node has been appended.
            items_.pop_back();
            throw;
        }
    }
    void erase(const std::string &key) {
        auto found = index_.find(key);
        if (found == index_.end())
            return;
        items_.erase(found->second);
        index_.erase(found);
    }
    void clear() {
        index_.clear();
        items_.clear();
    }

  private:
    Storage items_;
    std::unordered_map<std::string, iterator> index_;
};

} // namespace zipper
