#pragma once
#include "algorithm.hpp"

namespace zipper_verify {
// Observe the same coordinate maps as work_expression, without constructing
// its concatenations. Keep each arithmetic operation: a composed map could
// conceal an intermediate overflow, including two cancelling reflections.
inline I work_expression_size(const J &work) {
    try {
    struct Map { bool reflection; I value; };
    const I rank=as_i(work.at("rank"));
    const std::string phase=work.at("phase").get<std::string>();
    const bool reflected=work.at("world_phi").get<bool>();
    std::vector<I> shifts;
    for (auto context=work.at("contexts").rbegin(); context!=work.at("contexts").rend(); ++context)
        shifts.push_back(as_i(context->at("size")));
    auto count_word=[](const J &word,const std::vector<Map> &maps) {
        if (!word.is_array()) throw ReplayError("word must be an array");
        if (!maps.empty()) for (const J &token:word) {
            const bool generator=is_int(token);
            if (!generator && !is_slider(token)) throw ReplayError("word contains an invalid braid or slider token");
            I coordinate=as_i(generator?token:token.at(2));
            const I size=generator?0:as_i(token.at(1));
            for (const Map &map:maps) {
                if (generator) {
                    const I magnitude=map.reflection?checked_sub(map.value,abs_i(coordinate)):
                                                          checked_add(abs_i(coordinate),map.value);
                    coordinate=checked_mul(coordinate>0?1:-1,magnitude);
                } else coordinate=map.reflection?checked_add(checked_sub(checked_sub(map.value,size),coordinate),2):
                                                     checked_add(coordinate,map.value);
            }
        }
        return length(word);
    };
    std::vector<Map> maps;
    if (reflected) maps.push_back({true,rank});
    for (I shift:shifts) maps.push_back({false,shift});
    I count=0;
    auto add=[&](const J &word,const std::vector<Map> &mapping) { count=checked_add(count,count_word(word,mapping)); };
    if (phase=="single") {
        add(work.at("left"),maps);
        add(J::array({slider(1,as_i(work.at("position")))}),maps);
        add(work.at("processed"),maps);
    } else if (phase=="zip_up" || phase=="zip_down") {
        auto local_maps=maps;
        if (phase=="zip_down") local_maps.insert(local_maps.begin(),{true,rank});
        add(work.at("branch").at("word"),local_maps);
    } else add(work.at("word"),maps);
    add(work.at("boundary_left"),maps);
    add(work.at("boundary_right"),maps);
    // Context routes enter after their own shift, then receive outer shifts.
    std::size_t index=0;
    for (auto context=work.at("contexts").rbegin(); context!=work.at("contexts").rend(); ++context,++index) {
        std::vector<Map> outer;
        for (std::size_t i=index+1;i<shifts.size();++i) outer.push_back({false,shifts[i]});
        for (const char *field:{"top_route","component_word","bottom_route"}) add(context->at(field),outer);
    }
    return count;
    } catch (const Error &) {
        // Preserve which original check refuses first when a malformed state
        // violates several conditions. The cold refusal path may allocate.
        return length(work_expression(work));
    } catch (const J::exception &) {
        return length(work_expression(work));
    }
}
} // namespace zipper_verify
