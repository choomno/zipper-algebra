#pragma once
#include "typed_state.hpp"

namespace zipper_verify::typed {
// Observe every coordinate map separately, including cancelling reflections
// and zero shifts, in the same order as the literal expression materializer.
inline I work_expression_size(const Work &work) {
    try {
        struct Map { bool reflection; I value; };
        std::vector<I> shifts;
        for(auto context=work.contexts.rbegin();context!=work.contexts.rend();++context)
            shifts.push_back(context->size);
        auto count_word=[](const Word &word,const std::vector<Map> &maps) {
            if(!maps.empty()) for(const Token &token:word) {
                const bool generator=is_int(token);
                I coordinate=generator?as_i(token):std::get<Slider>(token).start;
                const I size=generator?0:std::get<Slider>(token).size;
                for(const Map &map:maps) {
                    if(generator) {
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
        if(work.world_phi) maps.push_back({true,work.rank});
        for(I shift:shifts) maps.push_back({false,shift});
        I count=0;
        auto add=[&](const Word &word,const std::vector<Map> &mapping) {
            count=checked_add(count,count_word(word,mapping));
        };
        if(const auto *single=std::get_if<Single>(&work.phase)) {
            add(single->left,maps);
            add(Word{slider(1,single->position)},maps);
            add(single->processed,maps);
        } else if(const auto *up=std::get_if<ZipUp>(&work.phase)) add(up->branch.word,maps);
        else if(const auto *down=std::get_if<ZipDown>(&work.phase)) {
            auto local_maps=maps;
            local_maps.insert(local_maps.begin(),{true,work.rank});
            add(down->branch.word,local_maps);
        } else if(const auto *layer=std::get_if<Layer>(&work.phase)) add(layer->word,maps);
        else add(std::get<Extract>(work.phase).word,maps);
        add(work.boundary_left,maps);
        add(work.boundary_right,maps);
        // A context route enters after its own shift and receives outer shifts.
        std::size_t index=0;
        for(auto context=work.contexts.rbegin();context!=work.contexts.rend();++context,++index) {
            std::vector<Map> outer;
            for(std::size_t i=index+1;i<shifts.size();++i) outer.push_back({false,shifts[i]});
            add(context->top_route,outer);
            add(context->component_word,outer);
            add(context->bottom_route,outer);
        }
        return count;
    } catch(const Error &) {
        // The materializer determines the first refusal when multiple maps fail.
        return length(work_expression(work));
    }
}
}
