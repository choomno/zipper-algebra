#include "typed_canonical.hpp"
namespace zipper_verify::typed {
namespace {
struct Writer {
    std::string text;
    struct Object {
        Writer &writer;
        bool first = true;
        explicit Object(Writer &w) : writer(w) {
            writer.text += '{';
        }
        ~Object() {
            writer.text += '}';
        }
        template <class F> void field(const char *name, F emit) {
            if (!first)
                writer.text += ',';
            first = false;
            writer.string(name);
            writer.text += ':';
            emit();
        }
    };
    void string(const std::string &s) {
        text += J(s).dump(-1, ' ', true);
    }
    void integer(I n) {
        text += std::to_string(n);
    }
    void boolean(bool b) {
        text += b ? "true" : "false";
    }
    template <class V, class F> void array(const V &values, F emit) {
        text += '[';
        bool first = true;
        for (const auto &value : values) {
            if (!first)
                text += ',';
            first = false;
            emit(value);
        }
        text += ']';
    }
    void labels(const Labels &v) {
        array(v, [&](I i) { integer(i); });
    }
    void word(const Word &v) {
        array(v, [&](const Token &token) {
            if (const auto *i = std::get_if<I>(&token))
                integer(*i);
            else {
                const auto &s = std::get<Slider>(token);
                text += "[\"slider\",";
                integer(s.size);
                text += ',';
                integer(s.start);
                text += ',';
                boolean(s.mirrored);
                text += ',';
                boolean(s.reversed);
                text += ',';
                boolean(s.barred);
                text += ']';
            }
        });
    }
    void scalar(const Laurent &v) {
        text += '[';
        bool first = true;
        for (const auto &[q, c] : v)
            if (c != 0) {
                if (!first)
                    text += ',';
                first = false;
                text += '[';
                integer(q);
                text += ',';
                text += c.get_str();
                text += ']';
            }
        text += ']';
    }
    void record(const ComponentRecord &r) {
        Object o(*this);
        o.field("bottom_ports", [&] { labels(r.bottom_ports); });
        o.field("bottom_route", [&] { word(r.bottom_route); });
        o.field("canonical_bottom_route", [&] { word(r.canonical_bottom_route); });
        o.field("canonical_parities", [&] {
            Object p(*this);
            p.field("lower_inverse", [&] { labels(r.lower_inverse_parities); });
            p.field("upper", [&] { labels(r.upper_parities); });
        });
        o.field("canonical_top_route", [&] { word(r.canonical_top_route); });
        if (r.inherited_bottom_route)
            o.field("inherited_bottom_route", [&] { word(*r.inherited_bottom_route); });
        if (r.inherited_top_route)
            o.field("inherited_top_route", [&] { word(*r.inherited_top_route); });
        if (r.input_hierarchy)
            o.field("input_hierarchy", [&] { labels(*r.input_hierarchy); });
        o.field("orientation", [&] { string(r.mirrored ? "mirrored" : "normal"); });
        o.field("original_bottom_ports", [&] { labels(r.original_bottom_ports); });
        o.field("original_top_ports", [&] { labels(r.original_top_ports); });
        if (r.output_hierarchy)
            o.field("output_hierarchy", [&] { labels(*r.output_hierarchy); });
        o.field("top_ports", [&] { labels(r.top_ports); });
        o.field("top_route", [&] { word(r.top_route); });
        o.field("twist_powers", [&] { labels(r.twist_powers); });
    }
    void context(const Context &c) {
        Object o(*this);
        o.field("bottom_route", [&] { word(c.bottom_route); });
        o.field("component_word", [&] { word(c.component_word); });
        o.field("rank", [&] { integer(c.rank); });
        o.field("record", [&] { record(c.record); });
        o.field("size", [&] { integer(c.size); });
        o.field("top_route", [&] { word(c.top_route); });
    }
    void contexts(const std::vector<Context> &v) {
        array(v, [&](const Context &c) { context(c); });
    }
    void factors(const std::vector<Factor> &v) {
        array(v, [&](const Factor &f) {
            Object o(*this);
            o.field("growth_size", [&] { integer(f.growth_size); });
            o.field("word", [&] { word(f.word); });
        });
    }
    void initial(const InitialSlider &s) {
        Object o(*this);
        o.field("input_power", [&] { integer(s.input_power); });
        o.field("j", [&] { integer(s.j); });
        o.field("k", [&] { integer(s.k); });
        o.field("mirror_policy", [&] { string(s.mirror_policy); });
        o.field("mirrored", [&] { boolean(s.mirrored); });
        o.field("output_power", [&] { integer(s.output_power); });
        o.field("processed_length", [&] { integer(s.processed_length); });
    }
    void event(const ProcessedEvent &e) {
        Object o(*this);
        o.field("mirror_policy", [&] { string(e.mirror_policy); });
        o.field("parameters", [&] {
            Object p(*this);
            for (const auto &[k, v] : e.parameters)
                p.field(k.c_str(), [&] { integer(v); });
        });
        o.field("rhs_index", [&] { integer(e.rhs_index); });
        o.field("rule_id", [&] { string(e.rule_id); });
    }
    void control(const Control &c) {
        Object o(*this);
        if (c.phase == ControlPhase::UnsupportedFrame) {
            o.field("message", [&] { string(c.message); });
            o.field("phase", [&] { string("unsupported_frame"); });
            o.field("rule_id", [&] { string(c.rule_id); });
            return;
        }
        if (c.phase == ControlPhase::TurnaroundComplete) {
            o.field("exterior_length", [&] { integer(c.exterior_length); });
            o.field("i", [&] { integer(c.i); });
            o.field("input_power", [&] { integer(c.input_power); });
            o.field("j", [&] { integer(c.j); });
            o.field("mirror_policy", [&] { string(c.mirror_policy); });
            o.field("mirrored", [&] { boolean(c.mirrored); });
            o.field("output_power", [&] { integer(c.output_power); });
            o.field("payload_length", [&] { integer(c.payload_length); });
            o.field("phase", [&] { string("turnaround_complete"); });
            o.field("processed_length", [&] { integer(c.processed_length); });
            o.field("routing_length", [&] { integer(c.routing_length); });
            o.field("tail_length", [&] { integer(c.tail_length); });
            o.field("transported_length", [&] { integer(c.transported_length); });
            return;
        }
        const bool framed =
            c.phase == ControlPhase::FramedUp || c.phase == ControlPhase::FramedComplete;
        o.field("i", [&] { integer(c.i); });
        if (framed && c.initial_slider)
            o.field("initial_slider", [&] { initial(*c.initial_slider); });
        if (framed)
            o.field("input_power", [&] { integer(c.input_power); });
        o.field("j", [&] { integer(c.j); });
        o.field("k", [&] { integer(c.k); });
        o.field("left_length", [&] { integer(c.left_length); });
        if (framed) {
            o.field("mirror_policy", [&] { string(c.mirror_policy); });
            o.field("mirrored", [&] { boolean(c.mirrored); });
            o.field("output_payload", [&] { word(c.output_payload); });
            o.field("output_power", [&] { integer(c.output_power); });
        }
        o.field("phase", [&] {
            switch (c.phase) {
            case ControlPhase::Up:
                string("zip_up");
                break;
            case ControlPhase::UpComplete:
                string("upper_pass_complete");
                break;
            case ControlPhase::FramedUp:
                string("framed_zip_up");
                break;
            default:
                string("framed_upper_pass_complete");
            }
        });
        if (framed)
            o.field("processed_events",
                    [&] { array(c.processed_events, [&](const ProcessedEvent &e) { event(e); }); });
        o.field("processed_length", [&] { integer(c.processed_length); });
        o.field("tail_length", [&] { integer(c.tail_length); });
    }
    void branch(const Branch &b, bool skip_scalar) {
        Object o(*this);
        if (!skip_scalar)
            o.field("coefficient", [&] { scalar(b.coefficient); });
        o.field("control", [&] { control(b.control); });
        o.field("word", [&] { word(b.word); });
    }
    void work(const Work &w, bool skip_scalar, bool pending) {
        const auto *s = std::get_if<Single>(&w.phase);
        const auto *u = std::get_if<ZipUp>(&w.phase);
        const auto *d = std::get_if<ZipDown>(&w.phase);
        const auto *l = std::get_if<Layer>(&w.phase);
        Object o(*this);
        o.field("bottom_labels", [&] { labels(w.bottom_labels); });
        o.field("boundary_left", [&] { word(w.boundary_left); });
        o.field("boundary_left_factors", [&] { factors(w.boundary_left_factors); });
        o.field("boundary_right", [&] { word(w.boundary_right); });
        o.field("boundary_right_factors", [&] { factors(w.boundary_right_factors); });
        if (u)
            o.field("branch", [&] { branch(u->branch, pending); });
        else if (d)
            o.field("branch", [&] { branch(d->branch, pending); });
        if (!skip_scalar)
            o.field("coefficient", [&] { scalar(w.coefficient); });
        if (l)
            o.field("completed_sides", [&] { integer(l->completed_sides); });
        o.field("contexts", [&] { contexts(w.contexts); });
        if (!pending)
            o.field("depth", [&] { integer(w.depth); });
        if (s)
            o.field("left", [&] { word(s->left); });
        if (l)
            o.field("payload_length", [&] { integer(l->payload_length); });
        o.field("phase", [&] { string(phase_name(w)); });
        if (s) {
            o.field("position", [&] { integer(s->position); });
            o.field("processed", [&] { word(s->processed); });
        }
        o.field("rank", [&] { integer(w.rank); });
        o.field("top_labels", [&] { labels(w.top_labels); });
        if (l)
            o.field("word", [&] { word(l->word); });
        else if (const auto *e = std::get_if<Extract>(&w.phase))
            o.field("word", [&] { word(e->word); });
        o.field("world_phi", [&] { boolean(w.world_phi); });
    }
};
} // namespace
std::string canonical_work(const Work &w) {
    Writer out;
    out.work(w, false, false);
    return std::move(out.text);
}
std::string continuation_key(const Work &w) {
    Writer out;
    out.work(w, true, true);
    return std::move(out.text);
}
std::string work_collection_key(const Work &w) {
    Writer out;
    out.work(w, true, false);
    return std::move(out.text);
}
std::string output_collection_key(const Output &v) {
    Writer out;
    {
        Writer::Object o(out);
        o.field("components", [&] { out.contexts(v.components); });
        o.field("word", [&] { out.word(v.word); });
    }
    return std::move(out.text);
}
} // namespace zipper_verify::typed
