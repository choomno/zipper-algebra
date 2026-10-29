"""Emit a standalone GAP check of a free-algebra cofactor identity over Q(q,t).

GAP performs the rational-function arithmetic and noncommutative expansion.
Input parsing and the defining presentation are shared with zipper-algebra.
The generated program does not compute a Groebner basis or sample parameters.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

from zipper_algebra import Cofactor, Presentation
from zipper_algebra.coefficients import SymbolicRing
from zipper_algebra.cofactor.cofactor_io import symbolic_target
from zipper_algebra.provenance import capture_provenance


def _coefficient(polynomial):
    return "(" + " + ".join(
        f"({value})*q^{i}*t^{j}" for (i, j), value in sorted(polynomial.items())
        if value) + ")" if any(polynomial.values()) else "(0)"


def _word(word):
    return "*".join(f"g[{index + 1}]" for index in word) or "One(A)"


def _polynomial(terms):
    return " + ".join(f"{_coefficient(value)}*{_word(word)}"
                      for word, value in sorted(terms.items()) if any(value.values())) or "Zero(A)"


def gap_script(cofactor, presentation, target=None, *, producer=None):
    """Return GAP source that exits 0 exactly when the supplied identity holds.

    Rational-function denominators are constructed from atom coefficient data,
    never atom display names. Words retain their original ordered occurrences.
    """
    cofactor.validate(presentation)
    target = cofactor.target if target is None else target
    if target is None or target == "":
        raise ValueError("an explicit target or a certificate target header is required")
    target_polynomial = symbolic_target(presentation, target)
    relations = presentation.relations(SymbolicRing())
    indices = cofactor.relation_indices(presentation)
    lines = ["# Exact free-algebra identity over Q(q,t), checked by GAP.",
             "# Parser and defining presentation: zipper-algebra; arithmetic: GAP.",
             "# Run: gap -A -q -r -T --quitonbreak check.g",
             "# claim: " + json.dumps({
                 "n": presentation.n, "untwist": presentation.untwist,
                 "target": str(target), "presentation_hash": presentation.content_hash(),
                 "terms": len(cofactor.terms), "generators": presentation.var_names(),
             }, sort_keys=True, ensure_ascii=True)]
    if producer is not None:
        lines.append("# producer: " + json.dumps(producer, sort_keys=True, ensure_ascii=True))
    lines += ["OnBreak := function() QUIT_GAP(1); end;;", "CheckCertificate := function()", "  local K, q, t, A, g, r, residual;",
              '  K := FunctionField(Rationals, ["q", "t"]);;',
              "  q := IndeterminatesOfFunctionField(K)[1];;",
              "  t := IndeterminatesOfFunctionField(K)[2];;",
              f"  A := FreeAssociativeAlgebraWithOne(K, {2 * (presentation.n - 1)});;",
              "  g := GeneratorsOfAlgebra(A);;",
              f"  if Length(g) <> {2 * (presentation.n - 1) + 1} or g[1] <> One(A) then",
              '    Error("Unexpected free-algebra generator ordering");', "  fi;", "  r := ["]
    for index, relation in enumerate(relations):
        suffix = "," if index + 1 < len(relations) else ""
        terms = {word: coefficient.v for word, coefficient in relation.terms.items()}
        label = json.dumps(presentation.labels()[index], ensure_ascii=True)
        lines.append("    " + _polynomial(terms) + suffix + " # " + label)
    lines += ["  ];;", "  residual := " + _polynomial(target_polynomial) + ";;"]
    for (left, _label, right, numerator, exponents), index in zip(cofactor.terms, indices, strict=True):
        denominator = "*".join(f"{_coefficient(atom.coefs)}^{power}"
                               for atom, power in zip(cofactor.atoms.atoms, exponents, strict=True)
                               if power) or "1"
        lines.append(f"  residual := residual - ({_coefficient(numerator)}/({denominator}))"
                     f"*{_word(left)}*r[{index + 1}]*{_word(right)};;")
    lines += ["  if residual = Zero(A) then", '    Print("PASS: exact cofactor identity over Q(q,t)\\n");',
              "    return true;", "  fi;", '  Print("FAIL: nonzero free-algebra residual\\n");',
              "  return false;", "end;;", "checked := CALL_WITH_CATCH(CheckCertificate, []);;",
              "if checked[1] = true and checked[2] = true then QUIT_GAP(0); fi;",
              "QUIT_GAP(1);", ""]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--n", required=True, type=int, help="presentation strand count")
    parser.add_argument("--target", help="target polynomial name; defaults to certificate header")
    parser.add_argument("--untwist", choices=("bigelow", "square-zero", "silly_question", "none"),
                        default="bigelow")
    args = parser.parse_args(argv)
    try:
        presentation = Presentation(args.n, untwist=args.untwist)
        raw = args.certificate.read_bytes()
        cofactor = Cofactor.from_text(raw.decode("utf-8"), pres=presentation)
        producer = capture_provenance()
        producer["input_certificate"] = {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}
        script = gap_script(cofactor, presentation, args.target, producer=producer)
    except (OSError, ValueError, TypeError, KeyError, UnicodeError) as error:
        parser.exit(2, f"Cannot export certificate: {error}\n")
    sys.stdout.write(script)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
