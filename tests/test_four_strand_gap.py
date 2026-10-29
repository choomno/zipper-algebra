"""Fast parser and proof-mutation checks for the optional independent GAP route.

The full 186-dimensional research computation is deliberately run locally.
"""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess

import pytest

from zipper_algebra import Presentation
from zipper_algebra.coefficients import SymbolicRing
from zipper_algebra.verification.groebner import GroebnerCertificate, ProofTerm


ROOT = Path(__file__).resolve().parents[1]
CHECKER = ROOT / "examples/verify_four_strand_square_zero.g"


@pytest.fixture(scope="module")
def gap():
    binary = shutil.which("gap")
    if binary is None:
        pytest.skip("GAP is an optional independent checker")
    result = subprocess.run(
        [binary, "-A", "-q", "-r", "--quitonbreak", "-c",
         'if LoadPackage("gbnp") <> true then QUIT_GAP(1); fi; QUIT_GAP(0);'],
        text=True, capture_output=True, timeout=15, check=False,
    )
    if result.returncode:
        pytest.skip("GBNP is an optional independent checker")
    return binary


def run_gap(gap, code):
    return subprocess.run(
        [gap, "-A", "-q", "-r", "--quitonbreak", "-c",
         "SquareZeroLibraryOnly:=true;;", str(CHECKER), "-c",
         'K:=FunctionField(Rationals,["q","t"]);;' + code + 'QUIT_GAP(0);'],
        text=True, capture_output=True, timeout=15, check=False,
    )


def tiny_proof(*, coefficient="1", reference="2", kind="1"):
    return (
        "META|31|1|2\n"
        + "".join(f"I|{i}|1|1,|0|\n" for i in range(1, 32))
        + "G|1|1|1,|0|\nG|2|1|1,|0|\n"
        + f"D|1|{coefficient}||{kind}|{reference}|\nD|2|1||0|1|\n"
        + "END|31|1|2\n"
    )


def check_data(gap, tmp_path, text, *, verify=True):
    path = tmp_path / "proof.txt"
    path.write_text(text)
    code = f"data:=SquareZero.read({json.dumps(str(path))},K);;"
    if verify:
        code += 'Print("VERIFIED|",SquareZero.verifyIdentities(data,K),"\\n");'
    return run_gap(gap, code)


def test_forward_reference_dag_is_checked_without_expansion(gap, tmp_path):
    result = check_data(gap, tmp_path, tiny_proof())
    assert result.returncode == 0, result.stdout + result.stderr
    assert "VERIFIED|2" in result.stdout


@pytest.mark.parametrize(("change", "message"), [
    ({"coefficient": "-1"}, "local identity 1"),
    ({"reference": "1", "coefficient": "0"}, "cyclic derivation references"),
    ({"reference": "3", "coefficient": "0"}, "derivation reference outside bounds"),
    ({"reference": "0", "coefficient": "0"}, "derivation reference outside bounds"),
    ({"kind": "2", "coefficient": "0"}, "derivation reference outside bounds"),
    ({"kind": "0", "reference": "31"}, "derivation reference outside bounds"),
    ({"kind": "0", "reference": "31", "coefficient": "0"}, "derivation reference outside bounds"),
])
def test_mutated_identities_and_raw_zero_coefficient_references_fail(gap, tmp_path, change, message):
    result = check_data(gap, tmp_path, tiny_proof(**change))
    assert result.returncode == 1
    assert message in result.stderr
    assert "VERIFIED" not in result.stdout


@pytest.mark.parametrize("change", [
    lambda s: s.replace("END|31|1|2\n", ""),
    lambda s: s.replace("END|31|1|2", "END|31|1|3"),
    lambda s: s + "G|1|1|1,|0|\n",
    lambda s: s.replace("G|1|1|1,|0|", "G|1|1|7,|0|"),
    lambda s: s.replace("G|1|1|1,|0|", "G|1|1|1,,|0|"),
    lambda s: s.replace("G|1|1|1,|0|\n", ""),
    lambda s: s.replace("D|1|1||1|2|\n", ""),
])
def test_malformed_missing_and_truncated_records_fail(gap, tmp_path, change):
    result = check_data(gap, tmp_path, change(tiny_proof()))
    assert result.returncode == 1
    assert "VERIFIED" not in result.stdout


def test_scalar_parser_cancels_exact_rational_functions(gap):
    result = run_gap(gap, 'z:=IndeterminatesOfFunctionField(K);;'
                     'if SquareZero.scalar("(q^2-t^2)/(q-t)",K) <> z[1]+z[2] then QUIT_GAP(1); fi;'
                     'if SquareZero.scalar("-q^2+2*t^-1",K) <> -z[1]^2+2/z[2] then QUIT_GAP(1); fi;')
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("text", ["1/0", "(q+t", "q;QUIT_GAP(0)", "q^t", "q**t"])
def test_scalar_parser_rejects_undefined_or_executable_text(gap, text):
    result = run_gap(gap, f"SquareZero.scalar({json.dumps(text)},K);")
    assert result.returncode == 1
    assert "FAIL:" in result.stderr


def test_independent_presentation_matches_all_31_symbolic_inputs(gap, tmp_path):
    lines = ["META|31|1|1"]
    for index, poly in enumerate(Presentation(4, "square-zero").build(SymbolicRing())[1], 1):
        for word, scalar in poly.terms.items():
            coefficient = "+".join(f"({c})*q^{q}*t^{t}" for (q, t), c in scalar.v.items())
            lines.append(f"I|{index}|{coefficient}|{''.join(f'{i},' for i in word)}|0|")
    lines += ["G|1|1|1,|0|", "D|1|1||0|1|", "END|31|1|1"]
    path = tmp_path / "inputs.txt"
    path.write_text("\n".join(lines) + "\n")
    result = run_gap(gap, f"data:=SquareZero.read({json.dumps(str(path))},K);;"
                     'if data.inputs <> SquareZero.presentation(K).inputs then QUIT_GAP(1); fi;')
    assert result.returncode == 0, result.stdout + result.stderr


def test_actual_producer_framed_proof_handoff_preserves_generic_fractions(gap, tmp_path):
    spec = importlib.util.spec_from_file_location("square_zero_example", ROOT / "examples/four_strand_square_zero.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    full = module.presentation()
    field = full.algebra.field
    q, t = (field.parameter(name) for name in ("q", "t"))
    coefficient = (q*q-t*t)/(2*(q-t))
    certificate = GroebnerCertificate(full.content_hash(), (coefficient*full.relations[0],),
                                     ((ProofTerm(coefficient, (), "input", 0, ()),),))
    result = check_data(gap, tmp_path, module.framed_proof(full, certificate, 1))
    assert result.returncode == 0, result.stdout + result.stderr
    assert "VERIFIED|1" in result.stdout


@pytest.mark.parametrize(("change", "message"), [
    ("", None),
    ("relations:=[x];;", "original relation has nonzero normal form"),
    ("basis:=[[[[1,1],[2]],[u,-u]],[[[1,2]],[u]]];;relations:=basis;;",
     "candidate is not a Groebner basis"),
    ("dimension:=3;;", "unexpected standard-word count"),
    ("x3:=x;;", "X3 is nonzero"),
    ("w:=x2;;", "w is zero"),
    ("w:=[[[]],[u]];;", "w*b*w is nonzero"),
])
def test_quotient_checks_on_tiny_square_zero_algebra(gap, change, message):
    # Q(q,t)<x>/(x^2) exercises the complete structural path in milliseconds.
    code = ('u:=One(K);;x:=[[[1]],[u]];;x2:=[[[1,1]],[u]];;'
            'basis:=[x2];;relations:=[x2];;x3:=x2;;w:=x;;dimension:=2;;'
            'multiply:=SquareZero.presentation(K).times;;' + change +
            'SquareZero.verifyQuotient(basis,relations,x3,w,dimension,1,multiply);;')
    result = run_gap(gap, code)
    assert result.returncode == int(message is not None), result.stdout + result.stderr
    if message is None:
        assert "PASS|complete_groebner_basis|1" in result.stdout
        assert "PASS|dimension_over_Q(q,t)|2" in result.stdout
        assert "square_zero_ideal_tests|2" in result.stdout
    else:
        assert message in result.stderr
