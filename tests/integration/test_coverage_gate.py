"""The release gate checks separate unrounded ratios and complete source inclusion."""
from contextlib import redirect_stdout
from copy import deepcopy
import io
from pathlib import Path
import pytest

from ci.check_coverage import check


class TestCoverageRequirements:
    @pytest.fixture(autouse=True)
    def coverage_data(self, tmp_path):
        self.root = tmp_path
        source = self.root / 'src/zipper_algebra/example.py'
        source.parent.mkdir(parents=True)
        source.write_text('value = 1\n')
        self.report = {'files': {'src/zipper_algebra/example.py': {
            'summary': {'num_statements': 1_000_000, 'covered_lines': 995_000,
                        'num_branches': 1_000_000, 'covered_branches': 995_000},
            'excluded_lines': [],
        }}}

    def run_gate(self, report):
        with redirect_stdout(io.StringIO()):
            return check(report, root=self.root)

    def test_the_exact_threshold_passes_for_both_measurements(self):
        result = self.run_gate(self.report)
        assert (result['covered_lines']) == (995000)
        assert (result['covered_branches']) == (995000)

    def test_display_rounding_cannot_admit_a_lower_line_ratio(self):
        report = deepcopy(self.report)
        report['files']['src/zipper_algebra/example.py']['summary']['covered_lines'] -= 1
        with pytest.raises(ValueError, match='before rounding'):
            self.run_gate(report)

    def test_line_coverage_cannot_compensate_for_a_lower_branch_ratio(self):
        report = deepcopy(self.report)
        summary = report['files']['src/zipper_algebra/example.py']['summary']
        summary['covered_lines'] = 1_000_000
        summary['covered_branches'] -= 1
        with pytest.raises(ValueError, match='before rounding'):
            self.run_gate(report)

    @pytest.mark.parametrize("kind", ["missing", "extra"])
    def test_missing_or_extra_source_files_are_rejected(self, kind):
        report = deepcopy(self.report)
        if kind == "missing":
            report["files"].clear()
        else:
            report["files"]["other.py"] = deepcopy(next(iter(report["files"].values())))
        with pytest.raises(ValueError, match="every maintained Python file"):
            self.run_gate(report)

    def test_excluded_source_lines_are_rejected(self):
        report = deepcopy(self.report)
        report['files']['src/zipper_algebra/example.py']['excluded_lines'] = [1]
        with pytest.raises(ValueError, match='no excluded lines'):
            self.run_gate(report)

    def test_branch_measurement_is_required(self):
        report = deepcopy(self.report)
        summary = report['files']['src/zipper_algebra/example.py']['summary']
        summary['num_branches'] = summary['covered_branches'] = 0
        with pytest.raises(ValueError, match='both required'):
            self.run_gate(report)
