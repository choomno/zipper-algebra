# Pull requests

Target `dev` for feature and fix work unless the maintainers specify another
branch. Use an imperative title that names the change, such as “Reject
certificates with mismatched presentations.” Keep each PR focused on a single
reviewable change, with tests and documentation where they belong. Separate
unrelated formatting and mechanical refactors from behavior changes. Open a
draft while the scope or validation is incomplete.

Apply a GitHub label for the primary type of change: `enhancement` for a new
feature or capability, `bug` for a fix, `refactor` for an internal restructuring,
`dependencies` for a dependency update, or `documentation` for a docs-only change.
Update the label if the PR changes scope.

For a stack, target each dependent PR at its immediate prerequisite branch.
Explain the merge order and final integration branch. As prerequisites land,
retarget the next PR and verify that its diff still contains only the intended
changes. Split large work at boundaries that can be understood and checked
independently; keep experimental backends separate when they are optional.

## Description format

Before writing or updating the description, verify the actual PR base and head
against current refs. Inspect `git diff <base>...<head>` using the immediate base
of a stacked PR. Repeat this after rebasing, retargeting, or changing scope.
Describe the final difference from that base. Include what the PR adds, changes,
and removes; leave out intermediate designs and changes already in the base.

Use the [PR template](../.github/pull_request_template.md). Keep **Summary**,
**Changes**, and **Validation**. Include **Dependencies** and **Migration** only
when applicable.

- **Summary:** State the concrete problem and resulting behavior. For a bug, give
  the trigger and what changes for that input.
- **Dependencies:** Link unmerged prerequisite PRs, with their merge order and
  final target branch. Identify external tools or revisions required to build,
  run, or review the change.
- **Changes:** List the meaningful changes relative to the actual base. Link the
  relevant entry points, implementations, mathematical statements, tests, or
  examples beside the claims they support. Explain compatibility where the diff
  affects it.
- **Validation:** Report checks actually run, their observed results, and gaps
  that affect review. Link tests, reproduction commands, or evidence. For a
  short change, one or two sentences can be enough.
- **Migration:** State what callers or operators must do, if anything. For
  changes to certificates, catalogs, or SQLite records, explain version
  compatibility, how to preserve existing results, and how to check the result
  after migration.

Pin source links to the full PR head commit using
`https://github.com/choomno/zipper-algebra/blob/<full-head-sha>/<path>`. Add a
line anchor for a specific definition or assertion. Check the path and anchor
at that commit and refresh them after changes. When discussing removed code,
link the diff or identify the base commit explicitly. Link directly to upstream
sources for external behavior.

## Evidence for zipper-algebra changes

Use the relevant checks from [Tests and continuous integration](TESTING.md).
State the scope of the evidence in terms a reviewer can reproduce:

- For an algebraic claim, give the presentation, coefficient domain,
  hypotheses, and conclusion. Distinguish a proof or exact certificate check
  from a specialization, sampled observation, or incomplete computation.
  Identify what an independent checker verifies when the change relies on one.
- For an implementation change, name the affected entry point and report the
  relevant regression, exact replay, native, or installed-package checks.
  State which optional backends or cases were exercised if it matters.
- For a performance or census claim, identify the source and build, input
  cohort, settings and limits, and what the timing includes. Use matched inputs
  for comparisons, report refusals and timeouts that affect the result, and
  distinguish computed candidates from independently certified results. State
  the enumeration scope and how stored results are preserved or reproduced.

Keep prompts, conversation transcripts, private working notes, and
machine-specific paths out of commits, PRs, and public release artifacts.
Link small, inspectable fixtures or portable research artifacts when they are
needed for review. See [Output provenance](PROVENANCE.md) for how the software
records the identities of its inputs and outputs.
