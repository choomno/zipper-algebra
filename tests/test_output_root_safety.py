"""Output preparation preserves checkout visibility and custom ignore rules."""

from pathlib import Path
import subprocess
import tempfile
import unittest

from experiments.output import prepare_output_dir
from zipper_algebra.execution.artifacts import ArtifactError, prepare_output_root


class TestOutputRootSafety(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def test_checkout_root_is_rejected_before_hiding_untracked_files(self):
        subprocess.run(["git", "init", "--quiet", str(self.root)], check=True)
        (self.root / "research.txt").write_text("untracked research")
        for prepare, error_type in ((prepare_output_root, ArtifactError),
                                    (prepare_output_dir, ValueError)):
            with self.subTest(prepare=prepare.__name__):
                with self.assertRaisesRegex(error_type, "Git checkout root"):
                    prepare(self.root)
                self.assertFalse((self.root / ".gitignore").exists())
                status = subprocess.check_output(
                    ["git", "-C", str(self.root), "status", "--porcelain"]
                ).decode()
                self.assertIn("?? research.txt", status)

    def test_worktree_git_file_is_rejected_with_existing_policy_preserved(self):
        (self.root / ".git").write_text("gitdir: /unused/worktree/metadata\n")
        policy = b"*.tmp\n"
        (self.root / ".gitignore").write_bytes(policy)
        for prepare in (prepare_output_root, prepare_output_dir):
            with self.subTest(prepare=prepare.__name__):
                with self.assertRaisesRegex(ValueError, "Git checkout root"):
                    prepare(self.root)
                self.assertEqual((self.root / ".gitignore").read_bytes(), policy)
                self.assertEqual(sorted(path.name for path in self.root.iterdir()),
                                 [".git", ".gitignore"])

    def test_nested_output_is_ignored_and_existing_custom_policy_survives(self):
        subprocess.run(["git", "init", "--quiet", str(self.root)], check=True)
        for index, prepare in enumerate((prepare_output_root, prepare_output_dir)):
            with self.subTest(prepare=prepare.__name__):
                output = self.root / "runs" / str(index)
                prepare(output)
                self.assertEqual((output / ".gitignore").read_text(), "*\n!.gitignore\n")
                (output / "certificate.txt").write_text("generated")
                ignored = subprocess.run(
                    ["git", "-C", str(self.root), "check-ignore", "--quiet",
                     str(output / "certificate.txt")], check=False,
                )
                self.assertEqual(ignored.returncode, 0)
                custom = b"# custom output rules\n*.txt\n"
                (output / ".gitignore").write_bytes(custom)
                prepare(output)
                self.assertEqual((output / ".gitignore").read_bytes(), custom)

    def test_existing_source_directory_preserves_visibility_of_new_files(self):
        subprocess.run(["git", "init", "--quiet", str(self.root)], check=True)
        for index, prepare in enumerate((prepare_output_root, prepare_output_dir)):
            with self.subTest(prepare=prepare.__name__):
                source = self.root / "examples" / str(index)
                source.mkdir(parents=True)
                (source / "README.md").write_text("example source")
                (source / "certificate.txt").write_text("example certificate")
                prepare(source)
                self.assertFalse((source / ".gitignore").exists())
                status = subprocess.check_output(
                    ["git", "-C", str(self.root), "status", "--porcelain",
                     "--untracked-files=all"],
                ).decode()
                self.assertIn(f"?? examples/{index}/README.md", status)
                self.assertIn(f"?? examples/{index}/certificate.txt", status)
                policy = b"*.tmp\n"
                (source / ".gitignore").write_bytes(policy)
                prepare(source)
                self.assertEqual((source / ".gitignore").read_bytes(), policy)

    def test_existing_generated_directories_get_missing_markers(self):
        subprocess.run(["git", "init", "--quiet", str(self.root)], check=True)
        for index, prepare in enumerate((prepare_output_root, prepare_output_dir)):
            for name in ("runs", ".research-data", ".cache"):
                with self.subTest(prepare=prepare.__name__, directory=name):
                    output = self.root / name / str(index)
                    output.mkdir(parents=True)
                    prepare(output)
                    self.assertEqual((output / ".gitignore").read_text(),
                                     "*\n!.gitignore\n")
                    (output / "generated.txt").write_text("generated")
                    ignored = subprocess.run(
                        ["git", "-C", str(self.root), "check-ignore", "--quiet",
                         str(output / "generated.txt")], check=False,
                    )
                    self.assertEqual(ignored.returncode, 0)

    def test_existing_output_outside_git_gets_ignore_marker(self):
        for index, prepare in enumerate((prepare_output_root, prepare_output_dir)):
            with self.subTest(prepare=prepare.__name__):
                output = self.root / str(index)
                output.mkdir()
                prepare(output)
                self.assertEqual((output / ".gitignore").read_text(),
                                 "*\n!.gitignore\n")

    def test_source_directory_in_worktree_preserves_missing_ignore_policy(self):
        (self.root / ".git").write_text("gitdir: /unused/worktree/metadata\n")
        for index, prepare in enumerate((prepare_output_root, prepare_output_dir)):
            with self.subTest(prepare=prepare.__name__):
                output = self.root / "examples" / str(index)
                output.mkdir(parents=True)
                prepare(output)
                self.assertFalse((output / ".gitignore").exists())

    def test_empty_git_named_directory_does_not_make_parent_a_checkout(self):
        (self.root / ".git").mkdir()
        for index, prepare in enumerate((prepare_output_root, prepare_output_dir)):
            with self.subTest(prepare=prepare.__name__):
                output = self.root / str(index)
                output.mkdir()
                prepare(output)
                self.assertEqual((output / ".gitignore").read_text(),
                                 "*\n!.gitignore\n")
