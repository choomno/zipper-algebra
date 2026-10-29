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
