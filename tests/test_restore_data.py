"""Exercise optional data restoration using small, local archives."""

import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest
from unittest.mock import patch

from experiments.restore_data import RestoreError, _verify_archive, restore_data


class TestRestoreData(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.output = self.root / "restored"
        self.manifest = self.root / "manifest.json"
        self.archive = self.root / "archive.tar.gz"
        self.contents = {"ms_inputs/example.ms": b"x*y - y*x\n",
                         "logs_and_traces/old.txt": b"historical output\n"}
        self.document = {
            "schema": "zipper-research-archive-v1",
            "source_revision": "a" * 40,
            "source_archive_url": "https://example.invalid/archive/" + "a" * 40 + ".tar.gz",
            "groups": {"inputs": "inputs", "x3-history": "history"},
            "files": [self.entry(name, data, group) for (name, data), group in
                      zip(self.contents.items(), ("inputs", "x3-history"))],
        }
        self.write_manifest()

    @staticmethod
    def entry(name, data, group):
        return {"path": name, "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(), "group": group}

    def write_manifest(self):
        self.manifest.write_text(json.dumps(self.document), encoding="utf-8")

    def make_archive(self, contents=None, *, prefix="", extra=()):
        contents = self.contents if contents is None else contents
        with tarfile.open(self.archive, "w:gz") as bundle:
            for name, data in contents.items():
                member = tarfile.TarInfo(prefix + name)
                member.size = len(data)
                bundle.addfile(member, io.BytesIO(data))
            for member in extra:
                bundle.addfile(member, io.BytesIO(b""))

    def restore(self, **kwargs):
        return restore_data(self.output, manifest=self.manifest, archive=self.archive, **kwargs)

    def assert_only_ignore_marker(self):
        self.assertEqual(list(self.output.iterdir()), [self.output / ".gitignore"])
        self.assertEqual((self.output / ".gitignore").read_text(), "*\n!.gitignore\n")

    def test_restores_checked_files_with_and_without_github_prefix(self):
        for prefix in ("", "zipper-algebra-abcdef/"):
            with self.subTest(prefix=prefix):
                self.make_archive(prefix=prefix)
                output = self.root / ("github" if prefix else "plain")
                paths = restore_data(output, manifest=self.manifest, archive=self.archive)
                self.assertEqual(len(paths), 2)
                self.assertEqual((output / ".gitignore").read_text(), "*\n!.gitignore\n")
                for name, data in self.contents.items():
                    self.assertEqual((output / name).read_bytes(), data)

    def test_group_selection_does_not_require_other_groups_or_restore_extra_files(self):
        first = next(iter(self.contents))
        self.make_archive({first: self.contents[first], "extra.txt": b"ignored"})
        self.assertEqual(len(self.restore(group="inputs")), 1)
        self.assertEqual((self.output / first).read_bytes(), self.contents[first])
        self.assertFalse((self.output / "extra.txt").exists())
        self.assertFalse((self.output / "logs_and_traces").exists())

    def test_missing_selected_file_prevents_any_installation(self):
        first = next(iter(self.contents))
        self.make_archive({first: self.contents[first]})
        with self.assertRaisesRegex(RestoreError, "missing selected files"):
            self.restore()
        self.assert_only_ignore_marker()

    def test_bad_hash_prevents_any_installation(self):
        changed = dict(self.contents)
        changed["logs_and_traces/old.txt"] = b"X" * len(changed["logs_and_traces/old.txt"])
        self.make_archive(changed)
        with self.assertRaisesRegex(RestoreError, "SHA256 mismatch"):
            self.restore()
        self.assert_only_ignore_marker()

    def test_ignore_marker_exists_before_archive_verification(self):
        self.make_archive()

        def inspect(bundle, entries, staging):
            self.assertEqual((self.output / ".gitignore").read_text(), "*\n!.gitignore\n")
            self.assertFalse((self.output / "ms_inputs/example.ms").exists())
            return _verify_archive(bundle, entries, staging)

        with patch("experiments.restore_data._verify_archive", side_effect=inspect) as verify:
            self.restore()
        verify.assert_called_once()

    def test_bad_size_is_rejected(self):
        self.make_archive({name: data + b"!" for name, data in self.contents.items()})
        with self.assertRaisesRegex(RestoreError, "Size mismatch"):
            self.restore()

    def test_unsafe_archive_member_paths_are_rejected(self):
        for name in ("../outside", "/outside", "root/../outside", "root\\outside", "C:/outside"):
            with self.subTest(name=name):
                self.make_archive(extra=[tarfile.TarInfo(name)])
                with self.assertRaisesRegex(RestoreError, "Unsafe archive path"):
                    self.restore()
                self.assert_only_ignore_marker()

    def test_links_and_special_members_are_rejected(self):
        for kind in (tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.FIFOTYPE):
            with self.subTest(kind=kind):
                member = tarfile.TarInfo("link")
                member.type = kind
                member.linkname = "ms_inputs/example.ms"
                self.make_archive(extra=[member])
                with self.assertRaisesRegex(RestoreError, "not a regular file or directory"):
                    self.restore()
                self.assert_only_ignore_marker()

    def test_selected_member_must_be_a_regular_file(self):
        directory = tarfile.TarInfo("ms_inputs/example.ms/")
        directory.type = tarfile.DIRTYPE
        self.make_archive({}, extra=[directory])
        with self.assertRaisesRegex(RestoreError, "Expected a regular file"):
            self.restore(group="inputs")

    def test_duplicate_matching_archive_members_are_rejected(self):
        first = next(iter(self.contents))
        duplicate = tarfile.TarInfo("github-prefix/" + first)
        duplicate.size = 0
        self.make_archive(extra=[duplicate])
        with self.assertRaisesRegex(RestoreError, "Duplicate archive member"):
            self.restore()

    def test_unsafe_and_duplicate_manifest_paths_are_rejected(self):
        for name in ("../outside", "/outside", "root\\outside", "root\nsecond-query",
                     "root\rsecond-query", "root\tname", "root\x7fname"):
            with self.subTest(name=name):
                self.document["files"][0]["path"] = name
                self.write_manifest()
                with self.assertRaisesRegex(RestoreError, "Unsafe archive path"):
                    self.restore()
        self.document["files"] = [self.document["files"][1]] * 2
        self.write_manifest()
        with self.assertRaisesRegex(RestoreError, "Duplicate manifest path"):
            self.restore()

    def test_existing_different_content_is_preserved(self):
        target = self.output / "logs_and_traces/old.txt"
        target.parent.mkdir(parents=True)
        target.write_bytes(b"local work")
        self.make_archive()
        with self.assertRaisesRegex(RestoreError, "Existing destination differs; preserved"):
            self.restore()
        self.assertEqual(target.read_bytes(), b"local work")
        self.assertFalse((self.output / "ms_inputs/example.ms").exists())
        self.assertFalse((self.output / ".gitignore").exists())

    def test_existing_identical_files_are_reused_without_reading_an_archive(self):
        for name, data in self.contents.items():
            target = self.output / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        self.assertFalse(self.archive.exists())
        self.assertEqual(len(self.restore()), 2)
        self.assertEqual((self.output / ".gitignore").read_text(), "*\n!.gitignore\n")

    def test_existing_custom_ignore_marker_is_preserved(self):
        self.output.mkdir()
        marker = self.output / ".gitignore"
        marker.write_text("# local rules\n*.custom\n", encoding="utf-8")
        self.make_archive()
        self.restore()
        self.assertEqual(marker.read_text(), "# local rules\n*.custom\n")
        self.restore()
        self.assertEqual(marker.read_text(), "# local rules\n*.custom\n")

    def test_existing_symlink_cannot_redirect_output(self):
        outside = self.root / "outside"
        outside.mkdir()
        self.output.mkdir()
        (self.output / "ms_inputs").symlink_to(outside, target_is_directory=True)
        self.make_archive()
        with self.assertRaisesRegex(RestoreError, "Destination contains a symbolic link"):
            self.restore()
        self.assertEqual(list(outside.iterdir()), [])

    def test_uses_raw_pinned_git_blobs_despite_autocrlf_and_attributes(self):
        repository = self.root / "repository"
        repository.mkdir()
        for name, data in self.contents.items():
            target = repository / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        def git(*args):
            return subprocess.run(["git", "-C", str(repository), *args], check=True,
                                  capture_output=True, text=True).stdout.strip()
        git("init", "--quiet")
        git("config", "core.autocrlf", "true")
        (repository / ".gitattributes").write_text("* text eol=crlf\n", encoding="utf-8")
        git("add", ".")
        git("-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid",
            "commit", "--quiet", "-m", "Fixture")
        self.document["source_revision"] = git("rev-parse", "HEAD")
        self.write_manifest()
        for name in self.contents:
            (repository / name).unlink()
        git("checkout-index", "--force", "--all")
        for name, data in self.contents.items():
            self.assertEqual((repository / name).read_bytes(), data.replace(b"\n", b"\r\n"))
        # Changes in the checkout must not change the archived revision.
        (repository / "ms_inputs/example.ms").unlink()
        with patch("experiments.restore_data.urlopen", side_effect=AssertionError("network")):
            paths = restore_data(self.output, manifest=self.manifest, source_root=repository)
        self.assertEqual(len(paths), 2)
        self.assertEqual((self.output / "ms_inputs/example.ms").read_bytes(),
                         self.contents["ms_inputs/example.ms"])

    def test_download_fallback_uses_the_pinned_url_and_checks_contents(self):
        self.make_archive(prefix="github-prefix/")
        with patch("experiments.restore_data.subprocess.run", side_effect=FileNotFoundError), \
                patch("experiments.restore_data.urlopen",
                      return_value=io.BytesIO(self.archive.read_bytes())) as download:
            paths = restore_data(self.output, manifest=self.manifest, source_root=self.root)
        self.assertEqual(len(paths), 2)
        self.assertEqual(download.call_args.args[0].full_url,
                         self.document["source_archive_url"])


if __name__ == "__main__":
    unittest.main()
