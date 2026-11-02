#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Test explicit temporary-artifact ownership on success and failure."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Tools.scratch import ScratchWorkspace


class ScratchTest(unittest.TestCase):
    def test_cleanup_preserves_unregistered_files(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            parent = Path(root)
            user_file = parent / 'input.vcf'
            user_file.touch()
            with ScratchWorkspace() as scratch:
                owned = scratch.register(parent / 'output.vcf')
                owned.touch()
                scratch.register(parent / 'missing.tbi')
                scratch.register(owned)
            self.assertFalse(owned.exists())
            self.assertTrue(user_file.exists())

    def test_keep_and_exception(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            for keep in (False, True):
                owned = Path(root) / str(keep)
                with self.assertRaisesRegex(ValueError, 'comparison failed'):
                    with ScratchWorkspace(keep=keep) as scratch:
                        scratch.register(owned).touch()
                        raise ValueError('comparison failed')
                self.assertEqual(owned.exists(), keep)

    def test_directory_and_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            parent = Path(root)
            target = parent / 'user-directory'
            target.mkdir()
            with ScratchWorkspace() as scratch:
                directory = scratch.register(parent / 'owned', directory=True)
                directory.mkdir()
                (directory / 'child').touch()
                link = parent / 'link'
                link.symlink_to(target, target_is_directory=True)
                scratch.register(link, directory=True)
            self.assertFalse(directory.exists())
            self.assertFalse(link.is_symlink())
            self.assertTrue(target.exists())


if __name__ == '__main__':
    unittest.main()
