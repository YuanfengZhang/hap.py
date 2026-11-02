#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fresh build directories avoid stale caches without hashing the source tree."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location('happy_installer', ROOT / 'install.py')
installer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(installer)


class BuildDirectoryTest(unittest.TestCase):
    def test_fresh_directories_preserve_existing_contents(self) -> None:
        with tempfile.TemporaryDirectory(dir=ROOT / '.work') as folder:
            base = Path(folder)
            first = installer.create_build_directory(ROOT, base / 'install', base / 'scratch')
            marker = first / 'marker'
            marker.write_text('preserve')
            second = installer.create_build_directory(ROOT, base / 'install', base / 'scratch')
            self.assertNotEqual(first, second)
            self.assertEqual(marker.read_text(), 'preserve')
            self.assertEqual(list(second.iterdir()), [])

    def test_unsafe_paths_rejected(self) -> None:
        with tempfile.TemporaryDirectory(dir=ROOT / '.work') as folder:
            target = Path(folder) / 'install'
            for scratch in (ROOT, target, target / 'nested'):
                with self.subTest(scratch=scratch), self.assertRaises(ValueError):
                    installer.create_build_directory(ROOT, target, scratch)


if __name__ == '__main__':
    unittest.main()
