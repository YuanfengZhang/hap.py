#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Preserve upstream notices and valid source-script interpreter directives."""
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[3]
NOTICE = ('# This is a modified version of the original software.\n'
          "# The original authors' notice follows.\n")


class SourceHeadersTest(unittest.TestCase):
    def test_modified_python_notices(self) -> None:
        for path in ROOT.rglob('*'):
            if any(part in ('.git', '.work', '.pixi', '__pycache__') for part in path.relative_to(ROOT).parts):
                continue
            if not path.is_file() or not (path.suffix == '.py' or path.name.endswith('.py.in')):
                continue
            text = path.read_text()
            for match in re.finditer(r'^# Copyright \(c\)[^\n]*Illumina', text, re.MULTILINE):
                with self.subTest(path=path):
                    self.assertTrue(text[:match.start()].endswith(NOTICE))

    def test_python_encoding_and_entrypoints(self) -> None:
        for path in ROOT.rglob('*'):
            if any(part in ('.git', '.work', '.pixi', '__pycache__') for part in path.relative_to(ROOT).parts):
                continue
            if not path.is_file() or not (path.suffix == '.py' or path.name.endswith('.py.in')):
                continue
            text = path.read_text()
            lines = text.splitlines()
            with self.subTest(path=path):
                index = 1 if lines[0].startswith('#!') else 0
                self.assertEqual(lines[index], '# -*- coding: utf-8 -*-')
                if re.search(r"if\s+__name__\s*==\s*['\"]__main__['\"]", text):
                    self.assertEqual(lines[0], '#!/usr/bin/env python3')

    def test_script_shebangs(self) -> None:
        for folder in ('src/python', 'src/sh', 'external'):
            for path in (ROOT / folder).rglob('*'):
                if not path.is_file() or path.suffix not in ('.py', '.sh'):
                    continue
                if any(part in ('.work', '.pixi', '__pycache__') for part in path.parts):
                    continue
                for index, line in enumerate(path.read_text().splitlines()):
                    if line.startswith('#!'):
                        with self.subTest(path=path, line=index + 1):
                            self.assertEqual(index, 0)
                            expected = '#!/usr/bin/env python3' if path.suffix == '.py' else None
                            if expected:
                                self.assertEqual(line, expected)
                            else:
                                self.assertIn(line, ('#!/bin/bash', '#!/bin/sh', '#!/usr/bin/env bash'))


if __name__ == '__main__':
    unittest.main()
