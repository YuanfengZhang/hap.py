# -*- coding: utf-8 -*-
"""Explicit ownership and cleanup of temporary comparison artifacts."""
import logging
import os
from pathlib import Path
import shutil
from types import TracebackType
from typing import Self


class ScratchWorkspace:
    """Clean registered artifacts only; preserve all artifacts when requested.

    Register only paths created by this run. Parent directories and unregistered
    files are never inferred or removed. Cleanup failures do not mask the
    comparison exception, but are reported rather than silently discarded.
    """

    def __init__(self, keep: bool = False) -> None:
        self.keep = keep
        self._paths: dict[Path, bool] = {}

    def register(self, path: str | os.PathLike[str], *, directory: bool = False) -> Path:
        """Take ownership of one file or an explicitly identified directory."""
        owned = Path(path)
        self._paths[owned] = directory
        return owned

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type: type[BaseException] | None,
                 exc: BaseException | None, traceback: TracebackType | None) -> None:
        if self.keep:
            logging.info('Scratch files kept: %s', list(map(str, self._paths)))
            return
        for path, directory in reversed(self._paths.items()):
            try:
                if directory and not path.is_symlink():
                    shutil.rmtree(path)
                else:
                    path.unlink()
            except FileNotFoundError:
                pass
            except OSError as error:
                logging.warning('Could not remove scratch artifact %s: %s', path, error)
