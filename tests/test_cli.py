# -*- coding: utf-8 -*-

"""
test_cli
----------------------------------

Tests for cli module.
"""

import pytest

from vcftoolz import cli


def test_error_on_empty_command_line():
    """Verify exception on empty command line."""
    with pytest.raises(SystemExit):
        cli.run_from_line("")
