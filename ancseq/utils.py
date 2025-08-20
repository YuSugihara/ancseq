"""Utility helpers for ancseq (logging, formatting)."""

import os
import sys
from datetime import datetime

__all__ = ["time_stamp", "clean_cmd", "call_log"]


def time_stamp():
    """Return a standardized timestamp prefix used in logs."""
    return f"[ancseq:{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]"


def clean_cmd(cmd):
    """Normalize whitespace in a command string."""
    return " ".join(cmd.split())


def call_log(name, cmd):
    print(time_stamp(), f"!!ERROR!! {cmd}\n", flush=True)
    print("please check below:\n")
    with open(name) as log:
        for line in log:
            print(line, end="")
