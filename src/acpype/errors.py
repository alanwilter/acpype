"""Exceptions ACPYPE raises on purpose.

A deliberate refusal carries a message written for the person running the tool. The
command line prints that message alone, where an unexpected exception gets the full
traceback.
"""


class AcpypeError(Exception):
    """A failure ACPYPE detected and explained itself; no traceback is needed."""


class UnsupportedTopologyError(AcpypeError):
    """The input topology uses a feature ACPYPE cannot convert yet."""
