"""Local web UI for p2d.  A shell over the library: no design logic lives here."""

from .server import serve

__all__ = ["serve"]
