"""Re-creation of the TED (The Encyclopedia of Domains) domain chopping and CATH labelling."""
import ctypes as _ctypes
import sys as _sys
from pathlib import Path as _Path

# The .venv python is the conda env "ted". Its pyarrow (pulled in by transformers -> sklearn when ESMFold
# is imported) needs a newer libstdc++ than /lib64 has. Whichever libstdc++ loads first wins, so load the
# env's copy globally before anything else imports a C++ extension.
_libstdcxx = _Path(_sys.executable).resolve().parent.parent / "lib" / "libstdc++.so.6"
if _libstdcxx.exists():
    try:
        _ctypes.CDLL(str(_libstdcxx), mode=_ctypes.RTLD_GLOBAL)
    except OSError:
        pass

from .chop import ted_chop, ted_chop_batch, ChopResult, Domain  # noqa: E402,F401
