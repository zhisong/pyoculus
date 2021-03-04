try:
    from importlib import metadata
except ImportError:
    # Running on pre-3.8 Python; use importlib-metadata package
    import importlib_metadata as metadata

__version__ = metadata.version('pyoculus')

from . import problems
from . import integrators
from . import irrationals
from . import solvers
