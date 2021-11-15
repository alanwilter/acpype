from __future__ import print_function, absolute_import
import os
from functools import wraps

# we duplicate code from .utils.check_and_assert here to avoid circular import


def register_pmap(f):
    @wraps(f)
    def inner(*args, **kwd):
        return f(*args, **kwd)

    inner._is_parallelizable = True
    return inner


def register_openmp(f):
    @wraps(f)
    def inner(*args, **kwd):
        return f(*args, **kwd)

    inner._openmp_capability = True
    return inner


def ensure_exist(f):
    @wraps(f)
    def inner(*args, **kwd):
        if 'filename' in kwd.keys() and not os.path.exists(kwd['filename']):
            raise RuntimeError('filename not exist')
        elif not os.path.exists(args[0]):
            raise RuntimeError('filename not exist')
        return f(*args, **kwd)

    return inner


def has_(lib):
    """check if having `lib` library

    Examples
    --------
    >>> has_("numpy")
    True
    """
    try:
        __import__(lib)
        return True
    except ImportError:
        return False


def makesureABC(classname):
    def inner(func):
        def _inner(self, *args, **kwd):
            if self.__class__.__name__ == classname:
                raise NotImplementedError("This is Abstract Base Class")
            else:
                return func(self, *args, **kwd)

        # update _inner doc for func
        _inner.__doc__ = func.__doc__
        _inner.__name__ = func.__name__
        return _inner

    return inner
