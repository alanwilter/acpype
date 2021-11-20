from __future__ import absolute_import
import os
import numbers


def eq(arr0, arr1):
    assert arr0 == arr1


def file_exist(filename):
    import os
    return os.path.isfile(filename)


def is_range(obj):
    return 'range' in obj.__class__.__name__


def is_array(obj):
    """check if a `word` is in obj.__class__.__name__
    """
    return 'array' in obj.__class__.__name__


def are_instance(obj_list, cls):
    """check if all elements have the same class `cls`"""
    for element in obj_list:
        if not isinstance(element, cls):
            return False
    return True


def is_generator(iter_obj):
    # use this method instead of `inspect` in python since this does not work with Cython
    # Reason: (I don't know)
    if iter_obj.__class__.__name__ == 'generator':
        return True
    else:
        return False


def is_frame_iter(iter_obj):
    """check if is frame_iter
    """
    if iter_obj.__class__.__name__ == 'FrameIterator':
        return True
    return False


def is_int(num):
    """wrapper class to check if `num` is int
    isinstance(nu, (int, long)) does not work with numpy.int64, so we use numbers.Integral
    """
    return isinstance(num, numbers.Integral)


def is_number(num):
    return isinstance(num, numbers.Number)


def ensure_exist(filename):
    '''
    >>> ensure_exist('xfdasfda33fe')
    Traceback (most recent call last):
        ...
    RuntimeError: can not find xfdasfda33fe
    '''
    if not os.path.exists(filename):
        txt = "can not find %s" % filename
        raise RuntimeError(txt)


def ensure_not_none_or_string(obj):
    '''
    >>> ensure_not_none_or_string(None)
    Traceback (most recent call last):
        ...
    ValueError: <None> is a wrong input. Can not use `None` or string type
    >>> ensure_not_none_or_string('test')
    Traceback (most recent call last):
        ...
    ValueError: <test> is a wrong input. Can not use `None` or string type
    '''
    name = obj.__str__()
    msg = "<%s> is a wrong input. Can not use `None` or string type" % name
    if obj is None or isinstance(obj, str):
        raise ValueError(msg)


def _import(modname):
    """has_numpy, np = _import('numpy')
    >>> has_np, np = _import('numpy')
    """
    has_module = False
    try:
        imported_mod = __import__(modname)
        has_module = True
        return (has_module, imported_mod)
    except ImportError:
        has_module = False
        return (has_module, None)


def has_(lib):
    """check if having `lib` library

    Examples
    --------
    >>> has_np = has_("numpy")
    """
    return _import(lib)[0]
