"""this file has commonly used actions such as rmsd calculation,
randomizeions, strip atoms, ..."""

from __future__ import print_function, absolute_import
import os
from glob import glob
from pytraj.core.c_options import set_world_silent

from .context import capture_stdout

try:
    from pytraj.externals.magic import from_file as file_type_info
except ImportError:
    file_type_info = None


def parallel_info(key=None):
    '''Parallel info

    Returns
    -------
    out: {str, dict}
        if key is None, return dict of all methods
        if key is not None, return specific method. key={'openmp', 'pmap'}
    '''
    import pytraj as pt
    from pytraj.analysis import matrix, vector, nmr
    from pytraj import cluster
    from itertools import chain
    methodlist_pmap = []
    methodlist_openmp = []

    for method_str in chain(
            dir(pt), dir(matrix), dir(vector), dir(nmr), dir(cluster)):
        try:
            method = getattr(pt, method_str)
            if hasattr(method,
                       '_is_parallelizable') and method._is_parallelizable:
                methodlist_pmap.append(method)
            if hasattr(method,
                       '_openmp_capability') and method._openmp_capability:
                methodlist_openmp.append(method)
        except AttributeError:
            pass

    pmap_ = []
    for method in set(methodlist_pmap):
        name = str(method).split()[1]
        # if 'calc_' in name:
        #    name = name.split('calc_')[-1]
        pmap_.append(name)
    supported_pmap_methods = sorted(pmap_)

    openmp_ = []
    for method in set(methodlist_openmp):
        name = str(method).split()[1]
        # if 'calc_' in name:
        #    name = name.split('calc_')[-1]
        openmp_.append(name)
    supported_openmp_methods = sorted(openmp_)
    d = {'pmap': supported_pmap_methods, 'openmp': supported_openmp_methods}
    return d if key is None else d[key]


def info(obj=None):  # pragma: no cover
    """get `help` for obj
    Useful for Actions and Analyses

    Since we use `set_worl_silent` to turn-off cpptraj' stdout, we need
    to turn on to use cpptraj's help methods
    """
    from pytraj import adict, analdict
    adict_keys = adict.keys()
    anal_keys = analdict.keys()

    if obj is None:
        print("action's keys", adict_keys)
        print("analysis' keys", anal_keys)
    else:
        if isinstance(obj, str):
            if obj == 'parallel':
                print(parallel_info())
                return
            elif obj in adict.keys():
                # make Action object
                _obj = adict[obj]
            elif obj in analdict.keys():
                # make Analysis object
                _obj = analdict[obj]
            else:
                raise ValueError("keyword must be an Action or Analysis")
        else:
            # assume `obj` hasattr `help`
            _obj = obj

        with capture_stdout() as (out, err):
            if hasattr(_obj, 'help'):
                set_world_silent(False)
                _obj.help()
                set_world_silent(True)
            elif hasattr(_obj, 'info'):
                set_world_silent(False)
                _obj.info()
                set_world_silent(True)
            elif 'calc_' in _obj.__name__:
                key = _obj.__name__.split("_")[-1]
                set_world_silent(False)
                adict[key].help()
                set_world_silent(True)
            elif hasattr(_obj, '__doc__'):
                print(_obj.__doc__)
            else:
                raise ValueError("object does not have `help` method")
        print(out.read())


def show_code(func, get_txt=False):  # pragma: no cover
    """show code of func or module"""
    import inspect
    txt = inspect.getsource(func)
    if not get_txt:
        print(txt)
    else:
        return txt


def get_atts(obj):  # pragma: no cover
    """get methods and atts from obj but excluding special methods __"""
    atts_dict = dir(obj)
    return [a for a in atts_dict if not a.startswith("__")]


def find_libcpptraj(**kwd):  # pragma: no cover
    '''
    '''
    return find_library('cpptraj', **kwd)


def find_library(libname, unique=False):  # pragma: no cover
    """return a list of all library files"""
    paths = os.environ.get('LD_LIBRARY_PATH', '').split(':')
    lib_path_list = []
    key = "lib" + libname + "*"

    for path in paths:
        path = path.strip()
        fnamelist = glob(os.path.join(path, key))
        for fname in fnamelist:
            if os.path.isfile(fname):
                lib_path_list.append(fname)

    if not lib_path_list:
        return None
    else:
        if unique:
            return set(lib_path_list)
        else:
            return lib_path_list
