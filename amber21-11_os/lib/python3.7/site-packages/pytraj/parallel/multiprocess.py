# do not use relative import here. Treat this module as a seperated package.
from __future__ import absolute_import
from functools import partial
from pytraj.core.c_options import info as compiled_info
from pytraj import matrix
from pytraj import mean_structure, volmap
from pytraj import ired_vector_and_matrix, rotation_matrix
from pytraj import NH_order_parameters
from multiprocessing import cpu_count
from pytraj.utils.tools import concat_dict
from pytraj.utils.get_common_objects import get_reference

from .base import worker_by_func
from .dataset import PmapDataset


def _pmap(func, traj, *args, **kwargs):
    '''use python's multiprocessing to accelerate calculation. Limited calculations.

    Parameters
    ----------
    func : a pytraj's methods or a list of string or simply as a cpptraj' text
    traj : pytraj.TrajectoryIterator
    n_cores : int, number of cores to be used, default 2. Specify n_cores=-1 to use all available cores
    iter_options : dict, default {}
        Specify trajectory iterating option. This will be done before calling ``func``.
    frame_indices : {None, array-like}, default None, optional
        if provided, pytraj will split this frame_indices into different chunks and let
        cpptraj perform calculation for specific indices.
        frame_indices must be pickable so is can be sent to different cores.

    *args, **kwargs: additional keywords

    Returns
    -------
    out : OrderedDict

    Notes
    -----
    - If you not sure about parallel's results, you should compare the output to serial run.

    - This is absolutely experimental. The syntax might be changed in future.

    Rule of thumbs: start with small number of frames (saying 10 frames), varying
    n_cores=1, 2, 3, 4 to see if the data makes sense or not.

    There are two modes in this method, use pytraj's methods (pytraj.rmsd, pytraj.radgyr,
    ...) or use cpptraj's command text syntax ('autoimage', 'rms', ...)

    If using cpptraj syntax::

        pytraj only supports limited cpptraj's Actions (not Analysis, checm Amber15 manual
        about Action and Analysis), say no  to 'matrix', 'atomicfluct', ... or any action
        that results output depending on the number of frames.


    This method only benifits you if your calculation is quite long (saying few minutes to
    few hours). For calculation that takes less than 1 minutes, you won't see the
    significant speed up (or even slower) since pytraj need to warm up and need to gather
    data when the calculation done.

    The parallel cacluation is very simple, trajectory will be split (almost equal) to
    different chunk (n_chunks = n_cores), pytraj/cpptraj perform calculation for each
    chunk in each core and then send data back to master. Note that we are using Python's
    built-in multiprocessing module, so you can use this method interactively in Ipython
    and ipython/jupyter notebook. This behavior is different from using MPI, in which you
    need to write a script, escaping ipython ession and type something like::

        mpirun -n 4 python my_script.py

    vs::

        In [1]: pt.pmap(pt.radgyr, traj, n_cores=4)
        Out[1]:
        OrderedDict([('RoG_00000',
                      array([ 18.91114428,  18.93654996,  18.84969884,  18.90449256,
                              18.8568644 ,  18.88917208,  18.9430491 ,  18.88878079,
                              18.91669565,  18.87069722]))])

    This is experimental method, you should expect its syntax, default output will be changed.

    When sending Topology to different cores, pytraj will reload Topology from
    traj.top.filename, so if you need to update Topology (in the fly), save it to disk and
    reload before using ``pytraj.pmap``

    Examples
    --------
    >>> import numpy as np
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')

    >>> # use iter_options
    >>> iter_options = {'autoimage': True, 'rmsfit': (0, '@CA')}
    >>> data = pt.pmap(pt.mean_structure, traj, iter_options=iter_options)

    >>> # cpptraj command style
    >>> data = pt.pmap(['distance :3 :7', 'vector mask :3 :12'], traj, n_cores=4)

    >>> # use reference. Need to explicitly use 'refindex', which is index of reflist
    >>> data = pt.pmap(['rms @CA refindex 0'], traj, ref=[traj[3],], n_cores=3)
    >>> data
    OrderedDict([('RMSD_00001', array([  2.68820312e-01,   3.11804885e-01,   2.58835452e-01,
             9.10475988e-08,   2.93310737e-01,   4.10197322e-01,
             3.96226694e-01,   3.66059215e-01,   3.90890362e-01,
             4.89180497e-01]))])

    >>> # use reference: if not want to use 'refindex', can use 'reference'
    >>> # the advantage is you can not specify a list of reference
    >>> data = pt.pmap(['rms @CA reference'], traj, ref=[traj[3],], n_cores=3)
    >>> data
    OrderedDict([('RMSD_00001', array([  2.68820312e-01,   3.11804885e-01,   2.58835452e-01,
             9.10475988e-08,   2.93310737e-01,   4.10197322e-01,
             3.96226694e-01,   3.66059215e-01,   3.90890362e-01,
             4.89180497e-01]))])

    >>> # use different references. Need to explicitly use 'refindex', which is index of reflist
    >>> # create a list of references
    >>> reflist = traj[3], traj[4]
    >>> # make sure to specify `refindex`
    >>> # `refindex 0` is equal to `reflist[0]`
    >>> # `refindex 1` is equal to `reflist[1]`
    >>> data = pt.pmap(['rms @CA refindex 0', 'rms !@H= refindex 1'], traj, ref=reflist, n_cores=2)
    >>> # convert to ndarray
    >>> data_arr = pt.tools.dict_to_ndarray(data)

    >>> # perform parallel calculation with given frame_indices
    >>> traj = pt.datafiles.load_tz2()
    >>> data = pt.pmap(pt.radgyr, traj, '@CA', frame_indices=range(10, 50), n_cores=4)
    >>> # serial version
    >>> data = pt.radgyr(traj, '@CA', frame_indices=range(10, 50))


    See also
    --------
    pytraj.pmap_mpi
    '''
    from multiprocessing import Pool
    from pytraj import TrajectoryIterator

    n_cores = kwargs.pop('n_cores') if 'n_cores' in kwargs else 2
    iter_options = kwargs.pop(
        'iter_options') if 'iter_options' in kwargs else {}
    apply = kwargs.pop('apply') if 'apply' in kwargs else None
    progress = kwargs.pop('progress') if 'progress' in kwargs else None
    progress_params = kwargs.pop(
        'progress_params') if 'progress_params' in kwargs else dict()

    if n_cores <= 0:
        # use all available cores
        n_cores = cpu_count()

    # update reference
    if 'ref' in kwargs:
        kwargs['ref'] = get_reference(traj, kwargs['ref'])

    if isinstance(func, (list, tuple, str)):
        # assume using _load_batch_pmap
        from pytraj.parallel.base import _load_batch_pmap
        #check_valid_command(func)
        data = _load_batch_pmap(
            n_cores=n_cores,
            traj=traj,
            lines=func,
            dtype='dict',
            root=0,
            mode='multiprocessing',
            **kwargs)
        data = concat_dict((x[0] for x in data))
        return data
    else:
        if not callable(func):
            raise ValueError('must callable argument')
        # pytraj's method
        if not hasattr(func, '_is_parallelizable'):
            raise ValueError("this method does not support parallel")
        elif not func._is_parallelizable:
            raise ValueError("this method does not support parallel")
        else:
            if hasattr(
                    func, '_openmp_capability'
            ) and func._openmp_capability and 'OPENMP' in compiled_info():
                raise RuntimeError(
                    "this method supports both openmp and pmap, but your cpptraj "
                    "version was installed with openmp. Should not use both openmp and pmap at the "
                    "same time. In this case, do not use pmap since openmp is more efficient"
                )

        if not isinstance(traj, TrajectoryIterator):
            raise ValueError('only support TrajectoryIterator')

        if 'dtype' not in kwargs and func not in [
                mean_structure,
                matrix.dist,
                matrix.idea,
                ired_vector_and_matrix,
                rotation_matrix,
                volmap,
        ]:
            kwargs['dtype'] = 'dict'

        # keyword
        if func is volmap:
            assert kwargs.get('size') is not None, 'must provide "size" value'

        p = Pool(n_cores)

        pfuncs = partial(
            worker_by_func,
            n_cores=n_cores,
            func=func,
            traj=traj,
            args=args,
            kwargs=kwargs,
            iter_options=iter_options,
            apply=apply,
            progress=progress,
            progress_params=progress_params)

        data = p.map(pfuncs, [rank for rank in range(n_cores)])
        p.close()

        dataset_processor = PmapDataset(
            data, func=func, kwargs=kwargs, traj=traj)
        return dataset_processor.process()


def pmap(func=None, traj=None, *args, **kwargs):
    if func != NH_order_parameters:
        return _pmap(func, traj, *args, **kwargs)
    else:
        if 'n_cores' in kwargs.keys():
            if kwargs['n_cores'] == 1:
                # use default n_cores=2 instead of 1
                kwargs['n_cores'] = 2
            if kwargs['n_cores'] <= 0:
                kwargs['n_cores'] = cpu_count()
        else:
            # use n_cores=2 for default value
            kwargs['n_cores'] = 2
        return NH_order_parameters(traj, *args, **kwargs)


pmap.__doc__ = _pmap.__doc__
