""""""
from __future__ import absolute_import
import inspect
from ...datasets.c_datasetlist import DatasetList as CpptrajDatasetList
from . import c_action
from ...utils.context import capture_stdout
from ...trajectory.shared_methods import iterframe_master
from io import StringIO


def do_action(traj, command, action_class, post_process=True, top=None):
    ''' For internal use

    Parameters
    ----------
    traj : Trajectory-like
    command : str
    action_class : derived class of c_action.Action
    '''
    assert inspect.isclass(
        action_class), 'must passing a derived class of c_action.Action'
    assert issubclass(action_class, c_action.Action)
    c_dslist = CpptrajDatasetList()
    top = traj.top if top is None else top
    act = action_class(command=command, top=top, dslist=c_dslist)
    with capture_stdout() as (out, _):
        for frame in iterframe_master(traj):
            act.compute(frame)
        if post_process:
            act.post_process()
    # make sure to use StringIO
    # If not, if not, there won't be output in the next function call
    return c_dslist, StringIO(out.read()).read()
