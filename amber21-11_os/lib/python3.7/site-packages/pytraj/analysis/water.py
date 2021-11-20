from ..utils.get_common_objects import get_data_from_dtype
from ..utils.decorators import register_pmap, register_openmp
from .c_action import c_action
from .c_action import do_action


@register_pmap
@register_openmp
def spam(traj, peak_file, command, dtype='dict'):
    ''' Limited support for SPAM. This method will be changed.

    TODO: update doc

    Parameters
    ----------
    traj : Trajectory-like
    peak_file : str
    command : str
        cpptraj spam command
    dtype : str
        return data type
    '''

    command_ = ' '.join((peak_file, command))
    c_dslist, _ = do_action(traj, command_, c_action.Action_Spam)
    return get_data_from_dtype(c_dslist, dtype=dtype)
