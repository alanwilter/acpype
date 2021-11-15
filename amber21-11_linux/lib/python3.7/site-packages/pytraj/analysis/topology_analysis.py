from functools import partial
from ..core.c_core import CpptrajState, Command
from ..utils.context import capture_stdout
from io import StringIO


def compare_topology(top0, top1):
    ''' top0, top1 are :class:`pytraj.Topology`

    Notes
    -----
    Experiment 
    '''

    # TODO : make copies of topologies?
    state = CpptrajState()
    t0 = state.data.add('topology', name='top0')
    t0.data = top0

    t1 = state.data.add('topology', name='top1')
    t1.data = top1

    with capture_stdout() as (out, _):
        with Command() as command:
            command.dispatch(state, 'comparetop top0 top1')
    content = StringIO(out.read())
    return content.read()


def _whatinfo(top, command, ref, task):
    ''' info for given Topology and Frame.
    This is not meant to be fast.

    Notes
    -----
    Experiment
    '''

    state = CpptrajState()
    top_set = state.data.add('topology')
    top_set.data = top

    ref_set = state.data.add('reference')
    ref_set.top = top
    ref_set.add_frame(ref)

    command_ = ' '.join((task, command, 'reference'))
    with capture_stdout() as (out, _):
        with Command() as cpp_command:
            cpp_command.dispatch(state, command_)
    content = StringIO(out.read())
    # note: need to flush the content by print
    # if not: calling more than two commands (e.g: 'bondinfo', 'angleinfo') will not
    # show anything.
    return content.read()
    # return content


atominfo = partial(_whatinfo, task='atominfo')
resinfo = partial(_whatinfo, task='resinfo')
bondinfo = partial(_whatinfo, task='bondinfo')
angleinfo = partial(_whatinfo, task='angleinfo')
dihedralinfo = partial(_whatinfo, task='dihedralinfo')
