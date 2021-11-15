from __future__ import absolute_import
from abc import ABCMeta, abstractmethod
from pytraj.externals.six import add_metaclass
from pytraj.shared_methods import my_str_method


@add_metaclass(ABCMeta)
class TrajectoryBaseIterator(object):
    @abstractmethod
    def __init__(self, filename=None, top=None):
        pass

    def __str__(self):
        return my_str_method(self)

    def __repr__(self):
        return self.__str__()

    @abstractmethod
    def __iter__(self):
        pass

    @abstractmethod
    def __getitem__(self, idx):
        pass

    @abstractmethod
    def __enter__(self):
        pass

    @abstractmethod
    def __exit__(self, *args):
        pass

    @property
    @abstractmethod
    def n_frames(self):
        pass

    @property
    @abstractmethod
    def n_atoms(self):
        pass

    @property
    @abstractmethod
    def xyz(self):
        pass

    @property
    @abstractmethod
    def filename(self):
        pass

    @property
    @abstractmethod
    def filelist(self):
        pass

    @property
    @abstractmethod
    def _frame_slice_list(self):
        pass

    @abstractmethod
    def iterframe(self, *args, **kwd):
        pass

    @abstractmethod
    def iterchunk(self, *args, **kwd):
        pass
