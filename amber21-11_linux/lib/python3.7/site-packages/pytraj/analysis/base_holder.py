from pytraj.datasets.datasetlist import DatasetList


class BaseDataHolder(object):
    def __init__(self, dslist=None):
        '''
        >>> holder = BaseDataHolder()
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')
        >>> dslist = pt.radgyr(traj, dtype='dataset')
        >>> holder = BaseDataHolder(dslist)

        >>> holder.data
        <pytraj.DatasetList with 1 datasets>
        RoG_00000
        [ 18.91114428  18.93654996  18.84969884  18.90449256  18.8568644
          18.88917208  18.9430491   18.88878079  18.91669565  18.87069722]

        >>> holder.to_dict()
        OrderedDict([('RoG_00000', array([ 18.91114428,  18.93654996,  18.84969884,  18.90449256,
                18.8568644 ,  18.88917208,  18.9430491 ,  18.88878079,
                18.91669565,  18.87069722]))])
        >>> for x in holder: pass
        '''
        if dslist is not None:
            self._dslist = dslist
        else:
            self._dslist = DatasetList()

    @property
    def data(self):
        '''return pytraj.DatasetList
        '''
        return self._dslist

    def to_dict(self):
        '''return OrderedDict
        '''
        return self._dslist.to_dict()

    def __getitem__(self, idx):
        return self.__class__(self._dslist[idx])

    def __iter__(self):
        return self._dslist.__iter__()

    @property
    def values(self):
        '''return raw numpy array
        '''
        return self._dslist.values
