import numpy as np

class Heap():
    """Heap is a class of binary min. heap.
       Attributes:
           data    -- a list of data
           size    -- the number of data in the list
           nlevels -- level of the heap
           indices -- a list of indices of each member in the heap
    """
    
    def __init__(self, datalist):
        self.data = datalist
        self.size = len(datalist)
        self.nlevels = int(np.floor(np.log2(self.size))) + 1
        self.indices = ran