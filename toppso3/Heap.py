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
        self.indices = range(self.size)
        self.BuildHeap()


    def BuildHeap(self):
        """
        BuildHeap heapifies from the second lowest level of the heap
        """
        for i in reversed(range(self.nlevels - 1)):
            for j in range(int(2**i)):
                self.Heapify(int(2**i - 1 