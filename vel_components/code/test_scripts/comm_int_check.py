from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    ids = np.arange(2787742702, 2787742722)
    pos = np.random.rand(20, 3) * 100.
    arr = np.zeros((20, 7))
    arr[:, 0] = ids
    arr[:, 1:4] = pos
    lst = np.split(arr, [10,])
    del arr

else:
    lst = None

arr = comm.scatter(lst, root=0)

print("Rank {}, arr:".format(rank), arr)
