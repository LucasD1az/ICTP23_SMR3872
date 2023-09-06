#!/usr/bin/env python

from mpi4py import MPI

comm = MPI.COMM_WORLD
me = comm.Get_rank()
ncpu = comm.Get_size()

# set rank of previous and next process
prev = me - 1
if prev < 0:
    prev = ncpu - 1

next = me + 1
if next == ncpu:
    next = 0

# at the beginning rank 0 will hold the message
mesgc = 0
mesgcc=0
if me == 0:
    mesgc = 1
    mesgcc=-1

for i in range(0,20):
    # post non-blocking receives first, then send, then wait for request to complete
    reqc = comm.irecv()
    reqcc = comm.irecv()
    s=comm.isend(mesgcc, prev)
    comm.send(mesgc, next)
    s.wait()
    #send.wait()
    mesg1 = reqc.wait()
    mesg2 = reqcc.wait()
    mesgcc=min(mesg1,mesg2)
    mesgc=max(mesg1,mesg2)
    # communication complete.
    if mesgcc:
        print("step %d: rank %d has counter-clockwise message" % (i, me))
    if mesgc:
        print("step %d: rank %d has clockwise message" % (i, me))
