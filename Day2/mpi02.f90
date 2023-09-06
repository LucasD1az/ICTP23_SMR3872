program bcast
implicit none
include 'mpif.h'

integer::myrank,ncpus,imesg,ierr
integer,parameter::comm=MPI_COMM_WORLD

call MPI_INIT(ierr)
call MPI_COMM_RANK(comm,myrank,ierr)
call MPI_COMM_SIZE(comm,ncpus,ierr)

imesg=myrank
print *, "Before Bcast operation I'm ",myrank, &
    "and my message content is ", imesg
call MPI_BCAST(imesg,1,MPI_INTEGER,0,comm,ierr)

print *, "After Bcast operation I'm ",myrank, &
    "and my message content is ", imesg

call MPI_FINALIZE(ierr)
end program bcast