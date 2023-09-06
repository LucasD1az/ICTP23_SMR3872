PROGRAM hello
    INCLUDE 'mpif.h'
    INTEGER :: ierr, rank, size
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    PRINT* , 'I am ', rank, ' of ', size
    CALL MPI_FINALIZE(ierr)
END