#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#define N 14
//#define PRINT_WITH_MPI_GATHER

#include <time.h>
#include <sys/time.h>

double seconds()
/* Returns elepsed seconds past from the last call to timer rest */
{
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *) 0 );
  sec = tmp.tv_sec + ( (double) tmp.tv_usec ) / 1000000.0;
  return sec;
}

void print_mat (double * mat, int n_loc){
    for(int i=0; i<n_loc; i++){
        for(int j=0; j<N; j++) fprintf(stdout, "%.3g ", mat[i*N+j]);
    fprintf(stdout, "\n");
    }
}

int main(int argc, char * argv[]){
    double * mat;
    int size=1,rank=0,n_loc,rest,offset, i_g,j_g;
    double t1 = 0.0, t2 = 0.0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!rank) fprintf(stdout, "Running with %d processes...\n", size);
    n_loc = N/size;
    rest=N%size;
    if(rank<rest) {
        n_loc++;
        offset=0;
    }
    else offset=rest;
    //if(!rank) mat = (double *) malloc(N*N*sizeof(double));
    //else mat = (double *) malloc(n_loc*N*sizeof(double));
    //t1=seconds();
    mat = (double *) malloc(n_loc*N*sizeof(double));
    for(int i=0; i<n_loc; i++){
        for(int j=0; j<N; j++){
            j_g=j;
            i_g=i+rank*n_loc+offset;
            if(i_g==j_g){
                mat[i*N+j] = 1.0;
            }
            else{
                mat[i*N+j] = 0.0;
            }
        //fprintf(stdout, "%.3g ", mat[i*N+j]);
        }
    //fprintf(stdout, "\n");
    }
#ifdef PRINT_WITH_MPI_GATHER
    double * printmat;
    if(!rank) printmat = (double *) malloc(N*N*sizeof(double));
    MPI_Gather(mat, n_loc*N, MPI_DOUBLE, printmat, n_loc*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //t2=seconds();
    if(!rank){
        fprintf(stdout, "Printing with MPI Gather\n");
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++) fprintf(stdout, "%.3g ", printmat[i*N+j]);
        fprintf(stdout, "\n");
        }
        free(printmat);
    }
    //fprintf(stdout, "Time to solution %.3g (sec.)\n", t2-t1);
#else
    if(!rank){
        fprintf(stdout, "Printing with send-receive\n");
        print_mat(mat, n_loc);
        for(int count=1;count<size;count++){
            MPI_Recv(mat, n_loc*N, MPI_DOUBLE, count, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //El status ignore es si no te importa
            if(count>=rest&&rest) print_mat(mat, n_loc-1);
            else print_mat(mat, n_loc);
            //print_mat(mat, n_loc);
        }
    }
    else{
        MPI_Send(mat, n_loc*N, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);
    }
#endif
    free(mat);
    MPI_Finalize();
    return 0;
}