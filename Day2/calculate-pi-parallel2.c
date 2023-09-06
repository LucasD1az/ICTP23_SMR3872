#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#define N 10000000

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

int main(int argc, char * argv[]){
    int n = N;
    double w=1.0/n;
    double x,fx,sum,tot_sum=0.0;


    double t1 = 0.0, t2 = 0.0;

    int size=1,rank=0,n_loc,start_i,end_i;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if( !rank ) fprintf( stdout, "Running with %d processes...\n", size);
    
    // N might be not divisible by size
    n_loc=n/size;

    t1 = seconds();
    
    for(int i = 1+rank; i<=n; i+=size){
        x = (i-0.5)*w;
        fx = 4.0/(1.0+x*x);
        sum += fx;
    }
    MPI_Reduce(&sum, &tot_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    t2 = seconds();
  
    if( !rank ) fprintf( stdout, "The value of PI is %.10g Vs %.10g\nTime to solution %.3g (sec.)\n", tot_sum * w, M_PI, t2 - t1 );
    //if(!rank) fprintf(stdout, "The value of PI is %.10g vs %.10g\n", tot_sum*w, M_PI);
    MPI_Finalize();
    return 0;
}