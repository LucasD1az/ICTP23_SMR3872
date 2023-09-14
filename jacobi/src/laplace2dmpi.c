#include <math.h>
#include <string.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>


/*** function declarations ***/

// save matrix to file
void save_gnuplot( double *M, size_t dim, int me, int ncpu, int nloc, int offset );

// evolve Jacobi
void evolve( double * A, double *Anew, size_t n, int me, int prev, int next, int ncpu, int nloc, int offset, double w, double c, double dt, size_t iter, size_t nc, double *s, double *r, double *l, double *snew, double *rnew, double *lnew );

// return the elapsed time
double seconds( void );

/*** end function declaration ***/


int main(int argc, char* argv[]) {

    // timing variables
    double t_start, t_end;

    //indexes for loops
    size_t i, j, iter;

    // initialize the matrices
    double  *A, *Anew, *Atmp,
            *s, *snew, *stmp,
            *r, *rnew, *rtmp,
            *l, *lnew, *ltmp,
            *mean;

    size_t n = 0, iter_max = 0, nc=0, mdim=0;

    int me, ncpu, prev, next, offset, nloc;

    /* set up MPI environment */
    i = MPI_Init(&argc,&argv);
    if (i != MPI_SUCCESS) {
        puts("problem initializing MPI");
        return 1;
    }

    MPI_Comm_size(MPI_COMM_WORLD,&ncpu);
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    
    prev = me - 1;
    next = me + 1;


    if ( argc != 3 ) {
        fprintf(stderr, "Usage: %s <n> <w>\n", argv[0]);
        exit(-1);
    }

    n  = atoi(argv[1]);
    double w = atof(argv[2]);

    // set MPI variables
    nloc = n/ncpu;
    int rest=n%ncpu;
    if (me==0||me==ncpu-1){
        nloc++;
        if(ncpu==1) nloc--;
    } else {
        nloc+=2;
    }
    if (me<rest){
        nloc++;
        offset=0;
    } else offset=rest;

    mdim=n*nloc;
    //int n = NN;
    //int m = NM;
    //int nc = NM/2;
    nc=n/2;
    //int mc = NM/2;
    double dx = 0.5 / (n + 1);
//    double w = 6000;

    double dt = 0.000000025;//dx / v * 0.5;
    iter_max=10/w/dt;    //itermax*dt*w=#ciclos = 10
    int d_iter=iter_max/50;
    double v =1700.0;//c*dx/dt;
    w *= 2*M_PI;
    double c = v * dt / dx;
    printf("c=%lf\n",c);

    //int iter = 0;
    //int i,j;
    //struct timeval temp_1, temp_2;
    //double ttime=0.;

    // allocate memory for the matrices
    A    = (double *) malloc( mdim * sizeof(double) );
    Anew = (double *) malloc( mdim * sizeof(double) );
    s    = (double *) malloc( mdim * sizeof(double) );
    snew = (double *) malloc( mdim * sizeof(double) );
    r    = (double *) malloc( mdim * sizeof(double) );
    rnew = (double *) malloc( mdim * sizeof(double) );
    l    = (double *) malloc( mdim * sizeof(double) );
    lnew = (double *) malloc( mdim * sizeof(double) );
    mean = (double *) malloc( mdim * sizeof(double) );

    memset(A, 0, mdim * sizeof(double));
    memset(Anew, 0, mdim * sizeof(double));
    memset(r, 0, mdim * sizeof(double));
    memset(rnew, 0, mdim * sizeof(double));
    memset(s, 0, mdim * sizeof(double));
    memset(snew, 0, mdim * sizeof(double));
    memset(l, 0, mdim * sizeof(double));
    memset(lnew, 0, mdim * sizeof(double));
    memset(mean, 0, mdim * sizeof(double));
    for( j = 1; j < nloc-1; j++) {
            for( i = 1; i < n-1; i++ ) {
                int i_g=i;
                int j_g=j+me*nloc+offset;
                if(j_g>=nc-2&&j<nc+2&&i_g>=nc-2&&i<nc+2)
                    s[j*n+i]=1000*w;
            }
    }
    printf("Chladni Plate Resonance Calculation: %ld x %ld mesh\n", n, n);

    //gettimeofday(&temp_1, (struct timezone*)0);

    t_start = seconds();

    for (iter=0; iter < iter_max; iter++ )
    {
        evolve( A, Anew, n, me, prev, next, ncpu, nloc, offset, w, c, dt, iter, nc, s, r, l, snew, rnew, lnew );
        
        //printf("I am %d and I've just evolved!\n",me);
        if(iter>0.8*iter_max){
            for(i=0;i<nloc;i++){
                for(j=0;j<n;j++){
                    mean[(i*n)+j]+=(A[(i*n)+j]*A[(i*n)+j]);
                }
            }
        }

        // swap pointers
        Atmp = A; A = Anew; Anew = Atmp;
        stmp = s; s = snew; snew = stmp;
        rtmp = r; r = rnew; rnew = rtmp;
        ltmp = l; l = lnew; lnew = ltmp;
    }
    printf("I am %d and I've just finished!\n",me);
    for(i=0;i<nloc;i++){
        for(j=0;j<n;j++){
            mean[(i*n)+j]=sqrt(mean[(i*n)+j]/(0.2*iter_max));
        }
    }
    t_end = seconds();
    printf("Elapsed time (s) = %lf\n", t_end - t_start);

    save_gnuplot( mean, n, me, ncpu, nloc, offset );
    printf("Stopped at iteration: %ld\n", iter);
    // // print the final grid to file
    // char filename[64];
    // snprintf(filename, 64, "modonormal.dat");
    // FILE* fp = fopen(filename, "w");
    // for (j = 0; j < n; j++)
    // {
    //     for (i = 0; i < m; i++) fprintf(fp, "%8.4lf ", mean[j][i]);
    //     fprintf(fp, "\n");
    // }
    // gettimeofday(&temp_2, (struct timezone*)0);
    // ttime = 0.000001*((temp_2.tv_sec-temp_1.tv_sec)*1.e6+(temp_2.tv_usec-temp_1 .tv_usec));

    // printf("Elapsed time (s) = %.2lf\n", ttime);
    // printf("Stopped at iteration: %u\n", iter);
    MPI_Finalize();
    free(A);
    free(Anew);
    free(s);
    free(snew);
    free(r);
    free(rnew);
    free(l);
    free(lnew);
    free(mean);

    return 0;
}

void evolve( double * A, double *Anew, size_t n, int me, int prev, int next, int ncpu, int nloc, int offset, double w, double c, double dt, size_t iter, size_t nc, double *s, double *r, double *l, double *snew, double *rnew, double *lnew ){
    size_t i,j;
    MPI_Request reqAp, reqsp, reqrp, reqlp, reqAn, reqsn, reqrn, reqln;
    if(prev>=0){
        MPI_Irecv( A , n , MPI_DOUBLE , prev , 0 , MPI_COMM_WORLD , &reqAp);
        MPI_Irecv( s , n , MPI_DOUBLE , prev , 0 , MPI_COMM_WORLD , &reqsp);
        MPI_Irecv( r , n , MPI_DOUBLE , prev , 0 , MPI_COMM_WORLD , &reqrp);
        MPI_Irecv( l , n , MPI_DOUBLE , prev , 0 , MPI_COMM_WORLD , &reqlp);
        MPI_Send( A+n , n , MPI_DOUBLE , prev , 0 , MPI_COMM_WORLD );
        MPI_Send( s+n , n , MPI_DOUBLE , prev , 0 , MPI_COMM_WORLD );
        MPI_Send( r+n , n , MPI_DOUBLE , prev , 0 , MPI_COMM_WORLD );
        MPI_Send( l+n , n , MPI_DOUBLE , prev , 0 , MPI_COMM_WORLD );
        //printf("I am %d communicating with %d\n",me,prev);
    }
    if(next<ncpu){
        MPI_Irecv( A + (nloc-1)*n, n , MPI_DOUBLE , next , 0 , MPI_COMM_WORLD , &reqAn);
        MPI_Irecv( s+ (nloc-1)*n , n , MPI_DOUBLE , next , 0 , MPI_COMM_WORLD , &reqsn);
        MPI_Irecv( r+ (nloc-1)*n , n , MPI_DOUBLE , next , 0 , MPI_COMM_WORLD , &reqrn);
        MPI_Irecv( l+ (nloc-1)*n , n , MPI_DOUBLE , next , 0 , MPI_COMM_WORLD , &reqln);
        MPI_Send( A+ (nloc-2)*n , n , MPI_DOUBLE , next , 0 , MPI_COMM_WORLD );
        MPI_Send( s+ (nloc-2)*n , n , MPI_DOUBLE , next , 0 , MPI_COMM_WORLD );
        MPI_Send( r+ (nloc-2)*n , n , MPI_DOUBLE , next , 0 , MPI_COMM_WORLD );
        MPI_Send( l+ (nloc-2)*n , n , MPI_DOUBLE , next , 0 , MPI_COMM_WORLD );
        //printf("I am %d communicating with %d\n",me,next);
    }
    if(prev>=0){
        MPI_Wait(&reqAp,MPI_STATUS_IGNORE);
        MPI_Wait(&reqsp,MPI_STATUS_IGNORE);
        MPI_Wait(&reqrp,MPI_STATUS_IGNORE);
        MPI_Wait(&reqlp,MPI_STATUS_IGNORE);
    }
    if(next<ncpu){
        MPI_Wait(&reqAn,MPI_STATUS_IGNORE);
        MPI_Wait(&reqsn,MPI_STATUS_IGNORE);
        MPI_Wait(&reqrn,MPI_STATUS_IGNORE);
        MPI_Wait(&reqln,MPI_STATUS_IGNORE);
    }
    int i_g, j_g;
    for( j = 1; j < nloc-1; j++) {
            for( i = 1; i < n-1; i++ ) {
                i_g=i;
                j_g=j+me*nloc+offset;
                rnew[(j*n)+i]=r[(j*n)+i]+c*(0.5*(s[(j+1)*n+i]-s[(j-1)*n+i])+0.5*c*(r[(j+1)*n+i]-2*r[(j*n)+i]+r[(j-1)*n+i]));
                lnew[(j*n)+i]=l[(j*n)+i]+c*(0.5*(s[(j*n)+i+1]-s[(j*n)+i-1])+0.5*c*(l[(j*n)+i+1]-2*l[(j*n)+i]+l[(j*n)+i-1])); 
                if(j_g>=nc-2&&j_g<nc+2&&i_g>=nc-2&&i_g<nc+2){      
                    Anew[(j*n)+i]=1000*sin(dt*iter*w);
                    snew[(j*n)+i]=1000*w*cos(dt*iter*w);
                } else {
                    snew[(j*n)+i]=s[(j*n)+i]+c*(0.5*(r[(j+1)*n+i]-r[(j-1)*n+i]+l[(j*n)+i+1]-l[(j*n)+i-1])+0.5*c*(s[(j+1)*n+i]-4*s[(j*n)+i]+s[(j-1)*n+i]+s[(j*n)+i+1]+s[(j*n)+i-1]));
                    Anew[(j*n)+i]=A[(j*n)+i]+dt*0.5*(snew[(j*n)+i]+s[(j*n)+i]);
                }
                
            }
        }
    if(me==0)
        for(i=1;i<n-1;i++){
            Anew[i]=Anew[n+i];
            snew[i]=snew[n+i];
            rnew[i]=0;
            lnew[i]=0;
        }
    if(me==ncpu-1)
        for(i=1;i<n-1;i++){
            Anew[(nloc-1)*n+i]=Anew[(nloc-2)*n+i];
            snew[(nloc-1)*n+i]=snew[(nloc-2)*n+i];
            rnew[(nloc-1)*n+i]=0;
            lnew[(nloc-1)*n+i]=0;
        }
        for(j=1;j<nloc-1;j++){
            /*for(i=0;i>=0;i--){
                Anew[j][i]=A[j][i+1];
                snew[j][i]=s[j][i+1];
                rnew[j][i]=0;
                lnew[j][i]=0;
            }
            for(i=m-1;i<m;i++){
                Anew[j][i]=A[j][i-1];
                snew[j][i]=s[j][i-1];
                rnew[j][i]=0;
                lnew[j][i]=0;
            }*/
            Anew[j*n]=Anew[j*n+1];
            Anew[j*n+n-1]=Anew[j*n+n-2];
            snew[j*n]=snew[j*n+1];
            snew[j*n+n-1]=snew[j*n+n-2];
            rnew[j*n]=0;
            rnew[j*n+n-1]=0;
            lnew[j*n]=0;
            lnew[j*n+n-1]=0;
        }
}

void save_gnuplot( double *M, size_t dimension, int me, int ncpu, int nloc, int offset ){
  size_t i , j;
  const double h = 0.1;
  FILE *file;
  int next_n;
  double *Mrecv;
  Mrecv = (double *) malloc( dimension * sizeof(double) );


if(!me){
  file = fopen( "solution.dat", "w" );
    for( i = 0; i < nloc; ++i ){
    for( j = 0; j < dimension; ++j ) fprintf(file, "%f ", M[ ( i * ( dimension ) ) + j ] );
    fprintf(file, "\n" );
  }
    printf("I am %d and I've just printed myself!\n",me);
  for(i=1;i<ncpu;i++){
    MPI_Recv( &next_n , 1 , MPI_INT , i , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE);
    for(j=1;j<next_n-1;j++){
        MPI_Recv( Mrecv , dimension , MPI_DOUBLE , i , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE);
        for( j = 0; j < dimension; ++j ) fprintf(file, "%f ", Mrecv[ j ] );
        fprintf(file, "\n" );
    }
    printf("I am %d and I've just printed %ld!\n",me,i);
  }
  fclose( file );
} else {
    MPI_Send( &nloc , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD );
    for( i = 1; i < nloc-1; ++i ){
        MPI_Send( M + ( i * ( dimension ) ) , dimension , MPI_DOUBLE , 0 , 0 , MPI_COMM_WORLD );
    }
}
  

}

// A Simple timer for measuring the walltime
double seconds(){

    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;
}