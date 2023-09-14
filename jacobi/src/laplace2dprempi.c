#include <math.h>
#include <string.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define NN 128
#define NM 128

/*** function declarations ***/

// save matrix to file
void save_gnuplot( double *M, size_t dim );

// evolve Jacobi
void evolve( double * A, double *Anew, size_t n, double w, double c, double dt, size_t iter, size_t nc, double *s, double *r, double *l, double *snew, double *rnew, double *lnew );

// return the elapsed time
double seconds( void );

/*** end function declaration ***/

double A[NN][NM];
double Anew[NN][NM];
double mean[NN][NM];
double r[NN][NM];
double rnew[NN][NM];
double s[NN][NM];
double snew[NN][NM];
double l[NN][NM];
double lnew[NN][NM];

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

    size_t n = 0, iter_max = 0, nc=0;
    
    if ( argc != 3 ) {
        fprintf(stderr, "Usage: %s <n> <w>\n", argv[0]);
        exit(-1);
    }

    n  = atoi(argv[1]);
    double w = atof(argv[2]);

    //int n = NN;
    //int m = NM;
    //int nc = NM/2;
    nc=n/2;
    //int mc = NM/2;
    double dx = 0.5 / (n - 1);
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
    A    = (double *) malloc( n * n * sizeof(double) );
    Anew = (double *) malloc( n * n * sizeof(double) );
    s    = (double *) malloc( n * n * sizeof(double) );
    snew = (double *) malloc( n * n * sizeof(double) );
    r    = (double *) malloc( n * n * sizeof(double) );
    rnew = (double *) malloc( n * n * sizeof(double) );
    l    = (double *) malloc( n * n * sizeof(double) );
    lnew = (double *) malloc( n * n * sizeof(double) );
    mean = (double *) malloc( n * n * sizeof(double) );

    memset(A, 0, n * n * sizeof(double));
    memset(Anew, 0, n * n * sizeof(double));
    memset(r, 0, n * n * sizeof(double));
    memset(rnew, 0, n * n * sizeof(double));
    memset(s, 0, n * n * sizeof(double));
    memset(snew, 0, n * n * sizeof(double));
    memset(l, 0, n * n * sizeof(double));
    memset(lnew, 0, n * n * sizeof(double));
    memset(mean, 0, n * n * sizeof(double));
    for(j=nc-2;j<nc+2;j++)
        for(i=nc-2;i<nc+2;i++)
            s[j*n+i]=1000*w;
    printf("Chladni Plate Resonance Calculation: %ld x %ld mesh\n", n, n);

    //gettimeofday(&temp_1, (struct timezone*)0);

    t_start = seconds();

    for (iter=0; iter < iter_max; iter++ )
    {
        evolve( A, Anew, n, w, c, dt, iter, nc, s, r, l, snew, rnew, lnew );
        
        if(iter>0.8*iter_max){
            for(i=0;i<n;i++){
                for(j=0;j<n;j++){
                    mean[(j*n)+i]+=(A[(j*n)+i]*A[(j*n)+i]);
                }
            }
        }

        // swap pointers
        Atmp = A; A = Anew; Anew = Atmp;
        stmp = s; s = snew; snew = stmp;
        rtmp = r; r = rnew; rnew = rtmp;
        ltmp = l; l = lnew; lnew = ltmp;
    }
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            mean[(j*n)+i]=sqrt(mean[(j*n)+i]/(0.2*iter_max));
        }
    }
    t_end = seconds();
    printf("Elapsed time (s) = %lf\n", t_end - t_start);

    save_gnuplot( mean, n );
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

    return 0;
}

void evolve( double * A, double *Anew, size_t n, double w, double c, double dt, size_t iter, size_t nc, double *s, double *r, double *l, double *snew, double *rnew, double *lnew ){
    size_t i,j;
    for( j = 1; j < n-1; j++) {
            for( i = 1; i < n-1; i++ ) {
                rnew[(j*n)+i]=r[(j*n)+i]+c*(0.5*(s[(j+1)*n+i]-s[(j-1)*n+i])+0.5*c*(r[(j+1)*n+i]-2*r[(j*n)+i]+r[(j-1)*n+i]));
                lnew[(j*n)+i]=l[(j*n)+i]+c*(0.5*(s[(j*n)+i+1]-s[(j*n)+i-1])+0.5*c*(l[(j*n)+i+1]-2*l[(j*n)+i]+l[(j*n)+i-1]));  
                if(j>=nc-2&&j<nc+2&&i>=nc-2&&i<nc+2){      
                    Anew[(j*n)+i]=1000*sin(dt*iter*w);
                    snew[(j*n)+i]=1000*w*cos(dt*iter*w);
                } else {
                    snew[(j*n)+i]=s[(j*n)+i]+c*(0.5*(r[(j+1)*n+i]-r[(j-1)*n+i]+l[(j*n)+i+1]-l[(j*n)+i-1])+0.5*c*(s[(j+1)*n+i]-4*s[(j*n)+i]+s[(j-1)*n+i]+s[(j*n)+i+1]+s[(j*n)+i-1]));
                    Anew[(j*n)+i]=A[(j*n)+i]+dt*0.5*(snew[(j*n)+i]+s[(j*n)+i]);
                }
                
            }
        }

        for(i=1;i<n-1;i++){
            Anew[i]=A[n+i];
            snew[i]=s[n+i];
            rnew[i]=0;
            lnew[i]=0;
            Anew[(n-1)*n+i]=A[(n-2)*n+i];
            snew[(n-1)*n+i]=s[(n-2)*n+i];
            rnew[(n-1)*n+i]=0;
            lnew[(n-1)*n+i]=0;
        }
        for(j=1;j<n-1;j++){
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

void save_gnuplot( double *M, size_t dimension ){
  
  size_t i , j;
  const double h = 0.1;
  FILE *file;

  file = fopen( "solution.dat", "w" );

  for( i = 0; i < dimension; ++i ){
    for( j = 0; j < dimension; ++j ) fprintf(file, "%f ", M[ ( i * ( dimension ) ) + j ] );
    fprintf(file, "\n" );
  }
  fclose( file );

}

// A Simple timer for measuring the walltime
double seconds(){

    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;
}