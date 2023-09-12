#include <math.h>
#include <string.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define NN 512
#define NM 512

double A[NN][NM];
double Aold[NN][NM];
double Anew[NN][NM];
double mean[NN][NM];

int main() {

    int n = NN;
    int m = NM;
    int nc = NM/2;
    int mc = NM/2;
//    int nv = 0; // Now the code forces it to be 0
    int nv = NN/3;
    int mv = NM/3;
    double v = 50; // velocity
    double dx = 0.1;//1.0 / (n - 1);
    double dt = 0.001;//dx / v * 0.5;
    double c = v * dt / dx;
    double c2=c*c;
    double k = (c-1.0)/(c+1.0);
    int iter_max = 10000;
    double tol = 1.0e-6;

    double error = 1.0;
    int iter = 0;
    int i,j;
    struct timeval temp_1, temp_2;
    double ttime=0.;

    memset(A, 0, n * m * sizeof(double));
    memset(Anew, 0, n * m * sizeof(double));
    memset(Aold, 0, n * m * sizeof(double));
    memset(mean, 0, n * m * sizeof(double));
    /*
    for (j = 0; j < n; j++)
    {
        A[j][0]    = 1.0;
        Anew[j][0] = 1.0;
    }
*/
    printf("Jacobi relaxation Calculation: %d x %d mesh\n", n, m);

    gettimeofday(&temp_1, (struct timezone*)0);

    #if defined(_OPENACC)
    printf("Using OpenACC\n");
    #pragma acc data copy(A), create(Anew)
    #endif
    while (iter < iter_max )
    {
        error = 0.0;
        #if defined(_OPENACC)
        #pragma acc parallel loop collapse(2) reduction(max:error)
        #endif
        for( j = 1; j < n-1; j++) {
            for( i = 1; i < m-1; i++ ) {
                if(j>=nc-5&&j<=nc+5&&i>=mc-5&&i<=mc+5){
                    Anew[j][i]=0;
                    error=0;
                // } else if(j==nv&&i==mv){
                //     Anew[j][i]=sin(dt*iter);
                //     error=fmax( error, fabs(Anew[j][i] - A[j][i]));
                } else {
                    Anew[j][i] = c2 *(0.5*A[j][i+1] + 0.5*A[j][i-1] + 0.5*A[j-1][i] + 0.5*A[j+1][i]
                                    +0.25*A[j+1][i+1]+0.25*A[j+1][i-1]+0.25*A[j-1][i+1]+0.25*A[j-1][i-1]
                                    -3*A[j][i])+2*A[j][i]-Aold[j][i];
                    //0.25 * ( A[j][i+1] + A[j][i-1] + A[j-1][i] + A[j+1][i]);
                    error = fmax( error, fabs(Anew[j][i] - A[j][i]));
                }
                
            }
        }
        
        //if(iter%200<100) A[0][mv]=100*sin(dt*iter*1000);
        //else A[0][mv]=0;
        
        // #if defined(_OPENACC)
        // #pragma acc parallel loop
        // #endif
        // for(i=1;i<m-1;i++){
        //     if(i>=mv-5&&i<=mv+5) A[0][i]=(iter%2==0)?0:1;
        //     else A[0][i]+=k*(Anew[1][i]-A[0][i]);
        //     A[n-1][i]+=k*(Anew[n-2][i]-A[n-1][i]);
        // }
        // #if defined(_OPENACC)
        // #pragma acc parallel loop
        // #endif
        // for(j=1;j<n-1;j++){
        //     A[j][0]+=k*(Anew[j][1]-A[j][0]);
        //     A[j][m-1]+=k*(Anew[j][m-2]-A[j][m-1]);
        // }
        #if defined(_OPENACC)
        #pragma acc parallel loop
        #endif
        for(i=1;i<m-1;i++){
            //if(i>=mv-5&&i<=mv+5) A[0][i]=(iter%200<100)?100*sin(dt*(iter%200)*1000):A[1][i];
            if(i==mv||i==m-mv) A[0][i]=100*sin(dt*iter*100);
            else A[0][i]=A[1][i];
            if(i==mv||i==m-mv) A[n-1][i]=100*sin(dt*iter*100);
            else A[n-1][i]=A[n-2][i];
        }
        #if defined(_OPENACC)
        #pragma acc parallel loop
        #endif
        for(j=1;j<n-1;j++){
            if(j==nv||j==n-nv) A[j][0]=100*sin(dt*iter*100);
            else A[j][0]=A[j][1];
            if(j==nv||j==n-nv) A[j][m-1]=100*sin(dt*iter*100);
            else A[j][m-1]=A[j][m-2];
        }
        //printf("%5d, %0.8f\n", iter, error);
        #if defined(_OPENACC)
        #pragma acc parallel loop collapse(2)
        #endif
        for( j = 1; j < n-1; j++) {
            for( i = 1; i < m-1; i++ ) {
                Aold[j][i] = A[j][i];
                A[j][i] = Anew[j][i];
            }
        }

        if(iter>0.8*iter_max){
            for(i=0;i<m;i++){
                for(j=0;j<n;j++){
                    mean[j][i]+=A[j][i]*A[j][i];
                }
            }
        }
        if(iter % 10 == 0&&iter<1000) //printf("%5d, %0.8lf\n", iter, error);
        {
            // print the grid to file laplace+iter for visualization
            char filename[64];
            snprintf(filename, 64, "laplace%d.dat", iter);
            FILE* fp = fopen(filename, "w");
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++) fprintf(fp, "%8.4lf ", A[j][i]);
                fprintf(fp, "\n");
            }
            fclose(fp);
            
        }
        iter++;
    }
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            mean[j][i]/=(0.2*iter_max);
        }
    }
    // print the final grid to file
    char filename[64];
    snprintf(filename, 64, "modonormal.dat");
    FILE* fp = fopen(filename, "w");
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < m; i++) fprintf(fp, "%8.4lf ", A[j][i]);
        fprintf(fp, "\n");
    }
    gettimeofday(&temp_2, (struct timezone*)0);
    ttime = 0.000001*((temp_2.tv_sec-temp_1.tv_sec)*1.e6+(temp_2.tv_usec-temp_1 .tv_usec));

    printf("Elapsed time (s) = %.2lf\n", ttime);
    printf("Stopped at iteration: %u\n", iter);

    return 0;

}
