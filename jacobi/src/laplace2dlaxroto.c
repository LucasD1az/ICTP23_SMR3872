#include <math.h>
#include <string.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define NN 128
#define NM 128

double A[NN][NM];
double Anew[NN][NM];
double Aold[NN][NM];
double mean[NN][NM];
double r[NN][NM];
double rnew[NN][NM];
double s[NN][NM];
double snew[NN][NM];
double l[NN][NM];
double lnew[NN][NM];

int main() {

    int n = NN;
    int m = NM;
    int nc = NM/2;
    int mc = NM/2;
//    int nv = 0; // Now the code forces it to be 0
    int nv = NN/3;
    int mv = NM/3;
    double dx = 0.5 / (n - 1);
    double lambda = 4.0;//dx*200.0;
    double dt = 0.00000001;//dx / v * 0.5;
    double v =1700.0;//c*dx/dt;
    double w = v/lambda;
//    double lambda = v/w;
    w *= 2*M_PI;
    double c = 0.005;//v * dt / dx;
    printf("c=%lf\n",c);
    double c2=c*c;
    double k = -(c-1.0)/(c+1.0);
    int iter_max = 100000;
    double tol = 1.0e-6;

    double error = 1.0;
    int iter = 0;
    int i,j;
    struct timeval temp_1, temp_2;
    double ttime=0.;

    memset(A, 0, n * m * sizeof(double));
    memset(Anew, 0, n * m * sizeof(double));
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
        for( j = 4; j < n-4; j++) {
            for( i = 4; i < m-4; i++ ) {
                if(j>=nc-2&&j<nc+2&&i>=mc-2&&i<mc+2){
                    // Anew[j][i]=0;
                    Anew[j][i]=1000*sin(dt*iter*w);
                    //snew[j][i]=1000*w*cos(dt*iter*w);
                    error=0;
                // } else if(j==nv&&i==mv){
                //     Anew[j][i]=sin(dt*iter);
                //     error=fmax( error, fabs(Anew[j][i] - A[j][i]));
                } else {
                            Anew[j][i]  = c2
                                    * ( -  1.0/560 * A[j][i-4]// u[1, 4:dimx-4, 0:dimy-8]  # c,   r-4
                                        +  8.0/315 * A[j][i-3]// u[1, 4:dimx-4, 1:dimy-7]  # c,   r-3
                                        -    1.0/5 * A[j][i-2]//u[1, 4:dimx-4, 2:dimy-6]  # c,   r-2
                                        +    8.0/5 * A[j][i-1]//u[1, 4:dimx-4, 3:dimy-5]  # c,   r-1

                                        - 1.0/560  * A[j-4][i]//u[1, 0:dimx-8, 4:dimy-4]  # c - 4, r
                                        + 8.0/315  * A[j-3][i]//u[1, 1:dimx-7, 4:dimy-4]  # c - 3, r
                                        -   1.0/5  * A[j-2][i]//u[1, 2:dimx-6, 4:dimy-4]  # c - 2, r
                                        +   8.0/5  * A[j-1][i]//u[1, 3:dimx-5, 4:dimy-4]  # c - 1, r
                                        - 410.0/72 * A[j][i]//u[1, 4:dimx-4, 4:dimy-4]  # c    , r
                                        +   8.0/5  * A[j+1][i]//u[1, 5:dimx-3, 4:dimy-4]  # c + 1, r
                                        -   1.0/5  * A[j+2][i]//u[1, 6:dimx-2, 4:dimy-4]  # c + 2, r
                                        + 8.0/315  * A[j+3][i]//u[1, 7:dimx-1, 4:dimy-4]  # c + 3, r
                                        - 1.0/560  * A[j+4][i]//u[1, 8:dimx  , 4:dimy-4]  # c + 4, r

                                        +    8.0/5 * A[j][i+1]//u[1, 4:dimx-4, 5:dimy-3]  # c  , r+1
                                        -    1.0/5 * A[j][i+2]//u[1, 4:dimx-4, 6:dimy-2]  # c  , r+2
                                        +  8.0/315 * A[j][i+3]//u[1, 4:dimx-4, 7:dimy-1]  # c  , r+3
                                        -  1.0/560 * A[j][i+4]//u[1, 4:dimx-4, 8:dimy  ]  # c  , r+4
                                        )
                                    + 2.0*A[j][i]//u[1, 4:dimx-4, 4:dimy-4] 
                                    -   Aold[j][i];//u[2, 4:dimx-4, 4:dimy-4]  
                    //Anew[j][i] = c2 *(0.5*A[j][i+1] + 0.5*A[j][i-1] + 0.5*A[j-1][i] + 0.5*A[j+1][i]
                    //                +0.25*A[j+1][i+1]+0.25*A[j+1][i-1]+0.25*A[j-1][i+1]+0.25*A[j-1][i-1]
                    //                -3*A[j][i])+2*A[j][i]-Aold[j][i];
                    //0.25 * ( A[j][i+1] + A[j][i-1] + A[j-1][i] + A[j+1][i]);
                    //Anew[j][i] = 0.995*Anew[j][i]+0.005*Aold[j][i];
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
            for(j=3;j>=0;j--){
                // if(((i>=mv-3&&i<=mv+3)||(i>=m-mv-3&&i<=m-mv+3))) {
                //     Anew[j][i]=1000*sin(dt*iter*w);
                //     snew[j][i]=1000*w*cos(dt*iter*w);
                // }
                // else 
                Anew[j][i]=A[j+1][i];//+k*(Anew[j+1][i]-A[j][i]);
                // Anew[j][i]=0;
                // snew[j][i]=0;

                // A[j][i]=A[j+1][i];
            }
            for(j=n-4;j<n;j++){
                // if(((i>=mv-3&&i<=mv+3)||(i>=m-mv-3&&i<=m-mv+3))) {
                //     Anew[j][i]=1000*sin(dt*iter*w);
                //     snew[j][i]=1000*w*cos(dt*iter*w);
                // }
                // else 
                Anew[j][i]=A[j-1][i];//+k*(Anew[j-1][i]-A[j][i]);
                // Anew[j][i]=0;
                // snew[j][i]=0;
                // A[j][i]=A[j-1][i];
            }
        }
        for(j=1;j<n-1;j++){
            for(i=3;i>=0;i--){
                // if(((j>=nv-3&&j<=nv+3)||(j>=n-nv-3&&j<=n-nv+3))) {
                //     Anew[j][i]=1000*sin(dt*iter*w);
                //     snew[j][i]=1000*w*cos(dt*iter*w);
                // }
                // else 
                Anew[j][i]=A[j][i+1];//+k*(Anew[j][i+1]-A[j][i]);
                // Anew[j][i]=0;
                // snew[j][i]=0;
            //    A[j][i]=A[j][i+1];
            }
            for(i=m-4;i<m;i++){
                // if(((j>=nv-3&&j<=nv+3)||(j>=n-nv-3&&j<=n-nv+3))) {
                //     Anew[j][i]=1000*sin(dt*iter*w);
                //     snew[j][i]=1000*w*cos(dt*iter*w);
                // }
                // else 
                Anew[j][i]=A[j][i-1];//+k*(Anew[j][i-1]-A[j][i]);
                // Anew[j][i]=0;
                // snew[j][i]=0;
            //    A[j][i]=A[j][i-1];
            }
        }
        // A[0][mv]=100*sin(dt*iter*100);
        // A[n-1][mv]=100*sin(dt*iter*100);
        // A[0][m-mv]=100*sin(dt*iter*100);
        // A[n-1][m-mv]=100*sin(dt*iter*100);
        // A[nv][0]=100*sin(dt*iter*100);
        // A[nv][m-1]=100*sin(dt*iter*100);
        // A[n-nv][0]=100*sin(dt*iter*100);
        // A[n-nv][m-1]=100*sin(dt*iter*100);

        // for(i=1;i<m-1;i++){
        //     //if(i>=mv-5&&i<=mv+5) A[0][i]=(iter%200<100)?100*sin(dt*(iter%200)*1000):A[1][i];
        //     if(i==mv||i==m-mv) A[0][i]=100*sin(dt*iter*100);
        //     else A[0][i]=A[1][i];
        //     if(i==mv||i==m-mv) A[n-1][i]=100*sin(dt*iter*100);
        //     else A[n-1][i]=A[n-2][i];
        // }
        // #if defined(_OPENACC)
        // #pragma acc parallel loop
        // #endif
        // for(j=1;j<n-1;j++){
        //     if(j==nv||j==n-nv) A[j][0]=100*sin(dt*iter*100);
        //     else A[j][0]=A[j][1];
        //     if(j==nv||j==n-nv) A[j][m-1]=100*sin(dt*iter*100);
        //     else A[j][m-1]=A[j][m-2];
        // }
        //printf("%5d, %0.8f\n", iter, error);
        #if defined(_OPENACC)
        #pragma acc parallel loop collapse(2)
        #endif
        for( j = 0; j < n; j++) {
            for( i = 0; i < m; i++ ) {
                A[j][i] = Anew[j][i];
                // r[j][i]=rnew[j][i];
                // l[j][i]=lnew[j][i];
                // s[j][i]=snew[j][i];
            }
        }

        if(iter>0.5*iter_max){
            for(i=0;i<m;i++){
                for(j=0;j<n;j++){
                    mean[j][i]+=A[j][i]*A[j][i];
                }
            }
        }
        if(iter % 1000 == 0) //printf("%5d, %0.8lf\n", iter, error);
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
