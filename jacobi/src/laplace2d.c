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
    double dx = 0.5 / (n - 1);
    double w = 5000;

    double dt = 0.000000025;//dx / v * 0.5;
    int iter_max=10/w/dt;    //itermax*dt*w=#ciclos = 10
    int d_iter=iter_max/50;
    double v =1700.0;//c*dx/dt;
    w *= 2*M_PI;
    double c = v * dt / dx;
    printf("c=%lf\n",c);

    int iter = 0;
    int i,j;
    struct timeval temp_1, temp_2;
    double ttime=0.;

    memset(A, 0, n * m * sizeof(double));
    memset(Anew, 0, n * m * sizeof(double));
    memset(mean, 0, n * m * sizeof(double));
    memset(r, 0, n * m * sizeof(double));
    memset(rnew, 0, n * m * sizeof(double));
    memset(s, 0, n * m * sizeof(double));
    memset(snew, 0, n * m * sizeof(double));
    memset(l, 0, n * m * sizeof(double));
    memset(lnew, 0, n * m * sizeof(double));
    for(j=nc-2;j<nc+2;j++)
        for(i=mc-2;i<mc+2;i++)
            s[j][i]=1000*w;
    printf("Chladni Plate Resonance Calculation: %d x %d mesh\n", n, m);

    gettimeofday(&temp_1, (struct timezone*)0);

    while (iter < iter_max )
    {
        for( j = 1; j < n-1; j++) {
            for( i = 1; i < m-1; i++ ) {
                rnew[j][i]=r[j][i]+c*(0.5*(s[j+1][i]-s[j-1][i])+0.5*c*(r[j+1][i]-2*r[j][i]+r[j-1][i]));
                lnew[j][i]=l[j][i]+c*(0.5*(s[j][i+1]-s[j][i-1])+0.5*c*(l[j][i+1]-2*l[j][i]+l[j][i-1]));  
                if(j>=nc-2&&j<nc+2&&i>=mc-2&&i<mc+2){      
                    Anew[j][i]=1000*sin(dt*iter*w);
                    snew[j][i]=1000*w*cos(dt*iter*w);
                } else {
                    snew[j][i]=s[j][i]+c*(0.5*(r[j+1][i]-r[j-1][i]+l[j][i+1]-l[j][i-1])+0.5*c*(s[j+1][i]-4*s[j][i]+s[j-1][i]+s[j][i+1]+s[j][i-1]));
                    Anew[j][i]=A[j][i]+dt*0.5*(snew[j][i]+s[j][i]);
                }
                
            }
        }

        for(i=1;i<m-1;i++){
            Anew[0][i]=A[1][i];
            snew[0][i]=s[1][i];
            rnew[0][i]=0;
            lnew[0][i]=0;
            Anew[n-1][i]=A[n-2][i];
            snew[n-1][i]=s[n-2][i];
            rnew[n-1][i]=0;
            lnew[n-1][i]=0;
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
            Anew[j][0]=Anew[j][1];
            Anew[j][m-1]=Anew[j][m-2];
            snew[j][0]=snew[j][1];
            snew[j][m-1]=snew[j][m-2];
            rnew[j][0]=0;
            rnew[j][m-1]=0;
            lnew[j][0]=0;
            lnew[j][m-1]=0;
        }

        if(iter>0.8*iter_max){
            for(i=0;i<m;i++){
                for(j=0;j<n;j++){
                    mean[j][i]+=(A[j][i]*A[j][i]);
                }
            }
        }

        for( j = 0; j < n; j++) {
            for( i = 0; i < m; i++ ) {
                A[j][i] = Anew[j][i];
                r[j][i]=rnew[j][i];
                l[j][i]=lnew[j][i];
                s[j][i]=snew[j][i];
            }
        }

/**/
        if(iter % d_iter == 0) //printf("%5d, %0.8lf\n", iter, error);
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
            
        }/**/
        iter++;
    }
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            mean[j][i]=sqrt(mean[j][i]/(0.2*iter_max));
        }
    }
    // print the final grid to file
    char filename[64];
    snprintf(filename, 64, "modonormal.dat");
    FILE* fp = fopen(filename, "w");
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < m; i++) fprintf(fp, "%8.4lf ", mean[j][i]);
        fprintf(fp, "\n");
    }
    gettimeofday(&temp_2, (struct timezone*)0);
    ttime = 0.000001*((temp_2.tv_sec-temp_1.tv_sec)*1.e6+(temp_2.tv_usec-temp_1 .tv_usec));

    printf("Elapsed time (s) = %.2lf\n", ttime);
    printf("Stopped at iteration: %u\n", iter);

    return 0;
}
