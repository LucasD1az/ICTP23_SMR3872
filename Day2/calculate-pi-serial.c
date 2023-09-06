#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 10000

int main(int argc, char * argv[]){
    int n = N;
    double w=1.0/n;
    double x,fx,sum=0.0;
    for(int i =1; i<=n; i++){
        x = (i-0.5)*w;
        fx = 4.0/(1.0+x*x);
        sum += fx;
    }
    fprintf(stdout, "The value of PI is %.10g vs %.10g\n", sum*w, M_PI);
    return 0;
}