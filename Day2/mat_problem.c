#include <stdlib.h>
#include <stdio.h>

#define N 16

int main(int argc, char * argv[]){
    double * mat;
    mat = (double *) malloc(N*N*sizeof(double));
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            if(i==j){
                mat[i*N+j] = 1.0;
            }
            else{
                mat[i*N+j] = 0.0;
            }
        fprintf(stdout, "%.3g ", mat[i*N+j]);
        }
    fprintf(stdout, "\n");
    }
    free(mat);
    return 0;
}