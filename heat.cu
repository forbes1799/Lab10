#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<math.h>

#define N 64
#define ITER_MAX 1000

int main(void){

    double tol = 1e-9;
    double error = tol + 1e-9;

    double * restrict A = malloc(N * N * sizeof(double));
    double * restrict Anew = malloc(N * N * sizeof(double));
	
	int x = 0;
	int radStart = floor((N-1)*0.3), radEnd = ceil((N-1)*0.7);

	int t0, t1;

	int pointx = floor((N-1)*0.5), pointy = floor((N-1)*0.5);
	
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
            int ij = i * N + j;
			Anew[ij] = 10;
            A[ij] = 10;
		}
	}

	for(int i = radStart; i <= radEnd; i++){
		Anew[i * N + (N-1)] = 100;
		A[i * N + (N-1)] = 100;
	}

    int iter = 0;

    double s_time = omp_get_wtime();
    
	while(error > tol && iter < ITER_MAX){
        error = 0.0;
		for(int i = 1; i < N-1; i++){
			for(int j = 1; j < N-1; j++){
                int ij = i * N + j;
                int ipj = (i + 1) * N + j;
                int imj = (i - 1) * N + j;
                int ijp = i * N + (j + 1);
                int ijm = i * N + (j - 1);
				Anew[ij] = 0.25 * (A[imj] + A[ipj] + A[ijm] + A[ijp]);
                error = fmax(error, fabs(A[ij] - Anew[ij]));
			}
		}
        for(int i = 0; i < N - 1; i++){
            for(int j = 0; j < N - 1; j++){
                int ij = i * N + j;
                A[ij] = Anew[ij];
            }
        }
        iter++;
    }


    double e_time = omp_get_wtime();
    
    printf("Time took %lf seconds", e_time - s_time); 
    

}
    


