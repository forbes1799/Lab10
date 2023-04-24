#include <stdio.h>
#include <stdlib.h>

int N =10;

__global__ void squareKernel(int N) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < N) {
    int square = i * i;
    printf("GPU thread %d of block %d (dim: %d): value %d, square %d\n", threadIdx.x, blockIdx.x, blockDim.x, i, square);
  }
}


int main(void) {
  int numGPUs;
  printf("hello from CPU\n");
  
  /* check have GPU else quit */
  cudaGetDeviceCount(&numGPUs);
  
  if (numGPUs >= 1 ) {
  /* Synchronize */
  cudaDeviceSynchronize();
  /* call GPU kernel  */
  squareKernel<<<4,3>>> (N);
  /* Synchronize */
  cudaDeviceSynchronize();
  }
  
  else {
   printf("no GPU present\n");
   }
   printf("hello again from CPU\n");
}


