#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <stdlib.h>


__global__ void hellofromGPU(void) {
	if (threadIdx.x == 5) {
		printf("Hello world from GPU thread %d!\n", threadIdx.x);
	}
}

int main()
{
	//hello world from CPU
	printf("Hello world from CPU!\n");

	hellofromGPU <<< 1, 10 >>> ();
	cudaDeviceSynchronize();
	return 0;
}