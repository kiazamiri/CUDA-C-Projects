#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

__global__ void sumArrayDevice(float *a, float *b, float *c, const int N) {
	int i = threadIdx.x;

	if (i < N) c[i] = a[i] + b[i];
}

void initialData(float *ip, int size) {
	time_t t;
	srand((unsigned int)time(&t));

	for (int i = 0; i < size; i++) {
		ip[i] = (float)(rand() & 0xFF) / 10.0f;
	}
}

int main(int argc, char **argv) {
	int nElem = 1024;
	size_t nBytes = nElem * sizeof(float);

	float *h_a, *h_b, *h_c;
	h_a = (float*)malloc(nBytes);
	h_b = (float*)malloc(nBytes);
	h_c = (float*)malloc(nBytes);

	float *d_a, *d_b, *d_c;
	cudaMalloc((float**)&d_a, nBytes);
	cudaMalloc((float**)&d_b, nBytes);
	cudaMalloc((float**)&d_c, nBytes);

	initialData(h_a, nElem);
	initialData(h_b, nElem);

	cudaMemcpy(d_a, h_a, nBytes, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, h_b, nBytes, cudaMemcpyHostToDevice);
	
	dim3 block(nElem);
	
	sumArrayDevice <<<1, block >>>(d_a, d_b, d_c, nElem);

	cudaMemcpy(h_c, d_c, nBytes, cudaMemcpyDeviceToHost);

	for (int i = 0; i < 10; i++) {
		printf("%f \n", h_c[i]);
	}
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);

	return(0);
}