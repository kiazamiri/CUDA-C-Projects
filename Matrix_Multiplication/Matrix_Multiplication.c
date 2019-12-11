// This program computes matrix multiplication
// By : Mohammadkia Zamiri
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <cassert>

using namespace std;

__global__ void matrixMul(int *a, int *b, int *c,int N) {
	// Calculate the global row and column for each thread
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	
	// Boundry check for our matrix
	if (row < N && col < N) {
		// Allocate partial results
		int temp = 0;
		for (int i = 0; i < N; i++) {
			temp += a[row * N + i] * b[i * N + col];
		}

		// Write back the result
		c[row * N + col] = temp;
	}

}
// Initilizes a square matrix with random numbers between 0-100
void init_matrix(int *m, int N) {
	for (int i = 0; i < N * N; i++) {
		m[i] = rand() % 100;
	}
}
// Verify the results on the CPU
void verify_result(int *a, int *b, int *c, int N) {
	int tmp;
	// For every row...
	for (int i = 0; i < N; i++) {
		// For every col....
		for (int j = 0; j < N; j++) {
			// For every element in row-col pair
			tmp = 0;
			for (int k = 0; k < N; k++) {
				tmp  += a[i * N + k] * b[k * N + j];
			}

			// Check each result
			assert(tmp == c[i * N + j]);
		}
	}

}
int main() {
	// Set our square matrix deminsion (2^10)*(2^10)
	int N = (1 << 10) + 1;
	size_t bytes = N * N * sizeof(int);
	// Allocate memory for our matrices
	int *a, *b, *c;
	cudaMallocManaged(&a, bytes);
	cudaMallocManaged(&b, bytes);
	cudaMallocManaged(&c, bytes);

	// Initialize our matrices
	init_matrix(a, N);
	init_matrix(b, N);

	// Set our block and Grid dimensions
	int threads = 16;
	int blocks = (N + threads - 1) / threads;

	// Set our kernel launch parameters 
	dim3 THREADS(threads, threads);
	dim3 BLOCKS(blocks, blocks); 

	// Launch our kernel
	matrixMul <<<BLOCKS, THREADS >>> (a,b,c,N);
	cudaDeviceSynchronize();


	//Verify the5 result
	verify_result(a, b, c, N);

	cout << "PROGRAM COMPLETED SUCCESSFULLY!" << endl;

	system("pause");
	return 0;
}