#include <stdio.h>
#include <cuda.h>
#include <cuComplex.h>
#include <time.h>
#include "g_kernel.h"
#include <complex>

__global__ 
void g_kernel(cuFloatComplex w, cuFloatComplex gamma, cuFloatComplex* quasi_k, cuFloatComplex* delta_k, cuFloatComplex* k_weight, cuFloatComplex* g_wk, const int N)
{
    int index=threadIdx.x+blockIdx.x*blockDim.x;
    if(index<N){
        cuFloatComplex wg = cuCaddf(w, gamma);
        cuFloatComplex t1 = cuCdivf(cuCmulf(delta_k[index], delta_k[index]), (cuCaddf(wg, quasi_k[index])) );
        g_wk[index] = cuCdivf(k_weight[index], cuCsubf(cuCsubf(wg,t1), quasi_k[index])  );
    }
}

//template <class T>
void test_nate(int N){
    printf("******** Calling test_nate *************");
}

// //double call_g_kernel(double real_w, complex<double> gamma, double* real_quasi_k, double* real_delta_k, double* real_k_weight, const int N)
// template void
// float call_g_kernel(float w, std::complex<float> complex_gamma, float* real_quasi_k, float* real_delta_k, float* real_k_weight, const int N);

// template void
// double call_g_kernel(double w, std::complex<double> complex_gamma, double* real_quasi_k, double* real_delta_k, double* real_k_weight, const int N);



double call_g_kernel(double real_w, std::complex<double> complex_gamma, double* real_quasi_k, double* real_delta_k, double* real_k_weight, const int N)
{
    cuFloatComplex w, gamma;
    cuFloatComplex *quasi_k, *delta_k, *k_weight, *g_wk;
    cuFloatComplex *dev_quasi_k, *dev_delta_k, *dev_k_weight, *dev_g_wk;

    w = make_cuFloatComplex(100.0,0);
    gamma = make_cuFloatComplex(10.0,0);
 

    quasi_k = (cuFloatComplex*)malloc( sizeof(cuFloatComplex)*N);
    delta_k = (cuFloatComplex*)malloc( sizeof(cuFloatComplex)*N);
    k_weight = (cuFloatComplex*)malloc( sizeof(cuFloatComplex)*N);
    g_wk = (cuFloatComplex*)malloc( sizeof(cuFloatComplex)*N);
    cudaMalloc((void**)&dev_quasi_k, sizeof(cuFloatComplex)*N);
    cudaMalloc((void**)&dev_delta_k, sizeof(cuFloatComplex)*N);
    cudaMalloc((void**)&dev_k_weight, sizeof(cuFloatComplex)*N);
    cudaMalloc((void**)&dev_g_wk, sizeof(cuFloatComplex)*N);

    init_complex_array(quasi_k, N);
    init_complex_array(delta_k, N);
    init_complex_array(k_weight, N);
    init_complex_array(g_wk, N);

    cudaMemcpy(dev_quasi_k, quasi_k, sizeof(cuFloatComplex)*N, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_delta_k, delta_k, sizeof(cuFloatComplex)*N, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_k_weight, k_weight, sizeof(cuFloatComplex)*N, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_g_wk, g_wk, sizeof(cuFloatComplex)*N, cudaMemcpyHostToDevice);
 
    g_kernel<<<N/256+1, 256>>>(w, gamma, dev_quasi_k, dev_delta_k, dev_k_weight, dev_g_wk, N);
    //printf("error code: %s\n",cudaGetErrorString(cudaGetLastError()));

    cudaMemcpy(g_wk, dev_g_wk, sizeof(cuFloatComplex)*N, cudaMemcpyDeviceToHost);

    cudaFree(dev_delta_k);
    cudaFree(dev_quasi_k);
    cudaFree(dev_k_weight);
    cudaFree(dev_g_wk);

    printf(">>>>>>>>>> final data:\n");
    print_complex_array(g_wk, N, "out-vector");

    return 0;
};

// template void
// float call_g_kernel(float w, std::complex<float> complex_gamma, float* real_quasi_k, float* real_delta_k, float* real_k_weight, const int N);

// template void
// double call_g_kernel(double w, std::complex<float> complex_gamma, double* real_quasi_k, double* real_delta_k, double* real_k_weight, const int N);


void init_complex_array(cuFloatComplex *a, const int N) {
    int i;
    int val = rand() % 4 + 1;
    for(i=0; i<N; i++)
        a[i] = make_cuFloatComplex(val,val);
}

void print_complex_array(cuFloatComplex *a, const int N, char *d){
    int i;
    for(i=0; i<N; i++)
        printf("\n%s[%d]: %f",d, i, cuCrealf(a[i]), cuCimagf(a[i]) );
    printf("\n");  
}
