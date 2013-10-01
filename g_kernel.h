// MyKernel.h
#ifndef g_kernel_h
#define g_kernel_h
#include <complex>
__global__
void g_kernel(cuFloatComplex w, cuFloatComplex gamma, cuFloatComplex* quasi_k, cuFloatComplex* delta_k, cuFloatComplex* k_weight, cuFloatComplex* g_wk, const int N);
void init_complex_array(cuFloatComplex *a, const int N);

void init_complex_array(cuFloatComplex *a, double* real_a, const int N);

void init_complex_array(cuFloatComplex *a, const float val, const int N);
void print_complex_array(cuFloatComplex *a, const int N, char *d);
cuFloatComplex sum_complex_array(cuFloatComplex *g_wk, const int N);

// template <class T>
// T call_g_kernel(T real_w, std::complex<T> complex_gamma, T* real_quasi_k, T* real_delta_k, T* real_k_weight, const int N);


// template <class T>
// T test_nate(int N);
void test_nate(int N);

std::complex<double> call_g_kernel(double real_w, std::complex<double> gamma, double* real_quasi_k, double* real_delta_k, double* real_k_weight, const int N);
#endif