#include <stdio.h>
#include <cusolverSp.h>
#include <cuda_runtime_api.h>
#include <stdlib.h>
#include <assert.h>


int main(){
    cusolverSpHandle_t cusolverH = NULL;
    cudaError_t cudaStatus;
    cudasparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCES;
    cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;

    csrqrInfo_t info = NULL;
    cusparseMatDescr_t descrA = NULL;

    int *d_csrRowPtrA = NULL;
    int reorder = 0;
    int singularity;
    double tol = 0.0000001;
    int *d_csrColIndA = NULL;
    double *d_csrValA = NULL;
    double *d_b = NULL;
    double *d_x = NULL;
    double *x = NULL;

    size_t size_qr = 0;
    size_t size_internal = 0;
    void *buffer_qr = NULL;


/*    | 1     0    0   0 |
      | 0     2    0   0 |
  A=  | 0     0    3   0 |
      | 0.1  0.1  0.1  4 |
  CSR of A is based-1

  b = [1 1 1 1 ]
*/

    const int m = 4; // row size
    const int n = 4; // nxn matrix A
    const int nnzA = 7; // number of non-zero elements in A
    const int csrRowPtrA[m+1] = {1, 2, 3, 4, 8};
    const int csrColIndA[nnzA] = {1, 2, 3, 1, 2, 3, 4};
    const double csrValA[nnzA] = {1.0, 2.0, 3.0, 0.1, 0.1, 0.1, 4.0};
    const double b[m] = {1.0, 1.0, 1.0, 1.0};

    x = (double *)malloc(m*sizeof(double));

    // Create cusolver handle, qr info and matrix descriptor
    cusolver_status = cuSolverSpCreate(&cusolverH);
    assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

    cusparse_status = cusparseCreateMatDescr(&descrA);
    assert(cusparse_status == CUSPARSE_STATUS_SUCCES);

    // copy A and b to device
    cudaStatus = cudaMalloc((void**)&d_csrValA, sizeof(double)*nnzA);
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaMalloc((void**)&d_csrColIndA, sizeof(int)*nnzA);
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaMalloc((void**)&d_csrRowPtrA, sizeof(int)*(m+1));
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaMalloc((void**)&d_b, sizeof(double)*m);
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaMalloc((void**)&d_x, sizeof(double)*m);
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaMemcpy(d_csrValA, csrValA, sizeof(double)*nnzA, cudaMemcpyHostToDevice);
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaMemcpy(d_csrColIndA, csrColIndA, sizeof(int)*nnzA, cudaMemcpyHostToDevice);
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaMemcpy(d_csrRowPtrA, csrRowPtrA, sizeof(int)*(m+1), cudaMemcpyHostToDevice);
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaMemcpy(d_b, b, sizeof(double)*m, cudaMemcpyHostToDevice);
    assert(cudaStatus == cudaSuccess);


    cusolver_status = cusolverSpDcsrlsvlu(cusolverH, n, nnzA, descrA, d_csrValA, d_csrRowPtrA,
            d_csrColIndA, d_b, tol, reorder, d_x, &singularity);

    cudaStatus = cudaMemcpy(x, d_x, sizeof(double)*m, cudaMemcpyDeviceToHost);
    assert(cudaStatus == cudaSuccess);

    return 0;

}
