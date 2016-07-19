#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <cusolverSp.h>

int assert2(bool t, char *msg){
    if(t){
        printf("todo bien %s\n", msg);
    }
    else{
        printf("fail %s\n",msg);
    }
    return 0;
}


int main(){
    cudaSetDevice(0);
    cudaStream_t stream1;
    cudaStreamCreate(&stream1);
    cusolverSpHandle_t cusolverH;
    cudaError_t cudaStatus;
    cusparseStatus_t cusparse_status;// = CUSPARSE_STATUS_SUCCESS;
    cusolverStatus_t cusolver_status;//= CUSOLVER_STATUS_SUCCESS;

    cusparseMatDescr_t descrA;

    int *d_csrRowPtrA;
    int reorder = 0;
    int singularity = 0;
    double tol = 0.00001;
    int *d_csrColIndA;
    double *d_csrValA;
    double *d_b;
    double *d_x;
    double *x;


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
    /*const int csrRowPtrA[m+1] = {1, 2, 3, 4, 8};
    const int csrColIndA[nnzA] = {1, 2, 3, 1, 2, 3, 4};
    const double csrValA[nnzA] = {1.0, 2.0, 3.0, 0.1, 0.1, 0.1, 4.0};
    const double b[m] = {1.0, 1.0, 1.0, 1.0};*/

    int csrRowPtrA[m+1];
    csrRowPtrA[0] = 1;
    csrRowPtrA[1] = 2;
    csrRowPtrA[2] = 3;
    csrRowPtrA[3] = 4;
    csrRowPtrA[4] = 8;

    int csrColIndA[nnzA];
    csrColIndA[0] = 1;
    csrColIndA[1] = 2;
    csrColIndA[2] = 3;
    csrColIndA[3] = 1;
    csrColIndA[4] = 2;
    csrColIndA[5] = 3;
    csrColIndA[6] = 4;

    double csrValA[nnzA];
    csrValA[0] = 1.0;
    csrValA[1] = 2.0;
    csrValA[2] = 3.0;
    csrValA[3] = 0.1;
    csrValA[4] = 0.1;
    csrValA[5] = 0.1;
    csrValA[6] = 4.0;

    double b[m];
    b[0] = 1.0;
    b[1] = 1.0;
    b[2] = 1.0;
    b[3] = 1.0;

    x = (double *)malloc(m*sizeof(double));

    // Create cusolver handle, qr info and matrix descriptor
    cusolver_status = cusolverSpCreate(&cusolverH);
    assert2(cusolver_status == CUSOLVER_STATUS_SUCCESS, "cusolverSpCreate");


    cusolver_status = cusolverSpSetStream(cusolverH, stream1);
    assert2(cusolver_status == CUSOLVER_STATUS_SUCCESS, "Assign cusolver stream");

    cusparse_status = cusparseCreateMatDescr(&descrA);
    assert2(cusparse_status == CUSPARSE_STATUS_SUCCESS, "cusparseCreateMatDescr");

    cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatDiagType(descrA, CUSPARSE_DIAG_TYPE_NON_UNIT);
    cudaDeviceSynchronize();



    // copy A and b to device
    cudaStatus = cudaMalloc((void**)&d_csrValA, sizeof(double)*nnzA);
    assert2(cudaStatus == cudaSuccess, "cudaMalloc d_csrValA");

    cudaStatus = cudaMalloc((void**)&d_csrColIndA, sizeof(int)*nnzA);
    assert2(cudaStatus == cudaSuccess, "cudaMalloc d_csrColIndA");

    cudaStatus = cudaMalloc((void**)&d_csrRowPtrA, sizeof(int)*(m+1));
    assert2(cudaStatus == cudaSuccess, "cudaMalloc d_csrRowPtrA");

    cudaStatus = cudaMalloc((void**)&d_b, sizeof(double)*m);
    assert2(cudaStatus == cudaSuccess,"cudaMalloc d_b");

    cudaStatus = cudaMalloc((void**)&d_x, sizeof(double)*m);
    assert2(cudaStatus == cudaSuccess,"cudaMalloc d_x");

    cudaStatus = cudaMemcpy(d_csrValA, csrValA, sizeof(double)*nnzA, cudaMemcpyHostToDevice);
    assert2(cudaStatus == cudaSuccess,"cudaMemcpy csrValA");

    cudaStatus = cudaMemcpy(d_csrColIndA, csrColIndA, sizeof(int)*nnzA, cudaMemcpyHostToDevice);
    assert2(cudaStatus == cudaSuccess,"cudaMemcpy csrColIndA");

    cudaStatus = cudaMemcpy(d_csrRowPtrA, csrRowPtrA, sizeof(int)*(m+1), cudaMemcpyHostToDevice);
    assert2(cudaStatus == cudaSuccess,"cudaMemcpy csrRowPtrA");

    cudaStatus = cudaMemcpy(d_b, b, sizeof(double)*m, cudaMemcpyHostToDevice);
    assert2(cudaStatus == cudaSuccess,"cudamemcpy b");


    cusolver_status = cusolverSpDcsrlsvqr(cusolverH, m , nnzA, descrA, d_csrValA,
            d_csrRowPtrA, d_csrColIndA, d_b, tol, reorder, d_x, &singularity);

    assert2(cusolver_status == CUSOLVER_STATUS_SUCCESS,"Call cusolverSpDcsrlsvqr");
    cudaDeviceSynchronize();

    cudaStatus = cudaMemcpy(x, d_x, sizeof(double)*m, cudaMemcpyDeviceToHost);
    assert2(cudaStatus == cudaSuccess,"cudaMemcpy d_x");

    cudaFree(d_csrValA);cudaFree(d_csrColIndA);cudaFree(d_csrRowPtrA);cudaFree(d_b);cudaFree(d_x);
    free(x);

    return 0;

}
