#include <bits/stdc++.h>
#include <armadillo>
#include "cell.cuh"
#include "arrayfire.h"
#include <af/cuda.h>

using namespace std;
using namespace arma;

#define db double
#define PI 3.14159265

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

void load_Matrix_A(mat &A, int Nx, int Ny,db Sx, db Sy){
  db diag = 2.0*Sx + 2.0*Sy + 1.0;
  int node,upper,lower,prev,next;
  //Begin

  for(int i=1; i<=Ny; i++){        // iterate over rows
    for(int j=1; j<=Nx; j++){    //iterate over cols
      node = (i-1)*Nx + (j-1);
      upper = node + (Nx);
      lower = node - (Nx);
      prev = node - 1;
      next = node + 1;
      A(node, node) =  diag;

      if(j==1 && i==1){         // bottom-left
        A(node,next) = -Sx;
        A(node,upper) = -Sy;
      }else if(j==1 && i==Ny){  // top-left
        A(node,next) = -Sx;
        A(node,lower) = -Sy;
      }else if(j==1){           // left-middle
        A(node,next) = -Sx;
        A(node,lower) = -Sy;
        A(node,upper) = -Sy;
      }else if(j==Nx && i==1){  // bottom-right
        A(node,prev) = -Sx;
        A(node,upper) = -Sy;
      }else if(j==Nx && i==Ny){ // top-right
        A(node,prev) = -Sx;
        A(node,lower) = -Sy;
      }else if(i==1){           // bottom-middle
        A(node,prev) = -Sx;
        A(node,next) = -Sx;
        A(node,upper) = -Sy;
      }else if(j==Nx){          // right_middle
        A(node,prev) = -Sx;
        A(node,lower) = -Sy;
        A(node,upper) = -Sy;
      }else if(i==Ny){          // top-middle
        A(node,prev) = -Sx;
        A(node,next) = -Sx;
        A(node,lower) = -Sy;
      }else{                    // central nodes
        A(node,prev) = -Sx;
        A(node,next) = -Sx;
        A(node,lower) = -Sy;
        A(node,upper) = -Sy;
      }
    }
  }
}


int create_voltage_file(db t, af::array &X,int nodesA, int iteration){
    char file_name[50];
    sprintf(file_name,"testParalelo%d.csv",iteration);
    ofstream myfile;
    myfile.open(file_name,ios::app);
    for (int i = 0; i < nodesA; i++) {
        myfile << t << ","<< X(i).host<double>()[0] << endl;
    }

    myfile.close();
    return 0;
}


__global__ void d_update_B(int Nx, int Ny, int dt, db Sx, db Sy, db Istim, db CurrStim, db aCm,
        db areaT, int nodes, int flag_stm, db begin_cell, Cell *cells ,db *B, db *prevV){
    int node = (blockIdx.x*blockDim.x + threadIdx.x)+(Nx+3);
    if(node<(nodes-(Nx+3))){
        int pos = 0, i, j, upper, lower, prev, next;
        db Iion, Jion;
        db BC = 0;       // boundary condition
        db rhs = 0;      // rigth hand side

        upper = node + (Nx+2);
        lower = node - (Nx+2);
        prev = node - 1;
        next = node + 1;
        j = node % (Nx+2);        //pos in x -> cols
        i = node / (Nx+2);        //pos in y -> rows
        // Estimulando toda una fila de celulas
        //printf("j=%d\n",j);
        if(!flag_stm && (node >= begin_cell && node <= begin_cell + Nx -1)){
            Istim = CurrStim;
        }
        else{
            Istim = 0.0;
        }

        if(j>0 && j<(Nx+1)){
            int t = node-(Nx+3);
            pos = ((t/(Nx+2)) *Nx)+j - 1;
            //printf("%lf\n",prevV[prev]);
            Iion = cells[node].getItot(dt);
            //printf("ECa=%lf\tENa=%lf\tEK=%lf\tENC=%lf\n", cells[node].ECa,cells[node].ENa,cells[node].EK,cells[node].ENC);
            if(j==1 && i==1){                           //bottom-left
                BC = Sx * prevV[prev] + Sy * prevV[lower];
            }else if(i==1 && j==Nx){                    //bottom-right
                BC = Sx*prevV[next] + Sy*prevV[lower];
            }else if(i==1){                             //bottom-middle
                BC = Sy*prevV[lower];
            }else if(i==Ny && j==1){                    //top-left
                BC = Sx*prevV[prev] + Sy*prevV[upper];
            }else if(i==Ny && j==Ny){                   //top-right
                BC = Sx*prevV[next] + Sy*prevV[upper];
            }else if(i==Ny){                            //top-middle
                BC = Sy*prevV[upper];
            }else if(j==1){                             //left-middle
                BC = Sx*prevV[prev];
            }else if(j==Nx){                            //right-middle
                BC = Sx*prevV[next];
            }

            rhs = Sx*prevV[prev] + (1.0-2.0*Sx-2.0*Sy)* prevV[node] + Sx*prevV[next] + Sy*prevV[lower] + Sy*prevV[upper];
            Jion = (Iion + Istim)/areaT;
            B[pos] = rhs + BC - (Jion*dt/aCm);
            //printf("j = %d \t node = %d \t node-(Nx+3) = %d \t pos = %d\n", j, node, node-(Nx+3), pos);
           // printf("B[%d] = %lf\n", pos,B[pos]);
        }
    }
}



__global__ void d_copy_voltage(Cell *cells, db *X, db *prevV, int Nx, int size){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    if(i<size){
        int idx = Nx;
        idx = (i%Nx==0)? idx+3: idx+1;
        cells[idx].V = X[i];
        prevV[idx] = X[i];
    }
}

int copy_voltage(vector<Cell> &cells,af::array &X, af::array &prevV, int Nx, int size){
    int idx = Nx;
    for (int i = 0; i < size; i++) {
        idx = (i%Nx==0)? idx+3 : idx+1;
        cells[idx].V = X(i).host<double>()[0];
        prevV(idx) = X(i);
    }
    return 0;
}

__global__ void testeando(Cell *cells,int nodes){
    //cells[0].init();
    cells = new Cell[nodes]();
}



__global__ void init_d_prevv(int nodes, db *d_prevV){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    if(i<nodes)
        d_prevV[i] = -81.2;
}

__global__ void init_d_B(int nodesA, db *d_B){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    if(i<nodesA)
        d_B[i] = 0.0;
}

int main(){
  db deltaX,deltaY;
  int Nx,Ny,nodes,nodesA;
  db Dx,Dy,Gx,Gy;
  db Sx,Sy;
  db nrepeat;     // numero de ciclos
  db tbegin;      // tiempo de inicio del primer estímulo
  db BCL;         // frecuencia de excitacion en mse
  db CI;          // intervalo de acoplamiento para el ultimo
  db dt;          // paso de tiempo
  db dtstim;      // duracion del estimulo
  db CurrStim;    // corriente de estimulo
  db cell_type;   // tipo de celula
  int nstp_prn;   // frecuencia con la que se imprimen los resultados
  db tend;
  db nstp;
  int cell_to_stim;
  int upper,lower,prev,next, pos,i,j;
  db Iion;
  db Jion;
  db cont_repeat = 0;
  db t = 0.0;
  int flag_stm = 1;
  db Istim = 0.0;
//-------------------------------------
  nrepeat = 1;   //60-> 1min, 600-> 10min
  tbegin = 50; //100; //50
  BCL =  600;//600;  //1000
  CI = 0;
  dtstim = 2;
  CurrStim = -8000;
  cell_type = 1;
  nstp_prn = 50;
  tend = tbegin+dtstim;
//-------------------------------------
  Nx = 10;
  Ny = 10;
  cell_to_stim = 47;   // 70 in plot
  db row_to_stim = 1;
  db begin_cell = row_to_stim*(Nx+2) + 1;

  dt = 0.02;
  deltaX = deltaY = 0.025;
  nstp = (tbegin+BCL*nrepeat+CI)/dt;
  nodes = (Nx+2)*(Ny+2);              // nodes including boundary conditions
  nodesA = Nx*Ny;                 //nodes calculated in matrix A, no boundary conditions.

  vector<Cell> cells(nodes);


  db areaT = PINUM*pow(RADIUSCELL,2);  // Capacitive membrane area
  db aCm = CAP / areaT;             // Capacitance per unit area pF/cm^2
  Dx = Dy = RADIUSCELL / (2.0*RI*aCm*1e-9); //D = 0.00217147 cm^2/ms

  Sx = (dt*Dx)/(2.0*pow(deltaX,2));
  Sy = (dt*Dy)/(2.0*pow(deltaY,2));

//-------------------------------------
  mat A = mat(nodesA,nodesA);         // A
  vec B = vec(nodesA);                // B
  vec X = vec(nodesA);                // X from AX=B;
  vec prevV = vec(nodes);             // Voltages of T time
//-------------------------------------
  prevV.fill(-81.2);
  load_Matrix_A(A, Nx, Ny, Sx, Sy);

  //Additional ArrayFire Code
  int device = 0;
  af::setDevice(device);


  Cell *d_cells, *h_cells;
  const size_t sz = nodes * sizeof(Cell);
  h_cells = new Cell[nodes]();
  gpuErrchk(cudaMalloc((void**)&d_cells, sz));
  gpuErrchk(cudaMemcpy(d_cells, h_cells, sz, cudaMemcpyHostToDevice))
  int af_id = af::getDevice();
  cudaStream_t af_stream = afcu::getStream(af_id);
  /*testeando<<<1,1,0,af_stream>>>(d_cells, nodes);
  gpuErrchk(cudaPeekAtLastError());
  gpuErrchk(cudaStreamSynchronize(af_stream));*/

  int blockSize = 1024;
  dim3 dimGrid(ceil(float(nodes)/float(blockSize)),1,1);
  dim3 dimBlock(blockSize,1,1);

  dim3 dimGridCopyV(ceil(float(nodesA)/float(blockSize)),1,1);


 // af::info();

  double *A_mem = (double*)malloc(nodesA*nodesA*sizeof(double));
  double *B_mem = (double*)malloc(nodesA*sizeof(double));
  double *X_mem = (double*)malloc(nodesA*sizeof(double));

  A_mem = A.memptr();
  B_mem = B.memptr();
  X_mem = X.memptr();

  af::array afA(nodesA,nodesA,A_mem);
  af::array afB(nodesA,f64);
  af::array afX(nodesA,f64);
  af::array afBC(1,1);
  af::array afrhs(1,1);
  db *d_B;
  db *d_prevV;
  db *d_x;

  gpuErrchk(cudaMalloc((void**)&d_B,nodesA*sizeof(db)));
  gpuErrchk(cudaMalloc((void**)&d_prevV,nodes*sizeof(db)));

  init_d_prevv<<<dimGrid,dimBlock,0,af_stream>>>(nodes, d_prevV);
  gpuErrchk(cudaPeekAtLastError());
  gpuErrchk(cudaStreamSynchronize(af_stream));

  init_d_B<<<dimGridCopyV,dimBlock,0,af_stream>>>(nodesA, d_B);
  gpuErrchk(cudaPeekAtLastError());
  gpuErrchk(cudaStreamSynchronize(af_stream));

  af::array afALU, pivot;
  af::lu(afALU,pivot,afA);

  af::array afPrevV = af::constant(-81.2,nodes);

 //var for printing only the last ncharts beats
  int ncharts = 4;
  int time_to_print = nstp- ((ncharts*BCL+tbegin)/dt);

  nstp=0;  // only for one iteration

  for(int k=0; k<nstp+2; k++,t+=dt){ //each time
    if(t>=tbegin && t<=tend){
      flag_stm = 0.0;
    }else{
      if(flag_stm==0.0){
        if(cont_repeat < nrepeat){
          tbegin=tbegin+BCL; //se establece el tiempo del próximo estimulo
        }else if(cont_repeat == nrepeat) tbegin=tbegin+CI;

        cont_repeat++;
        tend=tbegin+dtstim;
        flag_stm = 1.0;
      }
    }

    d_update_B<<<dimGrid, dimBlock,0,af_stream>>>(Nx,Ny,dt,Sx,Sy,Istim,CurrStim,aCm,areaT,nodes,
            flag_stm,begin_cell,d_cells,d_B,d_prevV);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaStreamSynchronize(af_stream));

    //af::array afB(nodesA,d_B,afDevice);
    afB.write(d_B,nodesA*sizeof(db),afDevice);

    ////Array Fire Solver
    //afX = af::solve(afA,afB);
    afX = af::solveLU(afALU, pivot, afB);

    //printf("hi\n");
    d_x = afX.device<db>();
    d_copy_voltage<<<dimGridCopyV,dimBlock,0,af_stream>>>(d_cells,d_x,d_prevV,Nx,nodesA);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaStreamSynchronize(af_stream));
    afX.unlock();
    //copy_voltage(cells,afX,d_prevV,Nx,nodesA);
    /*if(k%nstp_prn==0 && k>time_to_print) //use this for plot last beat
        create_voltage_file(t,afX,nodesA,k);*/

  }
  cudaFree(d_cells);cudaFree(d_prevV);cudaFree(d_B);
  af_print(afX);
  return 0;
}
