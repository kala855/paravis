#include<bits/stdc++.h>
#include<armadillo>
#include "utilities/cell.h"
#include "arrayfire.h"

using namespace std;
using namespace arma;

#define db double
#define PI 3.14159265

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

void copy_voltage(vector<Cell> &cells, vec &X, vec &prevV,int Nx,db dt){
  int idx = Nx;
  for(int i=0; i<X.size(); i++){
    idx =(i%Nx==0)? idx+3 : idx+1 ;
    cells[idx].dvdt = (X(i) - cells[idx].V)/dt;
    cells[idx].V = X(i);
    prevV(idx) = X(i);
  }
}

// Imprime los voltages calculados por DF en un tiempo determinado.
void print_solutions(vec &X,db t){
  cout<<t;
  for(int i=0; i<X.n_rows; i++)
    cout<<"  "<<X(i);
  cout<<endl;
}

void copy_prevV(vec &prevV,int Nx, vec X){
  for(int node=Nx+3,i=0; node<(prevV.n_rows-(Nx+3)); node++){
    if(node%(Nx+2)!=0 && (node+1)%(Nx+2)!=0){
      prevV(node) = X(i++);
    }
  }
}

void set_prevV(vec &prevV,int Nx, db value){
  for(int node=Nx+3,i=0; node<(prevV.n_rows-(Nx+3)); node++){
    if(node%(Nx+2)!=0 && (node+1)%(Nx+2)!=0){
      prevV(node) = value;
    }
  }
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
  nstp_prn = 20;
  tend = tbegin+dtstim;
//-------------------------------------
  Nx = 30;
  Ny = 30;
  cell_to_stim = 47;   // 70 in plot
  db row_to_stim = 4;
  db bengin_cell = row_to_stim*(Nx+2) + 1;

  dt = 0.02;
  deltaX = deltaY = 0.025;
  nstp = (tbegin+BCL*nrepeat+CI)/dt;
  nodes = (Nx+2)*(Ny+2);              // nodes including boundary conditions
  nodesA = Nx*Ny;                 //nodes calculated in matrix A, no boundary conditions.

  vector<Cell> cells(nodes);

  db areaT = cells[0].pi*pow(cells[0].a,2);  // Capacitive membrane area
  db aCm = cells[0].Cap / areaT;             // Capacitance per unit area pF/cm^2
  Dx = Dy = cells[0].a / (2.0*cells[0].Ri*aCm*1e-9); //D = 0.00217147 cm^2/ms

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

  //af::array afALU, pivot;
  //af::lu(afALU,pivot,afA);

  af::array afPrevV = af::constant(-81.2,nodes);

 //var for printing only the last ncharts beats
  int ncharts = 4;
  int time_to_print = nstp- ((ncharts*BCL+tbegin)/dt);

  nstp=-1;  // only for one iteration

  for(int k=0; k<1/*<nstp+2*/; k++,t+=dt){ //each time
    pos = 0;
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
    for(int node=Nx+3; node<(nodes-(Nx+3)); node++){
      db BC = 0;       // boundary condition
      db rhs = 0;      // rigth hand side

      upper = node + (Nx+2);
      lower = node - (Nx+2);
      prev = node - 1;
      next = node + 1;
      j = node % (Nx+2);        //pos in x -> cols
      i = node / (Nx+2);        //pos in y -> rows

      // Estimuando toda una fila de celulas
      if(!flag_stm && (node >= bengin_cell && node <= bengin_cell + Nx -1)){
        Istim = CurrStim;
      }
      else{
        Istim = 0.0;
      }

      if(j>0 && j<(Nx+1)){
        Iion = cells[node].getItot(dt);

        if(j==1 && i==1){                           //bottom-left
          afBC = Sx * afPrevV(prev) + Sy * afPrevV(lower);
        }else if(i==1 && j==Nx){                    //bottom-right
          afBC = Sx*afPrevV(next) + Sy*afPrevV(lower);
        }else if(i==1){                             //bottom-middle
          afBC = Sy*afPrevV(lower);
        }else if(i==Ny && j==1){                    //top-left
          afBC = Sx*afPrevV(prev) + Sy*afPrevV(upper);
        }else if(i==Ny && j==Ny){                   //top-right
          afBC = Sx*afPrevV(next) + Sy*afPrevV(upper);
        }else if(i==Ny){                            //top-middle
          afBC = Sy*afPrevV(upper);
        }else if(j==1){                             //left-middle
          afBC = Sx*afPrevV(prev);
        }else if(j==Nx){                            //right-middle
          afBC = Sx*afPrevV(next);
        }

        afrhs = Sx*afPrevV(prev) + (1.0-2.0*Sx-2.0*Sy)*afPrevV(node) + Sx*afPrevV(next) + Sy*afPrevV(lower) + Sy*afPrevV(upper);
        Jion = (Iion + Istim)/areaT;
        afB(pos++) = afrhs + afBC - (Jion*dt/aCm);
      }
    }
    ////Array Fire Solver
    afX = af::solve(afA,afB);
    //afX = af::solveLU(afALU,pivot,afB);
    afPrevV = afX;

    //copy_prevV(prevV,Nx,X);
    //copy_voltage(cells,X,prevV,Nx,dt);
    //if(k%nstp_prn==0)                //use this for plot all beats
  //  if(k%nstp_prn==0 && k>time_to_print)   //use this for plot last beat
    //  print_solutions(X,t);
  }
  //af_print(afX);
  return 0;
}
