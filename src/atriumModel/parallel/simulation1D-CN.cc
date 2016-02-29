#include<bits/stdc++.h>
#include<armadillo>
#include "utilities/cell.h"
#include "arrayfire.h"

using namespace std;
using namespace arma;
#define db double
#define PI 3.14159265

// Matriz A for AX=B.
void load_matriz_A(mat &A, int N, db r){
  for(int x=0; x<N; x++){
    // Presente en todos los tiempos.
    // Tiempo i=1, se cumple condicion de borde en i=0.
    A(x,x) = (1.0+2.0*r);
    if(x==0){
      A(x,x+1) = -r;
    }
    // Tiempo i=N, se cumple condicion de borde en i=N+1.
    else if(x==N-1){ //time N+1
      A(x,x-1) = -r;
    }
    // Diagonal para los otros tiempos.
    else{
      A(x,x-1) = -r;
      A(x,x+1) = -r;
    }
  }
}

void copy_voltage_af(vector<Cell> &cells, af::array afX, double *X_mem){
  X_mem = afX.host<double>();
  for (int i = 0; i < cells.size(); i++)
    cells[i].V = X_mem[i];
}

// Copiar los voltages calculados por DF a cada celula.
void copy_voltage(vector<Cell> &cells, vec X){
  for(int i=0; i<cells.size(); i++)
    cells[i].V = X(i);
}

// Imprime los voltages calculados por DF en un tiempo determinado.
void print_solutions(vec &X,db t){
  //cout<<t;
  for(int i=0; i<X.n_rows; i++)
    cout<<" "<<X(i);
  cout<<endl;
}

int main(int argc, char *argv[]){
  db deltaX;      // como determinar este valor? // Longitud de una celula?
  db deltaX2;
  db D;
  int N;          // Numero de celulas a simular
  db nrepeat;     // numero de ciclos
  db tbegin;      // tiempo de inicio del primer estímulo
  db BCL;         // frecuencia de excitacion en mse
  db CI;          // intervalo de acoplamiento para el ultimo
  db dt;          // paso de tiempo
  db dtstim;      // duracion del estimulo
  db CurrStim;    // corriente de estimulo
  db cell_type;   // tipo de celula
  db ICCT;        // condiciones iniciales para el tipo de celula
  int nstp_prn;   // frecuencia con la que se imprimen los resultados
  db t=0.0;
  db tend;
  bool flag_stm = 1;
  db Istim=0.0;
  db nstp;
  db cont_repeat = 0;
  db Iion;
  db Jion; // current density
  db beta;
  db cell_to_stim;  // numero de la celula a estimular

  //temporal init values
  N = 15000;
  nrepeat = 1;
  tbegin = 50; //100; //50
  BCL =  50;//600;  //1000
  CI = 0;
  dt = 0.5;//0.01;
  dtstim = 2;
  CurrStim = -14;
  cell_type = 1;
  ICCT = 1;
  nstp_prn = 1;
  cell_to_stim = 1;
  // ---------------------

  tend = tbegin+dtstim;
  nstp = (tbegin+BCL*nrepeat+CI)/dt;

  vector<Cell> cells(N+2);
  db areaT = cells[0].pi*pow(cells[0].a,2);  // Capacitive membrane area
  db aCm = cells[0].Cap / areaT;             // Capacitance per unit area pF/cm^2

  D = cells[0].a / (2.0*cells[0].Ri*aCm*1e-9); //D = 0.00217147 cm^2/ms

  deltaX = 0.025;
  deltaX2 = deltaX*deltaX;

  mat A = mat(N,N);          // A
  vec B = vec(N,1);          // B
  vec X = vec(N);            // X from AX=B;

  db r = (D*dt)/(2.0*deltaX2);

  double *A_mem = (double*)malloc(N*N*sizeof(double));
  double *B_mem = (double*)malloc(N*sizeof(double));
  double *X_mem = (double*)malloc(N*sizeof(double));

  A_mem = A.memptr();
  B_mem = B.memptr();
  X_mem = X.memptr();

  // Se llena la matriz A
  load_matriz_A(A, N, r);

  int device = 0;
  af::setDevice(device);
  af::info();

  //Versiones arrayfire de la matriz A y los vectores B y X
  af::array afA(N,N,A_mem);
  af::array afB(N,f64);//se deben especificar los tipos de datos a operar
  af::array afX(N,f64);

  vec preV = vec(N+2);
  preV.fill(-81.2);

  // Versión arrayfire de preV
  af::array afPrev = af::constant(-81.2,N+2);

  nstp = -1; // only for one iteration

  for(int k=0; k < nstp+2; k++,t+=dt){       //all nodes time
    if(t>=tbegin && t<=tend){
      flag_stm = 0;
    }
    else{
      if(flag_stm==0){
        if(cont_repeat < nrepeat){
          tbegin=tbegin+BCL; //se establece el tiempo del próximo estimulo
        }
        else if(cont_repeat == nrepeat) tbegin=tbegin+CI;
        cont_repeat++;
        tend=tbegin+dtstim;
        flag_stm = 1;
      }
    }
    // Metodo Crank Nicholson
    for(int x=1; x<N+1; x++){
      // La corriente de estimulo es aplicada en un periodo de tiempo.
      if(!flag_stm && x == cell_to_stim)
        Istim = CurrStim;
      else
        Istim = 0.0;

      Iion = cells[x].getItot(dt);
      if(x==1){
        beta = (1.0-2.0*r)*preV(x) + r*preV(x+1) +2*r*preV(0);
      }
      else if(x==N){
        beta = r*preV(x-1) + (1.0-2.0*r)*preV(x) + 2*r*preV(N+1);
      }
      else{
        beta= r*preV(x-1) + (1.0-2.0*r)*preV(x) + r*preV(x+1);
      }
      Jion = (Iion+Istim)/areaT;
      afB(x-1) = beta - (dt/aCm)*Jion;
    }
    //ArrayFire Solver
    afX = af::solve(afA,afB);
    afPrev = afX;

    //Copiar voltajes a las celulas.
    //copy_voltage(cells,X);
    // copy_voltage_af(cells, afX, X_mem);
    //        print_solutions(X,t);
  }
  return 0;
}
