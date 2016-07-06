#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "paravisAdaptor.h"


using namespace paravisAdaptor;


void delay(int milliseconds)
{
    long pause;
    clock_t now,then;

    pause = milliseconds*(CLOCKS_PER_SEC/1000);
    now = then = clock();
    while( (now-then) < pause )
        now = clock();
}

int updateData(double* data, int size){
    srand (time(NULL));
    for (int i = 0; i < size; i++) {
        data[i] = rand()*10;
    }

}

int main(int argc, char* argv[]){
    int numberOfTimeSteps = 100;
    double spacing[3];
    spacing[0] = 0.1;
    spacing[1] = 0.1;
    spacing[2] = 0.0;
    unsigned int sizeData = 150*150;
    double *data = (double*)malloc(sizeData*sizeof(double));
    //double data[100];
    unsigned int numPoints[3];
    numPoints[0] = 150;
    numPoints[1] = 150;
    numPoints[2] = 0;
//#ifdef CATALYST
    Initialize(argc, argv);
//#endif
    for (int timeStep=0; timeStep<numberOfTimeSteps;timeStep++)
    {
        updateData(data, sizeData);
        delay(1000);
//#ifdef CATALYST
        CoProcess(timeStep,timeStep*0.1,numPoints,spacing,data);
//#endif

    }
//#ifdef CATALYST
    Finalize();
//#endif
}
