#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int main(){
    int numberOfTimeSteps = 100;
    int data[100];
    srand (time(NULL));
    for (int i = 0; i < numberOfTimeSteps; i++) {
        data[i] = rand()*10;
    }

#ifdef CATALYST
    CatalystInit();
#endif
    for (int timeStep=0; timeStep<numberOfTimeSteps;timeStep++)
    {
#ifdef CATALYST
        CatalystCoProcess(timeStep, time,<grid info>,<field info>);
#endif

    }
#ifdef CATALYST
    CatalystFinalize();
#endif
}
