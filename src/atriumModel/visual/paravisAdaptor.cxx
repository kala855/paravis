#include <iostream>
#include "paravisAdaptor.h"

#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPointData.h>


vtkCPProcessor* Processor = NULL;

namespace paravisAdaptor{
    void Initialize(int numScripts, char* scripts[]){
        if(Processor == NULL){
            Processor = vtkCPProcessor::New();
            Processor->Initialize();
        }

        for (int i = 4; i < numScripts; i++) {
            vtkCPPythonScriptPipeline* pipeline = vtkCPPythonScriptPipeline::New();
            pipeline->Initialize(scripts[i]);
            Processor->AddPipeline(pipeline);
            pipeline->Delete();
        }
    }

    void Finalize(){
        if(Processor){
            Processor->Delete();
            Processor = NULL;
        }
    }


    void CoProcess(int timeStep, double time, unsigned int numPoints[3],
            double spacing[3], double *potential){

        vtkCPDataDescription *dataDescription = vtkCPDataDescription::New();
        dataDescription->AddInput("input");
        dataDescription->SetTimeData(time,timeStep);

        if(Processor->RequestDataDescription(dataDescription) != 0){
            //Catalyst needs to output data
            //Create an axis-aligned, uniform grid
            vtkImageData* grid = vtkImageData::New();
            grid->SetExtent(0,numPoints[0]-1, 0, numPoints[1]-1, 0, 0);
            dataDescription->GetInputDescriptionByName("input")->SetGrid(grid);
            grid->Delete();
            //Create a  potential associated with points
            vtkDoubleArray* array = vtkDoubleArray::New();
            array->SetName("Potential");
            array->SetArray(potential, grid->GetNumberOfPoints(),1);
            grid->GetPointData()->AddArray(array);
            array->Delete();
            Processor->CoProcess(dataDescription);
        }
        dataDescription->Delete();
    }
}
