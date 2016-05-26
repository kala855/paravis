#include <stdio.h>

void CatalystInit(int numScripts, char * scripts[])
{
    if (Processor == NULL){
        Processor = vtkCPProcessor::New();
        Processor->Initialize();
    }
    // Scripts are Passed in as command line arguments
    for (int i = 0; i < numScripts; i++) {
        vtkCPPythonScriptPipeline* pipeline = vtkCPPythonScriptPipeline::New();
        pipeline->Initialize(scripts[i]);
        Processor->AddPipeline(pipeline);
        pipeline->Delete();
    }

}


void CatalystFinalize(){
    if(Processor){
        Processor->Delete();
        Processor = NULL;
    }
}


//The grid is a uniform, rectilinear grid that can be specified
//with the number of points in each direction and the uniform
//spacing between points. There is only one fiel called
//potential which is specified over the points/nodes of the grid
void CatalysCoProcess(int timeStep, double time, unsigned int numPoints[3], double spacing[3], double *field){
    vtkCPDataDescription* dataDescription = vtkCPDataDescription::New();
    dataDescription->AddInput("input");
    dataDescription->SetTimeData(time, timeStep);

    if(Processor->RequestDataDescription(dataDescription)!=0){
        // Catalyst needs to output data
        // Create an axis-aligned, uniform grid
        vtkImageData* grid = vtkImageData::New();
        grid->SetExtents(0,numPoints[0]-1, 0, numPoints[1]-1,0, numPoints[2]-1);
        dataDescription->GetInputDescriptionByName("input")->SetGrid(grid);
        grid->Delete();
        // Create a field associated with points
        vtkDoubleArray* array = vtkDoubleArray::New();
        array->SetName("potential");
        array->SetArray(field, grid->GetNumberOfPoints(),1);
        grid->GetPointData()->AddArray(array);
        array->Delete();
        Processor->CoProcess(dataDescription);
    }
    dataDescription->Delete();
}
