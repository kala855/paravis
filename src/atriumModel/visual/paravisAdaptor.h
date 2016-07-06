

namespace paravisAdaptor
{
    void Initialize(int numScripts, char* scripts[]);

    void Finalize();

    void CoProcess(int timeStep, double time, unsigned int numPoints[3],
            double spacing[3], double *potential);
}
