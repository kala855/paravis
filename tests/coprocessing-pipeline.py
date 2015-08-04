
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 4.1.0-1257-gc87f9d4 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 4.1.0-1257-gc87f9d4

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView0 = CreateView('RenderView')
      renderView0.ViewSize = [868, 781]
      renderView0.CameraPosition = [0.0, 0.0, 66.92130429902464]
      renderView0.CameraParallelScale = 17.320508075688775
      renderView0.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView0,
                               filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=868, height=781, cinema={})

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'Wavelet'
      # create a producer from a simulation input
      wavelet0 = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Contour'
      contour0 = Contour(Input=wavelet0)
      contour0.ContourBy = ['POINTS', 'RTData']
      contour0.Isosurfaces = [157.0909652709961]
      contour0.PointMergeMethod = 'Uniform Binning'

      # create a new 'Parallel Image Data Writer'
      parallelImageDataWriter0 = servermanager.writers.XMLPImageDataWriter(Input=wavelet0)

      # register the writer with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the data, etc.
      coprocessor.RegisterWriter(parallelImageDataWriter0, filename='filename_%t.pvti', freq=1)

      # create a new 'Parallel PolyData Writer'
      parallelPolyDataWriter0 = servermanager.writers.XMLPPolyDataWriter(Input=contour0)

      # register the writer with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the data, etc.
      coprocessor.RegisterWriter(parallelPolyDataWriter0, filename='filename_%t.pvtp', freq=1)

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView0'
      # ----------------------------------------------------------------

      # show data from wavelet0
      wavelet0Display = Show(wavelet0, renderView0)
      # trace defaults for the display properties.
      wavelet0Display.Representation = 'Outline'
      wavelet0Display.ColorArrayName = ['POINTS', '']
      wavelet0Display.ScalarOpacityUnitDistance = 1.7320508075688779
      wavelet0Display.Slice = 10

      # show data from contour0
      contour0Display = Show(contour0, renderView0)
      # trace defaults for the display properties.
      contour0Display.ColorArrayName = [None, '']

    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1, 1, 1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor

#--------------------------------------------------------------
# Global variables that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView
coprocessor.EnableLiveVisualization(True, 1)


# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
