#### import the simple module from the paraview
from paraview.simple import *
import os, os.path
import random
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'

def loadSource():
    folderTam = '100x100BCL1000-2'
    folder='/home/john/Documents/Projects/paravis/src/atriumModel/parallel/outputdata/'+folderTam
    DIR = folder
    print len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))])
    list_of_files= os.listdir(DIR)
    names = []
    names2 = []
    i = 0
    for file_name in list_of_files:
        names.append(folder+'/'+file_name)
        names2.append(folder+'/testParalelo'+str(i)+'.csv')
        i = i + 20
    return names2


names = loadSource()
testParalelo = CSVReader(FileName=names)


# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on testParalelo
testParalelo.HaveHeaders = 0

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [952, 919]

# get layout
viewLayout1 = GetLayout()

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# place view in the layout
viewLayout1.AssignView(2, spreadSheetView1)

# show data in view
testParaleloDisplay = Show(testParalelo, spreadSheetView1)
# trace defaults for the display properties.
testParaleloDisplay.FieldAssociation = 'Row Data'

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

# close an empty frame
viewLayout1.Collapse(2)

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=testParalelo)
tableToPoints1.XColumn = 'Field 0'
tableToPoints1.YColumn = 'Field 0'
tableToPoints1.ZColumn = 'Field 0'

# Properties modified on tableToPoints1
tableToPoints1.YColumn = 'Field 1'
tableToPoints1.a2DPoints = 1

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1)
# trace defaults for the display properties.
tableToPoints1Display.ColorArrayName = [None, '']

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [49.5, 49.5, 10000.0]
renderView1.CameraFocalPoint = [49.5, 49.5, 0.0]

# set scalar coloring
ColorBy(tableToPoints1Display, ('POINTS', 'Field 2'))

# rescale color and/or opacity maps used to include current data range
tableToPoints1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
tableToPoints1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Field2'
field2LUT = GetColorTransferFunction('Field2')
field2LUT.RGBPoints = [-81.2, 0.231373, 0.298039, 0.752941, -81.199594, 0.865003, 0.865003, 0.865003, -81.199188, 0.705882, 0.0156863, 0.14902]
field2LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Field2'
field2PWF = GetOpacityTransferFunction('Field2')
field2PWF.Points = [-81.2, 0.0, 0.5, 0.0, -81.199188, 1.0, 0.5, 0.0]
field2PWF.ScalarRangeInitialized = 1

# Rescale transfer function
field2LUT.RescaleTransferFunction(-81.2, 20.0)

# Rescale transfer function
field2PWF.RescaleTransferFunction(-81.2, 20.0)

# create a new 'Delaunay 2D'
delaunay2D1 = Delaunay2D(Input=tableToPoints1)

# set active source
SetActiveSource(delaunay2D1)

# show data in view
delaunay2D1Display = Show(delaunay2D1, renderView1)
# trace defaults for the display properties.
delaunay2D1Display.ColorArrayName = ['POINTS', 'Field 2']
delaunay2D1Display.LookupTable = field2LUT

# show color bar/color legend
delaunay2D1Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(tableToPoints1, renderView1)

# show data in view
delaunay2D1Display = Show(delaunay2D1, renderView1)

# reset view to fit data
renderView1.ResetCamera()

# hide data in view
Hide(tableToPoints1, renderView1)

# show color bar/color legend
delaunay2D1Display.SetScalarBarVisibility(renderView1, True)

# get animation scene
animationScene1 = GetAnimationScene()

# Properties modified on animationScene1
animationScene1.PlayMode = 'Real Time'

# Properties modified on animationScene1
animationScene1.Duration = 1

animationScene1.Play()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [49.5, 49.5, 270.47302994931886]
renderView1.CameraFocalPoint = [49.5, 49.5, 0.0]
renderView1.CameraParallelScale = 70.0035713374682

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
