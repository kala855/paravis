#### import the simple module from the paraview
from paraview.simple import *
import os, os.path
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


def setDisplayProperties(plotSelectionDisplay):
    for plotDisplay in plotSelectionDisplay:
        plotDisplay.AttributeType = 'Row Data'
        plotDisplay.UseIndexForXAxis = 0
        plotDisplay.XArrayName = 'Time'
        plotDisplay.SeriesVisibility = ['Field 0 (Row Statistics)', 'Field 1 (Row Statistics)']
        plotDisplay.SeriesLabel = ['N (Row Statistics)', 'N (Row Statistics)', 'Field 0 (Row Statistics)', 'Field 0 (Row Statistics)', 'Field 1 (Row Statistics)', 'Field 1 (Row Statistics)', 'vtkOriginalRowIds (Row Statistics)', 'vtkOriginalRowIds (Row Statistics)', 'Time (Row Statistics)', 'Time (Row Statistics)']
        plotDisplay.SeriesColor = ['N (Row Statistics)', '0', '0', '0', 'Field 0 (Row Statistics)', '0.89', '0.1', '0.11', 'Field 1 (Row Statistics)', '0.22', '0.49', '0.72', 'vtkOriginalRowIds (Row Statistics)', '0.3', '0.69', '0.29', 'Time (Row Statistics)', '0.6', '0.31', '0.64']
        plotDisplay.SeriesPlotCorner = ['N (Row Statistics)', '0', 'Field 0 (Row Statistics)', '0', 'Field 1 (Row Statistics)', '0', 'vtkOriginalRowIds (Row Statistics)', '0', 'Time (Row Statistics)', '0']
        plotDisplay.SeriesLineStyle = ['N (Row Statistics)', '1', 'Field 0 (Row Statistics)', '1', 'Field 1 (Row Statistics)', '1', 'vtkOriginalRowIds (Row Statistics)', '1', 'Time (Row Statistics)', '1']
        plotDisplay.SeriesLineThickness = ['N (Row Statistics)', '2', 'Field 0 (Row Statistics)', '2', 'Field 1 (Row Statistics)', '2', 'vtkOriginalRowIds (Row Statistics)', '2', 'Time (Row Statistics)', '2']
        plotDisplay.SeriesMarkerStyle = ['N (Row Statistics)', '0', 'Field 0 (Row Statistics)', '0', 'Field 1 (Row Statistics)', '0', 'vtkOriginalRowIds (Row Statistics)', '0', 'Time (Row Statistics)', '0']
        # Properties modified on plotSelectionOverTime1Display
        plotDisplay.SeriesVisibility = ['Field 1 (Row Statistics)']
        plotDisplay.SeriesColor = ['N (Row Statistics)', '0', '0', '0', 'Field 0 (Row Statistics)', '0.889998', '0.100008', '0.110002', 'Field 1 (Row Statistics)', '0.220005', '0.489998', '0.719997', 'vtkOriginalRowIds (Row Statistics)', '0.300008', '0.689998', '0.289998', 'Time (Row Statistics)', '0.6', '0.310002', '0.639994']
        plotDisplay.SeriesPlotCorner = ['Field 0 (Row Statistics)', '0', 'Field 1 (Row Statistics)', '0', 'N (Row Statistics)', '0', 'Time (Row Statistics)', '0', 'vtkOriginalRowIds (Row Statistics)', '0']
        plotDisplay.SeriesLineStyle = ['Field 0 (Row Statistics)', '1', 'Field 1 (Row Statistics)', '1', 'N (Row Statistics)', '1', 'Time (Row Statistics)', '1', 'vtkOriginalRowIds (Row Statistics)', '1']
        plotDisplay.SeriesLineThickness = ['Field 0 (Row Statistics)', '2', 'Field 1 (Row Statistics)', '2', 'N (Row Statistics)', '2', 'Time (Row Statistics)', '2', 'vtkOriginalRowIds (Row Statistics)', '2']
        plotDisplay.SeriesMarkerStyle = ['Field 0 (Row Statistics)', '0', 'Field 1 (Row Statistics)', '0', 'N (Row Statistics)', '0', 'Time (Row Statistics)', '0', 'vtkOriginalRowIds (Row Statistics)', '0']



def showData(plotSelection,quartileChartView):
    i = 0
    plotSelectionDisplay = []
    for plots in plotSelection:
        plotSelectionDisplay.append(Show(plots, quartileChartView))
        i = i + 1
    return plotSelectionDisplay

def createPlotSelection():
    plotSelection = []
    for i in range (0,10):
        selection = IDSelectionSource()
        selection.IDs = [0L,i+60]
        selection.FieldType = 'ROW'
        # create a new 'Plot Selection Over Time'
        plotSelection.append(PlotSelectionOverTime(Input=testParalelo,Selection=selection))
    return plotSelection

def loadSource():
    folderTam = '10x10'
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
        i = i + 50
    return names2

# find source
names = loadSource()
testParalelo = CSVReader(FileName=names)
testParalelo.HaveHeaders = 0

plotSelection = createPlotSelection()

# Create a new 'Quartile Chart View'
quartileChartView1 = CreateView('QuartileChartView')
quartileChartView1.ViewSize = [777, 1168]
quartileChartView1.LeftAxisRangeMaximum = 6.66
quartileChartView1.BottomAxisRangeMaximum = 6.66

# get layout
viewLayout1 = GetLayout()

# place view in the layout
viewLayout1.AssignView(0, quartileChartView1)

# show data in view
plotSelectionDisplay = showData(plotSelection,quartileChartView1)


setDisplayProperties(plotSelectionDisplay)

# Properties modified on quartileChartView1
quartileChartView1.ShowLegend = 0

# set active source
SetActiveSource(testParalelo)

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
