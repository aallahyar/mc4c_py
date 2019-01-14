''' Present an interactive function explorer with slider widgets.
Scrub the sliders to change the properties of the ``sin`` curve, or
type into the title text box to update the title of the plot.
Use the ``bokeh serve`` command to run the example by executing:
	bokeh serve sliders.py
at your command prompt. Then navigate to the URL
	http://localhost:5006/sliders
in your browser.
'''
import numpy as np
import sys
from os import listdir
from os.path import isfile, join
from bokeh.io import curdoc
from bokeh.layouts import row, widgetbox, column
from bokeh.models import ColumnDataSource, Range1d,HoverTool, FactorRange, Spacer
from bokeh.models.widgets import MultiSelect, Slider, TextInput, Dropdown, DataTable, TableColumn, Button, RadioButtonGroup, RangeSlider
from bokeh.plotting import figure,ColumnDataSource
from bokeh.charts import Bar
from bokeh.charts.operations import blend
from bokeh.palettes import Colorblind, Greys, Set1

from time import sleep

from bokeh.charts.attributes import cat, color
import pandas


import operator as op

import dataTravel as dt

dataT = None

shownValues = []
shownLabels = []
shownIndexes = []

subSelection = []
subSelectionStr = []
subSelectionBox = []
subSelectionCount = 1

chromList = []

# Create a table to show statistics on the selected data
statisticsDict = dict(
    chromosomes=[],
    regions=[],
    selections=[],
    anchors=[]
)
tableSource = ColumnDataSource(statisticsDict)
statisticsCols = [
    TableColumn(field="chromosomes", title="Chromosome"),
    TableColumn(field="regions", title="Regions"),
    TableColumn(field="selections", title="In Selection"),
    TableColumn(field="anchors", title="Anchors"),
]
data_table = DataTable(source=tableSource, columns=statisticsCols, sortable=False)


# Load data upon selecting a dataset
mypath = '../pymc4c/data/'
outputDataFiles = [(f[:-8],join(mypath,f)) for f in listdir(mypath) if isfile(join(mypath, f)) and f.endswith('_out.npz')]
dropdown = Dropdown(label="Dataset", button_type="primary", menu=outputDataFiles)

def update_dataFile(new):
    dropdown.label='Loading...'
    dropdown.button_type='danger'
    global chromList, dataT
    dataT = dt.DataTraveller(new)

    chromList = [['',x,x] for x in dataT.byRegion.keys()]
    for x in chromList:
        if x[-1].startswith('chr'):
            x[-1]=x[-1][3:]
        if len(x[-1].split('_')[0]) == 1:
            try:
                int(x[-1][0])
                x[0] = '0'+x[-1]
            except ValueError:
                x[0] = x[-1]
        else:
            x[0] = x[-1]

    chromList.sort()
    #print chromList
    chromList = [x[1] for x in chromList]
    statisticsDict2 = dict(
        chromosomes=chromList,
        regions=[len(dataT.byRegion[i]) for i in chromList],
        selections=[len(dataT.byRegion[i]) for i in chromList],
        anchors=[0 for i in chromList]
    )
    tableSource.data = statisticsDict2

    buttonDrawPlot.label='Show selection'
    buttonDrawPlot.button_type='primary'
    buttonDrawPlot.on_click(update_dataShown)

    dropdown.label=new
    dropdown.button_type='default'

dropdown.on_click(update_dataFile)


def update_dataShown():

    # apply_subSelection()
    buttonDrawPlot.button_type='danger'

    global shownLabels,shownIndexes
    #for chromosome in byRegion:
    #	print chromosome, len(byRegion[chromosome])
    shownLabels = []
    shownValues = []

    shownData = dataT.getSelectedInfo(chromosomes=[chromList[x] for x in sorted(tableSource.selected['1d']['indices'])])
    shownIndexes = [(x[1],x[0]) for x in shownData]
    shownRange = range(len(shownData))
    shownLabels = [x[1]+':'+str(x[2])+'-'+str(x[3]) for x in shownData]
    shownVal = [x[4] for x in shownData] #shownValues#[7,9,5]
    shownSel = [x[5] for x in shownData]
    shownUnsel = [x[4]-x[5] for x in shownData]
    selData = [0]*len(shownVal)#[4,3,2]
    regStart = [x[2] for x in shownData]
    regEnd = [x[3] for x in shownData]


    df2=dict( \
        shownRange=shownRange,
        regionLabel=shownLabels, # This messes up if there is a : in a label...
        indexes=shownIndexes,
        chromosome=[x[1] for x in shownData],
        regStart = [x[2] for x in shownData],
        regEnd = [x[3] for x in shownData],
        selected=shownSel,
        unselected=shownUnsel,
        total=shownVal)

    source.data = df2

    global subSelectionBox
    subSelectionBox = []
    source.selected['1d']['indices']=[]

    update_viewFit()

    buttonDrawPlot.button_type='primary'


def update_viewFit(direction=None):
    start = 0
    end = len(source.data['selected'])
    if subSelectionBox != []:
        start=min(subSelectionBox)
        end=max(subSelectionBox)

    if direction == None or direction == 'x':
        bar2.x_range.start = start-0.5
        bar2.x_range.end = end+0.5


    #print start,end,max(source.data['selected'][start:end])
    #print source.data.keys()

    if direction == None or direction == 'y':
        bar2.y_range.start = 0.000001
        bar2.y_range.start = 0
        bar2.y_range.end = max(source.data['selected'][start:end])*1.2499999
        bar2.y_range.end = max(source.data['selected'][start:end])*1.25

def apply_subSelection():
    global subSelectionCount
    buttonApply.button_type='danger'
    if subSelection != []:
        dataT.applyAnchor(subSelection)
        selectSelHist.options.append('--- Group '+str(subSelectionCount)+' ---')
        selectSelHist.options.extend(selectSubSel.options)
        reset_subSelection()
    subSelectionCount+=1
    update_dataShown()
    #source.selected['1d']['indices']=[]

    tmpAnchors = dataT.getAppliedAnchors()
    tableSource.data['anchors']=[tmpAnchors[x] for x in tableSource.data['chromosomes']]

    tmpRegions = dataT.getChromosomeRegions()
    tableSource.data['selections']=[tmpRegions[x] for x in tableSource.data['chromosomes']]

    # print dataT.getChromosomeRegions()
    # print tableSource.data['chromosomes']
    buttonApply.button_type='primary'


def update_subSelection():
    subSelection.extend([shownIndexes[x] for x in subSelectionBox])
    subSelectionStr.extend([shownLabels[x] for x in subSelectionBox])
    selectSubSel.options=subSelectionStr


def reset_subSelection():
    global subSelection, subSelectionStr
    subSelection = []
    subSelectionStr = []
    selectSubSel.options=subSelectionStr

selectRange = RangeSlider(start=0, end=1, range=(0,1), step=1, title="Select")

buttonFitView = Button(label="Fit selection to view", button_type="default")
buttonFitView.on_click(update_viewFit)

buttonSubSelect = Button(label="Add Subselection", button_type="success")
buttonSubReset = Button(label="Reset Subselection", button_type="warning")
buttonApply = Button(label="Apply Subselection as Anchors", button_type="primary")

buttonSubSelect.on_click(update_subSelection)
buttonSubReset.on_click(reset_subSelection)
buttonApply.on_click(apply_subSelection)

buttonDrawPlot = Button(label="Load a dataset first", button_type="default")

selectSubSel = MultiSelect(title="Current Subselection:",
                           options=[])

selectSelHist = MultiSelect(title="Applied Anchors:",
                            options=[])

# Set up layouts and add to document
inputsData = widgetbox(dropdown,buttonDrawPlot,data_table,sizing_mode="scale_width")
inputsSelect = widgetbox(
    buttonFitView,
    buttonSubSelect,
    buttonSubReset,
    buttonApply,
    selectSubSel,
    selectSelHist,
    sizing_mode="scale_width")

#colorScheme = ['#DDDDDD','#E08787','#B3E087','#87E0E0','#B387E0']
colorScheme = Colorblind[8]
colorScheme = Set1[9]
greyScheme = Greys[3]

df=dict( \
    shownRange=[0],
    regionLabel=[0],
    indexes=[0],
    selected=[1],
    unselected=[2],
    total=[3])
data = pandas.DataFrame(df)#.assign(unselected=lambda df: df['total'] - df['selected'])
source = ColumnDataSource(data)
#print data['regionLabel'].tolist()

def selection_change(attrname, old, new):
    global subSelectionBox
    subSelectionBox = new['1d']['indices']
    # print 'dafuq'
    update_viewFit(direction='y')
    return

source.on_change('selected', selection_change)

tooltips = [
    ("Chromosome", "@chromosome"),
    ("Start", "@regStart"),
    ("End","@regEnd"),
    ("#Total", "@total"),
    ("#Selected", "@selected"),
    ("#Unselected", "@unselected"),
]

hover_bar = HoverTool(tooltips=tooltips)
TOOLS_bar = [hover_bar,'save','xbox_select']#,'xbox_zoom','xpan','reset']

bar2 = figure(width=800, height=150,
              x_range=(-1,1),#data['regionLabel'].tolist(),
              tools=TOOLS_bar,
              min_border=10,
              title='Welcome...')

bar2.xaxis.visible = False
bar2.xgrid.grid_line_color = None
bar2.xaxis.major_label_orientation = np.pi/3

pt1 = bar2.vbar(x='shownRange', bottom=0, top='selected', width=1.001,
                source=source,
                color=colorScheme[0],
                alpha=0.8,
                selection_color=colorScheme[2],nonselection_color=colorScheme[1],
                selection_line_color=colorScheme[2],nonselection_line_color=colorScheme[1],
                selection_alpha=0.8, nonselection_alpha=0.6)
pt2 = bar2.vbar(x='shownRange', bottom='selected', top='total', width=1.001,
                source=source,
                color=greyScheme[1],
                alpha=0.8,
                selection_color=greyScheme[0],nonselection_color=greyScheme[1],
                selection_line_color=None,nonselection_line_color=None,
                selection_alpha=0.5, nonselection_alpha=0.5)

ct1 = bar2.circle('shownRange', 1, source=source, color=None, selection_color=None, nonselection_color=None)

#bar2.y_range.end=20
#print dir(source.selected)

subBar1 = figure(width=800, height=50,
                 x_range=bar2.x_range,
                 tools=['xpan','xwheel_zoom'],
                 min_border=10)
subBar1.xaxis.visible = False
subBar1.xgrid.grid_line_color = None

subBar1.vbar(x='shownRange', bottom=0, top='selected', width=1,
             color=colorScheme[4], source=source,
             selection_color=colorScheme[5],nonselection_color=colorScheme[6],
             selection_line_color=None,nonselection_line_color=None,
             selection_alpha=0.8, nonselection_alpha=0.6)

# pt1.data_source.on_change('selected', selection_change)

# layout = column(row(bar2,sizing_mode='scale_width'), row(subBar1, Spacer(width=200, height=200),sizing_mode='scale_width'))
# curdoc().add_root(layout)
curdoc().add_root(row(bar2,sizing_mode="scale_width"))
curdoc().add_root(row(subBar1,sizing_mode="scale_width"))
curdoc().add_root(row(inputsData,inputsSelect,sizing_mode='scale_width'))
curdoc().title = "Sliders"
