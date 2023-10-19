#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.plotstoredtable Plot data from a SKIRT stored table file
#
# This module offer functions to create plots from the data in a SKIRT stored table file.
#

# -----------------------------------------------------------------

import logging
import numpy as np
import matplotlib.pyplot as plt
import pts.simulation as sm
import pts.storedtable as stab
import pts.utils as ut

# -----------------------------------------------------------------

## This function creates a plot of a particular 1D curve from the data in a specified SKIRT stored table file.
# The function accepts the following optional arguments to configure the plot:
#  - horAxis: zero-based index of the table axis on the horizontal plot axis (default = 0)
#  - verAxis: zero-based index of the table quantity on the vertical plot axis (default = 0)
#  - axis0, axis1, axis2, axis3, axis4: ordinate value for the table axis with the indicated zero-based index
#    (default = arithmetic or geometric mean of the axis range; ignored for axes not in the table)
#
# Thus, by default, the script plots the first table quantity as a function of the first table axis,
# with half-way values for the other axes, if any.
#
# The table file path is interpreted as described for the pts.utils.absPath() function. By default, the figure
# is saved as FigStoredTable.pdf in the current directory. This can be overridden with the out* arguments
# as described for the pts.utils.savePath() function. In interactive mode (see the pts.utils.interactive() function),
# the figure is not saved and it is left open so that is displayed in notebooks.
#
def plotStoredTableCurve(tableFilePath, horAxis=0, verAxis=0, *,
                         axis0=None, axis1=None, axis2=None, axis3=None, axis4=None,
                         outDirPath=None, outFileName=None, outFilePath=None, figSize=(8, 6), interactive=None):

    # load the complete stored table
    table = stab.readStoredTable(tableFilePath)

    # get info on horizontal axis
    horName = table['axisNames'][horAxis]
    horUnit = table['axisUnits'][horAxis]
    horScale = table['axisScales'][horAxis]
    horGrid = table[horName]

    # get info on vertical axis
    verName = table['quantityNames'][verAxis]
    verUnit = table['quantityUnits'][verAxis]
    verScale = table['quantityScales'][verAxis]

    # get the appropriate slice from the values hypercube
    index = []
    for axisName,axisScale,axisValue in zip(table['axisNames'], table['axisScales'], (axis0,axis1,axis2,axis3,axis4)):
        if axisName == horName:
            index.append(Ellipsis)
        else:
            axisGrid = table[axisName]
            if axisValue is None:
                if axisScale == 'log': axisValue = np.sqrt(axisGrid[0] * axisGrid[-1])
                else: axisValue = (axisGrid[0] + axisGrid[-1]) / 2
            index.append((np.abs(axisGrid - axisValue)).argmin())
    verValues = table[verName][tuple(index)]

    # create the plot
    plt.figure(figsize=figSize)
    if horScale == 'log': plt.xscale('log')
    if verScale == 'log': plt.yscale('log')
    plt.plot(horGrid, verValues)
    plt.xlabel(horName + sm.latexForUnit(horUnit))
    plt.ylabel(verName + sm.latexForUnit(verUnit))

    # if not in interactive mode, save the figure; otherwise leave it open
    if not ut.interactive(interactive):
        saveFilePath = ut.savePath("FigStoredTable.pdf", (".pdf",".png"),
                                   outDirPath=outDirPath, outFileName=outFileName, outFilePath=outFilePath)
        plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        logging.info("Created {}".format(saveFilePath))

# -----------------------------------------------------------------

## This function creates an interactive plot of a 1D curve from the data in a specified SKIRT stored table file.
# It is intended for use from a notebook environment such as the one offered by Jupyter notebook.
# The stored table file path and the figure size are passed as arguments to this function. Other aspects
# of the plot, such as what to show on the horizontal and vertical axis, can be configured interactively via
# a user interface created through the ipywidgets package using the metadata in the stored table file.
#
# The table file path is interpreted as described for the pts.utils.absPath() function.
# The generated figure is not saved and is left open so that is displayed in notebooks.
#
# The function returns the stored table dictionary loaded by pts.storedtable.io.readStoredTable so that
# the interactive user can inspect its contents in further detail if desired.
#
def plotStoredTableInteractive(tableFilePath, *, figSize=(8, 6)):

    # import this here so that the dependency is limited to this function
    import ipywidgets

    # this function plots a particular slice of the data table hypercube
    # according to the configuration created through the widgets.interact() function
    def plottable(**args):
        # get info on horizontal axis
        horName = args['horaxis']
        horIndx = table['axisNames'].index(horName)
        horUnit = table['axisUnits'][horIndx]
        horScale = table['axisScales'][horIndx]
        horGrid = table[horName].value

        # get info on vertical axis
        verName = args['veraxis']
        verIndx = table['quantityNames'].index(verName)
        verUnit = table['quantityUnits'][verIndx]
        verScale = table['quantityScales'][verIndx]

        # get the appropriate slice from the values hypercube
        index = []
        for axisName in table['axisNames']:
            if axisName == horName:
                index.append(Ellipsis)
            else:
                axisValue = args[axisName]
                axisGrid = table[axisName].value
                index.append((np.abs(axisGrid - axisValue)).argmin())
        verValues = table[verName][tuple(index)].value

        # create the plot
        plt.figure(figsize=figSize)
        if horScale == 'log': plt.xscale('log')
        if verScale == 'log': plt.yscale('log')
        plt.plot(horGrid, verValues)
        plt.vlines([args[horName]], verValues.min(), verValues.max(), linestyle='--')
        plt.xlabel(horName + sm.latexForUnit(horUnit))
        plt.ylabel(verName + sm.latexForUnit(verUnit))
        plt.show()

    # load the complete stored table
    table = stab.readStoredTable(tableFilePath)

    # build a dictionary with a slider for each axis in the data hypercube
    axisSliders = {}
    for axisName, axisScale in zip(table['axisNames'], table['axisScales']):
        axisGrid = table[axisName]
        minValue = axisGrid.min().value
        maxValue = axisGrid.max().value
        if axisScale != 'log':
            axisSliders[axisName] = ipywidgets.FloatSlider(min=minValue, max=maxValue,
                                        value=minValue, readout_format='.3g', continuous_update=False)
        else:
            if minValue<=0: minValue = axisGrid[1].value  # avoid zero values for logarithmic treatment
            axisSliders[axisName] = ipywidgets.FloatLogSlider(min=np.log10(minValue), max=np.log10(maxValue),
                                            value=minValue, readout_format='.3e', continuous_update=False)

    # create the interactive plot
    ipywidgets.interact(plottable, horaxis=table['axisNames'], veraxis=table['quantityNames'], **axisSliders)
    return table

# ----------------------------------------------------------------------
