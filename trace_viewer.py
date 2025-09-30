import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets
import sys
import numpy as np

def show_trace(time, values, title="Electrophysiological Trace Viewer"):
    """
    Displays a trace in a PyQtGraph window.

    Parameters:
        time (array-like): Time vector (x-axis).
        values (array-like): Signal values (y-axis).
        title (str): Window title.
    """
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)

    win = pg.GraphicsLayoutWidget(title=title)
    win.resize(800, 600)
    win.show()

    plot = win.addPlot(title="Trace Viewer")
    plot.plot(time, values, pen='y')

    plot.setMouseEnabled(x=True, y=True)
    plot.showGrid(x=True, y=True)

    plot.getViewBox().setMouseMode(pg.ViewBox.RectMode)
    plot.getViewBox().setLimits(xMin=min(time), xMax=max(time))

    app.exec_()


def show_two_traces(time, trace1, trace2, title="Two-Trace Viewer"):
    """
    Displays two traces in vertically stacked plots with synchronized time axis.

    Parameters:
        time (array-like): Shared time vector.
        trace1 (array-like): First trace (top).
        trace2 (array-like): Second trace (bottom).
        title (str): Window title.
    """
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)

    win = pg.GraphicsLayoutWidget(title=title)
    win.resize(1000, 600)
    win.show()

    # Ensure all arrays have the same length
    min_length = min(len(time), len(trace1), len(trace2))
    time = np.array(time)[:min_length]
    trace1 = np.array(trace1)[:min_length]
    trace2 = np.array(trace2)[:min_length]

    # First plot (top)
    p1 = win.addPlot(row=0, col=0, title="Trace 1")
    p1.plot(time, trace1, pen='c')
    p1.setMouseEnabled(x=True, y=True)
    p1.showGrid(x=True, y=True)

    # Second plot (bottom)
    p2 = win.addPlot(row=1, col=0, title="Trace 2")
    p2.plot(time, trace2, pen='m')
    p2.setMouseEnabled(x=True, y=True)
    p2.showGrid(x=True, y=True)

    # Link x-axes for synchronized zoom/pan
    p2.setXLink(p1)

    # Optional: set initial visible range
    for p in (p1, p2):
        p.getViewBox().setLimits(xMin=min(time), xMax=max(time))
        p.getViewBox().setMouseMode(pg.ViewBox.RectMode)  # rectangular zoom

    app.exec_()

def show_three_traces(time, trace1, trace2, trace3, title="Three-Trace Viewer"):
    """
        Displays three traces in vertically stacked plots with synchronized time axis.

        Parameters:
            time (array-like): Shared time vector.
            trace1 (array-like): First trace (top).
            trace2 (array-like): Second trace (middle).
            trace3 (array-like): Third trace (bottom).
            title (str): Window title.
    """
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)

    win = pg.GraphicsLayoutWidget(title=title)
    win.resize(1000, 800)
    win.show()

    # Ensure all arrays have the same length by trimming to the shortest
    min_length = min(len(time), len(trace1), len(trace2), len(trace3))
    time = np.array(time)[:min_length]
    trace1 = np.array(trace1)[:min_length]
    trace2 = np.array(trace2)[:min_length]
    trace3 = np.array(trace3)[:min_length]

    # First plot (top)
    p1 = win.addPlot(row=0, col=0, title="Trace 1")
    p1.plot(time, trace1, pen='c')
    p1.setMouseEnabled(x=True, y=True)
    p1.showGrid(x=True, y=True)

    # Second plot (middle)
    p2 = win.addPlot(row=1, col=0, title="Trace 2")
    p2.plot(time, trace2, pen='m')
    p2.setMouseEnabled(x=True, y=True)
    p2.showGrid(x=True, y=True)

    # Third plot (bottom)
    p3 = win.addPlot(row=2, col=0, title="Trace 3")
    p3.plot(time, trace3, pen='y')
    p3.setMouseEnabled(x=True, y=True)
    p3.showGrid(x=True, y=True)

    # Link x-axes for synchronized zoom/pan
    p2.setXLink(p1)
    p3.setXLink(p1)

    # Optional: set initial visible range
    for p in (p1, p2, p3):
        p.getViewBox().setLimits(xMin=min(time), xMax=max(time))
        p.getViewBox().setMouseMode(pg.ViewBox.RectMode)  # rectangular zoom

    app.exec_()