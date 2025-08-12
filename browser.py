import os, sys
import pyqtgraph as pg
import numpy as np
import heka_reader
from pyqtgraph.Qt import QtWidgets, QtCore
import pandas as pd
# Import the paths from main
from config import ROOT_FOLDER, IMPORT_FOLDER, METADATA_FILE
from config import analysis_points
import json


app = pg.mkQApp()

# Configure Qt GUI:

# Main window + splitters to let user resize panes
win = QtWidgets.QWidget()
layout = QtWidgets.QGridLayout()
win.setLayout(layout)
hsplit = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
layout.addWidget(hsplit, 0, 0)
vsplit = QtWidgets.QSplitter(QtCore.Qt.Vertical)
hsplit.addWidget(vsplit)
w1 = QtWidgets.QWidget()
w1l = QtWidgets.QGridLayout()
w1.setLayout(w1l)
vsplit.addWidget(w1)

# Button for loading .dat file
load_btn = QtWidgets.QPushButton("Load...")
w1l.addWidget(load_btn, 0, 0)

# Tree for displaying .pul structure
tree = QtWidgets.QTreeWidget()
tree.setHeaderLabels(['Node', 'Label'])
tree.setColumnWidth(0, 200)
w1l.addWidget(tree, 1, 0)

# Tree for displaying metadata for selected node
data_tree = pg.DataTreeWidget()
vsplit.addWidget(data_tree)

# Plot for displaying trace data
plot = pg.PlotWidget()
hsplit.addWidget(plot)

# Resize and show window
hsplit.setStretchFactor(0, 400)
hsplit.setStretchFactor(1, 600)
win.resize(1200, 800)
win.show()


def load_clicked():
    """Display a popup menu with available .dat files from metadata."""
    try:
        # Read metadata file
        metadata_df = pd.read_excel(METADATA_FILE)
        file_names = metadata_df['file_name'].tolist()

        # Create popup menu
        menu = QtWidgets.QMenu()
        for fname in file_names:
            action = menu.addAction(fname)

        # Show menu at button position
        selected_action = menu.exec_(load_btn.mapToGlobal(QtCore.QPoint(0, load_btn.height())))

        if selected_action:
            selected_file = selected_action.text()
            dat_path = os.path.join(IMPORT_FOLDER, selected_file)
            if os.path.exists(dat_path):
                load(dat_path)
            else:
                QtWidgets.QMessageBox.warning(
                    win,
                    "File Not Found",
                    f"The file {selected_file} was not found in the import folder."
                )
    except Exception as e:
        QtWidgets.QMessageBox.critical(
            win,
            "Error",
            f"Error loading metadata: {str(e)}"
        )


load_btn.clicked.connect(load_clicked)


def load(file_name):
    """Load a new .dat file into the browser."""
    global bundle, tree_items
    bundle = heka_reader.Bundle(file_name)

    # Clear and update tree
    tree.clear()
    update_tree(tree.invisibleRootItem(), [])
    replot()


def update_tree(root_item, index):
    """Recursively build tree structure."""
    global bundle
    root = bundle.pul
    node = root
    for i in index:
        node = node[i]
    node_type = node.__class__.__name__
    if node_type.endswith('Record'):
        node_type = node_type[:-6]
    try:
        node_type += str(getattr(node, node_type + 'Count'))
    except AttributeError:
        pass
    try:
        node_label = node.Label
    except AttributeError:
        node_label = ''
    item = QtWidgets.QTreeWidgetItem([node_type, node_label])
    root_item.addChild(item)
    item.node = node
    item.index = index
    if len(index) < 2:
        item.setExpanded(True)
    for i in range(len(node.children)):
        update_tree(item, index + [i])


def replot():
    """Update plot and data tree when user selects a trace."""
    plot.clear()
    data_tree.clear()

    selected = tree.selectedItems()
    if len(selected) < 1:
        return

    sel = selected[0]
    fields = sel.node.get_fields()
    data_tree.setData(fields)

    for sel in selected:
        index = sel.index
        if len(index) < 4:
            return

        # Extract indices from the selection
        group_id = index[0]
        series_id = index[1]
        sweep_id = index[2]
        trace_id = index[3]

        trace = sel.node
        plot.setLabel('bottom', trace.XUnit)
        plot.setLabel('left', trace.Label, units=trace.YUnit)
        data = bundle.data[index]
        time = np.linspace(trace.XStart, trace.XStart + trace.XInterval * (len(data) - 1), len(data))
        # Plot the main trace
        plot.plot(time, data)

        # Get the file name from the bundle
        file_name = os.path.basename(bundle.file_name)

        # Check if we have analysis points for this file and indices
        #print(f"Looking for points in file: {file_name}")
        #print(f"Existing analysis_points: {json.dumps(analysis_points, indent=2)}")
        #print(f"Points to plot: {points}")

        if (file_name in analysis_points and group_id in analysis_points[file_name] and series_id in analysis_points[file_name][group_id] and sweep_id in analysis_points[file_name][group_id][series_id] and trace_id in analysis_points[file_name][group_id][series_id][sweep_id]):
            points = analysis_points[file_name][group_id][series_id][sweep_id][trace_id]
            # Plot each type of point with different symbols and colors
            symbols = {
                'threshold': ('o', 'red', 'Threshold'),
                'half_duration_start': ('s', 'blue', 'Half Duration Start'),
                'half_duration_end': ('s', 'green', 'Half Duration End'),
                'peak': ('t', 'yellow', 'Peak'),
                'ahp': ('d', 'purple', 'AHP'),
                'dvdt_max': ('p', 'cyan', 'Max dV/dt')
            }

            for point_type, (symbol, color, label) in symbols.items():
                if points[point_type]:  # If we have points of this type
                    t_points, v_points = zip(*points[point_type])
                    scatter = pg.ScatterPlotItem(
                        x=t_points,
                        y=v_points,
                        symbol=symbol,
                        size=10,
                        pen=pg.mkPen(None),
                        brush=pg.mkBrush(color),
                        name=label
                    )
                    plot.addItem(scatter)

    # Add legend
    plot.addLegend()


tree.itemSelectionChanged.connect(replot)

# Optional: Load demo data if present
demo = 'DemoV9Bundle.dat'
if os.path.isfile(demo):
    load(demo)

if __name__ == '__main__':
    if sys.flags.interactive == 0:
        app.exec_()
