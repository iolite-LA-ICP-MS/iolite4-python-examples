# A python-based UI module for iolite 4 starts with some metadata
#/ Type: UI
#/ Name: UI Introduction
#/ Authors: Joe Petrus and Bence Paul
#/ Description: This is an example python-based UI plugin for iolite 4
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

"""
Qt imports can be done through 'iolite', e.g.
from iolite.QtGui import QLabel

Before functions are called a few additional objects are
added to the module:

data	an interface to iolite's C++ data. E.g. you can get
        existing time series data or selection groups with it
        as well as make new ones.

IoLog	an interface to iolite's logging facility. You can add
        messages with, e.g., IoLog.debug('My message')

ui		an interface to the PythonUI C++ class in iolite from 
		which an action associated with this plugin can be 
		added to the iolite interface either as a menu or 
		sidebar item
"""

from iolite.QtGui import QAction, QMessageBox
import numpy as np

def installUIHooks(window):
	"""
	This is the only method that must be defined in a python UI plugin.
	It will be automatically called to make the plugin available in the
	user interface.

	To make the plugin invokable from the user interface, you can
	use one of the following:
	
	ui.appendSidebarAction(action)
	ui.appendActionToMenu(menuStringList, action)
	ui.appendMenuToMenu(menuStringList, action)
	
	where the action is a QAction whose 'triggered' signal is connected
	to something useful, e.g.,
	action = QAction("Action name", window)
	action.triggered.connect(slot)		
	"""

	"""
	As an example, we'll add a menu item that when clicked pops up
	a message box summarizing the selections we've made.
	"""
	
	action = QAction("Summarize Selections", window)
	ui.appendActionToMenu(["Tools", "Examples"], action)
	action.triggered.connect(summarize)

def summarize():

	msg = ""

	for sg in data.selectionGroupList():
		durations = np.array([s.duration for s in sg.selections()])/1000

		msg += "<p><b>%s</b>: %i selections<br>"%(sg.name, len(durations))
		msg += "Average duration = %.3f sec.<br>"%(durations.mean())
		msg += "Shortest duration = %.3f sec.<br>"%(durations.min())
		msg += "Longest duration = %.3f sec.</p>"%(durations.max())


	QMessageBox.information(None, "Selection summary", msg)
