#/ Type: UI
#/ Name: Channel Calculator
#/ Authors: Joe Petrus and Bence Paul
#/ Description: This is an example that lets you create a channel from existing channels
#/ References: None
#/ Version: 1.0
#/ Contact: support@iolite-software.com

from iolite.QtGui import QAction, QLabel, QLineEdit, QCompleter, QDialog, QHBoxLayout, QVBoxLayout, QDialogButtonBox, QComboBox
from iolite.QtCore import Qt
import numpy as np
import re

class AutoCompleteEdit(QLineEdit):
	"""
	A slightly modified version of:
	http://blog.elentok.com/2011/08/autocomplete-textbox-for-multiple.html
	"""
	def __init__(self, model, separator = ' ', addSpaceAfterCompleting = True, parent = None):
		super().__init__(parent)
		self._separator = separator
		self._addSpaceAfterCompleting = addSpaceAfterCompleting
		self._completer = QCompleter(model)
		self._completer.setWidget(self)
		self._completer.activated.connect(self._insertCompletion)
		self._keysToIgnore = [Qt.Key_Enter,
							  Qt.Key_Return,
							  Qt.Key_Escape,
							  Qt.Key_Tab]

	def _insertCompletion(self, completion):
		extra = len(completion) - len(self._completer.completionPrefix)
		extra_text = completion[-extra:]
		if self._addSpaceAfterCompleting:
			extra_text += ' '
		self.setText(self.text + extra_text)

	def textUnderCursor(self):
		text = self.text
		textUnderCursor = ''
		i = self.cursorPosition - 1
		while i >=0 and text[i] != self._separator:
			textUnderCursor = text[i] + textUnderCursor
			i -= 1
		return textUnderCursor

	def keyPressEvent(self, event):
		if self._completer.popup().isVisible():
			if event.key() in self._keysToIgnore:
				event.ignore()
				return
		QLineEdit.keyPressEvent(self, event)
		completionPrefix = self.textUnderCursor()
		et = event.text()
		if completionPrefix != self._completer.completionPrefix:
			self._updateCompleterPopupItems(completionPrefix)
		if len(et) > 0 and len(completionPrefix) > 0:
			self._completer.complete()
		if len(completionPrefix) == 0:
			self._completer.popup().hide()

	def _updateCompleterPopupItems(self, completionPrefix):
		self._completer.setCompletionPrefix(completionPrefix)
		self._completer.popup().setCurrentIndex(self._completer.completionModel().index(0,0))

def createUIElements():
	action = QAction("Channel Calculator", ui)
	action.triggered.connect(calculate)
	ui.setAction(action)
	ui.setMenuName(['Examples'])
		

def calculate():

	d = QDialog()
	d.setWindowTitle('Channel Calculator')

	hl = QHBoxLayout()

	name_le = QLineEdit(d)
	hl.addWidget(name_le)

	hl.addWidget(QLabel('='))

	channel_names = [re.sub('\W|^(?=\d)','_', channel_name) for channel_name in data.timeSeriesNames()]
	eqn_le = AutoCompleteEdit(channel_names, parent=d)
	hl.addWidget(eqn_le)
	
	type_cb = QComboBox(d)
	channel_types = {'Input': data.Input, 'Intermediate': data.Intermediate, 'Output': data.Output}
	for tname, tv in channel_types.items():
		type_cb.addItem(tname, tv)
	hl.addWidget(type_cb)

	bb = QDialogButtonBox(QDialogButtonBox.Cancel | QDialogButtonBox.Ok, Qt.Horizontal, d)
	bb.accepted.connect(d.accept)
	bb.rejected.connect(d.reject)
	hl.addWidget(bb)

	d.setLayout(hl)

	if d.exec() == QDialog.Accepted:
		for channel_name in data.timeSeriesNames():
			exec('%s = data.timeSeries("%s").data()'%(re.sub('\W|^(?=\d)','_', channel_name), channel_name))

		new_data = eval('%s'%eqn_le.text)
		new_type = type_cb.currentData
		new_channel = data.createTimeSeries(name_le.text, new_type, None, new_data)
		new_channel.setProperty('Created by', 'Channel Calculator')
		data.dataChanged.emit()
