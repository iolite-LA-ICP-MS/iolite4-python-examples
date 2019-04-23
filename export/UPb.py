import numpy as np
import pandas as pd
from functools import partial
from xlwt import Workbook, easyxf

wb = Workbook()
ws = wb.add_sheet('Table', cell_overwrite_ok=True)

sels = [s for sg in data.selectionGroupList() for s in sg.selections()]

class ChannelNames:
	Pb206_cps = 'Pb206_CPS'
	Uppm = 'U ppm'
	Th_U = 'Th/U'
	Pb206_Pb204 = 'Final Pb206/Pb204'
	Pb206_U238 = 'Final Pb206/U238'
	Pb206_U238_age = 'Final Pb206/U238 age'
	Pb207_U235 = 'Final Pb207/U235'
	Pb207_U235_age = 'Final Pb207/U235 age'
	Pb208_Pb206 = 'Final Pb208/Pb206'
	U238_Pb206 = 'Final U238/Pb206'
	Pb207_Pb206 = 'Final Pb207/Pb206'
	Pb207_Pb206_age = 'Final Pb207/Pb206 age'
	Wetherill_rho = 'rho 206Pb/238U v 207Pb/235U'
	TW_rho = 'rho 207Pb/206Pb v 238U/206Pb'
	Pb208_Th232 = 'Final Pb208/Th232'
	Pb208_Th232_age = 'Final Pb208/Th232 age'

def channel_data(channel_name, selection):
	try:
		result = data.result(selection, data.timeSeries(channel_name))	
		u1s_pct = 0.5*100*result.uncertainty()/result.value()
		u2s_abs = result.uncertainty()
		pu2s_abs = result.propagatedUncertainty()
		return (result.value(), u1s_pct, u2s_abs, pu2s_abs)
	except ZeroDivisionError:
		return ('#N/A', '#N/A', '#N/A', '#N/A')

def selection_data(property_name, selection):
	return (selection.property(property_name),)

def associated_data(result_name, selection):
	return (data.associatedResult(selection, result_name).value(),)

def conc_pct(selection):
	try:
		age6_38 = data.result(selection, data.timeSeries(ChannelNames.Pb206_U238_age)).value()
		age7_6 = data.result(selection, data.timeSeries(ChannelNames.Pb207_Pb206_age)).value()	
		return (100*age6_38/age7_6,)
	except ZeroDivisionError:
		return ('#N/A', '#N/A', '#N/A', '#N/A')


def write_header():
	# write the PlasmaAge template header
	ws.write(0, 0, 'Header')

def write_borders():
	# add borders like the template
	ws.write(0, 0, None, style=easyxf('border: bottom thick'))

def write_footer(index, message):
	# write a footer
	ws.write(10, 0, 'Footer')

def write_column(column_index, data_func=None, uncert_types=[0], default_content='#N/A'):
	# iterate through selections and add column
	row_index = 3 # start at row 3 to leave room for header
	col_index_start = column_index

	for s in sels:
		column_index = col_index_start

		if data_func:
			ws.write(row_index, column_index, data_func(s)[0])
			for i in uncert_types:
				column_index += 1
				ws.write(row_index, column_index, data_func(s)[i])		
		else:
			ws.write(row_index, column_index, default_content)
	
		row_index += 1
	
write_header()
write_column(0, data_func = partial(selection_data, 'Name'))
write_column(1, data_func = partial(selection_data, 'Comment'))
write_column(2) # f206c
write_column(3, data_func = partial(channel_data, ChannelNames.Pb206_cps))
#write_column(4, data_func = partial(channel_data, ChannelNames.Uppm))
#write_column(5, data_func = partial(channel_data, ChannelNames.Th_U))
#write_column(6, data_func = partial(channel_data, ChannelNames.Pb206_Pb204), uncert_types=[1])
write_column(8, data_func = partial(channel_data, ChannelNames.U238_Pb206), uncert_types=[1])
write_column(10, data_func = partial(channel_data, ChannelNames.Pb207_Pb206), uncert_types=[1])
write_column(19, data_func = partial(associated_data, ChannelNames.TW_rho))
#write_column(13, data_func = partial(channel_data, ChannelNames.Pb208_Pb206), uncert_types=[1])
write_column(15, data_func = partial(channel_data, ChannelNames.Pb207_U235), uncert_types=[1])
write_column(17, data_func = partial(channel_data, ChannelNames.Pb206_U238), uncert_types=[1])
write_column(19, data_func = partial(associated_data, ChannelNames.Wetherill_rho))
write_column(20, data_func = partial(channel_data, ChannelNames.Pb208_Th232), uncert_types=[1])
write_column(22, data_func = partial(channel_data, ChannelNames.Pb207_Pb206_age), uncert_types=[2,3])
write_column(25, data_func = partial(channel_data, ChannelNames.Pb206_U238_age), uncert_types=[2,3])
write_column(28, data_func = partial(channel_data, ChannelNames.Pb207_U235_age), uncert_types=[2,3])
write_column(31, data_func = partial(channel_data, ChannelNames.Pb208_Th232_age), uncert_types=[2,3])
write_column(34, data_func = conc_pct)
write_borders()
write_footer(1, 'asdf')

wb.save(export_filepath)
