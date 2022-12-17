import os
import numpy as np
import matplotlib.pyplot as plt
from gprMax.gprMax import api
from tools.plot_Bscan import mpl_plot
from tools.plot_Ascan import mpl_plot as mpl_Ascan_plot

from tools.outputfiles_merge import get_output_data, merge_files

base_filename = 'user_models\cylinder_Bscan_Peplinski'
# base_filename = 'user_models\cylinder_Bscan_2D'
inp_filename = base_filename+'.in'




def gprmax_processing():
    return 0


api(inp_filename, n=2, geometry_only=False, gpu = {0})
# api(inp_filename, n=40, geometry_only=False)
# api(inp_filename, geometry_only=True)

merge_files(base_filename, removefiles=False)

filename_out_merge = base_filename+'_merged.out'
rxnumber = 1
rxcomponent = 'Ez'
output_data, dt = get_output_data(filename_out_merge, rxnumber, rxcomponent,)
np.savetxt('123.csv', output_data, delimiter=';', )

from gprMax.receivers import Rx
outputs = Rx.defaultoutputs
outputs = ['Ez']

print(outputs)
plt = mpl_plot(filename_out_merge, output_data, dt, rxnumber, rxcomponent)
plt.show()

print(dt)