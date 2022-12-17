import os
import numpy as np
import matplotlib.pyplot as plt
from gprMax.gprMax import api
from tools.outputfiles_merge import get_output_data, merge_files

# dmax = r".\GprmaxCode"
filename = '../user_models/cylinder_Bscan_GSSI_1500.in'
# filename = 'user_models\cylinder_Bscan_2D.in'


api(filename, n=1, geometry_only=False, gpu = {0})
# api(filename, n=10, geometry_only=False)
# # api(filename, n=100, geometry_only=False)
# merge_files(filename, removefiles=True)
# filename = 'user_models\cylinder_Ascan_2D.out'
#
filename_out_merge = filename[:-2]+'in_merged.out'
filename_out = filename[:-2]
rxnumber = 1
rxcomponent = 'Ez'
output_data, dt = get_output_data(filename_out_merge, rxnumber, rxcomponent,)
np.savetxt('123.txt', output_data, delimiter=' ')

from tools.plot_Ascan import mpl_plot
from tools.plot_Bscan import mpl_plot

from gprMax.receivers import Rx
outputs = Rx.defaultoutputs
outputs = ['Ez']

print(outputs)
plt = mpl_plot(filename_out_merge, outputs)
plt.show()