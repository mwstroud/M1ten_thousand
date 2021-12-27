from bmtools.cli.plugins.plot import plot as plt
from bmtools.cli.plugins.util import util
import matplotlib.pyplot as plot
import pandas as pd
import h5py
import numpy as np
import os
import sys
import time

plt.plot_inspikes('config.json')
#plt.plot_basic_cell_info('config.json')
#plt.plot_3d_positions(config='config.json',populations='cortex',group_by='pop_name',title='Cell Positions',save_file=False)#
##plt.plot_3d_positions(config='config.json',populations='shell',group_by='pop_name',title='Cell Positions',save_file=False)
#plot.show()