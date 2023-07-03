#!/usr/bin/env python3

import pandas as pd
import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib.transforms
import scipy.stats

font = {'weight' : 'normal',
        'size'   : 17}
matplotlib.rc('font', **font)
in_file = 'C:/Users/df391/OneDrive - University of Exeter/Post_Doc_Covex_Seascape/Shutler_Cross_Shelf_Transport/Results_500m_data.xlsx'
plot_location = 'plots/data_means/'
data = pd.read_excel(in_file)
print(data.columns)
