

import os
import sys

code_path = '/Users/mike/Documents/Research/Scripts/Python/Greenland Model Analysis/Fjord/Kangerlussuaq'

sys.path.insert(1,os.path.join(code_path,'Figures','Glacier'))
sys.path.insert(1,os.path.join(code_path,'Figures','Ocean'))
sys.path.insert(1,os.path.join(code_path,'Data','Ocean'))

# import plot_glacier_melt_rate as mr
# mr.plot_melt_rate_comparison()

# import collect_melange_melt_timeseries_to_nc as co
# co.collect_melange_melt()
# import plot_model_melange_melt_comparison as mc
# mc.plot_melange_comparison()

import plot_glacier_transect_mean


