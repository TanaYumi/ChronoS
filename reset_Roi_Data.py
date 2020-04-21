# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 20:41:46 2019

@author: gomak
"""

import pandas as pd

filepath = 'C:/Users/gomak/OneDrive - The University of Tokyo/Documents/sampledata_/ROI_Data.csv'
input_data = pd.read_csv(filepath)
return_data =  input_data.loc[:,:'408']
return_data.to_csv(filepath, index=None)

filepath = 'C:/Users/gomak/OneDrive - The University of Tokyo/Documents/sampledata_/time.csv'
input_data = pd.read_csv(filepath)
return_data =  input_data.loc[:,:'408']
return_data.to_csv(filepath, index=None)