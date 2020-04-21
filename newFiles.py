# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 23:04:41 2016

@author: tanakayumiko
"""

import pandas as pd

def newFiles(newData, pic, Roi_Data_file, status_file, time_file, n):

    data = pd.DataFrame([sorted(range(1,n+1)*4),[1,2,3,4]*n,newData]).T
    data.columns = ['Address', 'ROI', '1']
    data.to_csv(Roi_Data_file, index = False)

    status = pd.DataFrame([sorted(range(1,n+1)*4),[1,2,3,4]*n,[0,0,0,2]*n]).T
    status.columns = ['Address', 'ROI', 'Status']
    status.to_csv(status_file, index = False)
    
    time = pd.DataFrame([range(1,n+1),[i.metadata['t_ms'] for i in pic]]).T
    time.columns = ['Address', '1'] 
    time.to_csv(time_file, index = False)
    
    return data, status, time