# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 15:21:09 2016

@author: tanakayumiko
"""

import pandas as pd
from scipy.stats import linregress


def readFile(Roi_Data_file, status_file, time_file, rawdata, rawtime):
    status = pd.read_csv(status_file)
    rawdata = pd.pivot_table(rawdata, index = ['Address', 'ROI'])
    names = []
    for i in range(1,len(rawdata.columns)+1):
        names.append('%d'%i)
    rawdata = rawdata.reindex(columns = names)
    rawdata = rawdata.sub(rawdata['1'], axis = 0)
    data = rawdata.copy()
    
    if len(data)%4 != 0:
        raise ValueError('data length is not multiple of 4')
    for i in range(1,len(data)//4):
        back = int(status[status['Address']==i][status['Status']==2].ROI.values)
        for j in range(1,5):
            data.ix[i,j] = rawdata.ix[i,j] - rawdata.ix[i,back]
            
    time = pd.pivot_table(rawtime, index = ['Address'])
    status = pd.pivot_table(status, index = ['Address', 'ROI'])

    return data, status, time
    
def back_fit_suggest(scan, data, start, status):
    sigma = 15
    slope = []
    intercept = []
    stdev = []
    suggestion = []
    for i in range(1,1+data.shape[0]//4):
        for j in range(1,5):
            s, inter, r, p, stderr = \
            linregress(range(start+1,max(scan-2,start+2)), data.ix[i,j].values[start-1:max(scan-4,start)])
            slope.append(s)
            intercept.append(inter)
            stdev.append(stderr)
            line = data.ix[i,j].values
            if scan < 4:
                point1, point2, point3 = 0,0,0
            else:
                point1 = line[-1]>s*(scan)+inter+stderr*sigma
                point2 = line[-2]>s*(scan-1)+inter+stderr*sigma
                point3 = line[-3]>s*(scan-2)+inter+stderr*sigma
                
            if point1 and point2+point3:
                suggestion.append([i,j])
                
    status['slope']=slope
    status['intercept']=intercept
    status['stderr']=stdev
    return suggestion,status
            
