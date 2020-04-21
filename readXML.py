# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 12:52:42 2016

@author: tanakayumiko
"""

import xml.etree.ElementTree as ET


def readXML(position_file, n):
    tree = ET.parse(position_file)
    
    root=tree.getroot()
    
    tag = []
    for i in range(2,n):
        if root[0][i][0].attrib['value']=='true':
            tag.append(root[0][i].tag)
            
    return tag
        
position_file = '/Volumes/SINGLECELL/sampledata/160418.xml'