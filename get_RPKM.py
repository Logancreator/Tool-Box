# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 08:48:05 2021

@author: Jeff
"""

import pandas as pd
import numpy as np


data=pd.read_csv("C:/Users/Jeff/OneDrive/桌面/quail/RNAseq/featurecounts.result",sep="\t")

#print(data)


def get_RPKM():
    counts_sum_all={}
    sample_rpkm={}
    for i in data.keys()[2:]:
        counts_sum_all["%s" %i]=data[i].sum()
        sample_rpkm["%s" %i]=pd.Series((data[i]/(counts_sum_all[i]/1000000)))/(data["Length"]/1000)
    
    
    index=pd.Series(data["Geneid"])
    data_rpkm=pd.concat([index,sample_rpkm["control_1"],sample_rpkm["control_2"],sample_rpkm["ld7_1"],sample_rpkm["ld7_2"]],axis=1)
    
    data_rpkm.rename({0:'control_1',1:'control_2',2:'ld7_1',3:'ld7_2'},axis='columns',inplace=True)
    data_rpkm.loc['rpkm_sum']=data_rpkm.apply(lambda x:x.sum())
    return data_rpkm

data1=get_RPKM()
print(type(data1))

data1.to_csv("C:/Users/Jeff/OneDrive/桌面/quail/RNAseq/gene_rpkm.csv",sep="\t")