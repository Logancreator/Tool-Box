# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 17:07:54 2021

@author: Jeff
"""


import os
import numpy as np
import pandas as pd
  
motif_dir="C:/Users/Jeff/OneDrive/桌面/quail/next/cluster_1/test"
  
for i in os.listdir(motif_dir):
    df=pd.read_table(os.path.join(motif_dir, i),header=None)
    ATGC=["A","C","G","T"]
    df1=df.T
    symbol=["|"]*4
    df1.insert(0, 'a', ATGC)
    df1.insert(1, 'b', symbol)
    print(df1)
    df1.to_csv(os.path.join(motif_dir, i),sep="\t",index=False,header=False)

        