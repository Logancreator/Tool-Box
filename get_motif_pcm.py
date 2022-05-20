# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 17:07:54 2021

@author: Jeff
"""


import os
import numpy as np
import pandas as pd



# def get_motif_path(motif_dir):
#     files=[]
#     for i in os.listdir(motif_dir):
#         files.append(os.path.join(motif_dir, i))
#     return files
        
# def rename_motif_file(motif_dir):
#     motif_path=get_motif_path(motif_dir)
#     return motif_path
    
# temp=rename_motif_file(motif_dir)
# file_name=[]
# for i in range(0,len(temp),1):
    
#     with open(temp[i],'r') as f:
        
#         for j in f.readlines():
#             if j.startswith(">"):
#                 line=j.split()
                
#                 file_name.append(line[1])
                
#         f.close()
#     print(len(file_name))
#     print(len(temp))
    
#     os.rename(temp[i], str(motif_dir)+"/"+str(file_name[i]).replace("/","_").replace(":","_")+".pcm")
    
    
    
motif_dir="C:/Users/Jeff/OneDrive/桌面/quail/next/cluster_1/test"



# files=[]   
# for i in os.listdir(motif_dir):
#     files.append(os.path.join(motif_dir, i))
#     # with open(os.path.join(motif_dir, i)) as f:
#     #     for j in f.readlines():
#     #         if j.startswith(">")==False:
#     #             print(j)
                
                
#     with open(os.path.join(motif_dir, i),"r+") as f:
#         d = f.readlines()
#         f.seek(0)
#         for i in d:
#             if not i.startswith(">"):
#                 f.write(i)
#         f.truncate()                
for i in os.listdir(motif_dir):
    df=pd.read_table(os.path.join(motif_dir, i),header=None)
    ATGC=["A","C","G","T"]
    df1=df.T
    symbol=["|"]*4
    df1.insert(0, 'a', ATGC)
    df1.insert(1, 'b', symbol)
    print(df1)
    df1.to_csv(os.path.join(motif_dir, i),sep="\t",index=False,header=False)
    


        