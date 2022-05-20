# -*- coding: utf-8 -*-
"""
Created on Sun Aug 29 14:45:56 2021

@author: Jeff
"""



import pandas as pd

import numpy  as np

import re
# df=pd.read_table("C:\\Users\\Jeff\\OneDrive\\桌面\\quail\\peak\\idr_peakAnno.txt")


# #print(df)

# def get_gene_list(DEGfile):

#     gene_table=pd.read_csv(DEGfile,sep="\t")
#     #print(gene_table["Unnamed: 0"])
#     return gene_table["Unnamed: 0"].tolist()
    
    
    
# file1="C:\\Users\\Jeff\\OneDrive\\桌面\\quail\\RNAseq\\control_treat.DESeq2.up.txt"
# file2="C:\\Users\\Jeff\\OneDrive\\桌面\\quail\\RNAseq\\control_treat.DESeq2.down.txt"
# DEGs_up=get_gene_list(file1)
# print(len(DEGs_up))
# DEGs_down=get_gene_list(file2)
# print(len(DEGs_down))

# # def get_gene_length(gtf):
# #     gtf_file=pd.read_table(gtf,sep="\t",header= None)
    
# #     #if gtf_file.iloc[:,2] in DEGs_up:
# #     #    print(gtf_file)
# #     print(gtf_file[gtf_file.iloc[:,2].str.contains('gene')])



# #data = pd.read_excel("C:\\Users\\Jeff\\OneDrive\\桌面\\quail\\genomic.gtf.txt",header=None)
# #pd.read_table("C:\\Users\\Jeff\\OneDrive\\桌面\\quail\\genomic.gtf")


# #输入需要查找的文件
# gtf='C:\\Users\\Jeff\\OneDrive\\桌面\\quail\\genomic.gtf'
# def get_gene_length(gtf,DEG):
#     with open(gtf,'r') as gtf:
#         list=[]
    
#         for line in gtf:
#         #分成9列
#             line1=line.strip().split('\t',8)
#             #print(line1)
#     #解决IndexError: list index out of range问题
#             if line1[2]=="gene":
#                 if line1[8].split("\"")[3] in DEG:
#                      list.append(line1[0:9])
#     return list

# gene_length_up=get_gene_length(gtf,DEGs_up)
# print(len(gene_length_up))
# gene_length_dwon=get_gene_length(gtf,DEGs_down)
# print(len(gene_length_dwon))


'''
with open("C:/Users/Jeff/OneDrive/桌面/quail/peak/expand_down_list.txt") as f:
    for line in f.readlines():
        temp=line.split()
        #print(temp[0:3])
        with open("C:/Users/Jeff/OneDrive/桌面/quail/peak/idr_peakAnno.txt") as g:
            for peak in g.readlines():
                temp_peak=peak.split()
                #print(temp_peak)
                if temp_peak[0]==temp[0]:
                    if temp_peak[1]>=temp[1] and temp_peak[2]<=temp[2]:
                        #print(temp_peak)
                        with open("C:/Users/Jeff/OneDrive/桌面/quail/peak/expand_down_gene_+-50kb.txt",'a+') as h:
                            for line in temp_peak:
                                #print(temp_peak)
                                h.write(line+"\t")
                            h.write("\n")
                            
'''                            
with open("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose\\鹅基因组\\GeneName.txt") as g:
    for name in g.readlines():
        gene=name.split()
        #print(str(gene[0]))
        if len(gene)<=3:
            gene.append(str(gene[0])[0:5])
            #print(gene)


    with open("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose\\genes.gtf") as f:
        for line in f.readlines():
            print(line)
            new = line.replace("\""+str(gene[0])+"\"","\""+str(gene[2])+"\"")
            print(new)
