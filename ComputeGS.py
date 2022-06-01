#coding=utf-8
#python ComputeGS.py your_fasta
##Encode by Jianye   ##Use to get genome size for MACS2 
import sys
import os

def File_reader(path):
    with open(path,'r') as f:
        for line in f:
            line = line.strip()
            line = line.upper()
            if not line.startswith(">"):
                baseA = line.count("A")
                baseT = line.count("T")
                baseC = line.count("C")
                baseG = line.count("G")
                aList.extend([baseA, baseT, baseC, baseG])
    return aList,sum(aList)
def GC_content(list):
    baseC_count=baseG_count=[]

    base_C_num=2
    while base_C_num < len(list)-4:
        base_C_num+=4
        baseC_count.append(list[base_C_num])
    base_G_num=3
    while base_G_num < len(list)-4:
        base_G_num+=4
        baseG_count.append(list[base_G_num])


    return (sum(baseG_count)+sum(baseC_count))/(sum(list))

if __name__ == '__main__' :
    aList=[]

    fa_file = sys.argv[1]

    aList,res=File_reader(fa_file)
    
    print("Effective Genome Size :",res )

    GC=GC_content(aList)

    print("Genome GC Content :",round(GC*100,2),"%" )
