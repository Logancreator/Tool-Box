import pandas as pd
import numpy as np
import os

os.chdir("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\")

def read_bed(file):

    df=pd.read_csv(file,sep="\t",header=None)

    df.columns=['chr','start','end','strand','fold']  

    return df

def read_fai(file):

    df=pd.read_csv(file,sep="\t",header=None)

    df.columns=['chr','length','other1','other2','other3']

    return df

if __name__ == '__main__':
    bed_path="diffbind\\Lay_vs_Breed_deseq2_down.bed"

    bed_opth="diffbind\\Lay_vs_Breed_deseq2_down_reshape.bed"

    bed=read_bed(bed_path)

    fai=read_fai("GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.fa.fai")

    for row in bed.itertuples():
        if row.end >= fai[fai['chr'] == row.chr]["length"].iloc[0]:
            bed.loc[row.Index,"end"]=fai[fai['chr'] == row.chr]["length"].iloc[0]

        if row.start == 0 :
            bed.loc[row.Index,"start"]=1
    
    bed.to_csv(bed_opth,index=0,header=0,sep="\t")
