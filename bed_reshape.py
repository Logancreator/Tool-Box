import pandas as pd
import os
import argparse


def read_bed(file):

    df=pd.read_csv(file,sep="\t",header=None)

    df.columns=['chr','start','end','strand','fold']  

    return df

def read_fai(file):

    df=pd.read_csv(file,sep="\t",header=None)

    df.columns=['chr','length','other1','other2','other3']

    return df

def main():
    parser = argparse.ArgumentParser(description='bed reshape')
    parser.add_argument('--bed', type=str, required=True, help='bed file needed input')
    parser.add_argument('--fai', type=str, required=True,  help='fai file needed input')
    parser.add_argument('--out', type=str, required=True,  help='bed file output')
    args = parser.parse_args()

    bed_path=args.bed

    bed_opth=args.out

    fai_path=args.fai

    bed=read_bed(bed_path)

    fai=read_fai(fai_path)

    for row in bed.itertuples():
        if row.end >= fai[fai['chr'] == row.chr]["length"].iloc[0]:
            bed.loc[row.Index,"end"]=fai[fai['chr'] == row.chr]["length"].iloc[0]

        if row.start == 0 :
            bed.loc[row.Index,"start"]=1
    
    bed.to_csv(bed_opth,index=0,header=0,sep="\t")
    print("\n" * 5)

if __name__ == '__main__':
    try:
        main()
        
        print("bed reshape sucessfully!")

    except (ValueError, ArithmeticError):
        print("\n" * 5)
        print("A number format exception occurred in the program!")
    except :
        print("\n" * 5) 
        print("unknown abnormal")