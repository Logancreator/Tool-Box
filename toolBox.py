import pandas as pd
import os
import json
from pickle import NONE
import pandas as pd

def get_dict():

    data = pd.read_table("C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\0806Oryza_HJX74_top_level_v2.2.pep.fa",header=None,sep="\t")  

    data.columns = ['chr1','chr2',2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    data = data.drop_duplicates(subset='chr1', keep='first')

    dic=dict(zip(data['chr1'].tolist(),data['chr2'].tolist()))

    return dic

def get_csvFile(root):

    # 获取文件夹所有文件
    files_list = os.listdir(root)

    # 1.获取指定后缀（如txt）的文件
    filter_files_list = [fn for fn in files_list if fn.endswith("csv")]

    return filter_files_list

def merge(folder_name,result_file):
    # 2. 获取那个文件夹中所有的文件名字：
    file_names = os.listdir(folder_name)
    #3.创建一个输出表
    writer = pd.ExcelWriter(result_file)

    for file_name in file_names:
        if file_name.endswith("csv"):
            data = pd.read_csv(folder_name+"\\"+file_name,encoding='utf-8')
            data.to_excel(writer,file_name.split(".")[0],index=False)

    #4.保存，并关闭当前文件
    print('数据输出成功')
    writer.save()

def main():
    os.chdir("C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna\\sig")

    dic=get_dict()

    for file in get_csvFile("C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna\\sig"):
        os.chdir("C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna\\sig")
        sig=pd.read_csv(file)

        new_col1=[]
        
        new_col2=[]
        for chr in sig['Row.names'].tolist():

            new=dic.get(chr.replace('model.',''))

            new_col1.append(new)

            new_col2.append(str(new).split("-")[0])
            


        sig['rowname1']=new_col1
        
        sig['rowname2']=new_col2
        #print(sig.columns)
        #print(len(sig.columns))
        head = [0,30,31]
        behind = list(range(1,30,1))
        
        col=head+behind

        sig=sig.iloc[:,col]

        os.chdir("C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna\\sig_rename")

        sig.to_csv(file,sep=",",index=False)

    merge()
def run1():
        # main()
    #def merge():
    # 1. 获取一个要合并的文件夹的名称：
    folder_name = "C:\\Users\\Jeff\\OneDrive\\桌面\\scATAC\\P\\da_peaks_path"

    file_list = os.listdir(folder_name)

    # 2. 获取那个文件夹中所有的文件名字：
    filter_files_list = [fn for fn in file_list if fn.endswith("peaks.csv")]

    #3.创建一个输出表
    writer = pd.ExcelWriter('C:\\Users\\Jeff\\OneDrive\\桌面\\scATAC\\P\\da_peaks_path\\celltype_daPeak_peaks.xlsx')

    for file_name in filter_files_list:
        data = pd.read_csv('C:\\Users\\Jeff\\OneDrive\\桌面\\scATAC\\P\\da_peaks_path\\'+file_name,encoding='utf-8')
        data.to_excel(writer,file_name.split(".")[0],index=False)

    #4.保存，并关闭当前文件
    print('数据输出成功')
    writer.save()

def run2():
    os.chdir("C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna")
    dic = get_dict()

    file = get_csvFile("C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna")
    for file_name in file:
        data = pd.read_csv(file_name,encoding='utf-8')
        print(data)
        new_col1=[]
        
        new_col2=[]

        for chr in data['Unnamed: 0'].tolist():
            new=dic.get(chr.replace('model.',''))

            new_col1.append(new)

            new_col2.append(str(new).split("-")[0])

        data['rowname1']=new_col1
        
        data['rowname2']=new_col2
        print(data)
        data.to_csv("normCount_rename.csv")
def json_analysis(my_dir,outputName):
    def Read_json(file_path):
        with open(file_path) as f:
            return json.load(f)
        
    def Extract_info(col_name,col_info):
        value=[]
        key=['before_total_reads','before_total_bases','before_q20_bases','before_q30_bases','before_q20_rate','before_q30_rate','before_gc_content','before_duplication','after_total_reads','after_total_bases','after_q20_bases','after_q30_bases','after_q20_rate','after_q30_rate','after_gc_content','after_duplication']
        everyKey="summary"

        value.append(data.get(everyKey).get('before_filtering').get('total_reads'))
        value.append(data.get(everyKey).get('before_filtering').get('total_bases'))
        value.append(data.get(everyKey).get('before_filtering').get('q20_bases'))
        value.append(data.get(everyKey).get('before_filtering').get('q30_bases'))
        value.append(data.get(everyKey).get('before_filtering').get('q20_rate'))
        value.append(data.get(everyKey).get('before_filtering').get('q30_rate'))
        value.append(data.get(everyKey).get('before_filtering').get('gc_content'))
        value.append(data.get('duplication').get('rate'))

        value.append(data.get(everyKey).get('after_filtering').get('total_reads'))
        value.append(data.get(everyKey).get('after_filtering').get('total_bases'))
        value.append(data.get(everyKey).get('after_filtering').get('q20_bases'))
        value.append(data.get(everyKey).get('after_filtering').get('q30_bases'))
        value.append(data.get(everyKey).get('after_filtering').get('q20_rate'))
        value.append(data.get(everyKey).get('after_filtering').get('q30_rate'))
        value.append(data.get(everyKey).get('after_filtering').get('gc_content'))
        value.append(data.get('duplication').get('rate'))
        extracted_info = pd.DataFrame(list(zip(key,value)), columns = [col_name,col_info])
        return extracted_info

    def Df_Append():
        tmp=pd.concat([df, extracted_info], axis=1)
        return tmp

    def File_path(*args):
        list=[]
        for file in args:
            list.append(file)
        return list
    def Writer(things,csv_name):
        current_dir = os.path.abspath(my_dir)
        file_name = os.path.join(current_dir,csv_name)
        things.to_csv(file_name,index=None,index_label=None)
    #my_dir="C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\fastp"
    os.chdir(my_dir)
    json_path=[]
    for root, dirs, files in os.walk(my_dir):
        for file in files:
            if file.endswith(".json"):
                json_path.append(os.path.join(root, file))
    df = pd.DataFrame()
    for i in range(len(json_path)):
        data=Read_json(json_path[i])
        extracted_info=Extract_info(json_path[i].split("\\")[-1]+"_name",json_path[i].split("\\")[-1]+"_info")
        df=Df_Append()
        Writer(df,os.path.join(root,outputName))


if __name__ == '__main__':
    # folder_name = "C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna\\sig"
    # resule_file  =  "C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna\\sig\\sig_result.xlsx"
    # merge(folder_name,resule_file)

    # folder_name = "C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna\\up"
    # resule_file  =  "C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna\\up\\up_result.xlsx"
    # merge(folder_name,resule_file)

    # folder_name = "C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna\\down"
    # resule_file  =  "C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna\\down\\down_result.xlsx"
    # merge(folder_name,resule_file)

    # folder_name = "C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna\\all"
    # result_file  =  "C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna\\all\\all_result.xlsx"
    # merge(folder_name,result_file)
    my_dir = "C:\\Users\\Jeff\\OneDrive\\桌面\\hzh\\rna\\qc"
    outputName = "qc.csv"
    json_analysis(my_dir,outputName)