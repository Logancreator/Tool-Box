import json
import os
from pickle import NONE
import pandas as pd

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
    current_dir = os.path.abspath('E:\\github\\Interesting code\\Interesting\\json_analysis\\')
    file_name = os.path.join(current_dir,csv_name)
    things.to_csv(file_name,index=None,index_label=None)

if __name__ == '__main__':
    my_dir="E:\\github\\Interesting code\\Interesting\\json_analysis\\json\\"
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
        Writer(df,os.path.join(root,"cuttag_fastp_json.csv"))



