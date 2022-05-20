import json
import os
from pickle import NONE
import pandas as pd

def read_json(file_path):
    with open(file_path) as f:
        return json.load(f)
    
def extract_info():
    info=[]
    name=['total_reads','total_bases','q20_bases','q30_bases','q20_rate','q30_rate','gc_content','duplication']
    everyKey="summary"
    info.append(data.get(everyKey).get('after_filtering').get('total_reads'))
    info.append(data.get(everyKey).get('after_filtering').get('total_bases'))
    info.append(data.get(everyKey).get('after_filtering').get('q20_bases'))
    info.append(data.get(everyKey).get('after_filtering').get('q30_bases'))
    info.append(data.get(everyKey).get('after_filtering').get('q20_rate'))
    info.append(data.get(everyKey).get('after_filtering').get('q30_rate'))
    info.append(data.get(everyKey).get('after_filtering').get('gc_content'))
    info.append(data.get('duplication').get('rate'))
    df_append = pd.DataFrame(list(zip(name,info)), columns = ['name','info'])
    return df_append

def df_Append():
    tmp=pd.concat([df, df_append], axis=1)
    return tmp

def writer(csv_name):
    current_dir = os.path.abspath('C:\\Users\\Jeff\\OneDrive\\桌面\\')
    file_name = os.path.join(current_dir,csv_name)
    df.to_csv(file_name,index=None,index_label=None)

if __name__ == '__main__':
    json_path=['C:\\Users\\Jeff\\OneDrive\\桌面\\B1P1_L2_801D73_fastp.json',
        'C:\\Users\\Jeff\\OneDrive\\桌面\\B1P1_L2_801D73_fastp.json']
    csv_path="fastp_json_extract_info.csv"
    df = pd.DataFrame()
    for i in range(len(json_path)):
        data=read_json(json_path[i])
        df_append=extract_info()
        df=df_Append()
        writer(csv_path)
