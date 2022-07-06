import os
import re
from datetime import datetime
from io import StringIO
import glob
import pandas as pd
import pdfplumber

center=pd.read_excel("G:\\篮球赛发票\\unlock\\证明材料\\奖项统计v2.4.xlsx")
strlist1=center["姓名"].values.tolist()
strlist2=center["所属队伍"].values.tolist()
people_organ=dict(zip(strlist1,strlist2)) #{'k2': 'b', 'k1': 'a'}

class PDF:

    def __init__(self, path):
        self.path = path
        self.reader = pdfplumber.open(path)
        self.insurance = False
        # 检查pdf是否为保险
        if len(self.reader.pages) > 5:
            self.insurance = True
            print("yes")
        self.people = None
        self.start = None
        self.end = None

    def read_date(self):
        text = self.reader.pages[0].extract_text()
        start, end = re.findall(
            '(\d+年\d+月\d+日)',
            text,
            re.S | re.I
        )[1:]

        def x(s):
            return int(s[:4]), int(s[5:7]), int(s[8:10])

        start, end = x(start), x(end)
        iso = datetime
        self.start, self.end = iso(*start), iso(*end)
        return self.start, self.end

    def read_table(self, page_number=2):
        page = self.reader.pages[page_number]
        text = page.extract_text()
        buff = StringIO(text.strip('人员清单\n'))
        self.people = pd.read_csv(buff, sep=' ', index_col='序号')
        self.people['证件号码'] = self.people['证件号码'].astype(str)
        self.people['保险日期'] = self.start.strftime('%Y-%m-%d')
        self.people['所属队伍'] = self.people['被保险人'].map(lambda x: people_organ.get(x,x))
        self.people['签名'] = ''
        return self.people

    def to_csv(self, path, **kwargs):
        if not self.insurance:
            return
        self.people.to_csv(
            path,
            index=None,
            **kwargs
        )

    def __call__(self, *args, **kwargs):
        if self.insurance:
            self.read_date()
            self.read_table(*args, **kwargs)
        return self

    def __str__(self):
        return str(self.start)


def detect_pdfs(dirname, pdfs):
    pdfs += [
        i
        for i in
        glob.glob(os.path.join(dirname, '*pdf'))
    ]

    dirnames = [
        os.path.join(dirname, i)
        for i in os.listdir(dirname)
        if os.path.isdir(os.path.join(dirname, i))
    ]

    for dirname in dirnames:
        detect_pdfs(dirname, pdfs)


if __name__ == '__main__':
    # 输出的excel文件目录
    save_path = 'G:\\篮球赛发票\\unlock\\证明材料\\保险信息统计\\excel'
    insurance_path = 'G:\\篮球赛发票\\unlock\\证明材料\\保险信息统计\\'
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    pdfs = []
    detect_pdfs(insurance_path, pdfs)
    print(pdfs)
    for path in pdfs:
        pdf = PDF(path)()
        #print(pdf)
        if pdf.people is None:
            continue
        pdf.to_csv(
            os.path.join(
                save_path,
                pdf.start.strftime('%Y%m%d') +
                os.path.basename(path)
            ) + '.csv',
            columns = ['被保险人','性别','所属队伍','保险日期','签名']
        )
        print(path)