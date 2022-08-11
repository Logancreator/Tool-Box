# import pandas as pd
# import os
# os.chdir("C:\\Users\\Jeff\\OneDrive\\桌面\\")

# file_path = 'C:\\Users\\Jeff\\OneDrive\\桌面\\GCF_002166845.1_GooseV1.0_genomic.gtf'

# raw_gff = pd.read_table(file_path,header=None,sep="\t",comment="#")

# # print(raw_gff.head)

# # print(raw_gff.columns)

# raw_gff.columns = ['chr','source','type','start','end','other','strand','another','annotation']

# df = raw_gff["annotation"].str.split(';', expand=True)

# df.to_csv("annotation.csv")




# print(df.iloc[:,6])
# from cv2 import log
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.optimize import curve_fit
# import math
# #单个高斯模型，如果曲线有多个波峰，可以分段拟合
# def func(x, a,u, sig):
#     return a*np.exp(-(x - u) ** 2 / (2 * sig ** 2)) / (sig * math.sqrt(2 * math.pi))
# #混合高斯模型，多个高斯函数相加
# def func3(x, a1, a2, a3, m1, m2, m3, s1, s2, s3):
#     return a1 * np.exp(-((x - m1) / s1) ** 2) + a2 * np.exp(-((x - m2) / s2) ** 2) + a3 * np.exp(-((x - m3) / s3) ** 2)

#正弦函数拟合
#def fmax(x,a,b,c):
#    return a*np.sin(x*np.pi/6+b)+c
#fita,fitb=optimize.curve_fit(fmax,x,ymax,[1,1,1])
#非线性最小二乘法拟合
#def func(x, a, b,c):
#    return a*np.sqrt(x)*(b*np.square(x)+c)
#用3次多项式拟合，可推广到n次多项式，数学上可以证明，任意函数都可以表示为多项式形式
#f1 = np.polyfit(x, y, 3)
#p1 = np.poly1d(f1)
#yvals = p1(x)  #拟合y值
#也可使用yvals=np.polyval(f1, x)
#
#拟合，并对参数进行限制，bounds里面代表参数上下限，p0是初始范围，默认是[1,1,1]
# '''
# Author: CloudSir
# Date: 2021-08-01 13:40:50
# LastEditTime: 2021-08-02 09:41:54
# LastEditors: CloudSir
# Description: Python拟合多项式
# https://github.com/cloudsir
# '''
# import matplotlib.pyplot as plt
# import numpy as np
# import math
# import numpy

# x = [math.pow( 10, -11 ),math.pow( 10, -10 ),math.pow( 10, -9 ),math.pow( 10, -8 ),math.pow( 10, -7 )]
# x = -numpy.log10(x)
# y = [2275,1900,1000,640,391]
# z1 = np.polyfit(x, y, 10)  #用3次多项式拟合，输出系数从高到0
# p1 = np.poly1d(z1) #使用次数合成多项式
# y_pre = p1(x)
 
# plt.plot(x,y,'.')
# plt.plot(x,y_pre)
# plt.show()
