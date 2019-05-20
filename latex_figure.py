# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 12:23:58 2018

@author: J.Sun
"""

import re
import os
import shutil

def read_file(file_name):
    with open(file_name,'r') as data_file:
        return data_file.read()

def getFiles(path):
    fns=[]
    for root, dirs, files in os.walk( path ):
        for fn in files:
            fns.append( [ root , fn ] )
    return fns

def isTextFile(filename):
    b = False
    suffixs = ['txt','dat','log']
    for suffix in suffixs:
        if filename.split('.')[-1] == suffix:
            b = True
    return b

#latex_name = 'F:\\GitHub\\tgmf\\tgmf.tex'
#figures_directory = 'F:\\GitHub\\tgmf\\figures\\'

#latex_name = 'F:\\GitHub\\tmf\\tmf.tex'
#figures_directory = 'F:\\GitHub\\tmf\\figures\\'

#latex_name = 'F:\\GitHub\\fatigue\\fatigue.tex'
#figures_directory = 'F:\\GitHub\\fatigue\\figures\\'

latex_name = 'F:\\GitHub\\thermal\\thermal.tex'
figures_directory = 'F:\\GitHub\\thermal\\figures\\'
    
lines = read_file(latex_name)

"""
总结
^ 匹配字符串的开始。
$ 匹配字符串的结尾。
\b 匹配一个单词的边界。
\d 匹配任意数字。
\D 匹配任意非数字字符。
x? 匹配一个可选的 x 字符 (换言之，它匹配 1 次或者 0 次 x 字符)。
x* 匹配0次或者多次 x 字符。
x+ 匹配1次或者多次 x 字符。
x{n,m} 匹配 x 字符，至少 n 次，至多 m 次。
(a|b|c) 要么匹配 a，要么匹配 b，要么匹配 c。
(x) 一般情况下表示一个记忆组 (remembered group)。你可以利用 re.search 函数返回对象的 groups() 函数获取它的值。
正则表达式中的点号通常意味着 “匹配任意单字符”
"""

figure_list = []
figure_directory_list = []
suffix_list = ['.pdf','.png','.eps','.jpg']

for suffix in suffix_list:
#    regular = '{\w+-?\w?%s}' % (suffix)
    regular = '{[a-z|_|+|-|A-Z|0-9|-]+%s}' % (suffix)
    result_list = re.findall(regular, lines)
    for result in result_list:
        print result
        figure_list.append( result[1:-1] )    

regular = '\graphicspath{\S*}'
result_list = re.findall(regular, lines)
for result in result_list:
    regular = '{[a-z|A-Z|0-9|/|:]+}'
    for r in re.findall(regular, result):
        figure_directory_list.append( r[1:-1] )

figure_file_list = []
for figure_directory in figure_directory_list:
    for figure in getFiles(figure_directory):
        if figure[1] in figure_list:
            if figure[0][-1] <> '/':
                figure_file_list.append(figure[0]+'/'+figure[1])
            else:
                figure_file_list.append(figure[0]+figure[1])
                
figure_file_list = list(set(figure_file_list))
figure_file_list.sort()

if not os.path.isdir(figures_directory):
    os.makedirs(figures_directory)
    print 'Create new directory:',figures_directory
            
for figure_file in figure_file_list:
    shutil.copy(figure_file,figures_directory)