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

def calBlock(string):
    return string.count('{') - string.count('}')

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

latex_name_list = ['F:\\GitHub\\doctor\\chapter1\\chapter1.tex',
                   'F:\\GitHub\\doctor\\chapter3\\chapter3.tex',
                   'F:\\GitHub\\doctor\\chapter4\\chapter4.tex',
                   'F:\\GitHub\\doctor\\chapter5\\chapter5.tex',
                   'F:\\GitHub\\doctor\\chapter6\\chapter6.tex',
                   'F:\\GitHub\\doctor\\chapter7\\chapter7.tex',
                   'F:\\GitHub\\doctor\\chapter8\\chapter8.tex',
                   'F:\\GitHub\\doctor\\chapter9\\chapter9.tex',
#                   'F:\\GitHub\\doctor\\chapter10\\chapter10.tex',
                   ]
                   
bib_name = 'F:\\GitHub\\doctor\\mybib1.bib'

label_list = []
for latex_name in latex_name_list:
    latex_lines = read_file(latex_name)
    regular = '\cite{\S+?}'
    result_list = re.findall(regular, latex_lines)
    for result in result_list:
        for label in result[5:-1].strip().split(','):
            label_list.append(label)
            
#label_list = sorted(set(label_list),key=label_list.index)
label_list = sorted(set(label_list))
for label in label_list:
    print label
    
bib_lines = read_file(bib_name)

bib_list = []

for bib_line in bib_lines.split('@')[1:]:
    bib_dict = {}
    
    regular = '\w+{\S*,'
    result_list = re.findall(regular, bib_line)
    classification = result_list[0].split('{')[0]
    label = result_list[0].split('{')[1][:-1]
    bib_dict['label'] = label
    bib_dict['classification'] = classification
    
    regular = '\w+ *= *["|{][^=\f\r\t\v]+["|}]'
    result_list = re.findall(regular, bib_line)
    for result in result_list:
        item = result.split('=')[0].strip()
        tmp = result[len(result.split('=')[0])+1:].strip()
#        print item,tmp
        if len(tmp) > 3:
            value = re.findall('[^{"]+[\S|\s]*[^ ,}\n"]+', tmp)
            if len(value) == 1:
                value = value[0]
            else:
                value = 'error'
        else:
            value = tmp[1:-1]
#        print item, value
        bib_dict[item] = value
        
    bib_list.append(bib_dict)

bib_list = sorted(bib_list, key=lambda bib_list:(bib_list['label']), reverse=False)

#for bib in bib_list:
#    if bib.has_key('year'):
#        print int(bib['year'])

bibtexfile = open('mybib.bib', 'w')
for bib in bib_list:
    if bib.has_key('year') and bib['label'] in label_list:
        print >>bibtexfile, ('@%s{%s,' % (bib['classification'],bib['label']))
        for key in bib.keys():
#            if key not in ['classification','label','abstract','file','keywords','url']:
            if key not in ['classification','label','abstract']:
                if calBlock(bib[key]) == 0:
                    print >>bibtexfile, ('  %s={%s},' % (key,bib[key]))
                elif calBlock(bib[key]) == 1:
                    print >>bibtexfile, ('  %s={%s}},' % (key,bib[key]))
                elif calBlock(bib[key]) == -1:
                    print >>bibtexfile, ('  %s={{%s},' % (key,bib[key]))
        print >>bibtexfile, '}'
bibtexfile.close()