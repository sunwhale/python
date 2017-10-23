# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

import numpy as np
import os
from Data import *
from Constants import *

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

def FindFirstHeaderLineNumber(list1):
    Headers = ['Count','Diaplacement']
    for i in range(min(len(list1),20)):
        for Header in Headers:
            if Header in list1[i]:
                return i
                
def FindFirstUnitLineNumber(list1):
    Units = ['kN','mm/mm']
    for i in range(min(len(list1),20)):
        for Unit in Units:
            if Unit in list1[i]:
                return i

def transform_to_database(inpath,outpath,name):
    experiment_log = ExperimentLog(ExperimentLogFile)
    experiment_log.output(name)
    regular = r'\d+\.?\d*'
    d_in = float(experiment_log.obtainItem(name,'d_in',regular)[0]) * 1e-3
    d_out = float(experiment_log.obtainItem(name,'d_out',regular)[0]) * 1e-3
    gauge_length = float(experiment_log.obtainItem(name,'gauge_length',regular)[0]) * 1e-3
    
    r_out = d_out/2.0
    r_in = d_in/2.0
    area = np.pi * (r_out**2 - r_in**2)

    if not os.path.isdir(outpath):
        os.makedirs(outpath)
        print 'Create new directory:',outpath
    
    filenames = getFiles(inpath)
    print '=============================='
#    print filenames
    
    Units = []
    Headers = []
    
    for filename in filenames:
        fullname = filename[0] + '\\' + filename[1]
        outname = fullname.replace('\\','_')
        outname = outname.replace(':','')
        outname = outname[:-4] + '.csv'
        if ('Strain Test 3 - Data Acquisition 1 - (Delta Level)' in outname) \
        or ('Strain Test 3 - Data Acquisition 1 - (Timed)' in outname) \
        or ('Data Acquisition 1 - (Timed)' in outname):
            outname = name + '.csv'
        if isTextFile(filename[1]):
            print '------------------------------'
            print 'reading: ', filename
            print 
            infile = open(fullname,'r')
            list1 = infile.readlines()
            FirstHeaderLineNumber = FindFirstHeaderLineNumber(list1)
            FirstUnitLineNumber = None
            if FirstHeaderLineNumber <> None:
                FirstUnitLineNumber = FirstHeaderLineNumber + 1
                list2 = [i.strip() for i in list1[FirstHeaderLineNumber:]]
                for i in range(1,len(list2)):
                    list2[i] = list2[i].replace('        ',' ')
                    list2[i] = list2[i].replace('       ',' ')
                    list2[i] = list2[i].replace('      ',' ')
                    list2[i] = list2[i].replace('     ',' ')
                    list2[i] = list2[i].replace('    ',' ')
                    list2[i] = list2[i].replace('   ',' ')
                    list2[i] = list2[i].replace('  ',' ')
                    list2[i] = list2[i].replace(' ','\t')
                list3 = [i.split('\t') for i in list2[2:] if i <>'' and i[0].isdigit()]
                Header = list2[0].split('\t')
                Unit = list2[1].split('\t')
                if len(Unit) > 1:
                    for H in Header:
                        Headers.append(H)
                    for U in Unit:
                        Units.append(U)
                    print 'Header:'
                    print Header
                    print
                    print 'Unit:'
                    print Unit
                    print
                else:
                    print Unit
                    print fullname
                print

#                Headers = list(set(Headers))
#                Units = list(set(Units))
#                print Headers
#                print Units

                del list1
                del list2                    
                lines = [ [0,0,0,0,0,0,0,0,0,0,0,0,0,0] for i in range(len(list3)) ]
                for i in range(len(Header)):
                    if 'CycleCount' in Header[i]:
                        for j in range(len(list3)):
                            lines[j][0] = list3[j][i]
                    if 'Running Time ' == Header[i]:
                        for j in range(len(list3)):
                            lines[j][1] = str(float(list3[j][i]) - float(list3[0][i]) + 1)
                    if 'Specimen Temperature' in Header[i] or 'Lower Furnace Temperature' in Header[i]:
                        for j in range(len(list3)):
                            lines[j][2] = list3[j][i]
                    if 'Axial Displacement' in Header[i]:
                        for j in range(len(list3)):
                            lines[j][3] = list3[j][i]
                    if 'Axial Force' in Header[i]:
                        for j in range(len(list3)):
                            lines[j][4] = list3[j][i]
#                    if 'Axial Total Strain' in Header[i]:
#                        for j in range(len(list3)):
#                            lines[j][5] = list3[j][i]
                    if 'Axial Mech. Strain' in Header[i] or 'Axial Strain' in Header[i]:
                        for j in range(len(list3)):
                            lines[j][5] = list3[j][i]
                    if 'Torsional Rotation' in Header[i] or 'Torsional Angle' in Header[i]:
                        for j in range(len(list3)):
                            lines[j][7] = list3[j][i]
                    if 'Torsional Strain (ang)' in Header[i]:
                        for j in range(len(list3)):
                            lines[j][9] = list3[j][i]
                    if 'Torsional Torque' in Header[i]:
                        for j in range(len(list3)):
                            lines[j][8] = list3[j][i]
                    if 'Thermal Strain' in Header[i]:
                        for j in range(len(list3)):
                            lines[j][12] = list3[j][i]
                    if 'Axial Total Strain' in Header[i] or 'Axial Strain' in Header[i]:
                        for j in range(len(list3)):
                            lines[j][13] = list3[j][i]
                for i in range(len(lines)):
#                    力的单位为牛顿
                    force = float(lines[i][4]) * 1000.0
#                    转角单位为度
                    rotation_deg = np.rad2deg(float(lines[i][7]))
#                    角应变单位为度
                    theta = float(lines[i][9])
#                    扭矩单位为N.m
                    torque = float(lines[i][8])
#                    剪应变单位为1
                    gama = r_out * np.deg2rad(theta) / gauge_length
#                    薄壁管剪应力
                    shear_stress = 16.0*torque/(np.pi*(d_out**2-d_in**2)*(d_out+d_in))/1.0e6
#                    alpha表示
#                    alpha = d_in/d_out
#                    shear_stress = 16.0*torque/(np.pi*(d_out**3)*(1-alpha**4))/1.0e6
                    lines[i][4] = str(force)
                    lines[i][6] = str(force / area / 1.0e6) # 正应力单位为MPa
                    lines[i][7] = str(rotation_deg) # 转角单位为度
                    lines[i][9] = str(gama) # 剪应变单位为1
                    lines[i][10] = str(shear_stress) # 剪应力单位为MPa
                outfile = open(outpath + outname, 'w')
                print 'writing: ' + outpath + outname
                print
                outfile.writelines('Axial Segment Count,Running Time,Temperature,Axial Displacement,Axial Force,Axial Strain,Axial Stress,Rotation,Torque,Angle Strain,Shear Stress,Equivalent Plastic Strain,Thermal Strain,Axial Total Strain\n')
                outfile.writelines('cycles,sec,C,mm,N,mm/mm,MPa,deg,N*m,-,MPa,-,mm/mm,mm/mm\n')
                for l in lines:
                    line = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11],l[12],l[13])
                    outfile.writelines(line)
                del list3
                del lines

name = '7210'
inpath = 'F:\\Database\\IN718\\Sun Jingyu\\%s\\' % name
outpath = 'F:\\Database\\IN718\\Sun Jingyu\\'
transform_to_database(inpath,outpath,name)
print 'end'

#for name in experiment_type_dict['BIAXIAL']:
#    inpath = 'F:\\Database\\SS304\\MTS\\%s\\' % name
#    outpath = 'F:\\Database\\SS304\\MTS_OUT\\'
#    transform_to_database(inpath,outpath,name)