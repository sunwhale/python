# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

import zipfile
import os,sys


def unzip_dir(zipfilename, unzipdirname):
    fullzipfilename = os.path.abspath(zipfilename)
    fullunzipdirname = os.path.abspath(unzipdirname)
    print "Start to unzip file %s to folder %s ..." % (zipfilename, unzipdirname)
    #Check input ...
    if not os.path.exists(fullzipfilename):
        print "Dir/File %s is not exist, Press any key to quit..." % fullzipfilename
        inputStr = raw_input()
        return
    if not os.path.exists(fullunzipdirname):
        os.mkdir(fullunzipdirname)
    else:
        if os.path.isfile(fullunzipdirname):
            print "File %s is exist, are you sure to delet it first ? [Y/N]" % fullunzipdirname
            while 1:
                inputStr = raw_input()
                if inputStr == "N" or inputStr == "n":
                    return
                else:
                    if inputStr == "Y" or inuptStr == "y":
                        os.remove(fullunzipdirname)
                        print "Continue to unzip files ..."
                        break
            
    #Start extract files ...
    srcZip = zipfile.ZipFile(fullzipfilename, "r")
    for eachfile in srcZip.namelist():
        print "Unzip file %s ..." % eachfile
        eachfilename = os.path.normpath(os.path.join(fullunzipdirname, eachfile))
        eachdirname = os.path.dirname(eachfilename)
        if not os.path.exists(eachdirname):
            os.makedirs(eachdirname)
        fd=open(eachfilename, "wb")
        fd.write(srcZip.read(eachfile))
        fd.close()
    srcZip.close()
    print "Unzip file succeed"
    
def zip_dir(dirname, zipfilename):
    filelist = []
    #Check input ...
    fulldirname = os.path.abspath(dirname)
    fullzipfilename = os.path.abspath(zipfilename)
    print "Start to zip %s to %s ..." % (fulldirname, fullzipfilename)
    if not os.path.exists(fulldirname):
        print "Dir/File %s is not exist, Press any key to quit..." % fulldirname
        inputStr = raw_input()
        return
    if os.path.isdir(fullzipfilename):
        tmpbasename = os.path.basename(dirname)
        fullzipfilename = os.path.normpath(os.path.join(fullzipfilename, tmpbasename))
    if os.path.exists(fullzipfilename):    
        print "%s has already exist, are you sure to modify it ? [Y/N]" % fullzipfilename
        while 1:
            inputStr = raw_input()
            if inputStr == "N" or inputStr == "n" :
                return
            else:
                if inputStr == "Y" or inputStr == "y" :
                    print "Continue to zip files..."
                    break

    #Get file(s) to zip ...
    if os.path.isfile(dirname):
        filelist.append(dirname)
        dirname = os.path.dirname(dirname)
    else:
        #get all file in directory
        for root, dirlist, files in os.walk(dirname):
            for filename in files:
                filelist.append(os.path.join(root,filename))

    #Start to zip file ...
    destZip = zipfile.ZipFile(fullzipfilename, "w", zipfile.ZIP_DEFLATED)
    for eachfile in filelist:
        destfile = eachfile[len(dirname):]
        print "Zip file %s..." % destfile
        destZip.write(eachfile, destfile)
    destZip.close()
    print "Zip folder succeed!"
    
def print_help(toolname):
    print """
    This program can zip given folder to destination file, or unzip given zipped file to destination folder.
    Usage: %s [option] [arg]...
    -h: print this help message and exit (also --help)
    -u unzip: unzip given zipped file to destination folder,
        usage: %s -u/unzip zipfilename, unzipdirname
    -z zip: zip given folder to destination file
        usage: %s -z/zip dirname, zipfilename
    """  % (toolname, toolname, toolname)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print_help(sys.argv[0])
        sys.exit()
    if sys.argv[1].startswith('--'):
        option = sys.argv[1][2:]
        if option == 'version':
            print 'Get no version information'
        else:
            if option != 'help':
                print "Unknown option"
            print_help(sys.argv[0])
        sys.exit()
    if sys.argv[1].startswith('-'):
        if len(sys.argv) < 4 :
            print_help(sys.argv[0])
            sys.exit()
        option = sys.argv[1][1:]
        if option == 'u' or option == 'unzip':
            zipfilepath = sys.argv[2]
            unzipdirpath = sys.argv[3]
            unzip_dir(zipfilepath, unzipdirpath)
        else:
            if option == 'z' or option == 'zip':
                dirpath = sys.argv[2]
                zipfilepath = sys.argv[3]
                zip_dir(dirpath, zipfilepath)
            else:
                print_help(sys.argv[0])
                sys.exit()

#zip_dir('F:\\Temp\\SS304_Sim\\3001.csv', 'F:\\Temp\\SS304_Sim\\3001.zip')

#zip_dir('F:\\Temp\\SS304_Sim\\3001.csv', 'F:\\Temp\\SS304_Sim\\3001.zip')
#unzip_dir('F:\\Temp\\SS304_Sim\\3001.zip', 'F:\\Temp\\SS304_Sim\\')