############################################################################################
# Code to find coincidences within between different channels of the CAEN DAQ. 
# 
# This code iterates upon loadCAEN4.py that was created by Philippe.
#
# Created: 2024_05_01
# Author: Nick Swidinsky
###########################################################################################
#Math and basic python functionality packages
import os
import numpy as np
import pandas as pd
from scipy.stats import poisson, norm, truncnorm
from uncertainties import ufloat

#Graphing packages
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import collections #Used for bar graphs

#Other useful packages
import pathlib #Package for paths
import time #For timing
import inspect #Used to get function names
import xml.etree.ElementTree as ET # Parsing xml
import subprocess #used for C++ calls

##################################################
# Global Variables 
##################################################

ps_ns = 1.E-3 #Pico seconds to nano seconds
ps_s = 1.E-12 #Pico seconds to seconds
ms_s = 1.E-3 #ms seconds to seconds
s_h = 1./3600. #seconds to hours


class runInfo(object):
    def __init__(self, fullInFile=inpath, outpath=None, verbose=False):
        ''' A class that wrangles paths and file names, creates output directory if necessary.
        Ex: dir1/dir2/file.CSV
        Works for an arbitrary number of channels
        '''
        if verbose: print(locals())
            
        self.fullRoot = os.path.splitext(fullInFile)[0] # Store the full path to the input file. (example: dir1/dir2/file)
        self.directory, self.fileName = os.path.split(fullInFile) # Splits the full path of the input file. Returns home directory and file name. (ex: dir1/dir2 file.CSV)
        
        self.outPath = self.directory if outpath is None else outpath #Store the filepath to the output directory. If set to None then nothing is stored.
        if not os.path.exists(self.outPath): os.makedirs(self.outPath) # Creates the directory of the output file path if it doen't exist.
        print("Creating output directory " + str(self.outPath)) #Prints the output directory to terminal. 
        
        self.figTit = self.directory + "\n" + self.fileName
        self.baseName = os.path.splitext(self.fileName)[0] # file
        
        
        # Sorted data
        self.coinSortedNameCSV = self.baseName + "_coinSorted" + ".csv" # file_coinSorted.csv
        self.coinSortedName = self.baseName + "_coinSorted" + ".pd" # file_coinSorted.pd
        self.sortedName = self.baseName + "_sorted" + ".pd" # file_sorted.pd # May be deprecated
        self.fullSortedName = os.path.join(self.directory, self.sortedName) # dir1/dir2/file_sorted.pd  # May be deprecated
        self.fullCoinSortedNameCSV = os.path.join(self.directory, self.coinSortedNameCSV) # dir1/dir2/file_coinSorted.csv
        self.fullCoinSortedName = os.path.join(self.directory, self.coinSortedName) # dir1/dir2/file_coinSorted.pd
        
        self.figTitSorted = self.directory + "\n" + self.sortedName

        #self.sortedCleanName = self.baseName + "_sortedClean" + ".pd" # file_sortedClean.pd
        #self.fullSortedCleanName = os.path.join(self.directory, self.sortedCleanName) # dir1/dir2/file_sortedClean.pd
        
        
        
        # Setup xmlfile
        #import pathlib # For paths
        self.xmlFolder = pathlib.Path(self.directory).parent
        self.xmlName = "settings.xml"
        self.fullXmlName = os.path.join(self.xmlFolder, self.xmlName)
        with open(self.fullXmlName) as file:
            data = file.read()
        p1 = data.split('<timedRunDuration>')[1]
        p2 = p1.split('</timedRunDuration>')[0]
        self.runDuration_ms = int(p2)
        self.runDuration_s = self.runDuration_ms * ms_s
        self.runDuration_h = self.runDuration_s * s_h
        #print(self.runDuration_h)
        #<timedRunDuration>86400000</timedRunDuration>
        
        def getChList(strg, maxCh=20):
            ''' Since I don't know how to parse xml, brute force method to determine number of channels, and their numbers.
            strg: the xml file read as a string
            maxCh: the highest possible channel number
            returns the list of channels
            '''
            t2 = '            <values>\n'
            t3 = '                <entry>\n'
            t4 = '                    <key>SRV_PARAM_CH_ENABLED</key>\n'
            t5 = '                    <value xsi:type="xs:boolean" xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">true</value>\n'

            chList = []
            for i in range(maxCh):
                t1 = '        <channel>\n            <index>' + str(i) + '</index>\n'
                target = t1 + t2 + t3 + t4 + t5
                if strg.find(target) > 0: chList.append(i)

            return chList
        
        self.chList = getChList(data)



def setPars():  ####### Might be able to not put this in a function, but define it globally. 
    ###########################################################################################
    # Sets universal parameters for matplotlib graphs. 
    ###########################################################################################
    matplotlib.rcParams["font.size"] = 17
    matplotlib.rcParams["lines.linewidth"] = 3
    matplotlib.rcParams["mathtext.default"] = 'regular'
    matplotlib.rcParams['lines.markersize'] = 3
    matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["b", "r", "g", 'm'])



setUpFile = open('Coincidence_init.txt','r') #opens the init file
setUpLines = setUpFile.readlines() #returns a list of strings from the setUpFile. 
setUpFile.close() #Closes the setUpFile

try:
    fileLoc = setUpLines[0].split("\t")[1].split("\n") #Pulls out the filepath from the first line of the init file.
    filepath = f'{fileLoc[0]}'
    coincidence = setUpLines[1].split() #pulls out the coincidenc window from the init file. 
except:
    print("Information in the Coincidence_init text file are not separeted by a tab from the variable name.")
    
    
