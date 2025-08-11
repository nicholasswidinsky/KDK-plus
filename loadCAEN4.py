#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Code to load CAEN data
# Digitizer and Compass documentation, including error guide: https://www.caen.it/?downloadfile=7161
# P. Di Stefano, 2023 08 21
#
# Standard use:
#
# runInfo -> ready
# loadData -> ready
# sortCoins -> ready, writes to csv format E00_AU T00_ns Flag00 ... E15_AU T15_ns Flag15
# loadSortedData -> ready
#
#
# OK flags are 4000, 0100 and combination 4100


# In[2]:


# Export automatically to py script
try:
    get_ipython() # Ensures export only runs when running interactively
    get_ipython().system('jupyter nbconvert --to script "loadCAEN.ipynb"')
except: pass


# In[3]:


import os
import numpy as np
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import pathlib # For paths

import time # For timing
import inspect # To get function names
from deprecated import deprecated # Not sure this works

import subprocess # FOr c++ calls

import collections # For bar charts of text data

from scipy.stats import poisson, norm, truncnorm

import pandas as pd

from uncertainties import ufloat

from num2tex import num2tex

import xml.etree.ElementTree as ET # Parsing xml
import matplotlib.ticker as ticker

try:
    import sys  
    sys.path.insert(0, '../Utilities')
    import MatplotlibGraphSettings as mgs
    libLoad = True
except ImportError:
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MatplotlibGraphSettings")
    libLoad = False


# In[4]:


if False: print(matplotlib.rcParams.keys())

def setPars():
    if False: print(matplotlib.rcParams.keys())
    matplotlib.rcParams["font.size"] = 17
    matplotlib.rcParams["lines.linewidth"] = 3
    matplotlib.rcParams["mathtext.default"] = 'regular'
    matplotlib.rcParams['lines.markersize'] = 3
    matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["b", "r", "g", 'm'])

setPars()    

#%matplotlib inline
#plt.rcParams['lines.linewidth'] = 3  

baseSep = "============================================================================================================\n======================================================> "


# In[5]:


#inpath = rootPath + '/2SiPM_coincidence_90d_50d_1day_1/FILTERED/SDataF_2SiPM_coincidence_90d_50d_1day_1.CSV'


#inpath = rootPath + '/2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1200_200lsb_15h_raw/FILTERED/SDataF_2SiPM_coincidence_80d_50d_240nsBoard_240nsCoMPASS_CoincidenceWindow_1200_200lsb_15h_raw.CSV'

#inpath = rootPath + '/2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1400lsb_100ns_25pc_cfd_15h_raw/FILTERED/SDataR_2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1400lsb_100ns_25pc_cfd_15h_raw.CSV'
#inpath = rootPath + '/2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1400lsb_100ns_25pc_cfd_15h_raw/FILTERED/SDataF_2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1400lsb_100ns_25pc_cfd_15h_raw.CSV'

ps_ns = 1.E-3
ps_s = 1.E-12
ms_s = 1.E-3
s_h = 1./3600.


rootPath = '/home/nick/MSC/Plastic_scintillators/Compton_scattering_tests/Room_temperature_reflective_foil/Plastic_Scintillator_wRF_10CoarseGain_900HV_2024_05_01'
inpath = rootPath 

# In[6]:


class runInfo(object):
    def __init__(self, fullInFile=inpath, outpath=None, verbose=False):
        ''' A class that wrangles paths and file names, creates output directory if necessary.
        Ex: dir1/dir2/file.CSV
        Works for an arbitrary number of channels
        '''
        if verbose: print(locals())
            
        self.fullRoot = os.path.splitext(fullInFile)[0] # dir1/dir2/file
        self.directory, self.fileName = os.path.split(fullInFile) # dir1/dir2 file.CSV
        
        self.outPath = self.directory if outpath is None else outpath
        if not os.path.exists(self.outPath): os.makedirs(self.outPath); print("Creating output directory " + str(self.outPath))
            
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

        
        if False:
            tree = ET.parse(os.path.join(self.directory.rstrip('FILTERED'), self.xmlName))
            #print(tree)
            #root = tree.getroot()
            #print(root)
            #for child in root:
            #    print(child.tag, child.attrib)
            for node in tree.find('acquisitionMemento'):
                print("------------->", node)
                #print("ooooooooooooo>", node.find('timedRunDuration').text)
                for elt in node: 
                    print("******>", elt)
                    print("+++++++++++++++++++++++>", elt.tag, elt.attrib)
                    #print("ooooooooooooo>", elt.find('timedRunDuration').text)
            #for neighbor in root.iter('acquisitionMemento'):
            #    print(neighbor.attrib)
            #for neighbor in root.iter('timedRunDuration'):
            #    print(neighbor.attrib)
            #print("Run info", vars(self))
            #for page in tree.findall('acquisitionMemento'):
            #    print(page)

    
#if True: ri = runInfo()       


# In[7]:


def stabilityPlots(df, path, titExtra="NO OFFLINE CUTS", nameExtra="00"):
    ''' Basic plots, including energy distribution, stability, and flags
    Works for an arbitrary number of channels
    df: the dataframe read from DAQ disk output
    path: the path where the df is
    '''
    ri = runInfo(path)

    channels = set(df['CHANNEL'])
    for ch in channels: dfCh = {'C' + str(ch):  df[df['CHANNEL']==ch] for ch in channels}
    
    fig, ax = plt.subplots(4, len(dfCh), figsize=(7. * len(dfCh), 12)) 
    ax = ax.ravel()
    fig.suptitle(ri.directory + "\n" + str(ri.fileName) + ", " + str(len(df.index)) + " events, " + titExtra)

    for i, ch in enumerate(channels):
        ax[i].set_title('C' + str(ch) + ', ' + str(len(dfCh['C' + str(ch)].index)) + ' events')
        ax[i].hist(dfCh['C' + str(ch)]['ENERGY'], bins=100, histtype='step'); ax[i].set_xlabel('ENERGY (AU)')
        ax[i+len(dfCh)].hist2d(dfCh['C' + str(ch)]['TIMETAG'], dfCh['C' + str(ch)]['ENERGY'], bins=[50, 50]); ax[i+len(dfCh)].set_xlabel('TIMETAG (ps)'); ax[i+len(dfCh)].set_ylabel('ENERGY (AU)')
        ax[i+2*len(dfCh)].hist2d(dfCh['C' + str(ch)]['TIMETAG'], dfCh['C' + str(ch)]['ENERGY'].index.to_list(), bins=[50, 50]); ax[i+2*len(dfCh)].set_xlabel('TIMETAG (ps)'); ax[i+2*len(dfCh)].set_ylabel('INDEX')
        cnt = collections.Counter(dfCh['C' + str(ch)]['FLAGS'])
        ax[i+3*len(dfCh)].bar(cnt.keys(), cnt.values(), color='b', lw=3)
        ax[i+3*len(dfCh)].set_xlabel('FLAGS'); ax[i+3*len(dfCh)].set_yscale('log'); ax[i+3*len(dfCh)].set_ylabel('Counts'); ax[i+3*len(dfCh)].tick_params(axis='x', labelrotation = 45)

    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "stabilityPlots() loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_Stability_" + nameExtra + ext)
    plt.show(); plt.clf()
        
        
def flagPlots(df, path, log=True):
    ''' Effect of flags on spectra of various channels
    Works for an arbitrary number of channels
    df: the dataframe read from DAQ disk output
    path: the path where the df is
    '''
    ri = runInfo(path)

    channels = set(df['CHANNEL'])
    for ch in channels: dfCh = {'C' + str(ch):  df[df['CHANNEL']==ch] for ch in channels}

    flags = list(set(df['FLAGS']))
    
    # Make sure first flags are in standard order with '0x4000' first and '0x4100' second
    for fl in ['0x4100', '0x4000']:
        if fl in flags: flags.remove(fl); flags.insert(0, fl)

    fig, ax = plt.subplots(1, len(dfCh)+1, figsize=(5. * (len(dfCh) + 1), 7)) 
    ax = ax.ravel()
    fig.suptitle(ri.directory + "\n" + str(ri.fileName) + ", " + str(len(df.index)) + " events, NO OFFLINE CUTS")

    for ch, axx in zip(channels, ax[1:]):
        ddf = dfCh['C' + str(ch)]
        eMax_AU = max(ddf['ENERGY'])
        eBins_AU = np.linspace(0, eMax_AU, 101)
        axx.hist(ddf['ENERGY'], bins=eBins_AU, log=log, histtype='step', lw=6, ls=':', color='k', label='All')

        for fl, col in zip(flags, ['r', 'b', 'g', 'orange', 'm', 'c', 'pink', 'gray', 'yellow', 'maroon', 'violet', 'teal', 'linen']):
            axx.hist(ddf[ddf['FLAGS']==fl]['ENERGY'], bins=eBins_AU, log=log, histtype='step', lw=3, color=col, label=fl)
        axx.set_title('C' + str(ch) + ', ' + str(len(ddf.index)) + ' events'); axx.set_xlabel('ENERGY (AU)'); axx.legend()
    cnt = collections.Counter(df['FLAGS'])
    ax[0].bar(cnt.keys(), cnt.values(), color='k')
    ax[0].set_xlabel('FLAGS'); ax[0].set_yscale('log'); ax[0].set_ylabel('Counts'); ax[0].tick_params(axis='x', labelrotation = 45); ax[0].set_title('All ch')
    
    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "flagPlots() loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_Flags_00" + ext)
    plt.show(); plt.clf()
    
    
    
def loadData(inpath=inpath, verbose=False, plot=True, drop=True, nrows=None):
    ''' Loads CAEN data from .CSV to a pandas dataframe
    Optionnally makes data quality plots
    Does not apply any cuts, though has the option to drop certain columns
    Works for an arbitrary number of channels
    inpath: path to dataframe
    drop: get rid of colums BOARD, ENERGYSHORT
    Returns total dataframe, as well as df for each channel, and the input path
    '''
    print(baseSep + str(inspect.currentframe().f_code.co_name) + "\n" + str(locals()) + "\n" + baseSep)
    
    ri = runInfo(inpath)
    
    df = pd.read_csv(inpath, sep=";", nrows=nrows)
    if verbose: print("------------------ All data\n", df)
    
    print(df['CHANNEL'])
    channels = set(df['CHANNEL'])
    if verbose: print("Input file: " + str(inpath) + "\n\tFound headers: " + str(df.columns) + "\n\tFound channels: " + str(channels) + "\n\tTotal evts: " + str(len(df['CHANNEL'])))
    if drop:
        df.drop('BOARD', inplace=True, axis=1); df.drop('ENERGYSHORT', inplace=True, axis=1)
        if verbose: print("Dropping cols BOARD and ENERGY SHORT"); print("------------------ Cols dropped\n", df)
    
    dfCh = {'C' + str(ch): df[df['CHANNEL']==ch] for ch in channels}
    
    if plot: 
        stabilityPlots(df=df, path=inpath) # Channel plots
        flagPlots(df=df, path=inpath, log=True) # Effect of flags on spectra
    
    return df, dfCh, inpath


# In[8]:


@deprecated("Not useful")
def timeDiffs(df, dfCh, outpath=inpath, plot=False):
    ''' Find wait times for various channels
    NOT SO USEFUL
    '''
    ri = runInfo(outpath)
    
    deltat_ns = ps_ns * np.diff(df['TIMETAG'])
    deltatCh_ns = {'All': deltat_ns}
    for ch in dfCh: deltatCh_ns[ch] = ps_ns * np.diff(dfCh[ch]['TIMETAG'])

    if plot:
        setPars()
        fig, ax = plt.subplots(1, 1, figsize=(8, 6)) 
        fig.suptitle(ri.figTit)
        for ch in deltatCh_ns: ax.hist(deltatCh_ns[ch], log=True, histtype='step', lw=3 if ch=='All' else 5, alpha=1.0 if ch=='All' else 0.3, label=ch)
        ax.legend(); ax.set_xlabel('Time between consecutive events (ns)')

        fig.tight_layout()
        if libLoad: mgs.progWaterMark(plt.gcf(), "timeDiffs() loadCAEN.ipynb")
        for ext in ['.pdf']: plt.savefig(ri.outPath + "_TEST_02" + ext)
        plt.show(); plt.clf()


# In[41]:


@deprecated("Data format has changed with c++")        
def plotCoinsBasic(df, outpath=inpath, coinWin_ns=100., titExtra="NO OFFLINE CUTS", nameExtra="00"):
    ''' Plots information about coincidences without the reordered file
    '''
    ri = runInfo(outpath)
    
    setPars()
    
    fig, ax = plt.subplots(1, 4, figsize=(12, 6)) 
    fig.suptitle(ri.figTit + " " + str(titExtra))
    
    ax[0].hist(df['NUM_COIN'], log=True, histtype='step')
    ax[0].set_xlabel('Multiplicity of coins')
    
    ax[1].hist2d(df['CHANNEL'], df['NUM_COIN'], bins=[np.linspace(-0.5, 3.5, 5), np.linspace(-0.5, 3.5, 5)], norm=matplotlib.colors.LogNorm() if True else None, cmap='Greys', rasterized=True)
    ax[1].set_xlabel('Ch'); ax[1].set_ylabel('Multiplicity of coins')
    
    ax[2].hist(df['Self_coin'], bins=np.linspace(-0.5, 1.5, 3), histtype='step', color='b', log=True)
    ax[2].set_xlabel('Self coincidence'); ax[2].set_ylabel('Counts'); ax[2].set_ylim(bottom=0.1)
    
    ax[3].hist(df['COIN_WAITS_first_ns'], log=True, histtype='step') #; ax[1].axvline(coinWin_ns, color='g'); ax[1].axvline(-coinWin_ns, color='g')
    ax[3].axvspan(0, coinWin_ns, alpha=0.1, color='k')
    ax[3].set_xlabel('Wait time to first coin (ns)')
    
    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "plotCoinsBasic() loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_Coins_" + str(nameExtra) + ext)
    plt.show(); plt.clf()
    
    
@deprecated("Not useful")    
def plotCoins(dfc, outpath=inpath, coinWin_ns=100., titExtra="NO OFFLINE CUTS", nameExtra="02"):
    ''' To display ordered data
    0 is the main channel
    IGNORE
    '''
    ri = runInfo(outpath)
    
    setPars()
    
    fig, ax = plt.subplots(1, 4, figsize=(18, 6)) 
    fig.suptitle(ri.figTit + " " + str(titExtra))
    
    ax[0].hist(dfc['ENERGY'], histtype='step', log=True, lw=3, color='b', label='All')
    for i, col in zip([2, 3], ['r', 'g']): ax[0].hist(dfc[dfc['COIN_CH_first']==i]['ENERGY'], histtype='step', log=True, lw=3, color=col, label='Coin C'+str(i))
    ax[0].legend(); ax[0].set_xlabel('C0 ENERGY'); ax[0].set_ylabel('Counts')
    
    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "plotCoins() loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_Coins_" + str(nameExtra) + ext)
    plt.show(); plt.clf()
    
@deprecated("Too slow, replaced by sortCoins() which calls c++")        
def findCoins(df, outpath=inpath, coinWin_ns=100., back=True, verbose=False, plot=False, cStart=None, clean=False):
    ''' Finds coins
    Returns a df with extra columns which are the channels coincident to the original one of the line
    Coins are counted as many times as their multiplicity
    back: if should look back for coins as well; otherwise, just looks forward
    cStart: if should use a given channel like a TDC start
    Returns an ordered file, written to disk as _sorted.pd.  Events where no coins are found are kept, but listed with NaNs.
    '''
    print(baseSep + str(inspect.currentframe().f_code.co_name) + "\n" + str(locals()) + "\n" + baseSep)
    
    ri = runInfo(outpath)
    
    dfc = df.copy()
    def getEvt(i): return {'i': dfc.index[i], 'Ch': dfc['CHANNEL'][i], 'E': dfc['ENERGY'][i], 't_ps': dfc['TIMETAG'][i], 'flag': dfc['FLAGS'][i]}
    
    coinChAll, coinChTimeAll_ps, coinChDiffAll_ns, coinEAll, coinFlagAll, coinCh0, coinChTime0_ps, coinChDiff0_ns, coinE0, coinFlag0, numCoin = [], [], [], [], [], [], [], [], [], [], [] # Coin ch, time, diff for all coins with given evt, for all events
    tagged = [0]*len(dfc.index); selfCoin = [False]*len(dfc.index)
    #print(tagged)
    
    if verbose: print("------------------ Before looking for coins\n", dfc)
    start_time = time.time() # Timer
    for i in dfc.index:
        if ((i <= 30) & (i % 10 == 0)) or ((i <= 300) & (i % 100 == 0)) or ((i <= 3000) & (i % 1000 == 0)) or ((i <= 30000) & (i % 10000 == 0)) or ((i % 100000 == 0)): 
            dt_s = time.time() - start_time
            if dt_s < 300: print("Sorting " + str(i) + "/" + str(len(dfc.index)) + str(" in %s seconds." % round(dt_s, 2)))
            elif dt_s < 3600 : print("Sorting " + str(i) + "/" + str(len(dfc.index)) + str(" in %s minutes." % round(dt_s/60., 2)))
            else : print("Sorting " + str(i) + "/" + str(len(dfc.index)) + str(" in %s hours." % round(dt_s/3600., 2)))
                
        ei = getEvt(i)
        if verbose: print("======================================================================================\n" + f"{i=}" + "\n" + str(ei))
            
        coinCh, coinChTime_ps, coinChDiff_ns, coinE, coinFlag = [], [], [], [], [] # Coin ch, time, diff for all coins with given evt
        
        if back:
            if verbose: print("-----------------------------BACK")
            j = i-1
            dt_ns = 0.

            while (np.abs(dt_ns) <= coinWin_ns) and (j >= 0): # Go backwards
                ej = getEvt(j)
                dt_ns = ps_ns * (ej['t_ps'] - ei['t_ps'])
                if (np.abs(dt_ns) <= coinWin_ns): coinCh.append(ej['Ch']); coinChTime_ps.append(ej['t_ps']); coinChDiff_ns.append(dt_ns); coinE.append(ej['E']); coinFlag.append(ej['flag']); tagged[j] += 1;
                if verbose: print("\t" + f"{j=}" + "\t" + str(getEvt(j)) + "\t" + str(dt_ns))
                j = j-1
            
        if verbose: print("-----------------------------FORW")
        j = i+1
        dt_ns = 0.

        while (np.abs(dt_ns) <= coinWin_ns) and (j <= dfc.index[-1]): # Go forwards
            ej = getEvt(j)
            dt_ns = ps_ns * (ej['t_ps'] - ei['t_ps'])
            if (np.abs(dt_ns) <= coinWin_ns): coinCh.append(ej['Ch']); coinChTime_ps.append(ej['t_ps']); coinChDiff_ns.append(dt_ns); coinE.append(ej['E']); coinFlag.append(ej['flag']); tagged[j] += 1
            if verbose: print("\t" + f"{j=}" + "\t" + str(getEvt(j)) + "\t" + str(dt_ns))
            j = j+1
        
        if ei['Ch'] in coinCh: selfCoin[i] = True
            
        coinChAll.append(coinCh); coinChTimeAll_ps.append(coinChTime_ps); coinChDiffAll_ns.append(coinChDiff_ns); coinEAll.append(coinE); coinFlagAll.append(coinFlag)
        if len(coinCh) > 0: coinCh0.append(coinCh[0]); coinChTime0_ps.append(coinChTime_ps[0]); coinChDiff0_ns.append(coinChDiff_ns[0]); coinE0.append(coinE[0]); coinFlag0.append(coinFlag[0])
        else: coinCh0.append(np.NAN); coinChTime0_ps.append(np.NAN); coinChDiff0_ns.append(np.NAN); coinE0.append(np.NAN); coinFlag0.append(np.NAN)
        numCoin.append(len(coinChDiff_ns) + 1)
            
    dfc['NUM_COIN'] = numCoin
    dfc['COIN_CH'] = coinChAll
    dfc['COIN_TIMES_ps'] = coinChTimeAll_ps
    dfc['COIN_WAITS_ns'] = coinChDiffAll_ns
    dfc['COIN_CH_first'] = coinCh0; dfc = dfc.astype({'COIN_CH_first':'Int64'}) # Converts ints to floats because of NaN
    dfc['COIN_TIMES_first_ps'] = coinChTime0_ps
    dfc['COIN_WAITS_first_ns'] = coinChDiff0_ns
    dfc['ENERGY_first'] = coinE0
    dfc['Tagged'] = tagged
    dfc['Self_coin'] = selfCoin
    dfc['COIN_FLAGS'] = coinFlagAll
    dfc['COIN_FLAG_first'] = coinFlag0
    
    if verbose: 
        print("------------------ Data sorted for coins\n", dfc)
        print("------------------ Triple coins\n", dfc[dfc['NUM_COIN']==3])
    
    if not cStart is None:
        print("------------------ Only starts being kept: C" + str(cStart))
        dfc = dfc[dfc['CHANNEL']==cStart]
    
    if plot:
        extraTit = "SORTED, " + ("ALL STARTS" if cStart is None else "C" + str(cStart) + " STARTS ONLY")
        extraName = "SORTED_" + ("ALL_STARTS" if cStart is None else "C" + str(cStart) + "_STARTS_ONLY")                         
        #plotCoins(dfc, outpath=outpath, coinWin_ns=coinWin_ns, titExtra=extraTit, nameExtra=extraName + "_02")
        plotCoinsBasic(dfc, outpath=outpath, coinWin_ns=coinWin_ns, titExtra=extraTit, nameExtra=extraName + "_00")
        stabilityPlots(df=dfc, path=outpath, titExtra=extraTit, nameExtra=extraName)
        
    if clean:
        dfc.drop('COIN_CH', inplace=True, axis=1); dfc.drop('COIN_TIMES_ps', inplace=True, axis=1); dfc.drop('COIN_WAITS_ns', inplace=True, axis=1); dfc.drop('COIN_FLAGS', inplace=True, axis=1)
        if verbose: print("Dropping cols COIN_CH, COIN_TIMES_ps, COIN_FLAGS and COIN_WAITS_ns"); print("------------------ Dropping cols\n", df)
    
    dfc.to_pickle(ri.fullSortedName)
    return dfc


def plotCoinsSimple(df, ri, coinWin_ns=[-100., 100.], titExtra="W/O FLAG CUTS", nameExtra="000", chLst=None, okFlags=['0x4000', '0x4100']):
    ''' Plots information about coincidences without the reordered file
    Channel list needs to be input by hand
    '''
    #ri = runInfo(outpath)
    
    if chLst is None: chLst = ri.chList

    dfc = {'All': df}
    for ch in chLst: dfc['C' + str(ch)] = df[df['E'+str(ch) + '_AU'] > 0] # df for each specific channel in which there is an event

    shftCh = chLst[-1:] + chLst[:-1] # Put last ch in first position
    dfc2 = {'C'+str(ch) + '-C' + str(och): df[(df[('E' + str(ch) + '_AU')] > 0) & (df[('E' + str(och) + '_AU')] > 0)] for ch, och in zip(chLst, shftCh)} # Double coins
        
    okFlagsStr = "|".join(okFlags)
        
    #################
    # First set of plots
    fig, axx = plt.subplots(3, 1+len(chLst), figsize=(18, 10)) 
    fig.suptitle(ri.figTit + " " + str(titExtra))
    
    
    multMax = max(dfc['All']['Mult'])
    bins = np.linspace(-0.5, multMax+0.5, multMax+2)
    for k, col in zip(dfc, ['k', 'b', 'r', 'g']): axx[0][0].hist(dfc[k]['Mult'], bins=bins, histtype='step', lw=3, color=col, ls='-' if col=='k' else ':', log=True, label=k)
#    axx[0][0].set_xticks(range(multMax+1), labels=[str(i) for i in range(multMax+1)])
    axx[0][0].set_xlabel('Multiplicity'); axx[0][0].set_ylabel('Counts'); axx[0][0].legend()
    
    nBins = 100
    tBins_ns = np.linspace(1.2*coinWin_ns[0], 1.2*coinWin_ns[1], nBins)
    for i, ch in enumerate(chLst):
        axx[0][i+1].set_title('C' + str(ch))
        axx[0][i+1].hist(dfc['C' + str(ch)]['E' + str(ch) + '_AU'], bins=100, log=True, histtype='step', lw=6, color='k', label='C' + str(ch) + ' all')
        
        otherChLst = [c for c in chLst if not (c==ch)] # List of the other channels
        for och, col in zip(otherChLst, ['b', 'r', 'g']):
            dfcc = dfc['C' + str(ch)][dfc['C' + str(ch)]['E' + str(och) + '_AU']>0] # Coins between the two channels
            axx[0][i+1].hist(dfcc['E' + str(ch) + '_AU'], bins=nBins, log=True, histtype='step', lw=3, color=col, label='C'+str(ch) + '-' + 'C'+str(och))
            axx[1][i+1].hist((dfcc['T' + str(ch) + '_ps']-dfcc['T' + str(och) + '_ps'])*ps_ns, bins=tBins_ns, log=True, histtype='step', lw=3, color=col, label='C'+str(ch) + '-' + 'C'+str(och))
            axx[1][i+1].axvspan(coinWin_ns[0], coinWin_ns[1], alpha=0.1, color='k')
            
        axx[0][i+1].set_xlabel('C' + str(ch) + ' energy (AU)'); axx[0][i+1].legend()
        axx[1][i+1].set_xlabel('C' + str(ch) + ' wait times (ns)'); axx[1][i+1].legend()
        
    #print(dfc2)
    for ax, ch, och in zip(axx[-1][1:], chLst, shftCh):
        ax.hist2d(dfc2['C' + str(ch) + '-C' + str(och)]['E' + str(ch) + '_AU'], dfc2['C' + str(ch) + '-C' + str(och)]['E' + str(och)+'_AU'], bins=[nBins, nBins], norm=matplotlib.colors.LogNorm() if True else None, cmap='Greys', rasterized=True)
        ax.set_xlabel('C' + str(ch) + ' energy (AU)'); ax.set_ylabel('C' + str(och) + ' energy (AU)')
    
    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "plotCoinsSimple() loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_Coins_" + str(nameExtra) + ext)
    plt.show(); plt.clf()
    
    
    #################
    ### Biplots
    fig, ax = plt.subplots(len(chLst), len(chLst), figsize=(18, 18)) 
    fig.suptitle(ri.figTit + " " + str(titExtra))
        
    for i in range(len(chLst)):
        chi = chLst[i]
        
        dfi = df[df[('E' + str(chi) + '_AU')] > 0]
        dfiCut = dfi[dfi['Flag' + str(ch)].str.contains(okFlagsStr)]
        
        Eimax = max(dfi['E' + str(chi) + '_AU'])
        ax[i][i].hist(dfi['E' + str(chi) + '_AU'], bins=np.linspace(0, Eimax, nBins), log=True, histtype='step', lw=6, ls='-', color='k', label='E' + str(chi) + '>0')
        ax[i][i].hist(dfiCut['E' + str(chi) + '_AU'], bins=np.linspace(0, Eimax, nBins), log=True, histtype='step', lw=6, ls=':', color='k', alpha=0.3)
        
        for j, col in zip(range(len(chLst)), ['r', 'b', 'g', 'm', 'c']):
            chj = chLst[j]

            dfij = df[(df[('E' + str(chi) + '_AU')] > 0) & (df[('E' + str(chj) + '_AU')] > 0)] # Double coins
            dfijCut = dfij[(dfij['Flag' + str(chi)].str.contains(okFlagsStr)) & (dfij['Flag' + str(chj)].str.contains(okFlagsStr))]
            
            if j>i:
                ax[i][j].hist2d(dfij['E' + str(chi) + '_AU'], dfij['E' + str(chj) + '_AU'], bins=[nBins, nBins], norm=matplotlib.colors.LogNorm() if False else None, cmap='Greys', rasterized=True)
                ax[i][j].set_xlabel('E' + str(chi) + ' (AU)'); ax[i][j].set_ylabel('E' + str(chj) + ' (AU)')

                ax[j][i].hist((dfij['T' + str(chi) + '_ps']-dfij['T' + str(chj) + '_ps'])*ps_ns, bins=tBins_ns, log=True, histtype='step', lw=3, color='k', label='C'+str(chi) + '-' + 'C'+str(chj))
                ax[j][i].axvspan(coinWin_ns[0], coinWin_ns[1], alpha=0.1, color='k')
                ax[j][i].set_xlabel('T'+str(chi) + ' - ' + 'T'+str(chj) + ' (ns)')
                ax[j][i].xaxis.set_major_locator(ticker.LinearLocator(9))
                
            if j!=i: 
                ax[i][i].hist(dfij['E' + str(chi) + '_AU'], bins=np.linspace(0, Eimax, nBins), log=True, histtype='step', lw=2, ls='-', color=col, label='E' + str(chi) + '>0, E' + str(chj) + '>0')
                ax[i][i].hist(dfijCut['E' + str(chi) + '_AU'], bins=np.linspace(0, Eimax, nBins), log=True, histtype='step', lw=2, ls='-', color=col, alpha=0.3)
            
            #ax[i][j].legend(fontsize=10)
            
        ax[i][i].set_xlabel('E' + str(chi) + ' (AU)')    
        ax[i][i].legend(fontsize=10)
    
    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "plotCoinsSimple() loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_CoinsBiplots_" + str(nameExtra) + ext)
    plt.show(); plt.clf()
    
    
def sortCoins(inPath, coinWin_ns=[-100., 100.], plot=True):
    ''' Calls c++ code sortCAENdata.cpp to find coincidences in original RAW or FILTERED csv data.
    The c++ is a lot faster than the python implementation.
    Output format of csv is (tab-delimited):
    E00_AU T00_ns Flag00 ... E15_AU T15_ns Flag15
    Channels that do not participate in a coin are given values of 0.
    '''
    print(baseSep + str(inspect.currentframe().f_code.co_name) + "\n" + str(locals()) + "\n" + baseSep)
    
    # ri = runInfo(inPath)
    
    # compCmd = "g++ -o sortCAENdata sortCAENdata.cpp"
    # runCmd = "./sortCAENdata" + " " + inPath + " " + ri.fullCoinSortedNameCSV + " " + str(coinWin_ns[0]) + " " + str(coinWin_ns[1])

    # a = subprocess.run(compCmd.split(" "), capture_output=True) # Compile the c++
    # print(a)
    # b = subprocess.run(runCmd.split(" "), capture_output=True) # Run the c++ output
    # print(b)
    
    # Read the sorted data to df and write to disk as df after adding multiplicity
    #dfc = pd.read_csv(ri.fullCoinSortedNameCSV, sep="\t")
    print(f'C++ (sortCoin) Input file: {inPath}')
    dfc = pd.read_csv(inPath,delimiter="\t")
    print(dfc)
    
    eCols = [col for col in dfc.columns if '_AU' in col] # Get energy columns
    dfc['Mult'] = dfc[eCols].gt(0).sum(axis=1) # Add multiplicity to each line
    dfc.to_pickle(ri.fullCoinSortedName) # Write to disk
    print("Coincidence sorted data:\n", dfc)
    
    colLgths = [len(dfc[dfc['T' + str(i) + '_ps']>0]['T' + str(i) + '_ps']) for i in range(6)] # Check that number of events was conserved
    print("\tSorted events:", colLgths, np.sum(colLgths))
    dfi = pd.read_csv(inPath, sep=";")
    print("\tOriginal events:", str(len(dfi.index)))

    if plot: plotCoinsSimple(dfc, ri=ri, coinWin_ns=coinWin_ns, titExtra="NO OFFLINE CUTS", nameExtra="000")
        
    return dfc


# In[10]:


def getChLst(df):
    ''' Returns list of channels in a sorted df'''
    chLst = (list(set([int((col.replace("E","")).split("_AU")[0]) for col in df.columns if "_AU" in col]))); chLst.sort()
    return chLst

@deprecated("Use okFlagsStr = ' | '.join(okFlags) instead")
def flagToStrg(okFlags):
    ''' Takes a list of flags you want to keep, and transforms into a string with | separators (or)'''
    okFlagsStr = okFlags[0]
    for of in okFlags[1:]: okFlagsStr += "|" + of # For cuts
    return okFlagsStr
    
def plotSortedSpectraQuality(dfl, ri, okFlags=['0x4000', '0x4100'], log=False, eRanges_AU=None, tRanges_ns=None):
    ''' Spectra of various channels, w/o coins, and w/o cuts
    '''
    chLst = getChLst(dfl['All']) # Should use info in ri instead
    #print("TEST", dfl['All'].columns)
    
    okFlagsStr = "|".join(okFlags) #flagToStrg(okFlags) #okFlags[0]
    #okFlagsStr = flagToStrg(okFlags)
    #print(okFlagsStr)
    #for of in okFlags: okFlagsStr += "|" + of # For cuts
    
    fig, ax = plt.subplots(5, len(chLst), figsize=(18, 16)) 
    fig.suptitle(ri.figTitSorted + ", " + str(len(dfl['All'].index)) + " events\nOK flags: " + str(okFlags))
    
    for (i, ch), colc in zip(enumerate(chLst), ['b', 'r', 'yellow']): # 1 column per channel spectrum
        ax[0][i].set_title('C' + str(ch))
        dfc = dfl['C' + str(ch)] # Select the channel
        eBins_AU = np.linspace(0, max(dfc['E' + str(ch) + '_AU']), 71) if eRanges_AU is None else np.linspace(eRanges_AU[i][0], eRanges_AU[i][-1], 71) # Ensure the binning is the same
        tBins_ns = 100 if tRanges_ns is None else np.linspace(tRanges_ns[i][0], tRanges_ns[i][-1], 101)
        
        for j in [0, 1]: # Lines for different cuts
            if j>0: dfc = dfl['C' + str(ch)][dfl['C' + str(ch)]['Flag' + str(ch)].str.contains(okFlagsStr)]

            ax[j][i].hist(dfc['E' + str(ch) + '_AU'], bins=eBins_AU, log=True, histtype='step', color=colc, lw=6, label='E' + str(ch) + '>0' + ('' if (j==0) else ', cut'))

            for och in [c for c in chLst if not (c == ch)]: # Overlay coins with other channels
                clo , chi = min(ch, och), max(ch, och)
                
                sm = ch+och #; print('TEST', sm)
                col = 'Purple' if (sm == 2.) else 'Green' if (sm == 4.) else 'Orange' # Ugly hack
                
                dfcc = dfl['C'+ str(clo) + '-C' + str(chi)] # Select coin
                
                if j>0: dfcc = dfl['C'+ str(clo) + '-C' + str(chi)][dfl['C'+ str(clo) + '-C' + str(chi)]['Flag' + str(ch)].str.contains(okFlagsStr) & dfl['C'+ str(clo) + '-C' + str(chi)]['Flag' + str(och)].str.contains(okFlagsStr)]
                
                ax[j][i].hist(dfcc['E' + str(ch) + '_AU'], bins=eBins_AU, log=True, histtype='step', color=col, lw=3, label='E' + str(ch) + '>0, E' + str(och) + '>0' + ('' if (j==0) else ', cut'))

                if j>0: # 2D plots of energy
                    if ([ch, och] in [[0, 2], [2, 4], [4, 0]]): # Ugly hack
                        ax[2][i].hist2d(dfcc['E' + str(ch) + '_AU'], dfcc['E' + str(och) + '_AU'], bins=[eBins_AU, 50], norm=mpl.colors.LogNorm() if log else None, cmap=col+'s', rasterized=True)
                        ax[2][i].set_xlabel('E' + str(ch) + ' (AU)'); ax[2][i].set_ylabel('E' + str(och) + ' (AU)')
                        
                        ax[3][i].hist2d(dfcc['E' + str(ch) + '_AU'], (dfcc['T' + str(ch) + '_ps']- dfcc['T' + str(och) + '_ps']) * ps_ns, bins=[eBins_AU, tBins_ns], norm=mpl.colors.LogNorm() if log else None, cmap=col+'s', rasterized=True)
                        ax[3][i].set_xlabel('E' + str(ch) + ' (AU)'); ax[3][i].set_ylabel('T' + str(ch) + ' - T' + str(och) + ' (ns)')
                        
                if j>0: # Histo of wait times
                    ax[-1][i].hist((dfcc['T' + str(ch) + '_ps']- dfcc['T' + str(och) + '_ps']) * ps_ns, bins=tBins_ns, log=True, histtype='step', color=col, lw=3, label='E' + str(ch) + '>0, E' + str(och) + '>0' + ('' if (j==0) else ', cut'))
                    ax[-1][i].set_xlabel('Time difference (ns)'); ax[-1][i].set_ylabel('Counts per bin'); ax[-1][i].legend()
                        
            ax[j][i].set_xlabel('Energy (AU)'); ax[j][i].set_ylabel('Counts per bin'); ax[j][i].legend()
        
        
        
        #ax[2][i].set_xlabel('Flag'); ax[2][i].set_ylabel('Counts'); ax[2][i].set_yscale('log'); ax[2][i].tick_params(axis='x', labelrotation = 45); ax[2][i].legend()
        
    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "plotSortedSpectraQuality loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_CoinsQual_02" + ("" if eRanges_AU is None else "Zoom") + ext)
    plt.show(); plt.clf()
            
def plotSortedDataQuality(dfl, ri, eRng=None, log=False):
    ''' 
    dfl: list of sorted data frames
    ri: runinfo
    '''
    chLst = getChLst(dfl['All'])
    #print("TEST", dfl['All'].columns)
    
    fig, ax = plt.subplots(3, len(chLst), figsize=(18, 12)) 
    fig.suptitle(ri.figTitSorted + ", " + str(len(dfl['All'].index)) + " events")

    tBins_s, tStep_s = np.linspace(0, ri.runDuration_s, 100, retstep=True)
        
    for i, ch in enumerate(chLst):
        ax[0][i].set_title('C' + str(ch))
        dfc = dfl['C' + str(ch)]
        
        ax[0][i].hist2d(dfc['T' + str(ch) + '_ps']*ps_s, dfc['E' + str(ch) + '_AU'], bins=[50, 50] if eRng is None else [50, np.linspace(eRng[0][0], eRng[0][1], 50)], norm=mpl.colors.LogNorm() if log else None, cmap='Blues', rasterized=True, label='E' + str(ch) + '>0')
        ax[0][i].set_xlabel('Time tag (s)'); ax[0][i].set_ylabel('Energy (AU)'); ax[0][i].legend()
        
        ax[1][i].hist(dfc['T' + str(ch) + '_ps']*ps_s, bins=tBins_s, weights=[1.0/tStep_s]*len(dfc['T' + str(ch) + '_ps']), log=True, histtype='step', color='b', lw=3, rasterized=True, label='E' + str(ch) + '>0')
        cnt = collections.Counter(dfc['Flag' + str(ch)])
        ax[2][i].bar(cnt.keys(), cnt.values(), ec='b', fill=False, lw=3, label='E' + str(ch) + '>0')
        
        for och, col in zip(chLst[i+1:], ['r', 'g']):
            dfcc = dfl['C'+ str(ch) + '-C' + str(och)]
            ax[1][i].hist(dfcc['T' + str(ch) + '_ps']*ps_s, bins=tBins_s, weights=[1.0/tStep_s]*len(dfcc['T' + str(ch) + '_ps']), log=True, histtype='step', color=col, lw=3, rasterized=True, label='E' + str(ch) + '>0, E'+ str(och) + '>0')
            cnt = collections.Counter(dfcc['Flag' + str(ch)])
            ax[2][i].bar(cnt.keys(), cnt.values(), ec=col, fill=False, lw=3, label='E' + str(ch) + '>0, E'+ str(och) + '>0')
        
        ax[1][i].set_xlabel('Time since start of run (s)'); ax[1][i].set_ylabel('Rate (Hz)'); ax[1][i].legend()
        ax[2][i].set_xlabel('Flag'); ax[2][i].set_ylabel('Counts'); ax[2][i].set_yscale('log'); ax[2][i].tick_params(axis='x', labelrotation = 45); ax[2][i].legend()
        
    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "plotSortedDataQuality loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_CoinsQual_00" + ext)
    plt.show(); plt.clf()

@deprecated("Replaced by plotSortedDataQuality()")    
def plotSortedDataStability(df, ri, mainCh=0, eRng=None, log=False):
    ''' Assumes 3 channels (main and 2 others)
    eRng: [[e0Lo, e0Hi], [e1Lo, e1Hi], [e2Lo, e2Hi]]
    df: sorted data frame
    ri: runinfo
    '''
    #setPars()
    
    channels = list(set(df['CHANNEL'])); channels.extend(list(set(df['COIN_CH_first'])))
    
    fig, ax = plt.subplots(3, 3, figsize=(18, 12)) 
    ax = ax.ravel()
    fig.suptitle(ri.figTitSorted + ", " + str(len(df.index)) + " events")

    for i, ch in enumerate(channels):
        ax[i].set_title('C' + str(ch))
        #if ch == mainCh: ax[i].hist(df['ENERGY'], bins=100 if eRng is None else np.linspace(eRng[0][0], eRng[0][1], 100), histtype='step')
        #else: ax[i].hist(df[df['COIN_CH_first'] == ch]['ENERGY_first'], bins=100 if eRng is None else np.linspace(eRng[i][0], eRng[i][1], 100), histtype='step')
        #ax[i].set_xlabel('Energy (AU)')
        
        if ch == mainCh: ax[i].hist2d(df['TIMETAG']*ps_s, df['ENERGY'], bins=[50, 50] if eRng is None else [50, np.linspace(eRng[0][0], eRng[0][1], 50)], norm=mpl.colors.LogNorm() if log else None, cmap='Blues', rasterized=True)
        else: ax[i].hist2d(df[df['COIN_CH_first'] == ch]['TIMETAG']*ps_s, df[df['COIN_CH_first'] == ch]['ENERGY_first'], bins=[50, 50] if eRng is None else [50, np.linspace(eRng[i][0], eRng[i][1], 50)], norm=mpl.colors.LogNorm() if log else None, cmap='Blues', rasterized=True)
        ax[i].set_xlabel('Time tag (s)'); ax[i].set_ylabel('Energy (AU)')
        
        if ch == mainCh: ax[i+len(channels)].hist2d(df['TIMETAG']*ps_s, df['ENERGY'].index.to_list(), bins=[50, 50], norm=mpl.colors.LogNorm() if log else None, cmap='Blues', rasterized=True)
        else: ax[i+len(channels)].hist2d(df[df['COIN_CH_first'] == ch]['TIMETAG']*ps_s, df[df['COIN_CH_first'] == ch]['ENERGY'].index.to_list(), bins=[50, 50], norm=mpl.colors.LogNorm() if log else None, cmap='Blues', rasterized=True)
        ax[i+len(channels)].set_xlabel('Time tag (s)'); ax[i+len(channels)].set_ylabel('INDEX')

        if ch == mainCh: 
            #ax[i+2*len(channels)].hist(df[df['COIN_CH_first']==2]['FLAGS'], histtype='step', log=True, color='r', label='Coin with 2')
            cnt = collections.Counter(df[df['COIN_CH_first']==2]['FLAGS'])
            ax[i+2*len(channels)].bar(cnt.keys(), cnt.values(), ec='r', fill=False, lw=3)
            
            #ax[i+2*len(channels)].hist(df[df['COIN_CH_first']==4]['FLAGS'], histtype='step', log=True, color='g', label='Coin with 4')
            cnt = collections.Counter(df[df['COIN_CH_first']==4]['FLAGS'])
            ax[i+2*len(channels)].bar(cnt.keys(), cnt.values(), ec='g', fill=False, lw=3)
            
            ax[i+2*len(channels)].legend()
        else:
            #ax[i+2*len(channels)].hist(df[df['COIN_CH_first']==ch]['FLAGS'], histtype='step', log=True, color='b')
            cnt = collections.Counter(df[df['COIN_CH_first']==ch]['FLAGS'])
            ax[i+2*len(channels)].bar(cnt.keys(), cnt.values(), ec='b', fill=False, lw=3)
        ax[i+2*len(channels)].set_xlabel('Flag'); ax[i+2*len(channels)].set_ylabel('Counts'); ax[i+2*len(channels)].set_yscale('log'); ax[i+2*len(channels)].tick_params(axis='x', labelrotation = 45)
        
    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "plotSortedDataStability loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_Coins_03" + ext)
    plt.show(); plt.clf()

def plotSortedData(df, ri, mainCh=0, eRng=None):
    ''' Assumes 3 channels (main and 2 others)
    eRng: [[e0Lo, e0Hi], [e1Lo, e1Hi], [e2Lo, e2Hi]]
    df: sorted data frame
    ri: runinfo
    '''
    channels = list(set(df['CHANNEL'])); channels.extend(list(set(df['COIN_CH_first'])))
    
    #setPars()
    
    fig, ax = plt.subplots(4, 3, figsize=(18, 12)) 
    ax = ax.ravel()
    fig.suptitle(ri.figTitSorted + ", " + str(len(df.index)) + " events")

    for i, ch in enumerate(channels):
        ax[i].set_title('C' + str(ch))
        if ch == mainCh: ax[i].hist(df['ENERGY'], bins=100 if eRng is None else np.linspace(eRng[0][0], eRng[0][1], 100), histtype='step')
        else: ax[i].hist(df[df['COIN_CH_first'] == ch]['ENERGY_first'], bins=100 if eRng is None else np.linspace(eRng[i][0], eRng[i][1], 100), histtype='step')
        ax[i].set_xlabel('Energy (AU)')
        
        #if ch == mainCh: ax[i+2*len(channels)].hist2d(df['TIMETAG']*ps_s, df['ENERGY'], bins=[50, 50] if eRng is None else [50, np.linspace(eRng[0][0], eRng[0][1], 50)], norm=mpl.colors.LogNorm() if log else None, cmap='Blues', rasterized=True)
        #else: ax[i+2*len(channels)].hist2d(df[df['COIN_CH_first'] == ch]['TIMETAG']*ps_s, df[df['COIN_CH_first'] == ch]['ENERGY_first'], bins=[50, 50] if eRng is None else [50, np.linspace(eRng[i][0], eRng[i][1], 50)], norm=mpl.colors.LogNorm() if log else None, cmap='Blues', rasterized=True)
        #ax[i+2*len(channels)].set_xlabel('Time tag (s)'); ax[i+len(channels)].set_ylabel('Energy (AU)')
        
        #if ch == mainCh: ax[i+3*len(channels)].hist2d(df['TIMETAG']*ps_s, df['ENERGY'].index.to_list(), bins=[50, 50], norm=mpl.colors.LogNorm() if log else None, cmap='Blues', rasterized=True)
        #else: ax[i+3*len(channels)].hist2d(df[df['COIN_CH_first'] == ch]['TIMETAG']*ps_s, df[df['COIN_CH_first'] == ch]['ENERGY'].index.to_list(), bins=[50, 50], norm=mpl.colors.LogNorm() if log else None, cmap='Blues', rasterized=True)
        #ax[i+3*len(channels)].set_xlabel('Time tag (s)'); ax[i+3*len(channels)].set_ylabel('INDEX')
        
    for i, ch in enumerate(channels[1:]):
        ax[i+len(channels)+1].hist2d(df[df['COIN_CH_first'] == ch]['ENERGY'], df[df['COIN_CH_first'] == ch]['ENERGY_first'], bins=[50, 50] if eRng is None else [np.linspace(eRng[0][0], eRng[0][1], 50), np.linspace(eRng[i+1][0], eRng[i+1][1], 50)], norm=mpl.colors.LogNorm() if log else None, cmap='Blues', rasterized=True)
        ax[i+len(channels)+1].set_xlabel('Energy C' + str(mainCh) + " (AU)"); ax[i+len(channels)+1].set_ylabel('Energy C' + str(ch) + " (AU)")

        
        ax[i+2*len(channels)+1].hist2d(df[df['COIN_CH_first'] == ch]['ENERGY'], df[df['COIN_CH_first'] == ch]['COIN_WAITS_first_ns'], bins=[50, 50] if eRng is None else [np.linspace(eRng[0][0], eRng[0][1], 50), 50], norm=mpl.colors.LogNorm() if log else None, cmap='Blues', rasterized=True)
        ax[i+2*len(channels)+1].set_xlabel('Energy C' + str(mainCh) + " (AU)"); ax[i+2*len(channels)+1].set_ylabel('Wait C' + str(ch) + '-C' + str(mainCh) + " (ns)")
        
        ax[i+3*len(channels)+1].hist(df[df['COIN_CH_first'] == ch]['COIN_WAITS_first_ns'], log=True, bins=100, histtype='step')
        ax[i+3*len(channels)+1].set_xlabel('Wait C' + str(ch) + '-C' + str(mainCh) + " (ns)"); ax[i+3*len(channels)+1].set_ylabel('Counts')
        
    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "plotSortedData() loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_Coins_04" + ext)
    plt.show(); plt.clf()

@deprecated("Replaced by getSortedData()")
def loadSortedData(ri, plot=True, removeFlags=False):
    ''' Reads a coincident sorted df from disk
    Keeps only double coins
    Optiont to removes bad flags
    Returns cleaned data frame
    ri: runinfo
    '''
    print(baseSep + str(inspect.currentframe().f_code.co_name) + "\n" + str(locals()) + "\n" + baseSep)
    
    #print(runinfo)
    df = pd.read_pickle(ri.fullSortedName)
    print("------------------ Data sorted for coins\n", df)
    
    print("Filtering out non-coins")
    df = df[df['NUM_COIN'] == 2]
    print("------------------ Data sorted for coins, non-coins removed\n", df)
    
    print("Filtering out self-coins")
    df = df[df['Self_coin'] == False]    
    print("------------------ Data sorted for coins, non-coins and self-coins removed\n", df)
    
    if removeFlags:
        print("Removing bad flags")
        df = df[df['FLAGS'] == '0x4000']; df = df[df['COIN_FLAG_first'] == '0x4000']    
        print("------------------ Bad flags removed\n", df)
    
    if plot: 
        plotSortedDataStability(df, ri, mainCh=0, eRng=[[0., 600.], [0., 3000.], [0., 2000.]], log=False)
        plotSortedData(df, ri, mainCh=0, eRng=[[0., 600.], [0., 3000.], [0., 2000.]], log=False)
        
    return df


def getSortedData(ri, plot=True):
    ''' Reads a coincident sorted df from disk
    Keeps only double coins
    Not sure how to remove bad flags: if one ch is bad, remove the others in the coin as well?
    Returns cleaned data frame
    ri: runinfo
    '''
    print(baseSep + str(inspect.currentframe().f_code.co_name) + "\n" + str(locals()) + "\n" + baseSep)
    
    df = pd.read_pickle(ri.fullCoinSortedName)
    print("------------------ Data sorted for coins\n", df)
    
    dfc = df.loc[:, (df != 0).any(axis=0)] # Get rid of empty columns
    print("------------------ Empty columns removed\n", dfc)

    #if removeFlags:
    #    print("Removing bad flags")
    #    df = df[df['FLAGS'] == '0x4000']; df = df[df['COIN_FLAG_first'] == '0x4000']    
    #    print("------------------ Bad flags removed\n", df)
    
    
    # Get useful channel numbers
    #chLst = (list(set([((col.replace("E","")).split("_AU")[0]) for col in dfc.columns if "_AU" in col]))); chLst.sort()
    chLst = getChLst(dfc)
    print(chLst)
    #shftCh = chLst[-1:] + chLst[:-1]
    
    dfl = {'All': dfc} # All events
    for ch in chLst: dfl['C' + str(ch)] = dfc[dfc['E'+str(ch) + '_AU'] > 0] # df for each specific channel in which there is an event
    for i, ch in enumerate(chLst):
        for och in chLst[i+1:]: 
            dfl['C'+str(ch) + '-C' + str(och)] = dfc[(dfc[('E' + str(ch) + '_AU')] > 0) & (dfc[('E' + str(och) + '_AU')] > 0)] # Double coins

            
    #
    #
    if plot: 
        plotSortedDataQuality(dfl, ri, eRng=None)
        plotSortedSpectraQuality(dfl, ri, eRng=None)
    #    plotSortedDataStability(df, ri, mainCh=0, eRng=[[0., 600.], [0., 3000.], [0., 2000.]], log=False)
    #    plotSortedData(df, ri, mainCh=0, eRng=[[0., 600.], [0., 3000.], [0., 2000.]], log=False)
    
    return {'Data': dfc, 'Data coins': dfl, 'Ch lst': chLst}


# In[11]:


def fullLoop(inpath=inpath, reduce=False, coinWin_ns=100., removeFlags=False):
    ''' If reduce, does the whole coincidence search, otherwise just takes the coincident file form disk
    '''
    print(baseSep + str(inspect.currentframe().f_code.co_name) + "\n" + str(locals()) + "\n" + baseSep)
    
    ri = runInfo(fullInFile=inpath)
    
    ### ID coins
    if reduce:
        df, dfCh, inpath = loadData(inpath=inpath, verbose=True, plot=True, drop=True)
        timeDiffs(df, dfCh, outpath=inpath, plot=True)
        dfc = findCoins(df, outpath=inpath, coinWin_ns=coinWin_ns, back=False, verbose=False, plot=True, cStart=0, clean=True)
        plotCoinsBasic(dfc, outpath=inpath, coinWin_ns=coinWin_ns)
        plotCoins(dfc, outpath=inpath, coinWin_ns=coinWin_ns)
        
        
    ### Work with IDed data
    loadSortedData(ri, plot=True, removeFlags=removeFlags)


# In[12]:


def compareSpectra(inpath1, inpath2, plot=True):
    ''' Compares two sets of sorted spectra, with and without flag cut
    '''
    print(baseSep + str(inspect.currentframe().f_code.co_name) + "\n" + str(locals()) + "\n" + baseSep)

    ri1 = runInfo(fullInFile=inpath1); ri2 = runInfo(fullInFile=inpath2)
    df1, df1Cut = loadSortedData(ri1, plot=False, removeFlags=False), loadSortedData(ri1, plot=False, removeFlags=True)
    df2, df2Cut = loadSortedData(ri2, plot=False, removeFlags=False), loadSortedData(ri2, plot=False, removeFlags=True)
    
    for df in [df1, df1Cut, df2, df2Cut]: print('Flags present:'); print(list(set(df['FLAGS']))); print(list(set(df['COIN_FLAG_first']))) 
    
    if plot:
        fig, ax = plt.subplots(1, 3, figsize=(22, 10)) 
        ax = ax.ravel()
        fig.suptitle(ri1.figTitSorted + "\n" + ri2.figTitSorted)

        bins0 = np.linspace(0, 2000, 70)
        bins24 = np.linspace(0, 4200, 70)
        
        for df, lab, col, lw, alp, ls, htc in zip([df1, df1Cut, df2, df2Cut], ['Filt', 'Filt, Offline Cut', 'Raw', 'Raw, Offline Cut'], ["b", "b", "r", 'r'], [3, 6, 3, 6], [1., 0.3, 1., 0.3], ['-', ':', '-', ':'], ['', '/', '', '\\']):
            if False:
                ax[0].hist(df['ENERGY'], histtype='step', bins=bins0, log=True, linewidth=lw, linestyle=ls, alpha=alp, color=col, label=lab)
                ax[0].set_xlabel('C0 Energy (AU)'); ax[0].legend()

                for i in [2, 4]:
                    ax[int(i/2)].hist(df[df['COIN_CH_first']==i]['ENERGY_first'], histtype='step', bins=bins24, log=True, linewidth=lw, linestyle=ls, alpha=alp, color=col, label=lab)
                    ax[int(i/2)].set_xlabel('C' + str(i) + ' Energy (AU)'); ax[int(i/2)].legend()
            else:
                ax[0].hist(df['ENERGY'], bins=bins0, log=True, histtype='stepfilled', linewidth=lw, linestyle=ls, hatch=htc, facecolor="none", edgecolor=col, alpha=alp, color=col, label=lab)
                ax[0].set_xlabel('C0 Energy (AU)'); ax[0].legend()

                for i in [2, 4]:
                    ax[int(i/2)].hist(df[df['COIN_CH_first']==i]['ENERGY_first'], bins=bins24, log=True, histtype='stepfilled', linewidth=lw, linestyle=ls, hatch=htc, facecolor="none", edgecolor=col, alpha=alp, color=col, label=lab)
                    ax[int(i/2)].set_xlabel('C' + str(i) + ' Energy (AU)'); ax[int(i/2)].legend()

                
        fig.tight_layout()
        if libLoad: mgs.progWaterMark(plt.gcf(), "compareSpectra() loadCAEN.ipynb")
        for ext in ['.pdf']: plt.savefig(ri1.outPath + "_Compare_00" + ext)
        plt.show(); plt.clf()
        


# In[13]:


def threeChanAn(dfl, ri, okFlags=['0x4000', '0x4100'], log=False, eRanges_AU=None, tRanges_ns=None, verbose=True):
    ''' Treats the first channel as the main one
    '''
    chLst = getChLst(dfl['All']) if False else [0, 2, 4]
    chLstNotFirst = chLst[1:]
    chPairsLst = [[0, 2], [0, 4], [2, 4]]
    
    okFlagsStr = "|".join(okFlags)
    if verbose: print("TEST flags", okFlagsStr)
    
    fig, ax = plt.subplots(3, len(chLst), figsize=(18, 16)) 
    fig.suptitle(ri.figTitSorted + ", " + str(len(dfl['All'].index)) + " events\nOK flags: " + str(okFlags))
    
    #for (i, ch), colc in zip(enumerate(chLst), ['b', 'r', 'yellow']): # 1 column per channel spectrum
    #
    #fig.tight_layout()
    #if libLoad: mgs.progWaterMark(plt.gcf(), "compareSpectra() loadCAEN.ipynb")
    #for ext in ['.pdf']: plt.savefig(ri1.outPath + "_3ch_00" + ext)
    #plt.show(); plt.clf()
    


# In[42]:


if __name__ == '__main__':
    # Std trig, FILT and RAW

    ip10 = rootPath + '/RAW/SDataR_Plastic_Scintillator_wRF_10CoarseGain_900HV_2024_05_01_sorted.csv'
    ip11 = rootPath + ''

    # Lower HV, FILT and RAW
    subFolder = ''
    ip40 = rootPath + subFolder + '/FILTERED/SDataF_Setup2_UG_RF_1000V_80d_50d_240nsCoMPASS_1000_200lsb_15h.CSV'
    ip41 = rootPath + subFolder + '/RAW/SDataR_Setup2_UG_RF_1000V_80d_50d_240nsCoMPASS_1000_200lsb_15h.CSV'
    

    cw_ns = 400.
    nMax = None
    setPars()
    
    if True:
        for i in [ip10]:
            
            ### ID coins
            if True: ri = runInfo(fullInFile=i)
            if True: df, dfCh, path = loadData(inpath=rootPath + '/RAW/SDataR_Plastic_Scintillator_wRF_10CoarseGain_900HV_2024_05_01.CSV', verbose=True, plot=True, drop=True, nrows=nMax)
                
           # if True: dfc = pd.read_csv(f'{rootPath}/SDataR_plastic_woRF_sorted.csv', sep="\t") 
            print(f'Input file pre sortCoins: {i}')
            if True: dfc = sortCoins(inPath=i, coinWin_ns=[-cw_ns, cw_ns], plot=True)    
            if True: plotCoinsSimple(dfc, ri, coinWin_ns=[-cw_ns, cw_ns])

            ### Work with IDed data
            #if False: loadSortedData(ri, plot=True) # Deprecated
            if True: res = getSortedData(ri, plot=False)
            if True: plotSortedSpectraQuality(res['Data coins'], ri, okFlags=['0x4000', '0x4100'], log=False, eRanges_AU=None, tRanges_ns=None)
            if False: plotSortedSpectraQuality(res['Data coins'], ri, okFlags=['0x4000', '0x4100'], log=False, eRanges_AU=[[0, 1000], [0, 3000], [0, 3000]], tRanges_ns=[[-100, 0], [-100, 100], [-100, 100]])


        if False:
            fullLoop(inpath=inpath, reduce=False, coinWin_ns=cw_ns, removeFlags=True)
            
    if False:
        compareSpectra(ip10, ip11, plot=True)
        
    print("That's all folks!", '\a'); os.system("printf '\a'")


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


# nnn = 4

# fig, ax = plt.subplots(nnn, nnn, figsize=(18, 18)) 
    
# for i in range(nnn):
#     #for j in range(i+1, nnn):
#     for j in range(nnn):
#         chi, chj = 2*i, 2*j
#         print(i, j, chi, chj)
        

#         ax[i][j].set_xlabel('C' + str(chi) + ' energy (AU)'); ax[i][j].set_ylabel('C' + str(chj) + ' energy (AU)')
#         ax[i][j].text(0.5, 0.5, str(i) + ', ' + str(j),
#             verticalalignment='bottom', horizontalalignment='right',
#             transform=ax[i][j].transAxes,
#             color='green', fontsize=15)

# fig.tight_layout()
# plt.show(); plt.clf()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




# In[ ]:




