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
# findCoins
# loadSortedData
#
#
#


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

import time

import collections # For bar charts of text data

from scipy.stats import poisson, norm, truncnorm

import pandas as pd

from uncertainties import ufloat

from num2tex import num2tex

import xml.etree.ElementTree as ET # Parsing xml

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


# In[5]:


rootPath = '/Users/arnaudlemaire/ECL/Queens/Experiment/LS/ComptonNewDAQ/setup_2/80_50d'
#inpath = rootPath + '/2SiPM_coincidence_90d_50d_1day_1/FILTERED/SDataF_2SiPM_coincidence_90d_50d_1day_1.CSV'


inpath = rootPath + '/Setup2_UG_RF_1100V_80d_50d_240nsCoMPASS_1400_200lsb_15h/RAW/SDataR_Setup2_UG_RF_1100V_80d_50d_240nsCoMPASS_1400_200lsb_15h.CSV'
#inpath = rootPath + '/2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1200_200lsb_15h_raw/FILTERED/SDataF_2SiPM_coincidence_80d_50d_240nsBoard_240nsCoMPASS_CoincidenceWindow_1200_200lsb_15h_raw.CSV'

#inpath = rootPath + '/2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1400lsb_100ns_25pc_cfd_15h_raw/FILTERED/SDataR_2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1400lsb_100ns_25pc_cfd_15h_raw.CSV'
#inpath = rootPath + '/2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1400lsb_100ns_25pc_cfd_15h_raw/FILTERED/SDataF_2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1400lsb_100ns_25pc_cfd_15h_raw.CSV'

ps_ns = 1.E-3
ps_s = 1.E-12


# In[6]:


class runInfo(object):
    def __init__(self, fullInFile=inpath, outpath=None):
        ''' A class that wrangles paths and file names, creates output directory if necessary.
        Ex: dir1/dir2/file.CSV
        '''
        self.fullRoot = os.path.splitext(fullInFile)[0] # dir1/dir2/file
        self.directory, self.fileName = os.path.split(fullInFile) # dir1/dir2 file.CSV
        
        self.outPath = self.directory if outpath is None else outpath
        if not os.path.exists(self.outPath): os.makedirs(self.outPath); print("Creating output directory " + str(self.outPath))
            
        self.figTit = self.directory + "\n" + self.fileName
        self.baseName = os.path.splitext(self.fileName)[0] # file
        
        
        # Sorted data
        self.sortedName = self.baseName + "_sorted" + ".pd" # file_sorted.pd
        self.fullSortedName = os.path.join(self.directory, self.sortedName) # dir1/dir2/file_sorted.pd
        self.figTitSorted = self.directory + "\n" + self.sortedName

        #self.sortedCleanName = self.baseName + "_sortedClean" + ".pd" # file_sortedClean.pd
        #self.fullSortedCleanName = os.path.join(self.directory, self.sortedCleanName) # dir1/dir2/file_sortedClean.pd
        
        
        
        # Setup xmlfile
        self.xmlName = "settings.xml"
        
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

    
if False: ri = runInfo()       


# In[7]:


def stabilityPlots(df, path):
    ''' Basic plots, including energy distribution, stability, and flags
    df is the dataframe read from DAQ disk output
    path is the path where the df is
    '''
    ri = runInfo(path)

    channels = set(df['CHANNEL'])
    for ch in channels: dfCh = {'C' + str(ch):  df[df['CHANNEL']==ch] for ch in channels}
    
    fig, ax = plt.subplots(4, len(dfCh), figsize=(7. * len(dfCh), 12)) 
    ax = ax.ravel()
    fig.suptitle(ri.directory + "\n" + str(ri.fileName) + ", " + str(len(df.index)) + " events, NO OFFLINE CUTS")

    for i, ch in enumerate(channels):
        ax[i].set_title('C' + str(ch) + ', ' + str(len(dfCh['C' + str(ch)].index)) + ' events')
        ax[i].hist(dfCh['C' + str(ch)]['ENERGY'], bins=100, histtype='step'); ax[i].set_xlabel('ENERGY (AU)')
        ax[i+len(dfCh)].hist2d(dfCh['C' + str(ch)]['TIMETAG'], dfCh['C' + str(ch)]['ENERGY'], bins=[50, 50]); ax[i+len(dfCh)].set_xlabel('TIMETAG (ps)'); ax[i+len(dfCh)].set_ylabel('ENERGY (AU)')
        ax[i+2*len(dfCh)].hist2d(dfCh['C' + str(ch)]['TIMETAG'], dfCh['C' + str(ch)]['ENERGY'].index.to_list(), bins=[50, 50]); ax[i+2*len(dfCh)].set_xlabel('TIMETAG (ps)'); ax[i+2*len(dfCh)].set_ylabel('INDEX')
        cnt = collections.Counter(dfCh['C' + str(ch)]['FLAGS'])
        ax[i+3*len(dfCh)].bar(cnt.keys(), cnt.values(), color='b', lw=3)
        ax[i+3*len(dfCh)].set_xlabel('FLAGS'); ax[i+3*len(dfCh)].set_yscale('log'); ax[i+3*len(dfCh)].set_ylabel('Counts'); ax[i+3*len(dfCh)].tick_params(axis='x', labelrotation = 45)

    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "stabilityPlots(), loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_Stability_00" + ext)
    plt.show(); plt.clf()
        
        
def flagPlots(df, path, log=True):
    ''' Effect of flags on spectra of various channels
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
        axx.hist(ddf['ENERGY'], bins=100, log=log, histtype='step', lw=6, ls=':', color='k', label='All')

        for fl, col in zip(flags, ['r', 'b', 'g', 'orange', 'm', 'c', 'pink', 'gray', 'yellow', 'maroon', 'violet', 'teal', 'linen']):
            axx.hist(ddf[ddf['FLAGS']==fl]['ENERGY'], bins=100, log=log, histtype='step', lw=3, color=col, label=fl)
        axx.set_title('C' + str(ch) + ', ' + str(len(ddf.index)) + ' events'); axx.set_xlabel('ENERGY (AU)'); axx.legend()
    cnt = collections.Counter(df['FLAGS'])
    ax[0].bar(cnt.keys(), cnt.values(), color='k')
    ax[0].set_xlabel('FLAGS'); ax[0].set_yscale('log'); ax[0].set_ylabel('Counts'); ax[0].tick_params(axis='x', labelrotation = 45); ax[0].set_title('All ch')
    
    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "flagPlots(), loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_Flags_00" + ext)
    plt.show(); plt.clf()
    
    
    
def loadData(inpath=inpath, verbose=False, plot=True, drop=True):
    ''' Loads CAEN data from .CSV to a pandas dataframe
    Optionnally makes data quality plots
    Does not apply any cuts, though has the optio to drop certain columns
    inpath: path to dataframe
    drop: get rid of colums BOARD, ENERGYSHORT
    Returns total dataframe, as well as df for each channel
    '''
    ri = runInfo(inpath)
    
    df = pd.read_csv(inpath, sep=";")
    if verbose: print("------------------ All data\n", df)
    
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
        if libLoad: mgs.progWaterMark(plt.gcf(), "loadCAEN.ipynb")
        for ext in ['.pdf']: plt.savefig(ri.outPath + "_TEST_02" + ext)
        plt.show(); plt.clf()


# In[9]:


def plotCoinsBasic(df, outpath=inpath, coinWin_ns=100.):
    ''' Plots information about coincidences without the reordered file
    '''
    ri = runInfo(outpath)
    
    setPars()
    
    fig, ax = plt.subplots(1, 4, figsize=(12, 6)) 
    fig.suptitle(ri.figTit)
    
    ax[0].hist(df['NUM_COIN'], log=True, histtype='step')
    ax[0].set_xlabel('Multiplicity of coins')
    
    ax[1].hist2d(df['CHANNEL'], df['NUM_COIN'], bins=[np.linspace(-0.5, 3.5, 5), np.linspace(-0.5, 3.5, 5)], norm=matplotlib.colors.LogNorm() if False else None, cmap='Greys', rasterized=True)
    ax[1].set_xlabel('Ch'); ax[1].set_ylabel('Multiplicity of coins')
    
    ax[2].hist(df['Self_coin'], bins=np.linspace(-0.5, 1.5, 3), histtype='step', color='b', log=True)
    ax[2].set_xlabel('Self coincidence'); ax[2].set_ylabel('Counts'); ax[2].set_ylim(bottom=0.1)
    
    ax[3].hist(df['COIN_WAITS_first_ns'], log=True, histtype='step') #; ax[1].axvline(coinWin_ns, color='g'); ax[1].axvline(-coinWin_ns, color='g')
    ax[3].axvspan(0, coinWin_ns, alpha=0.1, color='k')
    ax[3].set_xlabel('Wait time to first coin (ns)')
    
    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_TEST_04" + ext)
    plt.show(); plt.clf()
    
    
def plotCoins(dfc, outpath=inpath, coinWin_ns=100.):
    ''' To display ordered data
    0 is the main channel
    '''
    ri = runInfo(outpath)
    
    setPars()
    
    fig, ax = plt.subplots(1, 4, figsize=(18, 6)) 
    fig.suptitle(ri.figTit)
    
    ax[0].hist(dfc['ENERGY'], histtype='step', log=True, lw=3, color='b', label='All')
    for i, col in zip([2, 3], ['r', 'g']): ax[0].hist(dfc[dfc['COIN_CH_first']==i]['ENERGY'], histtype='step', log=True, lw=3, color=col, label='Coin C'+str(i))
    ax[0].legend(); ax[0].set_xlabel('C0 ENERGY'); ax[0].set_ylabel('Counts')
    
    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_TEST_06" + ext)
    plt.show(); plt.clf()
    
    
def findCoins(df, outpath=inpath, coinWin_ns=100., back=True, verbose=False, plot=False, cStart=None, clean=False):
    ''' Finds coins
    Returns a df with extra columns which are the channels coincident to the original one of the line
    Coins are counted as many times as their multiplicity
    back: if should look back for coins as well; otherwise, just lookds forward
    cStart: if should use a given channel like a TDC start
    Returns an ordered file, written to disk as _sorted.pd.  Events where no coins are found are kept, but listed with NaNs.
    '''
    ri = runInfo(outpath)
    
    dfc = df.copy()
    def getEvt(i): return {'i': dfc.index[i], 'Ch': dfc['CHANNEL'][i], 'E': dfc['ENERGY'][i], 't_ps': dfc['TIMETAG'][i], 'flag': dfc['FLAGS'][i]}
    
    coinChAll, coinChTimeAll_ps, coinChDiffAll_ns, coinEAll, coinFlagAll, coinCh0, coinChTime0_ps, coinChDiff0_ns, coinE0, coinFlag0, numCoin = [], [], [], [], [], [], [], [], [], [], [] # Coin ch, time, diff for all coins with given evt, for all events
    tagged = [0]*len(dfc.index); selfCoin = [False]*len(dfc.index)
    #print(tagged)
    
    if verbose: print("------------------ Before looking for coins\n", dfc)
    start_time = time.time() # Timer
    for i in dfc.index:
        if ((i <= 30) & (i % 10 == 0)) or ((i <= 300) & (i % 100 == 0)) or ((i <= 3000) & (i % 1000 == 0)) or ((i <= 30000) & (i % 10000 == 0)) or ((i % 100000 == 0)): print("Sorting " + str(i) + "/" + str(len(dfc.index)) + str(" in %s seconds." % round(time.time() - start_time, 2)))
            
        
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
        plotCoins(dfc, outpath=outpath, coinWin_ns=coinWin_ns)
        plotCoinsBasic(dfc, outpath=outpath, coinWin_ns=coinWin_ns)
        stabilityPlots(df=dfc, path=outpath)
        
    if clean:
        dfc.drop('COIN_CH', inplace=True, axis=1); dfc.drop('COIN_TIMES_ps', inplace=True, axis=1); dfc.drop('COIN_WAITS_ns', inplace=True, axis=1); dfc.drop('COIN_FLAGS', inplace=True, axis=1)
        if verbose: print("Dropping cols COIN_CH, COIN_TIMES_ps, COIN_FLAGS and COIN_WAITS_ns"); print("------------------ Dropping cols\n", df)
    
    dfc.to_pickle(ri.fullSortedName)
    return dfc
        
        


# In[10]:


def plotSortedDataStability(df, ri, mainCh=0, eRng=None, log=False):
    ''' Assumes 3 channels (main and 2 others)
    eRng: [[e0Lo, e0Hi], [e1Lo, e1Hi], [e2Lo, e2Hi]]
    df: sorted data frame
    ri: runinfo
    '''
    setPars()
    
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
            ax[i+2*len(channels)].hist(df[df['COIN_CH_first']==2]['FLAGS'], histtype='step', log=True, color='r', label='Coin with 2')
            ax[i+2*len(channels)].hist(df[df['COIN_CH_first']==4]['FLAGS'], histtype='step', log=True, color='g', label='Coin with 4')
            ax[i+2*len(channels)].legend()
        else:
            ax[i+2*len(channels)].hist(df[df['COIN_CH_first']==ch]['FLAGS'], histtype='step', log=True, color='b')
        ax[i+2*len(channels)].set_xlabel('Flag'); ax[i+2*len(channels)].set_ylabel('Counts')
        #else: ax[i].hist2d(df[df['COIN_CH_first'] == ch]['TIMETAG']*ps_s, df[df['COIN_CH_first'] == ch]['ENERGY_first'], bins=[50, 50] if eRng is None else [50, np.linspace(eRng[i][0], eRng[i][1], 50)], norm=mpl.colors.LogNorm() if log else None, cmap='Blues', rasterized=True)
        #ax[i].set_xlabel('Time tag (s)'); ax[i].set_ylabel('Energy (AU)')
        #ax[i+3*len(dfCh)].hist(dfCh['C' + str(ch)]['FLAGS'], log=True, histtype='step'); ax[i+3*len(dfCh)].set_xlabel('FLAGS'); ax[i+3*len(dfCh)].set_ylabel('Counts')
        
        
    #for i, ch in enumerate(channels[1:]):
    #    ax[i+len(channels)+1].hist2d(df[df['COIN_CH_first'] == ch]['ENERGY'], df[df['COIN_CH_first'] == ch]['ENERGY_first'], bins=[50, 50] if eRng is None else [np.linspace(eRng[0][0], eRng[0][1], 50), np.linspace(eRng[i+1][0], eRng[i+1][1], 50)], norm=mpl.colors.LogNorm() if log else None, cmap='Blues', rasterized=True)
    #    ax[i+len(channels)+1].set_xlabel('Energy C' + str(mainCh) + " (AU)"); ax[i+len(channels)+1].set_ylabel('Energy C' + str(ch) + " (AU)")
#
    #    
    #    ax[i+4*len(channels)+1].hist2d(df[df['COIN_CH_first'] == ch]['ENERGY'], df[df['COIN_CH_first'] == ch]['COIN_WAITS_first_ns'], bins=[50, 50] if eRng is None else [np.linspace(eRng[0][0], eRng[0][1], 50), 50], norm=mpl.colors.LogNorm() if log else None, cmap='Blues', rasterized=True)
    #    ax[i+4*len(channels)+1].set_xlabel('Energy C' + str(mainCh) + " (AU)"); ax[i+4*len(channels)+1].set_ylabel('Wait C' + str(ch) + " (ns)")
        
    fig.tight_layout()
    if libLoad: mgs.progWaterMark(plt.gcf(), "plotSortedDataStability, loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_TEST_10" + ext)
    plt.show(); plt.clf()

def plotSortedData(df, ri, mainCh=0, eRng=None, log=False):
    ''' Assumes 3 channels (main and 2 others)
    eRng: [[e0Lo, e0Hi], [e1Lo, e1Hi], [e2Lo, e2Hi]]
    df: sorted data frame
    ri: runinfo
    '''
    channels = list(set(df['CHANNEL'])); channels.extend(list(set(df['COIN_CH_first'])))
    
    setPars()
    
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
    if libLoad: mgs.progWaterMark(plt.gcf(), "plotSortedData(), loadCAEN.ipynb")
    for ext in ['.pdf']: plt.savefig(ri.outPath + "_TEST_08" + ext)
    plt.show(); plt.clf()


def loadSortedData(ri, plot=True, removeFlags=False):
    '''
    ri: runinfo
    '''
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


# In[11]:


def fullLoop(inpath=inpath, reduce=False, coinWin_ns=100., removeFlags=False):
    ''' If reduce, does the whole coincidence search, otherwise just takes the coincident file form disk
    '''
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
    ''' Compares two set of sorted spectra, with and without flag cut
    '''
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
        
        for df, lab, col, lw, alp, ls, htc in zip([df1, df1Cut, df2, df2Cut], ['Filt', 'Filt Cut', 'Raw', 'Raw Cut'], ["b", "b", "r", 'r'], [3, 6, 3, 6], [1., 0.3, 1., 0.3], ['-', ':', '-', ':'], ['', '/', '', '\\']):
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
        if libLoad: mgs.progWaterMark(plt.gcf(), "compareSpectra(), loadCAEN.ipynb")
        for ext in ['.pdf']: plt.savefig(ri1.outPath + "_COMPARE_00" + ext)
        plt.show(); plt.clf()
        


# In[13]:


if __name__ == '__main__':
    # Std trig, FILT and RAW
    ip10 = rootPath + '/Setup2_UG_RF_1100V_80d_50d_240nsCoMPASS_1400_200lsb_15h/FILTERED/SDataF_Setup2_UG_RF_1100V_80d_50d_240nsCoMPASS_1400_200lsb_15h.CSV'
    ip11 = rootPath + '/Setup2_UG_RF_1100V_80d_50d_240nsCoMPASS_1400_200lsb_15h/RAW/SDataR_Setup2_UG_RF_1100V_80d_50d_240nsCoMPASS_1400_200lsb_15h.CSV'

    # CFD, FILT and RAW
    ip20 = rootPath + '/2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1400lsb_100ns_25pc_cfd_15h_raw/FILTERED/SDataF_2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1400lsb_100ns_25pc_cfd_15h_raw.CSV'
    ip21 = rootPath + '/2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1400lsb_100ns_25pc_cfd_15h_raw/RAW/SDataR_2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1400lsb_100ns_25pc_cfd_15h_raw.CSV'

    cw_ns = 300.

    setPars()
    
    for inpath in [ip10, ip11]:
        if True:
            ### ID coins
            if True: ri = runInfo(fullInFile=inpath)
            if True: df, dfCh, inpath = loadData(inpath=inpath, verbose=True, plot=True, drop=True)
            if False: timeDiffs(df, dfCh, outpath=inpath, plot=True)
            if False: dfc = findCoins(df, outpath=inpath, coinWin_ns=cw_ns, back=True, verbose=False, plot=True, cStart=0, clean=True)
            if False: plotCoinsBasic(dfc, outpath=inpath, coinWin_ns=cw_ns)
            if False: plotCoins(dfc, outpath=inpath, coinWin_ns=cw_ns)


            ### Work with IDed data
            if False: loadSortedData(ri, plot=True)

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





# In[ ]:


os.system("printf '\a'")


# In[ ]:


import matplotlib.pyplot as plt

dat = ['a', 'bb', 'a', 'c', 'bb', 'a', 'c', 'bb', 'a', 'c', 'c', 'c', 'c']
res = plt.hist(dat)
print(res)
plt.show(); plt.clf()

import collections
cnt = collections.Counter(dat)
plt.bar(cnt.keys(), cnt.values(), color='k')
print(cnt)
# np.histogram(dat) # UFuncTypeError error


# In[ ]:




