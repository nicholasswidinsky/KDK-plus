import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path
import statistics
import math
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import time 
import xml.etree.ElementTree as ET

startTime = time.time()

matplotlib.rcParams["font.size"] = 20
matplotlib.rcParams["lines.linewidth"] = 3
matplotlib.rcParams["mathtext.default"] = 'regular'
matplotlib.rcParams['lines.markersize'] = 3

# detectors = { #Global dictionary that converts channel number to the detector it is associated with, and a scale factor (KeV/ADC). Note the syntax for the channel #.
#     "Channel 4" : ["LS Left PMT",373.609/2104.961738528885],
#     "Channel 5" : ["LS Right PMT",373.609/1668.61105108102],
#     "Channel 8" : ["NaI 1"],
#     "Channel 10" : ["NaI 2"],
#     "Channel 12" : ["NaI 3"],
#     "Channel 14" : ["NaI 4"],
# }


# detectors = { #Global dictionary that converts channel number to the detector it is associated with, and a scale factor (KeV/ADC). Note the syntax for the channel #.
#     "Channel 10" : ["LS Left PMT",373.609/2104.961738528885],
#     "Channel 12" : ["LS Right PMT",373.609/1668.61105108102],
#     "Channel 14" : ["NaI"],
# }


class event:
    ############################################################################
    #   Class that stores all the details of an event.                         #
    #                                                                          #
    #   Events are within a single detector and contains information about the #
    #   Trigger time of the event, its energy, and any flags it may have.      #
    ############################################################################
    #   Variables:                                                             #
    #       -t (int): Trigger time for the event in ns (originally saved as ps)#
    #       - E (float): integral value of the event.                          #
    #       - flags (str): Flag given to the event by CoMPASS                  #
    #       - Ch (int): Channel number for the event                           #
    ############################################################################
    
    def __init__(self,t,E,flag,Ch,waveform = []):
        self.t = t
        self.E = E
        self.flag = flag
        self.Ch = Ch
        self.waveform = waveform
        
    
    
class Coincidence:
    #################################################################################################
    #   Class that stores a list of all coincident events. It also has a variable that is a binary  #
    #   representation of the channels involved in the coincidence. This will make it easier to     #
    #   parse the data later.                                                                       #
    #################################################################################################
    #   Variables:                                                                                  #
    #       - Events (array(event)): List of all the event objects in a given coincidence.          #
    #       - channels (array(int)): Array of 16 zeroes. A zero gets a 1 added to it if that channel#
    #         is in this set of coincidences.                                                       #
    #################################################################################################
    
    def __init__(self):
        self.Events = []
        self.Channels = [0]*16
        self.timeDiff = [[0]*16 for _ in range(16)] #Make a 16x16 array of 0's. These will get updated to be the time difference between channels n and m.
        
        
    def AddEvent(self,event,channel):
        #################################################################
        #   Function that adds a new event to this coincidence object.  #
        #   This function takes an input of an event object, and what   #
        #   channel that event is on.                                   #
        #################################################################
        self.Events.append(event)
        self.Channels[channel]+=1
        
        
    def calcTimeDiff(self):
        #################################################################
        #   Function that outputs an nxn array of the time differences  #
        #   between all the coincident events.                          #
        #################################################################
        #   Eg: for a triple coincidence:                               #
        #   Output: [[1-1,1-2,1-3],[2-1,2-2,2-3],[3-1,3-2,3-3]]         #
        #################################################################
        # print(self.timeDiff)
        
        for i in range(len(self.Events)):
            for j in range(len(self.Events)):
                self.timeDiff[self.Events[i].Ch][self.Events[j].Ch] += (self.Events[i].t - self.Events[j].t)
                

class Detector:
    ##########################################################################################
    #   Class that compiles all of the events within a specific detector. While this class   #
    #   currently only functions in coincidence mode, I do plan to make it work for          #
    #   non-coincidence data.                                                                #
    #       1. Coincidence mode:                                                             #
    #           This class is designed to work with the following class "totalCoincidence".  #
    #           That class will store a list of all the detector objects for a given         #
    #           coincidence. This class will store all of the events from a specific channel #
    #           in the coincidence. Thus for a triple coincidence between channesl x,y and z,#
    #           all of the data for channel x is stored in one "Detector object". This is    #
    #           the same for channels y and z.                                               #
    #       2. Non-coincidence mode:                                                         #
    #           This currently isn't set up but the functionality would be to store all the  #
    #           events from a given detector with no coincidence requirements.               #
    ##########################################################################################
    #   Variables:                                                                           #
    #       - events (list): List of all the events that are being compiled in this class.   #
    #       - t (list(int)): List of all the times (in ns) for each event.                   #
    #       - E (list(float))): List of all the integrals for each event.                    #
    #       - flag (list(str)): List of all the flags for each event.                        #
    #       - ch (int): Channel number that this object is for.                              #
    #       - tDiff (list(list(float))): 2D array of the time difference between the channel #
    #                                    and the other 16 channels of the digitizer.         #
    ##########################################################################################
    
    def __init__(self,Ch,dataName):
        self.events = []
        self.t = []
        self.E = []
        self.flag = []
        self.ch = Ch
        self.dataName = dataName
        self.tDiff = [[] for _ in range(16)]
        
    def AddEvent(self,event,tDiff):
        #################################################################
        #   Function that appends a new event to this class. Note that  #
        #   this class will always store events that have the same      #
        #   channel number.                                             #
        #################################################################
        
        self.events.append(event)
        self.t.append(event.t)
        self.E.append(event.E)
        self.flag.append(event.flag)
        for i in range(len(tDiff)):
            self.tDiff[i].append(tDiff[i])
            
    def EnergyHist1DPlot(self, ax, ChBins, BinRange, norm, log):
        #############################################################################
        #   This is a function that plots the energy spectra of a given event.      #
        #   The inputs of this function are the detector object itself, the axis    #
        #   that it will be plotted over, then the number of bins and the bin       #
        #   range. By default the spectra are all normalized to rate. I have no     #
        #   plans to change this right now, but maybe in the future it will be      #
        #   changed to be an option.                                                #
        #############################################################################
        #   Variables:                                                              #
        #       - self: The detector object that will be plotted.                   #
        #       - ax: The axis from plt.subplots that this spectra will be plotted  #
        #             on.                                                           #
        #       - ChBins (int): The number of bins that will be used for this       #
        #                       histogram.                                          #
        #       - BinRange (array(int)): The range that the data will be plotted    #
        #                                over.                                      #
        #       - norm (boolean): Normalizes the data to an area of 1.              #
        #       - log (boolean): Puts the y-axis to log scale.                      #
        #############################################################################
        
        Int, Bins = BinHistograms(self.E,BinRange,ChBins) #Calls the bin histogram function to bin the data. 
        ax.hist(Bins[:-1],bins =Bins,density = norm,weights=Int/(self.t[-1] - self.t[0]),histtype = 'step',label = f'{detectors[f'Channel {self.ch}'][0]} {self.dataName}', log = log,linewidth = 3) #Plot the data on the provided histogram. Note that this data is always rate normalized. 
        

        #Set various labels for the plots. 
        ax.set_title(detectors[f'Channel {self.ch}'][0])
        ax.set_xlabel('Integral (ADC Channel)')
        ax.set_ylabel('Counts/bin/ns')
        ax.legend(loc = 'best')
        
    def EnergyStabilityPlot(self, ax,fig, TimeBinRange,ChBinRange):
        
        h = ax.hist2d(self.t, self.E, bins = (TimeBinRange,ChBinRange))
        fig.colorbar(h[3],ax=ax)
        
        ax.set_title(detectors[f'Channel {self.ch}'][0])
        ax.set_xlabel('Trigger Time (ns)')
        ax.set_ylabel('Integral (ADC) Channel')
        
    def FlagStabilityPlot(self,ax):
        
        uniqueFlags = sorted(set(self.flag))
        
        flagDictionary = dict(zip(uniqueFlags, [0 for _ in range(len(self.flag))]))
        
        for i in (self.flag):
            flagDictionary[i] +=1
            
        ax.bar(range(len(flagDictionary)), list(flagDictionary.values()),align = 'center')
        ax.set_title(f'{detectors[f'Channel {self.ch}'][0]} flag distribution')
        ax.set_yscale('log')
        ax.set_xticks(range(len(flagDictionary)), list(flagDictionary.keys()),rotation = 70)
        ax.set_ylabel('Number of Events')
        

            
    
        
class totalCoinc:
    #####################################################################################
    #   This class stores a list of all coincident events between the same detectors.   #
    #   Events for this class are initialized externally, then parsed here into the     #
    #   detectors class. This class also sorts all the data in the coincidence events   #
    #   into the detector class that holds all the data for a given channel.            #
    #####################################################################################
    #   Variables                                                                       #
    #       - Coincidences (array(coincidence)): Array of coincident class objects. All #
    #         of these objects have the same coincidence channels.                      #
    #       - coincChannels (array(int)): 16-bit Binary array representing the          #
    #         coincidence channels.                                                     #
    #       - detectors (array(detectors)): Array of detectors objects. Each detector   #
    #         object contains all the data for one channel in the coincidence.          #
    #                                                                                   #
    #   Eg. coincidence between channesl 1, 2, and 4.                                   #
    #       - Every coincident event that has these (and only these) channels gets      #
    #         stored in the self.Coincidensces list.                                    #
    #       - Self.coincChannels stores the channels that are enabled for the           #
    #         coincidence.                                                              #
    #       - self.detectrs stores three detector class objects, one for detector 1,    #
    #         one for detector 2, and one for detector 4.                               #
    #####################################################################################
    
    def __init__(self,coincChannels, dataName):
        self.Coincidences = []
        self.coincChannels = coincChannels
        self.detectors = []
        self.dataName = dataName
        
    def SortChannels(self):
        #####################################################################################
        #   Function that looks through every coincident event in the list and sort them    #
        #   based on the channel that they have.                                            #
        #####################################################################################
        
        chInd = [] #Set and index to find the channel numbers that are in this coincidence. 
        for i, ch in enumerate(self.coincChannels):
            if ch != 0:
                chInd.append(i) #Decimal number for the channel -- Note that previously we have used binary, and this is now in decimal. 
        
        
        for i in chInd: #Loops over the number of detectors in the coincidence then makes 'n' blank detectors objects for 'n' channels. 
            self.detectors.append(Detector(i,self.dataName)) #Creates a blank list of detectors that correspond to the number of channels in coincidence. 
         
        #God I hope this doesn't run too slowly...
        for i, coinc in enumerate(self.Coincidences): #Loops through all of the events in the coincidence
            for k, event in enumerate(coinc.Events): #Loops through all of the events within that coincidence object. 
                for j, ch in enumerate(chInd): #Loops though every channel that is in this coincidence. 
                    if event.Ch == ch: #If the channel of the event matches the current channel, save that event.
                        # print(coinc.timeDiff[ch])
                        self.detectors[j].AddEvent(event,coinc.timeDiff[ch]) #Use the AddEvent method in the detectors class to add this event to that object. 
    
    def EnergyHist1D(self, Bins, Binrange, saveFilePath, fileName, norm, log, scale, additionalData = []):
        #################################################################
        #   This function plots the 1D spectra of the given spectra.    #
        #   Currently I have not fully implemented a system to handle   #
        #   plotting spectra from different runs at the same time.      #
        #################################################################
        #   Variables:                                                  #
        #       self: uses the current class object for this.           #
        #       additionalData: Array of other class objects to be      #
        #                       plotted alongside the "self" data.      #
        #       Bins: integer number of bins to use for each detector   #
        #       Binrange: This is the range that the bins will be       #
        #                 divided over.                                 #
        #       SaveFilePath: This is the filepath where the figures    #
        #                     will be saved.                            #
        #       fileName: This is the base name that these figures will #
        #                 have.                                         #
        #       norm: Boolean variable that turns normalization on or   #
        #             off.                                              #
        #       log: Boolean variable that turns log scale on or off.   #
        #       scale: Boolean variable that will use the scale factor  #
        #              to convert from ADC channels to KeV.             #
        #################################################################
        
        #############################################################################################################################
        #   Below is a lot of code that is designed to automatically set the number of subplots based on the number of detectors.   # 
        #   It looks confusing because it is, and I couldn't find a better way to automate the organization of up to 16 subplots.   #
        #   Sorry in advance if you are trying to read this.                                                                        #
        #############################################################################################################################
        
        root = np.sqrt(len(self.detectors)) #Take the square root of the number of detectors. If this number is a perfect square then make the subplots a n x n grid
        if root == int(root): #Check if the square root is an integer. 
            nrows = root #Set the number of rows and columns to the square root. 
            ncols = root
            
        
        else: #If the number of detectors isn't a perect square then we have to get more creative. 
            nrows = int(root) #set the number of rows to the floor of the square root. 
            ncols = len(self.detectors)//nrows #Set the number of columns equal to the number of detectors, integer divided by the number of rows. 
            if (ncols * nrows) != len(self.detectors): #If the number of rows times the number of columns doesn't equal the number of detectors, add a column. 
                ncols +=1
            
            
        nrows = int(nrows) #Ensure that the number of rows and columns is an integer and not a float (matplotlib panics if it is a float.)
        ncols = int(ncols)
        fig, ax = plt.subplots(nrows,ncols, figsize = (10*ncols,10*nrows)) #Initialize your subplot array based on the number of rows and columns that were provided. 
        fig.tight_layout() #Tight layout is always good. 
        plt.subplots_adjust(wspace = 0.2,hspace =0.3) #Adjust the spacing between the subplots. 
        
        if nrows == 1: #The way subplots works is annoying, and if there is only 1 row, then the index for ax is 1D. If there are n > 2 rows then its a 2D index. 
            
            for i,coinc in enumerate(self.detectors): #Loops through all the detectors. 
                coinc.EnergyHist1DPlot(ax = ax[i], ChBins =Bins[i], BinRange = Binrange[i], norm = norm, log =  log) #plot the 1D hist for the detector, passing the axis it is plotted on. 
                
                
                if len(additionalData) != 0: #if the length of the additional data is 0 (blank array) then plot the additional data. 
                    for j in additionalData: #Loop through all the different data series. 
                        for k, coinc2 in enumerate(j.detectors): #Loop through the detectors in the chosen totalCoinObject.
                            if coinc2.ch == coinc.ch: #Check if this detector has the same channel number as the series plotted by "self"
                                coinc2.EnergyHist1DPlot(ax = ax[i], ChBins =Bins[i], BinRange = Binrange[i], norm = norm, log =  log) #Plot the additional data on the same plot. 
                
                # # This is currently not functioning, so ignore. 
                # try:
                #     for j,coinc2 in enumerate(additionalData):
                #         if coinc2[i].ch == coinc.ch:
                #             coinc2[i].EnergyHist1DPlot(ax = ax[i], ChBins =Bins[i], BinRange = Binrange[i], norm = norm, log =  log)
                # except:
                #     pass
                
        #If there are more then 1 rows, then loop through all of them in a systematic manner. 
        else:
            row = 0
            col = 0
            
            for i, coinc in enumerate(self.detectors): #I think this is wrong, but I haven't tested it yet. 
                coinc.EnergyHist1DPlot(ax = ax[row], chBins =Bins[i], BinRange = Binrange[i], norm = norm, log =  log)
                
                col += 1
                if col == ncols:
                    col = 0
                    row +=1
         
        plt.savefig(saveFilePath / f'{fileName}_integral_1D_hist.png',bbox_inches='tight')
        
        plt.close()       
        
    def EnergyHist2D(self,ChBins,BinRange,saveFilePath,fileName,log,scale):
        #################################################################################
        #   2-dimensional histogram plotting the energy of two events against another.  #
        #   Histograms are plotted along the upper triangle in a square matrix. This is #
        #   because the plots along the y = -x line are plotting data against itself,   #
        #   and the lower triangle is just the inverse of the upper triangle plots.     #
        #################################################################################
        #   Variables:                                                                  #
        #       - self: the current totalCoinc object                                   #
        #       - ChBins (array(int)): Array for the number of bins per channel         #
        #       - BinRange (array(array(int))): 2D array containing the range overwhich #
        #                                       data will be plotted.                   #
        #       - saveFilePath (str): filepath to where it will be saved.               #
        #       - fileName (str): Name that the figure will be saved under.             #
        #       - log (bool): Boolean to set the y-axis to log scale.                   #
        #       - scale (bool): Boolean to scale the integral from ADC to KeV           #
        #################################################################################
        
        ChBinsArray = [] #Make a blank 2D array that will store the list of bins for each 2D array
        for i in range(len(self.detectors)): #Loop through each detector object to append their 
            ChBinsArray.append(np.linspace(BinRange[i][0],BinRange[i][1],ChBins[i]))
            
            
        fig,ax = plt.subplots(len(self.detectors),len(self.detectors),figsize = (10*len(self.detectors),10*len(self.detectors)))   
        plt.tight_layout() 
        plt.subplots_adjust(wspace=0.4,hspace=0.2)

        for i,det1 in enumerate(self.detectors):
            for j,det2 in enumerate(self.detectors):
                
                # if j < i: 
                #     ax[i,j].set_axis_off()
                #     continue
                
                if log:
                    h = ax[i,j].hist2d(det1.E,det2.E,bins = (ChBinsArray[i],ChBinsArray[j]),cmin = 1,norm=matplotlib.colors.LogNorm())
                else:
                    h = ax[i,j].hist2d(det1.E,det2.E,bins = (ChBinsArray[i],ChBinsArray[j]),cmin = 1)
                    
                fig.colorbar(h[3],ax=ax[i,j])
                
                if i ==0 and j ==1:
                    ax[i,j].plot(ChBinsArray[i],ChBinsArray[i], color = 'red', linestyle = '-.')
                    
                    
                ax[i,j].set_title(f'{detectors[f'Channel {det1.ch}'][0]}, {detectors[f'Channel {det2.ch}'][0]} 2D hist')
                if scale:
                    ax[i,j].set_xlabel(f'{detectors[f'Channel {det1.ch}'][0]} Integral (KeV)')
                    ax[i,j].set_ylabel(f'{detectors[f'Channel {det2.ch}'][0]} Integral (KeV)')
                else: 
                    ax[i,j].set_xlabel(f'{detectors[f'Channel {det1.ch}'][0]} Integral (ADC Channel)')
                    ax[i,j].set_ylabel(f'{detectors[f'Channel {det2.ch}'][0]} Integral (ADC Channel)')
                    
        if scale:
            plt.savefig(saveFilePath/ f'{fileName}_scaled_integral_2D_hist.png',bbox_inches='tight')
        else:
            plt.savefig(saveFilePath/ f'{fileName}_integral_2D_hist.png',bbox_inches='tight')
        plt.close()
    
    def TimeDiff1D(self,BinRange,saveFilePath,fileName,log):
                   
        timeBins = (BinRange[1] - BinRange[0])//8

        timebinsRange = np.linspace(BinRange[0],BinRange[1],timeBins)


        fig,ax = plt.subplots(len(self.detectors),len(self.detectors),figsize = (10*len(self.detectors),10*len(self.detectors)))
        
        for i,det1 in enumerate(self.detectors):
            for j,det2 in enumerate(self.detectors):
                
                ax[i,j].hist(det1.tDiff[det2.ch],bins = timebinsRange,log = log)
                ax[i,j].set_title(f'Time Difference {detectors[f'Channel {det1.ch}'][0]} - {detectors[f'Channel {det2.ch}'][0]}')
                ax[i,j].set_xlabel('Time Difference (s)')
                ax[i,j].set_ylabel('Counts')
                
        # print(saveFilePath/f'{fileName}_time_diff_1D_hist.png')
        plt.savefig(saveFilePath/f'{fileName}_time_diff_1D_hist.png',bbox_inches='tight')
        plt.close()
        
    def TimeDiff2D(self,ChBins,ChBinRange,TimeBinRange,saveFilePath,fileName,log,scale):
        #####################################################################################
        #   This function plots a 2D array of time difference between the two channels and  #
        #   the energy of one of the channels.                                              #
        #####################################################################################
        #   Variables:                                                                      #
        #       - 
        
        ChBinsArray = []
        for i in range(len(self.detectors)):
            ChBinsArray.append(np.linspace(ChBinRange[i][0],ChBinRange[i][1],ChBins[i]))
            
        timeBins = (TimeBinRange[1] - TimeBinRange[0])//2
        timebinsRange = np.linspace(TimeBinRange[0],TimeBinRange[1],timeBins)
        
        fig,ax = plt.subplots(len(self.detectors),len(self.detectors),figsize = (10*len(self.detectors),10*len(self.detectors)))
        
        for i,det1 in enumerate(self.detectors):
            for j,det2 in enumerate(self.detectors):
                
                if log:
                    h = ax[i,j].hist2d(det1.E,det1.tDiff[det2.ch],bins = (ChBinsArray[i],timebinsRange),cmin = 1,norm=matplotlib.colors.LogNorm())
                else:
                    h = ax[i,j].hist2d(det1.E,det1.tDiff[det2.ch],bins = (ChBinsArray[i],timebinsRange),cmin = 1)
                    
                
                fig.colorbar(h[3],ax=ax[i,j])
                    
                ax[i,j].set_title(f'{detectors[f'Channel {det1.ch}'][0]} Integral, {detectors[f'Channel {det1.ch}'][0]} - {detectors[f'Channel {det2.ch}'][0]} 2D hist')
                
                if scale:
                    ax[i,j].set_ylabel(f'{detectors[f'Channel {det1.ch}'][0]} - {detectors[f'Channel {det2.ch}'][0]} time difference (ns)')
                    ax[i,j].set_xlabel(f'{detectors[f'Channel {det1.ch}'][0]} Integral (KeV)')
                else: 
                    ax[i,j].set_ylabel(f'{detectors[f'Channel {det1.ch}'][0]} - {detectors[f'Channel {det2.ch}'][0]} time difference (ns)')
                    ax[i,j].set_xlabel(f'{detectors[f'Channel {det1.ch}'][0]} Integral (ADC Channel)')
   
        if scale:
            plt.savefig(saveFilePath/ f'{fileName}_scaled_integral_time_diff_2D_hist.png',bbox_inches='tight')
        else:
            plt.savefig(saveFilePath/f'{fileName}_integral_time_diff_2D_hist.png',bbox_inches='tight')
        plt.close()
                    
                    
    def Stabilityplots(self,ChBins,saveFilePath,fileName,scale):
        
        
        
        fig, ax = plt.subplots(3,len(self.detectors), figsize = (10*len(self.detectors),20))
        plt.tight_layout()
        plt.subplots_adjust(left = 0.06,wspace = 0.15,hspace = 0.25,top = 0.98,bottom = 0.04)
        
        for i,coinc in enumerate(self.detectors):
            Binrange = [0,4000]
            
            ChBinRange = np.linspace(Binrange[0],Binrange[1],ChBins[i])
            
            timebinsRange = np.linspace(min(coinc.t),max(coinc.t),100)
            
            coinc.EnergyHist1DPlot(ax = ax[0,i], ChBins =ChBins[i], BinRange = Binrange, norm = False, log = False) 
            
            coinc.EnergyStabilityPlot(ax[1,i],fig,TimeBinRange = timebinsRange, ChBinRange = ChBinRange)
            
            coinc.FlagStabilityPlot(ax[2,i])
            
        if scale:
            plt.savefig(saveFilePath/f'{fileName}_scaled_stability_plots.png',bbox_inches='tight')
        else:
            plt.savefig(saveFilePath/f'{fileName}_stability_plots.png',bbox_inches='tight')
            
# def EnergyStabilityPlot(self, ax,fig, TimeBinRange,ChBinRange):
                        
def readInFile(filepath,CoincWindow):
    #################################################################
    #   Reads in the csv file, and parses it into the coincidences. #
    #                                                               #
    #################################################################
    #   Variables:                                                  #
    #       - filepath (str): Filepath of the csv being read in.    #
    #       - CoincWindo (int): Coincidence window in ns for the    #
    #         given coincidence.                                    #
    #################################################################
    
    events = [] #creates blank arrays for all the event objects and the coincidence objects.
    coinc = []
    
    with open(filepath) as f: #Open the file
        next(f) #Skip the header

        for i,line in enumerate(f): #Read the file in line by line. This is lighter on the memory and faster overall. 
            data = line.split('\n')[0].split(';') #Take the line and split it first along the newline character, then the semi-colon delimeter. 
            events.append(event(t=int(data[2])/1e3,E=int(data[3]),flag=data[5],Ch = int(data[1]))) #Take the data from each line and create a new event object. 
            
            if i==0:  #The first time we loop through we need to grab the first event ot start the coincidences. 
                firstCoincEvent = events[0] #Set the save the event as the "firstCoincEvent" This is the event that is the first event in any arb. coincidence. 
                coinc.append(Coincidence()) #Make a new coincidence object.
            else: #Every other time check if the current event falls within the coinc window of the "firstCoincEvent". 
                dt = np.abs(firstCoincEvent.t - events[-1].t) #Calculate the time diff. 
                if dt > CoincWindow: #Check if the coinc is outisde the coinc window. 
                    firstCoincEvent = events[-1] #Set the firstCoincEvent to the latest event. 
                    coinc.append(Coincidence()) #Make a new coincidence object.

            coinc[-1].AddEvent(event = events[-1],channel = events[-1].Ch) #Add the latest event to the last coinc event in the list.
            
        #Once all the coincidences are read in, loop over them all and calculate the time differences. 
        for i in coinc:
            i.calcTimeDiff()
    return events,coinc #Return the list of events and the coincidences.       
          
def readInInit(initfilepath):
    ##########################################################
    #   This function reads in the init file to get some     #
    #   information regarding the analysis. The following is #
    #   included in the init file:                           #
    #       - data filepath                                  #
    #       - Coincidence window length (ns)                 #
    #       - Coincidence Channels (list of channels that    #
    #           are in coincidence, separated by tabs. )     #
    ##########################################################
    
    with open(initfilepath) as f: #Open the init file
        
        lines = f.readlines() #Read each line from the file

        dataFilePath = [Path(lines[0].split("\n")[0].split('\t')[1])] #Grab the data filepath from the first line 
        dataName = [lines[1].split("\n")[0].split("\t")[1]] #Grabs the name of the main data series you are plotting
        fileName = lines[0].split("\n")[0].split('\t')[1].split('/')[-1].split('.CSV')[0] #Splits the file name from the filepath
        coincEnabled = lines[2].split("\n")[0].split('\t')[1] #Boolean for enabling coincidences
        coincWindow = int(lines[3].split("\n")[0].split('\t')[1]) #Grab the coincidence window from the second line
        CoincChannelsList = lines[4].split("\n")[0].split('\t')[1:] # grab the list of all the coinc channels
        EnableAdditionalData = lines[5].split("\n")[0].split('\t')[1] #Boolean for plotting additional data on the 1D histograms. 
        ind = 6 #Start looping through the additional data. 
        while EnableAdditionalData == 'True': #If enableAdditionalData is true then read in the additional data. 
            try: #Read in lines until we get an exception then break.
                if ind == 6: #The first time we loop through the additional data will have text before it so we get rid of that. 
                    fileInd = 2
                else: 
                    fileInd = 1
                dataFilePath.append(lines[ind].split("\n")[0].split('\t')[fileInd]) #Read in the next additional data file path
                dataName.append(lines[ind].split("\n")[0].split('\t')[fileInd-1]) #Read in the next additional data name. 
                ind += 1
            except: #Break out of the wile loop when an error is encountered. 
                break
        

        
        CoincChannels = [] #Create a blank array for the final formatting of the coinc channel list. 
        for ch in CoincChannelsList: #Loop over all the different coinc channel pairs. Then separate them based into their own sub array (I might change this later but this works right now so I am keeping it. )
            CoincChannels.append(sorted(map(int,ch.split(",")))) #Append a sorted coinc channel lists to a general list. 
            
    mainData = zip(dataFilePath,dataName) #Zip the file paths and the filenames together. 
    
    coincChannelBin = [[0]*16 for _ in range(len(CoincChannels))] #Creates an array of 16 0's
    
    for i,ind in enumerate(coincChannelBin): #Create the binary list of coincidences. 
        for ch in CoincChannels[i]:
            ind[ch]+=1
            

    return list(mainData),fileName,coincEnabled,coincWindow,coincChannelBin

def ReadInChannelNames(settings):
    ##############################################################################################
    #   Function that reads in the settings.xml file produced by CoMPASS to grab the channel     #
    #   labels and make a dictionary for them. Note that I did use chatGPT to make this function #
    #   As much as it shames me, I had no idea how to read in this file so I cheated a bit.      #
    #   I have gone through the code, and the comments that ChatGPT left and it all makes sense. #
    #   I have also tested it to ensure that it is working properly.                             #
    ##############################################################################################    

    # Parse the XML file into a tree structure
    tree = ET.parse(settings)
    root = tree.getroot()

    # Initialize an empty dictionary to store channel index and label pairs
    channels = {}

    # Iterate over every <channel> element in the file
    for ch in root.iter("channel"):
        # Find the <index> tag within this <channel> element
        index_elem = ch.find("index")
        
        # Find the <values> tag, which contains parameter entries like labels
        values_elem = ch.find("values")
        
        # Continue only if both index and values are found
        if index_elem is not None and values_elem is not None:
            label = None  # Placeholder for the channel label text

            # Search through each <entry> element under <values>
            for entry in values_elem.findall("entry"):
                # Each entry has a <key> and <value>
                key_elem = entry.find("key")
                value_elem = entry.find("value")
                
                # Check if this entry is the label for the channel
                if key_elem is not None and key_elem.text == "SW_PARAMETER_CH_LABEL":
                    # If found, extract the label text (safely)
                    label = value_elem.text.strip() if value_elem is not None else None
                    break  # No need to check further entries for this channel

            # If a label was found, store it in the dictionary with its index
            if label is not None:
                channels[f'Channel {int(index_elem.text)}'] = [label]

    # Print the resulting dictionary: {channel_index: label, ...}
    return channels


def sortTotalCoincidences(coincidences,CoincCHList,name):
    #############################################################################################
    #   Sorts all of the coincidence events into the totalCoinc class, separating them based    #
    #   what channels are part of the coincidence.                                              #
    #############################################################################################
    #   Varibles:                                                                               #
    #       - coincidences (list(coincidence)): List of all coincidence objects that have been  #
    #         saved.                                                                            #
    #       - CincCHList (list(int)): Binary list of all coincidence channels that are being    #
    #         analyzed.                                                                         #
    #############################################################################################
    
    totalCoincList = [] #Make a blank array to store all the totalCoinc class objects. 
    for ch in CoincCHList: #Loop through all the coincidence channel pairs. 
        totalCoincList.append(totalCoinc(ch,name)) #Creates a new class object for the current coinc channel. 
        
    for coinc in coincidences: #Loops over all the coincident events
        for i,ch in enumerate(CoincCHList): #Loops over all the provided coincidence channel combinations. 
            if coinc.Channels == ch: #Checks to see if this event has a coincidence between the current channel coinc list. 
                totalCoincList[i].Coincidences.append(coinc) #Adds the coincident event to the totalcoincidence class. 
                
    for i in range(len(totalCoincList)): #loops through every totalCoinc object to sort them and create the detector objects. 
        totalCoincList[i].SortChannels() #Calls the SortChannels function which sorts the data into the detectors class. 
        
    return(totalCoincList)
            
def BinHistograms(data, Range, numBins):
    #########################################
    #   Small function to bin histograms    #
    #########################################
    
    binRange = np.linspace(Range[0],Range[1],numBins) #Sets the range to bin over. 
    hist,binedges = np.histogram(data,binRange) #Bins the histogram. 
    
    return hist, binedges
   
def CheckCoincChannels(coincList, coincChannels):
    #########################################################
    #   function that looks through all the coincidences    #
    #   saved by CoMPASS and compares it to the list of     #
    #   coincidences provided by the init file.             #
    #########################################################
    
    coincFound = 0 #Start a counter for whether a coincidence has been found
    for coinc in coincList: #Loop through all the channels provided by the data
        for ch in coincChannels: #Loop through all the coincidence channels provided by the init file
            if coinc == ch: #Check if they match
                coincFound +=1 #If they match, add 1 to the running tally. 
    if coincFound == 0: #If the tally is 0, return false, if it is greater than 0, then return true. 
        return False
    else: 
        return True
    
def bintodec(bin):
    #############################################
    #   Simple function that converts the       #
    #   binary representation of the channels   #
    #   to the decimal forms.                   #
    #############################################
    
    decimals = [] #make a blank array for the decimal numbers
    for i,ind in enumerate(bin): #Loop over all 16 channels in the bin array
        if ind != 0: #If the channel is not equal to 0, append the index of the channel to the decimal array
            decimals.append(i)
            
    return decimals
    




#Filepath for the CAEN_analysis_init.txt file
# initFilePath = '/home/nick/PhD/KDK+/code/CAEN_analysis_init.txt'
initFilePath = Path.cwd() / 'CAEN_analysis_init.txt'

#reads in the init file. 
mainData,fileName,coincEnabled,coincWindow,CoincChannels = readInInit(initFilePath)

#Define the filepath to the settings.xml file. 
settingsFilePath = mainData[0][0].parent.parent / 'settings.xml'
#Read in the detector names from the settings file. 
detectors = ReadInChannelNames(settingsFilePath)


#Read in the data from the csv file
events,coinc = readInFile(filepath=mainData[0][0],CoincWindow=coincWindow)



uniqueCoinc = []
for c in coinc:
    if c.Channels in uniqueCoinc:
        pass
    else:
        uniqueCoinc.append(c.Channels)
    

#Determines how long it took to read in the data. 
ReadTime = time.time()
print(f'Time to read in events: {ReadTime-startTime}')

if not coincEnabled:
    pass #Fix this later. 
else:
    if  not CheckCoincChannels(coincList=uniqueCoinc, coincChannels=CoincChannels):
        print("Coincidences provided aren't present in the data. \nThe following coincidences are present:")
        for i in uniqueCoinc:
            print(bintodec(i))
        
    else: 
        print("Coincidence's found, begin plotting:")


    #If coincidences are enebled then it will sort all the data based on the coincidence window. 
        if coincEnabled: #If true this will sort the events by the coincidence they have. 
            totalCoincList = sortTotalCoincidences(coincidences=coinc,CoincCHList=CoincChannels,name=mainData[0][1])
            
            additionalCoincidences = []
            if len(mainData) > 1: #If the list of all the data series is > 1 then read in the rest of the data as "Additional Coincidences".
                i = 1
                while i < len(mainData): #Loop through the filepaths from mainData, ignoring the first one. 
                    tempEvents,tempCoinc = readInFile(filepath=mainData[i][0],CoincWindow=coincWindow) #Reads in the additional data and sorts it based on coincidences
                    tempTotalCoinc = (sortTotalCoincidences(coincidences=tempCoinc,CoincCHList=CoincChannels,name=mainData[i][1])) #Sorts based on totalCoincidences
                    i +=1 
                    for j in tempTotalCoinc: #Loop through the data that was just parsed and append it to the additionalCoincides array. 
                        additionalCoincidences.append(j)
            
        else: #If false, then sort the events based on channel alone, no coincidence. 
            pass #Add in this functionality later. 

        #Prints the time it takes to sort the events. 
        sortTime = time.time()
        print(f'Time to sort events: {sortTime-ReadTime}')

        #A few miscellaneous settings for plotitng. 
        norm = False #Normalizes the histograms to have an area of 1. 
        log = False #Set the y-axis to log scaled. 
        scale = False #Scales the values for energy form ADC to KeV. 



        # bckfileName = bckfile.split('/')[-1].split('.CSV')[0]


        integralBins = 100

        timeBinRange = [-500,500]
        integralBinRange = [[0,4000],[0,4000],[0,4000]]
        data2 = []


        for i,coinc in enumerate(totalCoincList):
            chInd = [] #Set and index to find the channel numbers that are in this coincidence. 
            for j, ch in enumerate(coinc.coincChannels):
                if ch != 0:
                    chInd.append(j)
            Channels = ",".join(str(x) for x in chInd)

            # for j, coinc2 in enumerate(totalCoincList2):
            #     if coinc.coincChannels == coinc2.coincChannels:
            #         data2.append(coinc2)
            
            saveFilePath = mainData[0][0].with_suffix('') / 'figures' / f"Coincidence_Channels_{Channels}"
            Path(f"{saveFilePath}").mkdir(parents=True, exist_ok=True)
            plottingStartTime = time.time()
            
            coinc.EnergyHist1D(additionalData = additionalCoincidences, Bins = [integralBins,integralBins,integralBins], Binrange = integralBinRange, saveFilePath = saveFilePath, fileName = fileName, norm = norm, log = log, scale = scale)
            
            E1DTime = time.time()
            print(f'\t Time to Plot 1D Energy Hist: {E1DTime- plottingStartTime}')
            
            
            coinc.EnergyHist2D(ChBins = [integralBins,integralBins,integralBins], BinRange = integralBinRange,saveFilePath = saveFilePath, fileName = fileName, log = log, scale = scale)
            
            E2DTime = time.time()
            print(f'\t Time to Plot 2D Energy Hist: {E2DTime- E1DTime}')
            
            coinc.TimeDiff1D(BinRange = timeBinRange, saveFilePath = saveFilePath, fileName = fileName, log = True)
            
            T1DTime = time.time()
            print(f'\t Time to Plot 1D time Hist: {T1DTime- E2DTime}')
            
            coinc.TimeDiff2D(ChBins = [integralBins,integralBins,integralBins],ChBinRange = integralBinRange, TimeBinRange = timeBinRange,saveFilePath=saveFilePath,fileName=fileName,log=log,scale=scale)
            
            
            T2DTime = time.time()
            print(f'\t Time to Plot 2D Time Hist: {T2DTime- T1DTime}')
            
            
            coinc.Stabilityplots(ChBins = [integralBins,integralBins,integralBins],saveFilePath = saveFilePath, fileName = fileName, scale = scale)
            
            stabTiming = time.time()
            print(f'\t Time to Plot stability plots: {stabTiming- T2DTime}')
            
            coincTiming = time.time()
            print(f'Time to plot {i+1} coincidence: {coincTiming-plottingStartTime}')
            

print(f'Total run time: {time.time()-startTime}')

