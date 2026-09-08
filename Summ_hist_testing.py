import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import time 
import datetime
from matplotlib import cm
import matplotlib
from iminuit import Minuit
from iminuit.cost import ExtendedUnbinnedNLL, ExtendedBinnedNLL
from numba_stats import truncnorm, truncexpon
from scipy.stats import skewnorm
import scipy.stats as Stats
from ROOT import TCanvas, TH2D, TCutG,TProfile, TF1, kRed,TLegend,TH2F
import ROOT

ROOT.gStyle.SetLabelSize(0.05, "xyz")  # For axis labels
ROOT.gStyle.SetTitleSize(0.05, "xyz")  # For axis titles
ROOT.gStyle.SetTitleSize(0.1, "")     # For overall histogram/graph title


matplotlib.rcParams["font.size"] = 30
matplotlib.rcParams["lines.linewidth"] = 3
matplotlib.rcParams["mathtext.default"] = 'regular'
matplotlib.rcParams['lines.markersize'] = 3

class ChData:
    def __init__(self,ch,E):
        self.ch = int(ch)
        self.E = [E]
        # self.t = [t]
        
    def AddEvent(self,E):
        self.E.append(E)
        # self.t.append(t)
    
        
        
class coincData:
    def __init__(self,channel,Energy):
        self.chList = channel
        self.chData = {}
        for i,ch in enumerate(channel):
            self.chData.update({f'{ch}' : ChData(int(ch),Energy[i])})
        
        # print(self.chData['4'].ch)
        
    def AddEvent(self,channel,Energy):
        
        for i, ch in enumerate(channel):
            self.chData[f'{ch}'].AddEvent(Energy[i])

class Summed2DHist:
    def __init__(self,ch):
        self.Channels = ch
        self.cut = None 
        self.xInts = np.array([])
        self.yInts = np.array([])
        
        self.sumInts = np.array([])
        self.diffInts = np.array([])
        
    def addHistData(self,xData,yData):
        self.xInts = np.concatenate((self.xInts,xData))
        self.yInts = np.concatenate((self.yInts,yData))
        
        xData = np.asarray(xData,dtype='float64')
        yData = np.asarray(yData,dtype='float64')
        self.sumInts = np.concatenate((self.sumInts,(xData + yData)))# / np.sqrt(2)))
        self.diffInts = np.concatenate((self.diffInts,(xData - yData)))#/ np.sqrt(2)))
        
    def make2DHist(self,xRange,yRange,bins):
        self.hist2D = TH2D(f"Histogram_2D_ch_{self.Channels[0]}_ch_{self.Channels[1]}",f"Ch {self.Channels[0]} int. vs Ch {self.Channels[1]} int.",bins,xRange[0],xRange[1],bins,yRange[0],yRange[1])
        for x,y in zip(self.xInts,self.yInts):
            self.hist2D.Fill(x,y)
            
    def makeSumHist(self,xRange,bins):    
        self.sumIntsHist = ROOT.TH1D(f'Summed_int_hist_ch_{self.Channels[0]}_ch_{self.Channels[1]}', f"Ch {self.Channels[0]}, Ch {self.Channels[1]} summed integral", bins,xRange[0],xRange[1])
        
        for data in self.sumInts:
            self.sumIntsHist.Fill(data)
            
    def cutHist(self,cutEndPoints):
        self.cut = TCutG(f'Banana_cut_{self.Channels[0]}_vs_{self.Channels[1]}', len(cutEndPoints))
        for i in range(len(cutEndPoints)):
            self.cut.SetPoint(i,cutEndPoints[i][0], cutEndPoints[i][1])
        self.cut.SetLineColor(kRed+2)
        self.cut.SetLineWidth(3)
        
    def cutHistRotated(self,cutEndPoints):
        
        cutEndPointsRotated = [[point[0]-point[1],point[0]+point[1]] for point in cutEndPoints]
        
        self.cutRotated = TCutG(f'Banana_cut_{self.Channels[0]}_vs_{self.Channels[1]}_rotated', len(cutEndPoints))
        for i in range(len(cutEndPointsRotated)):
            self.cutRotated.SetPoint(i,cutEndPointsRotated[i][0], cutEndPointsRotated[i][1])
        self.cutRotated.SetLineColor(kRed+2)
        self.cutRotated.SetLineWidth(3)
            
    def make2DSumDiffHist(self,xRange,yRange,bins):
        self.sumDiff2DHist = TH2D(f'Sum_Diff_2D_hist_ch_{self.Channels[0]}_ch_{self.Channels[1]}', f"Ch {self.Channels[0]}, Ch {self.Channels[1]} sum vs diff integral", bins, xRange[0],xRange[1], bins, yRange[0], yRange[1])
        
        for sum,diff in zip(self.sumInts,self.diffInts):
            self.sumDiff2DHist.Fill(diff,sum)
            
    def rotateHistSlope(self, xRange, yRange, Bins):
        p1 = self.fit.GetParameter(1)
        
        self.slopeRotatedHist2D = TH2D(f'Slope_rot_hist_ch_{self.Channels[0]}_ch_{self.Channels[1]}',f"Ch {self.Channels[0]} Ch {self.Channels[1]} Slope Rot Hist", Bins, xRange[0],xRange[1],Bins,yRange[0],yRange[1])

        for x,y in zip(self.xInts,self.yInts):
            rotatedXData = (np.cos(np.arctan(-p1)) * x - np.sin(np.arctan(-p1))*y)
            rotatedYData = (np.sin(np.arctan(-p1)) * x + np.cos(np.arctan(-p1)) * y)
            
            self.slopeRotatedHist2D.Fill(rotatedXData,rotatedYData)
            
    def makeProfileHist(self,xRange,yRange,scaleFactor = 0.5):
        #Scale Factor is the constant that is used to cut the histogram before fitting. 

        self.profile = TProfile(f"profile_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}",f"profile Histogram Ch {self.Channels[0]} vs Ch {self.Channels[1]}", 100,xRange[0],xRange[1],yRange[0],yRange[1],"")
        
        self.cutHist2D = TH2D( f"cut_hist_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}",
        f"Cut Histogram Ch {self.Channels[0]} vs Ch {self.Channels[1]}",
        100, xRange[0], xRange[1], 100, yRange[0], yRange[1])
        
        
        maxBin = 0
        for bini in range(1, self.hist2D.GetNbinsX() + 1):
            for binj in range(1, self.hist2D.GetNbinsY() + 1): 
                x = self.hist2D.GetXaxis().GetBinCenter(bini)
                y = self.hist2D.GetYaxis().GetBinCenter(binj)
                z = self.hist2D.GetBinContent(bini,binj)
                
                if not self.cut.IsInside(x,y) and z > 0:
                    # self.profile.Fill(x,y,z)
                    if z > maxBin:
                        maxBin = z
        
        # cutThreshold = 0.5*np.sqrt(maxBin)
        cutThreshold = scaleFactor*maxBin
        
        for bini in range(1, self.hist2D.GetNbinsX() + 1):
            for binj in range(1, self.hist2D.GetNbinsY() + 1): 
                x = self.hist2D.GetXaxis().GetBinCenter(bini)
                y = self.hist2D.GetYaxis().GetBinCenter(binj)
                z = self.hist2D.GetBinContent(bini,binj)
                
                if not self.cut.IsInside(x,y) and z > cutThreshold:
                    self.profile.Fill(x,y,z)
                    self.cutHist2D.Fill(x,y,z)
        
        xBins = []
        for bini in range(1, self.cutHist2D.GetNbinsX() + 1):
            for binj in range(1, self.cutHist2D.GetNbinsY() + 1): 
                z = self.cutHist2D.GetBinContent(bini,binj)
                if z > 0:
                    xBins.append(self.cutHist2D.GetXaxis().GetBinCenter(bini))
            
        self.fitRange = [min(xBins),max(xBins)]
  
             
        self.fit = TF1(f"fit_{self.Channels[0]}_vs_{self.Channels[0]}","pol1",self.fitRange[0], self.fitRange[1])
        self.profile.Fit(self.fit, "RN")
        self.fit.SetLineColor(kRed)
        
    def makeProfileHistRotated(self,xRange,yRange,scaleFactor = 0.6):
        #Scale Factor is the constant that is used to cut the histogram before fitting. 

        self.profileRotated = TProfile(f"profile_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}_rotated",f"profile Histogram Ch {self.Channels[0]} vs Ch {self.Channels[1]}", 100,xRange[0],xRange[1],yRange[0],yRange[1],"")
        
        self.cutHist2DRotated = TH2D( f"cut_hist_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}_rotated",
        f"Cut Histogram Ch {self.Channels[0]} vs Ch {self.Channels[1]}",
        100, xRange[0], xRange[1], 100, yRange[0], yRange[1])
        
        
        maxBin = 0
        for bini in range(1, self.sumDiff2DHist.GetNbinsX() + 1):
            for binj in range(1, self.sumDiff2DHist.GetNbinsY() + 1): 
                x = self.sumDiff2DHist.GetXaxis().GetBinCenter(bini)
                y = self.sumDiff2DHist.GetYaxis().GetBinCenter(binj)
                z = self.sumDiff2DHist.GetBinContent(bini,binj)
                
                if not self.cutRotated.IsInside(x,y) and z > 0:
                    # self.profile.Fill(x,y,z)
                    if z > maxBin:
                        maxBin = z
        
        # cutThreshold = 0.5*np.sqrt(maxBin)
        cutThreshold = scaleFactor*maxBin
        
        for bini in range(1, self.sumDiff2DHist.GetNbinsX() + 1):
            for binj in range(1, self.sumDiff2DHist.GetNbinsY() + 1): 
                x = self.sumDiff2DHist.GetXaxis().GetBinCenter(bini)
                y = self.sumDiff2DHist.GetYaxis().GetBinCenter(binj)
                z = self.sumDiff2DHist.GetBinContent(bini,binj)
                
                if not self.cutRotated.IsInside(x,y) and z > cutThreshold:
                    self.profileRotated.Fill(x,y,z)
                    self.cutHist2DRotated.Fill(x,y,z)
        
        xBins = []
        for bini in range(1, self.cutHist2DRotated.GetNbinsX() + 1):
            for binj in range(1, self.cutHist2DRotated.GetNbinsY() + 1): 
                z = self.cutHist2DRotated.GetBinContent(bini,binj)
                if z > 0:
                    xBins.append(self.cutHist2DRotated.GetXaxis().GetBinCenter(bini))
            
        self.fitRange = [min(xBins),max(xBins)]
  
             
        self.fitRotated = TF1(f"fit_{self.Channels[0]}_vs_{self.Channels[0]}","pol1",self.fitRange[0], self.fitRange[1])
        self.profileRotated.Fit(self.fitRotated, "RN")
        self.fitRotated.SetLineColor(ROOT.kMagenta)
        
        # self.cutCanvasRotated = TCanvas(
        #     f"cut_canvas_ch{self.Channels[0]}_vs_ch{self.Channels[1]}",
        #     f"Cut Verification Ch {self.Channels[0]} vs Ch {self.Channels[1]}",
        #     1600, 1200
        # )
        # self.cutHist2DRotated.SetStats(0)
        # self.cutHist2DRotated.GetXaxis().SetTitle(f"Channel {self.Channels[0]} Energy (keV)")
        # self.cutHist2DRotated.GetYaxis().SetTitle(f"Channel {self.Channels[1]} Energy (keV)")
        # self.cutHist2DRotated.GetZaxis().SetTitle("Counts")
        # self.cutHist2DRotated.Draw("COLZ")
        # self.cutCanvasRotated.Update()
        
    def PlotProfileHist(self,canvas,padNumber,xRange,yRange):
        self.makeProfileHist(xRange,yRange)
        
        canvas.cd(padNumber)
        self.profile.Draw("same")
        self.fit.Draw("same")
        self.cut.Draw("same")
        canvas.Modified()
        canvas.Update()
        
    def PlotProfileHistRotated(self,canvas,padNumber,xRange,yRange):
        self.makeProfileHistRotated(xRange,yRange)
        
        canvas.cd(padNumber)
        self.profileRotated.Draw("same")
        self.fitRotated.Draw("same")
        self.cutRotated.Draw("same")
        canvas.Modified()
        canvas.Update()
        
    def PlotProfileHistRotatedSlope(self,canvas,padNumber,xRange,yRange,Bins):
        
        self.rotateHistSlope(xRange,yRange,Bins)
        
        canvas.cd(padNumber)
        
        self.slopeRotatedHist2D.Draw("COLZ")
        canvas.Update()
        
            
            
    def PlotRotatedFit(self,canvas,padNumber,fit):
        
        canvas.cd(padNumber)
        if fit == 'Rotated':
            p0 = self.fitRotated.GetParameter(0)
            p1 = self.fitRotated.GetParameter(1)
            xMin = self.fit.GetXmin()
            xMax = self.fit.GetXmax()
            
            slope = (p1*np.cos(np.arctan(p1)) + np.sin(np.arctan(p1))) / (-p1*np.sin(np.arctan(p1)) + np.cos(np.arctan(p1)))
            intercept = p0 / (-p1 * np.sin(np.arctan(p1)) + np.cos(np.arctan(p1)))
            
            self.fitRotatedTransposed = TF1(f"Transposed_rotated_fit_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}","pol1",xMin,xMax)      
            self.fitRotatedTransposed.SetParameters(intercept,slope)
            self.fitRotatedTransposed.SetLineColor(ROOT.kMagenta)
            self.fitRotatedTransposed.SetLineStyle(1) 
            
            self.fitRotatedTransposed.Draw("same")
            
            self.FitRotated = f"{intercept:.2f} + {slope:.4f}x"
            self.Fit = f"{self.fit.GetParameter(0):.2f} + {self.fit.GetParameter(1):.4f}x"

            
            
        elif fit == 'Raw':
            p0 = self.fit.GetParameter(0)
            p1 = self.fit.GetParameter(1)
            xMin = self.fitRotated.GetXmin()
            xMax = self.fitRotated.GetXmax()
            
            slope = (p1*np.cos(np.arctan(p1)) - np.sin(np.arctan(p1))) / (p1*np.sin(np.arctan(p1)) + np.cos(np.arctan(p1)))
            intercept = p0 / (p1 * np.sin(np.arctan(p1)) + np.cos(np.arctan(p1)))
            
            
            self.fitTransposed = TF1(f"Transposed_fit_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}","pol1",xMin,xMax)      
            self.fitTransposed.SetParameters(intercept,slope)
            self.fitTransposed.SetLineColor(ROOT.kRed)
            self.fitTransposed.SetLineStyle(1)   
            self.fitTransposed.Draw("same")
            
            self.Fit = f"{intercept:.2f} + {slope:.4f}x"
            self.FitRotated = f"{self.fitRotated.GetParameter(0):.2f} + {self.fitRotated.GetParameter(1):.4f}x"
            
            
        legend = TLegend(0.12, 0.65, 0.55, 0.88)     
        legend.AddEntry(self.fit, f'Raw Fit = {self.Fit}', "l")
        legend.AddEntry(self.fitRotated, f'Rotated Fit = {self.FitRotated}', "l")
        legend.SetBorderSize(1)
        legend.SetFillColorAlpha(0, 0.6)  # semi-transparent background
        legend.Draw()   
        if fit == 'Rotated':
            self.legendRotated = legend
        else:
            self.legendRaw = legend              # store to prevent garbage collection
            
        canvas.Modified()
        canvas.Update()
            
    def plot2DSumDiffHist(self,canvas,padNumber,xRange,yRange,bins,scaleData = False):
        self.make2DSumDiffHist(xRange=xRange,yRange=yRange,bins=bins)
        
        canvas.cd(padNumber)
        
        if scaleData:
            
            self.vLine = ROOT.TLine(xRange[0],662,xRange[1],662)
            self.vLine.SetLineColor(ROOT.kBlack)
            self.vLine.SetLineWidth(2)
            self.vLine.SetLineStyle(2)
            
            self.sumDiff2DHist.GetXaxis().SetTitle(f"Channel {self.Channels[0]} (keV) - Channel {self.Channels[1]} (keV)")
            self.sumDiff2DHist.GetYaxis().SetTitle(f"Channel {self.Channels[0]} (keV) + Channel {self.Channels[1]} (keV)")
        else:
            self.sumDiff2DHist.GetXaxis().SetTitle(f"Channel {self.Channels[0]} - Channel {self.Channels[1]}")
            self.sumDiff2DHist.GetYaxis().SetTitle(f"Channel {self.Channels[0]} + Channel {self.Channels[1]}")
            
            self.vLine = ROOT.TLine(xRange[0],1500,xRange[1],1500)
            self.vLine.SetLineColor(ROOT.kBlack)
            self.vLine.SetLineWidth(2)
            self.vLine.SetLineStyle(2)
        self.sumDiff2DHist.GetZaxis().SetTitle(f"Counts")
        
        ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)
        
        self.sumDiff2DHist.SetStats(0)
        self.sumDiff2DHist.Draw('COLZ')
        canvas.Update()
        self.vLine.Draw("SAME")
        canvas.Modified()
        canvas.Update()
        # canvas.Update()
        
            
    def plot2DHist(self,canvas,padNumber,xRange,yRange,bins,scaleData = False):
        
        self.make2DHist(xRange=xRange,yRange=yRange,bins=bins)
        
        canvas.cd(padNumber)
        
        if scaleData:        # stats = self.hist2D.GetListOfFunctions().FindObject("stats")
        
        # pad = 0.1
        # width = stats.GetX2NDC() - stats.GetX1NDC()
        # height = stats.GetY2NDC() - stats.GetY1NDC()
        
        # x1,x2,y1,y2 = 1-pad,1-pad-width,1-pad-height,1-pad
        
        # if stats:
        #     stats.SetX1NDC(x1)  # left edge
        #     stats.SetX2NDC(x2)  # right edge
        #     stats.SetY1NDC(y1)  # bottom edge
        #     stats.SetY2NDC(y2)  # top edge
        #     canvas.Modified()
        #     canvas.Update()
        #     # stats.SetTextSize(0.06) 
            self.hist2D.GetXaxis().SetTitle(f"Channel {self.Channels[0]} Integral (keV)")
            self.hist2D.GetYaxis().SetTitle(f"Channel {self.Channels[1]} Integral (keV)")
        else:
            self.hist2D.GetXaxis().SetTitle(f"Channel {self.Channels[0]} Integral (ADC)")
            self.hist2D.GetYaxis().SetTitle(f"Channel {self.Channels[1]} Integral (ADC)")
        self.hist2D.GetZaxis().SetTitle(f"Counts")
        self.hist2D.GetYaxis().SetLabelOffset(0.01)
        
        #KGreenREdViolet -> Nice colours, but I am not sure it has that much contrast. 
        #KOcean -> Good contrast, nice colours too.
        #KAvocado -> Very nice, and it is green
        #kBlueGreenYellow -> A darker shading of te default colours that I like.  
        ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)
        
        self.hist2D.SetStats(0)
        self.hist2D.Draw("COLZ")  
        canvas.Update()
        canvas.Modified()
        

        # stats = self.hist2D.GetListOfFunctions().FindObject("stats")
        
        # pad = 0.1
        # width = stats.GetX2NDC() - stats.GetX1NDC()
        # height = stats.GetY2NDC() - stats.GetY1NDC()
        
        # x1,x2,y1,y2 = 1-pad,1-pad-width,1-pad-height,1-pad
        
        # if stats:
        #     stats.SetX1NDC(x1)  # left edge
        #     stats.SetX2NDC(x2)  # right edge
        #     stats.SetY1NDC(y1)  # bottom edge
        #     stats.SetY2NDC(y2)  # top edge
        #     canvas.Modified()
        #     canvas.Update()
        #     # stats.SetTextSize(0.06) 
        
    def plot1DHist(self,canvas,padNumber,xRange,bins,scaleData = False):
        self.makeSumHist(xRange,bins)
        
        canvas.cd(padNumber)     
        
        self.sumIntsHist.GetXaxis().SetTitle(f"Ch {self.Channels[0]} int + Ch {self.Channels[1]} int")
        self.sumIntsHist.GetYaxis().SetTitle(f"Counts/bin")
        self.sumIntsHist.GetYaxis().SetLabelOffset(0.01)

        
        self.sumIntsHist.Draw("SAME")
        
        
class ScaleFactors:
    def __init__(self,ch):
        self.ch = ch
        self.scale = {}
        
    def AddData(self,date,scale):
        self.scale[date] = scale        
        
        
def readInData(file, scaleFactor = None):
    with open(file) as f: 
        for i in range(3): #skips the first 3 lines of the file. 
            next(f)
            
        date = str(file.parent.parent.parent.parent).split('/')[-1]
        date = date.replace('_','-')
        for i, line in enumerate(f):
            if i == 0: #Look at the third line to get the total number of channels in this coincidence. 
                numChannels = int(line.split("\n")[0].split(": ")[1])
            else:
                
                data = line.split("\n")[0].split(";")
                ch, E, t = [],[],[]
                for j in range(numChannels):
                    ch.append(int(data[6*j+1]))
                    # t.append(float(data[6*j+2])/1e3)
                    if scaleFactor is not None:
                        for s in scaleFactor:
                            if s.ch == int(data[6*j+1]):
                                scaleFac = s.scale[date]
                                break
                            else:
                                pass
                                
                        E.append(float(data[6*j+3])/scaleFac)
                        
                    else:
                        E.append(float(data[6*j+3]))
                if i == 1: #Use the first line to initialize the class objects using the channels. 
                    cData = coincData(ch,E)
                else:
                    cData.AddEvent(ch,E)
                    
    return cData

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

def BinHistograms(data, Range, numBins):
    #########################################
    #   Small function to bin histograms    #
    #########################################
    
    binRange = np.linspace(Range[0],Range[1],numBins) #Sets the range to bin over. 
    hist,binedges = np.histogram(data,binRange) #Bins the histogram. 
    
    return hist, binedges

def readInScaleFactors(filepath, channels,date):
    
    scale = []
    for ch in channels:
        scale.append(ScaleFactors(ch))
    with open(filepath) as f:
        next(f) #Skips the header
        
        for line in f:
            data = line.split(',')
            dataDate = data[0].split(' 00:00:00')[0]
            if dataDate == date: 
                for s in scale:
                    if s.ch == int(data[1]):
                        s.AddData(data[0].split(' 00:00:00')[0], float(data[3]))
            else:
                pass

    return scale




    
# rootfilePath = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/2026_06_17/2026_06_17_Daily_LSC_calibration_Cs137_coinc/RAW/coinc_sorted_500ns/') #Cs-137 Coinc Data
# # rootfilePath = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/2026_06_17/2026_06_17_Daily_LSC_calibration_bck_no_coinc_2/RAW/coinc_sorted_500ns') #Old Settings background
# LSCChannels = [4,5]
# NaIChannels = [8,10,12,14]

ScaleFactorFP = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/Results/Annulus_stability_data_average.txt')






rootfilePath = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/NaI_module_testing/2026_08_17/2026_08_17_Daily_LSC_calibration_Small_NaI_Module_testing/RAW/coinc_sorted_500ns/')
LSCChannels = [0,1]
NaIChannels = [2,3,4,5]

ScaleData = False
if ScaleData:
    scales = readInScaleFactors(ScaleFactorFP,LSCChannels + NaIChannels, '2026-06-17')

nBins = 100

pattern = re.compile(rf'_coinc_{LSCChannels[0]}_{LSCChannels[1]}_({NaIChannels[0]}|{NaIChannels[1]}|{NaIChannels[2]}|{NaIChannels[3]}).txt$')

coincFiles = sorted([
    f for f in rootfilePath.glob(f'**/*_coinc_{LSCChannels[0]}_{LSCChannels[1]}_*.txt')
    if pattern.search(f.name)
])
histogramData = []
for i in LSCChannels:
    for j in NaIChannels:
        histogramData.append(Summed2DHist([i,j]))

readInStartTime = time.time()
print('Reading in histogram data:')

for i,filePath in enumerate(coincFiles):
    date = filePath.parent.parent.parent.parent.stem
    settingsFilePath = filePath.parent.parent.parent / 'settings.xml'
    savefilepath = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/stability_figures/NaI_module_testing')
    savefilepath.mkdir(parents=True,exist_ok=True)
    if i == 0:
        Detectors = ReadInChannelNames(settingsFilePath)
    if ScaleData:
        cData = readInData(filePath,scales)
    else:
        cData = readInData(filePath)
    
    chPairs = [[cData.chList[0],cData.chList[2]], [cData.chList[1],cData.chList[2]]] #Makes a 2D array with the two channel pairs for the summed histogram data. 
    print(cData.chList)
    
    for j,hist in enumerate(histogramData):
            if hist.Channels == chPairs[0]:
                hist.addHistData(cData.chData[str(chPairs[0][0])].E,cData.chData[str(chPairs[0][1])].E)
            elif hist.Channels == chPairs[1]:
                hist.addHistData(cData.chData[str(chPairs[1][0])].E,cData.chData[str(chPairs[1][1])].E)
                

Hist2DCanvas = TCanvas("2D_Hist_Canvas","LSC Int vs NaI Int. Distribution", 2000,4000)
Hist2DCanvas.Divide(2,4)
Hist2DCanvas.SetLeftMargin(0.7)  
    
SumIntCanvas = TCanvas("Sum_int_canvas","Summed Channel Integrals",2000,4000)
SumIntCanvas.Divide(2,4)
SumIntCanvas.SetLeftMargin(0.7)  

SumDiffCanvas = TCanvas("Sum_diff_canvas", "Summed vs Difference Integrals", 2000,4000)
SumDiffCanvas.Divide(2,4)
SumDiffCanvas.SetLeftMargin(0.7)  

SlopeRotatedCanvas = TCanvas("Slope_rot_canvas","Slope Rotated Canvas", 2000,4000)
SlopeRotatedCanvas.Divide(2,4)

if ScaleData:
    plotRange2D = [[0,1000],[0,1000]]
    diffRange = [-800,800]
    plotRange1D = [0,1000]
    cutEndPoints = [[0,0], #[x coords, y coords]
                    [0,1000],
                    [20,1000],
                    [20,500],
                    [300,20],
                    [1000,20],
                    [1000,0],
                    [0,0]]
else:
    plotRange2D = [[0,2000],[0,2000]]
    diffRange = [-3000,3000]
    plotRange1D = [0,3000]
    cutEndPoints = [[0,0], #[x coords, y coords]
                    [0,4000],
                    [50,4000],
                    [50,800],
                    [800,50],
                    [4000,50],
                    [4000,0],
                    [0,0]]

for i,hist in enumerate(histogramData,start = 1):
    hist.plot2DHist(Hist2DCanvas,i,plotRange2D[0],plotRange2D[1],nBins,ScaleData)
    hist.cutHist(cutEndPoints)
    hist.PlotProfileHist(Hist2DCanvas,i,plotRange2D[0],plotRange2D[1])
    
    hist.plot2DSumDiffHist(SumDiffCanvas,i,diffRange,plotRange1D,nBins,ScaleData)
    hist.cutHistRotated(cutEndPoints)
    hist.PlotProfileHistRotated(SumDiffCanvas,i,diffRange,plotRange1D)
    
    hist.PlotRotatedFit(Hist2DCanvas,i,"Rotated")
    hist.PlotRotatedFit(SumDiffCanvas,i,"Raw")
    
    hist.PlotProfileHistRotatedSlope(SlopeRotatedCanvas,i,diffRange,plotRange1D,100)
    
    hist.plot1DHist(SumIntCanvas,i,plotRange1D,nBins,ScaleData)

if ScaleData:
    # Hist2DCanvas.SaveAs("/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/stability_figures/Testing_figures/Daily_cal_hist_2D_all_ch_scaled.png")    
    Hist2DCanvas.SaveAs(str(savefilepath / "Daily_cal_hist_2D_all_ch_scaled.png"))    
    # SumIntCanvas.SaveAs("/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/stability_figures/Testing_figures/Daily_cal_summed_int_1D_Hist_scaled.png")  
    SumIntCanvas.SaveAs(str(savefilepath / "Daily_cal_summed_int_1D_Hist_scaled.png"))  
    # SumDiffCanvas.SaveAs("/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/stability_figures/Testing_figures/Daily_cal_Sum_Diff_2D_Hist_scaled.png")
    SumDiffCanvas.SaveAs(str(savefilepath / "Daily_cal_Sum_Diff_2D_Hist_scaled.png"))    
else:
    # Hist2DCanvas.SaveAs("/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/stability_figures/Testing_figures/Daily_cal_hist_2D_all_ch.png")    
    # SumIntCanvas.SaveAs("/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/stability_figures/Testing_figures/Daily_cal_summed_int_1D_Hist.png")    
    # SumDiffCanvas.SaveAs("/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/stability_figures/Testing_figures/Daily_cal_Sum_Diff_2D_Hist.png")   
    Hist2DCanvas.SaveAs(str(savefilepath / "Daily_cal_hist_2D_all_ch.png")) 
    SumIntCanvas.SaveAs(str(savefilepath / "Daily_cal_summed_int_1D_Hist.png"))  
    SumDiffCanvas.SaveAs(str(savefilepath / "Daily_cal_Sum_Diff_2D_Hist.png"))    
    
    
    
