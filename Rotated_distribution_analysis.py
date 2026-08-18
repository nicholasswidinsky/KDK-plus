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
    def __init__(self,ch,E,t):
        self.ch = int(ch)
        self.E = [E]
        self.t = [t]
        
    def AddEvent(self,E,t):
        self.E.append(E)
        self.t.append(t)
    
        
        
class coincData:
    def __init__(self,channel,Energy,time):
        self.chList = channel
        self.chData = {}
        
        for i,ch in enumerate(channel):
            self.chData.update({f'{ch}' : ChData(int(ch),Energy[i],time[i])})
        
        # print(self.chData['4'].ch)
        
    def AddEvent(self,channel,Energy,time):
        
        for i, ch in enumerate(channel):
            self.chData[f'{ch}'].AddEvent(Energy[i],time[i])
            
            
class Summed2DHist:
    def __init__(self,ch):
        self.dates = []
        self.Channels = ch #List of the channels that are being summed (LSC is always along the x-axis. NaI are always along the y-axis)
        self.xInts = np.array([]) #Integrals along the x axis
        self.yInts = np.array([]) #Integrals along the y axis
        self.time = [] #Time for all the events. We are able to 
        
        self.SeparationPoints = []
        self.sumChannels = np.array([])
        
        self.slope = []
        self.slopeErr = []
        
        self.rotatedXData = None
        self.rotatedYData = None
        
        
    # def addSlopeData(self,slope):
    #     #The slopes for now are read in with the inverte histogram, so we need to invert them back
        
    #     self.slope = 1/slope
        
    def addHistData(self,xData,yData,time, date):
        runTime = time[-1] - time[0]
                

        xData = np.asarray(xData,dtype='float64')
        yData = np.asarray(yData,dtype='float64')
        
        self.xInts = np.concatenate((self.xInts,xData)) #Adds the latest data to the x and y arrays
        self.yInts = np.concatenate((self.yInts,yData))
        self.time.append(runTime)
        

        
        if date in self.dates:
            self.SeparationPoints[-1] + len(self.xInts)
        else:
            self.dates.append(date)
            self.SeparationPoints.append(len(self.xInts)) #Adds the separation point to note where each day stops. 
    
    def rotatedHistData(self,scaleData,slopeDict):
        self.rotatedXData = []
        self.rotatedYData = []
        
        
        index = 0
        if scaleData:
            for i,(x,y) in enumerate(zip(self.xInts,self.yInts)):
                xData = x-y
                yData = x+y
                self.rotatedXData.append(xData)
                self.rotatedYData.append(yData)
        else:        
            for i,(x,y) in enumerate(zip(self.xInts,self.yInts)):

                if i > self.SeparationPoints[index]:
                    index +=1
                    
                # slope = self.slope[index]
                slope = slopeDict[f'{self.Channels[0]},{self.Channels[1]}'][self.dates[0]][0]
                xData = np.cos(np.arctan(-slope)) * x - np.sin(np.arctan(-slope)) * y
                yData = np.sin(np.arctan( -slope)) *x + np.cos(np.arctan(-slope)) * y
                
                self.rotatedXData.append(xData)
                self.rotatedYData.append(yData)
            
            
    def plot2DHist(self,saveFilePath,scaleData,bins):
        # bins = 100
        
        if scaleData:
            binRangeX = [0,750]
            binRangeY = [0,750]
        else:
            binRangeX = [0,2250]
            binRangeY = [0,1750]
            
        self.Hist2DCanvas = ROOT.TCanvas(f"2D_Hist_Canvas_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}",f"2D Hist Canvas Ch {self.Channels[0]} vs ch {self.Channels[1]}",2000,4000)
        
        
        self.Hist2DCanvas.Divide(6,6)

        self.hist2D = []

        
        ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)
        
        # for i in range(36): #Loop over all 36 subplots.
        #     self.Hist2DCanvas.cd(i+1) #Switch the current pad to the next canvas. 
            
        #     self.hist2D.append(TH2D(f"Histogram_2D_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}_iteration_{i}",f"Histogram 2D Ch {self.Channels[0]} vs ch {self.Channels[1]} iteration {i}",bins,binRangeX[0],binRangeX[1],bins,binRangeY[0],binRangeY[1]))
            

        #     for x,y in zip(self.xInts[0:self.SeparationPoints[i]],self.yInts[0:self.SeparationPoints[i]]):
        #         self.hist2D[-1].Fill(x,y)
                

                
        #     # Hist2D.GetXaxis().SetTitle(f"")    
        #     self.hist2D[-1].SetStats(0)
        #     self.hist2D[-1].Draw('COLZ')
        #     self.Hist2DCanvas.Update()

        for i in range(36):
            self.Hist2DCanvas.cd(i+1)
            
            # Build THIS DAY's histogram only (not cumulative slicing)
            dayHist = TH2D(f"Day_Hist_2D_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}_iteration_{i}",
                            f"Day Hist 2D Ch {self.Channels[0]} vs ch {self.Channels[1]} iteration {i}",
                            bins, binRangeX[0], binRangeX[1], bins, binRangeY[0], binRangeY[1])
            dayHist.Sumw2()
            
            if i == 0:
                xSlice = self.xInts[0:self.SeparationPoints[i]]
                ySlice = self.yInts[0:self.SeparationPoints[i]]
            else:
                xSlice = self.xInts[self.SeparationPoints[i-1]:self.SeparationPoints[i]]
                ySlice = self.yInts[self.SeparationPoints[i-1]:self.SeparationPoints[i]]
            
            for x, y in zip(xSlice, ySlice):
                dayHist.Fill(x, y)
            
            dayHist.Scale(1.0 / self.time[i])   # normalize this day's counts by its own run time
            
            # Now build the CUMULATIVE histogram by adding this scaled day onto the running total
            if i == 0:
                self.hist2D.append(dayHist)
            else:
                cumHist = self.hist2D[-1].Clone(f"Histogram_2D_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}_iteration_{i}")
                cumHist.SetTitle(f"Histogram 2D Ch {self.Channels[0]} vs ch {self.Channels[1]} iteration {i}")
                cumHist.Add(dayHist)
                self.hist2D.append(cumHist)
            
            self.hist2D[-1].SetStats(0)
            self.hist2D[-1].Draw('COLZ')
            self.Hist2DCanvas.Update()            
            

        saveFilePath = saveFilePath / 'rotated_Histograms' / f'Channels_{self.Channels[0]}_{self.Channels[1]}'
        saveFilePath.mkdir(parents=True,exist_ok=True)
        if scaleData:
            saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_scaled_cummulative_2D_Hist.png'
        else:
            saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_cummulative_2D_Hist.png'

        self.Hist2DCanvas.SaveAs(str(saveFileName))
        
        # self.Hist2DCanvasIndividual = ROOT.TCanvas(f"Individual_2D_Hist_Canvas_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}",f"Invidual 2D Hist Canvas Ch {self.Channels[0]} vs ch {self.Channels[1]}",2000,4000)
        
        # self.Hist2DCanvasIndividual.Divide(6,6)
        
        # for i in range(36):
        #     if i == 0:
        #         for x,y in zip(self.xInts[0:self.SeparationPoints[i]],self.yInts[0:self.SeparationPoints[i]]):
        #             self.hist2DIndividual[-1].Fill(x,y)
        #     else:
        #         for x,y in zip(self.xInts[self.SeparationPoints[i-1]:self.SeparationPoints[i]],self.yInts[self.SeparationPoints[i-1]:self.SeparationPoints[i]]):
        #             self.hist2DIndividual[-1].Fill(x,y)
        
        # self.hist2DIndividual[-1].SetStats(0)
        # self.hist2DIndividual[-1].Draw('COLZ')
        # self.Hist2DCanvasIndividual.Update()
        

        # if scaleData:
        #     saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_scaled_non_cummulative_2D_Hist.png'
        # else:
        #     saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_noncummulative_2D_Hist.png'

        # self.Hist2DCanvasIndividual.SaveAs(str(saveFileName))
    # self.sumChannels = np.concatenate((self.sumChannels,(np.sin(np.arctan(slope))*xData + np.cos(np.arctan(slope))*yData)))
    
    def plotRotatedHistograms(self,saveFilePath,scaleData,slopeDict,threshold,bins):
        # bins = 100
        
        if scaleData:
            # threshold = 550
            binRangeX = [-800,800]
            binRangeY = [0,1000]
        else:
            # threshold = 900
            binRangeX = [-1500,1500]
            binRangeY = [0,1750]
            
        self.Hist2DRotCanvas = ROOT.TCanvas(f"2D_Hist_Canvas_Rot_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}",f"2D Hist Canvas Rotated Ch {self.Channels[0]} vs ch {self.Channels[1]}",2000,4000)
        self.Hist2DRotCanvas.Divide(6,6)
        
        self.hist2DRot = []
        
        ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)
        
        self.rotatedHistData(scaleData,slopeDict)
        
        for i in range(36):
        
            self.Hist2DRotCanvas.cd(i+1)
        
            dayHist = TH2D(f"Histogram_2D_Rot_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}_iteration_{i}",f"Histogram 2D Rotated Ch {self.Channels[0]} vs ch {self.Channels[1]} iteration {i}",bins,binRangeX[0],binRangeX[1],bins,binRangeY[0],binRangeY[1])
            
            dayHist.Sumw2()
            
            if i ==0:
                xSlice = self.rotatedXData[0:self.SeparationPoints[i]]
                ySlice = self.rotatedYData[0:self.SeparationPoints[i]]
            else:
                xSlice = self.rotatedXData[self.SeparationPoints[i-1]:self.SeparationPoints[i]]
                ySlice = self.rotatedYData[self.SeparationPoints[i-1]:self.SeparationPoints[i]]
                
            for x,y in zip(xSlice,ySlice):
                if y > threshold:
                    dayHist.Fill(x,y)
                
            dayHist.Scale(1.0/self.time[i])
            
            
            if i == 0:
                self.hist2DRot.append(dayHist)
            else:
                cumHist = self.hist2DRot[-1].Clone(f"Histogram_2D_Rot_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}_iteration_{i}")
                cumHist.SetTitle(f"Histogram 2D Rotated Ch {self.Channels[0]} vs ch {self.Channels[1]} iteration {i}")
                cumHist.Add(dayHist)
                self.hist2DRot.append(cumHist)
                
            self.hist2DRot[-1].SetStats(0)
            self.hist2DRot[-1].Draw('COLZ')
            self.Hist2DRotCanvas.Update()
            
        # for i in range(36):
        #     self.Hist2DRotCanvas.cd(i+1)
            
        #     self.hist2DRot.append(TH2D(f"Histogram_2D_Rot_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}_iteration_{i}",f"Histogram 2D Rotated Ch {self.Channels[0]} vs ch {self.Channels[1]} iteration {i}",bins,binRangeX[0],binRangeX[1],bins,binRangeY[0],binRangeY[1]))
                        
        #     for x,y in zip(self.rotatedXData[0:self.SeparationPoints[i]],self.rotatedYData[0:self.SeparationPoints[i]]):
        #         if y > threshold:
        #             self.hist2DRot[-1].Fill(x,y)
                
        #     # Hist2D.GetXaxis().SetTitle(f"")    
        #     self.hist2DRot[-1].SetStats(0)
        #     self.hist2DRot[-1].Draw('COLZ')
        #     self.Hist2DRotCanvas.Update()
            
        saveFilePath = saveFilePath / 'rotated_Histograms' / f'Channels_{self.Channels[0]}_{self.Channels[1]}'
        
        saveFilePath.mkdir(parents=True,exist_ok=True)
        if scaleData:
            saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_scaled_cummulative_2D_Hist_Rotated.png'
        else:
            saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_cummulative_2D_Hist_Rotated.png'
        self.Hist2DRotCanvas.SaveAs(str(saveFileName))
        
        
        self.Hist2DCanvasIndividual = ROOT.TCanvas(f"Individual_2D_Hist_Canvas_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}",f"Invidual 2D Hist Canvas Ch {self.Channels[0]} vs ch {self.Channels[1]}",2000,4000)
            
        self.Hist2DCanvasIndividual.Divide(6,6)
        self.hist2DIndividual = []
        
        for i in range(36):
            self.Hist2DCanvasIndividual.cd(i+1)
            self.hist2DIndividual.append(TH2D(f"Individual_Histogram_2D_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}_iteration_{i}",f"Individual_Histogram 2D Ch {self.Channels[0]} vs ch {self.Channels[1]} iteration {i}",bins,binRangeX[0],binRangeX[1],bins,binRangeY[0],binRangeY[1]))
            
            self.hist2DIndividual[-1].Sumw2()
            if i == 0:
                for x,y in zip(self.rotatedXData[0:self.SeparationPoints[i]],self.rotatedYData[0:self.SeparationPoints[i]]):
                    if y > threshold:
                        self.hist2DIndividual[-1].Fill(x,y)
            else:
                for x,y in zip(self.rotatedXData[self.SeparationPoints[i-1]:self.SeparationPoints[i]],self.rotatedYData[self.SeparationPoints[i-1]:self.SeparationPoints[i]]):
                    if y > threshold:
                        self.hist2DIndividual[-1].Fill(x,y)
            self.hist2DIndividual[-1].Scale(1.0/self.time[i])
            self.hist2DIndividual[-1].SetStats(0)
            self.hist2DIndividual[-1].Draw('COLZ')
            self.Hist2DCanvasIndividual.Update()
        

        if scaleData:
            saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_scaled_non_cummulative_2D_Hist.png'
        else:
            saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_noncummulative_2D_Hist.png'

        self.Hist2DCanvasIndividual.SaveAs(str(saveFileName))
        
    def PlotSummedIntegral(self,saveFilePath,scaleData,threshold,bins):
        if self.rotatedXData is None or self.rotatedYData is None:
            self.rotatedHistData(scaleData)
        
        if scaleData:
            # threshold = 550
            binRangeY = [0,1000]
        else:
            # threshold = 900
            binRangeY = [0,1750]
            
        # bins = 100
        
        fitRange = [threshold,binRangeY[1]]
            
        self.SummedIntHistCanvas = ROOT.TCanvas(f"Summed_int_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}",f"Summed Int Ch {self.Channels[0]} vs ch {self.Channels[1]}",4000,4000)
        self.SummedIntHistCanvas.Divide(6,6)
        
        self.summedIntHist = []
        self.vLine = []
        self.stdLLine = []
        self.stdHLine = []
        self.mean = []
        self.stdev = []
        self.fitMean,self.fitMeanErr = [],[]
        self.fitSigma, self.fitSigmaErr = [],[]
        self.fits = []
        self.legends = []
         
                
        for i in range(36):
            self.SummedIntHistCanvas.cd(i+1)
            
            dayHist = ROOT.TH1D(f"Summed_Int_hist_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}_iteration_{i}",f"Summed Int Hist Ch {self.Channels[0]} vs ch {self.Channels[1]} iteration {i}",bins,binRangeY[0],binRangeY[1])
            dayHist.Sumw2()
            
            if i == 0:
                ySlice = self.rotatedYData[0:self.SeparationPoints[i]]
            else:
                ySlice = self.rotatedYData[self.SeparationPoints[i-1]:self.SeparationPoints[i]]
            
            for y in ySlice:
                if y > threshold:
                    dayHist.Fill(y)
            
            dayHist.Scale(1.0 / self.time[i])
                
            if i == 0:
                self.summedIntHist.append(dayHist)
            else:
                cumHist = self.summedIntHist[-1].Clone(f"Summed_Int_hist_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}_iteration_{i}")
                cumHist.SetTitle(f"Summed Int Hist Ch {self.Channels[0]} vs ch {self.Channels[1]} iteration {i}")
                cumHist.Add(dayHist)
                self.summedIntHist.append(cumHist)
            
            # self.summedIntHist.append(ROOT.TH1D(f"Summed_Int_hist_ch_{self.Channels[0]}_vs_ch_{self.Channels[1]}_iteration_{i}",f"Summed Int Hist Ch {self.Channels[0]} vs ch {self.Channels[1]} iteration {i}",bins,binRangeY[0],binRangeY[1]))
                        
            # for y in self.rotatedYData[0:self.SeparationPoints[i]]:
            #     if y > threshold:
            #         self.summedIntHist[-1].Fill(y)
                
            # # Hist2D.GetXaxis().SetTitle(f"")    
            self.summedIntHist[-1].SetStats(0)
            self.summedIntHist[-1].Draw('HIST')
            
            yMax = self.summedIntHist[-1].GetMaximum() * 1.05
            
            
            ################## --- Truncated Gaussian fit --- ##################
            truncFormula = (
                f"[0]*TMath::Gaus(x,[1],[2],1)/"
                f"(ROOT::Math::normal_cdf({fitRange[1]},[2],[1]) - ROOT::Math::normal_cdf({fitRange[0]},[2],[1]))"
            )
            func = ROOT.TF1(f"truncGaus_{i}", truncFormula, fitRange[0], fitRange[1])
            
            # Initial parameter guesses
            ampGuess = self.summedIntHist[-1].Integral("width")
            func.SetParameters(ampGuess, self.summedIntHist[-1].GetMean(), self.summedIntHist[-1].GetStdDev())
            func.SetParNames("Amplitude", "Mean", "Sigma")
            ################## --- Truncated Gaussian fit --- ##################
            
            ################## --- Truncated Gaussian + Erf fit --- ##################
            # truncFormula = (
            #     f"[0]*TMath::Gaus(x,[1],[2],1)/"
            #     f"(ROOT::Math::normal_cdf({fitRange[1]},[2],[1]) - ROOT::Math::normal_cdf({fitRange[0]},[2],[1]))"
            #     f" + [3]*0.5*TMath::Erfc((x-[1])/([2]*sqrt(2)))"
            # )
            
            # ampGuess = self.summedIntHist[-1].Integral("width")
            # stepGuess = self.summedIntHist[-1].GetMean() - self.summedIntHist[-1].GetStdDev()
            # func = ROOT.TF1(f"truncGausErf_{i}", truncFormula, fitRange[0], fitRange[1])
            # func.SetParameters(ampGuess, self.summedIntHist[-1].GetMean(), self.summedIntHist[-1].GetStdDev(), stepGuess)
            # func.SetParNames("Amplitude", "Mean", "Sigma", "StepHeight")
            ################## --- Truncated Gaussian + Erf fit --- ##################
            
            
            fitResult = self.summedIntHist[-1].Fit(func, "RSQ")  # R=restrict to range, S=return result, Q=quiet
            
            self.fits.append(func)
            func.SetLineColor(ROOT.kGreen+2)
            func.Draw("SAME")
            
            self.fitMean.append(func.GetParameter(1))
            self.fitMeanErr.append(func.GetParError(1))
            self.fitSigma.append(func.GetParameter(2))
            self.fitSigmaErr.append(func.GetParError(2))    
            
            chi2     = func.GetChisquare()
            ndof     = func.GetNDF()
            p_value  = Stats.chi2.sf(chi2, ndof)  
            
            
            legend = ROOT.TLegend(0.17,0.65,0.60,0.88)
            legend.SetHeader("Fit Results", "C")
            legend.AddEntry(self.fits[-1],f"Mean: = {round(self.fitMean[-1],2)} #pm {round(self.fitMeanErr[-1],2)}", "l")
            legend.AddEntry((0), f"Sigma = {round(self.fitSigma[-1],)} #pm {round(self.fitSigmaErr[-1],2)}") 
            legend.AddEntry((0), f"#chi^{{2}} / ndof = {chi2:.2f} / {ndof} = {round(chi2/ndof,2)}", "")
            legend.AddEntry((0), f"p-value = {p_value:.4f}", "")
            
            print(f"Fit Status: {fitResult.Status()} \t Chi2: {chi2:.2f} \t ndof: {ndof} \t chi2/ndof: {chi2/ndof}")
            
            legend.SetBorderSize(1)
            legend.SetFillColorAlpha(0, 0)  # semi-transparent background
            self.legends.append(legend)               # store to prevent garbage collection
            self.legends[-1].Draw()
            # --- Truncated Gaussian fit ---
            # truncFormula = (
            #     f"[0]*TMath::Gaus(x,[1],[2],1)/"
            #     f"(ROOT::Math::normal_cdf({fitRange[1]},[2],[1]) - ROOT::Math::normal_cdf({fitRange[0]},[2],[1]))"
            # )
            # func = ROOT.TF1(f"truncGaus_{i}", truncFormula, fitRange[0], fitRange[1])
            
            # # Initial parameter guesses
            # ampGuess = self.summedIntHist[-1].Integral("width")
            # func.SetParameters(ampGuess, self.summedIntHist[-1].GetMean(), self.summedIntHist[-1].GetStdDev())
            # func.SetParNames("Amplitude", "Mean", "Sigma")
            
            # fitResult = self.summedIntHist[-1].Fit(func, "RSQ")  # R=restrict to range, S=return result, Q=quiet
            
            # self.fits.append(func)
            # func.SetLineColor(ROOT.kGreen+2)
            func.Draw("SAME")
            
            # self.fitMean.append(func.GetParameter(1))
            # self.fitMeanErr.append(func.GetParError(1))
            # self.fitSigma.append(func.GetParameter(2))
            # self.fitSigmaErr.append(func.GetParError(2))
            # print(f"Fit Mean: {self.fitMean}")
            # print(f"Fit STD: {self.fitSigma}")
            
            self.mean.append(self.summedIntHist[-1].GetMean())
            self.stdev.append(self.summedIntHist[-1].GetStdDev())
            
            # print(f"ROOT Mean: {self.mean}")
            
            self.vLine.append(ROOT.TLine(self.fitMean[-1],0,self.fitMean[-1],yMax))
            self.vLine[-1].SetLineColor(ROOT.kBlack)
            self.vLine[-1].SetLineWidth(2)
            self.vLine[-1].SetLineStyle(2)
            self.vLine[-1].Draw("SAME")
            
            self.stdLLine.append(ROOT.TLine(self.fitMean[-1]-self.fitSigma[-1],0,self.fitMean[-1]-self.fitSigma[-1],yMax))
            self.stdLLine[-1].SetLineColor(ROOT.kRed)
            self.stdLLine[-1].SetLineWidth(2)
            self.stdLLine[-1].SetLineStyle(2)
            self.stdLLine[-1].Draw("SAME")
            
            
            self.stdHLine.append(ROOT.TLine(self.fitMean[-1]+self.fitSigma[-1],0,self.fitMean[-1]+self.fitSigma[-1],yMax))
            self.stdHLine[-1].SetLineColor(ROOT.kRed)
            self.stdHLine[-1].SetLineWidth(2)
            self.stdHLine[-1].SetLineStyle(2)
            self.stdHLine[-1].Draw("SAME")
            
            
            self.SummedIntHistCanvas.Update()
            
        saveFilePath = saveFilePath / 'rotated_Histograms' / f'Channels_{self.Channels[0]}_{self.Channels[1]}'
        
        saveFilePath.mkdir(parents=True,exist_ok=True)
        if scaleData:
            saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_scaled_summed_int_dist.pdf'
        else:
            saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_summed_int_dist.pdf'
        self.SummedIntHistCanvas.SaveAs(str(saveFileName))
        
        if scaleData:
            saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_scaled_summed_int_dist.png'
        else:
            saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_summed_int_dist.png'
        self.SummedIntHistCanvas.SaveAs(str(saveFileName))
        
        
        for i,mean in enumerate(self.fitMean):
             print(f"Day Number: {i+1} \t Fit Mean: {mean} \t Fit Sigma: {self.fitSigma[i]}, \t ROOT mean: {self.mean[i]}")
        
        self.stdCanvas = TCanvas(f"Std_ch_{self.Channels[0]}_ch_{self.Channels[1]}_canvas", f"Std ch {self.Channels[0]} ch {self.Channels[1]} canvas",2400,600)
        # self.stdCanvas.SetLeftMargin(0.15)
        self.stdCanvas.Divide(3,1)
        
        
        stdDistIndex = np.array([i for i in range(len(self.summedIntHist))], dtype='float64')
        # stdDist = np.array([hist.GetStdDev() for hist in self.summedIntHist], dtype='float64')
        # stdDistErr = np.array([hist.GetStdDevError() for hist in self.summedIntHist], dtype='float64')
        stdDistIndexErr = np.zeros(len(self.summedIntHist), dtype='float64')
            
        stdDist = np.array([std for std in self.fitSigma], dtype='float64')    
        stdDistErr = np.array([stdErr for stdErr in self.fitSigmaErr], dtype='float64')
        
        meanDist = np.array([m for m in self.fitMean], dtype='float64')

        meanDistErr = np.array([mErr for mErr in self.fitMeanErr], dtype='float64')
        
        resolution = np.array([s/m for s,m in zip(stdDist,meanDist)], dtype='float64')
        resolutionErr = np.array([np.sqrt((s/m)**2 * ((sErr/s)**2 + (mErr/m)**2)) for s,sErr,m,mErr in zip(stdDist,stdDistErr,meanDist,meanDistErr)], dtype='float64')
        
        
        #######################################
        # Plot the standard deviation
        #######################################
        pad = self.stdCanvas.cd(1)
        self.stdDist = ROOT.TGraphErrors(len(self.summedIntHist),stdDistIndex,stdDist,stdDistIndexErr,stdDistErr)
        pad.SetLeftMargin(0.25)
        self.stdDist.SetMarkerStyle(20)
        self.stdDist.SetMarkerSize(1.0)
        self.stdDist.SetMarkerColor(ROOT.kBlue)
        self.stdDist.GetXaxis().SetTitle("Day Number")
        if scaleData:
            self.stdDist.SetTitle("Std of Histograms scaled to keV")
            saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_std_dist_scaled.png'
            self.stdDist.GetYaxis().SetTitle("Standard Deviation (keV)")
        else:
            self.stdDist.SetTitle("Std of Histograms in ADC")
            saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_std_dist.png'
            self.stdDist.GetYaxis().SetTitle("Standard Deviation (ADC)")
            
        self.stdDist.Draw("AP")
        self.stdCanvas.Update()
        

        
        #######################################
        # Plot the mean
        #######################################
        pad = self.stdCanvas.cd(2)
        self.meanDist = ROOT.TGraphErrors(len(self.summedIntHist),stdDistIndex,meanDist,stdDistIndexErr,meanDistErr)
        pad.SetLeftMargin(0.25)
        self.meanDist.SetMarkerStyle(20)
        self.meanDist.SetMarkerSize(1.0)
        self.meanDist.SetMarkerColor(ROOT.kRed)
        self.meanDist.GetXaxis().SetTitle("Day Number")
        if scaleData:
            self.meanDist.SetTitle("Mean of Histograms scaled to keV")
            # saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_std_dist_scaled.png'
            self.meanDist.GetYaxis().SetTitle("Mean (keV)")
        else:
            self.meanDist.SetTitle("Mean of Histograms in ADC")
            # saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_std_dist.png'
            self.meanDist.GetYaxis().SetTitle("Mean (ADC)")
        
        self.meanDist.Draw("AP")
        self.stdCanvas.Update()
        
        #######################################
        # Plot the resolution
        #######################################
        pad = self.stdCanvas.cd(3)
        self.resDist = ROOT.TGraphErrors(len(self.summedIntHist),stdDistIndex,resolution,stdDistIndexErr,resolutionErr)
        pad.SetLeftMargin(0.25)
        self.resDist.SetMarkerStyle(20)
        self.resDist.SetMarkerSize(1.0)
        self.resDist.SetMarkerColor(ROOT.kRed)
        self.resDist.GetXaxis().SetTitle("Day Number")

        self.resDist.SetTitle("Resolution of Histograms")
            # saveFileName = saveFilePath / f'Channels_{self.Channels[0]}_{self.Channels[1]}_std_dist_scaled.png'
        self.resDist.GetYaxis().SetTitle("Resolution")
        
        self.resDist.Draw("AP")
        self.stdCanvas.Update()
        
        self.stdCanvas.SaveAs(str(saveFileName))

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
                    t.append(float(data[6*j+2])/1e12) #Grab the time data for scaling purposes. Convert to s from ps.
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
                    cData = coincData(ch,E,t)
                else:
                    cData.AddEvent(ch,E,t)
                    
    return cData


def ReadInSlopes(file):
    slopeDict = {}
    
    with open(file) as f:
        next(f) #skip line 1
        
        for line in f:
            data = line.split('\n')[0].split(',')
            
            if int(data[1]) == 4 or int(data[1]) == 5:
                date = data[0].split(' 00:00:00')[0].replace("-","_")
                channels = [int(data[1]), int(data[11])]
                slope = float(data[9])
                slopeErr = float(data[10])
                
                print(f'Reading in slope data for:\n Slope Date {date} \tchannels {channels}')
                
                chstring = f"{channels[0]},{channels[1]}"
                if chstring in slopeDict:
                    slopeDict[chstring][date] = [1/slope,slopeErr/slope**2]
                else:
                    slopeDict[chstring] = {date: [1/slope,slopeErr/slope**2]}
                
                # for i,hist in enumerate(histData):
                #     if hist.Channels == channels:
                #         hist.slope.append(1/slope)
                #         hist.slopeErr.append(slopeErr/slope**2)
                #         print(f"Histogram Date: {hist.dates[i]}\t Channels: {hist.Channels}\t slope: {slope}")
                
            else: 
                pass
                
    return slopeDict     


def readInScaleFactors(filepath, channels):
    
    scale = []
    for ch in channels:
        scale.append(ScaleFactors(ch))
    with open(filepath) as f:
        next(f) #Skips the header
        
        for line in f:
            data = line.split(',')
            for s in scale:
                if s.ch == int(data[1]):
                    s.AddData(data[0].split(' 00:00:00')[0], float(data[5])) #index 3 is the scale factor for the X profile fits. Index 5 is the scale factor for the Y profile fits. 

    return scale

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




rootfilePath = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/') #filepath to the coinc sorted directory. 

#A list of data sets that are excluded from the data. This can be due to bad data or incorrect settings.
excludedDataFiles = [Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/2026_04_28/2026_04_28_Daily_LSC_calibration_Cs137_coinc'), #Data set was taken before the settings were finalized
                     ]

averageScaleFactorFP = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/Results/Annulus_stability_data_average.txt')

slopeFilePath = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/Results/Annulus_stability_data.txt')

savefilepath = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/stability_figures/')

LSCChannels = [4,5]
NaIChannels = [8,10,12,14] #Protects against the possibility of having a werid channel coincidence layout with incorrect channel numbers. 

# LSCChannels = [0,1]
# NaIChannels = [2,3,4,5]

scaleData = False
fitData = False
if scaleData:
    scaleFac = readInScaleFactors(averageScaleFactorFP, channels = LSCChannels + NaIChannels)

pattern = re.compile(rf'_coinc_{LSCChannels[0]}_{LSCChannels[1]}_({NaIChannels[0]}|{NaIChannels[1]}|{NaIChannels[2]}|{NaIChannels[3]}).txt$')

coincFiles = sorted([
    f for f in rootfilePath.glob(f'**/*_coinc_{LSCChannels[0]}_{LSCChannels[1]}_*.txt')
    if pattern.search(f.name)
])
histogramData = []
for i in LSCChannels:
    for j in NaIChannels:
        histogramData.append(Summed2DHist([i,j]))
        print(histogramData[-1].Channels)


readInStartTime = time.time()
print('Reading in histogram data:')

for i,filePath in enumerate(coincFiles):
    if filePath.parent.parent.parent in excludedDataFiles:
        print(f'Skipped file: {filePath}')
        pass
    else:
        date = filePath.parent.parent.parent.parent.stem
        settingsFilePath = filePath.parent.parent.parent / 'settings.xml'
        # saveFilePath = filePath.parent.parent / 'Daily_calibration_fits' / 'figures'
        # outputFileName = filePath.stem
        # outputFile = saveFilePath / f"{outputFileName}_Profile_hist_fit_results.txt"
        
        # if outputFile.exists():
        #     outputFile.unlink()

        # saveFilePath.mkdir(parents=True, exist_ok=True)

        if i == 0:
            Detectors = ReadInChannelNames(settingsFilePath)
        
        if scaleData:
            cData = readInData(filePath,scaleFac)
        else:
            cData = readInData(filePath)
            
        chPairs = [[cData.chList[0],cData.chList[2]], [cData.chList[1],cData.chList[2]]] #Makes a 2D array with the two channel pairs for the summed histogram data. 
        print(cData.chList)
        # print(len(cData.chData['4'].E))


        for j,hist in enumerate(histogramData):
            if hist.Channels == chPairs[0]:
                hist.addHistData(cData.chData[str(chPairs[0][0])].E,cData.chData[str(chPairs[0][1])].E,cData.chData[str(chPairs[0][0])].t,date)
            elif hist.Channels == chPairs[1]:
                hist.addHistData(cData.chData[str(chPairs[1][0])].E,cData.chData[str(chPairs[1][1])].E,cData.chData[str(chPairs[0][0])].t,date)



totalReadTime = time.time()
print(f'Total time to read in data: \t {totalReadTime - readInStartTime} s')

print("reading in Slope information")
slopeDict = ReadInSlopes(slopeFilePath)
slopeReadInTime = time.time()
print(f"Time to read in slope data: {slopeReadInTime - totalReadTime} s")

nbins = 100
if scaleData:
    threshold = 550
else:
    threshold = 900

for hist in histogramData:
    print(f'Plotting 2D Histogram from channels: {hist.Channels}')
    print(hist.dates)
    hist.plot2DHist(savefilepath,scaleData,nbins)
    hist.plotRotatedHistograms(savefilepath,scaleData,slopeDict,threshold,nbins)
    hist.PlotSummedIntegral(savefilepath,scaleData,threshold,nbins)
    break