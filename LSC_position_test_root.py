import numpy as np
import ROOT
from pathlib import Path
import xml.etree.ElementTree as ET
from itertools import combinations
import math


#########################################################################################
# Global Settings
#########################################################################################
# ROOT.gROOT.SetBatch(True)
# ROOT.EnableImplicitMT() #Enables implicit multithreading to speed up the code. 

ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)

ROOT.gStyle.SetLabelSize(0.05, "XYZ")   # axis tick-label size (numbers along the axis)
ROOT.gStyle.SetTitleSize(0.05, "XYZ")   # axis title size (e.g. "Energy (ADC)")
ROOT.gStyle.SetTitleFontSize(0.06)      # the histogram's own title (top of plot)
ROOT.gStyle.SetTextSize(0.05)           # general text (legends, TLatex, etc., unless overridden)
ROOT.gStyle.SetLegendTextSize(0.1)   

colourCycle = [ROOT.kBlue,ROOT.kOrange+2, ROOT.kGreen+2]
fitColourCycle = [ROOT.kRed,ROOT.kGreen+2]
#########################################################################################

class positionData:
    def __init__(self,file,channels,ncols,fitWindow,saveFP):
        self.file = ROOT.TFile.Open(str(file))
        # self.file.ls() #Prints the contents of the file 
        
        settingsFilePath = file.parent.parent.parent / 'settings.xml'
        
        self.saveFP = saveFP / "Position_Fits"
        self.saveFP.mkdir(parents = True, exist_ok = True)
        
        
        try:
            self.trialNumber = int(file.stem.split("_")[10]) -1 
        except:
            self.trialNumber = 0
            

        
        self.chLabels = ReadInChannelNames(settingsFilePath)
        
        self.filename = file.stem
        self.positionSTR = self.filename.split("_")[6]
        
        self.channels = channels
        
        if 'R' in self.positionSTR:
            self.position = int(self.positionSTR[:-1])
        elif 'L' in self.positionSTR:
            self.position = -1 * int(self.positionSTR[:-1])
        else:
            self.position = int(self.positionSTR)

        if "Left" in self.chLabels[f'Channel 0'][0]:
            self.Lch = 0
            self.Rch = 1
        else:
            self.Lch = 1
            self.Rch = 0
   
        self.NaIch = 2
        
        
        self.fitWindow = fitWindow
        
        ################################################################
        # Define the fits 
        ################################################################
        
        fitNames = ['Ngauss', 'mu', 'sigma','Nexp', 'tau', 'erfStep', ' erfMu', 'erfSigma']
        
        tag = f"{self.position}_{self.trialNumber}"
        
        #Left PMT fits
        #No erf
        self.fitFuncL = ROOT.TF1(
            f"fit_L_{tag}", 'gaus(0) + expo(3)',
            fitWindow[f'{self.position}'][0][0],
            fitWindow[f'{self.position}'][0][1])
        
        
        
        self.fitFuncL.SetLineColor(ROOT.kRed)
        self.fitFuncL.SetLineWidth(3)
        
        #erf
        self.fitFuncLerf = ROOT.TF1(
            f"fit_L_erf_{tag}", 'gaus(0) + expo(3) + [5]*TMath::Erfc((x-[6])/[7])',
            fitWindow[f'{self.position}'][0][0],
            fitWindow[f'{self.position}'][0][1])
        
        self.fitFuncLerf.SetLineColor(ROOT.kGreen+2)
        self.fitFuncLerf.SetLineWidth(3)
        
        #Right PMT fits
        #No erf
        self.fitFuncR = ROOT.TF1(
            f"fit_R_{tag}", 'gaus(0) + expo(3)',
            fitWindow[f'{self.position}'][1][0],
            fitWindow[f'{self.position}'][1][1])
        
        self.fitFuncR.SetLineColor(ROOT.kRed)
        self.fitFuncR.SetLineWidth(3)
        
        #erf
        self.fitFuncRerf = ROOT.TF1(
            f"fit_R_erf_{tag}", 'gaus(0) + expo(3) + [5]*TMath::Erfc((x-[6])/[7])',
            fitWindow[f'{self.position}'][1][0],
            fitWindow[f'{self.position}'][1][1])
        
        self.fitFuncRerf.SetLineColor(ROOT.kGreen+2)
        self.fitFuncRerf.SetLineWidth(3)
        
        
        for i,n in enumerate(fitNames):
            self.fitFuncL.SetParName(i,n)
            self.fitFuncLerf.SetParName(i,n)
            self.fitFuncR.SetParName(i,n)
            self.fitFuncRerf.SetParName(i,n)
        
        
        self.LpadInd = self.position + ncols//2 +1 #Sets the position on the pad for the left PMT
        self.RpadInd = self.position + ncols //2 + ncols +1 #Sets the position on the pad for the right PMT
        
        self.cols = []
        for ch in channels:
            ChCol = [f"Timestamp_{ch}",f"Board_{ch}_{self.trialNumber}",f"Energy_{ch}",f"EnergyShort_{ch}",f"Flags_{ch}"]
            self.cols += ChCol
            
        self.cPairsStr = '_'.join([str(i) for i in channels])
        self.tree = self.file.Get(f"coinc_ch_{self.cPairsStr}")
        
        self.df = ROOT.RDataFrame(self.tree)
        
        
    def calcLogRatio(self):
        logRatioCols = []
        
        # dfPos = self.df.Filter(f"Energy_{self.Lch} > 0 && Energy_{self.Rch} > 0", "positive energies")
        
        colName = f"LogRatio"
        self.df = self.df.Define(colName, f"std::log(static_cast<double>(Energy_{self.Rch} + 0.00000000001) / (Energy_{self.Lch}+ 0.00000000001))")
        
    def plotLogRatioHist(self):
        
        self.calcLogRatio()
        

        dfClean = self.df.Filter("std::isfinite(LogRatio)", "finite log ratio")
        dfClean.Report().Print()
        # dfClean = dfClean.Filter("LogRatio >= -0.05 && LogRatio <= 0.05")
        # self.logRatioHist = self.df.Histo1D(
        #     (f"Log_Ratio_Pos_{self.positionSTR}_iteration_{self.trialNumber}",
        #      f"Position {self.positionSTR} log Ratio;ln#left(#frac{{S_{{R}}}}{{S_{{L}}}}#right);Counts",
        #      100,-2,2),
        #     "LogRatio"
        # )
        
        self.logRatioHist = dfClean.Histo1D(
            (f"Log_Ratio_Pos_{self.positionSTR}_iteration_{self.trialNumber}",
             f"Position {self.positionSTR} log Ratio;ln#left(#frac{{S_{{R}}}}{{S_{{L}}}}#right);Counts",
             200,-2,2),
            "LogRatio"
        )
        
        # self.LogRatioEnergyL = self.df.Histo2D(
        #     (f"Pos_{self.positionSTR}_iter_{self.trialNumber}_log_Ratio_Left",
        #      f"Pos {self.positionSTR} iter {self.trialNumber} log Ratio Left;Integral Left (ADC);log Ratio;ln#left(#frac{{S_{{R}}}}{{S_{{L}}}}#right)",
        #      100,0,4000,
        #      100,-2,2),
        #     f"Energy_{self.Lch}",
        #     "LogRatio")
        
        # self.LogRatioEnergyR = self.df.Histo2D(
        #     (f"Pos_{self.positionSTR}_iter_{self.trialNumber}_log_Ratio_Right",
        #      f"Pos {self.positionSTR} iter {self.trialNumber} log Ratio Right;Integral Left (ADC);log Ratio;ln#left(#frac{{S_{{R}}}}{{S_{{L}}}}#right)",
        #      100,0,4000,
        #      100,-2,2),
        #     f"Energy_{self.Rch}",
        #     "LogRatio")
        
        self.LogRatioEnergyL = dfClean.Histo2D(
            (f"Pos_{self.positionSTR}_iter_{self.trialNumber}_log_Ratio_Left",
             f"Pos {self.positionSTR} iter {self.trialNumber} log Ratio Left;Integral Left (ADC);ln#left(#frac{{S_{{R}}}}{{S_{{L}}}}#right);Counts",
             100,0,4000,
             200,-2,2),
            f"Energy_{self.Lch}",
            "LogRatio")
        self.LogRatioEnergyL.SetStats(0)
        self.LogRatioEnergyR = dfClean.Histo2D(
            (f"Pos_{self.positionSTR}_iter_{self.trialNumber}_log_Ratio_Right",
             f"Pos {self.positionSTR} iter {self.trialNumber} log Ratio Right;Integral Left (ADC);ln#left(#frac{{S_{{R}}}}{{S_{{L}}}}#right);Counts",
             100,0,4000,
             200,-2,2),
            f"Energy_{self.Rch}",
            "LogRatio")
        self.LogRatioEnergyR.SetStats(0)
        
        
        self.EnergyLvR = dfClean.Histo2D(
            (f"Pos_{self.positionSTR}_iter_{self.trialNumber}_Left_vs_Right",
             f"Pos {self.positionSTR} iter {self.trialNumber} Left vs Right;Integral Left (ADC); Integral Right (ADC); Counts",
             100,0,4000,
             100,0,4000),
             f"Energy_{self.Lch}",
             f"Energy_{self.Rch}")

        
        self.EnergyLvR.SetStats(0)
        
        self.LogRatioEnergyL.GetXaxis().SetNdivisions(505) 
        self.LogRatioEnergyR.GetXaxis().SetNdivisions(505) 
        
        self.EnergyLvR.GetXaxis().SetNdivisions(505)
        self.EnergyLvR.GetYaxis().SetNdivisions(505)
        
        #####################################################################
        #Make the 1D Histograms
        self.LogRatioCanvas = ROOT.TCanvas(f"Log_Ratio_canvas_pos_{self.positionSTR}_iteration_{self.trialNumber}", f"Position {self.positionSTR} iteration {self.trialNumber} log ratio", 1400,800)
        
        pad = self.LogRatioCanvas.cd(1)
        canvasPadding(pad,0.15,0.1,0.15,0.2)
        
        self.logRatioHist.Draw()
        
        self.LogRatioCanvas.Update()
        
        self.LogRatioCanvas.SaveAs(str(self.saveFP / f"log_Ratio_{self.positionSTR}_{self.trialNumber}.pdf"))
        
        self.LogRatioCanvas.SaveAs(str(self.saveFP / f"log_Ratio_{self.positionSTR}_{self.trialNumber}.png"))
        
        ##################################################################
        # Make the 2D histograms
        
        
        self.LogRatioEnergyCanvas = ROOT.TCanvas(f"Pos_{self.positionSTR}_iter_{self.trialNumber}_log_vs_Energy", f"Pos {self.positionSTR} iter {self.trialNumber} log vs Energy", 1400,800)
        
        self.LogRatioEnergyCanvas.Divide(3,1)
        
        #Plot the Left PMT
        pad = self.LogRatioEnergyCanvas.cd(1)
        canvasPadding(pad,0.2,0.2,0.15,0.15)
        
        self.LogRatioEnergyL.Draw("COLZ")
        
        self.LogRatioEnergyCanvas.Update()
        
        
        #Plot the Right PMT
        pad = self.LogRatioEnergyCanvas.cd(2)
        canvasPadding(pad,0.2,0.2,0.15,0.15)
        
        self.LogRatioEnergyR.Draw("COLZ")
        
        self.LogRatioEnergyCanvas.Update()
        
        #Plot the Left Vs Right
        
        YeqX = ROOT.TH2F('Y_eq_x_line', "Y = X line", 100,0,4000,100,0,4000)
        self.YeqXLine = ROOT.TF1("YeqX", "x", 0,4000)
        self.YeqXLine.SetLineColor(ROOT.kRed)
        self.YeqXLine.SetLineWidth(2)
        self.YeqXLine.SetLineStyle(2)
        
        pad = self.LogRatioEnergyCanvas.cd(3)
        canvasPadding(pad,0.2,0.2,0.15,0.15)
        
        self.EnergyLvR.Draw("COLZ")
        self.YeqXLine.Draw("SAME")
        
        self.LogRatioEnergyCanvas.Update()
        

        self.LogRatioEnergyCanvas.SaveAs(str(self.saveFP / f"log_Ratio_Energy_{self.positionSTR}_{self.trialNumber}.png"))
        
        self.LogRatioEnergyCanvas.SaveAs(str(self.saveFP / f"log_Ratio_Energy_{self.positionSTR}_{self.trialNumber}.pdf"))

        
    def plot1DHist(self,canvas):
        
        #####################################################################
        # Create the two histograms 
        #####################################################################
        
        print(f"\n Plotting Position {self.positionSTR} \n")
        
        
        #Left PMT Hist
        self.LSCLHist1D = self.df.Histo1D(
            (f"Energy_Hist_{self.Lch}_{self.trialNumber}", 
             f"{self.chLabels[f'Channel {self.Lch}'][0]} Position {self.positionSTR};{self.chLabels[f'Channel {self.Lch}'][0]} Integral (ADC);Counts", 100,0,4000),
            f"Energy_{self.Lch}")
        
        #Right PMT Hist
        self.LSCRHist1D = self.df.Histo1D((f"Energy_Hist_{self.Rch}_{self.trialNumber}", f"{self.chLabels[f'Channel {self.Rch}'][0]} Position {self.positionSTR};{self.chLabels[f'Channel {self.Rch}'][0]} Integral (ADC);Counts", 100,0,4000),f"Energy_{self.Rch}")
        
        #Turn off stats
        self.LSCLHist1D.SetStats(0)
        self.LSCRHist1D.SetStats(0)
        
        
        #Swap to the first pad to plot the Left PMT data
        pad = canvas.cd(self.LpadInd)
        
        canvasPadding(pad,0.2,0.1,0.15,0.1)
        # self.LSCLHist1D.Draw()
        opt = "" if pad.GetListOfPrimitives().GetSize() == 0 else "SAME"
        self.LSCLHist1D.SetLineColor(colourCycle[self.trialNumber % len(colourCycle)])
        self.LSCLHist1D.Draw(opt)
        
        #Set initial parameters for the fits. 
        ampGuess = 500
        minRangeL, maxRangeL = self.fitWindow[f'{self.position}'][0]
        minRangeR, maxRangeR = self.fitWindow[f'{self.position}'][1]
        sigGuessL = (maxRangeL-minRangeL) / 10
        muGuessL = (maxRangeL+minRangeL)/2
        
    
        
        sigGuessR = (maxRangeR-minRangeR) / 10
        muGuessR = (maxRangeR+minRangeR)/2
        
        expAmpGuess = 5
        tau = 1e-10
        stepHeight = 30
        self.fitFuncL.SetParameters(ampGuess,muGuessL,sigGuessL,expAmpGuess,tau)
        self.fitFuncR.SetParameters(ampGuess,muGuessR,sigGuessR,expAmpGuess,tau)
        
        # self.fitFuncL.SetParLimits(1,minRangeL,maxRangeL)
        # self.fitFuncR.SetParLimits(1,minRangeR,maxRangeR)
        for i,fitFunc in enumerate([self.fitFuncL,self.fitFuncR]):
            fitFunc.SetParLimits(1, self.fitWindow[f'{self.position}'][i][0],self.fitWindow[f'{self.position}'][i][1]) # mean stays in the window
            fitFunc.SetParLimits(2, 1, (self.fitWindow[f'{self.position}'][i][1]-self.fitWindow[f'{self.position}'][i][0])) # sigma positive
            # fitFunc.SetParLimits(4, 0, 10) # tau positive

        
        self.fitFuncLerf.SetParameters(ampGuess,muGuessL,sigGuessL,expAmpGuess,tau,stepHeight,muGuessL,sigGuessL)
        self.fitFuncRerf.SetParameters(ampGuess,muGuessR,sigGuessR,expAmpGuess,tau,stepHeight,muGuessR,sigGuessR)
        
        for i,fitFunc in enumerate([self.fitFuncLerf,self.fitFuncRerf]):
            fitFunc.SetParLimits(1, self.fitWindow[f'{self.position}'][i][0],self.fitWindow[f'{self.position}'][i][1]) # mean stays in the window
            fitFunc.SetParLimits(2, 1, (self.fitWindow[f'{self.position}'][i][1]-self.fitWindow[f'{self.position}'][i][0])) # sigma positive
            fitFunc.SetParLimits(5, 0, 1e4)          # step height positive
            fitFunc.SetParLimits(4, 0, 10) # tau positive
            fitFunc.SetParLimits(6, self.fitWindow[f'{self.position}'][i][0],self.fitWindow[f'{self.position}'][i][1])          # mean stays in the window
            fitFunc.SetParLimits(7, 1, (self.fitWindow[f'{self.position}'][i][1]-self.fitWindow[f'{self.position}'][i][0]))    # sigma positive
            

        
        #Fit the data no erf
        self.resL = self.LSCLHist1D.Fit(self.fitFuncL, "RS")
        self.fitFuncL.Draw("SAME")
        
        #Fit the data erf
        self.resLErf = self.LSCLHist1D.Fit(self.fitFuncLerf, "RS")
        self.fitFuncLerf.Draw("SAME")
        self.LSCLHist1D.GetXaxis().SetNdivisions(505) 
       
        
        canvas.Update()
        
        pad = canvas.cd(self.RpadInd)
        canvasPadding(pad,0.2,0.1,0.10,0.15)
        # self.LSCRHist1D.Draw()
        self.LSCRHist1D.SetLineColor(colourCycle[self.trialNumber % len(colourCycle)])
        self.LSCRHist1D.Draw(opt)
        
        #Fit the data no erf
        self.resR = self.LSCRHist1D.Fit(self.fitFuncR, "RS")
        self.fitFuncR.Draw("SAME")
        
        #Fit the data erf
        self.resRErf = self.LSCRHist1D.Fit(self.fitFuncRerf, "RS")
        self.fitFuncRerf.Draw("SAME")
        self.LSCRHist1D.GetXaxis().SetNdivisions(505) 
        
        canvas.Update()
        
        
        #Create the plot for the individual position fits
        self.FitCanvas = ROOT.TCanvas(f"Individual_Fit_canvas_{self.positionSTR}_{self.trialNumber}", f"Position {self.positionSTR} fit", 2800,1600)
        
        self.FitCanvas.Divide(1,2)
        
        pad = self.FitCanvas.cd(1)
        canvasPaddingLegend(pad)
        
        leftMin = self.fitFuncL.GetXmin()
        leftMax = self.fitFuncL.GetXmax()
        
        self.LSCLHist1D.Draw()
        self.fitFuncL.Draw("SAME")
        self.fitFuncLerf.Draw("SAME")
        
        ROOT.gPad.Update()
        Leftymin = ROOT.gPad.GetUymin()
        Leftymax = ROOT.gPad.GetUymax()
        
        self.LeftMinLine = ROOT.TLine(leftMin,Leftymin,leftMin,Leftymax)
        self.LeftMaxLine = ROOT.TLine(leftMax,Leftymin,leftMax,Leftymax)
        
        for line in [self.LeftMinLine, self.LeftMaxLine]:
            line.SetLineColor(ROOT.kBlack)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw("same")
    
        self.LmeanLine = ROOT.TLine(self.resL.Parameter(1),Leftymin,self.resL.Parameter(1),Leftymax)
        self.LErfmeanLine = ROOT.TLine(self.resLErf.Parameter(1),Leftymin,self.resLErf.Parameter(1),Leftymax)
        
        
        
        for line in [self.LmeanLine]:
            line.SetLineColor(ROOT.kRed)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            
            
        for line in [self.LErfmeanLine]:
            line.SetLineColor(ROOT.kGreen+2)
            line.SetLineStyle(2)
            line.SetLineWidth(2)

        self.LmeanLine.Draw("same")
        self.LErfmeanLine.Draw("same")
        
        self.legendL = self.MakeLegend(self.resL,self.resLErf, self.fitFuncL, self.fitFuncLerf)
        self.legendL.Draw("same")
        
        self.FitCanvas.Update()
        
        pad = self.FitCanvas.cd(2)
        canvasPaddingLegend(pad)
        
       
        
        self.LSCRHist1D.Draw()
        self.fitFuncR.Draw("SAME")
        self.fitFuncRerf.Draw("SAME")
        
        rightMin = self.fitFuncR.GetXmin()
        rightMax = self.fitFuncR.GetXmax()
        
        ROOT.gPad.Update()
        
        Rightymin = ROOT.gPad.GetUymin()
        Rightymax = ROOT.gPad.GetUymax()
        
        self.RightMinLine = ROOT.TLine(rightMin,Rightymin,rightMin,Rightymax)
        self.RightMaxLine = ROOT.TLine(rightMax,Rightymin,rightMax,Rightymax)
        
        self.RmeanLine = ROOT.TLine(self.resR.Parameter(1),Rightymin,self.resR.Parameter(1),Rightymax)
        self.RErfmeanLine = ROOT.TLine(self.resRErf.Parameter(1),Rightymin,self.resRErf.Parameter(1),Rightymax)
        
        for line in [self.RightMinLine, self.RightMaxLine]:
            line.SetLineColor(ROOT.kBlack)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw("same")
            
        for line in [self.RmeanLine]:
            line.SetLineColor(ROOT.kRed)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            
            
        for line in [self.RErfmeanLine]:
            line.SetLineColor(ROOT.kGreen+2)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
        
        self.RmeanLine.Draw("same")
        self.RErfmeanLine.Draw("same")
        self.legendR = self.MakeLegend(self.resR,self.resRErf, self.fitFuncR, self.fitFuncRerf)
        self.legendR.Draw("same")
        
        self.FitCanvas.Update()
        
        self.FitCanvas.SaveAs(str(self.saveFP / f"Individual_fits_{self.positionSTR}_{self.trialNumber}.pdf"))
        self.FitCanvas.SaveAs(str(self.saveFP / f"Individual_fits_{self.positionSTR}_{self.trialNumber}.png"))
        
    # def MakeLegend(self,fitRes,fitResErf, fitFunc, fitFuncErf):
    #     ok = fitRes.IsValid() and fitRes.Status() == 0
    #     okErf = fitResErf.IsValid() and fitResErf.Status() == 0
        
    #     mu, muErr = fitRes.Parameter(1), fitRes.ParError(1)
    #     muErf, muErrErf = fitResErf.Parameter(1), fitResErf.ParError(1)
            
            
    #     # Build legend — positioned to the left of the colour bar
    #         # (x1, y1, x2, y2) in NDC coordinates (0-1)
    #     legend = ROOT.TLegend(0.02, 0.05, 0.98, 0.95)
        
    #     legend.SetNColumns(2)
    #     legend.SetBorderSize(1)
    #     legend.SetFillColorAlpha(0, 0)
    #     legend.SetTextSize(0)   
        
    #     legend.SetHeader("Fit Results", "C")
    #     legend.AddEntry(fitFunc, f"#mu = {mu:.1f} #pm {muErr:.1f}", "l")
    #     legend.AddEntry(ROOT.nullptr, f"#chi^{{2}}/ndf = {fitRes.Chi2():.1f} / {fitRes.Ndf()}", "")
    #     legend.AddEntry(ROOT.nullptr, f"p-value = {fitRes.Prob():.3f}", "")

    #     status = legend.AddEntry(ROOT.nullptr, "Fit OK" if ok else f"Fit FAILED (status {fitRes.Status()})", "")
    #     status.SetTextColor(ROOT.kGreen + 2 if ok else ROOT.kRed)
        
    #     #########################################################################
        
    #     legend.AddEntry(fitFuncErf, f"#mu = {muErf:.1f} #pm {muErrErf:.1f}", "l")
    #     legend.AddEntry(ROOT.nullptr, f"#chi^{{2}}/ndf = {fitResErf.Chi2():.1f} / {fitResErf.Ndf()}", "")
    #     legend.AddEntry(ROOT.nullptr, f"p-value = {fitResErf.Prob():.3f}", "")

    #     status = legend.AddEntry(ROOT.nullptr, "Fit OK" if okErf else f"Fit FAILED (status {fitResErf.Status()})", "")
    #     status.SetTextColor(ROOT.kGreen + 2 if okErf else ROOT.kRed)
        

    #     return legend
    
    def MakeLegend(self, fitRes, fitResErf, fitFunc, fitFuncErf):
        ok    = fitRes.IsValid()    and fitRes.Status() == 0
        okErf = fitResErf.IsValid() and fitResErf.Status() == 0

        leg = ROOT.TLegend(0.5, 0.2, 0.98, 0.9)
        leg.SetNColumns(2)
        leg.SetBorderSize(1)
        leg.SetFillColorAlpha(0, 0)
        leg.SetTextSize(0.07)          # 0 = auto-size text to fit the box

        # Row 0: column titles
        leg.AddEntry(ROOT.nullptr, "Gaus + exp", "")
        leg.AddEntry(ROOT.nullptr, "Gaus + exp + erf", "")

        # Row 1: mean (line sample uses each fit's own colour)
        leg.AddEntry(fitFunc,    f"#mu = {fitRes.Parameter(1):.1f} #pm {fitRes.ParError(1):.1f}", "l")
        leg.AddEntry(fitFuncErf, f"#mu = {fitResErf.Parameter(1):.1f} #pm {fitResErf.ParError(1):.1f}", "l")

        # Row 2: chi2/ndf
        leg.AddEntry(ROOT.nullptr, f"#chi^{{2}}/ndf = {fitRes.Chi2():.1f} / {fitRes.Ndf()}", "")
        leg.AddEntry(ROOT.nullptr, f"#chi^{{2}}/ndf = {fitResErf.Chi2():.1f} / {fitResErf.Ndf()}", "")

        # Row 3: p-value
        leg.AddEntry(ROOT.nullptr, f"p-value = {fitRes.Prob():.3f}", "")
        leg.AddEntry(ROOT.nullptr, f"p-value = {fitResErf.Prob():.3f}", "")

        # Row 4: status
        s1 = leg.AddEntry(ROOT.nullptr, "Fit OK" if ok else f"Fit FAILED (status {fitRes.Status()})", "")
        s1.SetTextColor(ROOT.kGreen + 2 if ok else ROOT.kRed)
        s2 = leg.AddEntry(ROOT.nullptr, "Fit OK" if okErf else f"Fit FAILED (status {fitResErf.Status()})", "")
        s2.SetTextColor(ROOT.kGreen + 2 if okErf else ROOT.kRed)

        return leg
        
class MeanPositionData:
    def __init__(self, norm, erf,saveFP):
        
        self.saveFP = saveFP
        #Store the means and their errors
        self.positions = []
        # self.positionsErr = [0.25 for i in self.positions]
        
        self.norm = norm
        self.erf = erf
        self.muL = []
        self.muLErr = []
        
        self.muR = []
        self.muRErr = []
        
        #Store the summed means
        self.summedMeans = []
        self.summedMeansErr = []
        
        self.sqrtProd = []
        self.sqrtProdErr = []
        
        #Store the log ratio
        self.logRatio = []
        self.logRatioErr = []
        
    def addData(self,muL,muLErr,muR,muRErr,position):
        
        self.positions.append(position)
        
        #Add the means and their errors        
        self.muL.append(muL)
        self.muLErr.append(muLErr)
        
        self.muR.append(muR)
        self.muRErr.append(muRErr)
        
        #Sum the data
        self.summedMeans.append((muL+muR)/2)
        self.summedMeansErr.append(muLErr + muRErr)
        
        self.sqrtProd.append(np.sqrt(muL*muR))
        self.sqrtProdErr.append(0.5*np.sqrt((muR/muL)*muLErr**2 + (muL/muR)*muRErr**2))
        
        #Calculate the log Ratio:
        self.logRatio.append(np.log(muR/muL))
        self.logRatioErr.append(np.sqrt((muLErr/muL)**2 + (muRErr/muR)**2))


    def plotResults(self):
        
        self.positionsErr = [0.25 for i in self.positions]
        
        if self.norm and not self.erf:
            resultName = "Position_test_no_outer_sleeve_mean_norm"
        elif self.erf and not self.norm:
            resultName = "Position_test_no_outer_sleeve_mean_erf"
        elif self.norm and self.erf:
            resultName = "Position_test_no_outer_sleeve_mean_norm_erf"
        else:
            resultName = "Position_test_no_outer_sleeve_mean"
            
            
        #Make the root TGraphErrors for the final results.     
        self.meanLPos = makeGraph(self.positions,self.muL,self.positionsErr,self.muLErr) #= ROOT.TGraphErrors(len(self.positions), self.positions, self.muL,self.positionsErr, self.muLErr)
        self.meanLPos.GetXaxis().SetTitle("Position (cm)")
        if self.norm:
            self.meanLPos.GetYaxis().SetTitle("mean position (normalized to position 0)")
        else:
            self.meanLPos.GetYaxis().SetTitle("mean position (ADC)")
        self.meanLPos.SetTitle("Mean Position Left PMT")
        self.meanLPos.SetMarkerSize(2.0)
        self.meanLPos.SetMarkerColor(ROOT.kBlue)
        self.meanLPos.SetMarkerStyle(22)
        
        
        
        self.meanRPos = makeGraph(self.positions,self.muR,self.positionsErr,self.muRErr) #ROOT.TGraphErrors(len(self.positions), self.positions, self.muR,self.positionsErr, self.muRErr)
        self.meanRPos.GetXaxis().SetTitle("Position (cm)")
        if self.norm:
            self.meanRPos.GetYaxis().SetTitle("mean position (normalized to position 0)")
        else:
            self.meanRPos.GetYaxis().SetTitle("mean position (ADC)")
        self.meanRPos.SetTitle("Mean Position Right PMT")
        self.meanRPos.SetMarkerSize(2.0)
        self.meanRPos.SetMarkerColor(ROOT.kRed)
        self.meanRPos.SetMarkerStyle(23)
        
        
        
        self.summedPosition = makeGraph(self.positions,self.summedMeans,self.positionsErr,self.summedMeansErr) #ROOT.TGraphErrors(len(self.positions), self.positions, self.summedMeans,self.positionsErr, self.summedMeansErr)
        self.summedPosition.GetXaxis().SetTitle("Position (cm)")
        if self.norm:
            self.summedPosition.GetYaxis().SetTitle("mean position (normalized to position 0)")
        else:
            self.summedPosition.GetYaxis().SetTitle("mean position (ADC)")
        self.summedPosition.SetTitle("Mean Position summed PMTs")
        self.summedPosition.SetMarkerSize(2.0)
        self.summedPosition.SetMarkerColor(ROOT.kGreen)
        self.summedPosition.SetMarkerStyle(21)
        
        
        
        self.sqrtProdPosition = makeGraph(self.positions,self.sqrtProd,self.positionsErr,self.sqrtProdErr) #ROOT.TGraphErrors(len(self.positions), self.positions, self.sqrtProd,self.positionsErr, self.sqrtProdErr)
        self.sqrtProdPosition.GetXaxis().SetTitle("Position (cm)")
        if self.norm:
            self.sqrtProdPosition.GetYaxis().SetTitle("mean position (normalized to position 0)")
        else:
            self.sqrtProdPosition.GetYaxis().SetTitle("mean position (ADC)")
        self.sqrtProdPosition.SetTitle("square root product")
        self.sqrtProdPosition.SetMarkerSize(2.0)
        self.sqrtProdPosition.SetMarkerColor(ROOT.kBlack)
        self.sqrtProdPosition.SetMarkerStyle(20)
        
        
        
        self.logRatioPosition = makeGraph(self.positions,self.logRatio,self.positionsErr,self.logRatioErr) #ROOT.TGraphErrors(len(self.positions),self.positions, self.logRatio, self.positionsErr,self.logRatioErr)  
        self.logRatioPosition.GetXaxis().SetTitle("Position (cm)")
        self.logRatioPosition.GetYaxis().SetTitle("ln(S2/S1)")
        self.logRatioPosition.SetTitle("Log Ratio")
        self.logRatioPosition.SetMarkerSize(2.0)
        self.logRatioPosition.SetMarkerColor(ROOT.kBlue)
        self.logRatioPosition.SetMarkerStyle(21)  
            
            
        #Make the results plot
        self.Canvas = ROOT.TCanvas(resultName,resultName, 1400,800)
        
        self.Canvas.Divide(1,2)
        
        pad = self.Canvas.cd(1)
        canvasPadding(pad,0.1,0.2,0.10,0.15)
        
        self.meanLPos.Draw("AP")
        self.meanRPos.Draw("P")
        self.summedPosition.Draw("P")
        self.sqrtProdPosition.Draw("P")
        
        self.leg = ROOT.TLegend(0.81, 0.60, 0.99, 0.88)
        self.leg.SetBorderSize(1)
        # self.leg.SetFillStyle(0)
        self.leg.SetTextSize(0.04)
        self.leg.AddEntry(self.meanLPos, "Left PMT", "pe")
        self.leg.AddEntry(self.meanRPos, "Right PMT", "pe")
        self.leg.AddEntry(self.summedPosition, "(#mu_{L} + #mu_{R}) / 2", "pe")
        self.leg.AddEntry(self.sqrtProdPosition, "#sqrt{#mu_{L}#mu_{R}}","pe")
        self.leg.Draw()
        
        
        self.Canvas.Update()
        
        pad = self.Canvas.cd(2)
        canvasPadding(pad,0.1,0.1,0.10,0.15)
        
        self.logRatioPosition.Draw("AP")
        self.Canvas.Update()
        
        self.Canvas.SaveAs(str(self.saveFP / f"{resultName}.pdf"))
        self.Canvas.SaveAs(str(self.saveFP / f"{resultName}.png"))
        
        
        
    
    
        
# def Fit1DData(df):


def makeResultPlots(data,centerData,saveFP):
    
    positions = np.linspace(-5,5,11)
    
    
    muCenterL = centerData.resL.Parameter(1)
    muCenterR = centerData.resR.Parameter(1)
    
    muCenterLErr = centerData.resL.ParError(1)
    muCenterRErr = centerData.resR.ParError(1)
    
    muCenterLErf = centerData.resLErf.Parameter(1)
    muCenterRErf = centerData.resRErf.Parameter(1)
    
    muCenterLErfErr = centerData.resLErf.ParError(1)
    muCenterRErfErr = centerData.resRErf.ParError(1)
    
    
    
    MeanData = MeanPositionData(norm = False,erf = False,saveFP = saveFP)
    MeanDataNorm = MeanPositionData(norm = True,erf = False,saveFP= saveFP)
    MeanDataErf = MeanPositionData(norm = False,erf = True,saveFP= saveFP)
    MeanDataNormErf = MeanPositionData(norm = True,erf = True,saveFP= saveFP)
    
    
    for d in data:
        #Grab the mu values and errors from each detector and the two fits. 
        muL = d.resL.Parameter(1)
        muLErr = d.resL.ParError(1)
        
        muR = d.resR.Parameter(1)
        muRErr = d.resR.ParError(1)
        
        muLErf = d.resLErf.Parameter(1)
        muLErfErr = d.resLErf.ParError(1)
        
        muRErf = d.resRErf.Parameter(1)
        muRErfErr = d.resRErf.ParError(1)
        
        position = d.position
        
        #Add the data to the class: MeanPositionData
        MeanData.addData(muL,muLErr,muR,muRErr,position)
        
        MeanDataNorm.addData(muL/muCenterL, 
                             muL/muCenterL * ((muLErr/muL)**2 + (muCenterLErr/muCenterL)**2),
                             muR/muCenterR,
                             muR/muCenterR * ((muRErr/muR)**2 + (muCenterLErr/muCenterL)**2),
                             position)
        
        MeanDataErf.addData(muLErf,muLErfErr,muRErf,muRErfErr,position)
        
        MeanDataNormErf.addData(muLErf/muCenterLErf, 
                             muLErf/muCenterLErf * ((muLErfErr/muLErf)**2 + (muCenterLErfErr/muCenterLErf)**2),
                             muRErf/muCenterRErf,
                             muRErf/muCenterRErf * ((muRErfErr/muRErf)**2 + (muCenterLErfErr/muCenterLErf)**2),
                             position)
        

        

    MeanData.plotResults()
    MeanDataNorm.plotResults()
    MeanDataErf.plotResults()
    MeanDataNormErf.plotResults()
    

def makeGraph(x, y, ex, ey):
    x, y, ex, ey = [np.asarray(a, dtype=np.float64) for a in (x, y, ex, ey)]
    return ROOT.TGraphErrors(len(x), x, y, ex, ey)

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

def ReadInNaIScale(file):
    with open(file) as f:
        next(f)
        
        for line in f:
            NaIScale = float(line.split(',')[4])
            
    return NaIScale 

def scaleNaIIntegrals(df,scale):
    df = df.Define("Energy_2_keV", f"Energy_2 / {scale}")
    
    return df

def canvasPadding(pad,L,R,T,B):
    pad.SetLeftMargin(L)
    pad.SetRightMargin(R)   # extra room for COLZ palette
    pad.SetTopMargin(T)
    pad.SetBottomMargin(B)
def canvasPaddingLegend(pad):
    pad.SetLeftMargin(0.05)
    pad.SetRightMargin(0.52)   # extra room for COLZ palette
    pad.SetTopMargin(0.10)
    pad.SetBottomMargin(0.15)
    

RootDirectory = Path("/home/nick/PhD/KDK+/Large_LSC_testing/Position_test_no_outer_sleeve/")
NaIScaleFile = Path("/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/Loaded_BNL_L_LSC_Small_NaI_LY_coinc_2/RAW/SDataR_Loaded_BNL_L_LSC_Small_NaI_LY_coinc_2/figures/SDataR_Loaded_BNL_L_LSC_Small_NaI_LY_coinc_2_scale_factors.txt")

settingsFilePath = Path("/home/nick/PhD/KDK+/Large_LSC_testing/Position_test_no_outer_sleeve/2026_11_03/Position_test_no_sleeve_position_0_NaI_lead_shield/settings.xml")

savefilepath = RootDirectory / 'results'

savefilepath.mkdir(parents = True, exist_ok = True)

NaIScale = ReadInNaIScale(NaIScaleFile)

#Read in the channel labels
files = sorted(RootDirectory.rglob("*_sorted_coinc*.root"))

# chLabels = ReadInChannelNames(settingsFilePath)

channels = [0,1,2]

fitWindow = {'-5': [[1800,3600],[900,1900]],
             '-4': [[1600,3600],[800,2200]],
             '-3': [[1400,3200],[1000,2200]],
             '-2': [[1400,3000],[1000,2500]],
             '-1': [[1200,2800],[1400,2600]],
             '0': [[1200,2500],[1500,2500]],
             '1': [[1400,2600],[1200,2800]],
             '2': [[1000,2500],[1400,3000]],
             '3': [[1000,2200],[1500,3300]],
             '4': [[800,2200],[1500,3500]],
             '5': [[800,2000],[1500,3500]]}

Hist1DFitsCanvas = ROOT.TCanvas("Hist_1D_Position_fits", "Position Histograms",25600,1600)
ncols = 11
Hist1DFitsCanvas.Divide(ncols,2)

PositionLogRatio2D = ROOT.TH2D(
    "Log_ratio_vs_position","Log Ratio Vs Position 2D Hist", 21,-5.25,5.25,50,-2,2
)


pData = []
for i,file in enumerate(files):
    pData.append(positionData(file,channels,ncols,fitWindow,savefilepath))  
    
    if pData[-1].position == 0 and pData[-1].trialNumber == 0:
        centralData = pData[-1]
    
    pData[-1].plot1DHist(Hist1DFitsCanvas)
    pData[-1].plotLogRatioHist()
    
    pos = pData[-1].position
    
    h1 = pData[-1].logRatioHist.GetValue()
    totalEntries = h1.GetEntries()
    for j in range(1,h1.GetNbinsX() + 1):
        nCounts = h1.GetBinContent(j)
        
        if nCounts > 0: 
            PositionLogRatio2D.Fill(pos,h1.GetBinCenter(j), nCounts/totalEntries)
    # for logRatio in pData[-1].logRatioHist:
    #     PositionLogRatio2D.Fill(pos,logRatio)
    
makeResultPlots(pData,centralData,savefilepath)
    
Hist1DFitsCanvas.SaveAs(str(savefilepath / 'Position_1D_hist_fits.pdf'))
Hist1DFitsCanvas.SaveAs(str(savefilepath / 'Position_1D_hist_fits.png'))


PositionLogRatio2D.GetXaxis().SetTitle('Position (cm)')
PositionLogRatio2D.GetYaxis().SetTitle('ln#left(#frac{S_{R}}{S_{L}}#right)')
PositionLogRatio2D.GetZaxis().SetTitle('Normalized Counts')
PositionLogRatio2D.GetZaxis().SetTitleOffset(1.5)
PositionLogRatio2D.SetStats(0)

logRatioPositionCanvas = ROOT.TCanvas(f"Log_ratio_Position_canvas", "Log Ratio Vs Position Canvas", 1400,800)

pad = logRatioPositionCanvas.cd(1)
# pad.SetLogz(1)
canvasPadding(pad,0.15,0.25,0.1,0.1)

PositionLogRatio2D.Draw("COLZ")

logRatioPositionCanvas.Update()
logRatioPositionCanvas.SaveAs(str(savefilepath / 'Log_Ratio_vs_position.png'))
logRatioPositionCanvas.SaveAs(str(savefilepath / 'Log_Ratio_vs_position.pdf'))

