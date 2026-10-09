import numpy as np
import ROOT
from pathlib import Path
import pandas as pd
import scipy.stats as Stats

# Colour cycle that uses distinct colours in root
# distinctColours = [
#     ROOT.kBlue, ROOT.kTeal+3, ROOT.kGreen+2, ROOT.kMagenta,
#     ROOT.kBlack, ROOT.kCyan+1, ROOT.kGreen+2, ROOT.kSpring+4,
#     ROOT.kPink+9, ROOT.kAzure+7, ROOT.kYellow+2, ROOT.kTeal+3,
#     ROOT.kGray+2, ROOT.kBlack
# ]

# distinctColours = [
#     # ROOT.TColor.GetColor("#332288"),
#     ROOT.kBlue,
#     ROOT.TColor.GetColor("#117733"),
#     ROOT.TColor.GetColor("#44AA99"),
#     ROOT.TColor.GetColor("#88CCEE"),
#     ROOT.TColor.GetColor("#DDCC77")
# ]

distinctColours = [
    ROOT.TColor.GetColor("#648FFF"),
    ROOT.TColor.GetColor("#785EF0"),
    ROOT.TColor.GetColor("#DC267F"),
    ROOT.TColor.GetColor("#FE6100"),
    ROOT.TColor.GetColor("#FFB000")
]
class sourceData:
    def __init__(self,source,file,fitWindows,nbins,channels):
        self.source = source
        self.file = file
        self.fitWindows = fitWindows
        self.fitFunc = None
        self.Ch = channels
        filetype = str(self.file).split(".")[-1]
       

        self.NaI1Data = ROOT.TH1D(f"NaI_1_Hist_data_{source}",f"NaI 1 Hist Data {source}",nbins,0,4000)
        self.NaI2Data = ROOT.TH1D(f"NaI_2_Hist_data_{source}",f"NaI 2 Hist Data {source}",nbins,0,4000)
        self.NaI3Data = ROOT.TH1D(f"NaI_3_Hist_data_{source}",f"NaI 3 Hist Data {source}",nbins,0,4000)
        self.NaI4Data = ROOT.TH1D(f"NaI_4_Hist_data_{source}",f"NaI 4 Hist Data {source}",nbins,0,4000)
        

        self.ChDict = {f'{self.Ch[0]}': [self.NaI1Data,self.fitWindows[0]],
                    f'{self.Ch[1]}': [self.NaI2Data,self.fitWindows[1]],
                    f'{self.Ch[2]}': [self.NaI3Data,self.fitWindows[2]],
                    f'{self.Ch[3]}': [self.NaI4Data,self.fitWindows[3]]}

            
    def ReadInData(self):
        #Reads in the data from the file and then saves it to the appropriate histogram
        df = pd.read_csv(self.file,sep = ";",skiprows = 1, header = None,
                         usecols = [1,3],names = ["ch", "Energy"], dtype = {"ch": "int32","Energy": "float64"},engine = "pyarrow")
        
        for ch,hist in self.ChDict.items():
            chValue = int(ch)
            Energy = df.loc[df["ch"] == chValue, "Energy"].to_numpy()
            hist[0].FillN(len(Energy),Energy,np.ones(len(Energy)))
            
    def ReadInDataROOT(self):
        

        self.RootFile = ROOT.TFile.Open(str(self.file))
        self.RootFile.ls()
        self.tree = self.RootFile.Get("Data_R")
        self.tree.Print()
        
        self.df = ROOT.RDataFrame(self.tree)
        
        
        results = {}
        for ch,data in self.ChDict.items():

            results[ch] = self.df.Filter(f"Channel == {ch}").Histo1D((f"Energy_hist_ch_{ch}", f"Ch {ch};Energy (ADC);Counts",100,0,4000), "Energy")
            # data[0] = h
            
        for ch, data in self.ChDict.items():
            h = results[ch].GetValue().Clone()
            h.SetDirectory(0)          # detach from the file so it survives file closure
            data[0] = h
            
        # self.h = self.df.Filter("Channel == 2").Histo1D(("Energy_hist_ch0", "Ch0;Energy;Counts",100,0,4000), "Energy")
        # self.c = ROOT.TCanvas("Test_canvas","Test Canvas", 1600,800)
        # self.h.Draw("Same")
        # self.c.Update()
        
    def PlotTotalData(self,canvas,colourInd,totalHists,legends):
        
        
        for i,ch in enumerate(self.Ch):
            pad = canvas.cd(i+1)
            pad.SetLogy(1)
            
            originalHist = self.ChDict[f"{ch}"][0]
            
            hist = originalHist.Clone(f"{originalHist.GetName()}_tota")
            hist.SetDirectory(0)
            hist.SetLineColor(distinctColours[colourInd])
            hist.SetStats(0)
            hist.SetLineWidth(3)
            
            legends[i].AddEntry(hist,self.source,"l")
            
            totalHists[i].append(hist)
            hist.Draw("SAME")
            
            
            frame = totalHists[i][0]
            yMax = max(x.GetBinContent(x.GetMaximumBin()) for x in totalHists[i])
            yMin = min(
                (x.GetBinContent(b) for x in totalHists[i]
                for b in range(1, x.GetNbinsX() + 1) if x.GetBinContent(b) > 0),
                default=1)
            frame.SetMaximum(3 * yMax)     # extra headroom since the axis is log
            frame.SetMinimum(0.5 * yMin)
            
            canvas.Update()
        
        
        
            
        
    def FitData(self,saveFP):    
        
        self.DataCanvas = ROOT.TCanvas(f"{self.source}_canvas",f"{self.source} Canvas", 3200,1600)
        self.DataCanvas.Divide(2,2)
        
        self.legends = []
        padInd = 1  
        for key,value in self.ChDict.items():
            
            hist = value[0]
            fitWindow = value[1]
            hist.SetStats(0)
            legend = ROOT.TLegend(0.1, 0.1, 0.9, 0.45)  #x1,y1,x2,y2
            legend.SetFillStyle(0)     # transparent background
            legend.SetTextSize(0.03)
            legend.SetNColumns(3)
            legend.SetHeader("Fit Results", "C")
            if self.source == "Co-60":
                # fitFormula = ("[0]*TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1) + [6]*TMath::Expo(x)")
                self.fitFunc =  ROOT.TF1(f"{self.source}_{key}_double_gaus","gaus(0) + gaus(3) + expo(6)", fitWindow[0], fitWindow[1])
                ampGuess = 5000
                sigGuess = fitWindow[1]-fitWindow[0]
                muGuess = (fitWindow[1]+fitWindow[0])/3
                muGuess2 = 2*(fitWindow[1]+fitWindow[0])/3
                
                expAmpGuess = 5
                tau = 1e-10
                self.fitFunc.SetParameters(ampGuess,muGuess,sigGuess,ampGuess,muGuess2,sigGuess,expAmpGuess,tau)
                self.fitFunc.SetLineColor(ROOT.kRed)
                self.fitFunc.SetLineWidth(3)
                
                
                # self.legends.append(legend)
                
            elif self.source == "Na-22" or self.source == "Zn-65":
                self.fitFunc1 =  ROOT.TF1(f"{self.source}_{key}_double_gaus_1","gaus(0) + expo(3)", fitWindow[0], fitWindow[1])
                self.fitFunc2 =  ROOT.TF1(f"{self.source}_{key}_double_gaus_2","gaus(0) + expo(3)", fitWindow[2], fitWindow[3])
                # self.fitFunc3 =  ROOT.TF1(f"{self.source}_{key}_double_gaus_3","gaus(0) + expo(3)", fitWindow[4], fitWindow[5])
                ampGuess = 1000
                expAmpGuess = 5
                sigGuess1 = fitWindow[1]-fitWindow[0]
                muGuess1 = (fitWindow[1]+fitWindow[0])/2
                
                sigGuess2 = fitWindow[3]-fitWindow[2]
                muGuess2 = (fitWindow[3]+fitWindow[2])/2
                
                # sigGuess3 = fitWindow[5]-fitWindow[4]
                # muGuess3 = (fitWindow[5]+fitWindow[4])/2
                tau = 1e-10
                self.fitFunc1.SetParameters(ampGuess,muGuess1,sigGuess1,expAmpGuess,tau)
                self.fitFunc2.SetParameters(ampGuess,muGuess2,sigGuess2,expAmpGuess,tau)
                # self.fitFunc3.SetParameters(ampGuess,muGuess3,sigGuess3,expAmpGuess,tau)
                
                self.fitFunc1.SetLineColor(ROOT.kRed)
                self.fitFunc2.SetLineColor(ROOT.kGreen+2)
                # self.fitFunc3.SetLineColor(ROOT.kRed)
                
                self.fitFunc1.SetLineWidth(3)
                self.fitFunc2.SetLineWidth(3)
                # self.fitFunc3.SetLineWidth(5)
                
            else:
                self.fitFunc =  ROOT.TF1(f"{self.source}_{key}_double_gaus","gaus(0) + expo(3)", fitWindow[0], fitWindow[1])
                ampGuess = 5000
                sigGuess = fitWindow[1]-fitWindow[0]
                muGuess = (fitWindow[1]+fitWindow[0])/2
                expAmpGuess = 5
                tau = 1e-10
                self.fitFunc.SetParameters(ampGuess,muGuess,sigGuess,expAmpGuess,tau)
                self.fitFunc.SetLineColor(ROOT.kRed)
                self.fitFunc.SetLineWidth(3)
                
            
            
            
            pad = self.DataCanvas.cd(padInd)
            pad.SetLogy(1)
            
            hist.Draw("Same")
            
            if self.fitFunc is not None:
                hist.Fit(self.fitFunc,"R")
                self.fitFunc.Draw("SAME")
                
                if self.source == "Co-60":
                    mu1 = self.fitFunc.GetParameter(1)
                    mu1Err = self.fitFunc.GetParError(1)
                    sig1 = self.fitFunc.GetParameter(2)
                    sig1Err = self.fitFunc.GetParError(2)
                    mu2 = self.fitFunc.GetParameter(4)
                    mu2Err = self.fitFunc.GetParError(4)
                    sig2 = self.fitFunc.GetParameter(5)
                    sig2Err = self.fitFunc.GetParError(5)
                    chi2 = self.fitFunc.GetChisquare()
                    ndof = self.fitFunc.GetNDF()
                    pVal = Stats.chi2.sf(chi2,ndof)
                    
                    legend.AddEntry(self.fitFunc, f"mu1 = {mu1:.2f} #pm {mu1Err:.2f}", "l")
                    legend.AddEntry(self.fitFunc, f"sig1 = {sig1:.2f} #pm {sig1Err:.2f}", "l")
                    legend.AddEntry(self.fitFunc, f"mu2 = {mu2:.2f} #pm {mu2Err:.2f}", "l")
                    legend.AddEntry(self.fitFunc, f"sig2 = {sig2:.2f} #pm {sig2Err:.2f}", "l")
                    legend.AddEntry(self.fitFunc, f"#chi^{{2}} / ndof = {chi2:.2f} / {ndof}", "l")
                    legend.AddEntry(self.fitFunc, f"p-value = {pVal:.4f}", "l")
                
                    self.legends.append(legend)
                else:
                    mu = self.fitFunc.GetParameter(1)
                    muErr = self.fitFunc.GetParError(1)
                    sig = self.fitFunc.GetParameter(2)
                    sigErr = self.fitFunc.GetParError(2)
                    chi2 = self.fitFunc.GetChisquare()
                    ndof = self.fitFunc.GetNDF()
                    pVal = Stats.chi2.sf(chi2,ndof)
                    
                    legend.AddEntry(self.fitFunc, f"mu = {mu:.2f} #pm {muErr:.2f}", "l")
                    legend.AddEntry(self.fitFunc, f"sig = {sig:.2f} #pm {sigErr:.2f}", "l")
                    legend.AddEntry(self.fitFunc, f"#chi^{{2}} / ndof = {chi2:.2f} / {ndof}", "l")
                    legend.AddEntry(self.fitFunc, f"p-value = {pVal:.4f}", "l")
                
                    self.legends.append(legend)
            else:
                hist.Fit(self.fitFunc1,"R")
                hist.Fit(self.fitFunc2,"R")
                # hist.Fit(self.fitFunc3,"R")
                
                
                self.fitFunc1.Draw("SAME")
                self.fitFunc2.Draw("SAME")
                
                mu1 = self.fitFunc1.GetParameter(1)
                mu1Err = self.fitFunc1.GetParError(1)
                sig1 = self.fitFunc1.GetParameter(2)
                sig1Err = self.fitFunc1.GetParError(2)
                mu2 = self.fitFunc2.GetParameter(4)
                mu2Err = self.fitFunc2.GetParError(4)
                sig2 = self.fitFunc2.GetParameter(5)
                sig2Err = self.fitFunc2.GetParError(5)
                chi2_1 = self.fitFunc1.GetChisquare()
                ndof1 = self.fitFunc1.GetNDF()
                pVal1 = Stats.chi2.sf(chi2_1,ndof1)
                
                chi2_2 = self.fitFunc2.GetChisquare()
                ndof2 = self.fitFunc2.GetNDF()
                pVal2 = Stats.chi2.sf(chi2_2,ndof2)
                
                legend.AddEntry(self.fitFunc1, f"mu1 = {mu1:.2f} #pm {mu1Err:.2f}", "l")
                legend.AddEntry(self.fitFunc1, f"sig1 = {sig1:.2f} #pm {sig1Err:.2f}", "l")

                legend.AddEntry(self.fitFunc1, f"#chi^{{2}} / ndof = {chi2_1:.2f} / {ndof1}", "l")
                legend.AddEntry(self.fitFunc1, f"p-value = {pVal1:.4f}", "l")
                
                legend.AddEntry(self.fitFunc2, f"mu2 = {mu2:.2f} #pm {mu2Err:.2f}", "l")
                legend.AddEntry(self.fitFunc2, f"sig2 = {sig2:.2f} #pm {sig2Err:.2f}", "l")
                
                legend.AddEntry(self.fitFunc2, f"#chi^{{2}} / ndof = {chi2_2:.2f} / {ndof2}", "l")
                legend.AddEntry(self.fitFunc2, f"p-value = {pVal2:.4f}", "l")
                self.legends.append(legend)
                # self.fitFunc3.Draw("SAME")
            
            try:
                self.legends[-1].Draw()
            except:
                pass
            self.DataCanvas.Update()

            
            saveFileName = saveFP / f"Annulus_cal_{self.source}_fit.png"
            self.DataCanvas.SaveAs(str(saveFileName))
            
            padInd +=1
        
                
                
def ReadInInitFile(file,nbins,channels):
    
    
    with open(file) as f:
        
        lines = f.readlines()
        
        saveFilePath = lines[0].split("\t")[1]
        
        data = []
        
        for i in range(5):
            source = lines[i+1].split("\t")[0]
            sourceFile = lines[i+1].split("\n")[0].split("\t")[1]
            fitWindowsStr = lines[i+6].split("\n")[0].split("\t")[1:]
            
            fitWindowsInt = [[int(x) for x in item.split(",")] for item in fitWindowsStr]

            data.append(sourceData(source,sourceFile,fitWindowsInt,nbins,channels))
            
    return data

        
        
            
# initFile = Path("/home/nick/PhD/KDK+/Annulus_Compton_scatter_V1/NaI_Annulus_calibration/CFD_lower_HV_Inside_Source.txt")   
initFile = Path("/home/nick/PhD/KDK+/Pilot_exp/Annulus_cal_data/Annulus_cal_init_v2.txt")
saveFP = initFile.parent / "Annulus_cal_figures"

saveFP.mkdir(parents= True, exist_ok = True)

fileType = "ROOT"
Channels = [2,3,4,5]

nbins = 100

sourceData = ReadInInitFile(initFile,nbins,Channels)

totalDataCanvas = ROOT.TCanvas("Total_data_canvas","Total Data", 3200,1600)
totalDataCanvas.Divide(2,2)
colourInd = 0
totalHists = [[] for _ in range(4)] 

legends = []
for i in range(4):
    legend = ROOT.TLegend(0.25, 0.15, 0.75, 0.35)   # NDC coords: x1,y1,x2,y2 - sits in the 0.20 right margin
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)     # transparent background
    legend.SetTextSize(0.06)
    legend.SetNColumns(3)
    
    legends.append(legend)

for data in sourceData:
    print(f"Reading in {data.source} data:")
    if fileType == "ROOT":
        data.ReadInDataROOT()
    else:
        data.ReadInData()
        
        
    data.PlotTotalData(totalDataCanvas,colourInd,totalHists,legends)
    colourInd+=1
        
    data.FitData(saveFP)
    
for i in range(4):
    pad = totalDataCanvas.cd(i+1)
    legends[i].Draw()
    
   
totalDataCanvas.SaveAs(str(saveFP / "Annulus_cal_all_sources.png"))
 
        

        
        