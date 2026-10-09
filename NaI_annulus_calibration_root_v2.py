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

SQRT2 = np.sqrt(2.0)

def SkewExpr(p):
    """TF1 expression for a skew-normal using parameters p..p+3 = [area, xi, omega, alpha]"""
    return (f"[{p}]*TMath::Gaus(x,[{p+1}],[{p+2}],1)"
            f"*(1+TMath::Erf([{p+3}]*(x-[{p+1}])/([{p+2}]*{SQRT2})))")

def MakeSkewFunc(name, lo, hi, nPeaks):
    expr = " + ".join(SkewExpr(4*i) for i in range(nPeaks)) + f" + expo({4*nPeaks})"
    return ROOT.TF1(name, expr, lo, hi)

def InitSkewPeak(f, p, area, mu, sigma, alpha=0.5):
    f.SetParameter(p,   area)
    f.SetParameter(p+1, mu)
    f.SetParameter(p+2, sigma)
    f.SetParameter(p+3, alpha)
    f.SetParLimits(p+2, 1e-3, 1e4)   # omega > 0
    f.SetParLimits(p+3, -20, 20)     # stops alpha running away when the peak is ~symmetric

def SkewMoments(xi, om, al):
    d = al / np.sqrt(1 + al**2)
    mean = xi + om * d * np.sqrt(2/np.pi)
    sd = om * np.sqrt(1 - 2*d**2/np.pi)
    return mean, sd

def SkewMomentsErr(res, p):
    """Mean, sd and their errors for the skew peak starting at parameter index p,
    propagated through the fit covariance matrix (numerical Jacobian)."""
    idx = [p+1, p+2, p+3]
    par = np.array([res.Parameter(i) for i in idx])
    cov = np.array([[res.CovMatrix(i, j) for j in idx] for i in idx])

    def g(v):
        return np.array(SkewMoments(*v))

    J = np.zeros((2, 3))
    for k in range(3):
        step = 1e-6 * max(abs(par[k]), 1.0)
        up, dn = par.copy(), par.copy()
        up[k] += step
        dn[k] -= step
        J[:, k] = (g(up) - g(dn)) / (2*step)

    outCov = J @ cov @ J.T
    mean, sd = g(par)
    return mean, np.sqrt(outCov[0, 0]), sd, np.sqrt(outCov[1, 1])


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
            
            hist = originalHist.Clone(f"{originalHist.GetName()}_total")
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
            
    def FitData(self, saveFP):
        self.DataCanvas = ROOT.TCanvas(f"{self.source}_canvas", f"{self.source} Canvas", 3200, 1600)
        self.DataCanvas.Divide(2, 2)
        self.legends = []
        self.funcs = []          # keep Python references so ROOT doesn't lose them
        self.results = {}        # ch -> list of (mean, meanErr, sd, sdErr)

        colours = [ROOT.kRed, ROOT.kGreen+2]

        for padInd, (key, (hist, w)) in enumerate(self.ChDict.items(), start=1):
            hist.SetStats(0)
            pad = self.DataCanvas.cd(padInd)
            pad.SetLogy(1)
            hist.Draw("SAME")

            legend = ROOT.TLegend(0.1, 0.1, 0.9, 0.45)
            legend.SetFillStyle(0)
            legend.SetTextSize(0.03)
            legend.SetNColumns(3)
            legend.SetHeader("Fit Results", "C")

            # Each entry: (lo, hi, nPeaks, [peak centre guesses])
            if self.source == "Co-60":
                fits = [(w[0], w[1], 2, [(w[0]+w[1])/3, 2*(w[0]+w[1])/3])]
            elif self.source in ("Na-22", "Zn-65"):
                fits = [(w[0], w[1], 1, [(w[0]+w[1])/2]),
                        (w[2], w[3], 1, [(w[2]+w[3])/2])]
            else:
                fits = [(w[0], w[1], 1, [(w[0]+w[1])/2])]

            self.results[key] = []
            for n, (lo, hi, nPeaks, centres) in enumerate(fits):
                f = MakeSkewFunc(f"{self.source}_{key}_skew_{n}", lo, hi, nPeaks)

                area = hist.Integral(hist.FindBin(lo), hist.FindBin(hi)) * hist.GetBinWidth(1) / nPeaks
                sigma = (hi - lo) / 6          # your original guess (full window width) was very wide
                for i, mu in enumerate(centres):
                    InitSkewPeak(f, 4*i, area, mu, sigma)
                f.SetParameter(4*nPeaks, 5)        # expo constant (log-amplitude)
                f.SetParameter(4*nPeaks + 1, 1e-10)  # expo slope

                f.SetLineColor(colours[n % 2])
                f.SetLineWidth(3)

                res = hist.Fit(f, "RS+")          # S -> returns TFitResult, + -> keep earlier fit functions
                f.Draw("SAME")
                self.funcs.append(f)

                chi2, ndof = f.GetChisquare(), f.GetNDF()
                pVal = Stats.chi2.sf(chi2, ndof) if ndof > 0 else float("nan")

                for i in range(nPeaks):
                    mean, meanErr, sd, sdErr = SkewMomentsErr(res, 4*i)
                    self.results[key].append((mean, meanErr, sd, sdErr))
                    tag = f"{n+1}" if len(fits) > 1 else (f"{i+1}" if nPeaks > 1 else "")
                    legend.AddEntry(f, f"mean{tag} = {mean:.2f} #pm {meanErr:.2f}", "l")
                    legend.AddEntry(f, f"sd{tag} = {sd:.2f} #pm {sdErr:.2f}", "l")
                    legend.AddEntry(f, f"#alpha{tag} = {f.GetParameter(4*i+3):.2f}", "l")

                legend.AddEntry(f, f"#chi^{{2}} / ndof = {chi2:.2f} / {ndof}", "l")
                legend.AddEntry(f, f"p-value = {pVal:.4f}", "l")

            self.legends.append(legend)
            legend.Draw()
            self.DataCanvas.Update()

        self.DataCanvas.SaveAs(str(saveFP / f"Annulus_cal_{self.source}_fit.png"))
        
        
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