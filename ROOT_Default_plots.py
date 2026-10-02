import ROOT
import numpy as np
from pathlib import Path
import xml.etree.ElementTree as ET
from itertools import combinations
import math
import matplotlib.pyplot as plt


#########################################################################################
# Global Settings
#########################################################################################
ROOT.gROOT.SetBatch(True)
ROOT.EnableImplicitMT() #Enables implicit multithreading to speed up the code. 

ROOT.gInterpreter.Declare("""
    std::string FlagsToHex(unsigned int flags) {
        char buf[16];
        snprintf(buf, sizeof(buf), "0x%04X", flags);
        return std::string(buf);
    }
    """)

ROOT.gStyle.SetLabelSize(0.05, "XYZ")   # axis tick-label size (numbers along the axis)
ROOT.gStyle.SetTitleSize(0.05, "XYZ")   # axis title size (e.g. "Energy (ADC)")
ROOT.gStyle.SetTitleFontSize(0.06)      # the histogram's own title (top of plot)
ROOT.gStyle.SetTextSize(0.05)           # general text (legends, TLatex, etc., unless overridden)
ROOT.gStyle.SetLegendTextSize(0.1)   
#########################################################################################


def SortData(df,cPairs,cols,tree):
    #This function takes a ROOT Data frame to creat the appropriate channels that are needed (flags -> hex, time difference etc)
    
    #Calculate the time difference for each ch pair
    tDiffCols = []
    for chi in cPairs:
        for chj in cPairs:
            colName = f"Time_diff_ch_{chi}_ch_{chj}"
            # df = df.Define(colName, f"(Timestamp_{chi} - Timestamp_{chj}) / 1000")
            df = df.Define(colName, f"((Long64_t)Timestamp_{chi} - (Long64_t)Timestamp_{chj})/1000")
            
            tDiffCols.append(colName)
    
    # tDiffCols = []
    # for chi, chj in combinations(cPairs, 2):   # only unique pairs, no self-pairs
    #     colName = f"Time_diff_ch_{chi}_ch_{chj}"
    #     df = df.Define(colName, f"Timestamp_{chi} - Timestamp_{chj}")
    #     tDiffCols.append(colName)
            
    #Convert the flasgs from decimal to hex
    
    dtCols = []
    flagsHexCols = []
    for ch in cPairs:
        colName = f"Flags_{ch}_hex"
        dtColName = f"Timestamp_dt_{ch}"
        
        startT = tree.GetLeaf(f"Timestamp_{ch}").GetValue(0)
        df = df.Define(colName, f"FlagsToHex(Flags_{ch})")
        df = df.Define(dtColName, f"((Long64_t)Timestamp_{ch} - (Long64_t){startT})/1e12")
        
        flagsHexCols.append(colName)
        dtCols.append(dtColName)
        
    cols += dtCols
    cols += tDiffCols
    cols += flagsHexCols
    return df, cols

def createFlagHist(data,ch):
    flagDict = {}
    
    # data = df.AsNumpy(columns = [f"Flags_{ch}_hex"])
    FlagList,flagCounts = np.unique(data[f'Flags_{ch}_hex'],return_counts = True)
    
    FlagList = [f.decode("utf-8") for f in FlagList] 
    
    return FlagList,flagCounts

def createFlagBarChart(data,ch,flagList,flagCounts,chLabels):
    
    nBins = len(flagCounts)
    hBar = ROOT.TH1D(f"Flags_bar_char_{ch}",f"{chLabels[f'Channel {ch}'][0]} Flag Counts;Flag;Counts",nBins,0,nBins)
    
    for i, (val,count) in enumerate(zip(flagList,flagCounts)):
        hBar.SetBinContent(i+1,int(count))
        hBar.GetXaxis().SetBinLabel(i+1,f"{val}")
        
    hBar.SetFillColor(ROOT.kAzure - 4)
    hBar.SetBarWidth(0.8)
    hBar.SetBarOffset(0.1)
    hBar.SetStats(0)
    hBar.GetXaxis().SetLabelSize(0.08)
    hBar.GetXaxis().LabelsOption("v")
    hBar.GetXaxis().SetTitleOffset(4.0)
    
    return hBar
    
def CreateFlagHists(df,ch,FlagList,chLabels):
    hists = {}
    
    for flag in FlagList:
        h = df.Filter(f'Flags_{ch}_hex == "{flag}"').Histo1D((f"Energy_ch_{ch}_flag_{flag}",f"{chLabels[f'Channel {ch}'][0]} Energy, Flag {flag};Energy (ADC);Counts",100, 0, 4200),f"Energy_{ch}")    
    
        hists[flag] = h
    return hists
    
def canvasPadding(pad):
    pad.SetLeftMargin(0.15)
    pad.SetRightMargin(0.20)   # extra room for COLZ palette
    pad.SetTopMargin(0.10)
    pad.SetBottomMargin(0.15)

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

file = Path("/home/nick/PhD/KDK+/Pilot_exp/2026_09_11/2026_09_11_Pilot_experiment_cal_data/RAW/coinc_sorted/SDataR_2026_09_11_Pilot_experiment_cal_data_sorted_coinc.root")
BaseSaveFilePath = file.parent / 'Default_plots'
settingsFilePath = file.parent.parent.parent / 'settings.xml'

chLabels = ReadInChannelNames(settingsFilePath) #Create a directory for the channel number to channel name 

fileName = file.stem #Grab the name of the file. 

coincPairs = [[0,1,2],
              [0,1,3],
              [0,1,4],
              [0,1,5]]

# coincPairs = [[0,1,2],
#               [0,1,3],
#               [0,1,4],
#               [0,1,5],
#               [0,1,2,3],
#               [0,1,2,4],
#               [0,1,2,5],
#               [0,1,3,4],
#               [0,1,3,5],
#               [0,1,4,5]]

f = ROOT.TFile.Open(str(file))
f.ls() #Print the tree contents of the file.

Hist1DCanvas = []
Hist2DCanvaslist = []
TDiff1DCanvas = []
TDiff2DCanvaslist = []
stabilityCanvas = []
# flagBar = []
# flagDist = []
# stackedFlagDist = []

#Create a distinct coloour pallet for use in the flag Histograms.
distinctColors = [
    ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta,
    ROOT.kOrange+1, ROOT.kCyan+1, ROOT.kViolet, ROOT.kSpring+4,
    ROOT.kPink+9, ROOT.kAzure+7, ROOT.kYellow+2, ROOT.kTeal+3,
    ROOT.kGray+2, ROOT.kBlack
]

for cPairs in coincPairs:
    
    ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)
    
    cols = []
    for ch in cPairs:
        ChCol = [f"Timestamp_{ch}",f"Board_{ch}",f"Energy_{ch}",f"EnergyShort_{ch}",f"Flags_{ch}"]
        # cols = np.concatenate((cols,LSCChCol))
        cols += ChCol
            
    
    cPairsStr = '_'.join([str(i) for i in cPairs])
    
    saveFilePath = BaseSaveFilePath / f"Coinc_{cPairsStr}_figures"
    saveFilePath.mkdir(parents = True, exist_ok = True)
    
    tree = f.Get(f"coinc_ch_{cPairsStr}")
    tree.Print() #Print an overview of the tree structure. 
    
    
    df = ROOT.RDataFrame(tree)
    df,cols = SortData(df,cPairs,cols,tree)
    
    
    Hist1D = []
    Hist2D = [[] for i in cPairs]
    TDiffHist = [[] for i in cPairs]
    TDiffHist2D = [[] for i in cPairs]
    TstabHist = []
    FlagLegend = []
    flagBar = []
    flagDist = []
    stackedFlagDist = []
    
    # Book all Max() requests together (still lazy - nothing runs yet)
    maxDtVals = {ch: df.Max(f'Timestamp_dt_{ch}') for ch in cPairs}

    # This first .GetValue() access triggers ONE combined pass covering every channel's Max at once
    maxDtResolved = {ch: maxDtVals[ch].GetValue() for ch in cPairs}
    
    
    for ch in cPairs:
        Hist1D.append(df.Histo1D((f"Energy_hist_ch_{ch}", f"{chLabels[f'Channel {ch}'][0]};{chLabels[f'Channel {ch}'][0]} Integral (ADC);Counts",100,0,4000), f"Energy_{ch}"))
        # Hist1D[-1].SetStats(0)
        
        # TstabHist.append(df.Histo2D((f"Time_stability_{chLabels[f'Channel {ch}'][0]}",f"{chLabels[f'Channel {ch}'][0]} integral vs Timestamp;TimeStamp (s);{chLabels[f'Channel {ch}'][0]} integral (ADC)",100,0,df.Max(f'Timestamp_dt_{ch}').GetValue(),100,0,4000),f"Timestamp_dt_{ch}",f"Energy_{ch}"))
        TstabHist.append(df.Histo2D((f"Time_stability_{chLabels[f'Channel {ch}'][0]}",f"{chLabels[f'Channel {ch}'][0]} integral vs Timestamp;TimeStamp (s);{chLabels[f'Channel {ch}'][0]} integral (ADC)",100,0,maxDtResolved[ch],100,0,4000),f"Timestamp_dt_{ch}",f"Energy_{ch}"))
        TstabHist[-1].SetStats(0)
        
    # for chi, chj in combinations(cPairs, 2):
    for i,chi in enumerate(cPairs):
        for chj in cPairs:
            Hist2D[i].append(df.Histo2D((f"Energy_hist_ch_{chi}_{chj}", f"{chLabels[f'Channel {chi}'][0]} integral vs {chLabels[f'Channel {chj}'][0]} integral;{chLabels[f'Channel {chi}'][0]} integral (ADC);{chLabels[f'Channel {chj}'][0]} integral (ADC)",100,0,4000,100,0,4000), f"Energy_{chi}",f"Energy_{chj}"))
            Hist2D[i][-1].SetStats(0)
            
            TDiffHist2D[i].append(df.Histo2D((f"TDiff_hist_ch_{chi}_{chj}_vs_Energy_{chi}", f"Tdiff {chLabels[f'Channel {chi}'][0]} - {chLabels[f'Channel {chj}'][0]} vs {chLabels[f'Channel {chi}'][0]} integral;Time Difference (ns);{chLabels[f'Channel {chi}'][0]} (ADC)",100,0,4000,100,-600,600),f"Energy_{chi}", f"Time_diff_ch_{chi}_ch_{chj}"))
            TDiffHist2D[i][-1].SetStats(0)
            
            TDiffHist[i].append(df.Histo1D((f"Time_Difference_{chLabels[f'Channel {chi}'][0]}-{chLabels[f'Channel {chj}'][0]}", f"Time Difference {chLabels[f'Channel {chi}'][0]} - {chLabels[f'Channel {chj}'][0]};Time Difference (ns);counts",100,-600,600),f"Time_diff_ch_{chi}_ch_{chj}") )
            TDiffHist[i][-1].SetStats(0)
        
        
       
        
    data = df.AsNumpy(columns = cols)
    
    # fig, ax = plt.subplots(1,1)
    # ax.hist(data['Time_diff_ch_0_ch_2'],bins = 100)
    
    # plt.show()
    
    flagOnly = df.AsNumpy(columns=[f"Flags_{ch}_hex" for ch in cPairs])

    flagListByCh = {}
    flagDist_ch = {}
     
    
    for ch in cPairs:
        # FlagList,flagCounts = createFlagHist(data,ch) 
        # flagBar.append(createFlagBarChart(data,ch,FlagList,flagCounts,chLabels))
        # flagDist.append(CreateFlagHists(df,ch,FlagList,chLabels))
        
        FlagList, flagCounts = createFlagHist(flagOnly, ch)
        flagBar.append(createFlagBarChart(flagOnly, ch, FlagList, flagCounts, chLabels))
        flagDist.append(CreateFlagHists(df, ch, FlagList, chLabels))
        
        
    H1DCanvas = ROOT.TCanvas(f"Energy_hist_ch_{cPairsStr}_canvas", f"Energy Hist ch {cPairsStr}",2400,1200)
    StabCanvas = ROOT.TCanvas(f"Stability_Plots_{cPairsStr}_canvas",f"Stability Plots ch {cPairsStr}",2400,1200)
    
    ncol = math.ceil(np.sqrt(len(cPairs)))
    nrow = math.ceil(len(cPairs)/ncol)
    
    H1DCanvas.Divide(ncol,nrow)
    StabCanvas.Divide(len(cPairs),4)
    
    for i,hist in enumerate(Hist1D):
        pad = H1DCanvas.cd(i+1)
        canvasPadding(pad)

        Hist1D[i].Draw()
        H1DCanvas.Update()
        
        pad = StabCanvas.cd(i+1)
        canvasPadding(pad)
        Hist1D[i].Draw()
        StabCanvas.Update()
        
        pad = StabCanvas.cd(len(cPairs)+i+1)
        canvasPadding(pad)
        TstabHist[i].Draw("colz")
        StabCanvas.Update()
        
        pad = StabCanvas.cd(2*len(cPairs) + i + 1)
        pad.SetLogy(1)
        pad.SetLeftMargin(0.15)
        pad.SetRightMargin(0.20)   # extra room for COLZ palette
        pad.SetTopMargin(0.10)
        pad.SetBottomMargin(0.35)
        # canvasPadding(pad)
        flagBar[i].Draw()
        StabCanvas.Update()
        
        
        ROOT.gStyle.SetPalette(ROOT.kRainBow)
        pad = StabCanvas.cd(3*len(cPairs)+i+1)
        pad.SetLogy(1)
        pad.SetLeftMargin(0.15)
        pad.SetRightMargin(0.40)   # extra room for COLZ palette
        pad.SetTopMargin(0.10)
        pad.SetBottomMargin(0.15)
        # canvasPadding(pad)
        
        sortedFlags = sorted(flagDist[i].items(), key=lambda item: item[1].GetEntries(), reverse = True)
        
        hs = ROOT.THStack(f"hs_{ch}",f"{chLabels[f'Channel {cPairs[i]}'][0]} Stacked Flag Histograms")
        for key,value in sortedFlags:
            # flagDist[i][f'{key}'].Draw("SAME")
            hs.Add(value.GetPtr())
  
        # hs.Draw("NOSTACK PLC")
        stackedFlagDist.append(hs)
        StabCanvas.Update()
        
        # Build the legend in the reserved right-margin space
        legend = ROOT.TLegend(0.62, 0.15, 0.99, 0.85)   # NDC coords: x1,y1,x2,y2 - sits in the 0.20 right margin
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)     # transparent background
        legend.SetTextSize(0.08)
        # legend.SetMargin(0.35)
        legend.SetNColumns(2)

        for idx,((key, value), hist) in enumerate(zip(sortedFlags, hs.GetHists())):
            colour = distinctColors[idx % len(distinctColors)]
            hist.SetLineColor(colour)
            hist.SetLineWidth(1)
            legend.AddEntry(hist, key, "l")   # "l" = line-style swatch, matches PLC line-colored hists

        hs.Draw("NOSTACK")
        legend.Draw()
        FlagLegend.append(legend)
        stackedFlagDist.append(hs)
        StabCanvas.Update()
        
    
    savefileName = saveFilePath / "Integral_dist_1D.pdf"
    H1DCanvas.SaveAs(str(savefileName))
    savefileName = saveFilePath / "Integral_dist_1D.png"
    H1DCanvas.SaveAs(str(savefileName))
    Hist1DCanvas.append(H1DCanvas)
    
    savefileName = saveFilePath / "Stability_plots.pdf"
    StabCanvas.SaveAs(str(savefileName))
    savefileName = saveFilePath / "Stability_plots.png"
    StabCanvas.SaveAs(str(savefileName))
    
    
    ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)
    
    TdiffCanvas = ROOT.TCanvas(f"TDiff_ch_{cPairsStr}_canvas", f"Tdiff ch {cPairsStr}",2400,1200)
    TdiffCanvas.Divide(len(cPairs),len(cPairs))
    
    H2DCanvas = ROOT.TCanvas(f"Energy_hist_2D_{cPairsStr}_canvas", f"Energy Hist ch {cPairsStr}", 2400,1200)
    H2DCanvas.Divide(len(cPairs),len(cPairs))
    
    TdiffECanvas = ROOT.TCanvas(f"TDiff_vs_Energy_ch_{cPairs}_canvas",f"Tdiff vs Integral ch {cPairsStr}",2400,1200)
    TdiffECanvas.Divide(len(cPairs),len(cPairs))
    
    padInd = 1
    # plotInd = 0
    for i,chi in enumerate(cPairs):
        for j,chj in enumerate(cPairs):
            if j < i:
                padInd +=1
            else: 
                pad = TdiffCanvas.cd(padInd)
                pad.SetLogy(1)
                canvasPadding(pad)
                
                TDiffHist[i][j].Draw()
                TdiffCanvas.Update()
                
                
                pad = H2DCanvas.cd(padInd)
                canvasPadding(pad)
                
                Hist2D[i][j].Draw("colz")
                H2DCanvas.Update()
                
                
                pad = TdiffECanvas.cd(padInd)
                canvasPadding(pad)
                
                TDiffHist2D[i][j].Draw("colz")
                TdiffECanvas.Update()
                
                
                
                
                padInd +=1
                
    TDiff1DCanvas.append(TdiffCanvas)
    Hist2DCanvaslist.append(H2DCanvas)
    TDiff2DCanvaslist.append(TdiffECanvas)
    
    TdiffFileName = saveFilePath / "Tdiff_dist_1D.pdf"
    TdiffCanvas.SaveAs(str(TdiffFileName))
    TdiffFileName = saveFilePath / "Tdiff_dist_1D.png"
    TdiffCanvas.SaveAs(str(TdiffFileName))
    
    
    H2DFileName = saveFilePath / "Hist_2D.pdf"
    H2DCanvas.SaveAs(str(H2DFileName))
    H2DFileName = saveFilePath / "Hist_2D.png"
    H2DCanvas.SaveAs(str(H2DFileName))
    
    
    TdiffEFileName = saveFilePath / "Tdiff_dist_2D.pdf"
    TdiffECanvas.SaveAs(str(TdiffEFileName))
    TdiffEFileName = saveFilePath / "Tdiff_dist_2D.png"
    TdiffECanvas.SaveAs(str(TdiffEFileName))
    
    
    