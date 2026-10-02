import ROOT
import numpy as np
from pathlib import Path

class CoicidenceData:
    #################################################################
    #   This class is a relic of my previous data reduction code.   #
    #   This class takes coincidence data then save it to a file    #
    #   while still reading in the data. This save on memory usage, #
    #   but makes it very slow code.                                #
    #################################################################
    
    def __init__(self,ch,t,data):
        self.ch = [ch]
        self.t = [t]
        self.data = [data]
        
    def addEvent(self,ch,t,data):
        self.ch.append(ch)
        self.t.append(t)
        self.data.append(data)
        
    def writeToDisk(self,filepath,filename):
        
        sortedCh = sorted(self.ch)
        self.t = [x for _,x in sorted(zip(self.ch,self.t),key=lambda pair:pair[0])]
        self.data = [x for _,x in sorted(zip(self.ch,self.data),key=lambda pair:pair[0])]
        self.ch = sortedCh
        
        coincChannelsStr = '_'.join([str(i) for i in sortedCh])
        
        coincfile = f"{filepath}/{filename}_coinc_{coincChannelsStr}.txt" 
        
        if coincfile in openedfiles:
            f = open(coincfile, 'a')
        else:
            coincList.append(sortedCh)
            openedfiles.append(coincfile)
            f = open(coincfile,"w+")
            f.write(f"Data reduction code output.\nEach line saves the information for every channel in the coincidence. Note that currently waveforms are not saved in this data reduction. \nEach channel has the following headers: BOARD;CHANNEL;TIMETAG;ENERGY;ENERGYSHORT;FLAGS;PROBE_CODE;SAMPLES\nNumber of channels in this coincidence: {len(self.ch)}\n")
        
        f.write(';'.join(self.data) + '\n')
        f.close()
        
def AddEvents(trees,buffers,chKeys,coincIndex,data):
    #Create a dictionary that determines the ROOT variable for the different types of ints/float/long/...
    rootTypeMap = {
        np.dtype("float64"): "D",
        np.dtype("float32"): "F",
        np.dtype("int64"):   "L",
        np.dtype("int32"):   "I",
        np.dtype("uint64"):  "l",
        np.dtype("uint32"):  "i",
        np.dtype("int16"):   "S",
        np.dtype("uint16"):  "s",
        np.dtype("int8"):    "B",
        np.dtype("uint8"):   "b",
    }
    
    perChannelCols = ["Timestamp", "Board", "Energy", "EnergyShort", "Flags"]
    
    if chKeys not in trees:
        name = "coinc_ch_" + "_".join(str(c) for c in chKeys)
        tree = ROOT.TTree(f"{name}", f"Coincidence channels {chKeys}")

        buf = {}
        for ch in chKeys:
            for col in perChannelCols:
                dt = data[col].dtype
                if dt not in rootTypeMap:
                    dt = np.dtype("float64")   # safe fallback for unsupported types
                arr = np.zeros(1, dtype=dt)
                buf[(col, ch)] = arr
                tree.Branch(f"{col}_{ch}", arr, f"{col}_{ch}/{rootTypeMap[dt]}")

        trees[chKeys] = tree
        buffers[chKeys] = buf

    tree = trees[chKeys]
    buf = buffers[chKeys]

    # --- fill: match each event's channel to its buffer slot ---
    eventChannels = data["Channel"][coincIndex]
    for idx, ch in zip(coincIndex, eventChannels):
        ch = int(ch)
        for col in perChannelCols:
            dt = buf[(col, ch)].dtype
            val = data[col][idx]
            # buf[(col, ch)][0] = data[col][idx].astype(dt) if dt != np.dtype("O") else data[col][idx]
            buf[(col,ch)][0] = dt.type(val)

    tree.Fill()
    
    # if chKeys in trees:
    #     pass
    # else:
    #     chStr = "_".join(str(c) for c in chKeys)
    #     name = "coinc_ch_" + "_".join(str(c) for c in chKeys)

    #     tree = ROOT.TTree(f"{name}_TTree",f"Coincidence Channels {chKeys} TTree")
        
    #     channelArr = data["Channel"][coincIndex]
    #     TimeArr =  data["Timestamp"][coincIndex]
    #     boardArr =  data["Board"][coincIndex]
    #     EnergyArr = data["Energy"][coincIndex]
    #     EnergyShortArr =  data["EnergyShort"][coincIndex]
    #     FlagsArr = data['Flags'][coincIndex]
    
    
    #     for i,ch in enumerate(chKeys):
    #         tree.Branch(f"Channel_{ch}",channelArr[i],f"Channel_{ch}/{rootTypeMap[channelArr.dtype]}")
        
        

def ReadInData(file,OutputFile,coincWindow):
    filename = Path(file).stem
    f = ROOT.TFile.Open(str(file))
    
    tree = f.Get("Data_R")
    tree.Print()
    tree.Show(0)
    
    cols = ["Channel","Timestamp","Board","Energy","EnergyShort","Flags"]
    
    df = ROOT.RDataFrame(tree) #Red in the TTree data into a RDataFrame for easier analysis. 
    data = df.AsNumpy(columns = cols) #Read in the data from the date frame into Numpy objects for easier analysis. 
    
    n = len(data["Channel"]) #Get the length of the data file that is being read in. 
    
    trees = {} #Create an empty dictionary that will store all of the different TTrees based on the coincidence channels
    buffers = {}
    
    outfile = ROOT.TFile.Open(f"{OutputFile}/{filename}_sorted_coinc.root","RECREATE")
    outfile.cd()
    
    coincCounter = 0
    i = 0
    while i < n: #loop over every event in the file. 
        startTime = data['Timestamp'][i] / 1e3
        j = i + 1
        while j < n:
            dt = (data['Timestamp'][j] / 1e3) - startTime #Calculate the time difference between this event and the first coinc event.

            if dt < coincWindow:
                j +=1
            else:
                break
            
        coincCounter +=1
        
        if coincCounter % 10000 == 0:
            print(f"Sorted {coincCounter} coincidence events")
        
        coincIndex = np.arange(i,j) #Make an array of the inices of all the events in the coincidence.
        channelKey = tuple(sorted(int(c) for c in data["Channel"][coincIndex])) #Make a tuple of all the coinc channels.

        AddEvents(trees,buffers,channelKey,coincIndex,data)
        i=j
    
        
    
    
    
    
    
    for t in trees.values():
        t.Write()
    outfile.Close()
    f.Close()

    
    
    # f.ls()

    # for key in f.GetListOfKeys():
    #     print(key.GetName(),key.GetClassName())
        
    # tree = f.Get("Data_R")
    # print("tree.Print command:")
    # tree.Print()
    # tree.Show(0)

    # df = ROOT.RDataFrame(tree)

    # h = df.Filter("Channel == 0").Histo1D(("H_energy", "Energy;Energy (ch);Counts", 100, 0, 4000), "Energy")

    # c = ROOT.TCanvas() 
    # c.Divide(2,1)
    # c.cd(1)
    # h.Draw()
    # c.Update()


    # EnergyHist = ROOT.TH1D('Energy_hist',"Energy Hist", 100,0,4000)

    # for i in range(tree.GetEntries()):
        
    #     nb = tree.GetEntry(i)
    #     if nb <= 0:
    #         continue
    #     flag = tree.Flags
    #     print(flag)
    #     ch = tree.Channel
    #     Energy = tree.Energy
    #     if ch == 0:
    #         EnergyHist.Fill(Energy)
    # c.cd(2)
    # EnergyHist.Draw()
    # c.Draw()



    # f.Print()
    

file = Path("/home/nick/PhD/KDK+/Pilot_exp/2026_09_11/2026_09_11_Pilot_experiment_cal_data/RAW/SDataR_2026_09_11_Pilot_experiment_cal_data.root")

OutputFileDir = file.parent / 'coinc_sorted'
OutputFileDir.mkdir(parents = True, exist_ok=True)

coincWindow = 500 #ns

ReadInData(file,OutputFileDir,coincWindow)