import numpy as np

def ReadInFileNaI(file, numChannel, channels):
    
    data = []
    for c in range(numChannel):
        data.append([ [] for _ in range(3)])
    # channels = [2*n for n in range(numChannel)]
    
    
    
    with open(file) as f:
        next(f)
        for line in f:
            lines = line.split('\n')[0].split(';')
            
            for c in range(numChannel):
                if int(lines[1]) == channels[c]:
                    data[c][0].append(float(lines[2]))
                    data[c][1].append(float(lines[3]))
                    data[c][2].append(lines[5])
    return data

def sumChannels(file):
    
    i = 0
    sumInt  = 0 
    sumTotals = []
    with open(file) as f:
        next(f)
        for line in f:
            lines = line.split('\n')[0].split(';')
            
            sumInt+=int(lines[1])
            
            i+=1
            
            if i == 4:
                sumTotals.append(sumInt)
                i = 0
                sumInt = 0
            
    return sumTotals

def DataCleanUp(file,channels,coincT):
    
    data = [[],[],[],[],[],[]] #Create empty arrays for the final data set.
    tempCh, tempT,tempE,tempF = [],[],[],[] #Create temporary arrays that will store the data of the 4 events while we verify that the 4 events have all been recorded properly. 
    board,energyS = [],[]
    
    i = 0 #Set a counter that gets reset once 4 events have been read. 
    counter= 0
    with open(file) as f:
        next(f)
        for line in f: 
            lines = line.split('\n')[0].split(';')
            
            board.append(int(lines[0]))
            tempCh.append(int(lines[1]))
            tempT.append(int(lines[2]))
            tempE.append(int(lines[3]))
            energyS.append(int(lines[4]))
            tempF.append(lines[5])
            
            
            
            i +=1 
            if i == len(channels): 
                
                sumCh = sum(tempCh)
                timeDiff = tempT[-1]-tempT[0]
                
                # print(sumCh)
                # print(sum(channels))
                # print(timeDiff)
                # print(coincT)
                
                if sumCh == sum(channels):# and timeDiff <= coincT:
                    # print('went to second if statement')
                    for j in range(len(channels)):
                        data[0].append(board[j])
                        data[1].append(tempCh[j])
                        data[2].append(tempT[j])
                        data[3].append(tempE[j])
                        data[4].append(energyS[j])
                        data[5].append(tempF[j])
                        
                    tempCh, tempT,tempE,tempF = [],[],[],[] #If there is no write error, then the temporary files get returned to empty. 
                    board,energyS = [],[]
                    
                    i =0 
                else:   
                    for j in range(len(channels)-1):
                        tempCh.pop(0)
                        tempT.pop(0)
                        tempE.pop(0)
                        tempF.pop(0)
                        board.pop(0) 
                        energyS.pop(0)  
  
                    i =1         
                
                
                
                
                
            counter +=1
            
    return data



filepath = '/home/nick/PhD/KDK+/Large_LSC_testing/Position_tests_outer_sleeve/2025_07_17/Large_LSC_vessel_Cs137_triple_coinc_Vertical_scatter_v4_central_position/RAW/'
filename = "SDataR_Large_LSC_vessel_Cs137_triple_coinc_Vertical_scatter_v4_central_position.CSV"

Channels = [0,2,6]


fileData = ReadInFileNaI(f'{filepath}/{filename}', len(Channels), Channels)


print("Number of events in each channel before clean-up:")
for i in range(len(Channels)):
    print(f"Ch {Channels[i]}: {len(fileData[i][0])}")

# print(f"Ch {Channels[0]}: {len(fileData[0][0])}")
# print(f"Ch {Channels[1]}: {len(fileData[1][0])}")
# print(f"Ch {Channels[2]}: {len(fileData[2][0])}")
# print(f"Ch {Channels[3]}: {len(fileData[3][0])}")



data = DataCleanUp(file = f"{filepath}/{filename}",channels = Channels,coincT= 700e3)



correctedFileName = filename.split('.')[0]

f = open(f'{filepath}/{correctedFileName}_corrected.CSV','w+')
f.write("BOARD;CHANNEL;TIMETAG;ENERGY;ENERGYSHORT;FLAGS\n")
for i in range(len(data[0])):
    f.write(f"{data[0][i]};{data[1][i]};{data[2][i]};{data[3][i]};{data[4][i]};{data[5][i]}\n")
    
f.close()

fileData = ReadInFileNaI(f'{filepath}/{correctedFileName}_corrected.CSV', len(Channels), Channels)

print("Number of events in each channel after clean-up:")
for i in range(len(Channels)):
    print(f"Ch {Channels[i]}: {len(fileData[i][0])}")
# print(f"Ch {Channels[0]}: {len(fileData[0][0])}")
# print(f"Ch {Channels[1]}: {len(fileData[1][0])}")
# print(f"Ch {Channels[2]}: {len(fileData[2][0])}")
# print(f"Ch {Channels[3]}: {len(fileData[3][0])}")

                        