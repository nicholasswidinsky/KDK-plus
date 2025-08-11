import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path
import statistics
import math
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

detectors = { #Global dictionary that converts channel number to the detector it is associated with, and a scale factor (KeV/ADC). Note the syntax for the channel #.
    "Channel 4" : ["LS Left PMT",373.609/2104.961738528885],
    "Channel 6" : ["LS Right PMT",373.609/1668.61105108102],
    "Channel 8" : "NaI 1",
    "Channel 10" : "NaI 2",
    "Channel 12" : "NaI 3",
    "Channel 14" : "NaI 4",
}

class event:
    ############################################################################
    #   Class that stores all the details of an event.                         #
    #                                                                          #
    #   Events are defined as a groupings of n detectors that have an 'and'    #
    #   coincidence between them. The number of detectors will be defined by   #
    #   the init.txt file.                                                     #
    ############################################################################
    #   Variables:                                                             #
    #       -t (array(int)): array of length n containing the trigger times    #
    #           for each channel. Note, this list is time sorted (due to the   #
    #           way that data is saved from compass). Sorting this array       #
    #           removes the correlation between this variable and any other    #
    #           variables in this object.                                      #
    #       - E (array(float)): Array of length n containing the energy        #
    #           (integral) of the pulse corresponding to the channel           #
    #       - flags (array(str)): Array of length n containing the flags       #
    #           provided by CoMPASS which corresponds to the other variables.  #
    #       - Ch (array(str)): Array of length n containing the channels of    #
    #           the events.                                                    #
    ############################################################################
    #   Notes:                                                                 #
    #       All of the variables here are correlated, that means that the ith  #
    #       index is the same event for all variables. Sorting any of these    #
    #       variables ruins this correlation, so don't sort them when          #
    #       initializing.                                                      #
    #                                                                          #
    #       eg. for an in dex of 1 its trigger time is self.t[1], energy is    #
    #           self.E[1], flag is self.flags[1], and Channel is self.ch[1]    #
    ############################################################################
    
    
    def __init__(self,t,E,flags,Ch):
        self.t = t
        self.E = E
        self.flags = flags
        self.Ch = Ch
        
    def coincCheck(self,coincRange,coincChannels): #Currently I have this in the class but I might move it outside the class. 
        #########################################################################################
        #   Function that checks whether a group of n coincidence channels are a true coinc.    #
        #   This function checks two things.                                                    #
        #       1. Checks whether the events fall within the coincidence window.                #
        #           - This is done by taking the min and max of the list of trigger times. If   #
        #             the difference is less than or equal to the coincidence window then these #
        #             are coinsidered coincidence events.                                       #
        #       2. Checks of the coincidence channels are one of the specified coincidences in  #
        #          in the init.txt file.                                                        #
        #########################################################################################
        #   Variables:                                                                          #
        #       - coincRange (int/float): Coincidence range in ns that is read in from the      #
        #         init.txt file.                                                                #
        #       - coincChannels (array(array(str))): Array that cointains arrays of all the     #
        #         coincident channels. These channels are strings that are read in from the     #
        #         inti.txt file.                                                                #
        ######################################################################################### 
        
                
        if max(self.t)-min(self.t) <= coincRange and sorted(self.Ch) in [sorted(x) for x in coincChannels]: 
            #Checks to see if the event is coincident within the correct time range, and if the event is one of the coincident channels I have specified. 
            return True
        else:
            return False
        
    # def chParse(self,Channel):
    #     #########################################################################################
    #     #   Function that checks if a specific channel is one of the coincident channels in     #
    #     #   the coincident event.                                                               #
    #     #########################################################################################
    #     #   Variables:                                                                          #
    #     #       - Channel (int): Integer of the channel that we are checking.                   #
    #     #########################################################################################
        
    #     if Channel in sorted(self.Ch):
    #         return True
    #     else: 
    #         return False
        
    
        
class Detector():
    #####################################################################################
    #   This class will further parse the data into individual detectors. This class    #
    #   will store the trigger times, integrals, and flags for every event, while also  #
    #   storing the channel number.                                                     #
    #####################################################################################
    #   Variables:                                                                      #
    #       - trigger (array(int)): array of the trigger times of every event for this  #
    #         channel number.                                                           #
    #       - integral (array(float)): Array of the integrals for every event for this  #
    #         channel number.                                                           #
    #       - flags (array(str)): Array of the flags for every event for this channel   #
    #         number.                                                                   #
    #       - Channel (int): Integer of the channel number. This can be used with the   #
    #         globally defined dictionary to conver to the detector name and the scale  #
    #         factor to convert from ADC units to KeV.                                  #
    #####################################################################################
    
    def __init__(self,trigger,integral,flags,channel):
        self.trigger = trigger
        self.integral = integral
        self.flags = flags
        self.channel = channel
        
def coincCheckReadIn(time,channels,coincRange,coincChannels): #Currently I have this in the class but I might move it outside the class. 
    #########################################################################################
    #   Function that checks whether a group of n coincidence channels are a true coinc.    #
    #   This function checks two things.                                                    #
    #       1. Checks whether the events fall within the coincidence window.                #
    #           - This is done by taking the min and max of the list of trigger times. If   #
    #             the difference is less than or equal to the coincidence window then these #
    #             are coinsidered coincidence events.                                       #
    #       2. Checks of the coincidence channels are one of the specified coincidences in  #
    #          in the init.txt file.                                                        #
    #########################################################################################
    #   Variables:                                                                          #
    #       - time (array(int)): Array of integers representing the trigger time for each   #
    #         channel.                                                                      #
    #       - channels (array(str)): Array of strings containing the channels for the       #
    #         associated time.                                                              #
    #       - coincRange (int/float): Coincidence range in ns that is read in from the      #
    #         init.txt file.                                                                #
    #       - coincChannels (array(array(str))): Array that cointains arrays of all the     #
    #         coincident channels. These channels are strings that are read in from the     #
    #         inti.txt file.                                                                #
    ######################################################################################### 
    
    if max(time)-min(time) <= coincRange and sorted(channels) in [sorted(x) for x in coincChannels]: 
        #Checks to see if the event is coincident within the correct time range, and if the event is one of the coincident channels I have specified. 
        return True
    else:
        return False
        
        
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
        
        dataFilePath = lines[0].split("\n")[0].split('\t')[1] #Grab the data filepath from the first line 
        coincWindow = int(lines[1].split("\n")[0].split('\t')[1]) #Grab the coincidence window from the second line
        CoincChannelsList = lines[2].split("\n")[0].split('\t')[1:] # grab the list of all the coinc channels
        
        CoincChannels = [] #Create a blank array for the final formatting of the coinc channel list. 
        for ch in CoincChannelsList: #Loop over all the different coinc channel pairs. Then separate them based into their own sub array (I might change this later but this works right now so I am keeping it. )
            CoincChannels.append(sorted(map(int,ch.split(",")))) #Append a sorted coinc channel lists to a general list. 


    return dataFilePath,coincWindow,CoincChannels

        
def readInData(file,CoincChannels,coincRange,scale = False,):
    
    totalChannels = [] #List of all the channels read in. This will be returned at the end as a sorted list of all the unique channels that are used. 
    tempT,tempCh,tempE,tempF = [],[],[],[] #Create temporary arrays to store data for an event. These events will then be used to initialize an event object. 
    events = [] #List of all the event objects that will be initialized during this process. 
    with open(file) as f: #Start reading in the file. 
        next(f) #Skip the header. 
        i = 1 #Start an index that will count up to the number of channels that are in coincidence. 
        j = 0
        for line in f: #Read each line one at a time. This will help with speed and memory.
            lines = line.split('\n')[0].split(';') #Split each line along semi-colons
            totalChannels.append(int(lines[1])) #Append the channel number to the toal list. 
            tempCh.append(int(lines[1])) #Append the channel number to a temporary list to check if this group of events is a proper coinicdent event. 
            tempT.append(int(lines[2])/1e3) #Append the trigger time to a temporary list. 
            tempE.append(int(lines[3])) #Append te integral to a temporary list
            tempF.append(lines[5]) #Append the flags to a temporary list. 
            
            if i == len(CoincChannels): #Check to see if the index matches the number of channels that are in coincidence. 
                if coincCheckReadIn(time = tempT,channels=tempCh,coincRange=coincRange,coincChannels=CoincChannels): #check to see if the coincidence is a proper coincidence (this will return true most of the time).
                    # if scale: #Do this later, make a global dictionary for the scale factors for each detector so I can just reference that. 
                        
                    events.append(event(t = tempT,E = tempE,flags = tempF,Ch = tempCh)) #initialize a new event object using the temporary arrays. 
                
                    tempT,tempCh,tempE,tempF = [],[],[],[] #Reset the temp arrays to being empty
                    i = 1 #Reset the index so it looks for 'n' new events       
                else:
                    for k in range(len(CoincChannels)-1): #If there was a write error for the event then we remove all other data for that event, only keeping the last data point. 
                        tempT.pop(0) #Removes the 0th element from the list. 
                        tempE.pop(0)
                        tempCh.pop(0)
                        tempF.pop(0)
                    i = 2 #Start the index with one item already in the list. 
            else:
                i +=1
                
    Channels = sorted(list(set(totalChannels))) #Sort, and remove dupliates of the total list of channels. 
    return events, Channels


def parseCh(Events,Channels):
    #####################################################################
    #   Function that parses all of the data in the event class and     #
    #   initializes the detectors class.                                #
    #####################################################################
    #   Variables:                                                      #
    #       - Events (array(event)): Array of all the event class       #
    #         objects.                                                  #
    #       - Channels (array(int)): Array of all the channels used     #
    #         for this experiment.                                      #
    #####################################################################
    
    ParsedData = [[] for _ in range(len(Channels))] #Create a large 2D array that will house all coincident events corresponding to each channel. 
    
    for i,event in enumerate(Events):
        a = 2
        
            
    
initFilePath = '/home/nick/PhD/KDK+/code/CAEN_analysis_init.txt'

dataFilePath,coincWindow,CoincChannels = readInInit(initFilePath)
events, channels = readInData(file = dataFilePath,CoincChannels=CoincChannels,coincRange=coincWindow)
print(len(events))


testEvent = event([10000,10100,10500],[10,20,30],['0x400','0x400','0x400'],[0,4,2])
# print(testEvent.coincCheck(500,[['Ch0','Ch2','Ch3'],['Ch4','Ch0','Ch2'],['Ch0','Ch2','Ch6'],['Ch0','Ch2','Ch7'],['Ch0','Ch2','Ch8']]))
print(testEvent.coincCheck(500,CoincChannels))
print(testEvent.Ch[0])




# print("Coincidence code ran successfully")
# print("""\

#                                            #######                                                  
#                                ####+     ##.......+#           ######                               
#                              #......+#  #...........# #...## #.......##                             
#                             #.........# #...........-#......#.........##                            
#                             #..........##............#......#..........#  ###...+#                  
#                    ##....#####..........#............#......#..........###........#                 
#                   #.........+#..........#######.....------######+.....-#...........#                
#                   #...........##.....##...............-------------##.#...........##                
#                   ##............#.#+...................---------------##..........#####             
#          ##+....+####...........#......................------------------#.#...##......+#           
#         #...........-#........#.........................-------------------#+##.........#           
#         #..............#....#...........................---------------------#+.........#           
#         ##..............###.............................-----------------------#-.....+#            
#          ##.............#...............................-------------------------#..####.....##     
#           ##.........##...............--###.............-----------###.-----------+#...........#    
#   ###.......##.......#...............--######...........----------#####.------------#..........#    
#  #.............##..#.................-+  ####..........----------+- ###.-------------#........#     
# #................##..................-+######..........----------+#####.--------------#.....##      
# #................#...................-+######.........-----------+#####.---------------#.#########  
#  #..............+-...................--######.........------------#####.----------------#.........##
#   #.............#....................--######........-------------#####.----------------#..........#
#     ##.........+-....................--######.......--------------####.-----------------#.........##
#      ######+...#......................--####.......---------------------------##--------#.......-#  
#   ##...........#..................................---------------------------##---------#....###    
#  #.............+-...............................----------------------------#----------########     
# #.........++....#.................#............---------------------------##----------#........##   
# #...............-+............####+..........-------------------------+##------------++..........#  
#  #...............#...............-#######..----------------------####.--------------#............#  
#    ##...........###-................------####################---------------------#...........-#   
#         #### ##.....#...............---------------------------------------------#-...#########     
#           ##.........#+---.....------------------------------------------------#+........+#         
#          #.......+.....##----------------------------------------------------#-............#        
#         #................-#-----------------------------------------------##.#..............#       
#         #..............##...##-----------------------------------------##......##..........##       
#         #............##.........##---------------------------------##-.#.........##########         
#         ##........###..........##....###---------------------###+.......##........#                 
#            ######   #..........#.....#.....#.-########+-.......#..........#.......#                 
#                     #.........#............#.........##........##.........###-..##                  
#                     #.......###...........##..........#........# #.........#                        
#                       ######  #...........##..........##......#   ##......#                         
#                                #.........# #.........#.-+#####       #####                          
#                                 #-.....##++#+.......##---+     +++.--------++                       
#                                    ##+...---+#.....#..---+  ++.---------+---+++++++++               
#                                  +---+..-------.+.-++...-+  .-----------+++.............++          
#                                  --+ +..-------------+.---+.+--.--------..................++        
#                                  +-++...--------------...--.--------.+...-+.................+       
#                                      ++...---------..--.-.-+--.----+..+......................+      
#                                      +......---------.+-.---+-.---+.-.........................+     
#                                   +++...............++-.+.-.-----++...........................+     
#                               ++.........................+.-.---+-............................-+    
#                             +-.................-.............++-+..............................+    
#                            +....................................+....................+++++....++    
#                            ++.....................................................++       +..+    .
#                            +.....................................+...............+          .-+     
#                           +......................................+.............++           -+      
#                           +......................................+............+        +++++        
#                           +.....................................+ +....+....+                       
#                           +......+++   +++.....................+    ++.++++                         
#                            +...+            +++............+++                                      
#                             +..+                 ++++++++++                                         
#                               ++++                                                                  
# """)

# print("""

#                                                                                                 .   
#                                                                                                 .   
#                                      :-.                                                        .   
#                                   :@@@@@@@#                                                     .   
#                                 :@@@@@@@@@@:                                                    .   
#                                  %@@@@@@@@@+                                                    .   
#                                  .@@@@@@@@@@                                                    .   
#                                   .@@@@@@@@@-                                                   .   
#                                    *@@@@@@@@@                                                   .   
#                                     @@@@@@@@@:                                                  .   
#                                     :@@@@@@@@#                                                  .   
#                                      @@@@@@@@@%*++=.                                            .   
#                                     .%@@@@@%*+======+-.                                         .   
#                                    =**#%%*+=------====+=.   .                                   .   
#                                   =++++#*==---::----======.                                     .   
#                                  :++==+*=----:::::::--=--===.                                   .   
#                                 .+=====+=---::::::::::------=+.                                 .   
#                                 -+======---:::::::::::::-----==+.                               .   
#                                .++======----:::::::::::::----===+#.                             .   
#                                -+=======----:::::::::::::::----=+#%.                            .   
#                               .++=======-----::::::::::::::--====*%+                            .   
#                               :++============-----:::::::::---====*%:                           .   
#                              .+++===++**********++=----::::::--===+#%.                          .   
#                              .*++********************++=---::--===++%-                          .   
#                              -**************************++=-----====*%                          .   
#                              +*****************************+==--====+%:                         .   
#                             .********************************+=======#+                         .   
#                             :*****++**++********+++**++******++=====++#.                        .   
#                             **+**++++*+++**++++++=+*++++******++=====+#.                        .   
#                            .***++++*+++++++++++++=+++++++++*****+====+*:                        .   
#                            :****@@@@@%+=+++==+++==+%@@@@#*++*****====++=                        .   
#                            :#**%@@@@@@%++++==++=+*@@@@*@@@@******+====+*                        .   
#                            .*++*@@@@@@%++=====+=+#@@@@@@@@*=++***+===++*.                       .   
#                             ***=*%@@#*##+++===+*##*#@@@@%*==+****+====+*.                       .   
#                             :##**++=++***+++++****++=====+*******=====+*:                       .   
#                              +##***+++=+*++++**++=++++++********+=====+*-                       .   
#                              .##**+++==+#***#%*+==++++++********======++=                       .   
#                               .#**++++==*%%%#*+=++++++++*******=--====++=                       .   
#                                +##**++***##%#*****++++++*****+=--=====+++                       .   
#                                .**##**##%%%%%%%##**+*******+=----=====+++                       .   
#                                 =**##%%################**+=-------====++=                       .   
#                                 :**+*##%%%######%#####*+=---------====++=                       .   
#                                 .*+++***###%%%%%##**++=----------=====+*-                       .   
#                                 .*+++++****##***+++===-----------=====+*:                       .   
#                                 .*+++++++++*+++====-------------====+++*.                       .   
#                                 .*++====+=++====---------------=====+++#.                       .   
#                                 :*+=======++===---------------======++**                        .   
#                                 :++=======++=----------------==+===+++*+.                       .   
#                                .++========+=----------------==++===+*+=+=.                      .   
#                              .+**+========+=---------------=+*#+===+*+====+-.                   .   
#                            .++***+========+=--------------=+*##+==+**++====+++.                 .   
#                          .++++*##+========+=---------::--==*##+===+*#**+++===+++.               .   
#                         :++++*#%#+========+--------------+*#*+===++*%@##*++=+=++*-.             .   
#                        -++**+*#.=*+======+=----------:-=+*#*====+++*#@. .++=++++**=             .   
#                       .++#%#**. :#+++===+*+==--::-----+***=====+++**#:  .++**+++***             .   
#                       -**@%  .  .%##*+=+##*+--------+*##*===+++*++*#=.   +*+**++***             .   
#                         .       .#@@#*+*%%*=-----=++**+====+++++**#+         -*+**=             .   
#                                  -%@@##@@%#+---=+***+===++++****##=            .                .   
#                                   *%@@@@@@%*++*#*+====+++*****#%%.                              .   
#                                   .#%@@@@@@@%%*=====+******##%##:                               .   
#                                    %#@@@@@@@#*++++++*****#%%#**#=                               .   
#                                   .##%@@@@@%#********##%%%#*****+                               .   
#                                   .###%%@@%####*#*####=.##******#                               .   
#                                   .####%%%%%.   ..  .   +#*******.                              .   
#                                  .+***###%%%-           :***+*****.                             .   
#                                  =******###%%           -**+++++**#.                            .   
#                                 .*+**+***###%.         .***++++*****                            .   
#                                  =*********-.          :***+++++*+**.                           .   
#                                                          ..:-:::..                              .   
#                                                                                                 .   
#                                                                                                 .   
#                     -@@       =@@=  -@@@@@@@@*        @@=       @@@                             .   
#                     -@@     .@@%    -@@     +@@@=     @@+     #@@:                              .   
#                     -@@    @@@.     -@@        @@%    @@+   .@@*        #@+                     .   
#                     -@@  =@@-       -@@         @@-   @@+  @@@          #@+                     .   
#                     -@@.@@#         -@@         *@#   @@+*@@        ----%@*----                 .   
#                     -@@@@@@#        -@@         *@#   @@@@#@@.     .@@@@@@@@@@@.                .   
#                     -@@-  *@@       -@@         @@=   @@@   @@#         *@+                     .   
#                     -@@    :@@+     -@@        #@@    @@+    #@@        *@+                     .   
#                     -@@      @@@    -@@     .*@@#     @@+     .@@-                              .   
#                     -@@       -@@.  -@@@@@@@@@-       @@+       @@@                             .   
#                                                                                                 .   
#                                                                                                 .   
#                                                                                                 .   
# """)

# print("""                                                                                        
#                                #%#.                                                                 
#                                %@#%%-                                                               
#                                 %#@#%%=   
#                                 #%%*%*#%*                                                           
#                                  :#%#*@*%%#                                                         
#                                  .%%*#%*#*#%#.                                                      
#                                    **##+%#%*+%%.                                                    
#                                    +%#=*#*%**#+%%:                                                  
#                                     =##%=+%**%#%*%%-                                                
#                                     .%%*+##*%#*%*%#@@-                                              
#                                      #@+*%*##+##*%#%#%%=                                            
#                                       *%##*##*%**%*%#@#@%+                                          
#                                        %%*+%**%##%#%#%#@%@@+                                        
#                                         +#%@##%+*##%%%%%%@#---                                      
#                                          :%%**%#+%#%%#@%@@*:..==                                    
#                                            =##%#*###@#%@@@*:...:++                                  
#                                             :@%*#%%#%@%@%%##....:-==                                
#                                              %@%##%##%%#@@@#*+..-:--#.                              
#                                              %%@%#*%@#@@%@*#*...::--===                             
#                                              -%##%##@@%%*-:....:....-++                             
#                                              :@%###%@@%::-:..:..=:..--=                             
#                                              :%%##%@@@+.........:.-=.:+                             
#                                                %##%%%*.....-:..::....+%@@@@@@@%@%.                  
#     .+%%#=------+#%@@%#*#=...:-*%%#%%%%#+-:..   :*+:... .   ...  ...:#@@@@@@#*@@@%%#                
#                           :%##:..  :*###%@@@@@-::....               -%%@#%%@@%@#@@@@@*.             
#                               .#=  .=**##%%#-...  ..                :#%%#@%%**=*#*#%.:              
#                               =@%= ..*@%@%-....   ..                 *###%*+***#%-                  
#                                 @@%%%@@@+:.:.      ...             ...*#*@#*#%@-                    
#                                 :@#=:=#--...-: : .....  ....  ... ....:*#@@@@@#*+                   
#                                  %#*==:=:.:-::......:.:.........:.::.::=+#@**++--+*.                
#                                    :#-+**@@@@@@@@@:    :%%@@@@@@@@@@@@@@@%%#+=-:..:=:               
#                                    ::#* .#@#%@@@%       .@@%%%%##@@%@%@@@@@@@#-::.:.=:              
#                                     -@*==.:%#%@+          .@%##%##%*#%#%%@@#@@%+:-.::=.             
#                                     *# #-===##.              .#%%*#*#%#%%%%%%@@##+.-:-+             
#                                       -%#. *%.                 @@%%#%*%#@%#%#@%@%%##+*+-            
#                                       @...-#                    *@%@%%*%%*#%#@@%@@@@%%@%.           
#                                        + -+                       -%-%%%%*%##%%#@@@@@@@%%           
#                                       +=-=                          +%%%**#*%#%*@@#@@@@%#-          
#                                      .#=.                             +@*%%*%#%*%##@%%@#%#          
#                                      %%                                *#%%*%*%#@####@#%#%#         
#                                    .%+                                 :%@##@#%*%*@*%##%*%#:        
#                                   :@-                                   . #%@*%*#*%*%*%*@%%#        
#                                  :%:                                      +#@*%*###*%#%*@#*#+       
#                                 +#                                         =*#%**%##%#%%****#.      
#                                %#                                           .%@*#%*%*##%#@*+#+      
#                               %-                                             .*##@+@#%*@##***#:     
#                             :%.                                                .%%*#*%*%**@++##     
#                            +=                                                    :@%#####+%#*+%-    
#                           -                                                        #@%%*@##@###%    
#                                                                                     .%@*%%*#@**#-   
#                                                                                       +%%@###%**#   
#                                                                                         *%###%##%   
#                                                                                           :%#%:%%.  
#                                                                                            .@@=:%:  
#                                                                                                 =   
# """)