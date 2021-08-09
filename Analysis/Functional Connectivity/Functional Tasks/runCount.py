import os
import pdb
i = 0
while i < 2500:
    file = r"/suphys/aloe8475/Documents/CODE/Analysis/Functional\ Connectivity/Functional\ Tasks/count.py"
    file2= r"/import/silo2/aloe8475/Documents/CODE/Data/Functional\ Connectivity/Leda\ Graphs/BetaSweep/BetaSweepLedaGraph" + str(i+1) +".gw"
    #pdb.set_trace()
    os.system(file+' ' + file2)    
    i += 1
