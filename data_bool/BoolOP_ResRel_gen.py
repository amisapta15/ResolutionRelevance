
import pandas as pd

from support_BoolOP import *

path="../data_extract/"
chunksize=20 #in minutes
rats=[20382,24101,21012,22295,20630,22098,23783,24116]
#rats=[22098]
for ratid in rats:
    dt=pd.read_json(path+"Rat_"+str(ratid)+"_RECdata_extracted.json");
    allD=[]
    print("for Rat",ratid,"\n")
    for nRec in range(dt.shape[0]):       
        #checking if the duration is better than chunk size
        if dt.iloc[nRec].REC_duration > (chunksize*60): #in seconds
                
                #segrecating CA1 and SUB units in the recording
                loci = {'CA1': [], 'SUB': [], 'others': []}
                for (s, v) in zip(dt.iloc[nRec].U_GIDs, dt.iloc[nRec].N_LOCs):
                        if 'CA1' in v:
                                loci['CA1'].append(s)
                        elif 'SUB' in v:
                                loci['SUB'].append(s)
                        else:
                                loci['others'].append(s)
                
                print(nRec,"->",loci)
                no_CA1=len(loci['CA1'])
                no_SUB=len(loci['SUB'])        
                
                #SUB neurons loop
                if no_SUB!= 0:
                        if no_SUB == 1:
                                allD=adding_rows(allD,dt,nRec,loci,chunksize,'SUB',False);
                        else:
                                allD=adding_rows(allD,dt,nRec,loci,chunksize,'SUB',True);
                
                #CA1 Loop
                if no_CA1 != 0:
                        if no_CA1 == 1:
                                allD=adding_rows(allD,dt,nRec,loci,chunksize,'CA1',False);
                        else:
                                allD=adding_rows(allD,dt,nRec,loci,chunksize,'CA1',True);
    rat_data=pd.DataFrame(allD)
    rat_data.to_json("Rat_"+str(ratid)+"_BOOLop_resrel_data.json",orient="records")
    del dt
    del rat_data
    
  
