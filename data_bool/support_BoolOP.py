import sys
sys.path.append('../lib')
from func_MSR import *
import numpy as np
from itertools import combinations

def empty_formatted(dt,nRec,u1,u2,op,nspike,LOC):

    rows={'RAT_ID': dt.iloc[nRec].RAT_ID,
        'Rec_GID' : dt.iloc[nRec].REC_GID, 
        'Rtask'   : dt.iloc[nRec].task,
        
        'U_LOC'   :LOC,
        'U1_GID':  u1,
        'OP': op,
        'U2_GID':u2,
 
        'Nspikes':nspike,
        'MSR':np.nan,
        'MHK':np.nan,
        'dt_MHK':np.nan,
        'OHK':np.nan,
        'OHS':np.nan,
        'dt_OHK':np.nan,
        'max_HSHK':np.nan
        };
    return rows;

def single_formatted(dt,nRec,chunksize,spike_train,st,u,LOC):
    nspikes,msr,OHKs,OHSs,maxHkpHs,Tohk,MHKs,MHSs,Tmhks=MSR(spike_train,st,st+(chunksize*60),200)
    
    rows={'RAT_ID': dt.iloc[nRec].RAT_ID,
        'Rec_GID' : dt.iloc[nRec].REC_GID, 
        'Rtask'   : dt.iloc[nRec].task,
        
        'U_LOC'   :LOC,
        'U1_GID':  u,
        'OP': 'NA',
        'U2_GID':u,
 
        'Nspikes':nspikes,
        'MSR':msr,
        'MHK':MHKs,
        'MHS':MHSs,
        'dt_MHK':Tmhks,
        'OHK':OHKs,
        'OHS':OHSs,
        'dt_OHK':Tohk,
        'max_HSHK':maxHkpHs
        };
    return rows;
def muli_formatted(dt,nRec,chunksize,spike_train,st,u1,u2,op,LOC):
    #for The spike train of length = chunksize
    nspikes,msr,OHKs,OHSs,maxHkpHs,Tohk,MHKs,MHSs,Tmhks=MSR(spike_train,st,st+(chunksize*60),200)
    
    rows={'RAT_ID': dt.iloc[nRec].RAT_ID,
        'Rec_GID' : dt.iloc[nRec].REC_GID, 
        'Rtask'   : dt.iloc[nRec].task,
        
        'U_LOC'   :LOC,
        'U1_GID':  u1,
        'OP': op,
        'U2_GID':u2,
 
        'Nspikes':nspikes,
        'MSR':msr,
        'MHK':MHKs,
        'MHS':MHSs,
        'dt_MHK':Tmhks,
        'OHK':OHKs,
        'OHS':OHSs,
        'dt_OHK':Tohk,
        'max_HSHK':maxHkpHs
        };
    return rows;


def adding_rows(datas,dt,nRec,loci,chunksize,location,mul_unit):
   bin_time=0.001 #1 milisec in seconds time resolution of binning for logical operations
   st,en=dt.iloc[nRec].REC_timerange
   edges = np.arange(st, en + bin_time, bin_time)
   # Binarization of spike train lambda function
   bined_spikes = lambda unit_id: np.array([False if x == 0 else True for x in np.histogram(dt.iloc[nRec].U_spiketimes[dt.iloc[nRec].U_GIDs.index(unit_id)], bins=edges)[0]])

   #For Single Units
   for items in loci[location]:
      sptr=np.array(dt.iloc[nRec].U_spiketimes[dt.iloc[nRec].U_GIDs.index(items)]); #spiketrain
      if len(np.where(sptr <= st+(chunksize*60))[0]) > 1:
         datas.append(single_formatted(dt,nRec,chunksize,sptr,st,items,location))
      else:
         datas.append(empty_formatted(dt,nRec,items,items,'NA',len(np.where(sptr <= st+(chunksize*60))[0]),location))
   if mul_unit:   
      #For combined Units
      for items in list(combinations(loci[location], 2)):
         uid1=items[0];uid2=items[1]; #unit id
         #print(uid1," ",uid2)
         
         #binarizing the spike train with 1 ms time resolution for logical operation
         bst1=bined_spikes(uid1);bst2=bined_spikes(uid2);
         # & is the element wise (bit wise) AND ; where retuns indices of the True after bool operation;
         #finally select from the time array from st to en in 1 ms resolution to create spike timestamps
         if sum(bst1)!=sum(bst2):
            # AND operation
            and_spike_train = np.arange(st,en+bin_time,bin_time)[np.where(bst1 & bst2)]
            if len(np.where(and_spike_train <= st+(chunksize*60))[0]) > 1:
               datas.append(muli_formatted(dt,nRec,chunksize,and_spike_train,st,uid1,uid2,'AND',location))
            else:
               datas.append(empty_formatted(dt,nRec,uid1,uid2,'AND',len(np.where(and_spike_train <= st+(chunksize*60))[0]),location))
            
            # OR operation
            or_spike_train = np.arange(st,en+bin_time,bin_time)[np.where(bst1 | bst2)]
            if len(np.where(or_spike_train <= st+(chunksize*60))[0]) > 1:
               datas.append(muli_formatted(dt,nRec,chunksize,or_spike_train,st,uid1,uid2,'OR',location))
            else:
               datas.append(empty_formatted(dt,nRec,uid1,uid2,'OR',len(np.where(or_spike_train <= st+(chunksize*60))[0]),location))
            
            # XOR operation
            xor_spike_train = np.arange(st,en+bin_time,bin_time)[np.where(bst1 ^ bst2)]
            if len(np.where(xor_spike_train <= st+(chunksize*60))[0]) > 1:
               datas.append(muli_formatted(dt,nRec,chunksize,xor_spike_train,st,uid1,uid2,'XOR',location))
            else:
               datas.append(empty_formatted(dt,nRec,uid1,uid2,'XOR',len(np.where(xor_spike_train <= st+(chunksize*60))[0]),location))    
         else:
            print("Duplicate! nREC=",nRec," loc: ",location," units: ",items)
     
            
   return datas 


  