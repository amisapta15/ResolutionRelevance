
import pandas as pd
from support_BoolOP import *
import multiprocessing as mp
from itertools import combinations

#path="../data_extract/"
path="../../Codes_5/data_extract/"
chunksize=20 #in minutes
rats=[20382,21012,22295,20630,22098,23783]
#rats=[24101,24116] #CA1 only Rats have single element Loc location

def worker_function(ratid):
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
    rat_data.to_json("Rat_"+str(ratid)+"_BOOLop_REC_resrel_data.json",orient="records")
    
    #Extra code to add the Neuronal Information to the DataFrame
    df2=pd.read_json(path+"Rat_"+str(ratid)+"_data_extracted.json");
    #For all Rats except 24101 and 24116, making the Loc as only CA1 or SUB
    #removing the other axis i.e. DISTAL or PROXIMAL
    df2['LOC']=df2['LOC'].apply(lambda x: x[1])
    
    
    dfP=pd.DataFrame()
    for (row_index1, row1), (row_index2, row2) in combinations(df2.iterrows(), 2):
        # Intreating over all possible unique Neuron Combinations
        for (ele_index1, ele1),(ele_index2, ele2) in zip(enumerate(row1['U_GID']),enumerate(row2['U_GID'])):
                # Intreating over all possible recordings of those Neuron pairs
                sdf=rat_data.query('U1_GID==' + str(ele1) + ' & U2_GID==' + str(ele2))
                if not sdf.empty: #if pair is found in the rat_data

                        # print("Neurons: ",row1['NeuID'],row1['N_DID'],row1['N_GID'],row1['LOC'], "--",row2['NeuID'],row2['N_DID'],row2['N_GID'],row2['LOC'])
                        # print("Units: ", ele1,ele2)
                        # print("Found")

                        #Making Sure the Locations are same coz for a recording session (because sdf comes from ratdata which is created from pairs recorded in same seesion and segregated by location ), any pair of neurons will have same location
                        assert row1['LOC']==row2['LOC'], "Locations are not same"+str(row1['LOC'])+"--"+str(row2['LOC'])
                        
                        #Adding the single recording Information of the pair to the DataFrame
                        sdf=pd.concat([sdf, 
                                        rat_data.query('U1_GID==' + str(ele1) + '& OP=="NA"'),
                                        rat_data.query('U1_GID==' + str(ele2) + '& OP=="NA"')], ignore_index=True)
                        #inserting the Neuron References to the DataFrame
                        sdf.insert(3, "N1_NeuID", np.repeat(row1['NeuID'],len(sdf)))
                        sdf.insert(4, "N1_DID", np.repeat(row1['N_DID'],len(sdf)))
                        sdf.insert(5, "N1_GID", np.repeat(row1['N_GID'],len(sdf)))
                        
                        sdf.insert(10, "N2_NeuID", np.repeat(row2['NeuID'],len(sdf)))
                        sdf.insert(11, "N2_DID", np.repeat(row2['N_DID'],len(sdf)))
                        sdf.insert(12, "N2_GID", np.repeat(row2['N_GID'],len(sdf)))
                        sdf.insert(13, "N_LOC", np.repeat(row1['LOC'],len(sdf)))
                        #print(sdf)
                        #Appending the Neuronal Information to the DataFrame
                        dfP=pd.concat([dfP,sdf], ignore_index=True)
    dfP.to_json("Rat_"+str(ratid)+"_BOOLop_Nresrel_data.json",orient="records")
    del dt
    del df2
    del dfP
    del rat_data

if __name__ == "__main__":
        pool = mp.Pool(mp.cpu_count()) 
        future_res = [pool.apply_async(worker_function, (param,)) for param in rats]
        # Close the pool and wait for the work to finish
        pool.close()
        pool.join()  