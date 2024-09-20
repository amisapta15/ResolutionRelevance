
import pandas as pd
from support_BoolOP import *
import multiprocessing as mp
from itertools import combinations,product

#path="../data_extract/"
path="../../Codes_5/data_extract/"
chunksize=20 #in minutes
rats=[20382,24101,21012,22295,20630,22098,23783,24116]
#rats=[24101,24116] #CA1 only Rats have single element Loc location



def worker_function(ratid,save_reclist=True): 
    dt=pd.read_json(path+"Rat_"+str(ratid)+"_RECdata_extracted.json");
    allD=[]
    #print("for Rat",ratid,"\n")
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
                
                #print(nRec,"->",loci)
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
    if save_reclist:
           rat_data.to_json("Rat_"+str(ratid)+"_BOOLop_REC_resrel_data.json",orient="records")

    del dt


    #Extra code to add the Neuronal Information to the DataFrame
    df2=pd.read_json(path+"Rat_"+str(ratid)+"_data_extracted.json");
    #For all Rats except 24101 and 24116, making the Loc as only CA1 or SUB
    #removing the other loc info i.e. DISTAL or PROXIMAL
    df2['LOC']=df2['LOC'].apply(lambda x: x[1] if len(x)>1 else x[0])
    
    
    dfP=pd.DataFrame()
    for (row_index1, row1), (row_index2, row2) in combinations(df2.iterrows(), 2):
        # Intreating over all possible unique Neuron Combinations
        if (row1['LOC']==row2['LOC']):
                #Making Sure the Locations are same coz for a recording session (because sdf comes from ratdata which is created from pairs recorded in same seesion and segregated by location ), any pair of neurons will have same location
                # assert row1['LOC']==row2['LOC'], "Locations are not same"+row1['LOC']+"--"+row2['LOC']
                for ele1,ele2 in set(tuple(sorted(t)) for t in product(row1['U_GID'], row2['U_GID']) if t[0] != t[1]):
                        if (row1['duration'][row1['U_GID'].index(ele1)]>=60*chunksize) and (row2['duration'][row2['U_GID'].index(ele2)]>=60*chunksize):  
        
                                # Intreating over all possible recordings of those Neuron pairs
                                sdf=rat_data.query('U1_GID==' + str(ele1) + ' & U2_GID==' + str(ele2))
                                if not sdf.empty: #if pair is found in the rat_data
                                        
                                        # print("Neurons: ",row1['NeuID'],row1['N_DID'],row1['N_GID'],row1['LOC'], "--",row2['NeuID'],row2['N_DID'],row2['N_GID'],row2['LOC'])
                                        # print("Units: ", ele1,ele2)
                                        # print("Found")
                                
                                        #Adding the single recording Information of the pair to the DataFrame
                                        sdf=pd.concat([sdf, 
                                                        rat_data.query('U1_GID==' + str(ele1) + '& OP=="NA"'),
                                                        rat_data.query('U1_GID==' + str(ele2) + '& OP=="NA"')], ignore_index=True)
                                        
                                        #inserting the Neuron References to the DataFrame
                                        sdf.insert(3, "N1_NeuID", row1['NeuID'])
                                        sdf.insert(4, "N1_DID", row1['N_DID'])
                                        sdf.insert(5, "N1_GID", row1['N_GID'])
                                        
                                        sdf.insert(10, "N2_NeuID", row2['NeuID'])
                                        sdf.insert(11, "N2_DID", row2['N_DID'])
                                        sdf.insert(12, "N2_GID", row2['N_GID'])
                                        
                                        sdf.insert(13, "N_LOC", row1['LOC'])
                                        #print(sdf.iloc[:,3:17])
                                        
                                        #Appending the Neuronal Information to the DataFrame
                                        dfP=pd.concat([dfP,sdf], ignore_index=True)
                                
    dfP.to_json("Rat_"+str(ratid)+"_BOOLop_Nresrel_data.json",orient="records")
    del dt
    del df2
    del dfP

                
if __name__ == "__main__":
        pool = mp.Pool(10) 
        future_res = [pool.apply_async(worker_function, (param,)) for param in rats]
        # Close the pool and wait for the work to finish
        pool.close()
        pool.join()  




# def Rec_List_Worker(ratid,save_reclist=True):
#     dt=pd.read_json(path+"Rat_"+str(ratid)+"_RECdata_extracted.json");
#     allD=[]
#     #print("for Rat",ratid,"\n")
#     for nRec in range(dt.shape[0]):       
#         #checking if the duration is better than chunk size
#         if dt.iloc[nRec].REC_duration > (chunksize*60): #in seconds
                
#                 #segrecating CA1 and SUB units in the recording
#                 loci = {'CA1': [], 'SUB': [], 'others': []}
#                 for (s, v) in zip(dt.iloc[nRec].U_GIDs, dt.iloc[nRec].N_LOCs):
#                         if 'CA1' in v:
#                                 loci['CA1'].append(s)
#                         elif 'SUB' in v:
#                                 loci['SUB'].append(s)
#                         else:
#                                 loci['others'].append(s)
                
#                 #print(nRec,"->",loci)
#                 no_CA1=len(loci['CA1'])
#                 no_SUB=len(loci['SUB'])        
                
#                 #SUB neurons loop
#                 if no_SUB!= 0:
#                         if no_SUB == 1:
#                                 allD=adding_rows(allD,dt,nRec,loci,chunksize,'SUB',False);
#                         else:
#                                 allD=adding_rows(allD,dt,nRec,loci,chunksize,'SUB',True);
                
#                 #CA1 Loop
#                 if no_CA1 != 0:
#                         if no_CA1 == 1:
#                                 allD=adding_rows(allD,dt,nRec,loci,chunksize,'CA1',False);
#                         else:
#                                 allD=adding_rows(allD,dt,nRec,loci,chunksize,'CA1',True);
#     rat_data=pd.DataFrame(allD)
#     if save_reclist:
#            rat_data.to_json("Rat_"+str(ratid)+"_BOOLop_REC_resrel_data.json",orient="records")
#     return rat_data
#     del dt
#     del rat_data

# def Neuron_ID_Worker(ratid,dt):  

#     #Extra code to add the Neuronal Information to the DataFrame
#     df2=pd.read_json(path+"Rat_"+str(ratid)+"_data_extracted.json");
#     #For all Rats except 24101 and 24116, making the Loc as only CA1 or SUB
#     #removing the other loc info i.e. DISTAL or PROXIMAL
#     df2['LOC']=df2['LOC'].apply(lambda x: x[1] if len(x)>1 else x[0])
    
    
#     dfP=pd.DataFrame()
#     for (row_index1, row1), (row_index2, row2) in combinations(df2.iterrows(), 2):
#         # Intreating over all possible unique Neuron Combinations
#         if (row1['LOC']==row2['LOC']):
#                 #Making Sure the Locations are same coz for a recording session (because sdf comes from ratdata which is created from pairs recorded in same seesion and segregated by location ), any pair of neurons will have same location
#                 # assert row1['LOC']==row2['LOC'], "Locations are not same"+row1['LOC']+"--"+row2['LOC']
#                 for ele1,ele2 in set(tuple(sorted(t)) for t in product(row1['U_GID'], row2['U_GID']) if t[0] != t[1]):
#                         if (row1['duration'][row1['U_GID'].index(ele1)]>=60*chunksize) and (row2['duration'][row2['U_GID'].index(ele2)]>=60*chunksize):  
        
#                                 # Intreating over all possible recordings of those Neuron pairs
#                                 sdf=dt.query('U1_GID==' + str(ele1) + ' & U2_GID==' + str(ele2))
#                                 if not sdf.empty: #if pair is found in the rat_data
                                        
#                                         # print("Neurons: ",row1['NeuID'],row1['N_DID'],row1['N_GID'],row1['LOC'], "--",row2['NeuID'],row2['N_DID'],row2['N_GID'],row2['LOC'])
#                                         # print("Units: ", ele1,ele2)
#                                         # print("Found")
                                
#                                         #Adding the single recording Information of the pair to the DataFrame
#                                         sdf=pd.concat([sdf, 
#                                                         dt.query('U1_GID==' + str(ele1) + '& OP=="NA"'),
#                                                         dt.query('U1_GID==' + str(ele2) + '& OP=="NA"')], ignore_index=True)
                                        
#                                         #inserting the Neuron References to the DataFrame
#                                         sdf.insert(3, "N1_NeuID", row1['NeuID'])
#                                         sdf.insert(4, "N1_DID", row1['N_DID'])
#                                         sdf.insert(5, "N1_GID", row1['N_GID'])
                                        
#                                         sdf.insert(10, "N2_NeuID", row2['NeuID'])
#                                         sdf.insert(11, "N2_DID", row2['N_DID'])
#                                         sdf.insert(12, "N2_GID", row2['N_GID'])
                                        
#                                         sdf.insert(13, "N_LOC", row1['LOC'])
#                                         #print(sdf.iloc[:,3:17])
                                        
#                                         #Appending the Neuronal Information to the DataFrame
#                                         dfP=pd.concat([dfP,sdf], ignore_index=True)
                                
#     dfP.to_json("Rat_"+str(ratid)+"_BOOLop_Nresrel_data.json",orient="records")
#     del dt
#     del df2
#     del dfP

# if __name__ == "__main__":
#         for rat in rats:
#                 print("Rat:",rat)
#                 allrat=Rec_List_Worker(rat,save_reclist=False)
#                 Neuron_ID_Worker(rat,allrat)
#                 del allrat
                
                
# # if __name__ == "__main__":
# #         pool = mp.Pool(10) 
# #         future_res = [pool.apply_async(worker_function, (param,)) for param in rats]
# #         # Close the pool and wait for the work to finish
# #         pool.close()
# #         pool.join()  
