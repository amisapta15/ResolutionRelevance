import sys
sys.path.append('../lib')
from func_MSR import *
from func_Info import *
import pandas as pd

path="../data_extract/"
chunksize=20 #in minutes
rats=[20382,24101,21012,22295,20630,22098,23783,24116]

for ratid in rats:
    df=pd.read_json(path+"Rat_"+str(ratid)+"_data_extracted.json");
    datas=[]
    print("for Rat",ratid,"\n")
    for nu in range(df.shape[0]):
        du=df.iloc[nu].duration
        for uni in range(len(du)):
                if du[uni]>(60*chunksize): ## if duration is more than 20 minutes 
                        print("Neuron ",nu," Unit ",uni," with NID ",df.iloc[nu].NeuID," Duration ",du[uni]/60," minutes")
                        
                        if du[uni]>(60*chunksize*2): ## if duration is more than 40 minutes
                                print("!! Neuron ",nu," Unit ",uni," with NID ",df.iloc[nu].NeuID," Duration ",du[uni]/60," minutes")
                                
                                
                        #poissonification of the entire spike train
                        p_st=poissonification(df.iloc[nu].u_spiketime[uni],df.iloc[nu].time_range[uni][0],df.iloc[nu].time_range[uni][-1],0.001) #wih 1ms time-bins
                        
                        #MSR
                        nspikes,msr,OHKs,OHSs,maxHkpHs,Tohk,MHKs,Tmhks=MSR(df.iloc[nu].u_spiketime[uni],df.iloc[nu].time_range[uni][0],df.iloc[nu].time_range[uni][0]+(chunksize*60),200) #with slices of chunk-size minutes
                        _,ps_msr,ps_OHKs,ps_OHSs,ps_maxHkpHs,ps_Tohk,ps_MHKs,ps_Tmhks=MSR(p_st,df.iloc[nu].time_range[uni][0],df.iloc[nu].time_range[uni][0]+(chunksize*60),200)  #with slices of chunk-size minutes
                        
        
                        #Information quantities          
                        pinfos_content,pinfos_rate=SPinfo(df.iloc[nu].u_spiketime[uni],df.iloc[nu].X[uni],df.iloc[nu].Y[uni],df.iloc[nu].t[uni],df.iloc[nu].time_range[uni][0],df.iloc[nu].time_range[uni][0]+(chunksize*60),0.01,30) #with slices of chunk-size minutes and 30x30 spatial bins for 150x150 cm box
                        ps_pinfos_rate,ps_pinfos_content=SPinfo(p_st,df.iloc[nu].X[uni],df.iloc[nu].Y[uni],df.iloc[nu].t[uni],df.iloc[nu].time_range[uni][0],df.iloc[nu].time_range[uni][0]+(chunksize*60),0.01,30) #with slices of chunk-size minutes and 30x30 spatial bins for 150x150 cm box
                        
                        
                        hinfos_content,hinfos_rate=HDinfo(df.iloc[nu].u_spiketime[uni],df.iloc[nu].HD[uni],df.iloc[nu].t[uni],df.iloc[nu].time_range[uni][0],df.iloc[nu].time_range[uni][0]+(chunksize*60),0.01,50)
                        #with slices of chunk-size minutes and ~7deg polar bins 
                        ps_hinfos_rate,ps_hinfos_content=HDinfo(p_st,df.iloc[nu].HD[uni],df.iloc[nu].t[uni],df.iloc[nu].time_range[uni][0],df.iloc[nu].time_range[uni][0]+(chunksize*60),0.01,50)
                        
                        
                        
                        dicts={'RAT_ID':df.iloc[nu].RAT_ID,
                                'NeuID':df.iloc[nu].NeuID, 
                                'N_GID': df.iloc[nu].N_GID,
                                'LOC':df.iloc[nu].LOC, 
                                'task':df.iloc[nu].task[uni],
                                'U_GID':df.iloc[nu].U_GID[uni], 
                                'time_range':df.iloc[nu].time_range[uni], 
                                'duration':df.iloc[nu].duration[uni], 
                                'REC_Date':df.iloc[nu].REC_Date[uni],
                                'extracted_pinfo':np.nan if df.iloc[nu].STAT == 'NA' else df.iloc[nu].STAT['pos']['info'],
                                'extracted_hinfo':np.nan if df.iloc[nu].STAT == 'NA' else df.iloc[nu].STAT['hd']['info'],
                                'Nspikes':nspikes,
                                'MSR':msr,
                                'MHK':MHKs,
                                'dt_MHK':Tmhks,
                                'OHK':OHKs,
                                'OHS':OHSs,
                                'dt_OHK':Tohk,
                                'max_HSHK':maxHkpHs,
                                'PInfo':{'rate':pinfos_rate,
                                        'content':pinfos_content
                                        },  
                                'HInfo':{'rate':hinfos_rate,
                                        'content':hinfos_content
                                        },  
                                'poss_MSR':ps_msr,
                                'poss_MHK':ps_MHKs,
                                'poss_dt_MHK':ps_Tmhks,
                                'poss_OHK':ps_OHKs,
                                'poss_OHS':ps_OHSs,
                                'poss_dt_OHK':ps_Tohk,
                                'poss_max_HSHK':ps_maxHkpHs,
                                'poss_PInfo':{'rate':ps_pinfos_rate,
                                        'content':ps_pinfos_content
                                        } , 
                                'poss_HInfo':{'rate':ps_hinfos_rate,
                                        'content':ps_hinfos_content
                                        },               
                                };
                        datas.append(dicts)
    rat_data=pd.DataFrame(datas)
    rat_data.to_json("Rat_"+str(ratid)+"_resrel_data.json",orient="records")
    del df
    del rat_data
