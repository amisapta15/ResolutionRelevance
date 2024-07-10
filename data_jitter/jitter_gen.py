import sys
sys.path.append('../lib')
from func_MSR import *
import pandas as pd
import numpy as np
from numpy.random import Generator, PCG64 ,MT19937, PCG64DXSM
from scipy.stats import truncnorm
import multiprocessing as mp

path="../data_extract/"
ratid=21012
chunksize=20 #in minutes
df=pd.read_json(path+"Rat_"+str(ratid)+"_data_extracted.json");

## Custom Imported from MATLAB Database
#D.animals(3).sessions(18) contains- id: 104
#D.animals(3).recordings(18).sessions.userData.containedNeurons.ids
neuron_Dids=[99,105,124,130,132,137,143,148,154,161,165,169,171,176,179,182,186,202,204,207,209,211,213,214,215,219,225,231,239,243,247,251,253,256,259,262]

unit_Gids=[4893,4942,5070,5122,5164,5262,5317,5374,5430,5486,5525,5562,5574,5650,5691,5735,5780,5941,5954,5967,5980,5990,6002,6010,6018,6041,6082,6122,6193,6240,6275,6295,6305,6327,6339,6353]

# Function to find the indices
def find_indices(df, column, items):
    results = []
    for item in items:
        for row_idx, row in df.iterrows():
            if item in row[column]:
                results.append((row_idx, row[column].index(item),item))
    return results
def truncated_normal_toss(mu, jitter, M, no_of_samples, seed=1244):
    # Initialize the random number generator
    #rng = Generator(PCG64DXSM(seed=seed))
    rng = Generator(PCG64DXSM())
    # Calculate the lower and upper bounds for the truncation
    sigma=jitter/(2*np.pi)
    lower, upper = (-jitter - mu) / sigma, (jitter - mu) / sigma
    truncated_normal = truncnorm(lower, upper, loc=mu, scale=sigma)
    # Generate the samples
    samples = truncated_normal.rvs(size=(M, no_of_samples), random_state=rng)
    return samples
def worker_function(idx):
    jtdat=df.iloc[nu].u_spiketime[uni] + samples[:,idx]
    return MSR(jtdat,df.iloc[nu].time_range[uni][0],df.iloc[nu].time_range[uni][0]+(chunksize*60),200)


if __name__ == "__main__":
    NuUnis = find_indices(df[df['N_DID'].isin(neuron_Dids)], 'U_GID', unit_Gids)
    no_of_samples = 10000

    for tt in [5,10,50,100,500,1000]: #5 10 50 100 500 1000 milisecs
        jitter=tt/1000 #5ms    
        datas=[]
        for (nu,uni,ugids) in NuUnis:
            assert df.iloc[nu].U_GID[uni]==ugids, "UIDs don't match"
            print(f"Calculating for neuron: %d, unit: %d, U_GID: %d" % (nu,uni,ugids))
            nj_results=np.array(MSR(df.iloc[nu].u_spiketime[uni],df.iloc[nu].time_range[uni][0],df.iloc[nu].time_range[uni][0]+(chunksize*60),200))
            # Output 9 dims
            #no-of-spikes->len(datas),msr,OHKs,OHSs,maxHkpHs,Tohks,MHKs,MHSs,Tmhks
            #samples from the truncated normal distribution for the jittering
            samples = truncated_normal_toss(0, jitter, int(nj_results[0]), no_of_samples)
            #Create pool of workers for jittering MSR calculation
            pool = mp.Pool(14) 
            future_res = [pool.apply_async(worker_function, (param,)) for param in range(no_of_samples)]
            res = [f.get() for f in future_res]
            # Close the pool and wait for the work to finish
            pool.close()
            pool.join()

            nspikes,msr,OHKs,OHSs,maxHkpHs,Tohk,MHKs,MHSs,Tmhks=nj_results
            #Calculating the displacement of the jittered MSR and other ResRel Quantities against unjittered MSR and other ResRel Quantities
            displacement=np.abs(res - nj_results) 
            _,mu_msrDSP,mu_OHKsDSP,mu_OHSsDSP,mu_maxHkpHsDSP,mu_TohkDSP,mu_MHKsDSP,mu_MHSsDSP,mu_TmhksDSP=np.mean(displacement,axis=0)

            row = {
                "ratid": ratid,
                "nu_index": nu,
                "unit_index": uni,
                "U_GID": ugids,
                "nspikes": nspikes,
                "msr": msr,
                "OHKs": OHKs,
                "OHSs": OHSs,
                "maxHkpHs": maxHkpHs,
                "Tohk": Tohk,
                "MHKs": MHKs,
                "MHSs": MHSs,
                "Tmhks": Tmhks,
                "mean_msrDSP": mu_msrDSP,
                "mean_OHKsDSP": mu_OHKsDSP,
                "mean_OHSsDSP": mu_OHSsDSP,
                "mean_maxHkpHsDSP": mu_maxHkpHsDSP,
                "mean_TohkDSP": mu_TohkDSP,
                "mean_MHKsDSP": mu_MHKsDSP,
                "mean_MHSsDSP": mu_MHSsDSP,
                "mean_TmhksDSP": mu_TmhksDSP
            }
            datas.append(row)

        jitter_data=pd.DataFrame(datas)
        jitter_data.to_json("Rat_"+str(ratid)+"_Rec_"+str(18)+"_jitter_"+str(int(jitter*1000))+"ms_data.json",orient="records")
    