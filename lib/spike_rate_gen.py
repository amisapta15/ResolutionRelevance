import numpy as np
import pandas as pd
from cmcrameri import cm



ratid=21012
Path="../../Codes_5/data_extract/"
df=pd.read_json(Path+"Rat_"+str(ratid)+"_data_extracted.json");
# df.keys()

#Following is manully set

## The DID is the unique ID of a particular neuron and U_GID is the unique ID of a recording of the neuron
Neu_DID=262;U_GID=6353; 

print(df.query('N_DID==@Neu_DID')['U_GID'])
# From the above output, we can see that 
# the nu in the index (row of the dataframe) 
# and uni is the index of the list containing recordings/data of that neuron

# For example, if we want to get the U_GID of the neuron with N_DID=262 and U_GID=6353
nu=169;uni=5; ## N_DID=262, UGID=6353;


print("Unit ID: ",df.iloc[nu].U_GID[uni]) # This will give us the U_GID of the neuron with N_DID=262 and U_GID=6353



print("Rec start: ",df.iloc[nu].time_range[uni][0])
print("Rec end: ",df.iloc[nu].time_range[uni][-1])
print("Rec Duration: ",df.iloc[nu].time_range[uni][-1] - df.iloc[nu].time_range[uni][0] , df.iloc[nu].duration[uni])



######################################################################################
from scipy.interpolate import UnivariateSpline

def interp1(xpt, ypt, X):
    curv = UnivariateSpline(xpt, ypt, w=np.ones(len(xpt)), k=1, s=0.0, ext=0)
    return curv(X)
    # You can uncomment the following line to also get the first derivative
    # return curv(X), curv.derivative()(X, 1)
    
#Scalling
def scales(array, lower, upper):
    scaled_array = (np.array(array) - np.min(array)) / np.ptp(array)
    return lower + (upper - lower) * scaled_array
def opSP(spikes_spliced, zpos_x, zpos_y, zpos_t, rec_start, rec_end, binning_time):
    zposx = scales(zpos_x, 0, 150)  # To convert from m to cm, and scaling between 0 to 150 cm
    zposy = scales(zpos_y, 0, 150)  # To convert from m to cm, and scaling between 0 to 150 cm
    zpos_t = np.array(zpos_t) 
    #######################################################################
    post = zpos_t[(rec_start <= zpos_t) & (zpos_t <= rec_end)]  # Splicing the time stamps
    posx_ = zposx[(rec_start <= zpos_t) & (zpos_t <= rec_end)]  # Spliced X-axis stamps
    posy_ = zposy[(rec_start <= zpos_t) & (zpos_t <= rec_end)]  # Spliced Y-axis stamps
    #######################################################################
    times = np.arange(rec_start, rec_end + binning_time, binning_time)
    # print("Projecting from", len(post),"to this",len(times))
    #######################################################################
    posx = interp1(post, posx_, times)
    posy = interp1(post, posy_, times)
    #######################################################################
    binned_spikes, _ = np.histogram(spikes_spliced, bins=times)
    return posx, posy, binned_spikes.astype(int)

def info_count(zdata,zpos_x,zpos_y,zpos_t,rec_start,rec_end,binning_time,nbins):
    #Box Size = 1.5m X 1.5m  #150x150cm
    #Spatial Bin =150/nbins cm x 150/nbins cm
    #Total spatial bins= nbins*nbins

    #if nbins=30
    #Spatial Bin = 0.05m X 0.05m #5cmx5cm
    #Total spatial bins= nbins*nbins=900

    zdata=np.array(zdata)
    data = zdata[(rec_start <= zdata) & (zdata <= rec_end)]  # Splicing the data
    NSpikes = len(data)
    posx, posy, binned_spikes = opSP(data, zpos_x, zpos_y, zpos_t, rec_start, rec_end, binning_time)

    bins = np.linspace(0, 150, nbins+1)
    spkX = np.repeat(posx[:-1], binned_spikes) 
    spkY = np.repeat(posy[:-1], binned_spikes)  # Repeat the position data as per the number of spikes in that bin
    
    return NSpikes,posx,posy,spkX,spkY
#######################################################################################




import matplotlib as mpl
import matplotlib.pyplot as plt

nps,psx,psy,spsx,spsy=info_count(df.iloc[nu].u_spiketime[uni],df.iloc[nu].X[uni],df.iloc[nu].Y[uni],df.iloc[nu].t[uni],df.iloc[nu].time_range[uni][0],df.iloc[nu].time_range[uni][-1],0.01,30)
#binned with 0.01 secs (or 10 mili seconds) and 30 bins, see the function info_count for more details


#cmap=cm.batlow.resampled(512)
cmap=mpl.colormaps['inferno']
#cmap.set_bad(color='black')

stt=np.array(df.iloc[nu].pmap[uni]['count'])
stts=stt.astype(np.float64)
stts[stt == -1] = np.nan

f,(ax1,ax2)=plt.subplots(1,2,figsize=(20,8))
ax1.imshow(stts,aspect='auto',cmap=cmap,origin='lower')
ax2.scatter(psx,psy,marker='o',s=1)
ax2.scatter(spsx,spsy,marker='o',c='red',s=3)

plt.suptitle("Rat_"+str(ratid)+"    "+str(df.iloc[nu]['LOC'][1])+" Neu DID: "+str(df.iloc[nu].N_DID)+"    U_GID: "+str(df.iloc[nu].U_GID[uni])+"    Duration: "+str(np.round((df.iloc[nu].duration[uni]/60),2))+" minutes")

plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap), ax=ax1, orientation='vertical',label='Spike Infos')

plt.savefig("Rat_"+str(ratid)+"_"+str(df.iloc[nu]['LOC'][1])+"_Neu_DID_"+str(df.iloc[nu].N_DID)+"_U_GID_"+str(df.iloc[nu].U_GID[uni])+".png",dpi=300)