import numpy as np
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

def opHD(spikes_spliced, zHD_, zpos_t, rec_start, rec_end, binning_time):
    zHD = scales(zHD_, 0, 2*np.pi) 
    zpos_t = np.array(zpos_t) 
    ######################################################################
    post = zpos_t[(rec_start<=zpos_t)&(zpos_t<=rec_end)] # Splicing the time stamps
    HDs_ = zHD[(rec_start<=zpos_t)&(zpos_t<=rec_end)] 
    #######################################################################
    times = np.arange(rec_start,rec_end+binning_time,binning_time);
    #print("Projecting from", len(post),"to this",len(times))
    #######################################################################
    HDs = interp1(post,HDs_, times);
    #######################################################################
    binned_spikes, _ = np.histogram(spikes_spliced, bins=times)
    return HDs,binned_spikes.astype(int)

def SPinfo(zdata,zpos_x,zpos_y,zpos_t,rec_start,rec_end,binning_time,nbins):
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

    # Fit a histogram
    hist, _, _ = np.histogram2d(spkX, spkY, bins=[bins, bins]) # count the number of spikes in each bin
    hist2, _, _ = np.histogram2d(posx, posy, bins=[bins, bins]) #count the visits at each bin divided by time-bin = time spent in that bin  

    # hist, _, _,_ = binned_statistic_2d(spkX, spkY, values=None, statistic='count', bins=[bins, bins])
    # hist2, _, _, _ = binned_statistic_2d(posx, posy, values=None, statistic='count', bins=[bins, bins])
    
    ZRaw = np.divide( hist, np.where((hist2 * binning_time) != 0, (hist2 * binning_time), np.nan))

    posPDF = hist2 / np.sum(hist2)   ## Occupation Probability; time spent in that bin / Total session time

    lbarP = NSpikes / (rec_end - rec_start)

    # log Argument
    Cp = ZRaw / lbarP  ## Rate MAp / mean Rate
    Cp[np.isnan(Cp)] = 1.0;   ## Removing NaNs
    Cp[Cp < 1.0] = 1.0;      ## Making All Values < 1 = 1
    Bp = np.log2(Cp)

    # Non Log Term
    Dp = ZRaw / lbarP
    Dp[np.isnan(Dp)] = 0.0;

    posIC = np.sum(posPDF * Dp * Bp)
    posIR = np.sum(posPDF * Dp * lbarP * Bp )
    return posIC,posIR

def HDinfo(zdata,zHD_,zpos_t,rec_start,rec_end,binning_time,hnbins):
    #Hd direction: 2pi - 360deg
    #angular bins= ~  360/hnbins deg per bin

    #if hnbins=50
    #angular bins= ~  7 deg per bin

    zdata=np.array(zdata)
    data = zdata[(rec_start <= zdata) & (zdata <= rec_end)]  # Splicing the data
    NSpikes = len(data)
    HDs, binned_spikes = opHD(data,zHD_, zpos_t, rec_start, rec_end, binning_time)

    hbins = np.linspace(0, 2*np.pi, hnbins+1)
    spkHDs = np.repeat(HDs[:-1], binned_spikes) 

    # Fit a histogram
    hist, _ = np.histogram(spkHDs, bins=hbins)
    hist2, _ = np.histogram(HDs,bins=hbins)

    # hist, _, _,_ = binned_statistic_2d(spkX, spkY, values=None, statistic='count', bins=[bins, bins])
    # hist2, _, _, _ = binned_statistic_2d(posx, posy, values=None, statistic='count', bins=[bins, bins])
    
    ZRaw = np.divide( hist, np.where((hist2 * binning_time) != 0, (hist2 * binning_time), np.nan))

    posPDF = hist2 / np.sum(hist2)   ## Occupation Probability; time spent in that bin / Total session time

    lbarP = NSpikes / (rec_end - rec_start)

    # log Argument
    Cp = ZRaw / lbarP  ## Rate MAp / mean Rate
    Cp[np.isnan(Cp)] = 1.0;   ## Removing NaNs
    Cp[Cp < 1.0] = 1.0;      ## Making All Values < 1 = 1
    Bp = np.log2(Cp)

    # Non Log Term
    Dp = ZRaw / lbarP
    Dp[np.isnan(Dp)] = 0.0;

    HDIC = np.sum(posPDF * Dp * Bp)
    HDIR = np.sum(posPDF * Dp * lbarP * Bp )
    return HDIC,HDIR
