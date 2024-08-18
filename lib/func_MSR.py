import numpy as np
from collections import Counter
from scipy.integrate import trapezoid

def msrc(rec_start, rec_end, data, binning_time):
    edges = np.arange(rec_start, rec_end + binning_time, binning_time)
    counts, _ = np.histogram(data, bins=edges)

    # Eliminating Zero-counts
    zero_loc = np.where(counts == 0)[0]
    counts = np.delete(counts, zero_loc)

    # Counting unique values and Frequencies
    o = Counter(counts)
    #vals, freq = zip(*o.items())
    vals, freq = np.array(list(o.keys())), np.array(list(o.values()))
    # Calculating total number of spikes
    M = sum(vals *freq)

    # Checking if total count and Freq/Values match or not
    assert len(data) == M, "Value/Frequency computation error!"

    # Calculating Hs & Hk
    hs = - np.sum(
        ((freq * vals) / M) * (np.log(vals/M) / np.log(M)) # Entropy calculation on Base M
    )

    hk = - np.sum(
        ((freq * vals) / M) * (np.log( (freq * vals) /M) / np.log(M)) # Entropy calculation on Base M
    )
    return [hs, hk]


def calculate_area(data_points):
    return np.sum(0.5 * (np.abs(data_points[:-1, 0] + data_points[1:, 0])) *
                   (np.abs(data_points[:-1, 1] - data_points[1:, 1])))

def MSR(zdata, rec_start, rec_end, num_points):
    duration=rec_end-rec_start;
    # Calculating the log-spaced time bin lengths from 10 milisec to chunk-size mins
    Tbin=np.logspace(np.log10(0.01), np.log10(duration), num=num_points, endpoint=True, base=10.0)
    #Spiliting the data into chunks of chunk-size mins
    data_=np.array(zdata)
    datas = data_[(rec_start <= data_) & (data_ <= rec_end)]
    #Number of Spikes in the chunk
    #NSpikes = len(datas)
    
    #Calculating the Hs and Hk for each time bin and taking transpose
    result = np.array([msrc(rec_start, rec_end, datas, tbins) for tbins in Tbin]).T
    HofS = np.abs(result[0])
    HofK = np.abs(result[1])
    
    
    # Adding the theoretical start and end points to the Hs and Hk arrays
    time_bins = np.concatenate(([0.00000000001], Tbin, [duration]))
    HofS = np.concatenate(([1.0], HofS, [0.0]))
    HofK = np.concatenate(([0.0], HofK, [0.0]))

    # Manually perform trapezoidal integration
    #msr = calculate_area(np.column_stack((HofK[np.argsort(HofS)], HofS[np.argsort(HofS)])))
    msr=trapezoid(HofK[np.argsort(HofS)], HofS[np.argsort(HofS)])
    
    idxxA = np.argmax(HofK) # Index of the maximum value of Hk
    idxxB = np.argmax(HofS + HofK) # Index of the maximum value of Hs + Hk i.e optimal point

    # Calculating the MHK and TOHK and other quantities
    OHKs = HofK[idxxB]
    OHSs = HofS[idxxB]
    maxHkpHs = OHKs + OHSs
    Tohks = time_bins[idxxB]
    MHKs = HofK[idxxA]
    MHSs=HofS[idxxA]
    Tmhks = time_bins[idxxA]

    return len(datas),msr,OHKs,OHSs,maxHkpHs,Tohks,MHKs,MHSs,Tmhks

def poissonification(zdata, rec_start, rec_end,dt):
    data_=np.array(zdata)
    datas = data_[(rec_start <= data_) & (data_ <= rec_end)]
    #Number of Spikes in the chunk
    #NSpikes = len(datas)

    times=np.arange(rec_start,rec_end,dt);
    binarized_datas=np.histogram(datas,bins=times)[0]
    rng = np.random.default_rng()

    avg_firing_rate = np.sum(binarized_datas) / (len(binarized_datas)*dt)
    poss_st = rng.random((len(binarized_datas),)) < (avg_firing_rate  *dt)

    return times[:-1][poss_st]

    #?? Random Shuffling
    # poss2_st=binarized_datas.astype(bool).copy(); 
    # for _ in range(1000):
    #     rng.shuffle(poss2_st)

    # return times[:-1][poss_st],times[:-1][poss2_st]
  
