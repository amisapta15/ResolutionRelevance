using StatsBase
using Dierckx
#Dropings Singleton Dimensions
squzee(x)= dropdims(x, dims = tuple(findall(size(x) .== 1)...))
repmat(x,v)=vcat(fill.(x, v)...)
#Interpolation
function interp1(xpt,ypt,X)
    curv=Spline1D(xpt, ypt; w=ones(length(xpt)), k=1, bc="nearest", s=0.0)
    return curv(X)
    #return curv(X),derivative(curv,X,1)
end
#Scalling
function scales(array,lower,upper)
    #scaling to 0 to 1
    array_=squzee((array .- minimum(array)) / (maximum(array) .- minimum(array)));
    #scaling to lower to upper
    scaled_array = lower .+ ((upper-lower) .* array_);
    return scaled_array
end


function opSP(spikes_spliced,zpos_x,zpos_y,zpos_t,rec_start,rec_end,binning_time)
    zposx=scales(zpos_x .* 100,0,150) #To convert from m to cm , and scaling between 0 to 150 cm
    zposy=scales(zpos_y .* 100,0,150) #To convert from m to cm , and scaling between 0 to 150 cm
    #######################################################################
    post=zpos_t[rec_start.<=zpos_t.<=rec_end] # Splicing the time stamps
    posx_=zposx[rec_start.<=zpos_t.<=rec_end] # Spliced X-axis stamps
    posy_=zposy[rec_start.<=zpos_t.<=rec_end] # Spliced Y-axis stamps
    #######################################################################
    times=rec_start:binning_time:rec_end+binning_time;
    #print("Projecting from", size(post),"to this",size(collect(times)))
    #######################################################################
    posx=interp1(post,posx_, times);
    posy=interp1(post,posy_, times);
    #######################################################################
   
    binned_spikes=fit(Histogram, spikes_spliced,  times, closed=:left).weights
    return posx,posy,binned_spikes
end

function opHD(spikes_spliced,zHD_,zpos_t,rec_start,rec_end,binning_time)
    zHD=scales(zHD_,0,2*pi) #To convert from m to cm , and scaling between 0 to 150 cm
    ######################################################################
    post=zpos_t[rec_start.<=zpos_t.<=rec_end] # Splicing the time stamps
    HD_=zHD[rec_start.<=zpos_t.<=rec_end] 
    #######################################################################
    times=rec_start:binning_time:rec_end+binning_time;
    #print("Projecting from", size(post),"to this",size(collect(times)))
    #######################################################################
    HD=interp1(post,HD_, times);
    #######################################################################
    binned_spikes=fit(Histogram, spikes_spliced,  times, closed=:left).weights
return HD,binned_spikes
end


function SPinfo(zdata,zpos_x,zpos_y,zpos_t,rec_start,rec_end,binning_time,nbins)
    #Box Size = 1.5m X 1.5m  #150x150cm
    #Spatial Bin =150/nbins cm x 150/nbins cm
    #Total spatial bins= nbins*nbins

    #if nbins=30
    #Spatial Bin = 0.05m X 0.05m #5cmx5cm
    #Total spatial bins= nbins*nbins=900

    data=zdata[rec_start.<=zdata.<=rec_end] # Spilcing the data
    NSpikes=length(data)
    posx,posy,binned_spikes=opSP(data,zpos_x,zpos_y,zpos_t,rec_start,rec_end,binning_time)

    bins=range(0,stop=150,length=nbins+1) 
    spkX=repmat(posx[1:end-1],binned_spikes);
    spkY=repmat(posy[1:end-1],binned_spikes);

    # Fit a histogram with StatsBase
    h1 = fit(Histogram, (spkX,spkY),(bins,bins)).weights;
    h2 = fit(Histogram, (posx,posy),(bins,bins)).weights;

    ZRaw=h1 ./ (h2 .* binning_time);  ## No of Spikes in a bin / Time Spent In that Bin

    posPDF=h2/sum(h2);   ##OCupation Probablity ; time spend in that bin / Total session time
    
    #Ap=ZRaw.*posPDF       ## Mean Firing Rate
    #Ap[isnan.(Ap)].=0      ## Removing NaNs
    #lbarP=sum(Ap);         ## Summing over all Bins
    lbarP=NSpikes/(rec_end-rec_start);
    #println(sum(h2).*binning_time)
    
    #log Argument
    Cp=ZRaw ./ lbarP       ## Rate MAp / mean Rate 
    Cp[isnan.(Cp)].=1.0;   ## Removing NaNs
    Cp[Cp.<1.0].=1.0;      ##Making All Values < 1 = 1
    Bp=log2.(Cp)
    
    #Non Log Term
    Dp= ZRaw ./ lbarP
    Dp[isnan.(Dp)].=0.0;
    
    posIC= sum(posPDF .* Dp .* Bp);
    posIR= sum(posPDF .* Dp.* lbarP .* Bp);
return posIC,posIR
end

function HDinfo(zdata,zHD_,zpos_t,rec_start,rec_end,binning_time,hnbins)
    #Hd direction: 2pi - 360deg
    #angular bins= ~  360/hnbins deg per bin

    #if hnbins=50
    #angular bins= ~  7 deg per bin
    
    data=zdata[rec_start.<=zdata.<=rec_end] # Spilcing the data
    NSpikes=length(data)
    HD,binned_spikes=opHD(data,zHD_,zpos_t,rec_start,rec_end,binning_time)

    #hnbins=50
    hbins=range(0,stop=2*pi,length=hnbins+1) #ox Size Approx
    spkHD=repmat(HD[1:end-1],binned_spikes);
    # Fit a histogram
    h1HD = fit(Histogram, (spkHD),(hbins)).weights;      # No of spikes at each Bin
    h2HD = fit(Histogram, (HD),(hbins)).weights;      # No of  (X,Y) at each bin
    
    HDRaw=h1HD ./ (h2HD .* binning_time);  ## No of Spikes in a bin / Time Spent In that Bin

    hdPDF=h2HD/sum(h2HD);   ##OCupation Probablity ; time spend in that bin / Total session time
    
    
    #Ah=HDRaw.*hdPDF       ## Mean Firing Rate
    #Ah[isnan.(Ah)].=0      ## Removing NaNs
    #lbarH=sum(Ah);         ## Summing over all Bins
    lbarH=NSpikes/(rec_end-rec_start);
    
    #log Argument
    Ch=HDRaw ./ lbarH       ## Rate MAp / mean Rate 
    Ch[isnan.(Ch)].=1.0;   ## Removing NaNs
    Ch[Ch.<1.0].=1.0;      ##Making All Values < 1 = 1
    Bh=log2.(Ch)
    
    #Non Log Term
    Dh= HDRaw ./ lbarH
    Dh[isnan.(Dh)].=0.0;
    
    hdIC= sum(hdPDF .* Dh .* Bh);
    hdIR= sum(hdPDF .* Dh.*lbarH .* Bh);
return hdIC,hdIR
end