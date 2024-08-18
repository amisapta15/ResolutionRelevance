
LogSpace(starting,ending,no_of_points::Int)=10 .^ range( log10(starting), log10(ending), length = no_of_points)

function msrc(rec_start,rec_end,data::Array,binning_time)
    
    edges=rec_start:binning_time:rec_end+binning_time;
    counts=fit(Histogram, data,  edges, closed=:left).weights;
    
    #println(sum(counts))
    #return counts

    #Eleminating Zero-counts
    zero_loc=findall(counts.==0); #Finding Zero Locations
    deleteat!(counts,zero_loc); #Eleminating Zeros

    #Counting unique values and Frequencies
    o=countmap(counts); # Couting Unique appearences
    vals=collect(keys(o));freq=collect(values(o)) # Values & their Freqiencies
    
    # #for elements in zip(vals,freq); println(elements); end # See

    #Checking if total count and Freq/Values matches or not
    length(data) == sum(vals.*freq)  || throw(AssertionError("Value/Frequency coumputation error!"))

    #Calculating Hs & Hk
    M=sum(vals.*freq)
    hs= - sum(  ((freq.*vals)./M) .* log.(M,(vals./M)) )
    hk= - sum(  ((freq.*vals)./M) .* log.(M,((freq.*vals)./M)) )
    
    return [hs,hk]
end

function calculate_area(data_points)
    # calculates integral using the trapezoid rule
    return sum(0.5*(abs.(data_points[:,1][1:end-1]+data_points[:,1][2:end]))
        .*(abs.(data_points[:,2][1:end-1]-data_points[:,2][2:end])))
end

function cal_msr(spike_times,timerange,num_points)
    data=spike_times[timerange[1].<=spike_times.<=timerange[2]] # Spilcing the data
    duration=timerange[2]-timerange[1]
    
    Tbin=LogSpace(0.01,duration,num_points);
    result=vcat([msrc(timerange[1],timerange[2],data,tbins) for tbins in Tbin]'...);
    
    HofS=abs.(result[:,1]);
    HofK=abs.(result[:,2]);
    
    Tbin=vcat(0.00000000001,Tbin)
    HofS=vcat(1.0,HofS)
    HofK=vcat(0.0,HofK)
    
    old_msr=calculate_area(hcat(HofK[sortperm(HofS)],HofS[sortperm(HofS)]))
    idxxB=argmax(HofS .+ HofK)
    #CRatio=HofK[idxxB]
    #MTbin=Tbin[idxxB]
    #OMSR=old_msr;

    return old_msr,HofK[idxxB],Tbin[idxxB],HofK[argmax(HofK)]
end
