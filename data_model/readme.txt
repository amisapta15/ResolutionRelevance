Results are for 10ms bin size and 40ms bin size.

Figures are saved both in jpg and fig format with 10/40ms indicated in name. 

I used subplots_main.m to plot the figures with 3 subplots.

Yasser suggestion: make a single plot with 6 subplots(2columns). First column is results of 10 ms and second column results of 40 ms. 


Data files are cells with 2 elements. Each sub element contains d number of doubles with size 1*20 which are results of 20 random realisations. d is the number of sub populations we considered which is 7 for CA1 and 8 for SUB. Each element of main array is results of recording defined by file NAMES.mat. 