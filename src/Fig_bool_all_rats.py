
import pandas as pd
import numpy as np
from itertools import combinations,product
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
mpl.rcParams.update(mpl.rcParamsDefault)

divs = lambda x, y: np.divide(np.array([x], dtype=float), np.array([y], dtype=float), out=np.full_like(np.array([x], dtype=float), np.nan), where=np.array([y], dtype=float)!=0)

# import warnings
# warnings.filterwarnings("ignore")

mss=50
fontssize=50
mpl.rcParams.update({
    'figure.figsize': (23.6,12.6),
    'font.family': 'serif',
    'font.serif': ['Liberation Sans'],  # Add or remove font names as needed
    'font.size': fontssize,  # Adjust as needed
    #'font.weight': 'bold',
    "svg.fonttype": 'none',
    'text.usetex': False,
    'axes.linewidth' : 2,
    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{amsfonts} \usepackage{cmbright}',
    'xtick.labelsize' : fontssize, # fontsize of the x tick labels
    'ytick.labelsize' : fontssize # fontsize of the y tick labels
})
ticks_font = font_manager.FontProperties(family='Liberation Sans', style='normal',
    size=fontssize, weight='bold', stretch='normal')


def avg_rows(df,operation,rat):
   op_tmps=df.query('OP=="'+operation+'"')[['N1_NeuID', 'N1_DID', 'N1_GID', 'N2_NeuID', 'N2_DID', 'N2_GID', 
         'Nspikes', 'MSR', 'MHK', 'MHS', 'dt_MHK', 'OHK', 'OHS', 'dt_OHK','max_HSHK']].mean().to_frame().T
   op_tmps=op_tmps.astype({'N1_NeuID':int, 'N1_DID':int, 'N1_GID':int, 'N2_NeuID':int, 'N2_DID':int, 'N2_GID':int})
   op_tmps.insert(0, 'Rat_ID', rat)
   assert np.array_equal(df.U_LOC.values,df.N_LOC.values)
   op_tmps.insert(1, 'LOC', df.U_LOC.values[0])
   sdd=df.query('OP=="NA"')[['MSR', 'MHK', 'MHS', 'dt_MHK', 'OHK', 'OHS', 'dt_OHK','max_HSHK']].max().to_frame().T
   sdd.rename(columns={col: 'Rmax_'+col for col in sdd.columns}, inplace=True)
   op_tmps = pd.concat([op_tmps,sdd], axis=1)
   return op_tmps



## Activate for Single recording plots
rats=[20382,24101,21012,22295,20630,22098,23783,24116]
path="../data_bool/"
chunksize=20 #in minutes

for quant in ['MSR','OHK','OHS']:
    if quant=='MSR':
        qmin=0.22;qmax=0.31;
        xqmax=0.31;xqmin=0.25;
        yqmax=0.32;yqmin=-0.01;
    elif quant=='OHK':
        qmin=0.1;qmax=0.6;
        xqmax=0.6;xqmin=0.1;
        yqmax=0.6;yqmin=-0.05;
    elif quant=='OHS':
        qmin=0.5;qmax=1.1;
        xqmax=1;xqmin=0.5;
        yqmax=1;yqmin=-0.05;
    
    
    for rat in rats:
        df=pd.read_json(path+"Rat_"+str(rat)+"_BOOLop_Nresrel_data.json")
        df = df.fillna(0)

        df_AND=pd.DataFrame(); df_OR=pd.DataFrame(); df_XOR=pd.DataFrame()
        # For Poisson spike trains
        df_possAND=pd.DataFrame(); df_possOR=pd.DataFrame(); df_possXOR=pd.DataFrame()
        for i,j in df.groupby(['N1_DID','N2_DID']).groups:
                tmps=df.query('N1_DID==' + str(i) + '& N2_DID==' + str(j))
                df_AND=pd.concat([df_AND,avg_rows(tmps,"AND",rat)], ignore_index=True)
                df_OR=pd.concat([df_OR,avg_rows(tmps,"OR",rat)], ignore_index=True)
                df_XOR=pd.concat([df_XOR,avg_rows(tmps,"XOR",rat)], ignore_index=True)
                # For Poisson spike trains
                df_possAND=pd.concat([df_possAND,avg_rows(tmps,"possAND",rat)], ignore_index=True)
                df_possOR=pd.concat([df_possOR,avg_rows(tmps,"possOR",rat)], ignore_index=True)
                df_possXOR=pd.concat([df_possXOR,avg_rows(tmps,"possXOR",rat)], ignore_index=True)

        # Specify the directory path
        directory = "../figures/bool_rats/"+str(rat)+"/"

        # Check if the directory already exists
        if not os.path.exists(directory):
            # Create the directory
            os.makedirs(directory)
            print("Directory created successfully!")
        else:
            print("Directory already exists!")
            
        assert np.array_equal(df_OR.N1_GID.values,df_XOR.N1_GID.values)
        assert np.array_equal(df_OR.N2_GID.values,df_XOR.N2_GID.values)
        # For Poisson spike trains
        assert np.array_equal(df_possOR.N1_GID.values,df_possXOR.N1_GID.values)
        assert np.array_equal(df_possOR.N2_GID.values,df_possXOR.N2_GID.values)

        print("Plotting for Rat: ",rat)
        ################################################################################
        fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='col',figsize=(50,30))

        ax1.scatter(df_OR.query('LOC=="CA1"')[quant].values,df_XOR.query('LOC=="CA1"')[quant].values, label='CA1 Neurons', marker='D', facecolors='none',edgecolors='r',s=mss)

        ax1.scatter(df_possOR.query('LOC=="CA1"')[quant].values,df_possXOR.query('LOC=="CA1"')[quant].values, label='poss CA1 Neurons',  marker='D', facecolors='none',edgecolors='grey',s=mss)

        ## Diagonal Line
        ax1.plot([0, qmax], [0, qmax],color="k",ls="-.", lw=3) 

        ax1.set_xlim([qmin,qmax]);ax1.set_ylim([qmin,qmax])
        # ax1.tick_params(axis='both', which='both', bottom=True, top=False, left=True, right=False, labelbottom=False,labelleft=True,rotation=0)

        ax1.set_xlabel(quant+r' - $OR\vee_{i,j}$');
        ax1.set_ylabel(quant+r' - $XOR\oplus_{i,j}$');

        ################################################################################

        ax3.scatter(df_OR.query('LOC=="SUB"')[quant].values,df_XOR.query('LOC=="SUB"')[quant].values,
                label='SUB Neurons', marker='o', facecolors='none',edgecolors='b',s=mss)

        #Poission Spike Trains
        ax3.scatter(df_possOR.query('LOC=="SUB"')[quant].values,df_possXOR.query('LOC=="SUB"')[quant].values, label='poss SUB Neurons',  marker='o', facecolors='none',edgecolors='grey',s=mss)


        ## Diagonal Line
        ax3.plot([0, qmax], [0, qmax],color="k",ls="-.", lw=3) 
        ax3.set_xlim([qmin,qmax]);ax1.set_ylim([qmin,qmax])

        # ax3.tick_params(axis='both', which='both', bottom=True, top=False, left=True, right=False, labelbottom=True,labelleft=True,rotation=0)

        ax3.set_xlabel(quant+r' - $OR\vee_{i,j}$');
        ax3.set_ylabel(quant+r' - $XOR\oplus_{i,j}$');

        ################################################################################

        ## Diagonal Line#
        ax2.plot([0, yqmax], [0, yqmax],color="k",ls="-.", lw=3) 

        clr='brown'
        ax2.set_xlim([xqmin,xqmax]);ax2.set_ylim([yqmin,yqmax]);

            
        ax2.set_xlabel(r'max $\{'+quant+r'_i,'+quant+r'_j\}$');

        ax2.scatter(df_AND.query('LOC=="CA1"')["Rmax_"+quant].values,df_AND.query('LOC=="CA1"')[quant].values,
                label='ANDed CA1 Neurons', marker='D', facecolors='none',edgecolors=clr,s=mss)

        ax2.scatter(df_possAND.query('LOC=="CA1"')["Rmax_"+quant].values,df_possAND.query('LOC=="CA1"')[quant].values,
                label='possANDed CA1 Neurons', marker='+',s=mss)


        ax2.set_ylabel(quant+r' - $AND\wedge_{i,j}$');
        ax2.tick_params(axis='y',labelcolor=clr)
        ax2.yaxis.label.set_color(clr)


        ax2_twin = ax2.twinx()
        clr='green'
            
        ax2_twin.scatter(df_OR.query('LOC=="CA1"')["Rmax_"+quant].values,df_OR.query('LOC=="CA1"')[quant].values,
                label='ORed CA1 Neurons', marker='D', color=clr,s=mss/2,alpha=0.3)

        ax2_twin.scatter(df_possOR.query('LOC=="CA1"')["Rmax_"+quant].values,df_possOR.query('LOC=="CA1"')[quant].values,label='possORed CA1 Neurons', marker='x',s=mss/2,alpha=0.3)


        ax2_twin.set_ylabel(quant+r' - $OR\vee_{i,j}$');
        ax2_twin.tick_params(axis='y', labelcolor=clr)
        ax2_twin.yaxis.label.set_color(clr)
        ax2_twin.set_ylim([yqmin,yqmax]);

        ## Activate if plotting for all recording 

        fracC_ad_or=divs(df_OR.query('LOC=="CA1" and '+quant+'>Rmax_'+quant).shape[0],df_OR.query('LOC=="CA1"').shape[0])
        fracC_bd_or=divs(df_OR.query('LOC=="CA1" and '+quant+'<Rmax_'+quant).shape[0],df_OR.query('LOC=="CA1"').shape[0])

        fracC_ad_and=divs(df_AND.query('LOC=="CA1" and '+quant+'>Rmax_'+quant).shape[0],df_AND.query('LOC=="CA1"').shape[0])
        fracC_bd_and=divs(df_AND.query('LOC=="CA1" and 0<'+quant+'<Rmax_'+quant).shape[0],df_AND.query('LOC=="CA1"').shape[0])
        fracC_zero=divs(df_AND.query('LOC=="CA1" and '+quant+'==0').shape[0],df_AND.query('LOC=="CA1"').shape[0])

        ax2.text(0.1, 0.5,   
                r'$OR^{ad}_\vee$ = %.2f%%' % (100*fracC_ad_or)+' '+r'$OR^{bd}_\vee$ = %.1f%%' %(100*fracC_bd_or)+'\n'+
                r'$AND^{ad}_\wedge$ = %.2f%%'% (100*fracC_ad_and)+' '+r'$AND^{bd}_\wedge$ = %.1f%%'%(100*fracC_bd_and)+'\n'+
                r'$AND^{0,na}_\wedge$ = %.2f%%'% (100*fracC_zero),
                transform=ax2.transAxes, fontsize=fontssize, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5, linewidth=0))
        ################################################################################


        ## Diagonal Line#
        ax4.plot([0, yqmax], [0, yqmax],color="k",ls="-.", lw=3) 

        clr='brown'
        ax4.set_xlim([xqmin,xqmax]);ax4.set_ylim([yqmin,yqmax]);

            
        ax4.set_xlabel(r'max $\{'+quant+r'_i,'+quant+r'_j\}$');

        ax4.scatter(df_AND.query('LOC=="SUB"')["Rmax_"+quant].values,df_AND.query('LOC=="SUB"')[quant].values,
                label='ANDed SUB Neurons', marker='D', facecolors='none',edgecolors=clr,s=mss)

        ax4.scatter(df_possAND.query('LOC=="SUB"')["Rmax_"+quant].values,df_possAND.query('LOC=="SUB"')[quant].values,
                label='possANDed SUB Neurons', marker='*',s=mss)

        ax4.set_ylabel(quant+r' - $AND\wedge_{i,j}$');
        ax4.tick_params(axis='y',labelcolor=clr)
        ax4.yaxis.label.set_color(clr)


        ax4_twin = ax4.twinx()
        clr='green'
            
        ax4_twin.scatter(df_OR.query('LOC=="SUB"')["Rmax_"+quant].values,df_OR.query('LOC=="SUB"')[quant].values,
                label='ORed SUB Neurons', marker='D', color=clr,s=mss/2,alpha=0.3)

        ax4_twin.scatter(df_possOR.query('LOC=="SUB"')["Rmax_"+quant].values,df_possOR.query('LOC=="SUB"')[quant].values,label='possORed SUB Neurons', marker='^',s=mss,alpha=0.3)

        ax4_twin.set_ylabel(quant+r' - $OR\vee_{i,j}$');
        ax4_twin.tick_params(axis='y', labelcolor=clr)
        ax4_twin.yaxis.label.set_color(clr)
        ax4_twin.set_ylim([yqmin,yqmax]);

        ## Activate if plotting for all recording 

        fracS_ad_or=divs(df_OR.query('LOC=="SUB" and '+quant+'>Rmax_'+quant).shape[0],df_OR.query('LOC=="SUB"').shape[0])
        fracS_bd_or=divs(df_OR.query('LOC=="SUB" and '+quant+'<Rmax_'+quant).shape[0],df_OR.query('LOC=="SUB"').shape[0])

        fracS_ad_and=divs(df_AND.query('LOC=="SUB" and '+quant+'>Rmax_'+quant).shape[0],df_AND.query('LOC=="SUB"').shape[0])
        fracS_bd_and=divs(df_AND.query('LOC=="SUB" and 0<'+quant+'<Rmax_'+quant).shape[0],df_AND.query('LOC=="SUB"').shape[0])
        fracS_zero=divs(df_AND.query('LOC=="SUB" and '+quant+'==0').shape[0],df_AND.query('LOC=="SUB"').shape[0])


        ax4.text(0.1, 0.5,   
                r'$OR^{ad}_\vee$ = %.2f%%' % (100*fracS_ad_or)+' '+r'$OR^{bd}_\vee$ = %.1f%%' %(100*fracS_bd_or)+'\n'+
                r'$AND^{ad}_\wedge$ = %.2f%%'% (100*fracS_ad_and)+' '+r'$AND^{bd}_\wedge$ = %.1f%%'%(100*fracS_bd_and)+'\n'+
                r'$AND^{0,na}_\wedge$ = %.2f%%'% (100*fracS_zero),
                transform=ax4.transAxes, fontsize=fontssize, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5, linewidth=0))


        handles, labels = [], []
        for ax in [ax1,ax2, ax2_twin,ax3,ax4,ax4_twin]:  # Add the axes that contain the plots you want in the legend
                for h, l in zip(*ax.get_legend_handles_labels()):
                    handles.append(h)
                    labels.append(l)

        # Create a single legend for the whole figure
        fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.99), ncol=4,fontsize=fontssize)

        fig.suptitle("Rat "+str(rat), x=0.1, y=1, horizontalalignment='left', verticalalignment='top', fontsize = 1.2*fontssize,color='red')
        fig.subplots_adjust(wspace=0.15,hspace=0.08)

            



        fig.savefig(directory+"Rat_"+str(rat)+"_"+quant+".png",bbox_inches='tight',dpi=300)

        plt.close(fig)

    