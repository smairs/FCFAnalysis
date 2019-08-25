import matplotlib.pyplot as plt
import numpy as np

def plot_by_tauxairmass_2lists(FCFP_850_before_newfilters1,\
                               FCFP_err_850_before_newfilters1,\
                               FCFP_850_after_newfilters1,\
                               FCFP_err_850_after_newfilters1,\
                               FCFA_850_before_newfilters1,\
                               FCFA_err_850_before_newfilters1,\
                               FCFA_850_after_newfilters1,\
                               FCFA_err_850_after_newfilters1,\
                               FCFP_450_before_SMUtrouble1,\
                               FCFP_err_450_before_SMUtrouble1,\
                               FCFP_450_after_SMUtrouble1,\
                               FCFP_err_450_after_SMUtrouble1,\
                               FCFA_450_before_SMUtrouble1,\
                               FCFA_err_450_before_SMUtrouble1,\
                               FCFA_450_after_SMUtrouble1,\
                               FCFA_err_450_after_SMUtrouble1,\
                               FCFP_850_before_newfilters2,\
                               FCFP_err_850_before_newfilters2,\
                               FCFP_850_after_newfilters2,\
                               FCFP_err_850_after_newfilters2,\
                               FCFA_850_before_newfilters2,\
                               FCFA_err_850_before_newfilters2,\
                               FCFA_850_after_newfilters2,\
                               FCFA_err_850_after_newfilters2,\
                               FCFP_450_before_SMUtrouble2,\
                               FCFP_err_450_before_SMUtrouble2,\
                               FCFP_450_after_SMUtrouble2,\
                               FCFP_err_450_after_SMUtrouble2,\
                               FCFA_450_before_SMUtrouble2,\
                               FCFA_err_450_before_SMUtrouble2,\
                               FCFA_450_after_SMUtrouble2,\
                               FCFA_err_450_after_SMUtrouble2,\
                               TauTimesAirmassBins,\
                               badsourcelists,\
                               constrained='no'):

    plt.clf()
    
    FCFP_850_nominal     = 537.0
    FCFP_850_nominal_err = 26.0/537.0
    
    FCFP_450_nominal     = 491.0
    FCFP_450_nominal_err = 67.0/491.0
    
    FCFA_850_nominal     = 2.34
    FCFA_850_nominal_err = 0.08/2.34
    
    FCFA_450_nominal     = 4.71
    FCFA_450_nominal_err = 0.5/4.71
    
    bincolours             = ['r','y','g','b','violet']
    TautimesAirmass_before = []
    TautimesAirmass_after  = []
    for i in range(len(TauTimesAirmassBins)-1):
        TautimesAirmass_before.append(TauTimesAirmassBins[i]+(TauTimesAirmassBins[i+1]-TauTimesAirmassBins[i])/4)
        TautimesAirmass_after.append(TauTimesAirmassBins[i]+(TauTimesAirmassBins[i+1]-TauTimesAirmassBins[i])/1.3333333)
    
    TautimesAirmass_before = np.array(TautimesAirmass_before)
    TautimesAirmass_after  = np.array(TautimesAirmass_after)

    # FCF P 850 #
    TauTimesAirmassDummy = 0
    for i in range(len(TauTimesAirmassBins)-1):
        plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours[TauTimesAirmassDummy])
        TauTimesAirmassDummy = TauTimesAirmassDummy+1
   
    split_up_markers = 0.005
 
    plt.errorbar(TautimesAirmass_before-split_up_markers,FCFP_850_before_newfilters1,yerr=FCFP_850_before_newfilters1*FCFP_err_850_before_newfilters1,linestyle='none',marker='o',color='gray',label='Before New Filters (Set 1)',markersize=10)
    plt.errorbar(TautimesAirmass_after-split_up_markers,FCFP_850_after_newfilters1,yerr=FCFP_850_after_newfilters1*FCFP_err_850_after_newfilters1,linestyle='none',marker='^',color='k',label='After New Filters (Set 1)',markersize=10)
    plt.errorbar(TautimesAirmass_before+split_up_markers,FCFP_850_before_newfilters2,yerr=FCFP_850_before_newfilters2*FCFP_err_850_before_newfilters2,linestyle='none',marker='s',label='Before New Filters (Set 2)',color='gray',markersize=10)
    plt.errorbar(TautimesAirmass_after+split_up_markers,FCFP_850_after_newfilters2,yerr=FCFP_850_after_newfilters2*FCFP_err_850_after_newfilters2,linestyle='none',marker='x',label='After New Filters (Set 2)',color='k',markersize=10)
    plt.xticks(TauTimesAirmassBins)
    plt.legend(loc='upper left')
    plt.ylabel('Peak FCF, 850 microns (Jy/beam/pW)')
    plt.xlabel('Tau225 x Airmass')
    plt.suptitle("Set 1: No {source}".format(', '.join(source for source in badsourcelist[0]))+"; Set 2:  No {source}".format(', '.join(source for source in badsourcelist[1])),fontsize=10)
    if constrained == 'no':
        plt.savefig('Figures/FCFP_by_TauxAM_850_SOURCESETS.png',format='png',dpi=300)
    if constrained == 'yes':
        plt.savefig('Figures/FCFP_by_TauxAM_850_constrained_SOURCESETS.png',format='png',dpi=300)
    plt.clf()
    
    
    # FCF A 850 #
    
    TauTimesAirmassDummy = 0
    for i in range(len(TauTimesAirmassBins)-1):
        plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours[TauTimesAirmassDummy])
        TauTimesAirmassDummy = TauTimesAirmassDummy+1
    
    plt.errorbar(TautimesAirmass_before-split_up_markers,FCFA_850_before_newfilters1,yerr=FCFA_850_before_newfilters1*FCFA_err_850_before_newfilters1,linestyle='none',marker='o',color='gray',label='Before New Filters (Set 1)',markersize=10)
    plt.errorbar(TautimesAirmass_after-split_up_markers,FCFA_850_after_newfilters1,yerr=FCFA_850_after_newfilters1*FCFA_err_850_after_newfilters1,linestyle='none',marker='^',color='k',label='After New Filters (Set 1)',markersize=10)
    plt.errorbar(TautimesAirmass_before+split_up_markers,FCFA_850_before_newfilters2,yerr=FCFA_850_before_newfilters2*FCFA_err_850_before_newfilters2,linestyle='none',marker='s',label='Before New Filters (Set 2)',color='gray',markersize=10)
    plt.errorbar(TautimesAirmass_after+split_up_markers,FCFA_850_after_newfilters2,yerr=FCFA_850_after_newfilters2*FCFA_err_850_after_newfilters2,linestyle='none',marker='x',label='After New Filters (Set 2)',color='k',markersize=10)
    plt.xticks(TauTimesAirmassBins)
    plt.legend(loc='upper left')
    plt.ylabel('Arcsec FCF, 850 microns (Jy/arcsec^2/pW)')
    plt.xlabel('Tau225 x Airmass')
    plt.suptitle("Set 1: No {source}".format(', '.join(source for source in badsourcelist[0]))+"; Set 2:  No {source}".format(', '.join(source for source in badsourcelist[1])),fontsize=10)
    if constrained == 'no':
        plt.savefig('Figures/FCFA_by_TauxAM_850_SOURCESETS.png',format='png',dpi=300)
    if constrained == 'yes':
        plt.savefig('Figures/FCFA_by_TauxAM_850_constrained_SOURCESETS.png',format='png',dpi=300)
    plt.clf()
    
    # FCF P 450 #
    
    TauTimesAirmassDummy = 0
    for i in range(len(TauTimesAirmassBins)-1):
        plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours[TauTimesAirmassDummy])
        TauTimesAirmassDummy = TauTimesAirmassDummy+1
    
    plt.errorbar(TautimesAirmass_before-split_up_markers,FCFP_450_before_SMUtrouble1,yerr=FCFP_450_before_SMUtrouble1*FCFP_err_450_before_SMUtrouble1,linestyle='none',marker='o',color='gray',label='Before SMU Fix (Set 1)',markersize=10)
    plt.errorbar(TautimesAirmass_after-split_up_markers,FCFP_450_after_SMUtrouble1,yerr=FCFP_450_after_SMUtrouble1*FCFP_err_450_after_SMUtrouble1,linestyle='none',marker='^',color='k',label='After SMU Fix (Set 1)',markersize=10)
    plt.errorbar(TautimesAirmass_before+split_up_markers,FCFP_450_before_SMUtrouble2,yerr=FCFP_450_before_SMUtrouble2*FCFP_err_450_before_SMUtrouble2,linestyle='none',marker='s',label='Before SMU Fix (Set 2)',color='gray',markersize=10)
    plt.errorbar(TautimesAirmass_after+split_up_markers,FCFP_450_after_SMUtrouble2,yerr=FCFP_450_after_SMUtrouble2*FCFP_err_450_after_SMUtrouble2,linestyle='none',marker='x',label='After SMU Fix (Set 2)',color='k',markersize=10)
    plt.xticks(TauTimesAirmassBins)
    plt.legend(loc='upper left')
    plt.ylabel('Peak FCF, 450 microns (Jy/beam/pW)')
    plt.xlabel('Tau225 x Airmass')
    plt.suptitle("Set 1: No {source}".format(', '.join(source for source in badsourcelist[0]))+"; Set 2:  No {source}".format(', '.join(source for source in badsourcelist[1])),fontsize=10)
    if constrained == 'no':
        plt.savefig('Figures/FCFP_by_TauxAM_450_SOURCESETS.png',format='png',dpi=300)
    if constrained == 'yes':
        plt.savefig('Figures/FCFP_by_TauxAM_450_constrained_SOURCESETS.png',format='png',dpi=300)
    plt.clf()
    
    # FCF A 450 #
    
    TauTimesAirmassDummy = 0
    for i in range(len(TauTimesAirmassBins)-1):
        plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours[TauTimesAirmassDummy])
        TauTimesAirmassDummy = TauTimesAirmassDummy+1
    
    plt.errorbar(TautimesAirmass_before-split_up_markers,FCFA_450_before_SMUtrouble1,yerr=FCFA_450_before_SMUtrouble1*FCFA_err_450_before_SMUtrouble1,linestyle='none',marker='o',color='gray',label='Before SMU Fix (Set 1)',markersize=10)
    plt.errorbar(TautimesAirmass_after-split_up_markers,FCFA_450_after_SMUtrouble1,yerr=FCFA_450_after_SMUtrouble1*FCFA_err_450_after_SMUtrouble1,linestyle='none',marker='^',color='k',label='After SMU Fix (Set 1)',markersize=10)
    plt.errorbar(TautimesAirmass_before+split_up_markers,FCFA_450_before_SMUtrouble2,yerr=FCFA_450_before_SMUtrouble2*FCFA_err_450_before_SMUtrouble2,linestyle='none',marker='s',label='Before SMU Fix (Set 2)',color='gray',markersize=10)
    plt.errorbar(TautimesAirmass_after+split_up_markers,FCFA_450_after_SMUtrouble2,yerr=FCFA_450_after_SMUtrouble2*FCFA_err_450_after_SMUtrouble2,linestyle='none',marker='x',label='After SMU Fix (Set 2)',color='k',markersize=10)
    plt.xticks(TauTimesAirmassBins)
    plt.legend(loc='upper left')
    plt.ylabel('Arcsec FCF, 450 microns (Jy/arcsec^2/pW)')
    plt.xlabel('Tau225 x Airmass')
    plt.suptitle("Set 1: No {source}".format(', '.join(source for source in badsourcelist[0]))+"; Set 2:  No {source}".format(', '.join(source for source in badsourcelist[1])),fontsize=10)
    if constrained == 'no':
        plt.savefig('Figures/FCFA_by_TauxAM_450_SOURCESETS.png',format='png',dpi=300)
    if constrained == 'yes':
        plt.savefig('Figures/FCFA_by_TauxAM_450_constrained_SOURCESETS.png',format='png',dpi=300)
    plt.clf()

    if constrained == 'no':
        print('\n\nInformation To Plot -- SET 1, No ',badsourcelist[0],' -- NOT CONSTRAINED\n\n')
    if constrained == 'yes':
        print('\n\nInformation To Plot -- SET 1, No ',badsourcelist[0],' -- CONSTRAINED\n\n')
    print('x_before = ',list(TautimesAirmass_before-split_up_markers))
    print('x_after = ',list(TautimesAirmass_after-split_up_markers))
    print('\n\nFCF Arcsec\n\n')
    print('FCFA_450_before     = np.array(',list(FCFA_450_before_SMUtrouble1),')')
    print('FCFA_450_before_err = np.array(',list(FCFA_450_before_SMUtrouble1*FCFA_err_450_before_SMUtrouble1),')')
    print('FCFA_450_after      = np.array(',list(FCFA_450_after_SMUtrouble1),')')
    print('FCFA_450_after_err  = np.array(',list(FCFA_450_after_SMUtrouble1*FCFA_err_450_after_SMUtrouble1),')')
    print('FCFA_850_before     = np.array(',list(FCFA_850_before_newfilters1),')')
    print('FCFA_850_before_err = np.array(',list(FCFA_850_before_newfilters1*FCFA_err_850_before_newfilters1),')')
    print('FCFA_850_after      = np.array(',list(FCFA_850_after_newfilters1),')')
    print('FCFA_850_after_err  = np.array(',list(FCFA_850_after_newfilters1*FCFA_err_850_after_newfilters1),')')
    print('\n\nFCF Peak\n\n')
    print('FCFP_450_before     = np.array(',list(FCFP_450_before_SMUtrouble1),')')
    print('FCFP_450_before_err = np.array(',list(FCFP_450_before_SMUtrouble1*FCFP_err_450_before_SMUtrouble1),')')
    print('FCFP_450_after      = np.array(',list(FCFP_450_after_SMUtrouble1),')')
    print('FCFP_450_after_err  = np.array(',list(FCFP_450_after_SMUtrouble1*FCFP_err_450_after_SMUtrouble1),')')
    print('FCFP_850_before     = np.array(',list(FCFP_850_before_newfilters1),')')
    print('FCFP_850_before_err = np.array(',list(FCFP_850_before_newfilters1*FCFP_err_850_before_newfilters1),')')
    print('FCFP_850_after      = np.array(',list(FCFP_850_after_newfilters1),')')
    print('FCFP_850_after_err  = np.array(',list(FCFP_850_after_newfilters1*FCFP_err_850_after_newfilters1),')')
    print('\n\n')

    if constrained == 'no':
        print('\n\nInformation To Plot -- SET 2, No ',badsourcelist[1],' -- NOT CONSTRAINED\n\n')
    if constrained == 'yes':
        print('\n\nInformation To Plot -- SET 2, No ',badsourcelist[1],' -- CONSTRAINED\n\n')
    print('x_before = ',list(TautimesAirmass_before-split_up_markers))
    print('x_after = ',list(TautimesAirmass_after-split_up_markers))
    print('\n\nFCF Arcsec\n\n')
    print('FCFA_450_before     = np.array(',list(FCFA_450_before_SMUtrouble2),')')
    print('FCFA_450_before_err = np.array(',list(FCFA_450_before_SMUtrouble2*FCFA_err_450_before_SMUtrouble2),')')
    print('FCFA_450_after      = np.array(',list(FCFA_450_after_SMUtrouble2),')')
    print('FCFA_450_after_err  = np.array(',list(FCFA_450_after_SMUtrouble2*FCFA_err_450_after_SMUtrouble2),')')
    print('FCFA_850_before     = np.array(',list(FCFA_850_before_newfilters2),')')
    print('FCFA_850_before_err = np.array(',list(FCFA_850_before_newfilters2*FCFA_err_850_before_newfilters2),')')
    print('FCFA_850_after      = np.array(',list(FCFA_850_after_newfilters2),')')
    print('FCFA_850_after_err  = np.array(',list(FCFA_850_after_newfilters2*FCFA_err_850_after_newfilters2),')')
    print('\n\nFCF Peak\n\n')
    print('FCFP_450_before     = np.array(',list(FCFP_450_before_SMUtrouble2),')')
    print('FCFP_450_before_err = np.array(',list(FCFP_450_before_SMUtrouble2*FCFP_err_450_before_SMUtrouble2),')')
    print('FCFP_450_after      = np.array(',list(FCFP_450_after_SMUtrouble2),')')
    print('FCFP_450_after_err  = np.array(',list(FCFP_450_after_SMUtrouble2*FCFP_err_450_after_SMUtrouble2),')')
    print('FCFP_850_before     = np.array(',list(FCFP_850_before_newfilters2),')')
    print('FCFP_850_before_err = np.array(',list(FCFP_850_before_newfilters2*FCFP_err_850_before_newfilters2),')')
    print('FCFP_850_after      = np.array(',list(FCFP_850_after_newfilters2),')')
    print('FCFP_850_after_err  = np.array(',list(FCFP_850_after_newfilters2*FCFP_err_850_after_newfilters2),')')
    print('\n\n')

