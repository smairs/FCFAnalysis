import matplotlib.pyplot as plt
import numpy as np

def plot_by_tauxairmass_2lists(FCFP_850_before_newfilters1,FCFP_err_850_before_newfilters1,FCFP_850_after_newfilters1,FCFP_err_850_after_newfilters1,FCFA_850_before_newfilters1,FCFA_err_850_before_newfilters1,FCFA_850_after_newfilters1,FCFA_err_850_after_newfilters1,FCFP_450_before_SMUtrouble1,FCFP_err_450_before_SMUtrouble1,FCFP_450_after_SMUtrouble1,FCFP_err_450_after_SMUtrouble1,FCFA_450_before_SMUtrouble1,FCFA_err_450_before_SMUtrouble1,FCFA_450_after_SMUtrouble1,FCFA_err_450_after_SMUtrouble1,FCFP_850_before_newfilters2,FCFP_err_850_before_newfilters2,FCFP_850_after_newfilters2,FCFP_err_850_after_newfilters2,FCFA_850_before_newfilters2,FCFA_err_850_before_newfilters2,FCFA_850_after_newfilters2,FCFA_err_850_after_newfilters2,FCFP_450_before_SMUtrouble2,FCFP_err_450_before_SMUtrouble2,FCFP_450_after_SMUtrouble2,FCFP_err_450_after_SMUtrouble2,FCFA_450_before_SMUtrouble2,FCFA_err_450_before_SMUtrouble2,FCFA_450_after_SMUtrouble2,FCFA_err_450_after_SMUtrouble2,TauTimesAirmassBins,constrained='no'):

    plt.clf()
    
    FCFP_850_nominal     = 537.0
    FCFP_850_nominal_err = 26.0/537.0
    
    FCFP_450_nominal     = 491.0
    FCFP_450_nominal_err = 67.0/491.0
    
    FCFA_850_nominal     = 2.34
    FCFA_850_nominal_err = 0.08/2.34
    
    FCFA_450_nominal     = 4.71
    FCFA_450_nominal_err = 0.5/4.71
    
    #FCFP_850_before_newfilters     = np.array([534.33,552.64,559.49,564.00,568.33])
    #FCFP_err_850_before_newfilters = np.array([0.0712,0.0754,0.0806,0.0707,0.0829])
    #FCFP_850_after_newfilters      = np.array([526.17,530.82,552.14,552.98,559.89])
    #FCFP_err_850_after_newfilters  = np.array([0.0452,0.0598,0.0680,0.0623,0.0914])
    
    #FCFA_850_before_newfilters     = np.array([2.2600,2.3200,2.3200,2.3600,2.3700])
    #FCFA_err_850_before_newfilters = np.array([0.0416,0.0433,0.0413,0.0443,0.0503])
    #FCFA_850_after_newfilters      = np.array([2.1200,2.1400,2.1900,2.2100,2.1400])
    #FCFA_err_850_after_newfilters  = np.array([0.0425,0.0394,0.0423,0.0337,0.0647])
    
    #FCFP_450_before_SMUtrouble     = np.array([565.15,589.67,530.61,519.35,np.nan])
    #FCFP_err_450_before_SMUtrouble = np.array([0.1921,0.2025,0.2337,0.1678,np.nan])
    #FCFP_450_after_SMUtrouble      = np.array([535.22,534.61,478.20,449.30,np.nan])
    #FCFP_err_450_after_SMUtrouble  = np.array([0.1335,0.1558,0.1967,0.1514,np.nan])
    
    #FCFA_450_before_SMUtrouble     = np.array([4.5800,4.7500,4.5700,4.4900,np.nan])
    #FCFA_err_450_before_SMUtrouble = np.array([0.1090,0.1419,0.2034,0.1884,np.nan])
    #FCFA_450_after_SMUtrouble      = np.array([3.9400,3.7500,3.7800,3.3800,np.nan])
    #FCFA_err_450_after_SMUtrouble  = np.array([0.0943,0.1269,0.1877,0.1444,np.nan])
    
    
    #TauTimesAirmassBins    = [0.03,0.084,0.138,0.192,0.246,0.32]
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
 
    plt.errorbar(TautimesAirmass_before-split_up_markers,FCFP_850_before_newfilters1,yerr=FCFP_850_before_newfilters1*FCFP_err_850_before_newfilters1,linestyle='none',marker='o',color='gray',label='Before New Filters (3 Cals)',markersize=10) # 3 Cals stands for "with CRL618" 
    plt.errorbar(TautimesAirmass_after-split_up_markers,FCFP_850_after_newfilters1,yerr=FCFP_850_after_newfilters1*FCFP_err_850_after_newfilters1,linestyle='none',marker='^',color='k',label='After New Filters (3 Cals)',markersize=10)
    plt.errorbar(TautimesAirmass_before+split_up_markers,FCFP_850_before_newfilters2,yerr=FCFP_850_before_newfilters2*FCFP_err_850_before_newfilters2,linestyle='none',marker='s',label='Before New Filters (2 Cals)',color='gray',markersize=10)
    plt.errorbar(TautimesAirmass_after+split_up_markers,FCFP_850_after_newfilters2,yerr=FCFP_850_after_newfilters2*FCFP_err_850_after_newfilters2,linestyle='none',marker='x',label='After New Filters (2 Cals)',color='k',markersize=10)
    plt.xticks(TauTimesAirmassBins)
    plt.legend(loc='upper left')
    plt.ylabel('Peak FCF, 850 microns (Jy/beam/pW)')
    plt.xlabel('Tau225 x Airmass')
    if constrained == 'no':
        plt.savefig('FCFP_by_TauxAM_850.png',format='png',dpi=300)
    if constrained == 'yes':
        plt.savefig('FCFP_by_TauxAM_850_constrained.png',format='png',dpi=300)
    plt.clf()
    
    
    # FCF A 850 #
    
    TauTimesAirmassDummy = 0
    for i in range(len(TauTimesAirmassBins)-1):
        plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours[TauTimesAirmassDummy])
        TauTimesAirmassDummy = TauTimesAirmassDummy+1
    
    plt.errorbar(TautimesAirmass_before-split_up_markers,FCFA_850_before_newfilters1,yerr=FCFA_850_before_newfilters1*FCFA_err_850_before_newfilters1,linestyle='none',marker='o',color='gray',label='Before New Filters (3 Cals)',markersize=10)
    plt.errorbar(TautimesAirmass_after-split_up_markers,FCFA_850_after_newfilters1,yerr=FCFA_850_after_newfilters1*FCFA_err_850_after_newfilters1,linestyle='none',marker='^',color='k',label='After New Filters (3 Cals)',markersize=10)
    plt.errorbar(TautimesAirmass_before+split_up_markers,FCFA_850_before_newfilters2,yerr=FCFA_850_before_newfilters2*FCFA_err_850_before_newfilters2,linestyle='none',marker='s',label='Before New Filters (2 Cals)',color='gray',markersize=10)
    plt.errorbar(TautimesAirmass_after+split_up_markers,FCFA_850_after_newfilters2,yerr=FCFA_850_after_newfilters2*FCFA_err_850_after_newfilters2,linestyle='none',marker='x',label='After New Filters (2 Cals)',color='k',markersize=10)
    plt.xticks(TauTimesAirmassBins)
    plt.legend(loc='upper left')
    plt.ylabel('Arcsec FCF, 850 microns (Jy/arcsec^2/pW)')
    plt.xlabel('Tau225 x Airmass')
    if constrained == 'no':
        plt.savefig('FCFA_by_TauxAM_850.png',format='png',dpi=300)
    if constrained == 'yes':
        plt.savefig('FCFA_by_TauxAM_850_constrained.png',format='png',dpi=300)
    plt.clf()
    
    # FCF P 450 #
    
    TauTimesAirmassDummy = 0
    for i in range(len(TauTimesAirmassBins)-1):
        plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours[TauTimesAirmassDummy])
        TauTimesAirmassDummy = TauTimesAirmassDummy+1
    
    plt.errorbar(TautimesAirmass_before-split_up_markers,FCFP_450_before_SMUtrouble1,yerr=FCFP_450_before_SMUtrouble1*FCFP_err_450_before_SMUtrouble1,linestyle='none',marker='o',color='gray',label='Before SMU Fix (3 Cals)',markersize=10)
    plt.errorbar(TautimesAirmass_after-split_up_markers,FCFP_450_after_SMUtrouble1,yerr=FCFP_450_after_SMUtrouble1*FCFP_err_450_after_SMUtrouble1,linestyle='none',marker='^',color='k',label='After SMU Fix (3 Cals)',markersize=10)
    plt.errorbar(TautimesAirmass_before+split_up_markers,FCFP_450_before_SMUtrouble2,yerr=FCFP_450_before_SMUtrouble2*FCFP_err_450_before_SMUtrouble2,linestyle='none',marker='s',label='Before SMU Fix (2 Cals)',color='gray',markersize=10)
    plt.errorbar(TautimesAirmass_after+split_up_markers,FCFP_450_after_SMUtrouble2,yerr=FCFP_450_after_SMUtrouble2*FCFP_err_450_after_SMUtrouble2,linestyle='none',marker='x',label='After SMU Fix (2 Cals)',color='k',markersize=10)
    plt.xticks(TauTimesAirmassBins)
    plt.legend(loc='upper left')
    plt.ylabel('Peak FCF, 450 microns (Jy/beam/pW)')
    plt.xlabel('Tau225 x Airmass')
    if constrained == 'no':
        plt.savefig('FCFP_by_TauxAM_450.png',format='png',dpi=300)
    if constrained == 'yes':
        plt.savefig('FCFP_by_TauxAM_450_constrained.png',format='png',dpi=300)
    plt.clf()
    
    # FCF A 450 #
    
    TauTimesAirmassDummy = 0
    for i in range(len(TauTimesAirmassBins)-1):
        plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours[TauTimesAirmassDummy])
        TauTimesAirmassDummy = TauTimesAirmassDummy+1
    
    plt.errorbar(TautimesAirmass_before-split_up_markers,FCFA_450_before_SMUtrouble1,yerr=FCFA_450_before_SMUtrouble1*FCFA_err_450_before_SMUtrouble1,linestyle='none',marker='o',color='gray',label='Before SMU Fix (3 Cals)',markersize=10)
    plt.errorbar(TautimesAirmass_after-split_up_markers,FCFA_450_after_SMUtrouble1,yerr=FCFA_450_after_SMUtrouble1*FCFA_err_450_after_SMUtrouble1,linestyle='none',marker='^',color='k',label='After SMU Fix (3 Cals)',markersize=10)
    plt.errorbar(TautimesAirmass_before+split_up_markers,FCFA_450_before_SMUtrouble2,yerr=FCFA_450_before_SMUtrouble2*FCFA_err_450_before_SMUtrouble2,linestyle='none',marker='s',label='Before SMU Fix (2 Cals)',color='gray',markersize=10)
    plt.errorbar(TautimesAirmass_after+split_up_markers,FCFA_450_after_SMUtrouble2,yerr=FCFA_450_after_SMUtrouble2*FCFA_err_450_after_SMUtrouble2,linestyle='none',marker='x',label='After SMU Fix (2 Cals)',color='k',markersize=10)
    plt.xticks(TauTimesAirmassBins)
    plt.legend(loc='upper left')
    plt.ylabel('Arcsec FCF, 450 microns (Jy/arcsec^2/pW)')
    plt.xlabel('Tau225 x Airmass')
    if constrained == 'no':
        plt.savefig('FCFA_by_TauxAM_450.png',format='png',dpi=300)
    if constrained == 'yes':
        plt.savefig('FCFA_by_TauxAM_450_constrained.png',format='png',dpi=300)
    plt.clf()

    if constrained == 'no':
        print('\n\nInformation To Plot -- NOT CONSTRAINED\n\n')
    if constrained == 'yes':
        print('\n\nInformation To Plot -- CONSTRAINED\n\n')
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
        print('\n\nInformation To Plot NO CRL618!!! -- NOT CONSTRAINED\n\n')
    if constrained == 'yes':
        print('\n\nInformation To Plot NO CRL618!!! -- CONSTRAINED\n\n')
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

