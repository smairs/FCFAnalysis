import numpy as np
import matplotlib.pyplot as plt
import pickle

def summarise_epochs(infofile,Tau225TimesAirmass_bins,badsourcelist):
    
    #infofile = 'FCF_AR_PER_EPOCH_450_EE.bin'
    #infofile = 'FCF_AR_PER_EPOCH_850_EE.bin'
    
    #infofile = 'FCF_AR_PER_EPOCH_450_SS.bin'
    #infofile = 'FCF_AR_PER_EPOCH_850_SS.bin'
    
    #infofile = 'FCF_AR_PER_EPOCH_450_EO.bin'
    #infofile = 'FCF_AR_PER_EPOCH_850_EO.bin'
    
    #infofile = 'FCF_AR_PER_EPOCH_450.bin'
    #infofile = 'FCF_AR_PER_EPOCH_850.bin'
    
    #badsourcelist = ['Arp220','NEPTUNE','MARS', 'CRL618']
    
    #Tau225TimesAirmass_bins = np.array([0.03,0.084])
    #Tau225TimesAirmass_bins = np.array([0.03 , 0.084, 0.138, 0.192, 0.246, 0.32  ])
    
    
    normalise  = False
    use_median = True
    show_details = False
    
    info = pickle.load(open(infofile,'rb'))
    
    last_part_of_file =  infofile.split('_')[-1].split('.bin')[0]
    if last_part_of_file[1]=='5':
        wave = int(infofile.split('_')[-1].split('.bin')[0])
    elif infofile.split('_')[-2][1]=='5':
        wave = int(infofile.split('_')[-2])
    else:
        wave = int(infofile.split('_')[-3]) 
    
    if wave== 850:
    
        nominal_FCF_peak      = 537.0
        nominal_FCF_peak_err  = 26.0/537.0 
        nominal_FCF_arcs      = 2.34
        nominal_FCF_arcs_err  = 0.08/2.34
        nominal_beam_2com     = 14.6
        nominal_beam_eff      = 14.1
    
    if wave == 450:
        nominal_FCF_peak      = 491.0
        nominal_FCF_peak_err  = 67.0/491.0
        nominal_FCF_arcs      = 4.71
        nominal_FCF_arcs_err  = 0.5/4.71
        nominal_beam_2com     = 9.8
        nominal_beam_eff      = 9.6
    # eff = sqrt([FCF_p/FCF_a]/1.133)
    
    markers                     = ['x','o','s','d','+','^']
    colours                     = ['#885FCD','#003593','#20b2aa','#800000','#daa520','#f6546a']
    labels                      = ['Arp220','CRL 2688','CRL 618','Mars','Neptune','Uranus']
    shortepochs                 = ['Pub','Sil1','No1','Sil2','No2','Bla','New F','Mem Off','Mem On','SMU Issue','Temp Fix','Perm Fix','Current']
    
    epoch_ind = 0
    
    x_axis        = []
    
    P_FCF_fig,P_FCF_ax       = plt.subplots()
    A_FCF_fig,A_FCF_ax       = plt.subplots()
    AR_fig,AR_ax             = plt.subplots()
    FWHM_eff_fig,FWHM_eff_ax = plt.subplots()
    
    
    weigh_avg_before_interest_FCFP = []
    weigh_avg_after_interest_FCFP  = []
    weigh_avg_before_interest_FCFA = []
    weigh_avg_after_interest_FCFA  = []
    weigh_std_before_interest_FCFP = []
    weigh_std_after_interest_FCFP  = []
    weigh_std_before_interest_FCFA = []
    weigh_std_after_interest_FCFA  = []
    
    
    if wave == 850:
        interest = 'New Filters'
    else:
        interest = 'SMU Gain Fix'
    
    epochlist = np.array(['Published','Silver WVM 1','No WVM 1','Silver WVM 2','No WVM 2','Black WVM','New Filters','Membrane Off','Membrane On','SMU Trouble','SMU Gain Fix','SMU HW Fix','Helco'])
    if wave == 850:
        badepochlist = ['No WVM 1','No WVM 2','Membrane Off','SMU Trouble']
    if wave == 450:
        badepochlist = ['No WVM 1','No WVM 2','Membrane Off','Membrane On','SMU Trouble']
    
    
    interesting_ind = np.where(epochlist==interest)[0][0]
    
    epoch_dummy = -1
    dummy1 = 0
    for eachepoch in epochlist:
        if show_details:
            print('\n############\n',eachepoch,'\n############\n')
        epoch_dummy = epoch_dummy+1
        Tau225TimesAirmassDummy = 0
    
        source_numbers_div_err2_FCFp_final     = []
        source_numbers_div_err2_FCFa_final     = []
        source_numbers_div_err2_Beam_eff_final = []
        source_numbers_div_err2_AR_final       = []
        
        FCFps_final          = []
        FCFas_final          = []
        ARs_final            = []
        Beam_effs_final      = []
        source_numbers_final = []
        source_numbers_indiv = {}
        FCFps_err_final      = []
        FCFas_err_final      = []
        ARs_err_final        = []
        Beam_effs_err_final  = []
    
        for eachTautimesAMbin in range(len(Tau225TimesAirmass_bins)-1):
            source_numbers_indiv[eachTautimesAMbin] = {}
    
            epoch_ind               = 0.5+dummy1
            epoch_ind_step          = 0.166666
            epoch_ind_step_Tau225AM = Tau225TimesAirmassDummy*0.166666/(len(Tau225TimesAirmass_bins)-1)
            dummy=0
        
            source_numbers_div_err2_FCFp         = []     
            source_numbers_div_err2_FCFa         = [] 
            source_numbers_div_err2_Beam_eff     = [] 
            source_numbers_div_err2_AR           = [] 
            FCFps                                = [] 
            FCFas                                = [] 
            ARs                                  = [] 
            Beam_effs                            = [] 
            source_numbers                       = [] 
            FCFps_err                            = [] 
            FCFas_err                            = [] 
            ARs_err                              = [] 
            Beam_effs_err                        = [] 
        
            for eachsource in sorted(list(info[eachepoch].keys())):
                epoch_ind_step_this = (epoch_ind_step*dummy)+epoch_ind_step_Tau225AM
                TauTimesAM_ind = np.where(np.logical_and(info[eachepoch][eachsource]['TAU225']*info[eachepoch][eachsource]['AM']>=Tau225TimesAirmass_bins[eachTautimesAMbin],info[eachepoch][eachsource]['TAU225']*info[eachepoch][eachsource]['AM']<=Tau225TimesAirmass_bins[eachTautimesAMbin+1]))
                source_numbers_indiv[eachTautimesAMbin][eachsource] = len(info[eachepoch][eachsource]['PEAK_FLUXES'][TauTimesAM_ind])
                if eachsource not in badsourcelist:
                    if show_details:
                        print(eachsource,' INCLUDED IN WEIGHTED AVERAGE')
                    #if not np.isnan(info[eachepoch][eachsource]['mean_Peak_FCF']):
                    if len(info[eachepoch][eachsource]['PEAK_FLUXES'][TauTimesAM_ind])>0:
                        source_numbers_final.append(len(info[eachepoch][eachsource]['PEAK_FLUXES'][TauTimesAM_ind]))
    
                        # This section is all just for the weighted average:
    
                        source_numbers.append(len(info[eachepoch][eachsource]['PEAK_FLUXES'][TauTimesAM_ind]))
                        source_numbers_div_err2_FCFp.append(len(info[eachepoch][eachsource]['PEAK_FCF'][TauTimesAM_ind])/(np.std(info[eachepoch][eachsource]['PEAK_FCF'][TauTimesAM_ind],ddof=1)**2))
                        source_numbers_div_err2_FCFa.append(len(info[eachepoch][eachsource]['ARCSEC_FCF'][TauTimesAM_ind])/(np.std(info[eachepoch][eachsource]['ARCSEC_FCF'][TauTimesAM_ind],ddof=1)**2))
                        source_numbers_div_err2_AR.append(len(info[eachepoch][eachsource]['ARs'][TauTimesAM_ind])/(np.std(info[eachepoch][eachsource]['ARs'][TauTimesAM_ind],ddof=1))**2)
                        source_numbers_div_err2_Beam_eff.append(len(info[eachepoch][eachsource]['Beam_Eff'][TauTimesAM_ind])/(np.std(info[eachepoch][eachsource]['Beam_Eff'][TauTimesAM_ind],ddof=1))**2)
                        if use_median:
                            FCFps.append(np.median(info[eachepoch][eachsource]['PEAK_FCF'][TauTimesAM_ind]))
                            FCFas.append(np.median(info[eachepoch][eachsource]['ARCSEC_FCF'][TauTimesAM_ind]))
                            ARs.append(np.median(info[eachepoch][eachsource]['ARs'][TauTimesAM_ind]))
                            Beam_effs.append(np.median(info[eachepoch][eachsource]['Beam_Eff'][TauTimesAM_ind])) 
                        else:
                            FCFps.append(np.mean(info[eachepoch][eachsource]['PEAK_FCF'][TauTimesAM_ind]))
                            FCFas.append(np.mean(info[eachepoch][eachsource]['ARCSEC_FCF'][TauTimesAM_ind]))
                            ARs.append(np.mean(info[eachepoch][eachsource]['ARs'][TauTimesAM_ind]))
                            Beam_effs.append(np.mean(info[eachepoch][eachsource]['Beam_Eff'][TauTimesAM_ind]))
                        FCFps_err.append(np.std(info[eachepoch][eachsource]['PEAK_FCF'][TauTimesAM_ind],ddof=1))
                        FCFas_err.append(np.std(info[eachepoch][eachsource]['ARCSEC_FCF'][TauTimesAM_ind],ddof=1))
                        ARs_err.append(np.std(info[eachepoch][eachsource]['ARs'][TauTimesAM_ind],ddof=1))
                        Beam_effs_err.append(np.std(info[eachepoch][eachsource]['Beam_Eff'][TauTimesAM_ind],ddof=1))
    
                        # End section for weighted average
    
                if normalise == False:
                    if eachepoch=='Published' and eachTautimesAMbin == 0:
                        x_axis.append(epoch_ind+epoch_ind_step_this)
                        P_FCF_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['PEAK_FCF'][TauTimesAM_ind])],yerr=[np.std(info[eachepoch][eachsource]['PEAK_FCF'][TauTimesAM_ind],ddof=1)],color=colours[dummy],label=labels[dummy],marker=markers[dummy])
                        A_FCF_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['ARCSEC_FCF'][TauTimesAM_ind])],yerr=[np.std(info[eachepoch][eachsource]['ARCSEC_FCF'][TauTimesAM_ind],ddof=1)],color=colours[dummy],label=labels[dummy],marker=markers[dummy])
                        AR_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['ARs'][TauTimesAM_ind])],yerr=[np.std(info[eachepoch][eachsource]['ARs'][TauTimesAM_ind],ddof=1)],color=colours[dummy],label=labels[dummy],marker=markers[dummy])
                        FWHM_eff_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['Beam_Eff'][TauTimesAM_ind])],yerr=[np.std(info[eachepoch][eachsource]['Beam_Eff'][TauTimesAM_ind],ddof=1)],color=colours[dummy],label=labels[dummy],marker=markers[dummy])
                    else:
                        P_FCF_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['PEAK_FCF'][TauTimesAM_ind])],yerr=[np.std(info[eachepoch][eachsource]['PEAK_FCF'][TauTimesAM_ind],ddof=1)],color=colours[dummy],marker=markers[dummy])
                        A_FCF_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['ARCSEC_FCF'][TauTimesAM_ind])],yerr=[np.std(info[eachepoch][eachsource]['ARCSEC_FCF'][TauTimesAM_ind],ddof=1)],color=colours[dummy],marker=markers[dummy])
                        AR_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['ARs'][TauTimesAM_ind])],yerr=[np.std(info[eachepoch][eachsource]['ARs'][TauTimesAM_ind],ddof=1)],color=colours[dummy],marker=markers[dummy])
                        FWHM_eff_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['Beam_Eff'][TauTimesAM_ind])],yerr=[np.std(info[eachepoch][eachsource]['Beam_Eff'][TauTimesAM_ind],ddof=1)],color=colours[dummy],marker=markers[dummy])
                else:
                    if eachepoch=='Published' and eachTautimesAMbin == 0:
                        x_axis.append(epoch_ind+epoch_ind_step_this)
                        P_FCF_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['PEAK_FCF'][TauTimesAM_ind])/nominal_FCF_peak],yerr=[np.std(info[eachepoch][eachsource]['PEAK_FCF'][TauTimesAM_ind],ddof=1)/nominal_FCF_peak],color=colours[dummy],label=labels[dummy],marker=markers[dummy])
                        A_FCF_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['ARCSEC_FCF'][TauTimesAM_ind])/nominal_FCF_arcs],yerr=[np.std(info[eachepoch][eachsource]['ARCSEC_FCF'][TauTimesAM_ind],ddof=1)/nominal_FCF_arcs],color=colours[dummy],label=labels[dummy],marker=markers[dummy])
                        AR_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['ARs'][TauTimesAM_ind])],yerr=[np.std(info[eachepoch][eachsource]['ARs'][TauTimesAM_ind],ddof=1)],color=colours[dummy],label=labels[dummy],marker=markers[dummy])
                        FWHM_eff_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['Beam_Eff'][TauTimesAM_ind])/nominal_beam_eff],yerr=[np.std(info[eachepoch][eachsource]['Beam_Eff'][TauTimesAM_ind],ddof=1)/nominal_beam_eff],color=colours[dummy],label=labels[dummy],marker=markers[dummy])
                    else:
                        P_FCF_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['PEAK_FCF'][TauTimesAM_ind])/nominal_FCF_peak],yerr=[np.std(info[eachepoch][eachsource]['PEAK_FCF'][TauTimesAM_ind],ddof=1)/nominal_FCF_peak],color=colours[dummy],marker=markers[dummy])
                        A_FCF_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['ARCSEC_FCF'][TauTimesAM_ind])/nominal_FCF_arcs],yerr=[np.std(info[eachepoch][eachsource]['ARCSEC_FCF'][TauTimesAM_ind],ddof=1)/nominal_FCF_arcs],color=colours[dummy],marker=markers[dummy])
                        AR_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['ARs'][TauTimesAM_ind])],yerr=[np.std(info[eachepoch][eachsource]['ARs'][TauTimesAM_ind],ddof=1)],color=colours[dummy],marker=markers[dummy])
                        FWHM_eff_ax.errorbar([epoch_ind+epoch_ind_step_this],[np.mean(info[eachepoch][eachsource]['Beam_Eff'][TauTimesAM_ind])/nominal_beam_eff],yerr=[np.std(info[eachepoch][eachsource]['Beam_Eff'][TauTimesAM_ind],ddof=1)/nominal_beam_eff],color=colours[dummy],marker=markers[dummy])
                
                dummy=dummy+1
        
            norm_weight_FCFp     = np.array(source_numbers_div_err2_FCFp)/sum(source_numbers)
            norm_weight_FCFa     = np.array(source_numbers_div_err2_FCFa)/sum(source_numbers)
            norm_weight_AR       = np.array(source_numbers_div_err2_AR)/sum(source_numbers)
            norm_weight_Beam_eff = np.array(source_numbers_div_err2_Beam_eff)/sum(source_numbers)
    
            weigh_avg_Peak_FCFs_num         = sum(norm_weight_FCFp*np.array(FCFps)) 
            weigh_avg_Arcsec_FCFs_num       = sum(norm_weight_FCFa*np.array(FCFas)) 
            weigh_avg_ARs_num               = sum(norm_weight_AR*np.array(ARs))
            weigh_avg_Beam_eff_num          = sum(norm_weight_Beam_eff*np.array(Beam_effs))
        
            weigh_avg_Peak_FCFs_denom       = sum(norm_weight_FCFp)
            weigh_avg_Arcsec_FCFs_denom     = sum(norm_weight_FCFa)
            weigh_avg_ARs_denom             = sum(norm_weight_AR)
            weigh_avg_Beam_eff_denom        = sum(norm_weight_Beam_eff)
        
            weigh_avg_Peak_FCF_err_num      = sum(np.array(source_numbers)*np.array(FCFps_err))
            weigh_avg_Arcsec_FCF_err_num    = sum(np.array(source_numbers)*np.array(FCFas_err))
            weigh_avg_AR_err_num            = sum(np.array(source_numbers)*np.array(ARs_err))
            weigh_avg_Beam_eff_err_num      = sum(np.array(source_numbers)*np.array(Beam_effs_err))
        
            weigh_avg_err_denom             = sum(source_numbers)
        
        
            simple_average_Peak_FCFs_num    = np.mean(FCFps)
            simple_average_Arcsec_FCFs_num  = np.mean(FCFas)
            simple_average_ARs_num          = np.mean(ARs)
            simple_average_Beam_eff_num     = np.mean(Beam_effs)
        
            simple_average_Peak_FCF_err     = np.std(FCFps,ddof=1)
            simple_average_Arcsec_FCF_err   = np.std(FCFas,ddof=1)
            simple_average_AR_err           = np.std(ARs,ddof=1)
            simple_average_Beam_eff_err     = np.std(Beam_effs,ddof=1)
        
        
            if len(FCFps)>0:
                if normalise == False:
                    if not np.isnan(weigh_avg_Peak_FCFs_num/weigh_avg_Peak_FCFs_denom):
                        P_FCF_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],weigh_avg_Peak_FCFs_num/weigh_avg_Peak_FCFs_denom,yerr=weigh_avg_Peak_FCF_err_num/weigh_avg_err_denom,color='k',marker='o',ms=10)
                        A_FCF_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],weigh_avg_Arcsec_FCFs_num/weigh_avg_Arcsec_FCFs_denom,yerr=weigh_avg_Arcsec_FCF_err_num/weigh_avg_err_denom,color='k',marker='o',ms=10)
                        AR_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],weigh_avg_ARs_num/weigh_avg_ARs_denom,yerr=weigh_avg_AR_err_num/weigh_avg_err_denom,color='k',marker='o',ms=10)
                        FWHM_eff_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],weigh_avg_Beam_eff_num/weigh_avg_Beam_eff_denom,yerr=weigh_avg_Beam_eff_err_num/weigh_avg_err_denom,color='k',marker='o',ms=10) 
                    else:
                        P_FCF_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],simple_average_Peak_FCFs_num,yerr=simple_average_Peak_FCF_err,markeredgecolor='k',markerfacecolor='none',ecolor='k',marker='o',ms=10)
                        A_FCF_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],simple_average_Arcsec_FCFs_num,yerr=simple_average_Arcsec_FCF_err,markeredgecolor='k',markerfacecolor='none',ecolor='k',marker='o',ms=10)
                        AR_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],simple_average_ARs_num,yerr=simple_average_AR_err,markeredgecolor='k',markerfacecolor='none',ecolor='k',marker='o',ms=10)
                        FWHM_eff_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],simple_average_Beam_eff_num,yerr=simple_average_Beam_eff_err,markeredgecolor='k',markerfacecolor='none',ecolor='k',marker='o',ms=10)
                else:
                    if not np.isnan(weigh_avg_Peak_FCFs_num/weigh_avg_Peak_FCFs_denom):
                        P_FCF_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],weigh_avg_Peak_FCFs_num/weigh_avg_Peak_FCFs_denom/nominal_FCF_peak,yerr=weigh_avg_Peak_FCF_err_num/weigh_avg_err_denom/nominal_FCF_peak,color='k',marker='o',ms=10)
                        A_FCF_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],weigh_avg_Arcsec_FCFs_num/weigh_avg_Arcsec_FCFs_denom/nominal_FCF_arcs,yerr=weigh_avg_Arcsec_FCF_err_num/weigh_avg_err_denom/nominal_FCF_arcs,color='k',marker='o',ms=10)
                        AR_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],weigh_avg_ARs_num/weigh_avg_ARs_denom,yerr=weigh_avg_AR_err_num/weigh_avg_err_denom,color='k',marker='o',ms=10)
                        FWHM_eff_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],weigh_avg_Beam_eff_num/weigh_avg_Beam_eff_denom/nominal_beam_eff,yerr=weigh_avg_Beam_eff_err_num/weigh_avg_err_denom/nominal_beam_eff,color='k',marker='o',ms=10)
                    else:
                        P_FCF_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],simple_average_Peak_FCFs_num/nominal_FCF_peak,yerr=simple_average_Peak_FCF_err/nominal_FCF_peak,markeredgecolor='k',markerfacecolor='none',ecolor='k',marker='o',ms=10)
                        A_FCF_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],simple_average_Arcsec_FCFs_num/nominal_FCF_arcs,yerr=simple_average_Arcsec_FCF_err/nominal_FCF_arcs,markeredgecolor='k',markerfacecolor='none',ecolor='k',marker='o',ms=10)
                        AR_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],simple_average_ARs_num,yerr=simple_average_AR_err,markeredgecolor='k',markerfacecolor='none',ecolor='k',marker='o',ms=10)
                        FWHM_eff_ax.errorbar([epoch_ind+(epoch_ind_step+epoch_ind_step_Tau225AM)*2.5],simple_average_Beam_eff_num/nominal_beam_eff,yerr=simple_average_Beam_eff_err/nominal_beam_eff,markeredgecolor='k',markerfacecolor='none',ecolor='k',marker='o',ms=10)
    
                if show_details:
                    print('\nPeak FCF (weighted Avg): ',round(weigh_avg_Peak_FCFs_num/weigh_avg_Peak_FCFs_denom,2),'+/-',round(weigh_avg_Peak_FCF_err_num/weigh_avg_err_denom,2))
                    print('Peak FCF (simple Avg)  : ',round(simple_average_Peak_FCFs_num,2),'+/-',round(simple_average_Peak_FCF_err,2))
                    print('Arcs FCF (weighted Avg): ',round(weigh_avg_Arcsec_FCFs_num/weigh_avg_Arcsec_FCFs_denom,2),'+/-',round(weigh_avg_Arcsec_FCF_err_num/weigh_avg_err_denom,2))
                    print('Arcs FCF (simple Avg)  : ',round(simple_average_Arcsec_FCFs_num,2),'+/-',round(simple_average_Arcsec_FCF_err,2))
                    print('Beam eff (weighted Avg): ',round(weigh_avg_Beam_eff_num/weigh_avg_Beam_eff_denom,2),'+/-',round(weigh_avg_Beam_eff_err_num/weigh_avg_err_denom,2))
                    print('Beam eff (simple Avg)  : ',round(simple_average_Beam_eff_num,2),'+/-',round(simple_average_Beam_eff_err,2),'\n')
               
                if epoch_dummy < interesting_ind:
                    if eachepoch not in badepochlist:
                         weigh_avg_before_interest_FCFP.append(weigh_avg_Peak_FCFs_num/weigh_avg_Peak_FCFs_denom)
                         weigh_std_before_interest_FCFP.append(weigh_avg_Peak_FCF_err_num/weigh_avg_err_denom)
                         weigh_avg_before_interest_FCFA.append(weigh_avg_Arcsec_FCFs_num/weigh_avg_Arcsec_FCFs_denom)
                         weigh_std_before_interest_FCFA.append(weigh_avg_Arcsec_FCF_err_num/weigh_avg_err_denom)
                else:
                    if eachepoch not in badepochlist:
                         if wave==850 and eachepoch=='Membrane On':
                             weigh_avg_after_interest_FCFP.append(np.nan)
                             weigh_std_after_interest_FCFP.append(np.nan)
                             weigh_avg_after_interest_FCFA.append(np.nan)
                             weigh_std_after_interest_FCFA.append(np.nan)
                         else:
                             weigh_avg_after_interest_FCFP.append(weigh_avg_Peak_FCFs_num/weigh_avg_Peak_FCFs_denom)
                             weigh_std_after_interest_FCFP.append(weigh_avg_Peak_FCF_err_num/weigh_avg_err_denom)
                             weigh_avg_after_interest_FCFA.append(weigh_avg_Arcsec_FCFs_num/weigh_avg_Arcsec_FCFs_denom)
                             weigh_std_after_interest_FCFA.append(weigh_avg_Arcsec_FCF_err_num/weigh_avg_err_denom)
                     
    
                source_numbers_div_err2_FCFp_final.append(sum(source_numbers)/(weigh_avg_Peak_FCF_err_num/weigh_avg_err_denom)**2.0)
                source_numbers_div_err2_FCFa_final.append(sum(source_numbers)/(weigh_avg_Arcsec_FCF_err_num/weigh_avg_err_denom)**2.0)
                source_numbers_div_err2_Beam_eff_final.append(sum(source_numbers)/(weigh_avg_Beam_eff_err_num/weigh_avg_err_denom)**2.0)
                source_numbers_div_err2_AR_final.append(sum(source_numbers)/(weigh_avg_AR_err_num/weigh_avg_err_denom)**2.0)
    
                FCFps_final.append(weigh_avg_Peak_FCFs_num/weigh_avg_Peak_FCFs_denom)
                FCFas_final.append(weigh_avg_Arcsec_FCFs_num/weigh_avg_Arcsec_FCFs_denom)
                Beam_effs_final.append(weigh_avg_Beam_eff_num/weigh_avg_Beam_eff_denom)
                ARs_final.append(weigh_avg_ARs_num/weigh_avg_ARs_denom)
    
                FCFps_err_final.append(weigh_avg_Peak_FCF_err_num/weigh_avg_err_denom)
                FCFas_err_final.append(weigh_avg_Arcsec_FCF_err_num/weigh_avg_err_denom)
                Beam_effs_err_final.append(weigh_avg_Beam_eff_err_num/weigh_avg_err_denom)
                ARs_err_final.append(weigh_avg_AR_err_num/weigh_avg_err_denom)
    
            else:
                if show_details:
                    print('\n Number of data points < 3 \n')
        
            P_FCF_ax.axvline(x=epoch_ind+epoch_ind_step*(dummy)-0.08,color='k',linestyle='dotted',linewidth=1)
            A_FCF_ax.axvline(x=epoch_ind+epoch_ind_step*(dummy)-0.08,color='k',linestyle='dotted',linewidth=1)
            AR_ax.axvline(x=epoch_ind+epoch_ind_step*(dummy)-0.08,color='k',linestyle='dotted',linewidth=1)
            FWHM_eff_ax.axvline(x=epoch_ind+epoch_ind_step*(dummy)-0.08,color='k',linestyle='dotted',linewidth=1)
            if epoch_ind == 5.5:
                if wave == 850:
                    P_FCF_ax.axvline(x=epoch_ind+epoch_ind_step*(dummy)-0.08,color='k',linestyle='solid',linewidth=2)
                    A_FCF_ax.axvline(x=epoch_ind+epoch_ind_step*(dummy)-0.08,color='k',linestyle='solid',linewidth=2)
                    AR_ax.axvline(x=epoch_ind+epoch_ind_step*(dummy)-0.08,color='k',linestyle='solid',linewidth=2)
                    FWHM_eff_ax.axvline(x=epoch_ind+epoch_ind_step*(dummy)-0.08,color='k',linestyle='solid',linewidth=2)
            if epoch_ind == 9.5:
                if wave == 450:
                    P_FCF_ax.axvline(x=epoch_ind+epoch_ind_step*(dummy)-0.08,color='k',linestyle='solid',linewidth=2)
                    A_FCF_ax.axvline(x=epoch_ind+epoch_ind_step*(dummy)-0.08,color='k',linestyle='solid',linewidth=2)
                    AR_ax.axvline(x=epoch_ind+epoch_ind_step*(dummy)-0.08,color='k',linestyle='solid',linewidth=2)
                    FWHM_eff_ax.axvline(x=epoch_ind+epoch_ind_step*(dummy)-0.08,color='k',linestyle='solid',linewidth=2)
    
    
            Tau225TimesAirmassDummy = Tau225TimesAirmassDummy+1
    
        if sum(source_numbers_final)>0: 
            FCFp_golden_star = sum((np.array(source_numbers_div_err2_FCFp_final)/sum(source_numbers_final))*np.array(FCFps_final))/sum(np.array(source_numbers_div_err2_FCFp_final)/sum(source_numbers_final))
            FCFp_golden_star_err = sum((np.array(source_numbers_div_err2_FCFp_final)/sum(source_numbers_final))*np.array(FCFps_err_final))/sum(np.array(source_numbers_div_err2_FCFp_final)/sum(source_numbers_final))
            FCFa_golden_star = sum((np.array(source_numbers_div_err2_FCFa_final)/sum(source_numbers_final))*np.array(FCFas_final))/sum(np.array(source_numbers_div_err2_FCFa_final)/sum(source_numbers_final))
            FCFa_golden_star_err = sum((np.array(source_numbers_div_err2_FCFa_final)/sum(source_numbers_final))*np.array(FCFas_err_final))/sum(np.array(source_numbers_div_err2_FCFa_final)/sum(source_numbers_final))
            Beam_eff_golden_star = sum((np.array(source_numbers_div_err2_Beam_eff_final)/sum(source_numbers_final))*np.array(Beam_effs_final))/sum(np.array(source_numbers_div_err2_Beam_eff_final)/sum(source_numbers_final))
            Beam_eff_golden_star_err = sum((np.array(source_numbers_div_err2_Beam_eff_final)/sum(source_numbers_final))*np.array(Beam_effs_err_final))/sum(np.array(source_numbers_div_err2_Beam_eff_final)/sum(source_numbers_final))
            AR_golden_star = sum((np.array(source_numbers_div_err2_AR_final)/sum(source_numbers_final))*np.array(ARs_final))/sum(np.array(source_numbers_div_err2_AR_final)/sum(source_numbers_final))
            AR_golden_star_err = sum((np.array(source_numbers_div_err2_AR_final)/sum(source_numbers_final))*np.array(ARs_err_final))/sum(np.array(source_numbers_div_err2_AR_final)/sum(source_numbers_final))
    
            P_FCF_ax.errorbar([epoch_ind+(epoch_ind_step)*2.5],[FCFp_golden_star],yerr=[FCFp_golden_star_err],markeredgecolor='gold',markerfacecolor='gold',ecolor='gold',marker='*',ms=8)
            A_FCF_ax.errorbar([epoch_ind+(epoch_ind_step)*2.5],[FCFa_golden_star],yerr=[FCFa_golden_star_err],markeredgecolor='gold',markerfacecolor='gold',ecolor='gold',marker='*',ms=8)
            FWHM_eff_ax.errorbar([epoch_ind+(epoch_ind_step)*2.5],[Beam_eff_golden_star],yerr=[Beam_eff_golden_star_err],markeredgecolor='gold',markerfacecolor='gold',ecolor='gold',marker='*',ms=8)
            AR_ax.errorbar([epoch_ind+(epoch_ind_step)*2.5],[AR_golden_star],yerr=[AR_golden_star_err],markeredgecolor='gold',markerfacecolor='gold',ecolor='gold',marker='*',ms=8)    
    
        dummy1 = dummy1+1
    
    
    P_FCF_ax.set_xticks(np.arange(1,13.1,1))
    P_FCF_ax.set_xticklabels(shortepochs,rotation=30)
    
    if normalise == False:
        P_FCF_ax.set_ylabel('Peak FCF (Jy/beam/pW)')
        P_FCF_ax.axhline(y=nominal_FCF_peak+nominal_FCF_peak_err*nominal_FCF_peak,linestyle='dashed',color='k')
        P_FCF_ax.axhline(y=nominal_FCF_peak-nominal_FCF_peak_err*nominal_FCF_peak,linestyle='dashed',color='k')
    else:
        P_FCF_ax.set_ylabel('Peak FCF WRT Nominal')
        P_FCF_ax.axhline(y=1+nominal_FCF_peak_err,linestyle='dashed',color='k')
        P_FCF_ax.axhline(y=1-nominal_FCF_peak_err,linestyle='dashed',color='k')
    P_FCF_ax.set_xlabel('Epochs')    
    P_FCF_ax.legend(loc='upper right')
    
    A_FCF_ax.set_xticks(np.arange(1,13.1,1))
    A_FCF_ax.set_xticklabels(shortepochs,rotation=30)
    if normalise == False:
        A_FCF_ax.set_ylabel('Arcsec FCF (Jy/arcsec^2/pW)')
        A_FCF_ax.axhline(y=nominal_FCF_arcs+nominal_FCF_arcs_err*nominal_FCF_arcs,linestyle='dashed',color='k')
        A_FCF_ax.axhline(y=nominal_FCF_arcs-nominal_FCF_arcs_err*nominal_FCF_arcs,linestyle='dashed',color='k')
    else:
        A_FCF_ax.set_ylabel('Arcsec FCF WRT Nominal')
        A_FCF_ax.axhline(y=1+nominal_FCF_arcs_err,linestyle='dashed',color='k')
        A_FCF_ax.axhline(y=1-nominal_FCF_arcs_err,linestyle='dashed',color='k')
    A_FCF_ax.set_xlabel('Epochs')
    if wave == 450:
        A_FCF_ax.set_ylim(ymin=2.5)
    A_FCF_ax.legend(loc='upper right')
    
    AR_ax.set_xticks(np.arange(1,13.1,1))
    AR_ax.set_xticklabels(shortepochs,rotation=30)
    AR_ax.set_ylabel('Aspect Ratio')
    AR_ax.set_xlabel('Epochs')
    AR_ax.legend(loc='upper right')
    
    FWHM_eff_ax.set_xticks(np.arange(1,13.1,1))
    FWHM_eff_ax.set_xticklabels(shortepochs,rotation=30)
    if normalise == False:
        FWHM_eff_ax.set_ylabel('FWHM_Eff (")')
        FWHM_eff_ax.axhline(y=nominal_beam_eff,linestyle='dotted',color='k')
        FWHM_eff_ax.axhline(y=nominal_beam_2com,linestyle='dashdot',color='k')
    else:
        FWHM_eff_ax.set_ylabel('FWHM_Eff WRT Nominal')
    FWHM_eff_ax.set_xlabel('Epochs')
    FWHM_eff_ax.legend(loc='upper right')
    
    #plt.show() 
    
    #P_FCF_fig.savefig('Peak_FCFs_epochs_'+str(wave)+'.png',dpi=300)
    #A_FCF_fig.savefig('Arcsec_FCFs_epochs_'+str(wave)+'.png',dpi=300)   
    
    P_FCF_fig.savefig('Peak_FCFs_epochs_'+str(wave)+'.pdf',format='pdf')
    A_FCF_fig.savefig('Arcsec_FCFs_epochs_'+str(wave)+'.pdf',format='pdf')
    
    #AR_fig.savefig('AR_epochs_'+str(wave)+'.png',dpi=300)
    #FWHM_eff_fig.savefig('FWHM_eff_epochs_'+str(wave)+'.png',dpi=300)
    
    weigh_avg_before_interest_FCFP = np.array(weigh_avg_before_interest_FCFP)
    weigh_avg_before_interest_FCFA = np.array(weigh_avg_before_interest_FCFA)
    weigh_avg_after_interest_FCFP = np.array(weigh_avg_after_interest_FCFP)
    weigh_avg_after_interest_FCFA = np.array(weigh_avg_after_interest_FCFA)
    
    weigh_std_before_interest_FCFP = np.array(weigh_std_before_interest_FCFP)
    weigh_std_before_interest_FCFA = np.array(weigh_std_before_interest_FCFA)
    weigh_std_after_interest_FCFP = np.array(weigh_std_after_interest_FCFP)
    weigh_std_after_interest_FCFA = np.array(weigh_std_after_interest_FCFA)
    
    print('\n\n################\n################\nFINAL RESULTS\n################\n################\n')
    print('Wavelength = ',wave)
    print('MINIMUM TAU225 X AM = ',Tau225TimesAirmass_bins[0],'\n')
    print('MAXIMUM TAU225 X AM = ',Tau225TimesAirmass_bins[-1],'\n')
    print('NOMINAL FCFP , FCFA = ',round(nominal_FCF_peak,2),'+/-',round(100*nominal_FCF_peak_err,2),'% ,',round(nominal_FCF_arcs,2),'+/-',round(100*nominal_FCF_arcs_err,2),'%','\n\n')
    print('FCFP BEFORE '+interest+' | \tFCFP AFTER '+interest+' | \tFCFA BEFORE '+interest+' | \tFCFA AFTER '+interest,'\n')
    print(round(np.nanmean(weigh_avg_before_interest_FCFP),2),'+/-',round(100*np.nanmean(weigh_std_before_interest_FCFP/weigh_avg_before_interest_FCFP),2),'% | \t',round(np.nanmean(weigh_avg_after_interest_FCFP),2),'+/-',round(100*np.nanmean(weigh_std_after_interest_FCFP/weigh_avg_after_interest_FCFP),2),'% | \t',round(np.nanmean(weigh_avg_before_interest_FCFA),2),'+/-',round(100*np.nanmean(weigh_std_before_interest_FCFA/weigh_avg_before_interest_FCFA),2),'% | \t',round(np.nanmean(weigh_avg_after_interest_FCFA),2),'+/-',round(100*np.nanmean(weigh_std_after_interest_FCFA/weigh_avg_after_interest_FCFA),2),'%')
    print('\n')
    P_FCF_fig.clf()
    A_FCF_fig.clf()
    AR_fig.clf()
    FWHM_eff_fig.clf()
    plt.clf()

    pickle.dump(source_numbers_indiv,open('source_numbers_'+str(Tau225TimesAirmass_bins[0])+'_'+str(Tau225TimesAirmass_bins[1])+'.bin','wb'))
    return(np.nanmean(weigh_avg_before_interest_FCFP),np.nanmean(weigh_std_before_interest_FCFP/weigh_avg_before_interest_FCFP),np.nanmean(weigh_avg_after_interest_FCFP),np.nanmean(weigh_std_after_interest_FCFP/weigh_avg_after_interest_FCFP),np.nanmean(weigh_avg_before_interest_FCFA),np.nanmean(weigh_std_before_interest_FCFA/weigh_avg_before_interest_FCFA),np.nanmean(weigh_avg_after_interest_FCFA),np.nanmean(weigh_std_after_interest_FCFA/weigh_avg_after_interest_FCFA)) 
