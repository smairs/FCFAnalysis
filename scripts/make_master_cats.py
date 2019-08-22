def total_flux_cor(aperture_diameter,wavelength):
    '''This function performs the 
       total flux correction based on Dempsey 2013's
       Curve of growth (adjusting for aperture diamter if not 1arcmin)
    '''

    import numpy as np

    diams = np.arange(20,105.1,5)
    factors_450 = [0.72,0.79,0.85,0.88,0.91,0.94,0.97,0.98,1.00,1.01,1.02,1.02,1.02,1.03,1.03,1.03,1.03,1.03]
    factors_850 = [0.69,0.79,0.85,0.88,0.91,0.94,0.97,0.99,1.00,1.01,1.02,1.03,1.04,1.05,1.05,1.06,1.07,1.08]

    if wavelength=='450':
        p =np.polyfit(diams,factors_450,10)
    else:
        p =np.polyfit(diams,factors_850,10)

    z = np.poly1d(p)

    return (z(aperture_diameter))

def SecCalLookup(sourcename):
    '''
    A lookup table for commonly used
    Secondary Calbrator fluxes
    '''
    SecCalDict = {}
    SecCalDict['CRL2688'] = {}
    SecCalDict['CRL618']  = {}
    SecCalDict['Arp220']  = {}
 
    SecCalDict['CRL2688']['fbeam_450']  = 24.9
    SecCalDict['CRL2688']['fbeam_850']  = 5.64
    SecCalDict['CRL2688']['ftotal_450'] = 29.1
    SecCalDict['CRL2688']['ftotal_850'] = 6.13

    SecCalDict['CRL618']['fbeam_450']   = 11.5
    SecCalDict['CRL618']['fbeam_850']   = 4.89
    SecCalDict['CRL618']['ftotal_450']  = 12.1
    SecCalDict['CRL618']['ftotal_850']  = 5.00

    SecCalDict['Arp220']['fbeam_450']   = 5.2
    SecCalDict['Arp220']['fbeam_850']   = 0.79
    SecCalDict['Arp220']['ftotal_450']  = 5.4
    SecCalDict['Arp220']['ftotal_850']  = 0.81

    return(SecCalDict[sourcename]['fbeam_850'],SecCalDict[sourcename]['ftotal_850'],SecCalDict[sourcename]['fbeam_450'],SecCalDict[sourcename]['ftotal_450'])


def make_cats(BINFILES=['TauRelPipeline_FullResults_CRL2688_26p0to26p5_0p012to0p012_20162017.bin','TauRelPipeline_FullResults_CRL618_26p0to26p5_0p012to0p012_20162017.bin','TauRelPipeline_FullResults_URANUS_26p0to26p5_0p012to0p012_20162017.bin'],BESTKEYS=['Run_0','Run_0','Run_0'],OBSSTLIM=7,OBSENLIM=16,LOWERTAULIM=0.0,UPPERTAULIM=0.32,wavelength='450',FWHMLIM=11.0,AMLIM=5.0,PIX_SIZE=1.0,present_epoch='20180501',present_epoch_only=0,INCLUDESOURCES=['CRL2688','CRL618','URANUS','MARS','NEPTUNE'],aperture_diam=60.0):

    '''
    IDENTIFY_OUTLIERS.PY                                                                         
                                                                                                 
    This program takes the results of TauRelPipeline                                             
    (the dictionary of all the peak flux, total flux, transmission,                              
    weather, metadata etc for all calibrator observations for all tau relations                 
    attempted in TauRelPipeline) and applies constraints to the data                            
    (time of night, size of the object, opacity, airmass)                                        
    to determine robust FCFs. Multiple dictionaries produced by TauRelPipeline                   
    can be used and each one can assume a different opacity relation.                            
                                                                                                 
                                                                                                 
    The program also displays histograms of the weather/source properties/metadata               
    and colours those which lie outside the standard deviation of an FCF vs MJD                  
    plot for both FCF_beam and FCF_total - this way, I can track which properties                
    outliers all seem to have in common and apply appropriate constraints.                       
                                                                                                 
    To Run:                                                                                             
   
    %python3
    >>>from identify_outliers import *
    >>>find_outliers(args*,kwargs*)                
   
    The BINFILES refer to the pickled files output by TauRelPipeline_NOEXT_faster.py
    The BESTKEYS refer to different opacity relations - the program will print out what opacity
    relation is being assumed for each source.                                                 
                                                                                               
    The rest of the keywords are constraints that I have modified to eliminate scatter. For instance,
    every observation with an airmass higher than 1.9 was an outlier in the data. I limit the        
    observation times occassionally to avoid dish deformation in the early evening and late morning, etc.
    '''

    #### Import Necessary Packages ####
    import pickle
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.time import Time
    import os
    from starlink import fluxes 

    # Any FCFs above these thresholds will be thrown out
    # Testing over the past year has shown 8.0 for Arcsec
    # and 900 for Peak to be good baseline values to use

    FCFarcsecupperlim = 8.0
    FCFbeamupperlim   = 900

    # Define the aperture correction in case you are using
    # an aperture to measure total flux that is not 1 arcmin:
    aperture_cor = total_flux_cor(aperture_diam,wavelength)

    # Define the Calibrator sources to be included in this analysis
    # And define a second array with corresponding indices

    CALIBRATOR_SOURCES   = np.array(['Arp220','CRL2688','CRL618','MARS','NEPTUNE','URANUS'])
    CALIBRATOR_SOURCE_ID = np.array([0,1,2,3,4,5])
 
    ##########################################################################
    ##########################################################################
    # Now we want to initialise empty lists for every measurement and metadata
    # parameter that we are going to record - these ones marked *_all
    # will include data from all calibrator sources 

    # "goodind" is going to be defined to be the calibrator observations that
    # meet a set of constraints (Airmass limit, Observation time limit, etc)
    # "nogoodind" represents the full dataset (minus the outliers)   
    ###

    # Lists of source names to match the data for all the other lists
    # This will act as indices for the other lists, to select only a specific
    # Calibrator when necessary. Also, generate corresponding lists of source IDs
    sourceind_all        = []
    sourceindID_all      = []

    # MJDs corresponding to every observation:
    MJD_dates_all        = []
  
    # Peak flux and total flux measurements
    peakfluxes_all       = []
    peakfluxes_all_norm  = []
    totalfluxes_all      = []
    totalfluxes_all_norm = []

    # Peak and Arcsec FCFs
    peak_fcfs_all        = []
    arcsec_fcfs_all      = []
  
    # Metadata
    airmasses_all        = []
    zentrans_all         = []
    transwitham_all      = []
    ATSTART_all          = []
    ATEND_all            = []
    AMSTART_all          = []
    AMEND_all            = []
    TAUST_all            = []
    TAUEN_all            = []
    TAUST_T_all          = []
    TAUEN_T_all          = []
    FWHM1s_all           = []
    FWHM2s_all           = []
    areas_all            = []
    delta_trans_all      = []
    OBSST_all            = []
    OBSEN_all            = []
    fnames_all           = []
    OBSNUM_all           = []
    UTDATE_all           = []
    AZSTART_all          = []
    AZEND_all            = []
    ELSTART_all          = []
    ELEND_all            = []
    HUMSTART_all         = []
    HUMEND_all           = []
    BPSTART_all          = []
    BPEND_all            = []
    WNDSPDST_all         = []
    WNDSPDEN_all         = []
    WNDDIRST_all         = []
    WNDDIREN_all         = []
    TAU225ST_all         = []
    TAU225EN_all         = []
    TAUDATST_all         = []
    TAUDATEN_all         = []
    SEEINGST_all         = []
    SEEINGEN_all         = []
    SEEDATST_all         = []
    SEEDATEN_all         = []
    FRLEGTST_all         = []
    FRLEGTEN_all         = []
    BKLEGTST_all         = []
    BKLEGTEN_all         = []
    MJD_OBS_all          = []
    MJD_END_all          = []
    ELAPTIME_all         = []
 
    ####################################################################################
    ####################################################################################

    # If you want to produce a catalogue for one epoch only,
    # call the function with the epoch defined and this kwarg==1
    if present_epoch_only == 1:
        # name the output catalogue with the epoch
        BINFILES_temp = []
        for i in BINFILES:
            if i.split('_')[4]==present_epoch:
                BINFILES_temp.append(i)
        BINFILES = BINFILES_temp

 
    ############################################3
    # Now, loop through every pipeline result file you've included
    # and gather all mesurements and metadata from that dictionary

    for eachresultsfile in range(len(BINFILES)):
        results = pickle.load(open(BINFILES[eachresultsfile],'rb'))

        date_to_print = BINFILES[eachresultsfile].split('_')[-2]
        
        # The best key refers to the best opacity relation - the one you want to use
        # to perform the FCF analysis
        BESTKEY = BESTKEYS[eachresultsfile]

        # What source are we dealing with:

        CurrentSource = str(results[BESTKEY]['source'])        

        # Skip if it is a weird source we haven't dealt with yet
        if CurrentSource in INCLUDESOURCES:

            # Print out the information we are gathering just so we can check and make sure it is right:

            print('\n'+date_to_print+': OPACITY RELATION ('+str(results[BESTKEY]['source'])+'): tau_'+wavelength+' = '+str(results[BESTKEY]['coeff1'])+'(tau225 - '+str(results[BESTKEY]['coeff2'])+') x Airmass')

            # When gathering this information, only use the data where a reliable peak flux measurement was made
            peakfluxes1 = np.array(results[BESTKEY]['peak_fluxes'])

            peak_fluxes_not_nan = np.where(~np.isnan(peakfluxes1))
    
            # Gather information into arrays and throw away the bad measurements (peak flux = nan)
            peakfluxes  = np.array(results[BESTKEY]['peak_fluxes'])[peak_fluxes_not_nan]
            airmasses   = (np.array(results[BESTKEY]['AMSTART'])+np.array(results[BESTKEY]['AMEND']))/2.0
            airmasses   = airmasses[peak_fluxes_not_nan]
            zentrans    = np.array(results[BESTKEY]['transmissions'])[peak_fluxes_not_nan]
            # Ensure the total fluxes are corrected if using an aperture to measure them 
            # that does not have a diameter of 1arcmin
            totalfluxes = np.array(results[BESTKEY]['total_fluxes'])[peak_fluxes_not_nan]/aperture_cor
            a           = results[BESTKEY]['coeff1'] 
            b           = results[BESTKEY]['coeff2']  
            trans_witham= zentrans**airmasses                                           
            ATSTART     = np.array(results[BESTKEY]['ATSTART'])[peak_fluxes_not_nan]
            ATEND       = np.array(results[BESTKEY]['ATEND'])[peak_fluxes_not_nan]
            AMSTART     = np.array(results[BESTKEY]['AMSTART'])[peak_fluxes_not_nan]
            AMEND       = np.array(results[BESTKEY]['AMEND'])[peak_fluxes_not_nan]
            TAUST       = np.array(results[BESTKEY]['WVMTAUST'])[peak_fluxes_not_nan]
            TAUEN       = np.array(results[BESTKEY]['WVMTAUEN'])[peak_fluxes_not_nan]
            TAUST_T     = np.array(results[BESTKEY]['WVMTAUST_TIME'])[peak_fluxes_not_nan]
            TAUEN_T     = np.array(results[BESTKEY]['WVMTAUEN_TIME'])[peak_fluxes_not_nan]
            FWHM1s      = np.array(results[BESTKEY]['FWHM1s'])[peak_fluxes_not_nan]
            FWHM2s      = np.array(results[BESTKEY]['FWHM2s'])[peak_fluxes_not_nan]
            areas       = np.array(results[BESTKEY]['areas'])[peak_fluxes_not_nan]
            delta_trans = np.array(results[BESTKEY]['delta_trans'])[peak_fluxes_not_nan]
            fnames      = np.array(results[BESTKEY]['fnames'])[peak_fluxes_not_nan]
            OBSNUM      = np.array(results[BESTKEY]['OBSNUM'])[peak_fluxes_not_nan]
            UTDATE      = np.array(results[BESTKEY]['UTDATE'])[peak_fluxes_not_nan]
            AZSTART     = np.array(results[BESTKEY]['AZSTART'])[peak_fluxes_not_nan]
            AZEND       = np.array(results[BESTKEY]['AZEND'])[peak_fluxes_not_nan]
            ELSTART     = np.array(results[BESTKEY]['ELSTART'])[peak_fluxes_not_nan]
            ELEND       = np.array(results[BESTKEY]['ELEND'])[peak_fluxes_not_nan]
            HUMSTART    = np.array(results[BESTKEY]['HUMSTART'])[peak_fluxes_not_nan]
            HUMEND      = np.array(results[BESTKEY]['HUMEND'])[peak_fluxes_not_nan]
            BPSTART     = np.array(results[BESTKEY]['BPSTART'])[peak_fluxes_not_nan]
            BPEND       = np.array(results[BESTKEY]['BPEND'])[peak_fluxes_not_nan]
            WNDSPDST    = np.array(results[BESTKEY]['WNDSPDST'])[peak_fluxes_not_nan]
            WNDSPDEN    = np.array(results[BESTKEY]['WNDSPDEN'])[peak_fluxes_not_nan]
            WNDDIRST    = np.array(results[BESTKEY]['WNDDIRST'])[peak_fluxes_not_nan]
            WNDDIREN    = np.array(results[BESTKEY]['WNDDIREN'])[peak_fluxes_not_nan]
            TAU225ST    = np.array(results[BESTKEY]['TAU225ST'])[peak_fluxes_not_nan]
            TAU225EN    = np.array(results[BESTKEY]['TAU225EN'])[peak_fluxes_not_nan]
            TAUDATST    = np.array(results[BESTKEY]['TAUDATST'])[peak_fluxes_not_nan]
            TAUDATEN    = np.array(results[BESTKEY]['TAUDATEN'])[peak_fluxes_not_nan]
            SEEINGST    = np.array(results[BESTKEY]['SEEINGST'])[peak_fluxes_not_nan]
            SEEINGEN    = np.array(results[BESTKEY]['SEEINGEN'])[peak_fluxes_not_nan]
            SEEDATST    = np.array(results[BESTKEY]['SEEDATST'])[peak_fluxes_not_nan]
            SEEDATEN    = np.array(results[BESTKEY]['SEEDATEN'])[peak_fluxes_not_nan]
            FRLEGTST    = np.array(results[BESTKEY]['FRLEGTST'])[peak_fluxes_not_nan]
            FRLEGTEN    = np.array(results[BESTKEY]['FRLEGTEN'])[peak_fluxes_not_nan]
            BKLEGTST    = np.array(results[BESTKEY]['BKLEGTST'])[peak_fluxes_not_nan]
            BKLEGTEN    = np.array(results[BESTKEY]['BKLEGTEN'])[peak_fluxes_not_nan]
            MJD_OBS     = np.array(results[BESTKEY]['MJD-OBS'])[peak_fluxes_not_nan]
            MJD_END     = np.array(results[BESTKEY]['MJD-END'])[peak_fluxes_not_nan]
            ELAPTIME    = np.array(results[BESTKEY]['ELAPTIME'])[peak_fluxes_not_nan]

            OBSST = []
            OBSEN = []

            for obstart,obsend in zip(np.array(results[BESTKEY]['OBSSTART'])[peak_fluxes_not_nan],np.array(results[BESTKEY]['OBSEND'])[peak_fluxes_not_nan]):
                # Just get the hour
                OBSST.append(int(obstart.split('T')[1].split(':')[0]))
                OBSEN.append(int(obsend.split('T')[1].split(':')[0]))

            OBSST=np.array(OBSST)
            OBSEN=np.array(OBSEN)

            # Now that the bad measurements are thrown away - include information in the
            # "ALL CALIBRATOR" lists  
            
            # MJD_dates and readable_time were used in deprecated part of the code
            # MJD_dates     = []
            # readable_time = []
            for i in range(len(ATSTART)):
                peakfluxes_all.append(peakfluxes[i])
                peakfluxes_all_norm.append(peakfluxes[i]/np.average(peakfluxes))
                totalfluxes_all.append(totalfluxes[i])
                totalfluxes_all_norm.append(totalfluxes[i]/np.average(totalfluxes))
                zentrans_all.append(zentrans[i])
                transwitham_all.append(trans_witham[i])
                airmasses_all.append(airmasses[i])   
                ATSTART_all.append(ATSTART[i])
                ATEND_all.append(ATEND[i])
                AMSTART_all.append(AMSTART[i])
                AMEND_all.append(AMEND[i])
                TAUST_all.append(TAUST[i])
                TAUEN_all.append(TAUEN[i])
                TAUST_T_all.append(TAUST_T[i])
                TAUEN_T_all.append(TAUEN_T[i])
                FWHM1s_all.append(FWHM1s[i])
                FWHM2s_all.append(FWHM2s[i])
                areas_all.append(areas[i])
                delta_trans_all.append(delta_trans[i])
                OBSST_all.append(OBSST[i])
                OBSEN_all.append(OBSEN[i])
                fnames_all.append(fnames[i])   
                OBSNUM_all.append(OBSNUM[i])      
                UTDATE_all.append(UTDATE[i])     
                AZSTART_all.append(AZSTART[i])      
                AZEND_all.append(AZEND[i])      
                ELSTART_all.append(ELSTART[i])       
                ELEND_all.append(ELEND[i])     
                HUMSTART_all.append(HUMSTART[i])      
                HUMEND_all.append(HUMEND[i])     
                BPSTART_all.append(BPSTART[i])     
                BPEND_all.append(BPEND[i])   
                WNDSPDST_all.append(WNDSPDST[i]) 
                WNDSPDEN_all.append(WNDSPDEN[i]) 
                WNDDIRST_all.append(WNDDIRST[i]) 
                WNDDIREN_all.append(WNDDIREN[i])  
                TAU225ST_all.append(TAU225ST[i])  
                TAU225EN_all.append(TAU225EN[i]) 
                TAUDATST_all.append(TAUDATST[i]) 
                TAUDATEN_all.append(TAUDATEN[i])
                # Weed Out Bad Seeing Data
                if len(str(SEEINGST[i]))>0 and str(SEEINGST[i])[0]=="<":
                    SEEINGST_all.append(np.nan)
                else:
                    SEEINGST_all.append(SEEINGST[i])
                if len(str(SEEINGEN[i])) and str(SEEINGEN[i])[0]=="<":
                    SEEINGEN_all.append(np.nan)
                else:
                    SEEINGEN_all.append(SEEINGEN[i])   
                if len(str(SEEDATST[i])) and str(SEEDATST[i])[0]=="<":
                    SEEDATST_all.append(np.nan)
                else:
                    SEEDATST_all.append(SEEDATST[i])  
                if len(str(SEEDATEN[i])) and str(SEEDATEN[i])[0]=="<":
                    SEEDATEN_all.append(np.nan)
                else:
                    SEEDATEN_all.append(SEEDATEN[i]) 
                FRLEGTST_all.append(FRLEGTST[i])     
                FRLEGTEN_all.append(FRLEGTEN[i])    
                BKLEGTST_all.append(BKLEGTST[i])   
                BKLEGTEN_all.append(BKLEGTEN[i])  
                MJD_OBS_all.append(MJD_OBS[i])   
                MJD_END_all.append(MJD_END[i])        
                ELAPTIME_all.append(ELAPTIME[i])     

                # Now let's get MJDs so we can plot light curves
                date  = fnames[i].split('_')[0].split('s')[1]
                year  = date[0:4]
                month = date[4:6]
                day   = date[6:8]
                datestr = year+'-'+month+'-'+day+'T'+'00:00:00'
                t = Time(datestr, format='isot', scale='utc')
                #MJD_dates.append(t.mjd)
                MJD_dates_all.append(t.mjd)
                #readable_time.append(date)
                sourceind_all.append(CurrentSource)
                sourceindID_all.append(CALIBRATOR_SOURCE_ID[np.where(CALIBRATOR_SOURCES==CurrentSource)][0])

            ###########################################
            ###########################################
            # Now we need to calculate the FCFs for each source
            # To do this, we need to know the actual peak and total fluxes
            # So, we need to call Starlink's Flux for the planets and
            # Refer to Dempsey 2013 for the secondary calibrators

            # Let's start with the planets
            if CurrentSource in ['URANUS','MARS','NEPTUNE']:

                # Include all the data
                peak850_actual   = []
                arcsec850_actual = []
                peak450_actual   = []
                arcsec450_actual = []

                # Peak flux = nan already culled!

                for eachfile in range(len(fnames)):
                    if len(np.array(results[BESTKEY]['OBSSTART'])[peak_fluxes_not_nan][eachfile].split('.'))>1:
                        datestr_long     = np.array(results[BESTKEY]['OBSSTART'])[peak_fluxes_not_nan][eachfile]
                    else:
                        datestr_long     = np.array(results[BESTKEY]['OBSSTART'])[peak_fluxes_not_nan][eachfile]+'.00'
                    peak850_actual.append(fluxes.get_flux(CurrentSource,datestr_long,filter_=850).f_beam)
                    arcsec850_actual.append(fluxes.get_flux(CurrentSource,datestr_long,filter_=850).f_total)
                    peak450_actual.append(fluxes.get_flux(CurrentSource,datestr_long,filter_=450).f_beam)
                    arcsec450_actual.append(fluxes.get_flux(CurrentSource,datestr_long,filter_=450).f_total)

                peak850_actual   = np.array(peak850_actual)
                arcsec850_actual = np.array(arcsec850_actual)
                peak450_actual   = np.array(peak450_actual)
                arcsec450_actual = np.array(arcsec450_actual)

                if wavelength == '850':
                    peak_fcfs   = peak850_actual/peakfluxes
                    arcsec_fcfs = (arcsec850_actual/((totalfluxes)*PIX_SIZE**2.0))
                if wavelength == '450':
                    peak_fcfs   = peak450_actual/peakfluxes
                    arcsec_fcfs = (arcsec450_actual/((totalfluxes)*PIX_SIZE**2.0))


            #################
            # Now do Secondary Calibrators

            if CurrentSource in ['CRL2688','CRL618','Arp220']:
                
                peak850_actual,arcsec850_actual,peak450_actual,arcsec450_actual = SecCalLookup(CurrentSource)

                if wavelength == '850':
                    peak_fcfs   = peak850_actual/peakfluxes
                    arcsec_fcfs = (arcsec850_actual/((totalfluxes)*PIX_SIZE**2.0))
                if wavelength == '450':
                    peak_fcfs   = peak450_actual/peakfluxes
                    arcsec_fcfs = (arcsec450_actual/((totalfluxes)*PIX_SIZE**2.0))


            ###################


            #########################
            #########################
            # Add FCF calculations to the "ALL Calibrators" lists

            for eachpeakfcf in peak_fcfs:
                peak_fcfs_all.append(eachpeakfcf)
            for eacharcsecfcf in arcsec_fcfs:
                arcsec_fcfs_all.append(eacharcsecfcf)

    # Now, we have all the metadata, measurements, peak flux calculations and everything from all epochs
    # With the only constraint being Peak Flux = Nan observations removed
    # Additionally, we have 2 lists: sourceind and sourceindID that match the length and order of all these
    # other lists - so we can pick out information about individual sources using these indexes

    peak_fcfs                  = np.array(peak_fcfs_all)
    arcsec_fcfs                = np.array(arcsec_fcfs_all)
    MJD_dates                  = np.array(MJD_dates_all)
    sourceind                  = np.array(sourceind_all)
    sourceindID                = np.array(sourceindID_all)
    zentrans                   = np.array(zentrans_all)
    trans_witham               = np.array(transwitham_all)
    airmasses                  = np.array(airmasses_all)
    ATSTART                    = np.array(ATSTART_all)    
    ATEND                      = np.array(ATEND_all)      
    AMSTART                    = np.array(AMSTART_all)    
    AMEND                      = np.array(AMEND_all)      
    TAUST                      = np.array(TAUST_all)      
    TAUEN                      = np.array(TAUEN_all)      
    TAUST_T                    = np.array(TAUST_T_all)    
    TAUEN_T                    = np.array(TAUEN_T_all)    
    FWHM1s                     = np.array(FWHM1s_all)     
    FWHM2s                     = np.array(FWHM2s_all)     
    areas                      = np.array(areas_all)      
    delta_trans                = np.array(delta_trans_all)
    OBSST                      = np.array(OBSST_all)      
    OBSEN                      = np.array(OBSEN_all)      
    peakfluxes                 = np.array(peakfluxes_all)
    peakfluxes_norm            = np.array(peakfluxes_all_norm)
    totalfluxes                = np.array(totalfluxes_all) 
    totalfluxes_norm           = np.array(totalfluxes_all_norm)
    fnames                     = np.array(fnames_all)
    OBSNUM                     = np.array(OBSNUM_all)
    UTDATE                     = np.array(UTDATE_all)
    AZSTART                    = np.array(AZSTART_all)
    AZEND                      = np.array(AZEND_all)
    ELSTART                    = np.array(ELSTART_all)
    ELEND                      = np.array(ELEND_all)
    HUMSTART                   = np.array(HUMSTART_all)
    HUMEND                     = np.array(HUMEND_all)
    BPSTART                    = np.array(BPSTART_all)
    BPEND                      = np.array(BPEND_all)
    WNDSPDST                   = np.array(WNDSPDST_all) 
    WNDSPDEN                   = np.array(WNDSPDEN_all) 
    WNDDIRST                   = np.array(WNDDIRST_all) 
    WNDDIREN                   = np.array(WNDDIREN_all) 
    TAU225ST                   = np.array(TAU225ST_all) 
    TAU225EN                   = np.array(TAU225EN_all) 
    TAUDATST                   = np.array(TAUDATST_all) 
    TAUDATEN                   = np.array(TAUDATEN_all) 
    SEEINGST                   = np.array(SEEINGST_all) 
    SEEINGEN                   = np.array(SEEINGEN_all) 
    SEEDATST                   = np.array(SEEDATST_all) 
    SEEDATEN                   = np.array(SEEDATEN_all) 
    FRLEGTST                   = np.array(FRLEGTST_all) 
    FRLEGTEN                   = np.array(FRLEGTEN_all) 
    BKLEGTST                   = np.array(BKLEGTST_all) 
    BKLEGTEN                   = np.array(BKLEGTEN_all) 
    MJD_OBS                    = np.array(MJD_OBS_all) 
    MJD_END                    = np.array(MJD_END_all) 
    ELAPTIME                   = np.array(ELAPTIME_all) 

    # NOW MAKE THE MASTER CATALOGUES WITH ALL THIS INFORMATION
    # In the name, a and b refer to the opacity relation coefficients being used

    if present_epoch_only == 1:
        catalogue_file=open('master_catalogue_'+wavelength+'_'+str(a)+'_'+str(b)+'_'+present_epoch+'_only.txt','w')
    else:
        catalogue_file=open('master_catalogue_'+wavelength+'_'+str(a)+'_'+str(b)+'.txt','w')

        catalogue_file.write('#SOURCE\tSOURCE_ID\tUTDATE\tOBSNUM\tPEAK_FLUX\tNORM_PEAK_FLUX\tTOTAL_FLUX\tNORM_TOTAL_FLUX\tFWHM1\tFWHM2\tAREA\tPEAK_FCF\tARCSEC_FCF\tTRANS\tZEN_TRANS\tMJDST\tMJDEN\tOBSST\tOBSEN\tAMSTART\tAMEND\tATSTART\tATEND\tFRLEGTST\tFRLEGTEN\tBKLEGTST\tBKLEGTEN\tAZSTART\tAZEND\tELSTART\tELEND\tHUMSTART\tHUMEND\tBPSTART\tBPEND\tWNDSPDST\tWNDSPDEN\tWNDDIRST\tWNDDIREN\tWVMTAUST\tWVMTAUEN\tWVMTAUST_TIME\tWVMTAUEN_TIME\tSMATAU225ST\tSMATAU225EN\tSMATAU225ST_TIME\tSMATAU225EN_TIME\tSEEINGST\tSEEINGEN\tSEEDATST\tSEEDATEN\tELAPTIME\tFNAME\tDIFF_TRANSOBS_TRANSMODEL\n')

    for eachfile in range(len(ATSTART)):
   
        catalogue_file.write(str(sourceind[eachfile])+'\t'+str(sourceindID[eachfile])+'\t'+str(UTDATE[eachfile])+'\t'+str(OBSNUM[eachfile])+'\t'+str(peakfluxes[eachfile])+'\t'+str(peakfluxes_norm[eachfile])+'\t'+str(totalfluxes[eachfile])+'\t'+str(totalfluxes_norm[eachfile])+'\t'+str(FWHM1s[eachfile])+'\t'+str(FWHM2s[eachfile])+'\t'+str(areas[eachfile])+'\t'+str(peak_fcfs[eachfile])+'\t'+str(arcsec_fcfs[eachfile])+'\t'+str(trans_witham[eachfile])+'\t'+str(zentrans[eachfile])+'\t'+str(MJD_OBS[eachfile])+'\t'+str(MJD_END[eachfile])+'\t'+str(OBSST[eachfile])+'\t'+str(OBSEN[eachfile])+'\t'+str(AMSTART[eachfile])+'\t'+str(AMEND[eachfile])+'\t'+str(ATSTART[eachfile])+'\t'+str(ATEND[eachfile])+'\t'+str(FRLEGTST[eachfile])+'\t'+str(FRLEGTEN[eachfile])+'\t'+str(BKLEGTST[eachfile])+'\t'+str(BKLEGTEN[eachfile])+'\t'+str(AZSTART[eachfile])+'\t'+str(AZEND[eachfile])+'\t'+str(ELSTART[eachfile])+'\t'+str(ELEND[eachfile])+'\t'+str(HUMSTART[eachfile])+'\t'+str(HUMEND[eachfile])+'\t'+str(BPSTART[eachfile])+'\t'+str(BPEND[eachfile])+'\t'+str(WNDSPDST[eachfile])+'\t'+str(WNDSPDEN[eachfile])+'\t'+str(WNDDIRST[eachfile])+'\t'+str(WNDDIREN[eachfile])+'\t'+str(TAUST[eachfile])+'\t'+str(TAUEN[eachfile])+'\t'+str(TAUST_T[eachfile])+'\t'+str(TAUEN_T[eachfile])+'\t'+str(TAU225ST[eachfile])+'\t'+str(TAU225EN[eachfile])+'\t'+str(TAUDATST[eachfile])+'\t'+str(TAUDATEN[eachfile])+'\t'+str(SEEINGST[eachfile])+'\t'+str(SEEINGEN[eachfile])+'\t'+str(SEEDATST[eachfile])+'\t'+str(SEEDATEN[eachfile])+'\t'+str(ELAPTIME[eachfile])+'\t'+str(fnames[eachfile])+'\t'+str(delta_trans[eachfile])+'\n')

    catalogue_file.close()



    # Now, apply additional constraints based on Airmass, FWHM, Observing time etc.. 

    ##################################################################################
    ############  DEFINE FURTHER CONSTRAINTS ON OBSERVING CONIDITONS #################
    ###########               TO BE USED LATER                       #################
    goodind_all     = np.where(np.logical_and(airmasses<AMLIM,np.logical_and(np.logical_and(FWHM1s<FWHMLIM,FWHM2s<FWHMLIM),np.logical_and(np.logical_and((TAUST+TAUEN)/2.0>LOWERTAULIM,(TAUST+TAUEN)/2.0<UPPERTAULIM),np.logical_and(OBSST>=OBSSTLIM,OBSEN<OBSENLIM)))))
    ##################################################################################

    peak_fcfs_goodind        = np.array(peak_fcfs_all)[goodind_all] 
    arcsec_fcfs_goodind      = np.array(arcsec_fcfs_all)[goodind_all] 
    MJD_dates_goodind        = np.array(MJD_dates_all)[goodind_all] 
    sourceind_goodind        = np.array(sourceind_all)[goodind_all] 
    sourceindID_goodind      = np.array(sourceindID_all)[goodind_all] 
    totalfluxes_goodind      = np.array(totalfluxes_all)[goodind_all]
    totalfluxes_norm_goodind = np.array(totalfluxes_all_norm)[goodind_all] 
    zentrans_goodind         = np.array(zentrans_all)[goodind_all]
    trans_witham_goodind     = np.array(transwitham_all)[goodind_all]            
    airmasses_goodind        = np.array(airmasses_all)[goodind_all]
    ATSTART_goodind          = np.array(ATSTART_all)[goodind_all]
    ATEND_goodind            = np.array(ATEND_all)[goodind_all]    
    AMSTART_goodind          = np.array(AMSTART_all)[goodind_all]   
    AMEND_goodind            = np.array(AMEND_all)[goodind_all]
    TAUST_goodind            = np.array(TAUST_all)[goodind_all]
    TAUEN_goodind            = np.array(TAUEN_all)[goodind_all] 
    TAUST_T_goodind          = np.array(TAUST_T_all)[goodind_all] 
    TAUEN_T_goodind          = np.array(TAUEN_T_all)[goodind_all]
    FWHM1s_goodind           = np.array(FWHM1s_all)[goodind_all] 
    FWHM2s_goodind           = np.array(FWHM2s_all)[goodind_all]
    areas_goodind            = np.array(areas_all)[goodind_all]
    delta_trans_goodind      = np.array(delta_trans_all)[goodind_all]  
    OBSST_goodind            = np.array(OBSST_all)[goodind_all]
    OBSEN_goodind            = np.array(OBSEN_all)[goodind_all]
    peakfluxes_goodind       = np.array(peakfluxes_all)[goodind_all]
    peakfluxes_norm_goodind  = np.array(peakfluxes_all_norm)[goodind_all]
    fnames_goodind           = np.array(fnames_all)[goodind_all] 
    OBSNUM_goodind           = np.array(OBSNUM_all)[goodind_all] 
    UTDATE_goodind           = np.array(UTDATE_all)[goodind_all] 
    AZSTART_goodind          = np.array(AZSTART_all)[goodind_all] 
    AZEND_goodind            = np.array(AZEND_all)[goodind_all] 
    ELSTART_goodind          = np.array(ELSTART_all)[goodind_all] 
    ELEND_goodind            = np.array(ELEND_all)[goodind_all] 
    HUMSTART_goodind         = np.array(HUMSTART_all)[goodind_all] 
    HUMEND_goodind           = np.array(HUMEND_all)[goodind_all] 
    BPSTART_goodind          = np.array(BPSTART_all)[goodind_all] 
    BPEND_goodind            = np.array(BPEND_all)[goodind_all]  
    WNDSPDST_goodind         = np.array(WNDSPDST_all)[goodind_all] 
    WNDSPDEN_goodind         = np.array(WNDSPDEN_all)[goodind_all] 
    WNDDIRST_goodind         = np.array(WNDDIRST_all)[goodind_all] 
    WNDDIREN_goodind         = np.array(WNDDIREN_all)[goodind_all] 
    TAU225ST_goodind         = np.array(TAU225ST_all)[goodind_all] 
    TAU225EN_goodind         = np.array(TAU225EN_all)[goodind_all] 
    TAUDATST_goodind         = np.array(TAUDATST_all)[goodind_all] 
    TAUDATEN_goodind         = np.array(TAUDATEN_all)[goodind_all] 
    SEEINGST_goodind         = np.array(SEEINGST_all)[goodind_all] 
    SEEINGEN_goodind         = np.array(SEEINGEN_all)[goodind_all] 
    SEEDATST_goodind         = np.array(SEEDATST_all)[goodind_all] 
    SEEDATEN_goodind         = np.array(SEEDATEN_all)[goodind_all] 
    FRLEGTST_goodind         = np.array(FRLEGTST_all)[goodind_all] 
    FRLEGTEN_goodind         = np.array(FRLEGTEN_all)[goodind_all] 
    BKLEGTST_goodind         = np.array(BKLEGTST_all)[goodind_all] 
    BKLEGTEN_goodind         = np.array(BKLEGTEN_all)[goodind_all] 
    MJD_OBS_goodind          = np.array(MJD_OBS_all)[goodind_all] 
    MJD_END_goodind          = np.array(MJD_END_all)[goodind_all] 
    ELAPTIME_goodind         = np.array(ELAPTIME_all)[goodind_all] 


    if present_epoch_only == 1:
        catalogue_file_goodind = open('master_catalogue_constrained_'+wavelength+'_'+str(results[BESTKEY]['coeff1'])+'_'+str(results[BESTKEY]['coeff2'])+'_'+present_epoch+'_only.txt','w')
    else:
        catalogue_file_goodind = open('master_catalogue_constrained_'+wavelength+'_'+str(results[BESTKEY]['coeff1'])+'_'+str(results[BESTKEY]['coeff2'])+'.txt','w')

    # Print the constraints at the top of the file so we know what has been done!
    catalogue_file_goodind.write('# Constraints:\n')
    catalogue_file_goodind.write('# (TAUST+TAUEN)/2.0 > '+str(LOWERTAULIM)+'\n')
    catalogue_file_goodind.write('# (TAUST+TAUEN)/2.0 < '+str(UPPERTAULIM)+'\n')
    if OBSSTLIM < 10:
        catalogue_file_goodind.write('# OBSST > UT 0'+str(OBSSTLIM)+':00\n')
    else:
        catalogue_file_goodind.write('# OBSST > UT '+str(OBSSTLIM)+':00\n')
    if OBSENLIM < 10:
        catalogue_file_goodind.write('# OBSEN < UT 0'+str(OBSENLIM)+':00\n')
    else:
        catalogue_file_goodind.write('# OBSEN < UT '+str(OBSENLIM)+':00\n')
    catalogue_file_goodind.write('# (AMSTART+AMSEND)/2.0 < '+str(AMLIM)+'\n')
    catalogue_file_goodind.write('# FWHM1 and FWHM2 < '+str(FWHMLIM)+'\n')
    catalogue_file_goodind.write('#SOURCE\tSOURCE_ID\tUTDATE\tOBSNUM\tPEAK_FLUX\tNORM_PEAK_FLUX\tTOTAL_FLUX\tNORM_TOTAL_FLUX\tFWHM1\tFWHM2\tAREA\tPEAK_FCF\tARCSEC_FCF\tTRANS\tZEN_TRANS\tMJDST\tMJDEN\tOBSST\tOBSEN\tAMSTART\tAMEND\tATSTART\tATEND\tFRLEGTST\tFRLEGTEN\tBKLEGTST\tBKLEGTEN\tAZSTART\tAZEND\tELSTART\tELEND\tHUMSTART\tHUMEND\tBPSTART\tBPEND\tWNDSPDST\tWNDSPDEN\tWNDDIRST\tWNDDIREN\tWVMTAUST\tWVMTAUEN\tWVMTAUST_TIME\tWVMTAUEN_TIME\tSMATAU225ST\tSMATAU225EN\tSMATAU225ST_TIME\tSMATAU225EN_TIME\tSEEINGST\tSEEINGEN\tSEEDATST\tSEEDATEN\tELAPTIME\tFNAME\tDIFF_TRANSOBS_TRANSMODEL\n')

    for eachfile in range(len(ATSTART_goodind)):

        catalogue_file_goodind.write(str(sourceind_goodind[eachfile])+'\t'+str(sourceindID_goodind[eachfile])+'\t'+str(UTDATE_goodind[eachfile])+'\t'+str(OBSNUM_goodind[eachfile])+'\t'+str(peakfluxes_goodind[eachfile])+'\t'+str(peakfluxes_norm_goodind[eachfile])+'\t'+str(totalfluxes[eachfile])+'\t'+str(totalfluxes_norm_goodind[eachfile])+'\t'+str(FWHM1s_goodind[eachfile])+'\t'+str(FWHM2s_goodind[eachfile])+'\t'+str(areas_goodind[eachfile])+'\t'+str(peak_fcfs_goodind[eachfile])+'\t'+str(arcsec_fcfs_goodind[eachfile])+'\t'+str(trans_witham_goodind[eachfile])+'\t'+str(zentrans_goodind[eachfile])+'\t'+str(MJD_OBS_goodind[eachfile])+'\t'+str(MJD_END_goodind[eachfile])+'\t'+str(OBSST_goodind[eachfile])+'\t'+str(OBSEN_goodind[eachfile])+'\t'+str(AMSTART_goodind[eachfile])+'\t'+str(AMEND_goodind[eachfile])+'\t'+str(ATSTART_goodind[eachfile])+'\t'+str(ATEND_goodind[eachfile])+'\t'+str(FRLEGTST_goodind[eachfile])+'\t'+str(FRLEGTEN_goodind[eachfile])+'\t'+str(BKLEGTST_goodind[eachfile])+'\t'+str(BKLEGTEN_goodind[eachfile])+'\t'+str(AZSTART_goodind[eachfile])+'\t'+str(AZEND_goodind[eachfile])+'\t'+str(ELSTART_goodind[eachfile])+'\t'+str(ELEND_goodind[eachfile])+'\t'+str(HUMSTART_goodind[eachfile])+'\t'+str(HUMEND_goodind[eachfile])+'\t'+str(BPSTART_goodind[eachfile])+'\t'+str(BPEND_goodind[eachfile])+'\t'+str(WNDSPDST_goodind[eachfile])+'\t'+str(WNDSPDEN_goodind[eachfile])+'\t'+str(WNDDIRST_goodind[eachfile])+'\t'+str(WNDDIREN_goodind[eachfile])+'\t'+str(TAUST_goodind[eachfile])+'\t'+str(TAUEN_goodind[eachfile])+'\t'+str(TAUST_T_goodind[eachfile])+'\t'+str(TAUEN_T_goodind[eachfile])+'\t'+str(TAU225ST_goodind[eachfile])+'\t'+str(TAU225EN_goodind[eachfile])+'\t'+str(TAUDATST_goodind[eachfile])+'\t'+str(TAUDATEN_goodind[eachfile])+'\t'+str(SEEINGST_goodind[eachfile])+'\t'+str(SEEINGEN_goodind[eachfile])+'\t'+str(SEEDATST_goodind[eachfile])+'\t'+str(SEEDATEN_goodind[eachfile])+'\t'+str(ELAPTIME_goodind[eachfile])+'\t'+str(fnames_goodind[eachfile])+'\t'+str(delta_trans_goodind[eachfile])+'\n')

    catalogue_file_goodind.close()
 
    if not os.path.exists('catalogues'): os.system('mkdir catalogues')
    os.system('mv master*txt catalogues/')

    print('\n\nMaster Catalogues now created and stored in directory "catalogues/"\n\n')

################# MASTER CATALOGUES NOW DONE READ OUT AND SORTED ################################


#####
#
#
# 20190821: The rest of this code is now deprecated and may be deleted entirely in the future
# The reason it is still here is to harvest useful code and to replot trends in
# the metadata along with FCF in case I want to revisit the "good ind" that I apply
# Really, the analysis has moved beyond this now, where I have several other codes
# to plot the useful information from the master catalogues generated above.
#
#
#####




















##################
##################
## Now begin second part of the code which will probably be moved to its own separate code (or several)
##################
##################
#
#    unique_sources = []
#    for source in sourceind:
#        if source not in unique_sources:
#            unique_sources.append(source)
#
#    def get_spaced_colors(n):
#        max_value = 16581375 #256**3
#        interval = int(max_value / n)
#        colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]
#        return [(int(i[:2], 16)/255.0, int(i[2:4], 16)/255.0, int(i[4:], 16)/255.0) for i in colors]
#
#    colours = []
#    for i in range(len(unique_sources)):
#        colours.append(get_spaced_colors(len(unique_sources))[i])
#
#    sd_num=[]
#    for i in peak_fcfs:
#        sd_num.append((i-np.average(peak_fcfs))**2.0)
#    sd = np.sqrt(sum(sd_num)/(len(sd_num)-1))
#
#    FCFoutlierind = np.where(np.logical_or(peak_fcfs>(np.average(peak_fcfs)+sd),peak_fcfs<(np.average(peak_fcfs)-sd)))
#    FCFgoodind    = np.where(np.logical_and(peak_fcfs<(np.average(peak_fcfs)+sd),peak_fcfs>(np.average(peak_fcfs)-sd)))
#
#    for eachsource in range(len(unique_sources)):
#        plt.scatter(MJD_dates[np.where(sourceind==unique_sources[eachsource])],peak_fcfs[np.where(sourceind==unique_sources[eachsource])],color=colours[eachsource],label=unique_sources[eachsource])
#        #plt.scatter(np.array(MJD_dates)[FCFoutlierind],peak_fcfs[FCFoutlierind],color='r',label='FCF_p > 1 SD')
#        if unique_sources[eachsource] == 'URANUS':
#            URANUS_MJDS_FCF_PEAK = MJD_dates[np.where(sourceind==unique_sources[eachsource])]
#            URANUS_FCF_PEAK = peak_fcfs[np.where(sourceind==unique_sources[eachsource])]
#            URANUS_AM_FCF_peak = (AMSTART[np.where(sourceind==unique_sources[eachsource])]+AMEND[np.where(sourceind==unique_sources[eachsource])])/2.0
#        if unique_sources[eachsource] == 'CRL2688':
#            CRL2688_MJDS_FCF_PEAK = MJD_dates[np.where(sourceind==unique_sources[eachsource])]
#            CRL2688_FCF_PEAK = peak_fcfs[np.where(sourceind==unique_sources[eachsource])]
#            CRL2688_AM_FCF_peak = (AMSTART[np.where(sourceind==unique_sources[eachsource])]+AMEND[np.where(sourceind==unique_sources[eachsource])])/2.0
#        if unique_sources[eachsource] == 'CRL618':
#            CRL618_MJDS_FCF_PEAK = MJD_dates[np.where(sourceind==unique_sources[eachsource])]
#            CRL618_FCF_PEAK = peak_fcfs[np.where(sourceind==unique_sources[eachsource])]
#            CRL618_AM_FCF_peak = (AMSTART[np.where(sourceind==unique_sources[eachsource])]+AMEND[np.where(sourceind==unique_sources[eachsource])])/2.0
#        if unique_sources[eachsource] == 'MARS':
#            MARS_MJDS_FCF_PEAK = MJD_dates[np.where(sourceind==unique_sources[eachsource])]
#            MARS_FCF_PEAK = peak_fcfs[np.where(sourceind==unique_sources[eachsource])]
#            MARS_AM_FCF_peak = (AMSTART[np.where(sourceind==unique_sources[eachsource])]+AMEND[np.where(sourceind==unique_sources[eachsource])])/2.0
#        if unique_sources[eachsource] == 'NEPTUNE':
#            NEPTUNE_MJDS_FCF_PEAK = MJD_dates[np.where(sourceind==unique_sources[eachsource])]
#            NEPTUNE_FCF_PEAK = peak_fcfs[np.where(sourceind==unique_sources[eachsource])]
#            NEPTUNE_AM_FCF_peak = (AMSTART[np.where(sourceind==unique_sources[eachsource])]+AMEND[np.where(sourceind==unique_sources[eachsource])])/2.0
#
#    plt.axhline(y=np.average(peak_fcfs)+sd,color='k',linestyle='dashed')
#    plt.axhline(y=np.average(peak_fcfs)-sd,color='k',linestyle='dashed')
#    plt.xlabel('MJD')
#    plt.ylabel('FCF_peak')
#    plt.legend(loc='upper right')
#    #plt.show()
#
#    plt.clf()
#
#    sd_num=[]
#    for i in peakfluxes[goodind]:
#        sd_num.append((i-np.average(peakfluxes[goodind]))**2.0)
#    sdPF = np.sqrt(sum(sd_num)/(len(sd_num)-1))
#
#    for eachsource in range(len(unique_sources)):
#        plt.scatter(trans_witham[goodind][np.where(sourceind==unique_sources[eachsource])],peakfluxes[goodind][np.where(sourceind==unique_sources[eachsource])],color=colours[eachsource],label=unique_sources[eachsource])
#        #plt.scatter(trans_witham[goodind][FCFoutlierind],peakfluxes[goodind][FCFoutlierind],color='r')
#    plt.axhline(y=np.average(peakfluxes[goodind])+sdPF,color='k',linestyle='dashed')
#    plt.axhline(y=np.average(peakfluxes[goodind])-sdPF,color='k',linestyle='dashed')
#    plt.xlabel('exp(-tau*airmass)')
#    plt.ylabel('peak_flux')
#    #plt.show()
#
#
#    plt.clf()
#
## PLOT FCF versus Transmission
#
#    sd_num=[]
#    for i in peak_fcfs:
#        if i<FCFbeamupperlim:
#            sd_num.append((i-np.average(np.array(peak_fcfs)[np.where(np.array(peak_fcfs)<FCFbeamupperlim)]))**2.0)
#    sd = np.sqrt(sum(sd_num)/(len(sd_num)-1))
#
#    for eachsource in range(len(unique_sources)):
#        plt.scatter(trans_witham[goodind][np.where(sourceind==unique_sources[eachsource])],peak_fcfs[np.where(sourceind==unique_sources[eachsource])],color=colours[eachsource],label=unique_sources[eachsource])
#        #plt.scatter(trans_witham[goodind][FCFoutlierind],peakfluxes[goodind][FCFoutlierind],color='r')
#    plt.axhline(y=np.average(peak_fcfs)+sd,color='k',linestyle='dashed')
#    plt.axhline(y=np.average(peak_fcfs)-sd,color='k',linestyle='dashed')
#    plt.xlabel('exp(-tau*airmass)')
#    plt.ylabel('FCF_peak')
#    #plt.show()
#
#
#
#    for eachsource in range(len(unique_sources)):
#        plt.scatter(MJD_dates[np.where(sourceind==unique_sources[eachsource])],peakfluxes[goodind][np.where(sourceind==unique_sources[eachsource])],color=colours[eachsource],label=unique_sources[eachsource])
#    #plt.scatter(np.array(MJD_dates)[FCFoutlierind],peakfluxes[goodind][FCFoutlierind],color='r',label='FCF_p > 1SD from avg')
#    #for i in range(len(np.array(MJD_dates)[np.where(peak_fcfs>600)])):
#    #    plt.text(np.array(MJD_dates)[np.where(peak_fcfs>600)][i],peakfluxes[goodind][np.where(peak_fcfs>600)][i],np.array(readable_time)[np.where(peak_fcfs>600)][i])
#    plt.xlabel('MJD')
#    plt.ylabel('peak_flux')
#    #plt.show()
#
#    AverageFCF_peak = np.average(np.array(peak_fcfs)[np.where(np.array(peak_fcfs)<FCFbeamupperlim)])
#    AverageFCF_peak_within1SD = np.average(peak_fcfs[FCFgoodind])
#
#    #print ('\n\nNumber of data points included in the FCF_peak calculation: '+str(len(np.array(peak_fcfs)[np.where(np.array(peak_fcfs)<FCFbeamupperlim)]))+'\n')
#
#    #print ('Average FCF_peak (all data within constraints)          : {:6.2f} \pm {:3.1f}\n'.format(AverageFCF_peak,sd))
#
######################################
######################################
#
#    arcsec_fcfs = arcsec_fcfs[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]
#
#    sd_num=[]
#    for i in arcsec_fcfs[np.where(arcsec_fcfs<FCFarcsecupperlim)]:
#        sd_num.append((i-np.average(arcsec_fcfs[np.where(arcsec_fcfs<FCFarcsecupperlim)]))**2.0)
#    sd = np.sqrt(sum(sd_num)/(len(sd_num)-1))
#
#
#    FCFoutlierindarcsec = np.where(np.logical_or(arcsec_fcfs>(np.average(arcsec_fcfs)+sd),arcsec_fcfs<(np.average(arcsec_fcfs)-sd)))
#    FCFgoodind    = np.where(np.logical_and(arcsec_fcfs<(np.average(arcsec_fcfs)+sd),arcsec_fcfs>(np.average(arcsec_fcfs)-sd)))
#
#    outlier_in_both = [i for i in FCFoutlierindarcsec[0] if i in FCFoutlierind[0]]
#
#    for eachsource in range(len(unique_sources)):
#        plt.scatter(MJD_dates[np.where(np.isnan(arcsec_fcfs)*-1+1>0)][np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])],arcsec_fcfs[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])],color=colours[eachsource],label=unique_sources[eachsource])
#    #plt.scatter(np.array(MJD_dates)[np.where(np.isnan(arcsec_fcfs)*-1+1>0)][FCFoutlierindarcsec],arcsec_fcfs[FCFoutlierindarcsec],color='darkgoldenrod',label='FCF_a > 1SD')
#    #plt.scatter(np.array(MJD_dates)[np.where(np.isnan(arcsec_fcfs)*-1+1>0)][outlier_in_both],arcsec_fcfs[outlier_in_both],color='magenta', label = 'FCF_p+FCF_a > 1SD')
#
#        if unique_sources[eachsource] == 'URANUS':
#            URANUS_MJDS_FCF_ARCSEC = MJD_dates[np.where(np.isnan(arcsec_fcfs)*-1+1>0)][np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]
#            URANUS_FCF_ARCSEC = arcsec_fcfs[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]
#            URANUS_AM_FCF_arcsec = (AMSTART[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]+AMEND[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])])/2.0
#        if unique_sources[eachsource] == 'CRL2688':
#            CRL2688_MJDS_FCF_ARCSEC = MJD_dates[np.where(np.isnan(arcsec_fcfs)*-1+1>0)][np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]
#            CRL2688_FCF_ARCSEC = arcsec_fcfs[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]
#            CRL2688_AM_FCF_arcsec = (AMSTART[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]+AMEND[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])])/2.0
#        if unique_sources[eachsource] == 'CRL618':
#            CRL618_MJDS_FCF_ARCSEC = MJD_dates[np.where(np.isnan(arcsec_fcfs)*-1+1>0)][np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]
#            CRL618_FCF_ARCSEC = arcsec_fcfs[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]
#            CRL618_AM_FCF_arcsec = (AMSTART[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]+AMEND[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])])/2.0
#        if unique_sources[eachsource] == 'MARS':
#            MARS_MJDS_FCF_ARCSEC = MJD_dates[np.where(np.isnan(arcsec_fcfs)*-1+1>0)][np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]
#            MARS_FCF_ARCSEC = arcsec_fcfs[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]
#            MARS_AM_FCF_arcsec = (AMSTART[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]+AMEND[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])])/2.0
#        if unique_sources[eachsource] == 'NEPTUNE':
#            NEPTUNE_MJDS_FCF_ARCSEC = MJD_dates[np.where(np.isnan(arcsec_fcfs)*-1+1>0)][np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]
#            NEPTUNE_FCF_ARCSEC = arcsec_fcfs[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]
#            NEPTUNE_AM_FCF_arcsec = (AMSTART[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])]+AMEND[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])])/2.0
#
#    plt.axhline(y=np.average(arcsec_fcfs[np.where(arcsec_fcfs<FCFarcsecupperlim)])+sd,color='k',linestyle='dashed')
#    plt.axhline(y=np.average(arcsec_fcfs[np.where(arcsec_fcfs<FCFarcsecupperlim)])-sd,color='k',linestyle='dashed')
#    plt.xlabel('MJD')
#    plt.ylabel('FCF_arcsec')
#    #plt.show()
#
#    plt.clf()
#
#    sd_num=[]
#    for i in totalfluxes:
#        sd_num.append((i-np.average(totalfluxes))**2.0)
#    sdTF = np.sqrt(sum(sd_num)/(len(sd_num)-1))
#
#
#    for eachsource in range(len(unique_sources)):
#        plt.scatter(trans_witham[goodind][np.where(sourceind==unique_sources[eachsource])],totalfluxes[np.where(sourceind==unique_sources[eachsource])],color=colours[eachsource],label=unique_sources[eachsource])
#    #plt.scatter(trans_witham[goodind][FCFoutlierindarcsec],totalfluxes[goodind][FCFoutlierindarcsec],color='darkgoldenrod')
#    #plt.scatter(trans_witham[goodind][outlier_in_both],totalfluxes[goodind][outlier_in_both],color='magenta')
#    plt.axhline(y=np.average(totalfluxes)+sdTF,color='k',linestyle='dashed')
#    plt.axhline(y=np.average(totalfluxes)-sdTF,color='k',linestyle='dashed')
#    plt.xlabel('exp(-tau*airmass)')
#    plt.ylabel('total_flux')
#    #plt.show()
#
#
#    plt.clf()
#
## PLOT THE FCF AGAINST TRANSMISSION
#
#    sd_num=[]
#    for i in arcsec_fcfs[np.where(arcsec_fcfs<FCFarcsecupperlim)]:
#        sd_num.append((i-np.average(arcsec_fcfs[np.where(arcsec_fcfs<FCFarcsecupperlim)]))**2.0)
#    sd = np.sqrt(sum(sd_num)/(len(sd_num)-1))
#
#    for eachsource in range(len(unique_sources)):
#        plt.scatter(trans_witham[goodind][np.where(np.isnan(arcsec_fcfs)*-1+1>0)][np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])],arcsec_fcfs[np.where(sourceind[np.where(np.isnan(arcsec_fcfs)*-1+1>0)]==unique_sources[eachsource])],color=colours[eachsource],label=unique_sources[eachsource])
#    #plt.scatter(np.array(MJD_dates)[np.where(np.isnan(arcsec_fcfs)*-1+1>0)][FCFoutlierindarcsec],arcsec_fcfs[FCFoutlierindarcsec],color='darkgoldenrod',label='FCF_a > 1SD')
#    #plt.scatter(np.array(MJD_dates)[np.where(np.isnan(arcsec_fcfs)*-1+1>0)][outlier_in_both],arcsec_fcfs[outlier_in_both],color='magenta', label = 'FCF_p+FCF_a > 1SD')
#    plt.axhline(y=np.average(arcsec_fcfs[np.where(arcsec_fcfs<FCFarcsecupperlim)])+sd,color='k',linestyle='dashed')
#    plt.axhline(y=np.average(arcsec_fcfs[np.where(arcsec_fcfs<FCFarcsecupperlim)])-sd,color='k',linestyle='dashed')
#    plt.xlabel('Transmission')
#    plt.ylabel('FCF_arcsec')
#    #plt.show()
#
#
#    for eachsource in range(len(unique_sources)):
#        plt.scatter(MJD_dates[np.where(sourceind==unique_sources[eachsource])],totalfluxes[np.where(sourceind==unique_sources[eachsource])],color=colours[eachsource],label=unique_sources[eachsource])
#    #plt.scatter(np.array(MJD_dates)[FCFoutlierindarcsec],totalfluxes[goodind][FCFoutlierindarcsec],color='darkgoldenrod')
#    #plt.scatter(np.array(MJD_dates)[outlier_in_both],totalfluxes[goodind][outlier_in_both],color='magenta')
#    #for i in range(len(np.array(MJD_dates)[np.where(arcsec_fcfs>600)])):
#    #    plt.text(np.array(MJD_dates)[np.where(arcsec_fcfs>600)][i],totalfluxes[goodind][np.where(arcsec_fcfs>600)][i],np.array(readable_time)[np.where(arcsec_fcfs>600)][i])
#    plt.xlabel('MJD')
#    plt.ylabel('total_flux')
#    #plt.show()
#
#    AverageFCF_arcsec = np.average(arcsec_fcfs[np.where(arcsec_fcfs<FCFarcsecupperlim)])
#    AverageFCF_arcsec_within1SD = np.average(arcsec_fcfs[FCFgoodind])
#
#    #print ('\n\nNumber of data points included in the FCF_arcsec calculation (ONLY INCLUDING FCF_ARCSEC<8.0): '+str(len(arcsec_fcfs[np.where(arcsec_fcfs<FCFarcsecupperlim)]))+'\n')
#    #print ('Average FCF_arcsec (all data within constraints)        : {:6.2f} \pm {:5.3f}\n\n'.format(AverageFCF_arcsec,sd))
#
#    plt.clf()
#    #plt.xlabel('sqrt[(FCFbeam / FCFarcsec) / 1.133]')
#    #plt.ylabel('N')
#    #plt.suptitle('Effective Beam Area (arcsec)')
#    #plt.scatter(arcsec_fcfs,peak_fcfs,color='k')
#    #slopeArea,interceptArea = np.polyfit(arcsec_fcfs,peak_fcfs,1)
#    #plt.plot(arcsec_fcfs,arcsec_fcfs*slopeArea+interceptArea,linestyle='dashed',color='g',linewidth=2.0,label='Beam_Eff = '+str(round(np.sqrt(slopeArea/1.133),2))+'"')
#    #print('')
#    #print('Effective Beam Areas:')
#    #print('beam_area_eff = ',list(np.sqrt((peak_fcfs/arcsec_fcfs)/1.133)))
#    #print('')
#    #plt.legend(loc='upper left')
#    #plt.show()
#    #plt.clf()
#
#
#    #PLOT FCF PEAK AGAINST FCF ARCSEC - FIT LINEAR FUNCTION, SLOPE = Effective Area, FWHM_eff = sqrt(Effective Area/1.133).
#    #
#    # 1.133 comes from:
#    # FWHM = 2*sqrt(2ln(2))sigma
#    # Area = 2*pi*sigma^2
#    #
#    # FWHM = 2 * sqrt([Area*2*ln(2)]/[2*pi])
#    # FWHM = 2 * sqrt(Area*ln(2)/pi)
#    # FWHM = sqrt(Area/[pi/4*ln(2)]) = sqrt(Area/1.133) 
#
#    allpeakFCF_allsources_for_fit   = []
#    allarcsecFCF_allsources_for_fit = []
#    linefit_slope_eachsource        = []
#    linefit_weighted_FWHM           = []
#    linefit_intercept_eachsource    = []
#    linefit_FCFarcsec_eachsource    = []
#    linefit_FCFpeak_eachsource      = []
#    linefit_colours                 = []
#    numpoints_eachsource            = []
#    weightingfactors_eachsource     = []
#
#    for eachsource in range(len(unique_sources)):
#        good_arcsec_fcfs = []
#        good_peak_fcfs   = []
#        good_sourceind   = []
#        for eachind in range(len(arcsec_fcfs)):
#            if np.logical_and(arcsec_fcfs[eachind]<FCFarcsecupperlim,peak_fcfs[eachind]<FCFbeamupperlim):
#                good_arcsec_fcfs.append(arcsec_fcfs[eachind])
#                good_peak_fcfs.append(peak_fcfs[eachind])
#                good_sourceind.append(sourceind[eachind])
#        
#        good_arcsec_fcfs = np.array(good_arcsec_fcfs)
#        good_peak_fcfs   = np.array(good_peak_fcfs)
#        good_sourceind   = np.array(good_sourceind)
#
#        plt.scatter(good_arcsec_fcfs[np.where(good_sourceind==unique_sources[eachsource])],good_peak_fcfs[np.where(good_sourceind==unique_sources[eachsource])],color=colours[eachsource],label=unique_sources[eachsource])
#        for eachpoint in good_arcsec_fcfs[np.where(good_sourceind==unique_sources[eachsource])]:
#            allarcsecFCF_allsources_for_fit.append(eachpoint)
#        for eachpoint in good_peak_fcfs[np.where(good_sourceind==unique_sources[eachsource])]:
#            allpeakFCF_allsources_for_fit.append(eachpoint)
#
#        m,b = np.polyfit(good_arcsec_fcfs[np.where(good_sourceind==unique_sources[eachsource])],good_peak_fcfs[np.where(good_sourceind==unique_sources[eachsource])],1)
#        linefit_slope_eachsource.append(m)
#        linefit_intercept_eachsource.append(b)
#        linefit_FCFarcsec_eachsource.append(good_arcsec_fcfs[np.where(good_sourceind==unique_sources[eachsource])])
#        linefit_FCFpeak_eachsource.append(good_peak_fcfs[np.where(good_sourceind==unique_sources[eachsource])])
#        linefit_colours.append(colours[eachsource])
#        linefit_weighted_FWHM.append(len(good_arcsec_fcfs[np.where(good_sourceind==unique_sources[eachsource])])/float(len(good_arcsec_fcfs[np.where(good_sourceind)]))*np.sqrt(m/1.133))
#        weightingfactors_eachsource.append(len(good_arcsec_fcfs[np.where(good_sourceind==unique_sources[eachsource])])/float(len(good_arcsec_fcfs[np.where(good_sourceind)])))
#
#    m_all,b_all=np.polyfit(allarcsecFCF_allsources_for_fit,allpeakFCF_allsources_for_fit,1)
#
#    #counter=0
#    #for eachfit in range(len(linefit_slope_eachsource)):
#    #    plt.plot(linefit_FCFarcsec_eachsource[eachfit],linefit_slope_eachsource[eachfit]*np.array(linefit_FCFarcsec_eachsource[eachfit])+linefit_intercept_eachsource[eachfit],linestyle='dashed',color=linefit_colours[eachfit])
#    #    plt.text(2.0,680-20*counter,'Slope (A_E, '+str(unique_sources[eachfit])+') = '+str(round(linefit_slope_eachsource[eachfit],3))+', FWHM_E = '+str(round(np.sqrt(linefit_slope_eachsource[eachfit]/1.133),3)))
#    #    counter=counter+1
#    #    print (unique_sources[eachfit]+': Slope (A_E, '+str(unique_sources[eachfit])+') = '+str(round(linefit_slope_eachsource[eachfit],3))+', FWHM_E = '+str(round(np.sqrt(linefit_slope_eachsource[eachfit]/1.133),3)))
#
#    #plt.plot(allarcsecFCF_allsources_for_fit,m_all*np.array(allarcsecFCF_allsources_for_fit)+b_all,linestyle='solid',color='g')
#    #plt.text(2.0,700,'Slope (A_E, all) = '+str(round(m_all,3))+', FWHM_E = '+str(round(np.sqrt(m_all/1.133),3)))
#    #plt.text(2.0,720,'Weighted FWHM = '+str(round(sum(np.array(linefit_weighted_FWHM)[np.where(np.isnan(np.array(linefit_weighted_FWHM))*-1+1>0)])/sum(np.array(weightingfactors_eachsource)[np.where(np.isnan(np.array(linefit_weighted_FWHM))*-1+1>0)]),3))+'"')
#    #print ('Weighted FWHM = '+str(round(sum(np.array(linefit_weighted_FWHM)[np.where(np.isnan(np.array(linefit_weighted_FWHM))*-1+1>0)])/sum(np.array(weightingfactors_eachsource)[np.where(np.isnan(np.array(linefit_weighted_FWHM))*-1+1>0)]),3))+'"')
#    #plt.xlabel('FCF_arcsec')
#    #plt.ylabel('FCF_peak')
#    #plt.legend(loc='upper right')
#    #plt.show()
#
###################### NOW FIND OUT WHAT OUTLIERS ALL HAVE IN COMMON! ###############
#
#    # This is what we are working with 
#
#    #peakfluxes = np.array(results[BESTKEY]['peak_fluxes'])
#    #airmasses=(np.array(results[BESTKEY]['AMSTART'])+np.array(results[BESTKEY]['AMEND']))/2.0
#    #airmasses = airmasses[np.where(np.isnan(np.array(peakfluxes))*-1+1>0)]
#    #zentrans =np.array(results[BESTKEY]['transmissions'])[np.where(np.isnan(np.array(peakfluxes))*-1+1>0)]
#    #totalfluxes=np.array(results[BESTKEY]['total_fluxes'])[np.where(np.isnan(np.array(peakfluxes))*-1+1>0)]
#    #a = results[BESTKEY]['coeff1']
#    #b = results[BESTKEY]['coeff2']
#    #trans_witham = zentrans**airmasses
#    #ATSTART     = np.array(results[BESTKEY]['ATSTART'][np.where(np.isnan(np.array(peakfluxes))*-1+1>0)]
#    #ATEND       = np.array(results[BESTKEY]['ATEND'][np.where(np.isnan(np.array(peakfluxes))*-1+1>0)]
#    ###TAUST       = np.array(results[BESTKEY]['WVMTAUST'][np.where(np.isnan(np.array(peakfluxes))*-1+1>0)]
#    #TAUEN       = np.array(results[BESTKEY]['WVMTAUEN'][np.where(np.isnan(np.array(peakfluxes))*-1+1>0)]
#    #TAUST_T     = np.array(results[BESTKEY]['WVMTAUST_TIME'][np.where(np.isnan(np.array(peakfluxes))*-1+1>0)]
#    ##TAUEN_T     = np.array(results[BESTKEY]['WVMTAUEN_TIME'][np.where(np.isnan(np.array(peakfluxes))*-1+1>0)]
#    #FWHM1s      = np.array(results[BESTKEY]['FWHM1s'][np.where(np.isnan(np.array(peakfluxes))*-1+1>0)]
#    #FWHM2s      = np.array(results[BESTKEY]['FWHM2s'][np.where(np.isnan(np.array(peakfluxes))*-1+1>0)]
#    #areas       = np.array(results[BESTKEY]['areas'][np.where(np.isnan(np.array(peakfluxes))*-1+1>0)]
#    #delta_trans = np.array(results[BESTKEY]['delta_trans'][np.where(np.isnan(np.array(peakfluxes))*-1+1>0)]
#    #OBSST=np.array(OBSST)
#    #OBSEN=np.array(OBSEN)
#    
#    plt.clf()
#    n,bins,patches=plt.hist(FWHM1s[goodind],edgecolor='b',linewidth=2.0,facecolor='none',label='Tau, Time, AM, Size = Constr.')
#    plt.hist(FWHM1s[goodind][FCFoutlierind],edgecolor='r',linestyle='dotted',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_p_out')
#    plt.hist(FWHM1s[goodind][FCFoutlierindarcsec],edgecolor='darkgoldenrod',linestyle='dashed',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_a_out')
#    plt.hist(FWHM1s[goodind][outlier_in_both],edgecolor='magenta',linewidth=2.0,facecolor='magenta',bins=bins,label='Constr. + FCF_*_out',alpha=0.4)
#    plt.xlabel('FWHM1')
#    plt.legend(loc='upper right')
#    plt.show()
#    plt.clf()
#
#    n,bins,patches=plt.hist(FWHM2s[goodind],edgecolor='b',linewidth=2.0,facecolor='none',label='Tau, Time, AM, Size = Constr.')
#    plt.hist(FWHM2s[goodind][FCFoutlierind],edgecolor='r',linestyle='dotted',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_p_out')
#    plt.hist(FWHM2s[goodind][FCFoutlierindarcsec],edgecolor='darkgoldenrod',linestyle='dashed',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_a_out')
#    plt.hist(FWHM2s[goodind][outlier_in_both],edgecolor='magenta',linewidth=2.0,facecolor='magenta',bins=bins,label='Constr. + FCF_*_out',alpha=0.4)
#    plt.xlabel('FWHM2')
#    plt.legend(loc='upper right')
#    plt.show()
#    plt.clf()
#
#    n,bins,patches=plt.hist(areas[goodind],edgecolor='b',linewidth=2.0,facecolor='none',label='Tau, Time, AM, Size = Constr.')
#    plt.hist(areas[goodind][FCFoutlierind],edgecolor='r',linestyle='dotted',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_p_out')
#    plt.hist(areas[goodind][FCFoutlierindarcsec],edgecolor='darkgoldenrod',linestyle='dashed',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_a_out')
#    plt.hist(areas[goodind][outlier_in_both],edgecolor='magenta',linewidth=2.0,facecolor='magenta',bins=bins,label='Constr. + FCF_*_out',alpha=0.4)
#    plt.xlabel('area')
#    plt.legend(loc='upper right')
#    plt.show()
#    plt.clf()
#
#    n,bins,patches=plt.hist(ATSTART[goodind],edgecolor='b',linewidth=2.0,facecolor='none',label='Tau, Time, AM, Size = Constr.')
#    plt.hist(ATSTART[goodind][FCFoutlierind],edgecolor='r',linestyle='dotted',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_p_out')
#    plt.hist(ATSTART[goodind][FCFoutlierindarcsec],edgecolor='darkgoldenrod',linestyle='dashed',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_a_out')
#    plt.hist(ATSTART[goodind][outlier_in_both],edgecolor='magenta',linewidth=2.0,facecolor='magenta',bins=bins,label='Constr. + FCF_*_out',alpha=0.4)
#    plt.xlabel('ATSTART')
#    plt.legend(loc='upper right')
#    plt.show()
#    plt.clf()
#
#    n,bins,patches=plt.hist(ATEND[goodind],edgecolor='b',linewidth=2.0,facecolor='none',label='Tau, Time, AM, Size = Constr.')
#    plt.hist(ATEND[goodind][FCFoutlierind],edgecolor='r',linestyle='dotted',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_p_out')
#    plt.hist(ATEND[goodind][FCFoutlierindarcsec],edgecolor='darkgoldenrod',linestyle='dashed',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_a_out')
#    plt.hist(ATEND[goodind][outlier_in_both],edgecolor='magenta',linewidth=2.0,facecolor='magenta',bins=bins,label='Constr. + FCF_*_out',alpha=0.4)
#    plt.xlabel('ATEND')
#    plt.legend(loc='upper right')
#    plt.show()
#    plt.clf()
#
#    n,bins,patches=plt.hist(AMSTART[goodind],edgecolor='b',linewidth=2.0,facecolor='none',label='Tau, Time, AM, Size = Constr.')
#    #print ('\n\nMEDIAN AIRMASS START: '+str(np.median(AMSTART[goodind]))+'\n\n')
#    plt.hist(AMSTART[goodind][FCFoutlierind],edgecolor='r',linestyle='dotted',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_p_out')
#    plt.hist(AMSTART[goodind][FCFoutlierindarcsec],edgecolor='darkgoldenrod',linestyle='dashed',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_a_out')
#    plt.hist(AMSTART[goodind][outlier_in_both],edgecolor='magenta',linewidth=2.0,facecolor='magenta',bins=bins,label='Constr. + FCF_*_out',alpha=0.4)
#    plt.xlabel('AMSTART')
#    plt.legend(loc='upper right')
#    plt.show()
#    plt.clf()
#
#    n,bins,patches=plt.hist(AMEND[goodind],edgecolor='b',linewidth=2.0,facecolor='none',label='Tau, Time, AM, Size = Constr.')
#    #print ('\n\nMEDIAN AIRMASS END: '+str(np.median(AMEND[goodind]))+'\n\n')
#    plt.hist(AMEND[goodind][FCFoutlierind],edgecolor='r',linestyle='dotted',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_p_out')
#    plt.hist(AMEND[goodind][FCFoutlierindarcsec],edgecolor='darkgoldenrod',linestyle='dashed',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_a_out')
#    plt.hist(AMEND[goodind][outlier_in_both],edgecolor='magenta',linewidth=2.0,facecolor='magenta',bins=bins,label='Constr. + FCF_*_out',alpha=0.4)
#    plt.xlabel('AMEND')
#    plt.legend(loc='upper right')
#    plt.show()
#    plt.clf()
#
#    n,bins,patches=plt.hist(TAUST[goodind],edgecolor='b',linewidth=2.0,facecolor='none',label='Tau, Time, AM, Size = Constr.')
#    plt.hist(TAUST[goodind][FCFoutlierind],edgecolor='r',linestyle='dotted',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_p_out')
#    plt.hist(TAUST[goodind][FCFoutlierindarcsec],edgecolor='darkgoldenrod',linestyle='dashed',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_a_out')
#    plt.hist(TAUST[goodind][outlier_in_both],edgecolor='magenta',linewidth=2.0,facecolor='magenta',bins=bins,label='Constr. + FCF_*_out',alpha=0.4)
#    plt.xlabel('TAUST')
#    plt.legend(loc='upper right')
#    plt.show()
#    plt.clf()
#
#    n,bins,patches=plt.hist(TAUEN[goodind],edgecolor='b',linewidth=2.0,facecolor='none',label='Tau, Time, AM, Size = Constr.')
#    plt.hist(TAUEN[goodind][FCFoutlierind],edgecolor='r',linestyle='dotted',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_p_out')
#    plt.hist(TAUEN[goodind][FCFoutlierindarcsec],edgecolor='darkgoldenrod',linestyle='dashed',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_a_out')
#    plt.hist(TAUEN[goodind][outlier_in_both],edgecolor='magenta',linewidth=2.0,facecolor='magenta',bins=bins,label='Constr. + FCF_*_out',alpha=0.4)
#    plt.xlabel('TAUEN')
#    plt.legend(loc='upper right')
#    plt.show()
#    plt.clf()
#
#    n,bins,patches=plt.hist(zentrans[goodind]*airmasses[goodind],edgecolor='b',linewidth=2.0,facecolor='none',label='Tau, Time, AM, Size = Constr.')
#    plt.hist(zentrans[goodind][FCFoutlierind]*airmasses[goodind][FCFoutlierind],edgecolor='r',linestyle='dotted',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_p_out')
#    plt.hist(zentrans[goodind][FCFoutlierindarcsec]*airmasses[goodind][FCFoutlierindarcsec],edgecolor='darkgoldenrod',linestyle='dashed',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_a_out')
#    plt.hist(zentrans[goodind][outlier_in_both]*airmasses[outlier_in_both],edgecolor='magenta',linewidth=2.0,facecolor='magenta',bins=bins,label='Constr. + FCF_*_out',alpha=0.4)
#    plt.xlabel('tau*airmass')
#    plt.legend(loc='upper right')
#    plt.show()
#    plt.clf()
#
#    n,bins,patches=plt.hist(OBSST[goodind],edgecolor='b',linewidth=2.0,facecolor='none',label='Tau, Time, AM, Size = Constr.')
#    plt.hist(OBSST[goodind][FCFoutlierind],edgecolor='r',linestyle='dotted',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_p_out')
#    plt.hist(OBSST[goodind][FCFoutlierindarcsec],edgecolor='darkgoldenrod',linestyle='dashed',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_a_out')
#    plt.hist(OBSST[goodind][outlier_in_both],edgecolor='magenta',linewidth=2.0,facecolor='magenta',bins=bins,label='Constr. + FCF_*_out',alpha=0.4)
#    plt.xlabel('OBSST')
#    plt.legend(loc='upper right')
#    plt.show()
#    plt.clf()
#
#    n,bins,patches=plt.hist(OBSEN[goodind],edgecolor='b',linewidth=2.0,facecolor='none',label='Tau, Time, AM, Size = Constr.')
#    plt.hist(OBSEN[goodind][FCFoutlierind],edgecolor='r',linestyle='dotted',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_p_out')
#    plt.hist(OBSEN[goodind][FCFoutlierindarcsec],edgecolor='darkgoldenrod',linestyle='dashed',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_a_out')
#    plt.hist(OBSEN[goodind][outlier_in_both],edgecolor='magenta',linewidth=2.0,facecolor='magenta',bins=bins,label='Constr. + FCF_*_out',alpha=0.4)
#    plt.xlabel('OBSEN')
#    plt.legend(loc='upper right')
#    plt.show()
#    plt.clf()
#
#    n,bins,patches = plt.hist(delta_trans[goodind],edgecolor='b',linewidth=2.0,facecolor='none',label='Tau, Time, AM, Size = Constr.')
#    plt.hist(delta_trans[goodind][FCFoutlierind],edgecolor='r',linestyle='dotted',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_p_out')
#    plt.hist(delta_trans[goodind][FCFoutlierindarcsec],edgecolor='darkgoldenrod',linestyle='dashed',linewidth=2.0,facecolor='none',bins=bins,label='Constr. + FCF_a_out')
#    plt.hist(delta_trans[goodind][outlier_in_both],edgecolor='magenta',linewidth=2.0,facecolor='magenta',bins=bins,label='Constr. + FCF_*_out',alpha=0.4)
#    plt.xlabel('delta_trans')
#    plt.legend(loc='upper right')
#    plt.show()
#    plt.clf()
#
#    
#    allsources = np.arange(len(unique_sources))
#    
#    #for eachsource in range(len(unique_sources)):
#    #    plt.bar([allsources[eachsource]], [100*len(sourceind[FCFoutlierind][np.where(sourceind[FCFoutlierind]==unique_sources[eachsource])])/float(len(sourceind[FCFoutlierind]))], align='center', edgecolor='r',linestyle='dotted',linewidth=2.0,facecolor='none',label='Constr. + FCF_p_out')
#    #    plt.bar([allsources[eachsource]], [100*len(sourceind[FCFoutlierindarcsec][np.where(sourceind[FCFoutlierindarcsec]==unique_sources[eachsource])])/float(len(sourceind[FCFoutlierindarcsec]))], align='center',edgecolor='darkgoldenrod',linestyle='dashed',linewidth=2.0,facecolor='none',label='Constr. + FCF_a_out')
#    #    plt.bar([allsources[eachsource]], [100*len(sourceind[outlier_in_both][np.where(sourceind[outlier_in_both]==unique_sources[eachsource])])/float(len(sourceind[outlier_in_both]))], align='center',edgecolor='magenta',linewidth=2.0,facecolor='magenta',label='Constr. + FCF_*_out',alpha=0.4)
#    #plt.xticks(allsources, unique_sources)
#    #plt.ylabel('Number')
#    #plt.title('Peak FCF Outliers by Source')
#    #plt.show()
#    #plt.clf()
#
###########################################
###########################################
###########################################
#
#####################################################################################################
############### NOW PLOT THE FCFs AS A FUNCTION OF MJD AND TRANSMISSION #############################
#
#    
#    #if BINFILES[0].split('_')[4][0:4]=='2016':
#    #    startdate          = Time("2016-11-01T00:00:00.00",format='isot',scale='utc')
#    #    enddate            = Time("2018-02-07T00:00:00.00",format='isot',scale='utc')
#    #else:
#    #    startdate          = Time("2011-05-01T00:00:00.00",format='isot',scale='utc')
#    #    enddate            = Time("2012-06-01T00:00:00.00",format='isot',scale='utc')
#    
#    begin_dempsey_regime  = Time("2011-05-01T00:00:00.00",format='isot',scale='utc') # Stated in paper
#    end_dempsey_regime    = Time("2012-06-01T00:00:00.00",format='isot',scale='utc') # stated in paper
#    silverWVM_dead        = Time("2015-01-27T00:00:00.00",format='isot',scale='utc') # from: engarchive.eao.hawaii.edu
#    blackWVM_installed    = Time("2015-04-10T00:00:00.00",format='isot',scale='utc') # See Fault 20150321.006
#    SCUBA2_filters_down   = Time("2016-10-05T00:00:00.00",format='isot',scale='utc') # From Tonight page archive
#    new_SCUBA2_filters    = Time("2016-11-16T00:00:00.00",format='isot',scale='utc') # from /jcmtdata/raw/scuba2 dates
#    membrane_removed      = Time("2017-12-06T00:00:00.00",format='isot',scale='utc') # I was there on the day
#    membrane_replaced     = Time("2018-01-11T00:00:00.00",format='isot',scale='utc') # I was there on the day
#    may2018_shutdown      = Time("2018-05-01T00:00:00.00",format='isot',scale='utc') # I was there on the day
#    post_may2018_shutdown = Time("2018-05-23T00:00:00.00",format='isot',scale='utc') # I was there on the day
#    SMU_gain_fixed        = Time("2018-06-30T08:10:00.00",format='isot',scale='utc') # From obs log - Jim turned gain up for 1 observation after this at 13:07 UT. That data point is no good
#    
#
#    #startdate = Time(str(min(alldates))[0:4]+'-'+str(min(alldates))[4:6]+'-'+str(min(alldates))[6:8]+'T00:00:00.00',format='isot',scale='utc')
#    #enddate = Time(str(max(alldates))[0:4]+'-'+str(max(alldates))[4:6]+'-'+str(max(alldates))[6:8]+'T00:00:00.00',format='isot',scale='utc')
#
#    #startdate_mjd         = startdate.mjd
#    #enddate_mjd           = enddate.mjd
#
#    startdate_mjd = min(MJD_dates_all)
#    enddate_mjd   = max(MJD_dates_all)
#
#    begin_dempsey_regime_mjd  = begin_dempsey_regime.mjd
#    end_dempsey_regime_mjd    = end_dempsey_regime.mjd
#    silverWVM_dead_mjd        = silverWVM_dead.mjd
#    blackWVM_installed_mjd    = blackWVM_installed.mjd
#    SCUBA2_filters_down_mjd   = SCUBA2_filters_down.mjd 
#    new_SCUBA2_filters_mjd    = new_SCUBA2_filters.mjd
#    membrane_removed_mjd      = membrane_removed.mjd
#    membrane_replaced_mjd     = membrane_replaced.mjd
#    may2018_shutdown_mjd      = may2018_shutdown.mjd
#    post_may2018_shutdown_mjd = post_may2018_shutdown.mjd
#    SMU_gain_fixed_mjd        = SMU_gain_fixed.mjd
#
#    epochlist      = [begin_dempsey_regime_mjd,end_dempsey_regime_mjd,silverWVM_dead_mjd,blackWVM_installed_mjd,SCUBA2_filters_down_mjd,new_SCUBA2_filters_mjd,membrane_removed_mjd,membrane_replaced_mjd,may2018_shutdown_mjd,post_may2018_shutdown_mjd,SMU_gain_fixed.mjd,enddate_mjd]
#    epochlabellist = ['Beg. Pub.','End Pub.','Silver WVM Dead','Black WVM','Filter Ch.','New Filters','Mem. Rem.','Mem. Rep.','May Shutdown','Post May Shutdown','SMU Gain Fix','End'] 
#
#    for i in range(len(epochlist)-1):
#        if i == range(len(epochlist)-1)[-1]:
#            if np.logical_and(startdate_mjd>=epochlist[i],startdate_mjd<=epochlist[i+1]):
#                epochs_included_labels  = [epochlabellist[i]]
#                epochs_included         = [epochlist[i]]
#        else:
#            if np.logical_and(startdate_mjd>=epochlist[i],startdate_mjd<epochlist[i+1]):
#                epochs_included_labels  = [epochlabellist[i]]
#                epochs_included         = [epochlist[i]]
#
#    if np.logical_and(end_dempsey_regime_mjd>startdate_mjd,end_dempsey_regime_mjd<enddate_mjd):
#        epochs_included_labels.append('End Pub.')
#        epochs_included.append(end_dempsey_regime_mjd)
#    if np.logical_and(silverWVM_dead_mjd>startdate_mjd,silverWVM_dead_mjd<enddate_mjd):
#        epochs_included_labels.append('Silver WVM Dead')
#        epochs_included.append(silverWVM_dead_mjd)
#    if np.logical_and(blackWVM_installed_mjd>startdate_mjd,blackWVM_installed_mjd<enddate_mjd):
#        epochs_included_labels.append('Black WVM')
#        epochs_included.append(blackWVM_installed_mjd)
#    if np.logical_and(SCUBA2_filters_down_mjd>startdate_mjd,SCUBA2_filters_down_mjd<enddate_mjd):
#        epochs_included_labels.append('Filter Ch.')
#        epochs_included.append(SCUBA2_filters_down_mjd)
#    if np.logical_and(new_SCUBA2_filters_mjd>startdate_mjd,new_SCUBA2_filters_mjd<enddate_mjd):
#        epochs_included_labels.append('New Filters')
#        epochs_included.append(new_SCUBA2_filters_mjd)
#    if np.logical_and(membrane_removed_mjd>startdate_mjd,membrane_removed_mjd<enddate_mjd):
#        epochs_included_labels.append('Mem. Rem.')
#        epochs_included.append(membrane_removed_mjd)
#    if np.logical_and(membrane_replaced_mjd>startdate_mjd,membrane_replaced_mjd<enddate_mjd):
#        epochs_included_labels.append('Mem. Rep.')
#        epochs_included.append(membrane_replaced_mjd)
#    if np.logical_and(may2018_shutdown_mjd>startdate_mjd,may2018_shutdown_mjd<enddate_mjd):
#        epochs_included_labels.append('May SD')
#        epochs_included.append(may2018_shutdown_mjd)
#    if np.logical_and(post_may2018_shutdown_mjd>startdate_mjd,post_may2018_shutdown_mjd<enddate_mjd):
#        epochs_included_labels.append('Post May SD')
#        epochs_included.append(post_may2018_shutdown_mjd)
#    if np.logical_and(SMU_gain_fixed_mjd>startdate_mjd,SMU_gain_fixed_mjd<enddate_mjd):
#        epochs_included_labels.append('SMU Fix')
#        epochs_included.append(SMU_gain_fixed_mjd)
#
#    for i in range(len(epochlist)-1):
#        if np.logical_and(enddate_mjd>epochlist[i],enddate_mjd<=epochlist[i+1]):
#            epochs_included_labels.append(epochlabellist[i+1])
#            epochs_included.append(epochlist[i+1])
#
#    eachsource_colours = ['#3cb44b','#0082c8','#f58231','#911eb4','#800000','#f032e6','#e6194b','#d2f53c','#fabebe','#008080','#e6beff','#aa6e28','#fffac8','#46f0f0','#aaffc3','#808000','#ffd8b1','#000080','#808080','#ffe119','#000000']
#    epoch_marker = ['o','x','^','v','s','+','<','>','8','p','h']
#    
#    if wavelength == '850':
#        FCF_peak_low_cutoff    = 400
#        FCF_peak_high_cutoff   = 700
#        FCF_arcsec_low_cutoff  = 1.8
#        FCF_arcsec_high_cutoff = 2.5
#        minylimpeak            = 300
#        maxylimpeak            = 800
#        minylimarcsec          = 1.5
#        maxylimarcsec          = 3.5
#    else:
#       FCF_peak_low_cutoff    = 200
#       FCF_peak_high_cutoff   = 900
#       FCF_arcsec_low_cutoff  = 2
#       FCF_arcsec_high_cutoff = 7
#       minylimpeak            = 200
#       maxylimpeak            = 1200
#       minylimarcsec          = 2.1
#       maxylimarcsec          = 8.0
#
#    # Just in case we aren't running all sources at once:
#
#    try:
#        URANUS_FCF_PEAK
#    except NameError:
#        URANUS_FCF_PEAK=np.array([])
#        URANUS_MJDS_FCF_PEAK=np.array([])
#        URANUS_FCF_ARCSEC=np.array([])
#        URANUS_MJDS_FCF_ARCSEC=np.array([])
#
#    try:
#        MARS_FCF_PEAK
#    except NameError:
#        MARS_FCF_PEAK=np.array([])
#        MARS_MJDS_FCF_PEAK=np.array([])
#        MARS_FCF_ARCSEC=np.array([])
#        MARS_MJDS_FCF_ARCSEC=np.array([])
#
#    try:
#        NEPTUNE_FCF_PEAK
#    except NameError:
#        NEPTUNE_FCF_PEAK=np.array([])
#        NEPTUNE_MJDS_FCF_PEAK=np.array([])
#        NEPTUNE_FCF_ARCSEC=np.array([])
#        NEPTUNE_MJDS_FCF_ARCSEC=np.array([])
#
#    try:
#        CRL618_FCF_PEAK
#    except NameError:
#        CRL618_FCF_PEAK=np.array([])
#        CRL618_MJDS_FCF_PEAK=np.array([])
#        CRL618_FCF_ARCSEC=np.array([])
#        CRL618_MJDS_FCF_ARCSEC=np.array([])
#
#    try:
#        CRL2688_FCF_PEAK
#    except NameError:
#        CRL2688_FCF_PEAK=np.array([])
#        CRL2688_MJDS_FCF_PEAK=np.array([])
#        CRL2688_FCF_ARCSEC=np.array([])
#        CRL2688_MJDS_FCF_ARCSEC=np.array([])
#
#    try:
#        Arp220_FCF_PEAK
#    except NameError:
#        Arp220_FCF_PEAK=np.array([])
#        Arp220_MJDS_FCF_PEAK=np.array([])
#        Arp220_FCF_ARCSEC=np.array([])
#        Arp220_MJDS_FCF_ARCSEC=np.array([])
#
#    source_information_dict = {}
#    for eachsource in INCLUDESOURCES:
#        source_information_dict[eachsource]={}
#        source_information_dict[eachsource]['Epoch_Labels'] = epochs_included_labels
#    
#    if 'CRL2688' in INCLUDESOURCES:
#        source_information_dict['CRL2688']['FCF_peaks']       = CRL2688_FCF_PEAK
#        source_information_dict['CRL2688']['FCF_peaks_MJD']   = CRL2688_MJDS_FCF_PEAK
#        source_information_dict['CRL2688']['FCF_arcsecs']     = CRL2688_FCF_ARCSEC
#        source_information_dict['CRL2688']['FCF_arcsecs_MJD'] = CRL2688_MJDS_FCF_ARCSEC
#
#    if 'CRL618' in INCLUDESOURCES:
#        source_information_dict['CRL618']['FCF_peaks']       = CRL618_FCF_PEAK
#        source_information_dict['CRL618']['FCF_peaks_MJD']   = CRL618_MJDS_FCF_PEAK
#        source_information_dict['CRL618']['FCF_arcsecs']     = CRL618_FCF_ARCSEC
#        source_information_dict['CRL618']['FCF_arcsecs_MJD'] = CRL618_MJDS_FCF_ARCSEC
#
#    if 'URANUS' in INCLUDESOURCES:
#        source_information_dict['URANUS']['FCF_peaks']       = URANUS_FCF_PEAK
#        source_information_dict['URANUS']['FCF_peaks_MJD']   = URANUS_MJDS_FCF_PEAK
#        source_information_dict['URANUS']['FCF_arcsecs']     = URANUS_FCF_ARCSEC
#        source_information_dict['URANUS']['FCF_arcsecs_MJD'] = URANUS_MJDS_FCF_ARCSEC
#
#    if 'MARS' in INCLUDESOURCES:
#        source_information_dict['MARS']['FCF_peaks']       = MARS_FCF_PEAK
#        source_information_dict['MARS']['FCF_peaks_MJD']   = MARS_MJDS_FCF_PEAK
#        source_information_dict['MARS']['FCF_arcsecs']     = MARS_FCF_ARCSEC
#        source_information_dict['MARS']['FCF_arcsecs_MJD'] = MARS_MJDS_FCF_ARCSEC
#
#    if 'NEPTUNE' in INCLUDESOURCES:
#        source_information_dict['NEPTUNE']['FCF_peaks']       = NEPTUNE_FCF_PEAK
#        source_information_dict['NEPTUNE']['FCF_peaks_MJD']   = NEPTUNE_MJDS_FCF_PEAK
#        source_information_dict['NEPTUNE']['FCF_arcsecs']     = NEPTUNE_FCF_ARCSEC
#        source_information_dict['NEPTUNE']['FCF_arcsecs_MJD'] = NEPTUNE_MJDS_FCF_ARCSEC
#
#    if 'Arp220' in INCLUDESOURCES:
#        source_information_dict['Arp220']['FCF_peaks']       = Arp220_FCF_PEAK
#        source_information_dict['Arp220']['FCF_peaks_MJD']   = Arp220_MJDS_FCF_PEAK
#        source_information_dict['Arp220']['FCF_arcsecs']     = Arp220_FCF_ARCSEC
#        source_information_dict['Arp220']['FCF_arcsecs_MJD'] = Arp220_MJDS_FCF_ARCSEC 
#
#    source_information_dict_by_epoch = {}
#    for eachsource in INCLUDESOURCES:
#        source_information_dict_by_epoch[eachsource]={}
#        # Don't include the last epoch which is "END" - there are no data points after that!
#        for eachepoch in range(len(epochs_included_labels)-1):
#            source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]          = []
#            source_information_dict_by_epoch[eachsource]['FCF_peaks_MJD_'+str(eachepoch)]      = [] 
#            source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]        = []
#            source_information_dict_by_epoch[eachsource]['FCF_arcsecs_MJD_'+str(eachepoch)]    = []
#
#    for eachsource in INCLUDESOURCES:
#        for i in range(len(epochs_included_labels)-1):
#                for j in range(len(source_information_dict[eachsource]['FCF_peaks_MJD'])):
#                    if i == range(len(epochs_included_labels)-1)[-1]:
#                        if np.logical_and(source_information_dict[eachsource]['FCF_peaks_MJD'][j]>=epochs_included[i],source_information_dict[eachsource]['FCF_peaks_MJD'][j]<=epochs_included[i+1]):
#                            source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(i)].append(source_information_dict[eachsource]['FCF_peaks'][j])
#                            source_information_dict_by_epoch[eachsource]['FCF_peaks_MJD_'+str(i)].append(source_information_dict[eachsource]['FCF_peaks_MJD'][j])
#                            source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(i)].append(source_information_dict[eachsource]['FCF_arcsecs'][j])
#                            source_information_dict_by_epoch[eachsource]['FCF_arcsecs_MJD_'+str(i)].append(source_information_dict[eachsource]['FCF_arcsecs_MJD'][j])
#                    else:
#                        if np.logical_and(source_information_dict[eachsource]['FCF_peaks_MJD'][j]>=epochs_included[i],source_information_dict[eachsource]['FCF_peaks_MJD'][j]<epochs_included[i+1]):
#                            source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(i)].append(source_information_dict[eachsource]['FCF_peaks'][j])
#                            source_information_dict_by_epoch[eachsource]['FCF_peaks_MJD_'+str(i)].append(source_information_dict[eachsource]['FCF_peaks_MJD'][j])
#                            source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(i)].append(source_information_dict[eachsource]['FCF_arcsecs'][j])
#                            source_information_dict_by_epoch[eachsource]['FCF_arcsecs_MJD_'+str(i)].append(source_information_dict[eachsource]['FCF_arcsecs_MJD'][j])
#
#    for eachsource in INCLUDESOURCES:
#        for eachlist in source_information_dict_by_epoch[eachsource].keys():
#            source_information_dict_by_epoch[eachsource][eachlist] = np.array(source_information_dict_by_epoch[eachsource][eachlist])
#
# 
#    # FCF_peak plot!
#    
#    plt.clf()
#
#    colourdummy = 0
#    for eachsource in sorted(INCLUDESOURCES):
#        already_labelled = 0
#        each_source_colour = eachsource_colours[colourdummy]
#        colourdummy=colourdummy+1
#        for eachepoch in range(len(epochs_included_labels)-1):
#            # If there are data points in this epoch:
#            if len(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)])>0:
#                if already_labelled==0:
#                    plt.scatter(source_information_dict_by_epoch[eachsource]['FCF_peaks_MJD_'+str(eachepoch)],source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)],color=each_source_colour,marker=epoch_marker[eachepoch],label=eachsource)
#                    already_labelled = 1
#                else:
#                    plt.scatter(source_information_dict_by_epoch[eachsource]['FCF_peaks_MJD_'+str(eachepoch)],source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)],color=each_source_colour,marker=epoch_marker[eachepoch])
#
#    for eachepoch in range(len(epochs_included)):
#        plt.axvline(x=epochs_included[eachepoch],linestyle='dashed',color='k',linewidth=2)
#        plt.text(epochs_included[eachepoch],maxylimpeak-(0.04*(maxylimpeak-minylimpeak))*(eachepoch+1),epochs_included_labels[eachepoch],clip_on=True)
#
#    plt.legend(loc='upper left')
#    plt.xlim(xmin=startdate_mjd-20,xmax=enddate_mjd+20)
#    plt.ylim(ymin=minylimpeak,ymax=maxylimpeak)
#    plt.xlabel('MJD',fontsize=15)
#    plt.ylabel('FCF_Peak (Jy/pW)',fontsize=15)
#    plt.show()
#    plt.clf()
#    
#
#    ######## PRINT FCF INFORMATION ############
#
#
#    for eachsource in INCLUDESOURCES:
#        print('\n\n#############\n\n'+eachsource+':\n\n')
#        for eachepoch in range(len(epochs_included)-1): # Avoid the epoch: "End"
#            if len(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))])==1:
#                print(epochs_included_labels[eachepoch]+' -> '+epochs_included_labels[eachepoch+1]+', FCF Peak: ',np.average(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))]),'+/- Only 1 data point! No SD can be calculated')
#            elif len(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))])==1:
#                print('No data points in this epoch')
#            else:
#                print(epochs_included_labels[eachepoch]+' -> '+epochs_included_labels[eachepoch+1]+', FCF Peak: ',np.average(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))]),'+/-',np.sqrt(sum((source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))]-np.average(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))]))**2.0)/(len(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))])-1)),' = ',100*(np.sqrt(sum((source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))]-np.average(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))]))**2.0)/(len(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))])-1)))/(np.average(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))])))
#
#    
#    # FCF_arcsec plot!
#
#    
#    plt.clf()
#
#    colourdummy=0
#    for eachsource in sorted(INCLUDESOURCES):
#        already_labelled = 0
#        each_source_colour = eachsource_colours[colourdummy]
#        colourdummy=colourdummy+1
#        for eachepoch in range(len(epochs_included_labels)-1): # Avoid the epoch: "End"
#            if len(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)])>0:
#                if already_labelled==0:
#                    plt.scatter(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_MJD_'+str(eachepoch)],source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)],color=each_source_colour,marker=epoch_marker[eachepoch],label=eachsource)
#                    already_labelled = 1
#                else:
#                    plt.scatter(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_MJD_'+str(eachepoch)],source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)],color=each_source_colour,marker=epoch_marker[eachepoch])
#
#    for eachepoch in range(len(epochs_included)):
#        plt.axvline(x=epochs_included[eachepoch],linestyle='dashed',color='k',linewidth=2)
#        plt.text(epochs_included[eachepoch],maxylimarcsec-(0.04*(maxylimarcsec-minylimarcsec))*(eachepoch+1),epochs_included_labels[eachepoch],clip_on=True)
#
#    plt.legend(loc='upper left')
#    plt.xlim(xmin=startdate_mjd-20,xmax=enddate_mjd+20)
#    plt.ylim(ymin=minylimarcsec,ymax=maxylimarcsec)
#    plt.xlabel('MJD',fontsize=15)
#    plt.ylabel('FCF_arcsec (Jy/pW/arcsec^2)',fontsize=15)
#    plt.show()
#    plt.clf()
#
#
#    ######## PRINT FCF INFORMATION ############
#
#
#    for eachsource in INCLUDESOURCES:
#        print('\n\n#############\n\n'+eachsource+':\n\n')
#        for eachepoch in range(len(epochs_included)-1): # Avoid the epoch: "End"
#            if len(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))])==1:
#                print(epochs_included_labels[eachepoch]+' -> '+epochs_included_labels[eachepoch+1]+', FCF Arcsec: ',np.average(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))]),'+/- Only 1 data point! No SD can be calculated')
#            elif len(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))])==0:
#                print('No data points in this epoch')
#            else:
#                print(epochs_included_labels[eachepoch]+' -> '+epochs_included_labels[eachepoch+1]+', FCF Arcsec: ',np.average(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))]),'+/-',np.sqrt(sum((source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))]-np.average(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))]))**2.0)/(len(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))])-1)),' = ',100*(np.sqrt(sum((source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))]-np.average(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))]))**2.0)/(len(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))])-1)))/(np.average(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))])))
#
#    print('\n\n')
#
#
#
#    # Effective Beam size info
#    for eachepoch in range(len(epochs_included)-1):
#        print('\n\n#############')
#        print('Epoch: '+epochs_included_labels[eachepoch])
#        all_FCF_peaks_beamsize_thisepoch   = []
#        all_FCF_arcsecs_beamsize_thisepoch = []
#        for eachsource in INCLUDESOURCES:
#
#            FCF_goodind_for_beamsize = np.where(np.logical_and(np.logical_and(np.array(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)])>FCF_peak_low_cutoff,np.array(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)])<FCF_peak_high_cutoff),np.logical_and(np.array(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)])>FCF_arcsec_low_cutoff,np.array(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)])<FCF_arcsec_high_cutoff)))
#            if len(np.array(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)])[FCF_goodind_for_beamsize])>1:
#                plt.clf()
#                plt.scatter(np.array(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)])[FCF_goodind_for_beamsize],np.array(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)])[FCF_goodind_for_beamsize])
#                m,b = np.polyfit(np.array(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)])[FCF_goodind_for_beamsize],np.array(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][FCF_goodind_for_beamsize]),1)
#                plt.plot(np.array(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)])[FCF_goodind_for_beamsize],m*np.array(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)])[FCF_goodind_for_beamsize]+b,color='k',linestyle='dashed',linewidth=2.0)
#                plt.text(np.nanmean(np.array(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)])[FCF_goodind_for_beamsize])-np.nanstd(np.array(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)])[FCF_goodind_for_beamsize]),np.nanmean(np.array(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)])[FCF_goodind_for_beamsize])+2*np.nanstd(np.array(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)])[FCF_goodind_for_beamsize]),eachsource+'\nBeam_eff = '+str(round(np.sqrt(m/1.133),3))+'"',fontsize=15)
#                plt.xlabel('FCF_arcsec')
#                plt.ylabel('FCF_peak')
#                plt.show()
#                plt.clf()
#                for eachFCFpeak in np.array(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)])[FCF_goodind_for_beamsize]:
#                    all_FCF_peaks_beamsize_thisepoch.append(eachFCFpeak)
#                for eachFCFarcsec in np.array(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)])[FCF_goodind_for_beamsize]:
#                    all_FCF_arcsecs_beamsize_thisepoch.append(eachFCFarcsec) 
#            else:
#                print('%%%%%% No Data for this source in this Epoch %%%%%%')
#
#        if len(np.array(all_FCF_arcsecs_beamsize_thisepoch))>1:
#            plt.clf()
#            plt.scatter(np.array(all_FCF_arcsecs_beamsize_thisepoch),np.array(all_FCF_peaks_beamsize_thisepoch))
#            m,b = np.polyfit(np.array(all_FCF_arcsecs_beamsize_thisepoch),np.array(all_FCF_peaks_beamsize_thisepoch),1)
#            plt.plot(np.array(all_FCF_arcsecs_beamsize_thisepoch),m*np.array(all_FCF_arcsecs_beamsize_thisepoch)+b,color='k',linestyle='dashed',linewidth=2.0)
#            plt.text(np.nanmean(np.array(all_FCF_arcsecs_beamsize_thisepoch))-np.nanstd(np.array(all_FCF_arcsecs_beamsize_thisepoch)),np.nanmean(np.array(all_FCF_peaks_beamsize_thisepoch))+2*np.nanstd(np.array(all_FCF_peaks_beamsize_thisepoch)),'Total\nBeam_eff = '+str(round(np.sqrt(m/1.133),3))+'"',fontsize=15)
#            plt.xlabel('FCF_arcsec')
#            plt.ylabel('FCF_peak')
#            plt.show()
#            plt.clf()
#        else:
#            print('%%%%%% No Beam Size Calculated for this Epoch %%%%%%')
#    #############
#
#    FCF_dictionary = {}
#    for eachsource in INCLUDESOURCES:
#        FCF_dictionary[eachsource]={}
#        FCF_dictionary[eachsource]['N_measures_peaks']   = []
#        FCF_dictionary[eachsource]['N_measures_arcsecs'] = []
#        FCF_dictionary[eachsource]['FCF_peak']           = []
#        FCF_dictionary[eachsource]['FCF_peak_err']       = []
#        FCF_dictionary[eachsource]['FCF_arcsec']         = []
#        FCF_dictionary[eachsource]['FCF_arcsec_err']     = []
#        FCF_dictionary[eachsource]['Epoch_Labels']       = epochs_included_labels
#        for eachepoch in range(len(epochs_included)-1):
#            FCF_dictionary[eachsource]['N_measures_peaks'].append(len(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))]))
#            FCF_dictionary[eachsource]['FCF_peak'].append(np.average(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))]))
#            if len(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))])<2:
#                FCF_dictionary[eachsource]['FCF_peak_err'].append(np.nan)
#            else:
#                FCF_dictionary[eachsource]['FCF_peak_err'].append(np.sqrt(sum((source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))]-np.average(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))]))**2.0)/(len(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]>FCF_peak_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_peaks_'+str(eachepoch)]<FCF_peak_high_cutoff))])-1)))
#            FCF_dictionary[eachsource]['N_measures_arcsecs'].append(len(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))]))
#            FCF_dictionary[eachsource]['FCF_arcsec'].append(np.average(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))]))
#            if len(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))])<2:
#                FCF_dictionary[eachsource]['FCF_arcsec_err'].append(np.nan)
#            else:
#                FCF_dictionary[eachsource]['FCF_arcsec_err'].append(np.sqrt(sum((source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))]-np.average(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))]))**2.0)/(len(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)][np.where(np.logical_and(source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]>FCF_arcsec_low_cutoff,source_information_dict_by_epoch[eachsource]['FCF_arcsecs_'+str(eachepoch)]<FCF_arcsec_high_cutoff))])-1)))
#
#    return(FCF_dictionary)
