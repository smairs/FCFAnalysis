# Now, with this TauRelResults code, the goal is to load in all the information 
# derived above for different calibrator sources (primarily Uranus and CRL2688)
# and:
#
# 1. Find the flattest, physical (defined using CSO transmission model) peak flux 
#    versus transmission slope
#
# 2. Derive a new, robust opacity relation to use for all calibrator maps
#
# 3. Use this new opacity relation to convert pW to Jy/beam and Jy/arcsec
#    based on the known fluxes of the sources.
#
# 4. Report the new FCFs and see how they have changed over time
#
################################################
# Steve Mairs - EAO Support Scientist - 2017.
################################################     
################################################

def TauRelResults(resultsdict,TransvsPWV,wave,pix_size,unphys_limit,tau_cutoff,evening_cutoff,morning_cutoff,airmass_cutoff,FWHM_cutoff):

    '''
    For instance, 

    resultsdict = "TauRelPipeline_FullResults_CRL618_26p0to26p5_0p012to0p012_20162017.bin"

    '''

    # Import all the necessary modules
    
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import astropy
    import pickle
    import datetime
    import os
    from astropy.time import Time
    
    # Set up the format of the plots so they look pretty
    
    matplotlib.rcParams['figure.figsize']=(9,9)
    matplotlib.rcParams['xtick.labelsize'] = 15
    matplotlib.rcParams['ytick.labelsize'] = 15
    matplotlib.rcParams['axes.labelsize'] = 15
    matplotlib.rcParams['lines.linewidth'] = 2
    
    
    ###############################################
    ###############################################
    ############## Make the dictionary ############
    ############## Easier to work with ############
    ###############################################
    ###############################################
    
    # The first thing I'll do is define a class for ease of data handling. Then, I will be able to access
    # all data and metadata easily - also print out descriptions of any specifc combination of Opacity
    # Relation coefficients. When I say coefficients throughout this document, I am referring to "a" and "b"
    # in the tau relations: tau_wavelength = a(tau_225 - b).
    
    def friendlyresults(fname):
        
        # The results were saved in binary format using pickle
        import pickle
        from copy import deepcopy
    
        class TauRelResults:
            
            # Define everything that is in the dictionary
            
            # "a" = coeff1 and "b" = coeff2 (see above)
            coeff1=0
            coeff2=0
            
            # The Calibrator we
            # are currently analysing
            source = ""
            
            # The peak fluxes (in pW)
            # of the calibrators after
            # the new extinction correction
            # has been applied
            peak_fluxes=[]
            
            # The sum (in picowatts)
            # of all the pixels contained
            # within the Guassian out to
            # a level of 0.5*rms after
            # the new extinction correction
            # has been applied
            total_fluxes=[]
            
            # The derived (GaussClumps)
            # FWHM values of the 1st dimension
            # of the 2D Gaussian
            FWHM1s=[]
            
            # The average FWHM
            # of the first dimension
            # in the 2D Gaussian
            # fit performed by GaussClumps
            averageFWHM1=0
            
            # The derived (GaussClumps)
            # FWHM values of the second 
            # dimension
            FWHM2s=[]
            
            # The average FWHM in
            # the second dimension
            averageFWHM2=0
    
            # The Area of the aperture
            # Used to derive total flux
            areas=[]
    
            # The file names of the 
            # calibrator reductions
            fnames=[]
    
            # The absolute value of the
            # Calculated and Expected Transmission
            delta_trans=[]
    
            # The calculated Zenith transmission
            # Values
            transmissions=[]
            
            # The percentage of data points
            # that have an "unphysical" transmission
            # Here, "unphysical" is defined as
            # 5% different than the CSO Transmission
            # Curve at this wavelength
            transunphysical_percent=0
    
            # 225GHz Tau value
            # derived from 183 GHz WVM
            WVMTAUST=[]
            WVMTAUEN=[]
            
            # Time of PWV measurments
            # sometimes these occur
            # far outside a representative time
            # of the observation
            WVMTAUST_TIME=[]
            WVMTAUEN_TIME=[]
            
            # Airmass
            AMSTART=[]
            AMEND=[]
    
            # UTC observation start time
            OBSSTART=[]
            OBSEND=[]
            
            # Air temperature
            ATSTART=[]
            ATEND=[]
            
            # Now make a function that will find the slope of
            # a peak flux versus transmission plot. Remember
            # To include the airmass!
            
            def fitPFvsT(self):
                import numpy as np
                m,b = np.polyfit(self.transmissions**((self.amstart+self.amend)/2.0),self.peak_fluxes,1)
                return [m,b]
    
            # Now make a function that will find the slope of
            # a total flux versus transmission plot. Remember
    	# To include the airmass!
            
            def fitTFvsT(self):
                import numpy as np
                m,b = np.polyfit(self.transmissions**((self.amstart+self.amend)/2.0),self.total_fluxes,1)
                return [m,b]
            
            # Make a function that will give a description
            # of the results of running TauRelPipeline with
            # this specific set of coefficients.
            
            def description(self):
    
                 print('\n In this run, \n\n              a = %s and b = %s.\n\n              There are %.1f observations of %s included in this calculation.\n\n              The best fitting line for the Peak Flux verus Transmission function is:\n\n              Peak Flux = %.5f Transmission + (%.5f)\n\n              The average FWHM1 is %.3f"\n              The average FWHM2 is %.3f"\n\n              %.2f percent of the calculated transmissions are unphysical.\n\n'              % (self.coeff1,self.coeff2,len(self.peak_fluxes),
                 self.source,self.fitPFvsT()[0],
                 self.fitPFvsT()[1],self.averageFWHM1,
                 self.averageFWHM2,self.transunphysical_percent))
    
        # Now, load in the results dictionary file generated by TauRelPipeline.py
        # The dictionary is first organised by 'Run_0', 'Run_1', etc - one run for
        # each pairing of coefficients ("a" and "b")
        results_dict = pickle.load(open(fname,'rb'))
    
        # Now define each dictionary entry to be the above
        # class in a friendier, more accessible dictionary.
        # Note that it is impossible to pickle or json
        # a dictionary of classes.
        
        import numpy as np
        
        friendlydict={}
        for eachkey in results_dict:
            friendlydict[eachkey]=type('CopyOfTauRelResults', TauRelResults.__bases__, 
                                       dict(TauRelResults.__dict__))
            if isinstance(results_dict[eachkey]['coeff1'],list):
                friendlydict[eachkey].coeff1                  = results_dict[eachkey]['coeff1'][0]
                friendlydict[eachkey].coeff2                  = results_dict[eachkey]['coeff2'][0]
            else:
                friendlydict[eachkey].coeff1                  = results_dict[eachkey]['coeff1']
                friendlydict[eachkey].coeff2                  = results_dict[eachkey]['coeff2']
            friendlydict[eachkey].source                  = results_dict[eachkey]['source']
            friendlydict[eachkey].total_fluxes            = np.array(results_dict[eachkey]['total_fluxes'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].FWHM1s                  = np.array(results_dict[eachkey]['FWHM1s'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].FWHM2s                  = np.array(results_dict[eachkey]['FWHM2s'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].areas                   = np.array(results_dict[eachkey]['areas'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].averageFWHM1            = results_dict[eachkey]['averageFWHM1']
            friendlydict[eachkey].averageFWHM2            = results_dict[eachkey]['averageFWHM2']
            friendlydict[eachkey].transmissions           = np.array(results_dict[eachkey]['transmissions'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].transunphysical_percent = results_dict[eachkey]['trans_unphys_per']
            friendlydict[eachkey].delta_trans             = np.array(results_dict[eachkey]['delta_trans'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].wvmtaust                = np.array(results_dict[eachkey]['WVMTAUST'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].wvmtauen                = np.array(results_dict[eachkey]['WVMTAUEN'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].wvmtaust_time           = np.array(results_dict[eachkey]['WVMTAUST_TIME'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].wvmtauen_time           = np.array(results_dict[eachkey]['WVMTAUEN_TIME'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].amstart                 = np.array(results_dict[eachkey]['AMSTART'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].amend                   = np.array(results_dict[eachkey]['AMEND'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].obsstart                = np.array(results_dict[eachkey]['OBSSTART'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].obsend                  = np.array(results_dict[eachkey]['OBSEND'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].atstart                 = np.array(results_dict[eachkey]['ATSTART'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].atend                   = np.array(results_dict[eachkey]['ATEND'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].fnames                  = np.array(results_dict[eachkey]['fnames'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            friendlydict[eachkey].peak_fluxes             = np.array(results_dict[eachkey]['peak_fluxes'])[np.where(np.isnan(np.array(results_dict[eachkey]['peak_fluxes']))*-1+1>0)]
            
        return friendlydict               
    
    ######################################
    ######################################
    ######################################
    # Now, we can proceed with the analysis
    ######################################
    ######################################
    ######################################
    
    # Create a friendly dictionary of classes
    # full of all the TauRelPipeline results:
    calinfo = friendlyresults(resultsdict)
    
    import matplotlib.pyplot as plt
    import numpy as np
    from datetime import datetime
    from TauRelAnalysis_20171215 import FitCSOTransvsPWV, CSOtrans
    
    #TransvsPWV = FitCSOTransvsPWV(int(wave))
     
    ############
    ############
    ############
    # Check to see if there are WVM values that were taken several minutes away from the observation start/end time
        
    for eachkey in ['Run_0']:
            wvmtausttimes=[]
            obssttimes   =[]
            wvmtauentimes=[]
            obsentimes   =[]
            for i in range(len(calinfo[eachkey].wvmtaust_time)):
                now = datetime.now()
                # There is a weird error in the WVM raw data sometimes where the
		# Time is set to 0.00000000e+00 or something like 1858-11-17T00:00:00.00
		# For example, see raw data for 20130401, Scan 00060.
                if len(calinfo[eachkey].wvmtaust_time[i].split('T'))>1:
                    if int(calinfo[eachkey].wvmtaust_time[i].split('T')[0].split('-')[0])<2000:
                        print('BAD WVM TIMING DATA: ',calinfo[eachkey].fnames[i])
                    wvmtaust_time_str = calinfo[eachkey].wvmtaust_time[i].split('T')[1]
                    wvmtausthr  = int(float(wvmtaust_time_str.split(':')[0]))
                    wvmtaustmin = int(float(wvmtaust_time_str.split(':')[1]))
                    wvmtaustsec = int(float(wvmtaust_time_str.split(':')[2]))
                    wvmtaust_sec_since_midnight = (now.replace(hour=0, minute=0, second=0, microsecond=0) - 
                                               now.replace(hour=wvmtausthr, minute=wvmtaustmin, 
                                                           second=wvmtaustsec, microsecond=0)).total_seconds()
                    wvmtausttimes.append(wvmtaust_sec_since_midnight)
            
                    wvmtauen_time_str = calinfo[eachkey].wvmtauen_time[i].split('T')[1]
                    wvmtauenhr  = int(float(wvmtauen_time_str.split(':')[0]))
                    wvmtauenmin = int(float(wvmtauen_time_str.split(':')[1]))
                    wvmtauensec = int(float(wvmtauen_time_str.split(':')[2]))
                    wvmtauen_sec_since_midnight = (now.replace(hour=0, minute=0, second=0, microsecond=0) - 
                                               now.replace(hour=wvmtauenhr, minute=wvmtauenmin, 
                                                           second=wvmtauensec, microsecond=0)).total_seconds()
                    wvmtauentimes.append(wvmtauen_sec_since_midnight)

                else:
                    wvmtausttimes.append(-50000)
                    wvmtauentimes.append(-50000)
                    print('BAD WVM TIMING DATA: ',calinfo[eachkey].fnames[i])
                
                obsst_time_str = calinfo[eachkey].obsstart[i].split('T')[1]
                obssthr  = int(float(obsst_time_str.split(':')[0]))
                obsstmin = int(float(obsst_time_str.split(':')[1]))
                obsstsec = int(float(obsst_time_str.split(':')[2]))
                obsst_sec_since_midnight = (now.replace(hour=0, minute=0, second=0, microsecond=0) - 
                                           now.replace(hour=obssthr, minute=obsstmin, 
                                                       second=obsstsec, microsecond=0)).total_seconds()
                obssttimes.append(obsst_sec_since_midnight)
             
                obsen_time_str = calinfo[eachkey].obsend[i].split('T')[1]
                obsenhr  = int(float(obsen_time_str.split(':')[0]))
                obsenmin = int(float(obsen_time_str.split(':')[1]))
                obsensec = int(float(obsen_time_str.split(':')[2]))
                obsen_sec_since_midnight = (now.replace(hour=0, minute=0, second=0, microsecond=0) - 
                                           now.replace(hour=obsenhr, minute=obsenmin, 
                                                       second=obsensec, microsecond=0)).total_seconds()
                obsentimes.append(obsen_sec_since_midnight)
    
    
            print ('\nThere are '+str(len(abs(np.array(obssttimes)-np.array(wvmtausttimes))[np.where(abs(np.array(obssttimes)-np.array(wvmtausttimes))>60.0)]))+' WVM data points taken more than 1 minute away from the start of the observation\n')
            print ('\nThere are '+str(len(abs(np.array(obsentimes)-np.array(wvmtauentimes))[np.where(abs(np.array(obsentimes)-np.array(wvmtauentimes))>60.0)]))+' WVM data points taken more than 1 minute away from the end of the observation\n')
    
    # Uncomment below to plot!
    
            #plt.hist(abs(np.array(obssttimes)-np.array(wvmtausttimes)))
            #plt.xlabel('Difference Between OBSSTART and WVMTAUST (sec)')
            #plt.show()
            
            #plt.hist(abs(np.array(obsentimes)-np.array(wvmtauentimes)))
            #plt.xlabel('Difference Between OBSEND and WVMTAUEN (sec)')
            #plt.show()
    
    ##############
    ##############
    ##############
    
    #########################################################################################
    #########################################################################################
    #########################################################################################
    # Next, we apply the series of  user-defined constraints so that we are sure we are dealing with
    # the stable part of the night:
    #########################################################################################
    #########################################################################################
    #########################################################################################
    
    # So get a list of TAU225 and OBSstart times:
    # And define observations to fall within the constraints
    # to be in a "good" index: goodind
    tau225_for_goodind  = []
    obsst_for_goodind   = []
    airmass_for_goodind = []
    FWHM1_for_goodind   = []
    FWHM2_for_goodind   = []
    expected_transmissions = {} #CSOtrans(int(wave),tau225,TransvsPWV)
    for eachkey in calinfo.keys():
        expected_transmissions[eachkey]=[]
    
    
    for eachkey in ['Run_0']:
        for i in range(len(calinfo[eachkey].wvmtaust)):
            tau225_for_goodind.append((calinfo[eachkey].wvmtaust[i]+calinfo[eachkey].wvmtauen[i])/2.0)
            FWHM1_for_goodind.append(calinfo[eachkey].FWHM1s[i])
            FWHM2_for_goodind.append(calinfo[eachkey].FWHM2s[i])
            obsst_time_str = calinfo[eachkey].obsstart[i].split('T')[1]
            obssthr  = int(float(obsst_time_str.split(':')[0]))
            obsstmin = int(float(obsst_time_str.split(':')[1]))
            obsstsec = int(float(obsst_time_str.split(':')[2]))
            obsst_for_goodind.append(obssthr)
            airmass_for_goodind.append((calinfo[eachkey].amend[i]+calinfo[eachkey].amstart[i])/2.0)
        
    tau225_for_goodind            = np.array(tau225_for_goodind)
    obsst_for_goodind             = np.array(obsst_for_goodind)
    airmass_for_goodind           = np.array(airmass_for_goodind) 
    FWHM1_for_goodind             = np.array(FWHM1_for_goodind)
    FWHM2_for_goodind             = np.array(FWHM2_for_goodind)
    # These 2 are from before, same order
    obsst_minus_wvmst_for_goodind = np.array(obssttimes)-np.array(wvmtausttimes)
    obsen_minus_wvmen_for_goodind = np.array(obsentimes)-np.array(wvmtauentimes)

    goodind = np.where(np.logical_and(np.logical_and(np.logical_and(FWHM1_for_goodind<FWHM_cutoff,FWHM2_for_goodind<FWHM_cutoff),np.logical_and(airmass_for_goodind<airmass_cutoff,np.logical_and(tau225_for_goodind<tau_cutoff,np.logical_and(obsst_for_goodind>=evening_cutoff,obsst_for_goodind<morning_cutoff)))),np.logical_and(obsst_minus_wvmst_for_goodind<60.0,obsen_minus_wvmen_for_goodind<60.0)))
    
    print ('\n\n '+str(len(goodind[0]))+' observations survived the constraints!\n\n')

    if len(goodind[0])>1:
    
        from scipy.stats import norm
        
        mu,std = norm.fit(calinfo[eachkey].peak_fluxes[goodind])
        good_peakflux_ind = np.where(np.logical_and(calinfo[eachkey].peak_fluxes[goodind]<mu+3*std,calinfo[eachkey].peak_fluxes[goodind]>mu-3*std))
        bad_peakflux_ind = np.where(np.logical_or(calinfo[eachkey].peak_fluxes[goodind]>mu+3*std,calinfo[eachkey].peak_fluxes[goodind]<mu-3*std))
        
        # Throw away everything that does not conform to the constraints!
        for eachkey in calinfo.keys():
            for everygoodfile in range(len(calinfo[eachkey].wvmtaust[goodind][good_peakflux_ind])):
                expected_transmissions[eachkey].append(CSOtrans(int(wave),(calinfo[eachkey].wvmtaust[goodind][good_peakflux_ind][everygoodfile]+calinfo[eachkey].wvmtauen[goodind][good_peakflux_ind][everygoodfile])/2.0,TransvsPWV))
            expected_transmissions[eachkey]  = np.array(expected_transmissions[eachkey])
            calinfo[eachkey].transmissions   = calinfo[eachkey].transmissions[goodind][good_peakflux_ind]
            calinfo[eachkey].FWHM1s          = calinfo[eachkey].FWHM1s[goodind][good_peakflux_ind]
            calinfo[eachkey].FWHM2s          = calinfo[eachkey].FWHM2s[goodind][good_peakflux_ind]
            calinfo[eachkey].areas           = calinfo[eachkey].areas[goodind][good_peakflux_ind]
            calinfo[eachkey].total_fluxes    = calinfo[eachkey].total_fluxes[goodind][good_peakflux_ind]
            calinfo[eachkey].peak_fluxes     = calinfo[eachkey].peak_fluxes[goodind][good_peakflux_ind]
            calinfo[eachkey].wvmtaust        = calinfo[eachkey].wvmtaust[goodind][good_peakflux_ind]
            calinfo[eachkey].wvmtauen        = calinfo[eachkey].wvmtauen[goodind][good_peakflux_ind]
            calinfo[eachkey].wvmtaust_time   = calinfo[eachkey].wvmtaust_time[goodind][good_peakflux_ind]
            calinfo[eachkey].wvmtauen_time   = calinfo[eachkey].wvmtauen_time[goodind][good_peakflux_ind]
            calinfo[eachkey].amstart         = calinfo[eachkey].amstart[goodind][good_peakflux_ind]
            calinfo[eachkey].amend           = calinfo[eachkey].amend[goodind][good_peakflux_ind]
            calinfo[eachkey].obsstart        = calinfo[eachkey].obsstart[goodind][good_peakflux_ind]
            calinfo[eachkey].obsend          = calinfo[eachkey].obsend[goodind][good_peakflux_ind]
            calinfo[eachkey].atstart         = calinfo[eachkey].atstart[goodind][good_peakflux_ind]
            calinfo[eachkey].atend           = calinfo[eachkey].atend[goodind][good_peakflux_ind]
            calinfo[eachkey].fnames          = calinfo[eachkey].fnames[goodind][good_peakflux_ind]
            calinfo[eachkey].delta_trans     = calinfo[eachkey].delta_trans[goodind][good_peakflux_ind]
        
        
        # Get lists of the coefficients that were used in each run of
        # TauRelPipeline.py
        
        coeff1s=[]
        coeff2s=[]
        for eachkey in calinfo.keys():
            coeff1s.append(calinfo[eachkey].coeff1)
            coeff2s.append(calinfo[eachkey].coeff2)
        
        
        ###############################
        ###############################
        ###############################
        ###### Begin Plotting #########
        ###############################
        ###############################
        ###############################
        
        ################################
        ################################
        # Now Plot the Peak Flux Versus Transmission Plots - REMEMBER TO INCLUDE AIRMASS!
        
        for eachkey in calinfo.keys():
            plt.plot(calinfo[eachkey].transmissions**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0),
                     (calinfo[eachkey]().fitPFvsT()[0]*np.array(calinfo[eachkey].transmissions**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0))+
                     calinfo[eachkey]().fitPFvsT()[1]))
        
            # Write which coefficients produced each line!
            plt.text(min(calinfo[eachkey].transmissions**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0)),min((calinfo[eachkey]().fitPFvsT()[0]*
                     np.array(calinfo[eachkey].transmissions**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0))+
                     calinfo[eachkey]().fitPFvsT()[1]))
                     +0.001,str(calinfo[eachkey].coeff1)+','+str(calinfo[eachkey].coeff2),size=15)
        
        plt.xlabel('Transmission= exp(-tau*airmass) (%)')
        plt.ylabel('Peak Flux (pW)')
        plt.show()
        plt.clf()
        
        
        
        ########################
        ########################
        # Now show a plot of the percentage of
        # unphysical values with respect to COEFF1
        
        #for eachkey in calinfo.keys():
        #    delta_trans_num = len(calinfo[eachkey].wvmtaust)
        #    if calinfo[eachkey].coeff2 == min(coeff2s):
        #        plt.scatter(calinfo[eachkey].coeff1,len(calinfo[eachkey].delta_trans[np.where(
        #            calinfo[eachkey].delta_trans>unphys_limit)])/len(calinfo[eachkey].delta_trans),color='r')
        #    else:
        #        plt.scatter(calinfo[eachkey].coeff1,len(calinfo[eachkey].delta_trans[np.where(
        #            calinfo[eachkey].delta_trans>unphys_limit)])/len(calinfo[eachkey].delta_trans),color='b')
        #        
        #plt.xlabel('COEFF1')
        #plt.ylabel('% Unphysical (|Trans_calc - Trans_exp| > '+str(unphys_limit)+')')
        #plt.show()
        #plt.clf()
        #
        #
        ########################
        ########################
        # Now show a plot of the percentage of
        # unphysical values with respect to COEFF2
        
        #for eachkey in calinfo.keys():
        #    if calinfo[eachkey].coeff2 == min(coeff2s):
        #        plt.scatter(calinfo[eachkey].coeff2,len(calinfo[eachkey].delta_trans[np.where(
        #            calinfo[eachkey].delta_trans>unphys_limit)])/len(calinfo[eachkey].delta_trans),color='r')
        #    else:
        #        plt.scatter(calinfo[eachkey].coeff2,len(calinfo[eachkey].delta_trans[np.where(
        #            calinfo[eachkey].delta_trans>unphys_limit)])/len(calinfo[eachkey].delta_trans),color='b')
        #        
        #plt.xlabel('COEFF2')
        #plt.ylabel('% Unphysical (|Trans_calc - Trans_exp| > '+str(unphys_limit)+')')
        #plt.show()
        #plt.clf()
        
        ##############################################
        ##############################################
        # Plot a histogram of the absolute value of the
        # difference between the observed transmission and the expected transmission.
        #
        
        for eachkey in calinfo.keys():
            plt.hist(calinfo[eachkey].delta_trans,bins=5)
            plt.suptitle('COEFF1 = '+str(calinfo[eachkey].coeff1)+', COEFF2 = '+str(calinfo[eachkey].coeff2),size=18)
            plt.axvline(x=0.05,linewidth=2.0,color='k')
            plt.text(0.052,50,'Unphysical -->',size=20)
            plt.ylabel('Count')
            plt.xlabel('|obs_trans - exp_trans|')
            #plt.show()
            plt.clf()
        
        ###############################################
        ###############################################
        # Plot the measured slope against the percentage of points that are unphysical
        
        for eachkey in calinfo.keys():
            transdiff = expected_transmissions[eachkey]**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0)-calinfo[eachkey].transmissions**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0)
            #plt.scatter(calinfo[eachkey]().fitPFvsT()[0],len(calinfo[eachkey].delta_trans[np.where(
            #            calinfo[eachkey].delta_trans>unphys_limit)])/len(calinfo[eachkey].delta_trans),color='b')
        
            plt.scatter(calinfo[eachkey]().fitPFvsT()[0],len(transdiff[np.where(abs(transdiff)>unphys_limit)])/len(transdiff),color='b')
        
            #plt.text(calinfo[eachkey]().fitPFvsT()[0],len(calinfo[eachkey].delta_trans[np.where(
            #            calinfo[eachkey].delta_trans>unphys_limit)])/len(calinfo[eachkey].delta_trans),str(calinfo[eachkey].coeff1)+','+str(
            #                                                                                            calinfo[eachkey].coeff2))
        
            plt.text(calinfo[eachkey]().fitPFvsT()[0],len(transdiff[np.where(abs(transdiff)>unphys_limit)])/len(transdiff),str(calinfo[eachkey].coeff1)+','+str(
                                                                                                                calinfo[eachkey].coeff2))
            
            plt.axhline(y=unphys_limit,linewidth=2.0,linestyle='dashed')
        
            #print(str(calinfo[eachkey].coeff1)+','+str(calinfo[eachkey].coeff2)+'	PF: '+str(calinfo[eachkey]().fitPFvsT()[0]))
        
        plt.ylabel('% Obs Unphysical')
        plt.xlabel('Peak Flux vs Transmission Slope (pW/%)')
        plt.show()
        plt.clf()
    
        tautimesairmass= []
        transdiffs     = []
        for eachkey in calinfo.keys():
            plt.scatter(((calinfo[eachkey].wvmtaust+calinfo[eachkey].wvmtauen)/2.0)*(calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0,expected_transmissions[eachkey]**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0)-calinfo[eachkey].transmissions**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0))
            for i in range(len((calinfo[eachkey].wvmtaust+calinfo[eachkey].wvmtauen)/2.0)):
                plt.text((((calinfo[eachkey].wvmtaust+calinfo[eachkey].wvmtauen)/2.0)*(calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0)[i],(expected_transmissions[eachkey]**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0)-calinfo[eachkey].transmissions**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0))[i],str((((calinfo[eachkey].wvmtaust+calinfo[eachkey].wvmtauen)/2.0)[i]))+','+str(((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0)[i]))
            plt.ylabel('transdiff (expected minus measured)')
            plt.xlabel('tau*airmass')
            plt.suptitle(str(calinfo[eachkey].coeff1)+','+str(calinfo[eachkey].coeff2))
            plt.show()
            plt.clf()
        ##################################################
        ##################################################
        ##################################################
        ### Now, out of the coeff1+coeff2 combinations
        ### that have the minimum number of unphysical
        ### points (0 in the ideal case we are going for)
        ### - Find the one with the smallest slope
        ### and choose it as the "Best Combination"
        ### for the FCF_peak opacity relation
        ##################################################
        ##################################################
        ##################################################
        
        slopes       = []
        unphys       = []
        coeff1s_temp = []
        coeff2s_temp = []
        keys         = []
        for eachkey in calinfo.keys():
            slopes.append(calinfo[eachkey]().fitPFvsT()[0])
            unphys.append(len(calinfo[eachkey].delta_trans[np.where(
                    calinfo[eachkey].delta_trans>unphys_limit)])/len(calinfo[eachkey].delta_trans))
            coeff1s_temp.append(calinfo[eachkey].coeff1)
            coeff2s_temp.append(calinfo[eachkey].coeff2)
            keys.append(eachkey)
        slopes       = np.array(slopes)
        unphys       = np.array(unphys)
        coeff1s_temp = np.array(coeff1s_temp)
        coeff2s_temp = np.array(coeff2s_temp)
        keys         = np.array(keys)
        
        unphysical = 1
        while unphysical == 1:
            best_ind = np.argmin(abs(slopes))
            percentage_unphysical = unphys[best_ind]
            if percentage_unphysical > min(unphys):
                coeff1s_temp = np.delete(coeff1s_temp,best_ind)
                coeff2s_temp = np.delete(coeff2s_temp,best_ind)
                slopes       = np.delete(slopes,best_ind)
                unphys       = np.delete(unphys,best_ind)
                keys         = np.delete(keys,best_ind)
            #else:
            #    if slopes[best_ind]<0:
            #        coeff1s_temp = np.delete(coeff1s_temp,best_ind)
            #        coeff2s_temp = np.delete(coeff2s_temp,best_ind)
            #        slopes       = np.delete(slopes,best_ind)
            #        unphys       = np.delete(unphys,best_ind)
            #        keys         = np.delete(keys,best_ind)
            else:
                unphysical = 0
                
        #print ('\n\nBest Coefficients: \n\n')
        #print ('\ntau_'+wave+' = '+str(round(coeff1s_temp[best_ind],2))+' x (tau_225 -'+str(round(coeff2s_temp[best_ind],5))+')\n\n')
        #print ('\n\nDictionary key: '+keys[best_ind]+'\n\n')
        
        a_peak = coeff1s_temp[best_ind]
        b_peak = coeff2s_temp[best_ind]
        bestkey_peak = keys[best_ind]
       
        source = calinfo[bestkey_peak].source
    
        #######
        # Make a plot of how changing the tau relation from the nominal coefficients
        # to these new coeffcients would change the FCF
        #######
        
        airmassmodel = np.arange(1.0,2.5,0.0006)
        tau225model  = np.arange(0.0,0.25,0.0001)
        
        if wave=='850':
            transmissionmodel_nominal = np.exp(-4.6*(tau225model-0.0043)*airmassmodel)
            transmissionmodel_new     = np.exp(-1*calinfo[bestkey_peak].coeff1*(tau225model-calinfo[bestkey_peak].coeff2)*airmassmodel)
        else:
            transmissionmodel_nominal = np.exp(-26.0*(tau225model-0.012)*airmassmodel)
            transmissionmodel_new     = np.exp(-1*calinfo[bestkey_peak].coeff1*(tau225model-calinfo[bestkey_peak].coeff2)*airmassmodel)
        
        plt.scatter(tau225model*airmassmodel,transmissionmodel_nominal,color='b',label='Nominal')
        plt.scatter(tau225model*airmassmodel,transmissionmodel_nominal,color='r',label='New')
        plt.xlabel('Tau225*Airmass')
        plt.ylabel('Transmission')
        #plt.show()
        plt.clf()
        
        plt.scatter(tau225model*airmassmodel,transmissionmodel_new/transmissionmodel_nominal,color='r')
        plt.axvline(x=0.065*1.2,label='Typical Band 2, Airmass = 1.2',color='k',linestyle='dashed')
        plt.xlabel('Tau225*Airmass')
        plt.ylabel('FCFnew/FCFnominal')
        plt.legend(loc='upper right')
        #plt.show()
        plt.clf()
        
            
        #################################################
        #################################################
        #################################################
        
        #################################################
        #################################################
        ########### Now Plot the Total Flux #############
        ###########   Against Transmission  #############
        ###########  For the FCF_arcsec     #############
        ########### Does it agree with peak?#############
        #################################################
        #################################################
        
        for eachkey in calinfo.keys():
            plt.plot(calinfo[eachkey].transmissions**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0),
                     (calinfo[eachkey]().fitTFvsT()[0]*np.array(calinfo[eachkey].transmissions**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0))+
                     calinfo[eachkey]().fitTFvsT()[1]))
        
            # Write the COEFF1+COEFF2 pair over the curve
            plt.text(min(calinfo[eachkey].transmissions**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0)),min((
                     calinfo[eachkey]().fitTFvsT()[0]*np.array(calinfo[eachkey].transmissions**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0))+
                     calinfo[eachkey]().fitTFvsT()[1]))
                     +0.002,str(calinfo[eachkey].coeff1)+','+str(calinfo[eachkey].coeff2),size=15)
        
            #print(str(calinfo[eachkey].coeff1)+','+str(calinfo[eachkey].coeff2)+'       TF: '+str(calinfo[eachkey]().fitTFvsT()[0]))
        
        plt.xlabel('Transmission (%)')
        plt.ylabel('Total Flux (pW)')
        #plt.show()
        plt.clf()
        
        ###################################################
        ###################################################
        ###################################################
        # Plot the measured slope against the percentage of points that are unphysical
        
        for eachkey in calinfo.keys():
            #plt.scatter(calinfo[eachkey]().fitTFvsT()[0],len(calinfo[eachkey].delta_trans[np.where(
            #            calinfo[eachkey].delta_trans>unphys_limit)])/len(calinfo[eachkey].delta_trans),color='b')
        
            #plt.text(calinfo[eachkey]().fitTFvsT()[0],len(calinfo[eachkey].delta_trans[np.where(
            #            calinfo[eachkey].delta_trans>unphys_limit)])/len(calinfo[eachkey].delta_trans),str(calinfo[eachkey].coeff1)+','+str(
            #                                                                                            calinfo[eachkey].coeff2))
        
            transdiff = expected_transmissions[eachkey]**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0)-calinfo[eachkey].transmissions**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0)
        
            plt.scatter(calinfo[eachkey]().fitTFvsT()[0],len(transdiff[np.where(abs(transdiff)>unphys_limit)])/len(transdiff),color='b')
        
            plt.text(calinfo[eachkey]().fitTFvsT()[0],len(transdiff[np.where(abs(transdiff)>unphys_limit)])/len(transdiff),str(calinfo[eachkey].coeff1)+','+str(
                                                                                                                calinfo[eachkey].coeff2))
        
            plt.axhline(y=unphys_limit,linewidth=2.0,linestyle='dashed')
        
            #print(str(calinfo[eachkey].coeff1)+','+str(calinfo[eachkey].coeff2)+'       PF: '+str(calinfo[eachkey]().fitPFvsT()[0]))
        
        plt.ylabel('% Obs Unphysical')
        plt.xlabel('Total Flux vs Transmission Slope (pW/%)')
        plt.show()
        plt.clf()
        
        ##################################################
        ##################################################
        ##################################################
        ### Now, out of the coeff1+coeff2 combinations
        ### that have the minimum number of unphysical
        ### points (0 in the ideal case we are going for)
        ### - Find the one with the smallest slope
        ### and choose it as the "Best Combination"
        ### for the FCF_arcsec opacity relation
        ##################################################
        ##################################################
        ##################################################
        
        slopes       = []
        unphys       = []
        coeff1s_temp = []
        coeff2s_temp = []
        keys         = []
        for eachkey in calinfo.keys():
            slopes.append(calinfo[eachkey]().fitTFvsT()[0])
            unphys.append(len(calinfo[eachkey].delta_trans[np.where(
                    calinfo[eachkey].delta_trans>unphys_limit)])/len(calinfo[eachkey].delta_trans))
            coeff1s_temp.append(calinfo[eachkey].coeff1)
            coeff2s_temp.append(calinfo[eachkey].coeff2)
            keys.append(eachkey)
        
        slopes       = np.array(slopes)
        unphys       = np.array(unphys)
        coeff1s_temp = np.array(coeff1s_temp)
        coeff2s_temp = np.array(coeff2s_temp)
        keys         = np.array(keys)
        
        unphysical = 1
        while unphysical ==1:
            best_ind = np.argmin(abs(slopes))
            percentage_unphysical = unphys[best_ind]
            if percentage_unphysical > min(unphys):
                coeff1s_temp = np.delete(coeff1s_temp,best_ind)
                coeff2s_temp = np.delete(coeff2s_temp,best_ind)
                slopes       = np.delete(slopes,best_ind)
                unphys       = np.delete(unphys,best_ind)
                keys         = np.delete(keys,best_ind)
            #else:
            #    if slopes[best_ind]<0:
            #        coeff1s_temp = np.delete(coeff1s_temp,best_ind)
            #        coeff2s_temp = np.delete(coeff2s_temp,best_ind)
            #        slopes       = np.delete(slopes,best_ind)
            #        unphys       = np.delete(unphys,best_ind)
            #        keys         = np.delete(keys,best_ind)
            else:
                unphysical = 0
        
        
        #print ('\n\nBest Coefficients: \n\n')
        #print ('\ntau_'+wave+' = '+str(round(coeff1s_temp[best_ind],2))+' x (tau_225 -'+str(round(coeff2s_temp[best_ind],5))+')\n\n')
        #print ('\n\nDictionary key: '+keys[best_ind]+'\n\n')
        
        a_total = coeff1s_temp[best_ind]
        b_total = coeff2s_temp[best_ind]
        bestkey_total = keys[best_ind]
       
        return(source,[a_peak,b_peak],bestkey_peak,[a_total,b_total],bestkey_total)

    else:
        source = calinfo[sorted(list(calinfo.keys()))[0]].source 
        print('\n\nLESS THAN 2 OBSERVATIONS SURVIVED CONSTRAINTS - CANNOT MEASURE A SLOPE OR DETERMINE BEST OPACITY RELATION!\n\n')
        return(source,[np.nan,np.nan],'Less than 2 obs',[np.nan,np.nan],'Less than 2 obs')


    ########### IDENTIFY_OUTLIERS_FIND_FCFs DOES THIS NOW - IGNORE BELOW EXCEPT FOR HISTORICAL REASONS #####################
#
#    #####################################################
#    #####################################################
#    #####################################################
#    
#    # Now we want to determine an FCF For EVERY observation and plot that over time and with Transmission
#    # But note that the total flux is assumed to be taken from a 60" aperture - we used Gaussclumps
#    # which never performs the integration over exactly 60 arcseconds - so we need to correct for this
#    # when we are calculating the FCF_arcsec value
#    
#    def total_flux_cor(aperture_diameter,wavelength):
#    
#    
#        import numpy as np
#    
#        diams = np.arange(20,105.1,5)
#        factors_450 = [0.72,0.79,0.85,0.88,0.91,0.94,0.97,0.98,1.00,1.01,1.02,1.02,1.02,1.03,1.03,1.03,1.03,1.03]
#        factors_850 = [0.69,0.79,0.85,0.88,0.91,0.94,0.97,0.99,1.00,1.01,1.02,1.03,1.04,1.05,1.05,1.06,1.07,1.08]
#    
#        if wavelength=='450':
#            p =np.polyfit(diams,factors_450,10)
#        else:
#            p =np.polyfit(diams,factors_850,10)
#    
#        z = np.poly1d(p)
#    
#        return (z(aperture_diameter))
#    
#    
#    
#    #Output From "fluxes":
#    #URANUS
#    #South pole is Earth-facing; Inclination Angle =  51.53 degrees
#    #Semi-diameter =  1.85 arcsecs    Solid angle =  2.54E-10 sterads
#    #Filter    Centre   Filter   Total    Flux in    Brightness         HPBW
#    #Wavel.     Freq     Width    Flux     beam      Temperature       assumed
#    #micron     (GHz)    (GHz)    (Jy)      (Jy)         (K)          (arcsecs)
#    #  850      350.0    30.0     71.32     69.38     82.7 +-  1.0       13.0
#    #  450      677.0    30.0    191.17    178.02     68.4 +-  0.5        7.9
#    
#    #From: Dempsey / 
#    #http://www.eaobservatory.org/jcmt/instrumentation/continuum/scuba-2/calibration/calibrators/#Secondary_Calibrators-2
#    #name       RA              DEC             Peak 850        Int. 850     Peak 450    Int. 450
#    #CRL 618 	04:42:53.67 	+36:06:53.17 	4.89 ± 0.24 	5.0 ± 0.2 	11.5 ± 1.4 	12.1 ± 1.05
#    #CRL 2688 	21:02:18.27 	+36:41:37.00 	5.64 ± 0.27 	6.13 ± 0.2 	24.9 ± 2.9 	29.1 ± 2.5
#    #Arp 220 	15:34:57.27 	+23:30:10.48 	0.79 ± 0.09 	0.81 ± 0.07 	5.2 ± 0.8 	5.4 ± 0.7
#    #HD 169142 	18:24:29.78 	-29:46:49.37 	0.52 ± 0.03 	0.58 ± 0.02 	2.21 ± 0.25 	2.78 ± 0.24
#    #HD 135344B 	15:15:48.44 	-37:09:16.02 	0.46 ± 0.03 	0.53 ± 0.02 	1.66 ± 0.27 	3.30 ± 0.34
#    #MWC 349 	20:32:45.53 	40:39:36.6 	2.21 ± 0.11 	2.19 ± 0.08 	3.4 ± 0.26 	3.2 ± 0.25
#    #V883 Ori 	05:38:19 	-07:02:22 	1.55 ± 0.09 	2.0 ± 0.07 	7.8 ± 1.0 	11.0 ± 0.94
#    #HL Tau 	04:31:38.44 	+18:13:57.65 	2.32 ± 0.01 	2.42 ± 0.08 	8.3 ± 1.03 	10.3 ± 0.86
#    
#    from starlink import fluxes # Thanks Sarah Graves!
#    
#    print (calinfo[bestkey_peak].source)
#    
#    if calinfo[bestkey_peak].source == 'URANUS':
#    
#        peak850_actual   = []
#        arcsec850_actual = []
#        peak450_actual   = []
#        arcsec450_actual = []
#        aperture_cor     = []
#    
#        for eachfile in range(len(calinfo['Run_0'].fnames)):
#            datestr_long     = calinfo['Run_0'].obsstart[eachfile]+'.00'
#            peak850_actual.append(fluxes.get_flux('URANUS',datestr_long,filter_=850).f_beam)
#            arcsec850_actual.append(fluxes.get_flux('URANUS',datestr_long,filter_=850).f_total) 
#            peak450_actual.append(fluxes.get_flux('URANUS',datestr_long,filter_=450).f_beam)     
#            arcsec450_actual.append(fluxes.get_flux('URANUS',datestr_long,filter_=450).f_total) 
#            aperture_cor.append(total_flux_cor(2*np.sqrt(calinfo[bestkey_total].areas[eachfile]/np.pi),wave))
#    
#        aperture_cor     = np.array(aperture_cor)
#        peak850_actual   = np.array(peak850_actual)    
#        arcsec850_actual = np.array(arcsec850_actual)
#        peak450_actual   = np.array(peak450_actual)
#        arcsec450_actual = np.array(arcsec450_actual)
#        
#        if wave == '850':    
#            peak_fcfs   = peak850_actual/calinfo[bestkey_peak].peak_fluxes
#            arcsec_fcfs = (arcsec850_actual/((calinfo[bestkey_total].total_fluxes/aperture_cor)*pix_size**2.0))
#        if wave == '450':
#            peak_fcfs   = peak450_actual/calinfo[bestkey_peak].peak_fluxes
#            arcsec_fcfs = (arcsec450_actual/((calinfo[bestkey_total].total_fluxes/aperture_cor)*pix_size**2.0))
#    
#    if calinfo[bestkey_peak].source == 'CRL2688':
#        CRL2688_peak_850   = 5.64
#        CRL2688_arcsec_850 = 6.13
#        CRL2688_peak_450   = 24.9
#        CRL2688_arcsec_450 = 29.1
#    
#    
#        aperture_cor     = []
#        for eachfile in range(len(calinfo['Run_0'].fnames)):
#            aperture_cor.append(total_flux_cor(2*np.sqrt(calinfo[bestkey_total].areas[eachfile]/np.pi),wave))
#        aperture_cor     = np.array(aperture_cor)
#    
#    
#        if wave == '850':    
#            peak_fcfs   = CRL2688_peak_850/calinfo[bestkey_peak].peak_fluxes
#            arcsec_fcfs = (CRL2688_arcsec_850/((calinfo[bestkey_total].total_fluxes/aperture_cor)*pix_size**2.0))
#        if wave == '450':
#            peak_fcfs   = CRL2688_peak_450/calinfo[bestkey_peak].peak_fluxes
#            arcsec_fcfs = (CRL2688_arcsec_450/((calinfo[bestkey_total].total_fluxes/aperture_cor)*pix_size**2.0))
#    
#    
#    if calinfo[bestkey_peak].source == 'CRL618':
#        CRL618_peak_850   = 4.89
#        CRL618_arcsec_850 = 5.00
#        CRL618_peak_450   = 11.5
#        CRL618_arcsec_450 = 12.1
#    
#    
#        aperture_cor     = []
#        for eachfile in range(len(calinfo['Run_0'].fnames)):
#            aperture_cor.append(total_flux_cor(2*np.sqrt(calinfo[bestkey_total].areas[eachfile]/np.pi),wave))
#        aperture_cor     = np.array(aperture_cor)
#    
#    
#        if wave == '850':
#            peak_fcfs   = CRL618_peak_850/calinfo[bestkey_peak].peak_fluxes
#            arcsec_fcfs = (CRL618_arcsec_850/((calinfo[bestkey_total].total_fluxes/aperture_cor)*pix_size**2.0))
#        if wave == '450':
#            peak_fcfs   = CRL618_peak_450/calinfo[bestkey_peak].peak_fluxes
#            arcsec_fcfs = (CRL618_arcsec_450/((calinfo[bestkey_total].total_fluxes/aperture_cor)*pix_size**2.0))
#    
#    
#    # In[ ]:
#    
#    # Get a list of dates to plot the FCFs against
#    
#    
#    # In[ ]:
#    
#    MJD_dates = []
#    for eachfile in calinfo['Run_0'].fnames:
#        date  = eachfile.split('_')[0].split('s')[1]
#        year  = date[0:4]
#        month = date[4:6]
#        day   = date[6:8]
#        datestr = year+'-'+month+'-'+day+'T'+'00:00:00'
#        t = Time(datestr, format='isot', scale='utc')
#        MJD_dates.append(t.mjd)
#    
#    MJD_dates = np.array(MJD_dates)
#    
#    plt.scatter(MJD_dates,calinfo[bestkey_peak].peak_fluxes)
#    plt.xlabel('MJD')
#    plt.ylabel('Peak Fluxes')
#    plt.show()
#    plt.clf()
#    
#    plt.scatter(calinfo[bestkey_peak].transmissions,calinfo[bestkey_peak].peak_fluxes)
#    plt.xlabel('Transmission (%)')
#    plt.ylabel('Peak Fluxes')
#    plt.show()
#    plt.clf()
#    
#    for eachkey in calinfo.keys():
#    
#        plt.scatter(expected_transmissions[eachkey]**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0),calinfo[eachkey].transmissions**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0))
#        plt.plot(expected_transmissions[eachkey]**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0),expected_transmissions[eachkey]**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0),linestyle='dashed',linewidth=2,color='k')
#        plt.xlabel('Expected Transmission from CSO Atmospheric Plotter (%)')
#        plt.ylabel('Actual Transmission (%)')
#        plt.suptitle(str(calinfo[eachkey].coeff1)+','+str(calinfo[eachkey].coeff2))
#        plt.show()
#        plt.clf()
#    
#        plt.hist(expected_transmissions[eachkey]**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0)-calinfo[eachkey].transmissions**((calinfo[eachkey].amstart+calinfo[eachkey].amend)/2.0))
#        plt.xlabel('Expected Transmission - Measured Transmission')
#        plt.ylabel('Number')
#        plt.suptitle(str(calinfo[eachkey].coeff1)+','+str(calinfo[eachkey].coeff2))
#        plt.show()
#        plt.clf()
#    
#    
#    sd_num=[]
#    for i in peak_fcfs:
#        sd_num.append((i-np.average(peak_fcfs))**2.0)
#    sd = np.sqrt(sum(sd_num)/(len(sd_num)-1))
#    
#    #print('')
#    #print(np.average(peak_fcfs[np.where(np.logical_and(peak_fcfs<np.average(peak_fcfs)+3*sd,peak_fcfs>np.average(peak_fcfs)-3*sd))]),sd)
#    #print('')
#    
#    plt.scatter(MJD_dates,peak_fcfs)
#    plt.axhline(y=np.average(peak_fcfs[np.where(np.logical_and(peak_fcfs<np.average(peak_fcfs)+3*sd,peak_fcfs>np.average(peak_fcfs)-3*sd))])+sd,color='k',linestyle='dashed')
#    plt.axhline(y=np.average(peak_fcfs[np.where(np.logical_and(peak_fcfs<np.average(peak_fcfs)+3*sd,peak_fcfs>np.average(peak_fcfs)-3*sd))])-sd,color='k',linestyle='dashed')
#    plt.axhline(y=np.average(peak_fcfs[np.where(np.logical_and(peak_fcfs<np.average(peak_fcfs)+3*sd,peak_fcfs>np.average(peak_fcfs)-3*sd))]),color='k',linestyle='dashed')
#    plt.ylabel('FCF_peak '+wave+' um')
#    plt.xlabel('MJD')
#    plt.show()
#    plt.clf()
#    
#    sd_num=[]
#    for i in arcsec_fcfs:
#        sd_num.append((i-np.average(arcsec_fcfs))**2.0)
#    sd = np.sqrt(sum(sd_num)/(len(sd_num)-1))
#    
#    #print('')
#    #print(np.average(arcsec_fcfs[np.where(np.logical_and(arcsec_fcfs<np.average(arcsec_fcfs)+3*sd,arcsec_fcfs>np.average(arcsec_fcfs)-3*sd))]),sd)
#    #print('')
#    
#    plt.scatter(calinfo[bestkey_total].transmissions,calinfo[bestkey_total].total_fluxes)
#    plt.xlabel('Transmission (%)')
#    plt.ylabel('Total Fluxes')
#    #plt.show()
#    plt.clf()
#    
#    plt.scatter(MJD_dates,arcsec_fcfs)
#    plt.axhline(y=np.average(arcsec_fcfs[np.where(np.logical_and(arcsec_fcfs<np.average(arcsec_fcfs)+3*sd,arcsec_fcfs>np.average(arcsec_fcfs)-3*sd))])+sd,linestyle='dashed',color='k')
#    plt.axhline(y=np.average(arcsec_fcfs[np.where(np.logical_and(arcsec_fcfs<np.average(arcsec_fcfs)+3*sd,arcsec_fcfs>np.average(arcsec_fcfs)-3*sd))])-sd,linestyle='dashed',color='k')
#    plt.axhline(y=np.average(arcsec_fcfs[np.where(np.logical_and(arcsec_fcfs<np.average(arcsec_fcfs)+3*sd,arcsec_fcfs>np.average(arcsec_fcfs)-3*sd))]),color='k',linestyle='dashed')
#    plt.ylabel('FCF_arcsec '+wave+' um')
#    plt.xlabel('MJD')
#    #plt.show()
#    plt.clf()
#    
#    
#    # Get a list of calculated transmissions to plot the FCFs against
#    
#    
#    
#    #for eachkey in calinfo.keys():
#    #    obssttimes   = []
#    #    obsentimes   = []
#    #    delta_transvals  = []
#    #    for i in range(len(calinfo[eachkey].wvmtaust_time)):
#    #            now = datetime.now()
#    #            
#    #            obsst_time_str = calinfo[eachkey].obsstart[i].split('T')[1]
#    #            obssthr  = int(obsst_time_str.split(':')[0])
#    #            obsstmin = int(obsst_time_str.split(':')[1])
#    #            obsstsec = int(obsst_time_str.split(':')[2])
#    #            obsst_sec_since_midnight = (now.replace(hour=0, minute=0, second=0, microsecond=0) - 
#    #                                       now.replace(hour=obssthr, minute=obsstmin, 
#    #                                                   second=obsstsec, microsecond=0)).total_seconds()
#    #            obssttimes.append(obsst_sec_since_midnight)
#    #         
#    #            obsen_time_str = calinfo[eachkey].obsend[i].split('T')[1]
#    #            obsenhr  = int(obsen_time_str.split(':')[0])
#    #            obsenmin = int(obsen_time_str.split(':')[1])
#    #            obsensec = int(obsen_time_str.split(':')[2])
#    #            obsen_sec_since_midnight = (now.replace(hour=0, minute=0, second=0, microsecond=0) - 
#    #                                       now.replace(hour=obsenhr, minute=obsenmin, 
#    #                                                   second=obsensec, microsecond=0)).total_seconds()
#    #            obsentimes.append(obsen_sec_since_midnight)
#    #            delta_transvals.append(calinfo[eachkey].delta_trans[i])
#    #
#    #    plt.scatter(delta_transvals,obssttimes,color='b')
#    #    plt.scatter(np.array(delta_transvals)[np.where(np.array(delta_transvals)>0.05)],
#    #                                                   np.array(obssttimes)[np.where(np.array(delta_transvals)>0.05)],
#    #                color='r')
#    #    plt.xlabel('\Delta_trans')
#    #    plt.ylabel('OBSSTART')
#    #    plt.show()
#    #
#    #for eachkey in calinfo.keys():
#    #    unphysical_ind = np.where(calinfo[eachkey].delta_trans>0.05)
#    #    plt.scatter(calinfo[eachkey].delta_trans,calinfo[eachkey].amstart,color='b')
#    #    plt.scatter(calinfo[eachkey].delta_trans[unphysical_ind],calinfo[eachkey].amstart[unphysical_ind],color='r')
#    #    plt.xlabel('\Delta Trans')
#    #    plt.ylabel('Airmass')
#    #    plt.show()
#    #
#    #
#    #for eachkey in calinfo.keys():
#    #    unphysical_ind = np.where(calinfo[eachkey].delta_trans>0.05)
#    #    plt.scatter(calinfo[eachkey].delta_trans,calinfo[eachkey].peak_fluxes,color='b')
#    #    plt.scatter(calinfo[eachkey].delta_trans[unphysical_ind],calinfo[eachkey].peak_fluxes[unphysical_ind],color='r')
#    #    plt.xlabel('\Delta Trans')
#    #    plt.ylabel('Peak Flux')
#    #    plt.show()
#    #
#    #
#    #for eachkey in calinfo.keys():
#    #    unphysical_ind = np.where(calinfo[eachkey].delta_trans>0.05)
#    #    plt.scatter(calinfo[eachkey].delta_trans,calinfo[eachkey].total_fluxes,color='b')
#    #    plt.scatter(calinfo[eachkey].delta_trans[unphysical_ind],calinfo[eachkey].total_fluxes[unphysical_ind],color='r')
#    #    plt.xlabel('\Delta Trans')
#    #    plt.ylabel('Total Flux')
#    #    plt.show()
#    
#    
#    
