def TauRelPipeline(source,OR_coeffs,wave,mindate,maxdate,aperture_diam=0.01666666667,physical_thresh=0.05):
    '''
    This pipeline relies on a variety of programs found in the 
    
    `TauRelPrepFunctions.py` and 
    `TauRelAnalysis_20171215.py` packages 
    
    (also `get_noext_reductions_from_kamaka.py`, coded with Graham Bell's help). 
    
    The first has functions that crop images, find files, etc while the second one has functions 
    which can read fits headers using .sdf file format and perform general analysis like running KAPPA's beamfit.

    This code relies on a consistent reduction being performed nightly. Currently, there is a kamaka code that takes 
    all the calibrator data and reduced it on 1 arcsecond pixels with NO EXTINCTION CORRECTION. The data is in pW and it has no extinction correction applied.
    These files can be found by running get_noext_reductions_from_kamaka.py.

    This is the order of operations for `TauRelPipeline.py`:

    1. Gather the reduced data that has no extinction correction or FCF factor applied
    2. Select, from those, the source you are currently interested in
    3. Construct the Transmission versus PWV function as before
    4. For each pair of physical coefficients (a,b) and each calibrator observation:
    
        a. Get the tau225 from the header by averaging WVMTAUST and WVMTAUEN (these WVM values may be wrong because we did not apply an EXT model to the data! Check with the raw WVM information)
        b. Get the airmass from the header, again by averaging over the short observation
        c. Calculate the transmission and compare it to the expected transmission
        d. Apply the new extinction correction (transmission) using CMULT
        e. Crop the image in order to run GaussClumps
        f. Fit the peak flux with GaussClumps and obtain the size & area of the calibrator
        g. Lay down a 1 arcminute diameter aperture centered on the calibrator with beamfit and measure the total flux

    5. Save all the information for every (a,b) pair (size, location, peak flux, total flux, etc.)
    6. Plot the Peak Flux versus Transmission for every calibrator observation and for every (a,b) extinction correction and fit a linear regression
    7. Plot the Total Flux versus Transmission for every calibrator observation and for every (a,b) extinction correction and fit a linear regression
    8. Find the shallowest slopes = the optimal (a,b)


    ######
    To call:

    TauRelPipeline('CRL618',[(26.0, 0.012),(26.5,0.012)],'450','20161101','25001231')

    source          = STR: name of source from fits header
    OR_coeffs       = LIST: a list of tuples (a,b) representing opacity relation coefficients you would like to test:

                      Im = I0 exp(-tau x airmass)

                      So to correct for the extinction, we apply a version of: exp(-tau x airmass)

                      This requires an opacity relation for the wavelengths we care about.
 
                                       tau_wave = a(tau_225 - b)
 
                      The current data reduction uses: 
   
                                       tau_850 = 4.6(tau_225 - 0.0043) 
                                       tau_450 = 26.0(tau_225 - 0.012)
 

    wave            = STR: Wavelength in microns
    mindate         = STR: Minimum date to consider: e.g. '20131231' -- WILL START ONE DAY AFTER
    maxdate         = STR: Maximum date to consider: e.g  '20131231' -- WILL STOP ONE DAY BEFORE
    aperture_diam   = FLOAT: Diameter of aperture to measure total flux
    physical_thresh = FLOAT: The percentage/100.0 that the calculated transmission can differ from the CSO model and still be considered physical

    '''

    print('\n##############\n##############\n##############')
    print('### Beginning ###')
    print('##############\n##############\n##############')
    
    import time
    start_time = time.time()
    import os
    from os import sys
    import subprocess
    import numpy as np
    import matplotlib.pyplot as plt
    import pickle
    import glob
    from starlink import kappa
    
    trans_minus_expectedtrans_thresh = physical_thresh #historical     
    coeffs_from_user = np.array(OR_coeffs)
    
    coeff1s_from_user = []
    coeff2s_from_user = []
    
    for eachcoeffpair in coeffs_from_user:
        coeff1s_from_user.append(eachcoeffpair[0])
        coeff2s_from_user.append(eachcoeffpair[1])
    
    ###############################################
    ###############################################
    ###############################################
    
    
    
    ###########
    ########### Section 1: Preparation
    ###########
    ###################
    # The files we are accessing
    # are uncalibrated and have no extinction
    # correction applied through the makemap
    # procedure - so lets build a list
    # of directories and files
    ###################
    ###########
    ###########
    ###########
    
    os.system('mkdir '+source)
    
    from TauRelPrepFunctions import crop_img
    from TauRelAnalysis_20171215 import getsdfhdr
    
    print('\n##############\n##############\n##############')
    print('### "Gathering list of NOEXT, UNCALIBRATED Files" ###')
    print('##############\n##############\n##############')
    
    getdir_command = "python scripts/get_noext_reductions_from_kamaka.py"
    process = subprocess.Popen(getdir_command.split(), stdout=subprocess.PIPE)
    alldirsbytes, error = process.communicate()
    alldirs = str(alldirsbytes).split('\\n')[0:-1]   
 
    reduced_cal_noext_files=[]
    
    if source+'_'+str(wave)+'_NOEXTfiles_'+str(mindate)+'.txt' not in os.listdir('FileLists/'):
        previous_file_exists = False
        noextfilelist = open('FileLists/'+source+'_'+str(wave)+'_NOEXTfiles_'+str(mindate)+'.txt','w')
        for eachdir in alldirs:
            if eachdir[0]=='b':
                directoryname = eachdir.split("b'")[-1]
            else:
                directoryname = eachdir
            for eachfile in os.listdir(directoryname):
                if eachfile.split('_')[-1]=='reduced.sdf':
                    if eachfile.split('_')[-2]==wave:
                        if np.logical_and(int(eachfile.split('_')[0].split('s')[-1])>=int(mindate),int(eachfile.split('_')[0].split('s')[-1])<int(maxdate)):
                            if kappa.fitsval(directoryname+'/'+eachfile,'OBJECT').value==source:
                                reduced_cal_noext_files.append(directoryname+'/'+eachfile)
                                noextfilelist.write(directoryname+'/'+eachfile+'\n')
        noextfilelist.close()
    else:
        previous_file_exists = True
        dates_of_files_added = []
        scans_of_files_added = []
        for eachfile in open('FileLists/'+source+'_'+str(wave)+'_NOEXTfiles_'+str(mindate)+'.txt','r'):
            print(eachfile,' already done!')
            dates_of_files_added.append(float(eachfile.split('/')[-1].split('_')[0].split('s')[-1]))
            scans_of_files_added.append(float(eachfile.split('/')[-1].split('_')[1]))
        if len(dates_of_files_added)>0:
            date_of_last_file_added = int(max(dates_of_files_added))
            scan_of_last_file_added = int(max(scans_of_files_added))
        else:
            date_of_last_file_added = int(mindate)
            scan_of_last_file_added = -1.0
        noextfilelist = open('FileLists/'+source+'_'+str(wave)+'_NOEXTfiles_'+str(mindate)+'.txt','a')
        for eachdir in alldirs:
            if eachdir[0]=='b':
                directoryname = eachdir.split("b'")[-1]
            else:
                directoryname = eachdir
            for eachfile in os.listdir(directoryname):
                if eachfile.split('_')[-1]=='reduced.sdf':
                    if eachfile.split('_')[-2]==wave:
                        if float(eachfile.split('_')[0].split('s')[-1]) > date_of_last_file_added:
                            if float(eachfile.split('_')[0].split('s')[-1]) < int(maxdate):
                                if kappa.fitsval(directoryname+'/'+eachfile,'OBJECT').value==source:
                                    reduced_cal_noext_files.append(directoryname+'/'+eachfile)
                                    noextfilelist.write(directoryname+'/'+eachfile+'\n')
                        elif float(eachfile.split('_')[0].split('s')[-1]) == date_of_last_file_added:
                            if float(eachfile.split('_')[0].split('s')[-1]) < int(maxdate):
                                if float(eachfile.split('_')[1]) > scan_of_last_file_added:
                                     if kappa.fitsval(directoryname+'/'+eachfile,'OBJECT').value==source:
                                         reduced_cal_noext_files.append(directoryname+'/'+eachfile)
                                         noextfilelist.write(directoryname+'/'+eachfile+'\n')
    
    print('\n##############\n##############\n##############')
    print('### "List of files has been gathered!" ###')
    print('##############\n##############\n##############')
    
    
    ############
    ############ Section 2: Analysis
    ############
    ####################
    # Now, let's define a range of
    # a and b values to step through
    # depending on the wavelength
    # so we can apply a new extinction
    # correction, fit the calibrator
    # and see which slope is the 
    # flattest on a Flux versus Transmission plot.
    ####################
    #############
    #############
    #############
    
    from TauRelAnalysis_20171215 import FitCSOTransvsPWVFast,CSOtrans,fitcal,fitcalBF
    from starlink.ndfpack import Ndf
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from starlink import convert
    import astropy.io.fits as apfits

    # Construct a functional form of Transmission vs PWV
    
    print('\n##############\n##############\n##############')
    print('### Building Transmission vs PWV function')
    print('##############\n##############\n##############')
    
    TransvsPWV = FitCSOTransvsPWVFast(int(wave))
    
    coeff1_range = coeff1s_from_user
    coeff2_range = coeff2s_from_user
    
    slopes            = []
    intercepts        = []
    coeff1s           = []
    coeff2s           = []
    trans_unphys_perc = []
    
    print('\n##############\n##############\n##############')
    print('### Testing Pairs of Coefficients for the Tau Relation')
    print('##############\n##############\n##############')
    
    results_dict = {}
    if previous_file_exists == True:    
        previous_results_dict = pickle.load(open(sorted(glob.glob('results/TauRelPipeline_FullResults_'+wave+'_'+source+'_'+str(mindate)+'*bin'))[-1],'rb'))
    
    paircount = 0
    
    numfiles = len(reduced_cal_noext_files)
    
    
    for eachpair in range(len(coeff1_range)):
        #for coeff2 in coeff2_range:
    
            paircount = paircount + 1
            print ('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n######\n# Pair '+str(paircount)+' out of '+str(len(coeff1_range))+'\n######\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
    
            run_num = paircount-1
    
            run_num_str = "Run_"+str(run_num)
    
            # Values to store for each file
            transphysical = [] 
            transmissions = []
            peak_fluxes   = []
            map_peak_fluxes = []
            BFPeak_fluxes = []
            total_fluxes  = []
            areas         = []
            FWHM1s        = []
            FWHM2s        = []
            fnames        = []
            WVMST         = []
            WVMEN         = []
            WVMST_TIME    = []
            WVMEN_TIME    = []
            AMSTART       = []
            AMEND         = []
            obsstart      = []
            obsend        = []
            ATSTART       = []
            ATEND         = []
            delta_trans_list  = [] #abs(trans - expected_trans)
            # New Params
            OBSNUM        = []
            UTDATE        = [] 
            AZSTART       = [] 
            AZEND         = [] 
            ELSTART       = [] 
            ELEND         = [] 
            HUMSTART      = [] 
            HUMEND        = [] 
            BPSTART       = [] 
            BPEND         = [] 
            WNDSPDST      = [] 
            WNDSPDEN      = [] 
            WNDDIRST      = [] 
            WNDDIREN      = [] 
            TAU225ST      = [] 
            TAU225EN      = [] 
            TAUDATST      = [] 
            TAUDATEN      = [] 
            SEEINGST      = [] 
            SEEINGEN      = [] 
            SEEDATST      = [] 
            SEEDATEN      = [] 
            FRLEGTST      = [] 
            FRLEGTEN      = [] 
            BKLEGTST      = [] 
            BKLEGTEN      = [] 
            MJD_OBS       = [] #In FITS HEADER - hyphen
            MJD_END       = [] #IN FITS HEADER - hyphen
            ELAPTIME      = [] 
            
    
            filecount     = 0 
    
            for eachfile in reduced_cal_noext_files:
                    eachfile = eachfile.split('\n')[0]
                    filecount = filecount + 1
                    print ('\n\n File: '+eachfile+' ('+str(filecount)+' out of '+str(numfiles)+')\n\n')
                    # Get the header information to calculate the Transmission 
                    convert.ndf2fits(eachfile,eachfile.split('/')[-1].split('.sdf')[0]+'.fits')
                    hdr = apfits.getheader(eachfile.split('/')[-1].split('.sdf')[0]+'.fits')
                    list_of_header_values = []
                    for eachcard in range(len(hdr.cards)):
                        list_of_header_values.append(hdr.cards[eachcard][0])
                    os.system('rm -f '+eachfile.split('/')[-1].split('.sdf')[0]+'.fits')
                    #hdr = getsdfhdr(eachfile) 
                    
    
                    tau225       = (hdr['WVMTAUST']+hdr['WVMTAUEN'])/2.0
    
                    tau_jcmt     = coeff1_range[eachpair]*(tau225-coeff2_range[eachpair])
    
                    airmass      = (hdr['AMSTART']+hdr['AMEND'])/2.0
    
                    #transmission = np.exp(-1.0*tau_jcmt*airmass)
    
                    # The Tau225 value is the zenith value! So don't modify it by the arimass!!
    		# We are comparing it to the zenith transmission!
    
                    transmission = np.exp(-1.0*tau_jcmt)
    
                    expected_transmission = CSOtrans(int(wave),tau225,TransvsPWV)
    
                    delta_trans = abs(transmission - expected_transmission)
                    delta_trans_list.append(delta_trans)
    
                    if delta_trans > trans_minus_expectedtrans_thresh:
                        
                        transphysical.append('no')
    
                    else:
    
                        transphysical.append('yes')
    
                    # Apply the new extinction correction to the file - YOU NEED THE AIRMASS HERE - THE ORIGINAL EXT CORRECTION INCLUDES THE AIRMASS
                    
                    new_ext_cor    = np.exp(-1.0*tau_jcmt*airmass) 
    
                    extcor_command = '$KAPPA_DIR/cdiv in='+eachfile
                    extcor_command += ' scalar='+str(new_ext_cor)+' out='+source+'/'+eachfile.split('/')[-1].split('.sdf')[0]+'_extcor.sdf'
                    subprocess.call(extcor_command, shell=True)
    
                    #Crop the image!
                    crop_img(source+'/'+eachfile.split('/')[-1].split('.sdf')[0]+'_extcor.sdf','./',CROP_METHOD='CIRCLE',MAP_RADIUS=200)
    
                    # Now measure the peak and extended structure of the calibrator with Gaussclumps, return these values
                    # peak_flux,total_flux, FWHM1, FWHM2 -- need to compare Beamfit to Gaussclumps - benefit of GC: Total and Peak at the same time
    		    # so run gaussclumps, grab the brightest source - and output the info - also tells you FWHM_maj and FWHM_min
                    peak_flux,total_fluxbad,FWHM1,FWHM2,area,peakX,peakY = fitcal(eachfile.split('/')[-1].split('.sdf')[0]+'_extcor_crop.sdf',"GCParms/GCParmsm.txt")
                    BFPeak,BFPeakunc,BFTotal,BFmajFWHM,BFmajFWHMunc,BFminFWHM,BFminFWHMunc = fitcalBF(eachfile.split('/')[-1].split('.sdf')[0]+'_extcor_crop.sdf',hdr,source)
                    if not np.isnan(peak_flux): #Sometimes Gaussclumps does not fit a calibrator- testing with beamfit currently 20170920
   
                        if float(peakX.split('h ')[1].split('m ')[1].split('s')[0])>=60:
                            peakX                                                = peakX.split('h ')[0]+':'+str(int(float(peakX.split('h ')[1].split('m ')[0])+1))+':0'+"{:1.2f}".format(float(peakX.split('h ')[1].split('m ')[1].split('s')[0]) % 60)
                        else:
                            peakX                                                = peakX.split('h ')[0]+':'+peakX.split('h ')[1].split('m ')[0]+':'+peakX.split('h ')[1].split('m ')[1].split('s')[0]
                        coord = SkyCoord(peakX,peakY,unit=(u.hourangle, u.deg))
                        ndf = Ndf(eachfile.split('/')[-1].split('.sdf')[0]+'_extcor_crop.sdf')
                        get_pixel_coord = ndf.wcs.tran([[coord.ra.radian], [coord.dec.radian],[float(wave) / 1000000]],False)
                        (xgrid, ygrid, zgrid) = (get_pixel_coord - 1).flatten().round().astype(int)
                        MAPPeak = ndf.data[0, ygrid, xgrid].item()
                        total_flux                                           = kappa.aperadd(eachfile.split('/')[-1].split('.sdf')[0]+'_extcor_crop.sdf',centre='"'+peakX+','+peakY+'"',diam=aperture_diam).total
                        annulus_ARD                                          = open('annulus.ard','w')
                        annulus_ARD.write('COFRAME(SKY,SYSTEM=FK5,EQUINOX=2000)')

##################################
                        # THE FOLLOWING LINE IS INCORRECT - 20180423 - 2 things wrong here, the background should be calculated between 90 and 120 arcsecond diameter apertures. These are 60 and 90 arcseconds and in annulus_ARD, they should be given in readius! Not Diameter! See the fix one line below.
                        #annulus_ARD.write('CIRCLE('+peakX+','+peakY+',0.025) .AND. .NOT. CIRCLE('+peakX+','+peakY+',0.01666666)') 
                        #
                        # THE FOLLOWING LINE IS CORRECT   - fixed at 12:00pm on 20180423 - everything produced since then will have this fix
                        annulus_ARD.write('CIRCLE('+peakX+','+peakY+',0.016666) .AND. .NOT. CIRCLE('+peakX+','+peakY+',0.0125)')
###################################
                        annulus_ARD.close()
                        kappa.ardmask(eachfile.split('/')[-1].split('.sdf')[0]+'_extcor_crop.sdf','annulus.ard',out='ardmask_bkgnd',INSIDE=False)
                        total_flux_bckgnd                                = kappa.stats('ardmask_bkgnd').mean
                        os.system('rm -f annulus.ard ardmask_bkgnd*')
    
                        peak_fluxes.append(peak_flux)
                        map_peak_fluxes.append(MAPPeak)
                        BFPeak_fluxes.append(BFPeak)
                        total_fluxes.append(total_flux-total_flux_bckgnd)
    
                        areas.append(area)
                        FWHM1s.append(FWHM1)
                        FWHM2s.append(FWHM2)
                        transmissions.append(transmission) # THIS IS THE ZENITH TRANSMISSION - NO AIRMASS
                        fnames.append(source+'/'+eachfile.split('/')[-1].split('.sdf')[0]+'_extcor.sdf')
                        WVMST.append(hdr['WVMTAUST'])
                        WVMEN.append(hdr['WVMTAUEN'])
                        WVMST_TIME.append(hdr['WVMDATST'])
                        WVMEN_TIME.append(hdr['WVMDATEN'])
                        AMSTART.append(hdr['AMSTART'])
                        AMEND.append(hdr['AMEND'])
                        obsstart.append(hdr['DATE-OBS'])
                        obsend.append(hdr['DATE-END'])
                        ATSTART.append(hdr['ATSTART'])
                        ATEND.append(hdr['ATEND'])
                        OBSNUM.append(hdr['OBSNUM']) 
                        UTDATE.append(hdr['UTDATE']) 
                        AZSTART.append(hdr['AZSTART']) 
                        AZEND.append(hdr['AZEND']) 
                        ELSTART.append(hdr['ELSTART']) 
                        ELEND.append(hdr['ELEND']) 
                        HUMSTART.append(hdr['HUMSTART']) 
                        HUMEND.append(hdr['HUMEND']) 
                        BPSTART.append(hdr['BPSTART']) 
                        BPEND.append(hdr['BPEND']) 
                        WNDSPDST.append(hdr['WNDSPDST']) 
                        WNDSPDEN.append(hdr['WNDSPDEN']) 
                        WNDDIRST.append(hdr['WNDDIRST']) 
                        WNDDIREN.append(hdr['WNDDIREN']) 
                        TAU225ST.append(hdr['TAU225ST']) 
                        TAU225EN.append(hdr['TAU225EN']) 
                        TAUDATST.append(hdr['TAUDATST']) 
                        TAUDATEN.append(hdr['TAUDATEN']) 
                        if 'SEEINGST' in list_of_header_values:
                            SEEINGST.append(hdr['SEEINGST'])
                        else:
                            SEEINGST.append(np.nan)
                        if 'SEEINGEN' in list_of_header_values:
                            SEEINGEN.append(hdr['SEEINGEN']) 
                        else:
                            SEEINGEN.append(np.nan)
                        if 'SEEDATST' in list_of_header_values:
                            SEEDATST.append(hdr['SEEDATST']) 
                        else:
                            SEEDATST.append(np.nan)
                        if 'SEEDATEN' in list_of_header_values:
                            SEEDATEN.append(hdr['SEEDATEN']) 
                        else:
                            SEEDATEN.append(np.nan)
                        FRLEGTST.append(hdr['FRLEGTST']) 
                        FRLEGTEN.append(hdr['FRLEGTEN']) 
                        BKLEGTST.append(hdr['BKLEGTST']) 
                        BKLEGTEN.append(hdr['BKLEGTEN']) 
                        MJD_OBS.append(hdr['MJD-OBS']) 
                        MJD_END.append(hdr['MJD-END']) 
                        ELAPTIME.append(hdr['ELAPTIME']) 


                    os.system('rm -f '+eachfile.split('.sdf')[0]+'_extcor_crop.sdf')
                    os.system('rm -f '+source+'/'+eachfile.split('/')[-1].split('.sdf')[0]+'_extcor*.sdf')
    
            if 'no' in transphysical:
                print ('\n For this pairing of coefficients:\n')
                print ('a = '+str(coeff1_range[eachpair]))
                print ('b = '+str(coeff2_range[eachpair]))
                print ('\nTranmission Unphysical for '+str(100*float(len(np.where(np.array(transphysical)=='no')[0]))/float(len(np.array(transphysical))))+'% of measurements \n')
                trans_unphys_perc.append(float(100*len(np.where(np.array(transphysical)=='no')[0]))/float(len(np.array(transphysical))))
    
            else:
                trans_unphys_perc.append(0.0)
    
            coeff1s.append(coeff1_range[eachpair])
            coeff2s.append(coeff2_range[eachpair])
            if len(peak_fluxes)>1:
                m,b = np.polyfit(transmissions,peak_fluxes,1)
                slopes.append(m)
                intercepts.append(b)
            else:
                slopes.append(np.nan)
                intercepts.append(np.nan)

            if previous_file_exists == True:
                 previous_results_dict[run_num_str]['fnames'].extend(list(fnames))
                 previous_results_dict[run_num_str]['peak_fluxes'].extend(list(peak_fluxes))
                 previous_results_dict[run_num_str]['BFPeak'].extend(list(BFPeak_fluxes))
                 previous_results_dict[run_num_str]['MapPeak'].extend(list(map_peak_fluxes))
                 previous_results_dict[run_num_str]['total_fluxes'].extend(list(total_fluxes))
                 previous_results_dict[run_num_str]['areas'].extend(list(areas))
                 previous_results_dict[run_num_str]['FWHM1s'].extend(list(FWHM1s))
                 previous_results_dict[run_num_str]['FWHM2s'].extend(list(FWHM2s))
                 previous_results_dict[run_num_str]['transmissions'].extend(list(transmissions))
                 previous_results_dict[run_num_str]['transmissions'].extend(list(delta_trans_list))
                 previous_results_dict[run_num_str]['WVMTAUST'].extend(list(WVMST))
                 previous_results_dict[run_num_str]['WVMTAUEN'].extend(list(WVMEN))
                 previous_results_dict[run_num_str]['WVMTAUST_TIME'].extend(list(WVMST_TIME))
                 previous_results_dict[run_num_str]['WVMTAUEN_TIME'].extend(list(WVMEN_TIME))
                 previous_results_dict[run_num_str]['AMSTART'].extend(list(AMSTART))
                 previous_results_dict[run_num_str]['AMEND'].extend(list(AMEND))
                 previous_results_dict[run_num_str]['OBSSTART'].extend(list(obsstart))
                 previous_results_dict[run_num_str]['OBSEND'].extend(list(obsend))
                 previous_results_dict[run_num_str]['ATSTART'].extend(list(ATSTART))
                 previous_results_dict[run_num_str]['ATEND'].extend(list(ATEND))
                 previous_results_dict[run_num_str]['OBSNUM'].extend(list(OBSNUM))
                 previous_results_dict[run_num_str]['UTDATE'].extend(list(UTDATE))
                 previous_results_dict[run_num_str]['AZSTART'].extend(list(AZSTART))
                 previous_results_dict[run_num_str]['AZEND'].extend(list(AZEND))
                 previous_results_dict[run_num_str]['ELSTART'].extend(list(ELSTART))
                 previous_results_dict[run_num_str]['ELEND'].extend(list(ELEND))
                 previous_results_dict[run_num_str]['HUMSTART'].extend(list(HUMSTART))
                 previous_results_dict[run_num_str]['HUMEND'].extend(list(HUMEND))
                 previous_results_dict[run_num_str]['BPSTART'].extend(list(BPSTART))
                 previous_results_dict[run_num_str]['BPEND'].extend(list(BPEND))
                 previous_results_dict[run_num_str]['WNDSPDST'].extend(list(WNDSPDST))
                 previous_results_dict[run_num_str]['WNDSPDEN'].extend(list(WNDSPDEN))
                 previous_results_dict[run_num_str]['WNDDIRST'].extend(list(WNDDIRST))
                 previous_results_dict[run_num_str]['WNDDIREN'].extend(list(WNDDIREN))
                 previous_results_dict[run_num_str]['TAU225ST'].extend(list(TAU225ST))
                 previous_results_dict[run_num_str]['TAU225EN'].extend(list(TAU225EN))
                 previous_results_dict[run_num_str]['TAUDATST'].extend(list(TAUDATST))
                 previous_results_dict[run_num_str]['TAUDATEN'].extend(list(TAUDATEN))
                 previous_results_dict[run_num_str]['SEEINGST'].extend(list(SEEINGST))
                 previous_results_dict[run_num_str]['SEEINGEN'].extend(list(SEEINGEN))
                 previous_results_dict[run_num_str]['SEEDATST'].extend(list(SEEDATST))
                 previous_results_dict[run_num_str]['SEEDATEN'].extend(list(SEEDATEN))
                 previous_results_dict[run_num_str]['FRLEGTST'].extend(list(FRLEGTST))
                 previous_results_dict[run_num_str]['FRLEGTEN'].extend(list(FRLEGTEN))
                 previous_results_dict[run_num_str]['BKLEGTST'].extend(list(BKLEGTST))
                 previous_results_dict[run_num_str]['BKLEGTEN'].extend(list(BKLEGTEN))
                 previous_results_dict[run_num_str]['MJD-OBS'].extend(list(MJD_OBS))
                 previous_results_dict[run_num_str]['MJD-END'].extend(list(MJD_END))
                 previous_results_dict[run_num_str]['ELAPTIME'].extend(list(ELAPTIME))

            results_dict[run_num_str]={}
            results_dict[run_num_str]['fnames']                = fnames
            results_dict[run_num_str]['coeff1']                = coeff1_range[eachpair]
            results_dict[run_num_str]['coeff2']                = coeff2_range[eachpair]
            results_dict[run_num_str]['source']                = source
            results_dict[run_num_str]['peak_fluxes']           = peak_fluxes
            results_dict[run_num_str]['BFPeak']                = BFPeak_fluxes
            results_dict[run_num_str]['MapPeak']               = map_peak_fluxes
            results_dict[run_num_str]['total_fluxes']          = total_fluxes
            results_dict[run_num_str]['areas']                 = areas
            results_dict[run_num_str]['FWHM1s']                = FWHM1s
            results_dict[run_num_str]['FWHM2s']                = FWHM2s
            results_dict[run_num_str]['averageFWHM1']          = np.average(FWHM1s)
            results_dict[run_num_str]['averageFWHM2']          = np.average(FWHM2s)
            results_dict[run_num_str]['transmissions']         = transmissions
            results_dict[run_num_str]['delta_trans']           = delta_trans_list
            results_dict[run_num_str]['WVMTAUST']              = WVMST
            results_dict[run_num_str]['WVMTAUEN']              = WVMEN
            results_dict[run_num_str]['WVMTAUST_TIME']         = WVMST_TIME
            results_dict[run_num_str]['WVMTAUEN_TIME']         = WVMEN_TIME
            results_dict[run_num_str]['AMSTART']               = AMSTART
            results_dict[run_num_str]['AMEND']                 = AMEND
            results_dict[run_num_str]['OBSSTART']              = obsstart
            results_dict[run_num_str]['OBSEND']                = obsend
            results_dict[run_num_str]['ATSTART']               = ATSTART
            results_dict[run_num_str]['ATEND']                 = ATEND
            results_dict[run_num_str]['OBSNUM']                = OBSNUM
            results_dict[run_num_str]['UTDATE']                = UTDATE
            results_dict[run_num_str]['AZSTART']               = AZSTART
            results_dict[run_num_str]['AZEND']                 = AZEND
            results_dict[run_num_str]['ELSTART']               = ELSTART
            results_dict[run_num_str]['ELEND']                 = ELEND
            results_dict[run_num_str]['HUMSTART']              = HUMSTART
            results_dict[run_num_str]['HUMEND']                = HUMEND
            results_dict[run_num_str]['BPSTART']               = BPSTART
            results_dict[run_num_str]['BPEND']                 = BPEND
            results_dict[run_num_str]['WNDSPDST']              = WNDSPDST
            results_dict[run_num_str]['WNDSPDEN']              = WNDSPDEN
            results_dict[run_num_str]['WNDDIRST']              = WNDDIRST
            results_dict[run_num_str]['WNDDIREN']              = WNDDIREN
            results_dict[run_num_str]['TAU225ST']              = TAU225ST
            results_dict[run_num_str]['TAU225EN']              = TAU225EN
            results_dict[run_num_str]['TAUDATST']              = TAUDATST
            results_dict[run_num_str]['TAUDATEN']              = TAUDATEN
            results_dict[run_num_str]['SEEINGST']              = SEEINGST
            results_dict[run_num_str]['SEEINGEN']              = SEEINGEN
            results_dict[run_num_str]['SEEDATST']              = SEEDATST
            results_dict[run_num_str]['SEEDATEN']              = SEEDATEN
            results_dict[run_num_str]['FRLEGTST']              = FRLEGTST
            results_dict[run_num_str]['FRLEGTEN']              = FRLEGTEN
            results_dict[run_num_str]['BKLEGTST']              = BKLEGTST
            results_dict[run_num_str]['BKLEGTEN']              = BKLEGTEN
            results_dict[run_num_str]['MJD-OBS']               = MJD_OBS
            results_dict[run_num_str]['MJD-END']               = MJD_END
            results_dict[run_num_str]['ELAPTIME']              = ELAPTIME

            #if 'no' in transphysical:
            #    results_dict[run_num_str]['trans_unphys_per'] = float(100*len(np.where(np.array(transphysical)=='no')[0]))/float(len(np.array(transphysical)))
            #else:
            #    results_dict[run_num_str]['trans_unphys_per'] = 0.0 
    
    
    #plt.clf()
    #plt.hist(delta_trans_list)
    #plt.savefig('deltatrans_hist_'+str(time.localtime().tm_year)+str(time.localtime().tm_mon)+str(time.localtime().tm_mday)+str(time.localtime().tm_hour)+str(time.localtime().tm_min)+'.pdf',format='pdf')
    #plt.clf()
    
    coeff1s           = np.array(coeff1s)
    coeff2s           = np.array(coeff2s)
    slopes            = np.array(slopes)
    trans_unphys_perc = np.array(trans_unphys_perc)
    
    #plt.hist(trans_unphys_perc)
    #plt.savefig('trans_unphys_perc_'+str(time.localtime().tm_year)+str(time.localtime().tm_mon)+str(time.localtime().tm_mday)+str(time.localtime().tm_hour)+str(time.localtime().tm_min)+'.pdf',format='pdf')
    #plt.clf()
    
    #plt.scatter(np.arange(0,len(coeff1s),1),trans_unphys_perc)
    #plt.savefig('trans_unphys_perc_scatter_'+str(time.localtime().tm_year)+str(time.localtime().tm_mon)+str(time.localtime().tm_mday)+str(time.localtime().tm_hour)+str(time.localtime().tm_min)+'.pdf',format='pdf')
    #plt.clf()
    
    # Save the results dictionary
    YearForFilename  = str(time.localtime().tm_year) 
    if float(time.localtime().tm_mon)<10:
        MonthForFilename = '0'+str(time.localtime().tm_mon)
    else:
        MonthForFilename = str(time.localtime().tm_mon)
    if float(time.localtime().tm_mday)<10:
        DayForFilename   = '0'+str(time.localtime().tm_mday)
    else:
        DayForFilename   = str(time.localtime().tm_mday)
    if float(time.localtime().tm_hour)< 10:
        HourForFilename  = '0'+str(time.localtime().tm_hour)
    else:
        HourForFilename  = str(time.localtime().tm_hour)
    if float(time.localtime().tm_min)<10:
        MinForFilename   = '0'+str(time.localtime().tm_min)
    else:
        MinForFilename   = str(time.localtime().tm_min)

    pickle.dump(results_dict,open('TauRelPipeline_FullResults_'+wave+'_'+source+'_'+mindate+'_'+YearForFilename+MonthForFilename+DayForFilename+HourForFilename+MinForFilename+".bin",'wb'))

    if not os.path.exists('results'): os.system('mkdir results')
    if not os.path.exists('Figures'): os.system('mkdir Figures')
    os.system('rm -f *crop*')
    os.system('mv *FullResults*bin results/')
    #clean up
    os.system('rm -f disp.dat log.group rules.badobs s201*sdf')
    os.system('mv *pdf Figures/')
    
    
    #print('\n##############\n##############\n##############')
    #print('### Finding the best Tau Relation!')
    #print('##############\n##############\n##############')
    
    
    # Now find the best ind:
    unphysical = 1
    while unphysical == 1:
        best_ind              = np.argmin(abs(slopes))
        precentage_unphysical = trans_unphys_perc[best_ind]
        if precentage_unphysical>40:
            coeff1s           = np.delete(coeff1s,best_ind)
            coeff2s           = np.delete(coeff2s,best_ind)
            slopes            = np.delete(slopes,best_ind)
            trans_unphys_perc = np.delete(trans_unphys_perc,best_ind)
        else:
            unphysical = 0
    
    #print ('\n\nBest Coefficients: \n\n')
    #print ('a = '+str(coeff1s[best_ind]))
    #print ('b = '+str(coeff2s[best_ind]))
    #print ('\ntau_'+wave+' = '+str(round(coeff1s[best_ind],2))+' x (tau_225 - '+str(round(coeff2s[best_ind],5))+')\n\n')
    #print ('Percentage of points that have unphysical transmissions: '+str(trans_unphys_perc[best_ind])+'%\n\n')
    
    # Write out the results
    #results_file = open("TauRelPipeline_Results_"+str(time.localtime().tm_year)+str(time.localtime().tm_mon)+str(time.localtime().tm_mday)+str(time.localtime().tm_hour)+str(time.localtime().tm_min)+".txt","w")
    #results_file.write('\n\nBest Coefficients: \n\n')
    #results_file.write('a = '+str(coeff1s[best_ind])+'\n')
    #results_file.write('b = '+str(coeff2s[best_ind]))
    #results_file.write('\ntau_'+wave+' = '+str(round(coeff1s[best_ind],2))+' x (tau_225 - '+str(round(coeff2s[best_ind],5))+')\n\n')
    #results_file.write('Slope and intercept of Peak Flxu versus Transmission plot: \n')
    #results_file.write('Percentage of points that have unphysical transmissions: '+str(trans_unphys_perc[best_ind])+'%\n\n')
    #results_file.write('\n\n'+str(((time.time() - start_time)/60.0))+' minutes to complete.\n\n')
    #results_file.close()
    
    os.system('mv *Results*txt results/')
    os.system('rm -f '+source+'/*')
    os.system('rmdir '+source+'/')

    print ('\n\n'+str(((time.time() - start_time)/60.0))+' minutes to complete.\n\n')
    
