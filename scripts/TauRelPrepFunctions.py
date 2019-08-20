# This Library contains all of the necessary preparatory functions
# to determine the Tau Relations based on SCUBA-2 calibrators.
#
# Here, I have functions to copy all of the reduced calibrator files to a local spot
# (all reduced with "bright_compact" and 1 arcsecond pixels) and uncalibrate the data
# and remove the extinction correction, 
#

# check Project_description.txt for a more detailed explanation.


############################
############################
# First, check if we have all the necessary reduced 
# calibrator files stored locally, and if not - copy them over
############################
############################

def copy_calreducs(orig_reduced_cal_dir,local_reduced_cal_dir,wave):

    # Import the OS package to perform system operations

    import os

    # Define the main directory from where we will be copying the calstats logs

    calstats_main_dir = orig_reduced_cal_dir

    # For each date in the directory, above ...

    date_direcs       = os.listdir(calstats_main_dir)

    print ('\n\nChecking to see which reduced calibrator files need to be copied to the local drive...\n\n')

    list_of_calibrator_reduced_files = []

    for i in date_direcs:

        # if the reduced data file isn't already in the reduced directory...

        for eachfile in os.listdir(calstats_main_dir+i):

            if eachfile.split('.')[-1]=='sdf':

                if len(eachfile.split('_'))>2:

                    if eachfile.split('_')[-2]==wave:

                        if eachfile.split('_')[-1]=='reduced.sdf':

                            list_of_calibrator_reduced_files.append(eachfile)
        
    for eachfile in list_of_calibrator_reduced_files:
 
            if os.path.isfile(local_reduced_cal_dir+'/'+eachfile):

                check='good'

            else:

                date = eachfile.split('_')[-4].split('s')[1]
                # ...Copy the reduced calibrator data file over

                print ('\n\tCopying '+eachfile+' over to '+local_reduced_cal_dir+'\n')

                os.system('cp '+calstats_main_dir+date+'/*_850_reduced.sdf '+local_reduced_cal_dir)

############################
############################
#
# Then sort all of files copied over
# by source
#
############################
############################


def sort_direc_tree(parent_directory,calinfo):

    # parent_directory = 'reduced/'
    # calinfo          = 'sc2cal_info.txt'

    import os                                                                                                                                                                                                       
    import numpy as np                                                                                                                                                                                              

    ############
    # Read in the calibrator data! Note this was all
    # reduced with 1 arcsecond pixels - I will
    # very much want to redo this with 3" or 4"
    # at 850
    ############

    master_calinfo = np.loadtxt(calinfo,ndmin=1,dtype={'names':('UT','HST','obs','Source','Mode','Filter',
                                                                    'El','Airmass','Trans','Tau225','Tau','Flux_ap',
                                                                    'Err_ap','Noise','Peak_obs','Peak_fit','FCFasec',
                                                                    'FCFasec_err','FCFbeam','FCFbeam_err','FCFmatch',
                                                                    'FCFmatch_err','FWHMmain','Error_beam','ID','obsid'
                                                                    ,'obsstatus','project'),
                                                           'formats':('<U16','<U19','i4','<U10','S5','i4','i4','f4',
                                                                      'f4','f4','f4','f4','f4','f4','f4','f4','f4',
                                                                      'f4','f4','f4','f4','f4','f4','f4','f4',
                                                                      '<U27','i4','<U7')})


    # Get a short list of all targets found within the sc2cal database:

    sources_used_as_cals = []

    for eachtarget in master_calinfo['Source']:
        if eachtarget not in sources_used_as_cals:
            sources_used_as_cals.append(eachtarget)

    # Then, associate which source belongs to which date/scan combination

    for eachuniqtarget in sources_used_as_cals:
        # make a directory for the source:
        os.system('mkdir '+parent_directory+eachuniqtarget)
        source_ind                  = np.where(master_calinfo['Source'] == eachuniqtarget)
        date_associated_with_source = master_calinfo['UT'][source_ind]
        scan_associated_with_source = master_calinfo['obs'][source_ind]
        for eachobs in range(len(date_associated_with_source)):

            date_str_identifier = date_associated_with_source[eachobs].split('-')[0]+date_associated_with_source[eachobs].split('-')[1]+date_associated_with_source[eachobs].split('-')[2]
            if len(str(scan_associated_with_source[eachobs]))==5:
                scan_str_identifier = str(scan_associated_with_source[eachobs])
            if len(str(scan_associated_with_source[eachobs]))==4:
                scan_str_identifier = '0'+str(scan_associated_with_source[eachobs])
            if len(str(scan_associated_with_source[eachobs]))==3:
                scan_str_identifier = '00'+str(scan_associated_with_source[eachobs])
            if len(str(scan_associated_with_source[eachobs]))==2:
                scan_str_identifier = '000'+str(scan_associated_with_source[eachobs])
            if len(str(scan_associated_with_source[eachobs]))==1:
                scan_str_identifier = '0000'+str(scan_associated_with_source[eachobs])
            str_identifier=date_str_identifier+'_'+scan_str_identifier
            os.system('mv '+parent_directory+'*'+str_identifier+'*.sdf '+parent_directory+eachuniqtarget)



###############################
###############################
# Now, uncalibrate the data and remove the extinction correction
###############################
###############################

# the default a and b values are for 850 microns!
# for 450 microns: a = 26.0 and b = 0.012

def uncal(fname,local_reduced_cal_dir,a=4.6,b=0.0043):

    import subprocess
    import pdb
    import os

    ##### Image Preparation Functions #####

    # Perform the "uncalibration". This divides by the FCF used
    # and corrects the units.

    print ('\nPerforming the "Uncalibration" back to pW... File: '+fname+'\n')

    if os.path.isfile(local_reduced_cal_dir+fname.split('.sdf')[0]+'_uncal.sdf'):
        
        print ('\tAlready Uncalibrated')

    else:

        uncal_command = '${ORAC_DIR}/etc/picard_start.sh -log f -nodisplay UNCALIBRATE_SCUBA2_DATA '+local_reduced_cal_dir+fname

        print (uncal_command)

        subprocess.call(uncal_command, shell=True)

        # Perform some cleanup from running picard

        if os.path.isfile(fname.split('.sdf')[0]+'_uncal.sdf'):
            os.system('mv *uncal.sdf '+local_reduced_cal_dir)
        else:
            os.system('cp '+local_reduced_cal_dir+fname+' '+local_reduced_cal_dir+fname.split('.sdf')[0]+'_uncal.sdf')
        os.system('rm -f log.group disp.dat rules.badobs s201*sdf')

        print('\tUNCALIBRATE = DONE')

    print ('\nBeginning the Extinction Correction removal...\n')

    if os.path.isfile(local_reduced_cal_dir+fname.split('.sdf')[0]+'_uncalNOext.sdf'):

        print ('\tAlready undergone extinction removal')

    else:

        # Next, remove the extinction correction. Recall: Im = I0 exp(-Tau*Airmass)
        # So, the extinction correction originally DIVIDED the data by exp(-Tau*Airmass)
        # Therefore, we need to MULTIPLY the data by import subprocess exp(-Tau*Airmass)
        # Note that Tau = a(TAU225 - b), where a is 4.6 and b is 0.0043 ny default
    
        # First, find the average TAU225 measurement

        # This requires reading in a FITS header from an NDF file and we have a sneaky way of doing this

        # First, import the Ndf module from starlink.ndfpack

        from starlink.ndfpack import Ndf

        # Then load in the NDF file
        myndf    = Ndf(local_reduced_cal_dir+fname)
       
        # Make a line break separated string of the header
        myndf_string_list = []
        for i in myndf.head['FITS']:
            myndf_string_list.append(i.decode('UTF-8'))

        # myndf.head['FITS'] is a fits header that is given as a list of strings
        # so split up the strings into values and comments and load it into a
        # standard astropy fits header

        from astropy.io import fits

        hdr = fits.Header.fromstring('\n'.join(myndf_string_list), sep='\n')

        # Now get the TAU_850 or TAU_450  measurement: use the average between the START and END of
        # the observation

        tau225 = (hdr['WVMTAUST']+hdr['WVMTAUEN'])/2.0

        # Then, get the airmass: again use the average from the start and the end
        # of the observation
    
        airmass = (hdr['AMSTART']+hdr['AMEND'])/2.0

        # Compute the extinction correction - the default a and b values are for 850 microns!

        import numpy as np

        ext_cor = np.exp(-1.0*a*(tau225-b)*airmass)

        # Perform the removal of the extinction correction

        extcor_command = '$KAPPA_DIR/cmult in='+local_reduced_cal_dir+fname.split('.sdf')[0]+'_uncal.sdf'
        extcor_command += ' scalar='+str(ext_cor)+' out='+local_reduced_cal_dir+fname.split('.sdf')[0]
        extcor_command += '_uncalNOext.sdf'
        subprocess.call(extcor_command, shell=True)
        
        print ('Extinction Correction = Removed')

def crop_img(fname,local_reduced_cal_dir,CROP_METHOD='CIRCLE',MAP_RADIUS=200):
    #Make the cropping parameter file.

    import subprocess
    import os
    
    print ('\nBeginning the Cropping Procedure'+' Method = '+CROP_METHOD+' Radius = '+str(MAP_RADIUS)+'...\n')

    if os.path.isfile(local_reduced_cal_dir+fname[:-4]+'_crop.sdf'):

        print('\tFile already Cropped!')

    else:

        crop_parms = open('crop.ini', 'w')
        crop_parms.write('[CROP_SCUBA2_IMAGES]\n')
        crop_parms.write('CROP_METHOD = '+str(CROP_METHOD)+'\n')
        crop_parms.write('MAP_RADIUS = '+str(MAP_RADIUS)+'\n')
        crop_parms.close()

        #Perform the cropping.
        crop_command = '${ORAC_DIR}/etc/picard_start.sh CROP_SCUBA2_IMAGES '
        crop_command += '-log f -recpars crop.ini '+local_reduced_cal_dir+fname+' ${1+"$@"};'
        subprocess.call(crop_command, shell=True)
        os.system('mv '+fname[:-4]+'_crop.sdf '+local_reduced_cal_dir)

        print('\tCrop = Done')
