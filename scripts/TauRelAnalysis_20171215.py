########################
########################
#
# These functions are for the
# analysis of the calibrators
# after their extinction cor.
# has been removed and a new
# one has been supplied.
#
########################
########################

######
#
# This function
# extracts a header
# from an sdf file
#
######

def getsdfhdr(fname):
    # First, import the Ndf module from starlink.ndfpack

    from starlink.ndfpack import Ndf

    # Then load in the NDF file
    myndf    = Ndf(fname)

    # Make a line break separated string of the header
    myndf_string_list = []
    for i in myndf.head['FITS']:
        myndf_string_list.append(i.decode('UTF-8'))

    # myndf.head['FITS'] is a fits header that is given as a list of strings
    # so split up the strings into values and comments and load it into a
    # standard astropy fits header

    from astropy.io import fits

    hdr = fits.Header.fromstring('\n'.join(myndf_string_list), sep='\n')
 
    return (hdr)

####################                                                                                                                                                                                                      
# This program takes the                                                                                                                                                                                                  
# Maunakea Transmission                                                                                                                                                                                                   
# data from the CSO interactive                                                                                                                                                                                           
# plotter and interpolates                                                                                                                                                                                                
# over it. That way, we make a                                                                                                                                                                                            
# function of Transmission versus PWV.                                                                                                                                                                                    
# Then, we can                                                                                                                                                                                                            
# test to see if our transmission                                                                                                                                                                                         
# values at 450 and 850 microns                                                                                                                                                                                           
# are physical                                                                                                                                                                                                            
####################                                                                                                                                                                                                      

def FitCSOTransvsPWV(wave,MakeDict='no'): # Wave is: 850 or 450

    import numpy as np
    import scipy.interpolate
    import matplotlib.pyplot as plt
    import pickle

    # Calculate the frequency based on the supplied wavelength

    c       = 299792459 # m/s
    freq_in = ( c / (wave*1e-6) ) / 1e9 # GHz

    # The output of this program is a transmission versus PWV curve

    pwvs  = []
    trans = []

    # Read in the CSO data at a variety of PWVs
    # and save the transmission information that comes out

    transfile='CSOtrans/CSO_05mm_'+str(wave)+'.txt'
    cso_data = np.loadtxt(transfile,delimiter='    ',usecols=(1,3),dtype={'names':('freq','trans'),'formats':('f4','f4')})
    interp_cso_data = scipy.interpolate.interp1d(cso_data['freq'],cso_data['trans'],kind='cubic')                         

    pwvs.append(0.5)
    trans.append(interp_cso_data(float(freq_in)))

    transfile='CSOtrans/CSO_075mm_'+str(wave)+'.txt'
    cso_data = np.loadtxt(transfile,delimiter='    ',usecols=(1,3),dtype={'names':('freq','trans'),'formats':('f4','f4')})
    interp_cso_data = scipy.interpolate.interp1d(cso_data['freq'],cso_data['trans'],kind='cubic')                         

    pwvs.append(0.75)
    trans.append(interp_cso_data(float(freq_in)))

    transfile='CSOtrans/CSO_1mm_'+str(wave)+'.txt'
    cso_data = np.loadtxt(transfile,delimiter='    ',usecols=(1,3),dtype={'names':('freq','trans'),'formats':('f4','f4')})
    interp_cso_data = scipy.interpolate.interp1d(cso_data['freq'],cso_data['trans'],kind='cubic')                         

    pwvs.append(1.0)
    trans.append(interp_cso_data(float(freq_in)))

    transfile='CSOtrans/CSO_15mm_'+str(wave)+'.txt'
    cso_data = np.loadtxt(transfile,delimiter='    ',usecols=(1,3),dtype={'names':('freq','trans'),'formats':('f4','f4')})
    interp_cso_data = scipy.interpolate.interp1d(cso_data['freq'],cso_data['trans'],kind='cubic')                         

    pwvs.append(1.5)
    trans.append(interp_cso_data(float(freq_in)))

    transfile='CSOtrans/CSO_2mm_'+str(wave)+'.txt'
    cso_data = np.loadtxt(transfile,delimiter='    ',usecols=(1,3),dtype={'names':('freq','trans'),'formats':('f4','f4')})
    interp_cso_data = scipy.interpolate.interp1d(cso_data['freq'],cso_data['trans'],kind='cubic')                         

    pwvs.append(2.0)
    trans.append(interp_cso_data(float(freq_in)))

    transfile='CSOtrans/CSO_25mm_'+str(wave)+'.txt'
    cso_data = np.loadtxt(transfile,delimiter='    ',usecols=(1,3),dtype={'names':('freq','trans'),'formats':('f4','f4')})
    interp_cso_data = scipy.interpolate.interp1d(cso_data['freq'],cso_data['trans'],kind='cubic')                         

    pwvs.append(2.5)
    trans.append(interp_cso_data(float(freq_in)))

    transfile='CSOtrans/CSO_3mm_'+str(wave)+'.txt'
    cso_data = np.loadtxt(transfile,delimiter='    ',usecols=(1,3),dtype={'names':('freq','trans'),'formats':('f4','f4')})
    interp_cso_data = scipy.interpolate.interp1d(cso_data['freq'],cso_data['trans'],kind='cubic')

    pwvs.append(3.0)
    trans.append(interp_cso_data(float(freq_in)))

    transfile='CSOtrans/CSO_35mm_'+str(wave)+'.txt'
    cso_data = np.loadtxt(transfile,delimiter='    ',usecols=(1,3),dtype={'names':('freq','trans'),'formats':('f4','f4')})
    interp_cso_data = scipy.interpolate.interp1d(cso_data['freq'],cso_data['trans'],kind='cubic')

    pwvs.append(3.5)
    trans.append(interp_cso_data(float(freq_in)))

    transfile='CSOtrans/CSO_4mm_'+str(wave)+'.txt'
    cso_data = np.loadtxt(transfile,delimiter='    ',usecols=(1,3),dtype={'names':('freq','trans'),'formats':('f4','f4')})
    interp_cso_data = scipy.interpolate.interp1d(cso_data['freq'],cso_data['trans'],kind='cubic')

    pwvs.append(4.0)
    trans.append(interp_cso_data(float(freq_in)))

    transfile='CSOtrans/CSO_45mm_'+str(wave)+'.txt'
    cso_data = np.loadtxt(transfile,delimiter='    ',usecols=(1,3),dtype={'names':('freq','trans'),'formats':('f4','f4')})
    interp_cso_data = scipy.interpolate.interp1d(cso_data['freq'],cso_data['trans'],kind='cubic')

    pwvs.append(4.5)
    trans.append(interp_cso_data(float(freq_in)))

    transfile='CSOtrans/CSO_5mm_'+str(wave)+'.txt'
    cso_data = np.loadtxt(transfile,delimiter='    ',usecols=(1,3),dtype={'names':('freq','trans'),'formats':('f4','f4')})
    interp_cso_data = scipy.interpolate.interp1d(cso_data['freq'],cso_data['trans'],kind='cubic')

    pwvs.append(5.0)
    trans.append(interp_cso_data(float(freq_in)))

    # Make a Transmission versus PWV function

    p=np.poly1d(np.polyfit(pwvs,trans,3))

    if MakeDict == 'yes':
        pwv_trans_dict = {}
        pwv_trans_dict['pwv'] = pwvs
        pwv_trans_dict['trans'] = trans

        pickle.dump(pwv_trans_dict,open('pwv_trans_'+str(int(wave))+'.bin','wb'))

    return (p)


def FitCSOTransvsPWVFast(wave): # Wave is: 850 or 450

    import numpy as np
    import pickle
  
    info_dict = pickle.load(open('CSOtrans/pwvfunction/pwv_trans_'+str(int(wave))+'.bin','rb'))

    p=np.poly1d(np.polyfit(info_dict['pwv'],info_dict['trans'],3))

    return(p)

    

####################
# This program works in conjunction
# with FITCSOTransvsPWV, above.
# Given a wavelength 
# and a TAU225zen measurement,
# it returns the expected Transmission.
####################

def CSOtrans(wave,tau225zen,TransvsPWVfunc): # Wave is: 850 or 450, tau225 is zenith value

    import numpy as np
    import scipy.interpolate
    import matplotlib.pyplot as plt

    c       = 299792459 # m/s
    freq_in = ( c / (wave*1e-6) ) / 1e9 # GHz

    PWVzen = (tau225zen - 0.017) / 0.04

    transmission_at_this_PWVzen = TransvsPWVfunc(PWVzen)

    return (transmission_at_this_PWVzen)

######
#
#
# This program measures the peak and total fluxes of the calibrator with GaussClumps
#
#
######

def fitcal(fname,paramfile):
    import subprocess
    import math
    import numpy as np
    import astropy.io.fits as apfits
    import os
    from mairs3 import sexiges

    ### get_noise ###

    #Convert the .sdf to .fits.
    fits_name = fname[:-4]+'.fits'
    sdf_conv_command = '$CONVERT_DIR/ndf2fits '+fname+' !'+fits_name
    subprocess.call(sdf_conv_command, shell=True)

    #Read in the variance extension of the file.
    vardata = apfits.getdata(fits_name, 1)
    
    #Find the noise in the central portion of the map.
    vardata_rad = vardata.shape[1]/2
    #Pick a region in the middle of the map to find the median noise. The X and
    # Y ranges will be the same.
    vardata_lo = int((vardata_rad-(vardata_rad/2)))
    vardata_hi = int((vardata_rad+(vardata_rad/2)))
    vardata_mid = vardata[0, vardata_lo:vardata_hi, vardata_lo:vardata_hi]
    noise = math.sqrt(np.nanmedian(vardata_mid))
    os.system('rm '+fits_name)

    # Now run Gaussclumps

    #Default outputs are "fname"+suffix.extension
    outputfile = fname[:-4]+'_clumps.sdf'
    outcatfile = fname[:-4]+'.FIT'
    logfile = fname[:-4]+'.txt'
    method = 'GaussClumps'
    noise_str = str( ('{:.5e}').format( noise ) )

    #Write out the .sh script that will be executed.
    shellscript = open('./findclumpscript.sh', 'w')
    shellscript.write('#!/bin/bash\n')
    shellscript.write('{ . $CUPID_DIR/cupid.sh; }\n')
    shellscript.write("findclumps '"
    +fname+"' '"+outputfile+"' '"+outcatfile+"' '"+method+"' CONFIG=^'"
    +paramfile+"' LOGFILE='"+logfile+"' REPCONF=TRUE RMS="+noise_str+
    " DECONV=FALSE MSG_FILTER=3 SHAPE=Ellipse WCSPAR=TRUE BACKOFF=TRUE") #Set MSG_FILTER to 3
    shellscript.close()
    subprocess.call('source ./findclumpscript.sh', shell='TRUE')
    os.system('rm findclumpscript.sh') 

    # Now extract the peak, total flux, and (FWHM1, FWHM2)
    # From the output table

    if os.path.isfile(outcatfile):
        gaussian_data = apfits.getdata(outcatfile)
        peak_flux     = gaussian_data['Peak'][0]
        total_flux    = gaussian_data['Sum'][0]
        FWHM1         = gaussian_data['GCFWHM1'][0]
        FWHM2         = gaussian_data['GCFWHM2'][0]
        area          = gaussian_data['Volume'][0]
        peakX1        = gaussian_data['Peak1'][0]
        peakY1        = gaussian_data['Peak2'][0]
        if peakX1 < 0.0001:
            peakX         = sexiges(0.0,'HMS')
        else:
            peakX         = sexiges(gaussian_data['Peak1'][0],'HMS')
        if peakY1 < 0.0001:
            peakY         = sexiges(0.0,'DMS')
        else:
            peakY         = sexiges(gaussian_data['Peak2'][0],'DMS')
    else:
        peak_flux     = np.nan
        total_flux    = np.nan
        FWHM1         = np.nan
        FWHM2         = np.nan
        area          = np.nan
        peakX         = np.nan
        peakY         = np.nan

    return peak_flux,total_flux,FWHM1,FWHM2,area,peakX,peakY

#####################
# After some testing with GaussClumps,
# Lets also use BeamFit and see if
# we get the same results/see if
# the tau relations are the same but the FCFs
# are different
#####################

def fitcalBF(fname,hdr,source):
    import subprocess
    import math
    import numpy as np
    import astropy.io.fits as apfits
    import os
    from mairs3 import sexiges
    from starlink import kappa

    RA1   = sexiges((hdr['OBSRABL']+hdr['OBSRABR'])/2.0,'HMS')
    RAhr  = RA1.split('h ')[0]
    RAmin = RA1.split('h ')[1].split('m ')[0]
    RAsec = RA1.split('m ')[1].split('s')[0]

    if int(RAhr)<10:
        RAhr='0'+RAhr

    RAstr   = RAhr+':'+RAmin+':'+RAsec
    DECstr  = sexiges((hdr['OBSDECTL']+hdr['OBSDECBL'])/2.0,'DMS') 
    if int(DECstr.split(':')[0])<10:
        DECstr = '0'+DECstr

    if source == 'URANUS':
        RAstr  = '00:00:00.00'
        DECstr = '00:00:00.00'

    try:
        BFresults = kappa.beamfit(fname,MODE='i',POS='"'+RAstr+','+DECstr+'"')#,LOGFILE=fname.split('.sdf')[0]+'_BFlog.txt')
        if BFresults is None:
            return(np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
        else:
            return(BFresults.amp[0],BFresults.amp[1],BFresults.sum,BFresults.majfwhm[0]*205264.0,BFresults.majfwhm[1]*206264.0,BFresults.minfwhm[0]*206264.0,BFresults.minfwhm[1]*206264.0)
    except:
        return(np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
    #beamfit_command = "$KAPPA_DIR/beamfit "+fname+' interface 1 '+"'"+'"'+RAstr+','+DECstr+'"'+"'"+' LOGFILE = '+fname.split('.sdf')[0]+'_BFlog.txt'
    #subprocess.call(beamfit_command, shell=True)
    

