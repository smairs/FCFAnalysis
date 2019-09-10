import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import datetime
import pickle

#catstring = '../catalogues/master_catalogue_850_4.6_0.0043_20190904.txt'
catstring = '../catalogues/master_catalogue_450_26.0_0.012_20190904.txt'

plot_date_min = Time('2019-01-01T00:00:00.00',format='isot').mjd
plot_date_max = Time('2019-09-30T00:00:00.00',format='isot').mjd 

# 58705 = 2019/08/10 - Back on sky 2 days later (protests)
# 58253 = 2018/05/15 - SMU Trouble
# 57693 = New Filters (2016/11/01)
# 58620.416667 Gain Turned Down to 70%
vert_lines    = [Time('2019-08-10T00:00:00.00',format='isot').mjd]

wave = catstring.split('_')[2]

######
######
######
# Apply ARP220 Flux Fix - WILL NOT WORK UNLESS INIDIVIDUAL PLOTS SELECTED

# Original Values:

ARP220_peak_850   = 0.79
ARP220_arcsec_850 = 0.81
ARP220_peak_450   = 5.2
ARP220_arcsec_450 = 5.4

# Enter updated value below and a fix will be applied to the plotted results
######

if wave == '850':
    ARP220_True_Peak  = 0.87
    ARP220_True_Total = 0.87
    ARP220_Peak_fix   = ARP220_True_Peak/ARP220_peak_850
    ARP220_Total_fix  = ARP220_True_Total/ARP220_arcsec_850

else:
    ARP220_True_Peak  = 5.2
    ARP220_True_Total = 5.4
    ARP220_Peak_fix   = ARP220_True_Peak/ARP220_peak_450
    ARP220_Total_fix  = ARP220_True_Total/ARP220_arcsec_450

#####
#####
#####



a = float(catstring.split('_')[3])
b = float(catstring.split('_')[4])

cat                     = np.genfromtxt(catstring,dtype=None,names=True)
epoch                   = 'All' #'Published','Silver WVM 1','No WVM 1','Silver WVM 2','No WVM 2','Black WVM','New Filters','Membrane Off','Membrane On','SMU Trouble','SMU Gain Fix','SMU HW Fix','Helco','All'
epochs_for_comparison   = ['Helco']
#epochs_for_comparison   = ['Published','Silver WVM 1','No WVM 1','Silver WVM 2','No WVM 2','Black WVM','New Filters','Membrane Off','Membrane On','SMU Trouble','SMU Gain Fix','SMU HW Fix','Helco']

plot_nominal_epoch_mean   = False
individual_source_plots   = True
saveplots                 = True

date_constraint_Obs_after_UT  = datetime.date(2011,5,1)
date_constraint_Obs_before_UT = datetime.date(2030,12,31)

# ALL TIMES OF DAY AND NIGHT
Time_Constraint_Obs_after_UT  = datetime.time(0,0,0)    # midnight UT  
Time_Constraint_Obs_before_UT = datetime.time(23,59,59) # 11:59 pm UT  

# EXTENDED OBSERVING
#Time_Constraint_Obs_after_UT  = datetime.time(17,30,0) # 7:30:00 am HST
#Time_Constraint_Obs_before_UT = datetime.time(22,59,59)# 12:59:59 pm HST

# SWEET SPOT
#Time_Constraint_Obs_after_UT  = datetime.time(7,0,0)  # 21:00 HST 
#Time_Constraint_Obs_before_UT = datetime.time(13,0,0) # 03:00 HST 

# EARLY EVENING
#Time_Constraint_Obs_after_UT  = datetime.time(3,0,0)    # 17:00 HST
#Time_Constraint_Obs_before_UT = datetime.time(7,0,0) # 21:00 HST

tag = '' # For naming the output files

tick_spacing            = 72000*2 # min
#tick_spacing            = 7200*24 # min
tick_spacing_mjd        = tick_spacing/1440 # 1440 minutes per day

############################################################
############################################################
############################################################

if wave=='450':
    ARCSEC_UPPER_LIM        = 8.#10.0
    ARCSEC_LOWER_LIM        = 1.
    PEAK_UPPER_LIM          = 1100.0 #1200.0
    PEAK_LOWER_LIM          = 200.0
if wave=='850':
    ARCSEC_UPPER_LIM        = 3.0 #8.
    ARCSEC_LOWER_LIM        = 1.7 #1.
    PEAK_UPPER_LIM          = 800.0#1000.0
    PEAK_LOWER_LIM          = 200.0

if tag != '':
    ending                  = '_'+wave+'_'+epoch.replace(' ','_')+'_oldextcor'+'_'+tag+'.png'
else:
    ending                  = '_'+wave+'_'+epoch.replace(' ','_')+'_oldextcor.png'


# RxA Warmups - which can apparently affect SCUBA-2 Performance!

RxA_warmups   = [Time("2015-11-21T05:00:00.00",format='isot',scale='utc').mjd,Time("2016-08-08T00:00:00.00",format='isot',scale='utc').mjd,Time("2016-10-21T03:00:00.00",format='isot',scale='utc').mjd,Time("2016-11-02T02:00:00.00",format='isot',scale='utc').mjd,Time("2016-11-09T06:00:00.00",format='isot',scale='utc').mjd,Time("2017-01-19T10:00:00.00",format='isot',scale='utc').mjd,Time("2017-03-16T21:00:00.00",format='isot',scale='utc').mjd,Time("2017-05-23T16:00:00.00",format='isot',scale='utc').mjd,Time("2017-07-03T10:00:00.00",format='isot',scale='utc').mjd,Time("2018-02-21T22:00:00.00",format='isot',scale='utc').mjd,Time("2018-02-28T10:00:00.00",format='isot',scale='utc').mjd,Time("2018-06-26T10:00:00.00",format='isot',scale='utc').mjd]
RxA_cooldowns = [Time("2016-01-06T02:00:00.00",format='isot',scale='utc').mjd,Time("2016-08-12T02:00:00.00",format='isot',scale='utc').mjd,Time("2016-10-28T05:00:00.00",format='isot',scale='utc').mjd,Time("2016-11-03T02:00:00.00",format='isot',scale='utc').mjd,Time("2016-11-10T06:00:00.00",format='isot',scale='utc').mjd,Time("2017-01-30T10:00:00.00",format='isot',scale='utc').mjd,Time("2017-03-24T00:00:00.00",format='isot',scale='utc').mjd,Time("2017-06-07T10:00:00.00",format='isot',scale='utc').mjd,Time("2017-07-10T10:00:00.00",format='isot',scale='utc').mjd,Time("2018-02-22T22:00:00.00",format='isot',scale='utc').mjd,Time("2018-03-06T10:00:00.00",format='isot',scale='utc').mjd]

# EPOCHS - ALL INCLUSIVE!!!

PUBMJDstart        = Time("2011-05-01T00:00:00.00",format='isot',scale='utc').mjd
PUBMJDend          = Time("2012-05-31T00:00:00.00",format='isot',scale='utc').mjd
SILMJDstart        = Time("2012-06-01T00:00:00.00",format='isot',scale='utc').mjd
WVM_out_of_service = Time("2013-03-15T00:00:00.00",format='isot',scale='utc').mjd
WVM_back_in_service= Time("2013-04-09T00:00:00.00",format='isot',scale='utc').mjd
SILMJDend          = Time("2015-01-27T00:00:00.00",format='isot',scale='utc').mjd
BLAMJDstart        = Time("2015-04-10T00:00:00.00",format='isot',scale='utc').mjd
BLAMJDend          = Time("2016-10-05T00:00:00.00",format='isot',scale='utc').mjd
NEWFMJDstart       = Time("2016-10-06T00:00:00.00",format='isot',scale='utc').mjd
NEWFMJDend         = Time("2017-12-05T00:00:00.00",format='isot',scale='utc').mjd
MEMOFFMJDstart     = Time("2017-12-06T00:00:00.00",format='isot',scale='utc').mjd
MEMOFFMJDend       = Time("2018-01-10T00:00:00.00",format='isot',scale='utc').mjd
MEMONMJDstart      = Time("2018-01-11T00:00:00.00",format='isot',scale='utc').mjd
MEMONMJDend        = Time("2018-05-01T00:00:00.00",format='isot',scale='utc').mjd
SMUTMJDstart       = Time("2018-05-02T00:00:00.00",format='isot',scale='utc').mjd
SMUTMJDend         = Time("2018-06-30T08:10:00.00",format='isot',scale='utc').mjd
SMUGMJDstart       = Time("2018-06-30T08:11:00.00",format='isot',scale='utc').mjd
SMUGMJDend         = Time("2018-07-25T00:00:00.00",format='isot',scale='utc').mjd
SMUHWMJDstart      = Time("2018-07-28T00:00:00.00",format='isot',scale='utc').mjd
SMUHWMJDend        = Time("2018-11-26T00:00:00.00",format='isot',scale='utc').mjd
CoolDownHelcostart = Time("2018-11-26T00:00:01.00",format='isot',scale='utc').mjd
CoolDownHelcoend   = Time("2500-12-30T00:00:00.00",format='isot',scale='utc').mjd

epochs             = np.array(['Published','Silver WVM 1','No WVM 1','Silver WVM 2','No WVM 2','Black WVM','New Filters','Membrane Off','Membrane On','SMU Trouble','SMU Gain Fix','SMU HW Fix','Helco','All'])
#yearcodes          = np.array(['20110501','20120601,20140101','20120601','20120601,20140101','20150127','20150410','20161101','20171206','20180111','20180501','20180630','20180724','20181126','20110501,20120601,20140101,20150127,20150410,20161101,20171206,20180111,20180501,20180630,20180724,20181126'])

epoch_dates        = np.array([(PUBMJDstart,PUBMJDend),(SILMJDstart,WVM_out_of_service),(WVM_out_of_service+1,WVM_back_in_service-1),(WVM_back_in_service,SILMJDend),(SILMJDend+1,BLAMJDstart-1),(BLAMJDstart,BLAMJDend),(NEWFMJDstart,NEWFMJDend),(MEMOFFMJDstart,MEMOFFMJDend),(MEMONMJDstart,MEMONMJDend),(SMUTMJDstart,SMUTMJDend),(SMUGMJDstart,SMUGMJDend),(SMUHWMJDstart,SMUHWMJDend),(CoolDownHelcostart,CoolDownHelcoend),(PUBMJDstart,CoolDownHelcoend)])
# Approximate: Orange, Light Brown, Bright Green, Bright Pink, Gold, Purple, Blue-Green, Red, Dark Green, Soft Pink, Dark Purple, Gray, Blue, Black
comparison_colours = np.array(['#ff4f00','#a48662','#73e309','#ff00c1','#ffb400','#a33260','#0cecd8','r','#2d8938','#ccaabb','#502a72','#95a5a6','b','k'])


MJDstart       = epoch_dates[np.where(epochs == epoch)][0][0]
MJDend         = epoch_dates[np.where(epochs == epoch)][0][1] 


comparison_dict = {}

for eachcomparison in epochs_for_comparison:

    COMPAREMJDstart = epoch_dates[np.where(epochs == eachcomparison)][0][0]
    COMPAREMJDend   = epoch_dates[np.where(epochs == eachcomparison)][0][1]

    
    comparison_dict[eachcomparison] = {}
    if individual_source_plots:
        comparison_dict[eachcomparison]['CRL2688']    = {}
        comparison_dict[eachcomparison]['CRL618']     = {}
        comparison_dict[eachcomparison]['MARS']       = {}
        comparison_dict[eachcomparison]['NEPTUNE']    = {}
        comparison_dict[eachcomparison]['URANUS']     = {}
        comparison_dict[eachcomparison]['Arp220']     = {}
    else:
        comparison_dict[eachcomparison]['ALLSOURCES'] = {}

    comparison_AR         = []
    comparison_SOURCE     = []
    comparison_Peak_FCF   = []
    comparison_Arcsec_FCF = []
    comparison_weather    = []
    comparison_AM         = []
    comparison_TauTimesAM = []
    comparison_PEAK_FLUX  = []
    comparison_TOTAL_FLUX = []
    comparison_FWHM1      = []
    comparison_FWHM2      = [] 
    comparison_Beam_Eff   = []
    comparison_MJD        = []
    
    for i in range(len(cat['MJDST'])):
        if cat['MJDST'][i]>=COMPAREMJDstart and cat['MJDST'][i]<=COMPAREMJDend: 
#            if cat['MJDST'][i] < SILMJDend or cat['MJDST'][i] > BLAMJDstart: # No reliable WVM here
#                if cat['MJDST'][i]< WVM_out_of_service or cat['MJDST'][i] > WVM_back_in_service:
            if datetime.date(Time(cat['MJDST'][i],format='mjd').datetime.year,Time(cat['MJDST'][i],format='mjd').datetime.month,Time(cat['MJDST'][i],format='mjd').datetime.day) >= date_constraint_Obs_after_UT and datetime.date(Time(cat['MJDST'][i],format='mjd').datetime.year,Time(cat['MJDST'][i],format='mjd').datetime.month,Time(cat['MJDST'][i],format='mjd').datetime.day) <= date_constraint_Obs_before_UT:
                    if datetime.time(Time(cat['MJDST'][i],format='mjd').datetime.hour,Time(cat['MJDST'][i],format='mjd').datetime.minute,Time(cat['MJDST'][i],format='mjd').datetime.second) >= Time_Constraint_Obs_after_UT and datetime.time(Time(cat['MJDST'][i],format='mjd').datetime.hour,Time(cat['MJDST'][i],format='mjd').datetime.minute,Time(cat['MJDST'][i],format='mjd').datetime.second) <= Time_Constraint_Obs_before_UT:
                        if cat['ARCSEC_FCF'][i] < ARCSEC_UPPER_LIM and cat['ARCSEC_FCF'][i] > ARCSEC_LOWER_LIM:
                            if cat['PEAK_FCF'][i] < PEAK_UPPER_LIM and cat['PEAK_FCF'][i] > PEAK_LOWER_LIM:
                                comparison_AR.append(max(cat['FWHM1'][i],cat['FWHM2'][i])/min(cat['FWHM1'][i],cat['FWHM2'][i]))
                                comparison_Peak_FCF.append(cat['PEAK_FCF'][i])
                                comparison_SOURCE.append(cat['SOURCE_ID'][i]) 
                                comparison_Arcsec_FCF.append(cat['ARCSEC_FCF'][i])
                                comparison_weather.append((cat['WVMTAUST'][i]+cat['WVMTAUEN'][i])/2.0)
                                comparison_AM.append((cat['AMSTART'][i]+cat['AMEND'][i])/2.0)
                                comparison_TauTimesAM.append(((cat['WVMTAUST'][i]+cat['WVMTAUEN'][i])/2.0)*((cat['AMSTART'][i]+cat['AMEND'][i])/2.0))
                                comparison_PEAK_FLUX.append(cat['PEAK_FLUX'][i])
                                comparison_TOTAL_FLUX.append(cat['TOTAL_FLUX'][i])
                                comparison_FWHM1.append(cat['FWHM1'][i])
                                comparison_FWHM2.append(cat['FWHM2'][i])
                                comparison_Beam_Eff.append(np.sqrt((cat['PEAK_FCF'][i]/cat['ARCSEC_FCF'][i])/1.133))
                                comparison_MJD.append(cat['MJDST'][i])
    
    
    comparison_AR         = np.array(comparison_AR)
    comparison_SOURCE     = np.array(comparison_SOURCE)
    # Switch up Arp220 source IDs - I know this is stupid, but it is what it is.
    comparison_SOURCE[np.where(comparison_SOURCE==0)] = 6
    comparison_Peak_FCF   = np.array(comparison_Peak_FCF)
    comparison_Arcsec_FCF = np.array(comparison_Arcsec_FCF)
    comparison_weather    = np.array(comparison_weather)
    comparison_AM         = np.array(comparison_AM) 
    comparison_TauTimesAM = np.array(comparison_TauTimesAM)
    comparison_PEAK_FLUX  = np.array(comparison_PEAK_FLUX)
    comparison_TOTAL_FLUX = np.array(comparison_TOTAL_FLUX)
    comparison_FWHM1      = np.array(comparison_FWHM1)
    comparison_FWHM2      = np.array(comparison_FWHM2)
    comparison_Beam_Eff   = np.array(comparison_Beam_Eff)
    comparison_MJD        = np.array(comparison_MJD)
    
    source_1_ind    = np.where(comparison_SOURCE==1)[0]
    source_2_ind    = np.where(comparison_SOURCE==2)[0]
    source_3_ind    = np.where(comparison_SOURCE==3)[0]
    source_4_ind    = np.where(comparison_SOURCE==4)[0]
    source_5_ind    = np.where(comparison_SOURCE==5)[0]
    source_6_ind    = np.where(comparison_SOURCE==6)[0]

    if individual_source_plots:
        if len(source_1_ind)>3:
            source_1_mean_AR         = np.mean(comparison_AR[source_1_ind])
            source_1_AR_SD           = np.std(comparison_AR[source_1_ind],ddof=1)
            source_1_mean_Peak_FCF   = np.mean(comparison_Peak_FCF[source_1_ind])
            source_1_Peak_FCF_SD     = np.std(comparison_Peak_FCF[source_1_ind],ddof=1)
            source_1_mean_Arcsec_FCF = np.mean(comparison_Arcsec_FCF[source_1_ind])
            source_1_Arcsec_FCF_SD   = np.std(comparison_Arcsec_FCF[source_1_ind],ddof=1)
            source_1_mean_FWHM1         = np.mean(comparison_FWHM1[source_1_ind])
            source_1_FWHM1_SD           = np.std(comparison_FWHM1[source_1_ind],ddof=1)
            source_1_mean_FWHM2         = np.mean(comparison_FWHM2[source_1_ind])
            source_1_FWHM2_SD           = np.std(comparison_FWHM2[source_1_ind],ddof=1)
            source_1_mean_Beam_Eff      = np.mean(comparison_Beam_Eff[source_1_ind])
            source_1_Beam_Eff_SD        = np.std(comparison_Beam_Eff[source_1_ind],ddof=1)
            comparison_dict[eachcomparison]['CRL2688']['ARs']                = comparison_AR[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['Beam_Eff']           = comparison_Beam_Eff[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['ARCSEC_FCF']         = comparison_Arcsec_FCF[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['PEAK_FCF']           = comparison_Peak_FCF[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['TAU225']             = comparison_weather[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['AM']                 = comparison_AM[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['PEAK_FLUXES']        = comparison_PEAK_FLUX[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['TOTAL_FLUXES']       = comparison_TOTAL_FLUX[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['FWHM1']              = comparison_FWHM1[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['FWHM2']              = comparison_FWHM2[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['MJD']                = comparison_MJD[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['mean_Beam_eff']      = source_1_mean_Beam_Eff
            comparison_dict[eachcomparison]['CRL2688']['SD_Beam_eff']        = source_1_Beam_Eff_SD
            comparison_dict[eachcomparison]['CRL2688']['mean_FWHM1']         = source_1_mean_FWHM1
            comparison_dict[eachcomparison]['CRL2688']['SD_FWHM1']           = source_1_FWHM1_SD
            comparison_dict[eachcomparison]['CRL2688']['mean_FWHM2']         = source_1_mean_FWHM2
            comparison_dict[eachcomparison]['CRL2688']['SD_FWHM2']           = source_1_FWHM2_SD
            comparison_dict[eachcomparison]['CRL2688']['mean_AR']            = source_1_mean_AR
            comparison_dict[eachcomparison]['CRL2688']['SD_AR']              = source_1_AR_SD
            comparison_dict[eachcomparison]['CRL2688']['mean_Peak_FCF']      = source_1_mean_Peak_FCF
            comparison_dict[eachcomparison]['CRL2688']['SD_Peak_FCF']        = source_1_Peak_FCF_SD
            comparison_dict[eachcomparison]['CRL2688']['mean_Arcsec_FCF']    = source_1_mean_Arcsec_FCF
            comparison_dict[eachcomparison]['CRL2688']['SD_Arcsec_FCF']      = source_1_Arcsec_FCF_SD
        else:
            comparison_dict[eachcomparison]['CRL2688']['ARs']                = comparison_AR[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['Beam_Eff']           = comparison_Beam_Eff[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['ARCSEC_FCF']         = comparison_Arcsec_FCF[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['PEAK_FCF']           = comparison_Peak_FCF[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['TAU225']             = comparison_weather[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['AM']                 = comparison_AM[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['PEAK_FLUXES']        = comparison_PEAK_FLUX[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['TOTAL_FLUXES']       = comparison_TOTAL_FLUX[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['FWHM1']              = comparison_FWHM1[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['FWHM2']              = comparison_FWHM2[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['MJD']                = comparison_MJD[source_1_ind]
            comparison_dict[eachcomparison]['CRL2688']['mean_Beam_eff']      = np.nan
            comparison_dict[eachcomparison]['CRL2688']['SD_Beam_eff']        = np.nan
            comparison_dict[eachcomparison]['CRL2688']['mean_FWHM1']         = np.nan
            comparison_dict[eachcomparison]['CRL2688']['SD_FWHM1']           = np.nan
            comparison_dict[eachcomparison]['CRL2688']['mean_FWHM2']         = np.nan
            comparison_dict[eachcomparison]['CRL2688']['SD_FWHM2']           = np.nan
            comparison_dict[eachcomparison]['CRL2688']['mean_AR']            = np.nan 
            comparison_dict[eachcomparison]['CRL2688']['SD_AR']              = np.nan 
            comparison_dict[eachcomparison]['CRL2688']['mean_Peak_FCF']      = np.nan 
            comparison_dict[eachcomparison]['CRL2688']['SD_Peak_FCF']        = np.nan 
            comparison_dict[eachcomparison]['CRL2688']['mean_Arcsec_FCF']    = np.nan 
            comparison_dict[eachcomparison]['CRL2688']['SD_Arcsec_FCF']      = np.nan 
        
        if len(source_2_ind)>3:
            source_2_mean_AR         = np.mean(comparison_AR[source_2_ind])
            source_2_AR_SD           = np.std(comparison_AR[source_2_ind],ddof=1)
            source_2_mean_Peak_FCF   = np.mean(comparison_Peak_FCF[source_2_ind])
            source_2_Peak_FCF_SD     = np.std(comparison_Peak_FCF[source_2_ind],ddof=1)
            source_2_mean_Arcsec_FCF = np.mean(comparison_Arcsec_FCF[source_2_ind])
            source_2_Arcsec_FCF_SD   = np.std(comparison_Arcsec_FCF[source_2_ind],ddof=1)
            source_2_mean_FWHM1         = np.mean(comparison_FWHM1[source_2_ind])
            source_2_FWHM1_SD           = np.std(comparison_FWHM1[source_2_ind],ddof=1)
            source_2_mean_FWHM2         = np.mean(comparison_FWHM2[source_2_ind])
            source_2_FWHM2_SD           = np.std(comparison_FWHM2[source_2_ind],ddof=1)
            source_2_mean_Beam_Eff      = np.mean(comparison_Beam_Eff[source_2_ind])
            source_2_Beam_Eff_SD        = np.std(comparison_Beam_Eff[source_2_ind],ddof=1)
            comparison_dict[eachcomparison]['CRL618']['ARs']                = comparison_AR[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['Beam_Eff']           = comparison_Beam_Eff[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['ARCSEC_FCF']         = comparison_Arcsec_FCF[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['PEAK_FCF']           = comparison_Peak_FCF[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['TAU225']             = comparison_weather[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['AM']                 = comparison_AM[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['PEAK_FLUXES']        = comparison_PEAK_FLUX[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['TOTAL_FLUXES']       = comparison_TOTAL_FLUX[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['FWHM1']              = comparison_FWHM1[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['FWHM2']              = comparison_FWHM2[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['MJD']                = comparison_MJD[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['mean_Beam_eff']      = source_2_mean_Beam_Eff
            comparison_dict[eachcomparison]['CRL618']['SD_Beam_eff']        = source_2_Beam_Eff_SD
            comparison_dict[eachcomparison]['CRL618']['mean_FWHM1']         = source_2_mean_FWHM1
            comparison_dict[eachcomparison]['CRL618']['SD_FWHM1']           = source_2_FWHM1_SD
            comparison_dict[eachcomparison]['CRL618']['mean_FWHM2']         = source_2_mean_FWHM2
            comparison_dict[eachcomparison]['CRL618']['SD_FWHM2']           = source_2_FWHM2_SD
            comparison_dict[eachcomparison]['CRL618']['mean_AR']         = source_2_mean_AR
            comparison_dict[eachcomparison]['CRL618']['SD_AR']           = source_2_AR_SD
            comparison_dict[eachcomparison]['CRL618']['mean_Peak_FCF']   = source_2_mean_Peak_FCF
            comparison_dict[eachcomparison]['CRL618']['SD_Peak_FCF']     = source_2_Peak_FCF_SD
            comparison_dict[eachcomparison]['CRL618']['mean_Arcsec_FCF'] = source_2_mean_Arcsec_FCF
            comparison_dict[eachcomparison]['CRL618']['SD_Arcsec_FCF']   = source_2_Arcsec_FCF_SD
        else:
            comparison_dict[eachcomparison]['CRL618']['ARs']                = comparison_AR[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['Beam_Eff']           = comparison_Beam_Eff[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['ARCSEC_FCF']         = comparison_Arcsec_FCF[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['PEAK_FCF']           = comparison_Peak_FCF[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['TAU225']             = comparison_weather[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['AM']                 = comparison_AM[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['PEAK_FLUXES']        = comparison_PEAK_FLUX[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['TOTAL_FLUXES']       = comparison_TOTAL_FLUX[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['FWHM1']              = comparison_FWHM1[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['FWHM2']              = comparison_FWHM2[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['MJD']                = comparison_MJD[source_2_ind]
            comparison_dict[eachcomparison]['CRL618']['mean_Beam_eff']      = np.nan
            comparison_dict[eachcomparison]['CRL618']['SD_Beam_eff']        = np.nan
            comparison_dict[eachcomparison]['CRL618']['mean_FWHM1']         = np.nan
            comparison_dict[eachcomparison]['CRL618']['SD_FWHM1']           = np.nan
            comparison_dict[eachcomparison]['CRL618']['mean_FWHM2']         = np.nan
            comparison_dict[eachcomparison]['CRL618']['SD_FWHM2']           = np.nan
            comparison_dict[eachcomparison]['CRL618']['mean_AR']         = np.nan
            comparison_dict[eachcomparison]['CRL618']['SD_AR']           = np.nan
            comparison_dict[eachcomparison]['CRL618']['mean_Peak_FCF']   = np.nan
            comparison_dict[eachcomparison]['CRL618']['SD_Peak_FCF']     = np.nan
            comparison_dict[eachcomparison]['CRL618']['mean_Arcsec_FCF'] = np.nan
            comparison_dict[eachcomparison]['CRL618']['SD_Arcsec_FCF']   = np.nan
        
        if len(source_3_ind)>3:
            source_3_mean_AR         = np.mean(comparison_AR[source_3_ind])
            source_3_AR_SD           = np.std(comparison_AR[source_3_ind],ddof=1)
            source_3_mean_Peak_FCF   = np.mean(comparison_Peak_FCF[source_3_ind])
            source_3_Peak_FCF_SD     = np.std(comparison_Peak_FCF[source_3_ind],ddof=1)
            source_3_mean_Arcsec_FCF = np.mean(comparison_Arcsec_FCF[source_3_ind])
            source_3_Arcsec_FCF_SD   = np.std(comparison_Arcsec_FCF[source_3_ind],ddof=1)
            source_3_mean_FWHM1         = np.mean(comparison_FWHM1[source_3_ind])
            source_3_FWHM1_SD           = np.std(comparison_FWHM1[source_3_ind],ddof=1)
            source_3_mean_FWHM2         = np.mean(comparison_FWHM2[source_3_ind])
            source_3_FWHM2_SD           = np.std(comparison_FWHM2[source_3_ind],ddof=1)
            source_3_mean_Beam_Eff      = np.mean(comparison_Beam_Eff[source_3_ind])
            source_3_Beam_Eff_SD        = np.std(comparison_Beam_Eff[source_3_ind],ddof=1)
            comparison_dict[eachcomparison]['MARS']['ARs']                = comparison_AR[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['Beam_Eff']           = comparison_Beam_Eff[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['ARCSEC_FCF']         = comparison_Arcsec_FCF[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['PEAK_FCF']           = comparison_Peak_FCF[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['TAU225']             = comparison_weather[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['AM']                 = comparison_AM[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['PEAK_FLUXES']        = comparison_PEAK_FLUX[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['TOTAL_FLUXES']       = comparison_TOTAL_FLUX[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['FWHM1']              = comparison_FWHM1[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['FWHM2']              = comparison_FWHM2[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['MJD']                = comparison_MJD[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['mean_Beam_eff']      = source_3_mean_Beam_Eff
            comparison_dict[eachcomparison]['MARS']['SD_Beam_eff']        = source_3_Beam_Eff_SD
            comparison_dict[eachcomparison]['MARS']['mean_FWHM1']         = source_3_mean_FWHM1
            comparison_dict[eachcomparison]['MARS']['SD_FWHM1']           = source_3_FWHM1_SD
            comparison_dict[eachcomparison]['MARS']['mean_FWHM2']         = source_3_mean_FWHM2
            comparison_dict[eachcomparison]['MARS']['SD_FWHM2']           = source_3_FWHM2_SD
            comparison_dict[eachcomparison]['MARS']['mean_AR']         = source_3_mean_AR
            comparison_dict[eachcomparison]['MARS']['SD_AR']           = source_3_AR_SD
            comparison_dict[eachcomparison]['MARS']['mean_Peak_FCF']   = source_3_mean_Peak_FCF
            comparison_dict[eachcomparison]['MARS']['SD_Peak_FCF']     = source_3_Peak_FCF_SD
            comparison_dict[eachcomparison]['MARS']['mean_Arcsec_FCF'] = source_3_mean_Arcsec_FCF
            comparison_dict[eachcomparison]['MARS']['SD_Arcsec_FCF']   = source_3_Arcsec_FCF_SD    
        else:
            comparison_dict[eachcomparison]['MARS']['ARs']                = comparison_AR[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['Beam_Eff']           = comparison_Beam_Eff[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['ARCSEC_FCF']         = comparison_Arcsec_FCF[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['PEAK_FCF']           = comparison_Peak_FCF[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['TAU225']             = comparison_weather[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['AM']                 = comparison_AM[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['PEAK_FLUXES']        = comparison_PEAK_FLUX[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['TOTAL_FLUXES']       = comparison_TOTAL_FLUX[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['FWHM1']              = comparison_FWHM1[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['FWHM2']              = comparison_FWHM2[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['MJD']                = comparison_MJD[source_3_ind]
            comparison_dict[eachcomparison]['MARS']['mean_Beam_eff']      = np.nan
            comparison_dict[eachcomparison]['MARS']['SD_Beam_eff']        = np.nan
            comparison_dict[eachcomparison]['MARS']['mean_FWHM1']         = np.nan
            comparison_dict[eachcomparison]['MARS']['SD_FWHM1']           = np.nan
            comparison_dict[eachcomparison]['MARS']['mean_FWHM2']         = np.nan
            comparison_dict[eachcomparison]['MARS']['SD_FWHM2']           = np.nan
            comparison_dict[eachcomparison]['MARS']['mean_AR']         = np.nan
            comparison_dict[eachcomparison]['MARS']['SD_AR']           = np.nan
            comparison_dict[eachcomparison]['MARS']['mean_Peak_FCF']   = np.nan
            comparison_dict[eachcomparison]['MARS']['SD_Peak_FCF']     = np.nan
            comparison_dict[eachcomparison]['MARS']['mean_Arcsec_FCF'] = np.nan
            comparison_dict[eachcomparison]['MARS']['SD_Arcsec_FCF']   = np.nan
    
        if len(source_4_ind)>3:
            source_4_mean_AR         = np.mean(comparison_AR[source_4_ind])
            source_4_AR_SD           = np.std(comparison_AR[source_4_ind],ddof=1)
            source_4_mean_Peak_FCF   = np.mean(comparison_Peak_FCF[source_4_ind])
            source_4_Peak_FCF_SD     = np.std(comparison_Peak_FCF[source_4_ind],ddof=1)
            source_4_mean_Arcsec_FCF = np.mean(comparison_Arcsec_FCF[source_4_ind])
            source_4_Arcsec_FCF_SD   = np.std(comparison_Arcsec_FCF[source_4_ind],ddof=1)
            source_4_mean_FWHM1         = np.mean(comparison_FWHM1[source_4_ind])
            source_4_FWHM1_SD           = np.std(comparison_FWHM1[source_4_ind],ddof=1)
            source_4_mean_FWHM2         = np.mean(comparison_FWHM2[source_4_ind])
            source_4_FWHM2_SD           = np.std(comparison_FWHM2[source_4_ind],ddof=1)
            source_4_mean_Beam_Eff      = np.mean(comparison_Beam_Eff[source_4_ind])
            source_4_Beam_Eff_SD        = np.std(comparison_Beam_Eff[source_4_ind],ddof=1)
            comparison_dict[eachcomparison]['NEPTUNE']['ARs']                = comparison_AR[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['Beam_Eff']           = comparison_Beam_Eff[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['ARCSEC_FCF']         = comparison_Arcsec_FCF[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['PEAK_FCF']           = comparison_Peak_FCF[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['TAU225']             = comparison_weather[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['AM']                 = comparison_AM[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['PEAK_FLUXES']        = comparison_PEAK_FLUX[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['TOTAL_FLUXES']       = comparison_TOTAL_FLUX[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['FWHM1']              = comparison_FWHM1[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['FWHM2']              = comparison_FWHM2[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['MJD']                = comparison_MJD[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['mean_Beam_eff']      = source_4_mean_Beam_Eff
            comparison_dict[eachcomparison]['NEPTUNE']['SD_Beam_eff']        = source_4_Beam_Eff_SD
            comparison_dict[eachcomparison]['NEPTUNE']['mean_FWHM1']         = source_4_mean_FWHM1
            comparison_dict[eachcomparison]['NEPTUNE']['SD_FWHM1']           = source_4_FWHM1_SD
            comparison_dict[eachcomparison]['NEPTUNE']['mean_FWHM2']         = source_4_mean_FWHM2
            comparison_dict[eachcomparison]['NEPTUNE']['SD_FWHM2']           = source_4_FWHM2_SD
            comparison_dict[eachcomparison]['NEPTUNE']['mean_AR']         = source_4_mean_AR
            comparison_dict[eachcomparison]['NEPTUNE']['SD_AR']           = source_4_AR_SD
            comparison_dict[eachcomparison]['NEPTUNE']['mean_Peak_FCF']   = source_4_mean_Peak_FCF
            comparison_dict[eachcomparison]['NEPTUNE']['SD_Peak_FCF']     = source_4_Peak_FCF_SD
            comparison_dict[eachcomparison]['NEPTUNE']['mean_Arcsec_FCF'] = source_4_mean_Arcsec_FCF
            comparison_dict[eachcomparison]['NEPTUNE']['SD_Arcsec_FCF']   = source_4_Arcsec_FCF_SD
        else:
            comparison_dict[eachcomparison]['NEPTUNE']['ARs']                = comparison_AR[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['Beam_Eff']           = comparison_Beam_Eff[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['ARCSEC_FCF']         = comparison_Arcsec_FCF[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['PEAK_FCF']           = comparison_Peak_FCF[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['TAU225']             = comparison_weather[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['AM']                 = comparison_AM[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['PEAK_FLUXES']        = comparison_PEAK_FLUX[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['TOTAL_FLUXES']       = comparison_TOTAL_FLUX[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['FWHM1']              = comparison_FWHM1[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['FWHM2']              = comparison_FWHM2[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['MJD']                = comparison_MJD[source_4_ind]
            comparison_dict[eachcomparison]['NEPTUNE']['mean_Beam_eff']      = np.nan
            comparison_dict[eachcomparison]['NEPTUNE']['SD_Beam_eff']        = np.nan
            comparison_dict[eachcomparison]['NEPTUNE']['mean_FWHM1']         = np.nan
            comparison_dict[eachcomparison]['NEPTUNE']['SD_FWHM1']           = np.nan
            comparison_dict[eachcomparison]['NEPTUNE']['mean_FWHM2']         = np.nan
            comparison_dict[eachcomparison]['NEPTUNE']['SD_FWHM2']           = np.nan
            comparison_dict[eachcomparison]['NEPTUNE']['mean_AR']         = np.nan
            comparison_dict[eachcomparison]['NEPTUNE']['SD_AR']           = np.nan
            comparison_dict[eachcomparison]['NEPTUNE']['mean_Peak_FCF']   = np.nan
            comparison_dict[eachcomparison]['NEPTUNE']['SD_Peak_FCF']     = np.nan
            comparison_dict[eachcomparison]['NEPTUNE']['mean_Arcsec_FCF'] = np.nan
            comparison_dict[eachcomparison]['NEPTUNE']['SD_Arcsec_FCF']   = np.nan
        
        if len(source_5_ind)>3:
            source_5_mean_AR         = np.mean(comparison_AR[source_5_ind])
            source_5_AR_SD           = np.std(comparison_AR[source_5_ind],ddof=1)
            source_5_mean_Peak_FCF   = np.mean(comparison_Peak_FCF[source_5_ind])
            source_5_Peak_FCF_SD     = np.std(comparison_Peak_FCF[source_5_ind],ddof=1)
            source_5_mean_Arcsec_FCF = np.mean(comparison_Arcsec_FCF[source_5_ind])
            source_5_Arcsec_FCF_SD   = np.std(comparison_Arcsec_FCF[source_5_ind],ddof=1)
            source_5_mean_FWHM1         = np.mean(comparison_FWHM1[source_5_ind])
            source_5_FWHM1_SD           = np.std(comparison_FWHM1[source_5_ind],ddof=1)
            source_5_mean_FWHM2         = np.mean(comparison_FWHM2[source_5_ind])
            source_5_FWHM2_SD           = np.std(comparison_FWHM2[source_5_ind],ddof=1)
            source_5_mean_Beam_Eff      = np.mean(comparison_Beam_Eff[source_5_ind])
            source_5_Beam_Eff_SD        = np.std(comparison_Beam_Eff[source_5_ind],ddof=1)
            comparison_dict[eachcomparison]['URANUS']['ARs']                = comparison_AR[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['Beam_Eff']           = comparison_Beam_Eff[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['ARCSEC_FCF']         = comparison_Arcsec_FCF[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['PEAK_FCF']           = comparison_Peak_FCF[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['TAU225']             = comparison_weather[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['AM']                 = comparison_AM[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['PEAK_FLUXES']        = comparison_PEAK_FLUX[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['TOTAL_FLUXES']       = comparison_TOTAL_FLUX[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['FWHM1']              = comparison_FWHM1[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['FWHM2']              = comparison_FWHM2[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['MJD']                = comparison_MJD[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['mean_Beam_eff']      = source_5_mean_Beam_Eff
            comparison_dict[eachcomparison]['URANUS']['SD_Beam_eff']        = source_5_Beam_Eff_SD
            comparison_dict[eachcomparison]['URANUS']['mean_FWHM1']         = source_5_mean_FWHM1
            comparison_dict[eachcomparison]['URANUS']['SD_FWHM1']           = source_5_FWHM1_SD
            comparison_dict[eachcomparison]['URANUS']['mean_FWHM2']         = source_5_mean_FWHM2
            comparison_dict[eachcomparison]['URANUS']['SD_FWHM2']           = source_5_FWHM2_SD
            comparison_dict[eachcomparison]['URANUS']['mean_AR']         = source_5_mean_AR
            comparison_dict[eachcomparison]['URANUS']['SD_AR']           = source_5_AR_SD
            comparison_dict[eachcomparison]['URANUS']['mean_Peak_FCF']   = source_5_mean_Peak_FCF
            comparison_dict[eachcomparison]['URANUS']['SD_Peak_FCF']     = source_5_Peak_FCF_SD
            comparison_dict[eachcomparison]['URANUS']['mean_Arcsec_FCF'] = source_5_mean_Arcsec_FCF
            comparison_dict[eachcomparison]['URANUS']['SD_Arcsec_FCF']   = source_5_Arcsec_FCF_SD
        else:
            comparison_dict[eachcomparison]['URANUS']['ARs']                = comparison_AR[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['Beam_Eff']           = comparison_Beam_Eff[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['ARCSEC_FCF']         = comparison_Arcsec_FCF[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['PEAK_FCF']           = comparison_Peak_FCF[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['TAU225']             = comparison_weather[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['AM']                 = comparison_AM[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['PEAK_FLUXES']        = comparison_PEAK_FLUX[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['TOTAL_FLUXES']       = comparison_TOTAL_FLUX[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['FWHM1']              = comparison_FWHM1[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['FWHM2']              = comparison_FWHM2[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['MJD']                = comparison_MJD[source_5_ind]
            comparison_dict[eachcomparison]['URANUS']['mean_Beam_eff']      = np.nan
            comparison_dict[eachcomparison]['URANUS']['SD_Beam_eff']        = np.nan
            comparison_dict[eachcomparison]['URANUS']['mean_FWHM1']         = np.nan
            comparison_dict[eachcomparison]['URANUS']['SD_FWHM1']           = np.nan
            comparison_dict[eachcomparison]['URANUS']['mean_FWHM2']         = np.nan
            comparison_dict[eachcomparison]['URANUS']['SD_FWHM2']           = np.nan
            comparison_dict[eachcomparison]['URANUS']['mean_AR']         = np.nan
            comparison_dict[eachcomparison]['URANUS']['SD_AR']           = np.nan
            comparison_dict[eachcomparison]['URANUS']['mean_Peak_FCF']   = np.nan
            comparison_dict[eachcomparison]['URANUS']['SD_Peak_FCF']     = np.nan
            comparison_dict[eachcomparison]['URANUS']['mean_Arcsec_FCF'] = np.nan
            comparison_dict[eachcomparison]['URANUS']['SD_Arcsec_FCF']   = np.nan

        if len(source_6_ind)>3:
            source_6_mean_AR         = np.mean(comparison_AR[source_6_ind])
            source_6_AR_SD           = np.std(comparison_AR[source_6_ind],ddof=1)
            source_6_mean_Peak_FCF   = np.mean(comparison_Peak_FCF[source_6_ind])*ARP220_Peak_fix
            source_6_Peak_FCF_SD     = np.std(comparison_Peak_FCF[source_6_ind],ddof=1)*ARP220_Peak_fix
            source_6_mean_Arcsec_FCF = np.mean(comparison_Arcsec_FCF[source_6_ind])*ARP220_Total_fix
            source_6_Arcsec_FCF_SD   = np.std(comparison_Arcsec_FCF[source_6_ind],ddof=1)*ARP220_Total_fix
            source_6_mean_FWHM1         = np.mean(comparison_FWHM1[source_6_ind])
            source_6_FWHM1_SD           = np.std(comparison_FWHM1[source_6_ind],ddof=1)
            source_6_mean_FWHM2         = np.mean(comparison_FWHM2[source_6_ind])
            source_6_FWHM2_SD           = np.std(comparison_FWHM2[source_6_ind],ddof=1)
            source_6_mean_Beam_Eff      = np.mean(comparison_Beam_Eff[source_6_ind])
            source_6_Beam_Eff_SD        = np.std(comparison_Beam_Eff[source_6_ind],ddof=1)
            comparison_dict[eachcomparison]['Arp220']['ARs']                = comparison_AR[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['Beam_Eff']           = comparison_Beam_Eff[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['ARCSEC_FCF']         = comparison_Arcsec_FCF[source_6_ind]*ARP220_Total_fix
            comparison_dict[eachcomparison]['Arp220']['PEAK_FCF']           = comparison_Peak_FCF[source_6_ind]*ARP220_Peak_fix
            comparison_dict[eachcomparison]['Arp220']['TAU225']             = comparison_weather[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['AM']                 = comparison_AM[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['PEAK_FLUXES']        = comparison_PEAK_FLUX[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['TOTAL_FLUXES']       = comparison_TOTAL_FLUX[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['FWHM1']              = comparison_FWHM1[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['FWHM2']              = comparison_FWHM2[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['MJD']                = comparison_MJD[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['mean_Beam_eff']      = source_6_mean_Beam_Eff
            comparison_dict[eachcomparison]['Arp220']['SD_Beam_eff']        = source_6_Beam_Eff_SD
            comparison_dict[eachcomparison]['Arp220']['mean_FWHM1']         = source_6_mean_FWHM1
            comparison_dict[eachcomparison]['Arp220']['SD_FWHM1']           = source_6_FWHM1_SD
            comparison_dict[eachcomparison]['Arp220']['mean_FWHM2']         = source_6_mean_FWHM2
            comparison_dict[eachcomparison]['Arp220']['SD_FWHM2']           = source_6_FWHM2_SD
            comparison_dict[eachcomparison]['Arp220']['mean_AR']         = source_6_mean_AR
            comparison_dict[eachcomparison]['Arp220']['SD_AR']           = source_6_AR_SD
            comparison_dict[eachcomparison]['Arp220']['mean_Peak_FCF']   = source_6_mean_Peak_FCF
            comparison_dict[eachcomparison]['Arp220']['SD_Peak_FCF']     = source_6_Peak_FCF_SD
            comparison_dict[eachcomparison]['Arp220']['mean_Arcsec_FCF'] = source_6_mean_Arcsec_FCF
            comparison_dict[eachcomparison]['Arp220']['SD_Arcsec_FCF']   = source_6_Arcsec_FCF_SD
        else:
            comparison_dict[eachcomparison]['Arp220']['ARs']                = comparison_AR[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['Beam_Eff']           = comparison_Beam_Eff[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['ARCSEC_FCF']         = comparison_Arcsec_FCF[source_6_ind]*ARP220_Total_fix
            comparison_dict[eachcomparison]['Arp220']['PEAK_FCF']           = comparison_Peak_FCF[source_6_ind]*ARP220_Peak_fix
            comparison_dict[eachcomparison]['Arp220']['TAU225']             = comparison_weather[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['AM']                 = comparison_AM[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['PEAK_FLUXES']        = comparison_PEAK_FLUX[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['TOTAL_FLUXES']       = comparison_TOTAL_FLUX[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['FWHM1']              = comparison_FWHM1[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['FWHM2']              = comparison_FWHM2[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['MJD']                = comparison_MJD[source_6_ind]
            comparison_dict[eachcomparison]['Arp220']['mean_Beam_eff']      = np.nan
            comparison_dict[eachcomparison]['Arp220']['SD_Beam_eff']        = np.nan
            comparison_dict[eachcomparison]['Arp220']['mean_FWHM1']         = np.nan
            comparison_dict[eachcomparison]['Arp220']['SD_FWHM1']           = np.nan
            comparison_dict[eachcomparison]['Arp220']['mean_FWHM2']         = np.nan
            comparison_dict[eachcomparison]['Arp220']['SD_FWHM2']           = np.nan
            comparison_dict[eachcomparison]['Arp220']['mean_AR']         = np.nan
            comparison_dict[eachcomparison]['Arp220']['SD_AR']           = np.nan
            comparison_dict[eachcomparison]['Arp220']['mean_Peak_FCF']   = np.nan
            comparison_dict[eachcomparison]['Arp220']['SD_Peak_FCF']     = np.nan
            comparison_dict[eachcomparison]['Arp220']['mean_Arcsec_FCF'] = np.nan
            comparison_dict[eachcomparison]['Arp220']['SD_Arcsec_FCF']   = np.nan


    else:
        if len(comparison_Peak_FCF)>3:
            allsources_mean_AR         = np.mean(comparison_AR)
            allsources_AR_SD           = np.std(comparison_AR,ddof=1)
            allsources_mean_Peak_FCF   = np.mean(comparison_Peak_FCF)
            allsources_Peak_FCF_SD     = np.std(comparison_Peak_FCF,ddof=1)
            allsources_mean_Arcsec_FCF = np.mean(comparison_Arcsec_FCF)
            allsources_Arcsec_FCF_SD   = np.std(comparison_Arcsec_FCF,ddof=1)
            allsources_mean_FWHM1         = np.mean(comparison_FWHM1[allsources_ind])
            allsources_FWHM1_SD           = np.std(comparison_FWHM1[allsources_ind],ddof=1)
            allsources_mean_FWHM2         = np.mean(comparison_FWHM2[allsources_ind])
            allsources_FWHM2_SD           = np.std(comparison_FWHM2[allsources_ind],ddof=1)
            allsources_mean_Beam_Eff      = np.mean(comparison_Beam_Eff[allsources_ind])
            allsources_Beam_Eff_SD        = np.std(comparison_Beam_Eff[allsources_ind],ddof=1)
            comparison_dict[eachcomparison]['ALLSOURCES']['ARs']                = comparison_AR[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['Beam_Eff']           = comparison_Beam_Eff[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['ARCSEC_FCF']         = comparison_Arcsec_FCF[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['PEAK_FCF']           = comparison_Peak_FCF[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['TAU225']             = comparison_weather[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['AM']                 = comparison_AM[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['PEAK_FLUXES']        = comparison_PEAK_FLUX[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['TOTAL_FLUXES']       = comparison_TOTAL_FLUX[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['FWHM1']              = comparison_FWHM1[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['FWHM2']              = comparison_FWHM2[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['MJD']                = comparison_MJD[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['mean_Beam_eff']      = allsources_mean_Beam_Eff
            comparison_dict[eachcomparison]['ALLSOURCES']['SD_Beam_eff']        = allsources_Beam_Eff_SD
            comparison_dict[eachcomparison]['ALLSOURCES']['mean_FWHM1']         = allsources_mean_FWHM1
            comparison_dict[eachcomparison]['ALLSOURCES']['SD_FWHM1']           = allsources_FWHM1_SD
            comparison_dict[eachcomparison]['ALLSOURCES']['mean_FWHM2']         = allsources_mean_FWHM2
            comparison_dict[eachcomparison]['ALLSOURCES']['SD_FWHM2']           = allsources_FWHM2_SD
            comparison_dict[eachcomparison]['ALLSOURCES']['mean_AR']         = allsources_mean_AR
            comparison_dict[eachcomparison]['ALLSOURCES']['SD_AR']           = allsources_AR_SD
            comparison_dict[eachcomparison]['ALLSOURCES']['mean_Peak_FCF']   = allsources_mean_Peak_FCF
            comparison_dict[eachcomparison]['ALLSOURCES']['SD_Peak_FCF']     = allsources_Peak_FCF_SD
            comparison_dict[eachcomparison]['ALLSOURCES']['mean_Arcsec_FCF'] = allsources_mean_Arcsec_FCF
            comparison_dict[eachcomparison]['ALLSOURCES']['SD_Arcsec_FCF']   = allsources_Arcsec_FCF_SD
        else:
            comparison_dict[eachcomparison]['ALLSOURCES']['ARs']                = comparison_AR[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['Beam_Eff']           = comparison_Beam_Eff[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['TAU225']             = comparison_weather[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['AM']                 = comparison_AM[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['PEAK_FLUXES']        = comparison_PEAK_FLUX[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['TOTAL_FLUXES']       = comparison_TOTAL_FLUX[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['FWHM1']              = comparison_FWHM1[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['FWHM2']              = comparison_FWHM2[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['MJD']                = comparison_MJD[allsources_ind]
            comparison_dict[eachcomparison]['ALLSOURCES']['mean_Beam_eff']      = np.nan
            comparison_dict[eachcomparison]['ALLSOURCES']['SD_Beam_eff']        = np.nan
            comparison_dict[eachcomparison]['ALLSOURCES']['mean_FWHM1']         = np.nan
            comparison_dict[eachcomparison]['ALLSOURCES']['SD_FWHM1']           = np.nan
            comparison_dict[eachcomparison]['ALLSOURCES']['mean_FWHM2']         = np.nan
            comparison_dict[eachcomparison]['ALLSOURCES']['SD_FWHM2']           = np.nan
            comparison_dict[eachcomparison]['ALLSOURCES']['mean_AR']         = np.nan
            comparison_dict[eachcomparison]['ALLSOURCES']['SD_AR']           = np.nan
            comparison_dict[eachcomparison]['ALLSOURCES']['mean_Peak_FCF']   = np.nan
            comparison_dict[eachcomparison]['ALLSOURCES']['SD_Peak_FCF']     = np.nan
            comparison_dict[eachcomparison]['ALLSOURCES']['mean_Arcsec_FCF'] = np.nan
            comparison_dict[eachcomparison]['ALLSOURCES']['SD_Arcsec_FCF']   = np.nan

nominal_dict = {}
if individual_source_plots:
        nominal_dict['CRL2688']    = {}
        nominal_dict['CRL618']     = {}
        nominal_dict['MARS']       = {}
        nominal_dict['NEPTUNE']    = {}
        nominal_dict['URANUS']     = {}
        nominal_dict['Arp220']     = {}
else:
        nominal_dict['ALLSOURCES'] = {}

UT         = []
DATE       = []
MJD        = []
AR         = []
Peak_FCF   = []
Arcsec_FCF = []
SOURCE     = []
weather    = []
AM         = []
TauTimesAM = []
PEAK_FLUX  = []
TOTAL_FLUX = []
FWHM1      = []
FWHM2      = []

for i in range(len(cat['MJDST'])):
    if cat['MJDST'][i]>=MJDstart and cat['MJDST'][i]<=MJDend:
#        if cat['MJDST'][i] < SILMJDend or cat['MJDST'][i] > BLAMJDstart: # No reliable WVM here
#            if cat['MJDST'][i]< WVM_out_of_service or cat['MJDST'][i] > WVM_back_in_service:
                if datetime.date(Time(cat['MJDST'][i],format='mjd').datetime.year,Time(cat['MJDST'][i],format='mjd').datetime.month,Time(cat['MJDST'][i],format='mjd').datetime.day) >= date_constraint_Obs_after_UT and datetime.date(Time(cat['MJDST'][i],format='mjd').datetime.year,Time(cat['MJDST'][i],format='mjd').datetime.month,Time(cat['MJDST'][i],format='mjd').datetime.day) <= date_constraint_Obs_before_UT:
                    if datetime.time(Time(cat['MJDST'][i],format='mjd').datetime.hour,Time(cat['MJDST'][i],format='mjd').datetime.minute,Time(cat['MJDST'][i],format='mjd').datetime.second) >= Time_Constraint_Obs_after_UT and datetime.time(Time(cat['MJDST'][i],format='mjd').datetime.hour,Time(cat['MJDST'][i],format='mjd').datetime.minute,Time(cat['MJDST'][i],format='mjd').datetime.second) <= Time_Constraint_Obs_before_UT:
                        if cat['ARCSEC_FCF'][i] < ARCSEC_UPPER_LIM and cat['ARCSEC_FCF'][i] > ARCSEC_LOWER_LIM:
                            if cat['PEAK_FCF'][i] < PEAK_UPPER_LIM and cat['PEAK_FCF'][i] > PEAK_LOWER_LIM:
                                MJD.append(cat['MJDST'][i])
                                UT.append(str(Time(cat['MJDST'][i],format='mjd').isot).split('T')[-1])
                                DATE.append(cat['MJDST'][i])
                                AR.append(max(cat['FWHM1'][i],cat['FWHM2'][i])/min(cat['FWHM1'][i],cat['FWHM2'][i]))
                                Peak_FCF.append(cat['PEAK_FCF'][i])
                                Arcsec_FCF.append(cat['ARCSEC_FCF'][i])
                                SOURCE.append(cat['SOURCE_ID'][i])
                                weather.append((cat['WVMTAUST'][i]+cat['WVMTAUEN'][i])/2.0)
                                AM.append((cat['AMSTART'][i]+cat['AMEND'][i])/2.0)
                                TauTimesAM.append(((cat['WVMTAUST'][i]+cat['WVMTAUEN'][i])/2.0)*((cat['AMSTART'][i]+cat['AMEND'][i])/2.0))
                                PEAK_FLUX.append(cat['PEAK_FLUX'][i])
                                TOTAL_FLUX.append(cat['TOTAL_FLUX'][i])
                                FWHM1.append(cat['FWHM1'][i])
                                FWHM2.append(cat['FWHM2'][i])

HST = []
for eachUT in UT:
    HSThour = int(eachUT.split(':')[0])-10
    if HSThour<0:
        HSThour = HSThour+24
    HST.append(str(HSThour)+':'+eachUT.split(':')[1]+':'+eachUT.split(':')[2].split('.')[0])

UT             = np.array(UT)
DATE           = np.array(DATE)
MJD            = np.array(MJD)
AR             = np.array(AR)
Peak_FCF       = np.array(Peak_FCF)
Arcsec_FCF     = np.array(Arcsec_FCF)
SOURCE         = np.array(SOURCE)
# Switch up Arp220 source IDs - I know this is stupid, but it is what it is.
SOURCE[np.where(SOURCE==0)] = 6
weather        = np.array(weather)
AM             = np.array(AM)
HST            = np.array(HST)
TauTimesAM     = np.array(TauTimesAM)
COR_PEAK_FLUX  = np.array(PEAK_FLUX)
COR_TOTAL_FLUX = np.array(TOTAL_FLUX)
FWHM1          = np.array(FWHM1)
FWHM2          = np.array(FWHM2)

source_1_ind = np.where(SOURCE==1)[0]
source_2_ind = np.where(SOURCE==2)[0]
source_3_ind = np.where(SOURCE==3)[0]
source_4_ind = np.where(SOURCE==4)[0]
source_5_ind = np.where(SOURCE==5)[0]
source_6_ind = np.where(SOURCE==6)[0]

if individual_source_plots:
    if len(source_1_ind)>3:
        source_1_mean_AR         = np.mean(AR[source_1_ind])
        source_1_AR_SD           = np.std(AR[source_1_ind],ddof=1)
        source_1_mean_Peak_FCF   = np.mean(Peak_FCF[source_1_ind])
        source_1_Peak_FCF_SD     = np.std(Peak_FCF[source_1_ind],ddof=1)
        source_1_mean_Arcsec_FCF = np.mean(Arcsec_FCF[source_1_ind])
        source_1_Arcsec_FCF_SD   = np.std(Arcsec_FCF[source_1_ind],ddof=1)
        source_1_mean_FWHM1         = np.mean(FWHM1[source_1_ind])
        source_1_FWHM1_SD           = np.std(FWHM1[source_1_ind],ddof=1)
        source_1_mean_FWHM2         = np.mean(FWHM2[source_1_ind])
        source_1_FWHM2_SD           = np.std(FWHM2[source_1_ind],ddof=1)
        nominal_dict['CRL2688']['mean_FWHM1']         = source_1_mean_FWHM1
        nominal_dict['CRL2688']['SD_FWHM1']           = source_1_FWHM1_SD
        nominal_dict['CRL2688']['mean_FWHM2']         = source_1_mean_FWHM2
        nominal_dict['CRL2688']['SD_FWHM2']           = source_1_FWHM2_SD
        nominal_dict['CRL2688']['mean_AR']         = source_1_mean_AR
        nominal_dict['CRL2688']['SD_AR']           = source_1_AR_SD
        nominal_dict['CRL2688']['mean_Peak_FCF']   = source_1_mean_Peak_FCF
        nominal_dict['CRL2688']['SD_Peak_FCF']     = source_1_Peak_FCF_SD
        nominal_dict['CRL2688']['mean_Arcsec_FCF'] = source_1_mean_Arcsec_FCF
        nominal_dict['CRL2688']['SD_Arcsec_FCF']   = source_1_Arcsec_FCF_SD
    else:
        nominal_dict['CRL2688']['mean_FWHM1']         = np.nan 
        nominal_dict['CRL2688']['SD_FWHM1']           = np.nan
        nominal_dict['CRL2688']['mean_FWHM2']         = np.nan
        nominal_dict['CRL2688']['SD_FWHM2']           = np.nan
        nominal_dict['CRL2688']['mean_AR']         = np.nan
        nominal_dict['CRL2688']['SD_AR']           = np.nan
        nominal_dict['CRL2688']['mean_Peak_FCF']   = np.nan
        nominal_dict['CRL2688']['SD_Peak_FCF']     = np.nan
        nominal_dict['CRL2688']['mean_Arcsec_FCF'] = np.nan
        nominal_dict['CRL2688']['SD_Arcsec_FCF']   = np.nan

    if len(source_2_ind)>3:
        source_2_mean_AR         = np.mean(AR[source_2_ind])
        source_2_AR_SD           = np.std(AR[source_2_ind],ddof=1)
        source_2_mean_Peak_FCF   = np.mean(Peak_FCF[source_2_ind])
        source_2_Peak_FCF_SD     = np.std(Peak_FCF[source_2_ind],ddof=1)
        source_2_mean_Arcsec_FCF = np.mean(Arcsec_FCF[source_2_ind])
        source_2_Arcsec_FCF_SD   = np.std(Arcsec_FCF[source_2_ind],ddof=1)
        source_2_mean_FWHM1         = np.mean(FWHM1[source_2_ind])
        source_2_FWHM1_SD           = np.std(FWHM1[source_2_ind],ddof=1)
        source_2_mean_FWHM2         = np.mean(FWHM2[source_2_ind])
        source_2_FWHM2_SD           = np.std(FWHM2[source_2_ind],ddof=1)
        nominal_dict['CRL618']['mean_FWHM1']         = source_2_mean_FWHM1
        nominal_dict['CRL618']['SD_FWHM1']           = source_2_FWHM1_SD
        nominal_dict['CRL618']['mean_FWHM2']         = source_2_mean_FWHM2
        nominal_dict['CRL618']['SD_FWHM2']           = source_2_FWHM2_SD
        nominal_dict['CRL618']['mean_AR']         = source_2_mean_AR
        nominal_dict['CRL618']['SD_AR']           = source_2_AR_SD
        nominal_dict['CRL618']['mean_Peak_FCF']   = source_2_mean_Peak_FCF
        nominal_dict['CRL618']['SD_Peak_FCF']     = source_2_Peak_FCF_SD
        nominal_dict['CRL618']['mean_Arcsec_FCF'] = source_2_mean_Arcsec_FCF
        nominal_dict['CRL618']['SD_Arcsec_FCF']   = source_2_Arcsec_FCF_SD
    else:
        nominal_dict['CRL618']['mean_FWHM1']         = np.nan
        nominal_dict['CRL618']['SD_FWHM1']           = np.nan
        nominal_dict['CRL618']['mean_FWHM2']         = np.nan
        nominal_dict['CRL618']['SD_FWHM2']           = np.nan 
        nominal_dict['CRL618']['mean_AR']         = np.nan
        nominal_dict['CRL618']['SD_AR']           = np.nan
        nominal_dict['CRL618']['mean_Peak_FCF']   = np.nan
        nominal_dict['CRL618']['SD_Peak_FCF']     = np.nan
        nominal_dict['CRL618']['mean_Arcsec_FCF'] = np.nan
        nominal_dict['CRL618']['SD_Arcsec_FCF']   = np.nan    

    if len(source_3_ind)>3:
        source_3_mean_AR         = np.mean(AR[source_3_ind])
        source_3_AR_SD           = np.std(AR[source_3_ind],ddof=1)
        source_3_mean_Peak_FCF   = np.mean(Peak_FCF[source_3_ind])
        source_3_Peak_FCF_SD     = np.std(Peak_FCF[source_3_ind],ddof=1)
        source_3_mean_Arcsec_FCF = np.mean(Arcsec_FCF[source_3_ind])
        source_3_Arcsec_FCF_SD   = np.std(Arcsec_FCF[source_3_ind],ddof=1)
        source_3_mean_FWHM1         = np.mean(FWHM1[source_3_ind])
        source_3_FWHM1_SD           = np.std(FWHM1[source_3_ind],ddof=1)
        source_3_mean_FWHM2         = np.mean(FWHM2[source_3_ind])
        source_3_FWHM2_SD           = np.std(FWHM2[source_3_ind],ddof=1)
        nominal_dict['MARS']['mean_FWHM1']         = source_3_mean_FWHM1
        nominal_dict['MARS']['SD_FWHM1']           = source_3_FWHM1_SD
        nominal_dict['MARS']['mean_FWHM2']         = source_3_mean_FWHM2
        nominal_dict['MARS']['SD_FWHM2']           = source_3_FWHM2_SD
        nominal_dict['MARS']['mean_AR']         = source_3_mean_AR
        nominal_dict['MARS']['SD_AR']           = source_3_AR_SD
        nominal_dict['MARS']['mean_Peak_FCF']   = source_3_mean_Peak_FCF
        nominal_dict['MARS']['SD_Peak_FCF']     = source_3_Peak_FCF_SD
        nominal_dict['MARS']['mean_Arcsec_FCF'] = source_3_mean_Arcsec_FCF
        nominal_dict['MARS']['SD_Arcsec_FCF']   = source_3_Arcsec_FCF_SD
    else:
        nominal_dict['MARS']['mean_FWHM1']         = np.nan
        nominal_dict['MARS']['SD_FWHM1']           = np.nan
        nominal_dict['MARS']['mean_FWHM2']         = np.nan
        nominal_dict['MARS']['SD_FWHM2']           = np.nan 
        nominal_dict['MARS']['mean_AR']         = np.nan
        nominal_dict['MARS']['SD_AR']           = np.nan
        nominal_dict['MARS']['mean_Peak_FCF']   = np.nan
        nominal_dict['MARS']['SD_Peak_FCF']     = np.nan
        nominal_dict['MARS']['mean_Arcsec_FCF'] = np.nan
        nominal_dict['MARS']['SD_Arcsec_FCF']   = np.nan    

    if len(source_4_ind)>3:
        source_4_mean_AR         = np.mean(AR[source_4_ind])
        source_4_AR_SD           = np.std(AR[source_4_ind],ddof=1)
        source_4_mean_Peak_FCF   = np.mean(Peak_FCF[source_4_ind])
        source_4_Peak_FCF_SD     = np.std(Peak_FCF[source_4_ind],ddof=1)
        source_4_mean_Arcsec_FCF = np.mean(Arcsec_FCF[source_4_ind])
        source_4_Arcsec_FCF_SD   = np.std(Arcsec_FCF[source_4_ind],ddof=1)
        source_4_mean_FWHM1         = np.mean(FWHM1[source_4_ind])
        source_4_FWHM1_SD           = np.std(FWHM1[source_4_ind],ddof=1)
        source_4_mean_FWHM2         = np.mean(FWHM2[source_4_ind])
        source_4_FWHM2_SD           = np.std(FWHM2[source_4_ind],ddof=1)
        nominal_dict['NEPTUNE']['mean_FWHM1']         = source_4_mean_FWHM1
        nominal_dict['NEPTUNE']['SD_FWHM1']           = source_4_FWHM1_SD
        nominal_dict['NEPTUNE']['mean_FWHM2']         = source_4_mean_FWHM2
        nominal_dict['NEPTUNE']['SD_FWHM2']           = source_4_FWHM2_SD
        nominal_dict['NEPTUNE']['mean_AR']         = source_4_mean_AR
        nominal_dict['NEPTUNE']['SD_AR']           = source_4_AR_SD
        nominal_dict['NEPTUNE']['mean_Peak_FCF']   = source_4_mean_Peak_FCF
        nominal_dict['NEPTUNE']['SD_Peak_FCF']     = source_4_Peak_FCF_SD
        nominal_dict['NEPTUNE']['mean_Arcsec_FCF'] = source_4_mean_Arcsec_FCF
        nominal_dict['NEPTUNE']['SD_Arcsec_FCF']   = source_4_Arcsec_FCF_SD
    else:
        nominal_dict['NEPTUNE']['mean_FWHM1']         = np.nan
        nominal_dict['NEPTUNE']['SD_FWHM1']           = np.nan
        nominal_dict['NEPTUNE']['mean_FWHM2']         = np.nan
        nominal_dict['NEPTUNE']['SD_FWHM2']           = np.nan
        nominal_dict['NEPTUNE']['mean_AR']         = np.nan
        nominal_dict['NEPTUNE']['SD_AR']           = np.nan
        nominal_dict['NEPTUNE']['mean_Peak_FCF']   = np.nan
        nominal_dict['NEPTUNE']['SD_Peak_FCF']     = np.nan
        nominal_dict['NEPTUNE']['mean_Arcsec_FCF'] = np.nan
        nominal_dict['NEPTUNE']['SD_Arcsec_FCF']   = np.nan    

    if len(source_5_ind)>3:
        source_5_mean_AR         = np.mean(AR[source_5_ind])
        source_5_AR_SD           = np.std(AR[source_5_ind],ddof=1)
        source_5_mean_Peak_FCF   = np.mean(Peak_FCF[source_5_ind])
        source_5_Peak_FCF_SD     = np.std(Peak_FCF[source_5_ind],ddof=1)
        source_5_mean_Arcsec_FCF = np.mean(Arcsec_FCF[source_5_ind])
        source_5_Arcsec_FCF_SD   = np.std(Arcsec_FCF[source_5_ind],ddof=1)
        source_5_mean_FWHM1         = np.mean(FWHM1[source_5_ind])
        source_5_FWHM1_SD           = np.std(FWHM1[source_5_ind],ddof=1)
        source_5_mean_FWHM2         = np.mean(FWHM2[source_5_ind])
        source_5_FWHM2_SD           = np.std(FWHM2[source_5_ind],ddof=1)
        nominal_dict['URANUS']['mean_FWHM1']         = source_5_mean_FWHM1
        nominal_dict['URANUS']['SD_FWHM1']           = source_5_FWHM1_SD
        nominal_dict['URANUS']['mean_FWHM2']         = source_5_mean_FWHM2
        nominal_dict['URANUS']['SD_FWHM2']           = source_5_FWHM2_SD
        nominal_dict['URANUS']['mean_AR']         = source_5_mean_AR
        nominal_dict['URANUS']['SD_AR']           = source_5_AR_SD
        nominal_dict['URANUS']['mean_Peak_FCF']   = source_5_mean_Peak_FCF
        nominal_dict['URANUS']['SD_Peak_FCF']     = source_5_Peak_FCF_SD
        nominal_dict['URANUS']['mean_Arcsec_FCF'] = source_5_mean_Arcsec_FCF
        nominal_dict['URANUS']['SD_Arcsec_FCF']   = source_5_Arcsec_FCF_SD
    else:
        nominal_dict['URANUS']['mean_FWHM1']         = np.nan
        nominal_dict['URANUS']['SD_FWHM1']           = np.nan
        nominal_dict['URANUS']['mean_FWHM2']         = np.nan
        nominal_dict['URANUS']['SD_FWHM2']           = np.nan
        nominal_dict['URANUS']['mean_AR']         = np.nan
        nominal_dict['URANUS']['SD_AR']           = np.nan
        nominal_dict['URANUS']['mean_Peak_FCF']   = np.nan
        nominal_dict['URANUS']['SD_Peak_FCF']     = np.nan
        nominal_dict['URANUS']['mean_Arcsec_FCF'] = np.nan
        nominal_dict['URANUS']['SD_Arcsec_FCF']   = np.nan

    if len(source_6_ind)>3:
        source_6_mean_AR         = np.mean(AR[source_6_ind])
        source_6_AR_SD           = np.std(AR[source_6_ind],ddof=1)
        source_6_mean_Peak_FCF   = np.mean(Peak_FCF[source_6_ind])*ARP220_Peak_fix
        source_6_Peak_FCF_SD     = np.std(Peak_FCF[source_6_ind],ddof=1)*ARP220_Peak_fix
        source_6_mean_Arcsec_FCF = np.mean(Arcsec_FCF[source_6_ind])*ARP220_Total_fix
        source_6_Arcsec_FCF_SD   = np.std(Arcsec_FCF[source_6_ind],ddof=1)*ARP220_Total_fix
        source_6_mean_FWHM1         = np.mean(FWHM1[source_6_ind])
        source_6_FWHM1_SD           = np.std(FWHM1[source_6_ind],ddof=1)
        source_6_mean_FWHM2         = np.mean(FWHM2[source_6_ind])
        source_6_FWHM2_SD           = np.std(FWHM2[source_6_ind],ddof=1)
        nominal_dict['Arp220']['mean_FWHM1']         = source_6_mean_FWHM1
        nominal_dict['Arp220']['SD_FWHM1']           = source_6_FWHM1_SD
        nominal_dict['Arp220']['mean_FWHM2']         = source_6_mean_FWHM2
        nominal_dict['Arp220']['SD_FWHM2']           = source_6_FWHM2_SD
        nominal_dict['Arp220']['mean_AR']            = source_6_mean_AR
        nominal_dict['Arp220']['SD_AR']              = source_6_AR_SD
        nominal_dict['Arp220']['mean_Peak_FCF']      = source_6_mean_Peak_FCF
        nominal_dict['Arp220']['SD_Peak_FCF']        = source_6_Peak_FCF_SD
        nominal_dict['Arp220']['mean_Arcsec_FCF']    = source_6_mean_Arcsec_FCF
        nominal_dict['Arp220']['SD_Arcsec_FCF']      = source_6_Arcsec_FCF_SD
    else:
        nominal_dict['Arp220']['mean_FWHM1']         = np.nan
        nominal_dict['Arp220']['SD_FWHM1']           = np.nan
        nominal_dict['Arp220']['mean_FWHM2']         = np.nan
        nominal_dict['Arp220']['SD_FWHM2']           = np.nan
        nominal_dict['Arp220']['mean_AR']            = np.nan
        nominal_dict['Arp220']['SD_AR']              = np.nan
        nominal_dict['Arp220']['mean_Peak_FCF']      = np.nan
        nominal_dict['Arp220']['SD_Peak_FCF']        = np.nan
        nominal_dict['Arp220']['mean_Arcsec_FCF']    = np.nan
        nominal_dict['Arp220']['SD_Arcsec_FCF']      = np.nan

else:
    if len(Peak_FCF)>3:
        allsources_mean_AR         = np.mean(AR)
        allsources_AR_SD           = np.std(AR,ddof=1)
        allsources_mean_Peak_FCF   = np.mean(Peak_FCF)
        allsources_Peak_FCF_SD     = np.std(Peak_FCF,ddof=1)
        allsources_mean_Arcsec_FCF = np.mean(Arcsec_FCF)
        allsources_Arcsec_FCF_SD   = np.std(Arcsec_FCF,ddof=1)
        allsources_mean_FWHM1         = np.mean(FWHM1[allsources_ind])
        allsources_FWHM1_SD           = np.std(FWHM1[allsources_ind],ddof=1)
        allsources_mean_FWHM2         = np.mean(FWHM2[allsources_ind])
        allsources_FWHM2_SD           = np.std(FWHM2[allsources_ind],ddof=1)
        nominal_dict['ALLSOURCES']['mean_FWHM1']         = allsources_mean_FWHM1
        nominal_dict['ALLSOURCES']['SD_FWHM1']           = allsources_FWHM1_SD
        nominal_dict['ALLSOURCES']['mean_FWHM2']         = allsources_mean_FWHM2
        nominal_dict['ALLSOURCES']['SD_FWHM2']           = allsources_FWHM2_SD
        nominal_dict['ALLSOURCES']['mean_AR']         = allsources_mean_AR
        nominal_dict['ALLSOURCES']['SD_AR']           = allsources_AR_SD
        nominal_dict['ALLSOURCES']['mean_Peak_FCF']   = allsources_mean_Peak_FCF
        nominal_dict['ALLSOURCES']['SD_Peak_FCF']     = allsources_Peak_FCF_SD
        nominal_dict['ALLSOURCES']['mean_Arcsec_FCF'] = allsources_mean_Arcsec_FCF
        nominal_dict['ALLSOURCES']['SD_Arcsec_FCF']   = allsources_Arcsec_FCF_SD
    else:
        nominal_dict['ALLSOURCES']['mean_AR']         = np.nan
        nominal_dict['ALLSOURCES']['SD_AR']           = np.nan
        nominal_dict['ALLSOURCES']['mean_Peak_FCF']   = np.nan
        nominal_dict['ALLSOURCES']['SD_Peak_FCF']     = np.nan
        nominal_dict['ALLSOURCES']['mean_Arcsec_FCF'] = np.nan
        nominal_dict['ALLSOURCES']['SD_Arcsec_FCF']   = np.nan

unique_sources = []
source_names   = []
colours        = []
for i in SOURCE:
    if i not in unique_sources:
        unique_sources.append(i)
        if i == 1:
            source_names.append('CRL 2688')
            colours.append('#003593')
        if i == 2:
            source_names.append('CRL 618')
            colours.append('#20b2aa')
        if i == 3:
            source_names.append('Mars')
            colours.append('#800000')
        if i == 4:
            source_names.append('Neptune')
            colours.append('#daa520')
        if i == 5:
            source_names.append('Uranus')
            colours.append('#f6546a')
        if i == 6:
            source_names.append('Arp220')
            #colours.append('#FFDEAD')
            colours.append('#885FCD')


if tick_spacing < 300: # If tick spacing is less than 5 hours (300 minutes), show time in HST, otherwise, just show date
    if individual_source_plots:
        HSTticks_source1 = []
        HSTticks_source2 = []
        HSTticks_source3 = []
        HSTticks_source4 = []
        HSTticks_source5 = []
        HSTticks_source6 = []
        
        for eachsource in source_names:
            if eachsource == 'CRL 2688':
                for i in np.arange(min(MJD[np.where(SOURCE == 1)]),max(MJD[np.where(SOURCE == 1)])+tick_spacing_mjd/2.0,tick_spacing_mjd):
                    UTtime = str(Time(i,format='mjd').isot).split('T')[-1]
                    HSThour = int(UTtime.split(':')[0])-10
                    if HSThour<0:
                        HSThour = HSThour+24
                    HSTticks_source1.append(str(HSThour)+':'+UTtime.split(':')[1])    
            if eachsource == 'CRL 618':
                for i in np.arange(min(MJD[np.where(SOURCE == 2)]),max(MJD[np.where(SOURCE == 2)])+tick_spacing_mjd/2.0,tick_spacing_mjd):
                    UTtime = str(Time(i,format='mjd').isot).split('T')[-1]
                    HSThour = int(UTtime.split(':')[0])-10
                    if HSThour<0:
                        HSThour = HSThour+24
                    HSTticks_source2.append(str(HSThour)+':'+UTtime.split(':')[1])
            if eachsource == 'Mars':
                for i in np.arange(min(MJD[np.where(SOURCE == 3)]),max(MJD[np.where(SOURCE == 3)])+tick_spacing_mjd/2.0,tick_spacing_mjd):
                    UTtime = str(Time(i,format='mjd').isot).split('T')[-1]
                    HSThour = int(UTtime.split(':')[0])-10
                    if HSThour<0:
                        HSThour = HSThour+24
                    HSTticks_source3.append(str(HSThour)+':'+UTtime.split(':')[1])
            if eachsource == 'Neptune':
                for i in np.arange(min(MJD[np.where(SOURCE == 4)]),max(MJD[np.where(SOURCE == 4)])+tick_spacing_mjd/2.0,tick_spacing_mjd):
                    UTtime = str(Time(i,format='mjd').isot).split('T')[-1]
                    HSThour = int(UTtime.split(':')[0])-10
                    if HSThour<0:
                        HSThour = HSThour+24
                    HSTticks_source4.append(str(HSThour)+':'+UTtime.split(':')[1])
            if eachsource == 'Uranus':
                for i in np.arange(min(MJD[np.where(SOURCE == 5)]),max(MJD[np.where(SOURCE == 5)])+tick_spacing_mjd/2.0,tick_spacing_mjd):
                    UTtime = str(Time(i,format='mjd').isot).split('T')[-1]
                    HSThour = int(UTtime.split(':')[0])-10
                    if HSThour<0:
                        HSThour = HSThour+24
                    HSTticks_source5.append(str(HSThour)+':'+UTtime.split(':')[1])
            if eachsource == 'Arp220':
                for i in np.arange(min(MJD[np.where(SOURCE == 6)]),max(MJD[np.where(SOURCE == 6)])+tick_spacing_mjd/2.0,tick_spacing_mjd):
                    UTtime = str(Time(i,format='mjd').isot).split('T')[-1]
                    HSThour = int(UTtime.split(':')[0])-10
                    if HSThour<0:
                        HSThour = HSThour+24
                    HSTticks_source6.append(str(HSThour)+':'+UTtime.split(':')[1])
    else:
        HSTticks_allsources = []
        for i in np.arange(min(MJD),max(MJD)+tick_spacing_mjd/2.0,tick_spacing_mjd):
            UTtime = str(Time(i,format='mjd').isot).split('T')[-1]
            HSThour = int(UTtime.split(':')[0])-10
            if HSThour<0:
                HSThour = HSThour+24
            HSTticks_allsources.append(str(HSThour)+':'+UTtime.split(':')[1])

else:

    if individual_source_plots:
        HSTticks_source1 = []
        HSTticks_source2 = []
        HSTticks_source3 = []
        HSTticks_source4 = []
        HSTticks_source5 = []
        HSTticks_source6 = []
        for eachsource in source_names:
            if eachsource == 'CRL 2688':
                for i in np.arange(min(MJD[np.where(SOURCE == 1)]),max(MJD[np.where(SOURCE == 1)])+tick_spacing_mjd/2.0,tick_spacing_mjd):
                    UTdate = str(Time(i,format='mjd').isot).split('T')[0]
                    HSTticks_source1.append(UTdate)
        for eachsource in source_names:
            if eachsource == 'CRL 618':
                for i in np.arange(min(MJD[np.where(SOURCE == 2)]),max(MJD[np.where(SOURCE == 2)])+tick_spacing_mjd/2.0,tick_spacing_mjd):
                    UTdate = str(Time(i,format='mjd').isot).split('T')[0]
                    HSTticks_source2.append(UTdate)
        for eachsource in source_names:
            if eachsource == 'Mars':
                for i in np.arange(min(MJD[np.where(SOURCE == 3)]),max(MJD[np.where(SOURCE == 3)])+tick_spacing_mjd/2.0,tick_spacing_mjd):
                    UTdate = str(Time(i,format='mjd').isot).split('T')[0]
                    HSTticks_source3.append(UTdate)
        for eachsource in source_names:
            if eachsource == 'Neptune':
                for i in np.arange(min(MJD[np.where(SOURCE == 4)]),max(MJD[np.where(SOURCE == 4)])+tick_spacing_mjd/2.0,tick_spacing_mjd):
                    UTdate = str(Time(i,format='mjd').isot).split('T')[0]
                    HSTticks_source4.append(UTdate)
        for eachsource in source_names:
            if eachsource == 'Uranus':
                for i in np.arange(min(MJD[np.where(SOURCE == 5)]),max(MJD[np.where(SOURCE == 5)])+tick_spacing_mjd/2.0,tick_spacing_mjd):
                    UTdate = str(Time(i,format='mjd').isot).split('T')[0]
                    HSTticks_source5.append(UTdate)
        for eachsource in source_names:
            if eachsource == 'Arp220':
                for i in np.arange(min(MJD[np.where(SOURCE == 6)]),max(MJD[np.where(SOURCE == 6)])+tick_spacing_mjd/2.0,tick_spacing_mjd):
                    UTdate = str(Time(i,format='mjd').isot).split('T')[0]
                    HSTticks_source6.append(UTdate)
    else:
        HSTticks_allsources = []
        for i in np.arange(min(MJD),max(MJD)+tick_spacing_mjd/2.0,tick_spacing_mjd):
            UTdate = str(Time(i,format='mjd').isot).split('T')[0]
            HSTticks_allsources.append(UTdate)

if individual_source_plots:
    for i in range(len(unique_sources)):
        plt.scatter(MJD[np.where(SOURCE==unique_sources[i])],AR[np.where(SOURCE==unique_sources[i])],color=colours[i])
    
        if source_names[i] == 'CRL 2688':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['CRL2688']['mean_AR'],color=colours[i],linestyle='solid',label='CRL 2688 '+epoch)
                plt.axhline(y=nominal_dict['CRL2688']['mean_AR']+nominal_dict['CRL2688']['SD_AR'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['CRL2688']['mean_AR']-nominal_dict['CRL2688']['SD_AR'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['CRL2688']['mean_AR']-nominal_dict['CRL2688']['SD_AR'],nominal_dict['CRL2688']['mean_AR']+nominal_dict['CRL2688']['SD_AR'],alpha=0.5,color=colours[i])
            print('\nCRL2688 AR (mean,SD): ',nominal_dict['CRL2688']['mean_AR'],nominal_dict['CRL2688']['SD_AR'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison) 
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_AR']+comparison_dict[eachcomparison]['CRL2688']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_AR']-comparison_dict[eachcomparison]['CRL2688']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['CRL2688']['mean_AR']-comparison_dict[eachcomparison]['CRL2688']['SD_AR'],comparison_dict[eachcomparison]['CRL2688']['mean_AR']+comparison_dict[eachcomparison]['CRL2688']['SD_AR'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 1)]),max(MJD[np.where(SOURCE == 1)]),tick_spacing_mjd),HSTticks_source1,rotation='20')
            plt.suptitle(wave+' microns - Aspect Ratio, '+source_names[i])
    
        if source_names[i] == 'CRL 618':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['CRL618']['mean_AR'],color=colours[i],linestyle='solid',label='CRL 618 '+epoch)
                plt.axhline(y=nominal_dict['CRL618']['mean_AR']+nominal_dict['CRL618']['SD_AR'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['CRL618']['mean_AR']-nominal_dict['CRL618']['SD_AR'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['CRL618']['mean_AR']-nominal_dict['CRL618']['SD_AR'],nominal_dict['CRL618']['mean_AR']+nominal_dict['CRL618']['SD_AR'],alpha=0.5,color=colours[i])
            print('CRL618 AR (mean,SD): ',nominal_dict['CRL618']['mean_AR'],nominal_dict['CRL618']['SD_AR'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_AR']+comparison_dict[eachcomparison]['CRL618']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_AR']-comparison_dict[eachcomparison]['CRL618']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['CRL618']['mean_AR']-comparison_dict[eachcomparison]['CRL618']['SD_AR'],comparison_dict[eachcomparison]['CRL618']['mean_AR']+comparison_dict[eachcomparison]['CRL618']['SD_AR'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 2)]),max(MJD[np.where(SOURCE == 2)]),tick_spacing_mjd),HSTticks_source2,rotation='20')
            plt.suptitle(wave+' microns - Aspect Ratio, '+source_names[i])
    
        if source_names[i] == 'Mars':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['MARS']['mean_AR'],color=colours[i],linestyle='solid',label='Mars '+epoch)
                plt.axhline(y=nominal_dict['MARS']['mean_AR']+nominal_dict['MARS']['SD_AR'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['MARS']['mean_AR']-nominal_dict['MARS']['SD_AR'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['MARS']['mean_AR']-nominal_dict['MARS']['SD_AR'],nominal_dict['MARS']['mean_AR']+nominal_dict['MARS']['SD_AR'],alpha=0.5,color=colours[i])
            print('MARS AR (mean,SD): ',nominal_dict['MARS']['mean_AR'],nominal_dict['MARS']['SD_AR'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_AR']+comparison_dict[eachcomparison]['MARS']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_AR']-comparison_dict[eachcomparison]['MARS']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['MARS']['mean_AR']-comparison_dict[eachcomparison]['MARS']['SD_AR'],comparison_dict[eachcomparison]['MARS']['mean_AR']+comparison_dict[eachcomparison]['MARS']['SD_AR'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 3)]),max(MJD[np.where(SOURCE == 3)]),tick_spacing_mjd),HSTticks_source3,rotation='20')
            plt.suptitle(wave+' microns - Aspect Ratio, '+source_names[i])
    
        if source_names[i] == 'Neptune':
            try:
                if plot_nominal_epoch_mean:
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_AR'],color=colours[i],linestyle='solid',label='Neptune '+epoch)
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_AR']+nominal_dict['NEPTUNE']['SD_AR'],color=colours[i],linestyle='dotted')
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_AR']-nominal_dict['NEPTUNE']['SD_AR'],color=colours[i],linestyle='dotted')
                    plt.axhspan(nominal_dict['NEPTUNE']['mean_AR']-nominal_dict['NEPTUNE']['SD_AR'],nominal_dict['NEPTUNE']['mean_AR']+nominal_dict['NEPTUNE']['SD_AR'],alpha=0.5,color=colours[i])
                print('NEPTUNE AR (mean,SD): ',nominal_dict['NEPTUNE']['mean_AR'],nominal_dict['NEPTUNE']['SD_AR'],'\n')
                for eachcomparison in comparison_dict.keys():
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_AR']+comparison_dict[eachcomparison]['NEPTUNE']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_AR']-comparison_dict[eachcomparison]['NEPTUNE']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                    plt.axhspan(comparison_dict[eachcomparison]['NEPTUNE']['mean_AR']-comparison_dict[eachcomparison]['NEPTUNE']['SD_AR'],comparison_dict[eachcomparison]['NEPTUNE']['mean_AR']+comparison_dict[eachcomparison]['NEPTUNE']['SD_AR'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
                plt.xticks(np.arange(min(MJD[np.where(SOURCE == 4)]),max(MJD[np.where(SOURCE == 4)]),tick_spacing_mjd),HSTticks_source4,rotation='20')
                plt.suptitle(wave+' microns - Aspect Ratio, '+source_names[i])
            except NameError:
                print('\nNeptune has no previous data with which to compare!\n')
                plt.xticks(np.arange(min(MJD[np.where(SOURCE == 4)]),max(MJD[np.where(SOURCE == 4)]),tick_spacing_mjd),HSTticks_source4,rotation='20')
                plt.suptitle('Neptune')

        if source_names[i] == 'Uranus':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['URANUS']['mean_AR'],color=colours[i],linestyle='solid',label='Uranus '+epoch)
                plt.axhline(y=nominal_dict['URANUS']['mean_AR']+nominal_dict['URANUS']['SD_AR'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['URANUS']['mean_AR']-nominal_dict['URANUS']['SD_AR'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['URANUS']['mean_AR']-nominal_dict['URANUS']['SD_AR'],nominal_dict['URANUS']['mean_AR']+nominal_dict['URANUS']['SD_AR'],alpha=0.5,color=colours[i])
            print('URANUS AR (mean,SD): ',nominal_dict['URANUS']['mean_AR'],nominal_dict['URANUS']['SD_AR'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_AR']+comparison_dict[eachcomparison]['URANUS']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_AR']-comparison_dict[eachcomparison]['URANUS']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['URANUS']['mean_AR']-comparison_dict[eachcomparison]['URANUS']['SD_AR'],comparison_dict[eachcomparison]['URANUS']['mean_AR']+comparison_dict[eachcomparison]['URANUS']['SD_AR'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 5)]),max(MJD[np.where(SOURCE == 5)]),tick_spacing_mjd),HSTticks_source5,rotation='20')
            plt.suptitle(wave+' microns - Aspect Ratio, '+source_names[i])

        if source_names[i] == 'Arp220':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['Arp220']['mean_AR'],color=colours[i],linestyle='solid',label='Arp220 '+epoch)
                plt.axhline(y=nominal_dict['Arp220']['mean_AR']+nominal_dict['Arp220']['SD_AR'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['Arp220']['mean_AR']-nominal_dict['Arp220']['SD_AR'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['Arp220']['mean_AR']-nominal_dict['Arp220']['SD_AR'],nominal_dict['Arp220']['mean_AR']+nominal_dict['Arp220']['SD_AR'],alpha=0.5,color=colours[i])
            print('Arp220 AR (mean,SD): ',nominal_dict['Arp220']['mean_AR'],nominal_dict['Arp220']['SD_AR'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_AR']+comparison_dict[eachcomparison]['Arp220']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_AR']-comparison_dict[eachcomparison]['Arp220']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['Arp220']['mean_AR']-comparison_dict[eachcomparison]['Arp220']['SD_AR'],comparison_dict[eachcomparison]['Arp220']['mean_AR']+comparison_dict[eachcomparison]['Arp220']['SD_AR'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 6)]),max(MJD[np.where(SOURCE == 6)]),tick_spacing_mjd),HSTticks_source6,rotation='20')
            plt.suptitle(wave+' microns - Aspect Ratio, '+source_names[i])
    
        #plt.xticks(np.arange(min(MJD),max(MJD),tick_spacing_mjd),HSTticks,rotation='20')
        if tick_spacing<300:
            plt.xlabel('HST')
        else:
            plt.xlabel('UT Date')
        plt.ylabel('Aspect Ratio')
        plt.legend(loc='upper left')
        if saveplots == False:
            #plt.axvline(x=58448.0,color='k',linestyle='dashed',linewidth=2)
            plt.show()
        else:
            #plt.xlim(xmin=57693.0) # 57693.0 = New Filters
            plt.xlim(xmin=plt_date_min,xmax=plot_date_max) # 58484 = January 1st, 2019
            for eachvertline in vert_lines:
                plt.axvline(x=eachvertline,color='k',linestyle='solid',linewidth=2)
            plt.savefig('AR_vs_HST_'+source_names[i].replace(' ','_')+ending,dpi=300)
        plt.clf()

else:
    plt.scatter(MJD,AR,color='k')
    if plot_nominal_epoch_mean:
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_AR'],color='k',linestyle='solid',label='All Sources '+epoch)
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_AR']+nominal_dict['ALLSOURCES']['SD_AR'],color='k',linestyle='dotted')
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_AR']-nominal_dict['ALLSOURCES']['SD_AR'],color='k',linestyle='dotted')
        plt.axhspan(nominal_dict['ALLSOURCES']['mean_AR']-nominal_dict['ALLSOURCES']['SD_AR'],nominal_dict['ALLSOURCES']['mean_AR']+nominal_dict['ALLSOURCES']['SD_AR'],alpha=0.5,color='k')
    print('\nALLSOURCES AR (mean,SD): ',nominal_dict['ALLSOURCES']['mean_AR'],nominal_dict['ALLSOURCES']['SD_AR'])
    for eachcomparison in comparison_dict.keys():
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_AR']+comparison_dict[eachcomparison]['ALLSOURCES']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_AR']-comparison_dict[eachcomparison]['ALLSOURCES']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
        plt.axhspan(comparison_dict[eachcomparison]['ALLSOURCES']['mean_AR']-comparison_dict[eachcomparison]['ALLSOURCES']['SD_AR'],comparison_dict[eachcomparison]['ALLSOURCES']['mean_AR']+comparison_dict[eachcomparison]['ALLSOURCES']['SD_AR'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
    plt.xticks(np.arange(min(MJD),max(MJD),tick_spacing_mjd),HSTticks_allsources,rotation='20')
    plt.suptitle(wave+' microns - Aspect Ratio, All Sources')

    #plt.xticks(np.arange(min(MJD),max(MJD),tick_spacing_mjd),HSTticks,rotation='20')
    if tick_spacing<300:
        plt.xlabel('HST')
    else:
        plt.xlabel('UT Date')
    plt.ylabel('Aspect Ratio')
    plt.legend(loc='upper left')
    if saveplots == False:
        #plt.axvline(x=58448.0,color='k',linestyle='dashed',linewidth=2)
        plt.show()
    else:
        plt.savefig('AR_vs_HST_ALLSOURCES_'+ending,dpi=300)
    plt.clf()

plt.clf()


if individual_source_plots:

    for i in range(len(unique_sources)):
        if unique_sources[i]!=6:
            plt.scatter(MJD[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(Peak_FCF[np.where(SOURCE==unique_sources[i])]<PEAK_UPPER_LIM,Peak_FCF[np.where(SOURCE==unique_sources[i])]>PEAK_LOWER_LIM))],Peak_FCF[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(Peak_FCF[np.where(SOURCE==unique_sources[i])]<PEAK_UPPER_LIM,Peak_FCF[np.where(SOURCE==unique_sources[i])]>PEAK_LOWER_LIM))],color=colours[i])
        else: # APPLY ARP220 CORRECTION
            plt.scatter(MJD[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(Peak_FCF[np.where(SOURCE==unique_sources[i])]<PEAK_UPPER_LIM,Peak_FCF[np.where(SOURCE==unique_sources[i])]>PEAK_LOWER_LIM))],ARP220_Peak_fix*Peak_FCF[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(Peak_FCF[np.where(SOURCE==unique_sources[i])]<PEAK_UPPER_LIM,Peak_FCF[np.where(SOURCE==unique_sources[i])]>PEAK_LOWER_LIM))],color=colours[i])
    
        if source_names[i] == 'CRL 2688':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['CRL2688']['mean_Peak_FCF'],color=colours[i],linestyle='solid',label='CRL 2688 '+epoch)
                plt.axhline(y=nominal_dict['CRL2688']['mean_Peak_FCF']+nominal_dict['CRL2688']['SD_Peak_FCF'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['CRL2688']['mean_Peak_FCF']-nominal_dict['CRL2688']['SD_Peak_FCF'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['CRL2688']['mean_Peak_FCF']-nominal_dict['CRL2688']['SD_Peak_FCF'],nominal_dict['CRL2688']['mean_Peak_FCF']+nominal_dict['CRL2688']['SD_Peak_FCF'],alpha=0.5,color=colours[i])
            print('\nCRL2688 Peak_FCF (mean,SD): ',nominal_dict['CRL2688']['mean_Peak_FCF'],nominal_dict['CRL2688']['SD_Peak_FCF'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_Peak_FCF']+comparison_dict[eachcomparison]['CRL2688']['SD_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_Peak_FCF']-comparison_dict[eachcomparison]['CRL2688']['SD_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['CRL2688']['mean_Peak_FCF']-comparison_dict[eachcomparison]['CRL2688']['SD_Peak_FCF'],comparison_dict[eachcomparison]['CRL2688']['mean_Peak_FCF']+comparison_dict[eachcomparison]['CRL2688']['SD_Peak_FCF'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 1)]),max(MJD[np.where(SOURCE == 1)]),tick_spacing_mjd),HSTticks_source1,rotation='20')
            plt.suptitle(wave+' microns - Peak FCF, '+source_names[i])
    
        if source_names[i] == 'CRL 618':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['CRL618']['mean_Peak_FCF'],color=colours[i],linestyle='solid',label='CRL 618 '+epoch)
                plt.axhline(y=nominal_dict['CRL618']['mean_Peak_FCF']+nominal_dict['CRL618']['SD_Peak_FCF'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['CRL618']['mean_Peak_FCF']-nominal_dict['CRL618']['SD_Peak_FCF'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['CRL618']['mean_Peak_FCF']-nominal_dict['CRL618']['SD_Peak_FCF'],nominal_dict['CRL618']['mean_Peak_FCF']+nominal_dict['CRL618']['SD_Peak_FCF'],alpha=0.5,color=colours[i])
            print('CRL618 Peak_FCF (mean,SD): ',nominal_dict['CRL618']['mean_Peak_FCF'],nominal_dict['CRL618']['SD_Peak_FCF'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_Peak_FCF']+comparison_dict[eachcomparison]['CRL618']['SD_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_Peak_FCF']-comparison_dict[eachcomparison]['CRL618']['SD_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['CRL618']['mean_Peak_FCF']-comparison_dict[eachcomparison]['CRL618']['SD_Peak_FCF'],comparison_dict[eachcomparison]['CRL618']['mean_Peak_FCF']+comparison_dict[eachcomparison]['CRL618']['SD_Peak_FCF'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 2)]),max(MJD[np.where(SOURCE == 2)]),tick_spacing_mjd),HSTticks_source2,rotation='20')
            plt.suptitle(wave+' microns - Peak FCF, '+source_names[i])
    
        if source_names[i] == 'Mars':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['MARS']['mean_Peak_FCF'],color=colours[i],linestyle='solid',label='Mars '+epoch)
                plt.axhline(y=nominal_dict['MARS']['mean_Peak_FCF']+nominal_dict['MARS']['SD_Peak_FCF'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['MARS']['mean_Peak_FCF']-nominal_dict['MARS']['SD_Peak_FCF'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['MARS']['mean_Peak_FCF']-nominal_dict['MARS']['SD_Peak_FCF'],nominal_dict['MARS']['mean_Peak_FCF']+nominal_dict['MARS']['SD_Peak_FCF'],alpha=0.5,color=colours[i])
            print('MARS Peak_FCF (mean,SD): ',nominal_dict['MARS']['mean_Peak_FCF'],nominal_dict['MARS']['SD_Peak_FCF'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_Peak_FCF']+comparison_dict[eachcomparison]['MARS']['SD_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_Peak_FCF']-comparison_dict[eachcomparison]['MARS']['SD_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['MARS']['mean_Peak_FCF']-comparison_dict[eachcomparison]['MARS']['SD_Peak_FCF'],comparison_dict[eachcomparison]['MARS']['mean_Peak_FCF']+comparison_dict[eachcomparison]['MARS']['SD_Peak_FCF'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 3)]),max(MJD[np.where(SOURCE == 3)]),tick_spacing_mjd),HSTticks_source3,rotation='20')
            plt.suptitle(wave+' microns - Peak FCF, '+source_names[i])
    
        if source_names[i] == 'Neptune':
            try:
                if plot_nominal_epoch_mean:
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_Peak_FCF'],color=colours[i],linestyle='solid',label='Neptune '+epoch)
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_Peak_FCF']+nominal_dict['NEPTUNE']['SD_Peak_FCF'],color=colours[i],linestyle='dotted')
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_Peak_FCF']-nominal_dict['NEPTUNE']['SD_Peak_FCF'],color=colours[i],linestyle='dotted')
                    plt.axhspan(nominal_dict['NEPTUNE']['mean_Peak_FCF']-nominal_dict['NEPTUNE']['SD_Peak_FCF'],nominal_dict['NEPTUNE']['mean_Peak_FCF']+nominal_dict['NEPTUNE']['SD_Peak_FCF'],alpha=0.5,color=colours[i])
                print('NEPTUNE Peak_FCF (mean,SD): ',nominal_dict['NEPTUNE']['mean_Peak_FCF'],nominal_dict['NEPTUNE']['SD_Peak_FCF'],'\n')
                for eachcomparison in comparison_dict.keys():
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_Peak_FCF']+comparison_dict[eachcomparison]['NEPTUNE']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_Peak_FCF']-comparison_dict[eachcomparison]['NEPTUNE']['SD_AR'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                    plt.axhspan(comparison_dict[eachcomparison]['NEPTUNE']['mean_Peak_FCF']-comparison_dict[eachcomparison]['NEPTUNE']['SD_AR'],comparison_dict[eachcomparison]['NEPTUNE']['mean_Peak_FCF']+comparison_dict[eachcomparison]['NEPTUNE']['SD_AR'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
                plt.xticks(np.arange(min(MJD[np.where(SOURCE == 4)]),max(MJD[np.where(SOURCE == 4)]),tick_spacing_mjd),HSTticks_source4,rotation='20')
                plt.suptitle('Peak FCF, '+source_names[i])
            except NameError:
                print('\nNeptune has no previous data with which to compare!\n')
                plt.xticks(np.arange(min(MJD[np.where(SOURCE == 4)]),max(MJD[np.where(SOURCE == 4)]),tick_spacing_mjd),HSTticks_source4,rotation='20')
                plt.suptitle(wave+' microns - Peak FCF, '+source_names[i])
    
        if source_names[i] == 'Uranus':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['URANUS']['mean_Peak_FCF'],color=colours[i],linestyle='solid',label='Uranus '+epoch)
                plt.axhline(y=nominal_dict['URANUS']['mean_Peak_FCF']+nominal_dict['URANUS']['SD_Peak_FCF'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['URANUS']['mean_Peak_FCF']-nominal_dict['URANUS']['SD_Peak_FCF'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['URANUS']['mean_Peak_FCF']-nominal_dict['URANUS']['SD_Peak_FCF'],nominal_dict['URANUS']['mean_Peak_FCF']+nominal_dict['URANUS']['SD_Peak_FCF'],alpha=0.5,color=colours[i])
            print('URANUS Peak_FCF (mean,SD): ',nominal_dict['URANUS']['mean_Peak_FCF'],nominal_dict['URANUS']['SD_Peak_FCF'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_Peak_FCF']+comparison_dict[eachcomparison]['URANUS']['SD_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_Peak_FCF']-comparison_dict[eachcomparison]['URANUS']['SD_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['URANUS']['mean_Peak_FCF']-comparison_dict[eachcomparison]['URANUS']['SD_Peak_FCF'],comparison_dict[eachcomparison]['URANUS']['mean_Peak_FCF']+comparison_dict[eachcomparison]['URANUS']['SD_Peak_FCF'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 5)]),max(MJD[np.where(SOURCE == 5)]),tick_spacing_mjd),HSTticks_source5,rotation='20')
            plt.suptitle(wave+' microns - Peak FCF, '+source_names[i])

        if source_names[i] == 'Arp220':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['Arp220']['mean_Peak_FCF'],color=colours[i],linestyle='solid',label='Arp220 '+epoch)
                plt.axhline(y=nominal_dict['Arp220']['mean_Peak_FCF']+nominal_dict['Arp220']['SD_Peak_FCF'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['Arp220']['mean_Peak_FCF']-nominal_dict['Arp220']['SD_Peak_FCF'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['Arp220']['mean_Peak_FCF']-nominal_dict['Arp220']['SD_Peak_FCF'],nominal_dict['Arp220']['mean_Peak_FCF']+nominal_dict['Arp220']['SD_Peak_FCF'],alpha=0.5,color=colours[i])
            print('Arp220 Peak_FCF (mean,SD): ',nominal_dict['Arp220']['mean_Peak_FCF'],nominal_dict['Arp220']['SD_Peak_FCF'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_Peak_FCF']+comparison_dict[eachcomparison]['Arp220']['SD_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_Peak_FCF']-comparison_dict[eachcomparison]['Arp220']['SD_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['Arp220']['mean_Peak_FCF']-comparison_dict[eachcomparison]['Arp220']['SD_Peak_FCF'],comparison_dict[eachcomparison]['Arp220']['mean_Peak_FCF']+comparison_dict[eachcomparison]['Arp220']['SD_Peak_FCF'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 6)]),max(MJD[np.where(SOURCE == 6)]),tick_spacing_mjd),HSTticks_source6,rotation='20')
            plt.suptitle(wave+' microns - Peak FCF, '+source_names[i])
    
        if tick_spacing<300:
            plt.xlabel('HST')
        else:
            plt.xlabel('UT Date')    
    
        plt.ylabel('Peak FCF')
        #plt.xticks(np.arange(min(MJD),max(MJD),tick_spacing_mjd),HSTticks,rotation='20')
        plt.legend(loc='upper left')
        if saveplots == False:
            #plt.axvline(x=57693.0,color='k',linestyle='solid',linewidth=2) # 57693 = 2016/11/01 - New Filters
            #plt.axvline(x=58323.0,color='k',linestyle='solid',linewidth=2)  # 58323 = SMU HW Fix  (2018/07/24))
            plt.show()
        else:
            #plt.xlim(xmin=57693.0) # 58323 = 2016/11/01 - New Filters
            plt.xlim(xmin=plt_date_min,xmax=plot_date_max) # 58484 = January 1st, 2019
            #plt.xlim(xmin=57357.0,xmax=58703.0) #57357 = December 1 2015
            for eachvertline in vert_lines:
                plt.axvline(x=eachvertline,color='k',linestyle='solid',linewidth=2)
            plt.savefig('PeakFCF_vs_HST_'+source_names[i].replace(' ','_')+ending,dpi=300)   
        plt.clf()

else:
    plt.scatter(MJD[np.where(np.logical_and(Peak_FCF<PEAK_UPPER_LIM,Peak_FCF>PEAK_LOWER_LIM))],Peak_FCF[np.where(np.logical_and(Peak_FCF<PEAK_UPPER_LIM,Peak_FCF>PEAK_LOWER_LIM))],color='k')
    if plot_nominal_epoch_mean:
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_Peak_FCF'],color='k',linestyle='solid',label='All Sources '+epoch)
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_Peak_FCF']+nominal_dict['ALLSOURCES']['SD_Peak_FCF'],color='k',linestyle='dotted')
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_Peak_FCF']-nominal_dict['ALLSOURCES']['SD_Peak_FCF'],color='k',linestyle='dotted')
        plt.axhspan(nominal_dict['ALLSOURCES']['mean_Peak_FCF']-nominal_dict['ALLSOURCES']['SD_Peak_FCF'],nominal_dict['ALLSOURCES']['mean_Peak_FCF']+nominal_dict['ALLSOURCES']['SD_Peak_FCF'],alpha=0.5,color='k')
    print('ALLSOURCES Peak_FCF (mean,SD): ',nominal_dict['ALLSOURCES']['mean_Peak_FCF'],nominal_dict['ALLSOURCES']['SD_Peak_FCF'])
    for eachcomparison in comparison_dict.keys():
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_Peak_FCF']+comparison_dict[eachcomparison]['ALLSOURCES']['SD_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_Peak_FCF']-comparison_dict[eachcomparison]['ALLSOURCES']['SD_Peak_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
        plt.axhspan(comparison_dict[eachcomparison]['ALLSOURCES']['mean_Peak_FCF']-comparison_dict[eachcomparison]['ALLSOURCES']['SD_Peak_FCF'],comparison_dict[eachcomparison]['ALLSOURCES']['mean_Peak_FCF']+comparison_dict[eachcomparison]['ALLSOURCES']['SD_Peak_FCF'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
    plt.xticks(np.arange(min(MJD),max(MJD),tick_spacing_mjd),HSTticks_allsources,rotation='20')
    plt.suptitle(wave+' microns - Peak FCF, All Sources')

    if tick_spacing<300:
        plt.xlabel('HST')
    else:
        plt.xlabel('UT Date')

    plt.ylabel('Peak FCF')
    #plt.xticks(np.arange(min(MJD),max(MJD),tick_spacing_mjd),HSTticks,rotation='20')
    plt.legend(loc='upper left')
    if saveplots == False:
        #plt.axvline(x=57693.0,color='k',linestyle='solid',linewidth=2) # 57693 = 2016/11/01 - New Filters
        plt.show()
    else:
        plt.savefig('PeakFCF_vs_HST_ALLSOURCES_'+ending,dpi=300)
    plt.clf()


if individual_source_plots:
    for i in range(len(unique_sources)):
        if i != 6:
            plt.scatter(MJD[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM))],Arcsec_FCF[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM))],color=colours[i])
        else:
            plt.scatter(MJD[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM))],ARP220_Total_fix*Arcsec_FCF[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM))],color=colours[i])
    
    # Plot Colouring by weather:
    
        #plt.scatter(MJD[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM),np.logical_and(weather[np.where(SOURCE==unique_sources[i])]>0.12,weather[np.where(SOURCE==unique_sources[i])]<0.2)))],Arcsec_FCF[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM),np.logical_and(weather[np.where(SOURCE==unique_sources[i])]>=0.12,weather[np.where(SOURCE==unique_sources[i])]<0.2)))],color='y')
        #plt.scatter(MJD[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM),weather[np.where(SOURCE==unique_sources[i])]>=0.2))],Arcsec_FCF[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM),weather[np.where(SOURCE==unique_sources[i])]>=0.2))],color='r')
    
    
        if source_names[i] == 'CRL 2688':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['CRL2688']['mean_Arcsec_FCF'],color=colours[i],linestyle='solid',label='CRL 2688 '+epoch)
                plt.axhline(y=nominal_dict['CRL2688']['mean_Arcsec_FCF']+nominal_dict['CRL2688']['SD_Arcsec_FCF'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['CRL2688']['mean_Arcsec_FCF']-nominal_dict['CRL2688']['SD_Arcsec_FCF'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['CRL2688']['mean_Arcsec_FCF']-nominal_dict['CRL2688']['SD_Arcsec_FCF'],nominal_dict['CRL2688']['mean_Arcsec_FCF']+nominal_dict['CRL2688']['SD_Arcsec_FCF'],alpha=0.5,color=colours[i])
            print('\nCRL2688 Arcsec_FCF (mean,SD): ',nominal_dict['CRL2688']['mean_Arcsec_FCF'],nominal_dict['CRL2688']['SD_Arcsec_FCF'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_Arcsec_FCF']+comparison_dict[eachcomparison]['CRL2688']['SD_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_Arcsec_FCF']-comparison_dict[eachcomparison]['CRL2688']['SD_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['CRL2688']['mean_Arcsec_FCF']-comparison_dict[eachcomparison]['CRL2688']['SD_Arcsec_FCF'],comparison_dict[eachcomparison]['CRL2688']['mean_Arcsec_FCF']+comparison_dict[eachcomparison]['CRL2688']['SD_Arcsec_FCF'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 1)]),max(MJD[np.where(SOURCE == 1)]),tick_spacing_mjd),HSTticks_source1,rotation='20')
            plt.suptitle(wave+' microns - Arcsec FCF, '+source_names[i])
    
        if source_names[i] == 'CRL 618':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['CRL618']['mean_Arcsec_FCF'],color=colours[i],linestyle='solid',label='CRL 618 '+epoch)
                plt.axhline(y=nominal_dict['CRL618']['mean_Arcsec_FCF']+nominal_dict['CRL618']['SD_Arcsec_FCF'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['CRL618']['mean_Arcsec_FCF']-nominal_dict['CRL618']['SD_Arcsec_FCF'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['CRL618']['mean_Arcsec_FCF']-nominal_dict['CRL618']['SD_Arcsec_FCF'],nominal_dict['CRL618']['mean_Arcsec_FCF']+nominal_dict['CRL618']['SD_Arcsec_FCF'],alpha=0.5,color=colours[i])
            print('CRL618 Arcsec_FCF (mean,SD): ',nominal_dict['CRL618']['mean_Arcsec_FCF'],nominal_dict['CRL618']['SD_Arcsec_FCF'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_Arcsec_FCF']+comparison_dict[eachcomparison]['CRL618']['SD_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_Arcsec_FCF']-comparison_dict[eachcomparison]['CRL618']['SD_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['CRL618']['mean_Arcsec_FCF']-comparison_dict[eachcomparison]['CRL618']['SD_Arcsec_FCF'],comparison_dict[eachcomparison]['CRL618']['mean_Arcsec_FCF']+comparison_dict[eachcomparison]['CRL618']['SD_Arcsec_FCF'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 2)]),max(MJD[np.where(SOURCE == 2)]),tick_spacing_mjd),HSTticks_source2,rotation='20')
            plt.suptitle(wave+' microns - Arcsec FCF, '+source_names[i])
    
        if source_names[i] == 'Mars':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['MARS']['mean_Arcsec_FCF'],color=colours[i],linestyle='solid',label='Mars '+epoch)
                plt.axhline(y=nominal_dict['MARS']['mean_Arcsec_FCF']+nominal_dict['MARS']['SD_Arcsec_FCF'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['MARS']['mean_Arcsec_FCF']-nominal_dict['MARS']['SD_Arcsec_FCF'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['MARS']['mean_Arcsec_FCF']-nominal_dict['MARS']['SD_Arcsec_FCF'],nominal_dict['MARS']['mean_Arcsec_FCF']+nominal_dict['MARS']['SD_Arcsec_FCF'],alpha=0.5,color=colours[i])
            print('MARS Arcsec_FCF (mean,SD): ',nominal_dict['MARS']['mean_Arcsec_FCF'],nominal_dict['MARS']['SD_Arcsec_FCF'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_Arcsec_FCF']+comparison_dict[eachcomparison]['MARS']['SD_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_Arcsec_FCF']-comparison_dict[eachcomparison]['MARS']['SD_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['MARS']['mean_Arcsec_FCF']-comparison_dict[eachcomparison]['MARS']['SD_Arcsec_FCF'],comparison_dict[eachcomparison]['MARS']['mean_Arcsec_FCF']+comparison_dict[eachcomparison]['MARS']['SD_Arcsec_FCF'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 3)]),max(MJD[np.where(SOURCE == 3)]),tick_spacing_mjd),HSTticks_source3,rotation='20')
            plt.suptitle(wave+' microns - Arcsec FCF, '+source_names[i])
    
        if source_names[i] == 'Neptune':
            try:
                if plot_nominal_epoch_mean:
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_Arcsec_FCF'],color=colours[i],linestyle='solid',label='Neptune '+epoch)
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_Arcsec_FCF']+nominal_dict['NEPTUNE']['SD_Arcsec_FCF'],color=colours[i],linestyle='dotted')
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_Arcsec_FCF']-nominal_dict['NEPTUNE']['SD_Arcsec_FCF'],color=colours[i],linestyle='dotted')
                    plt.axhspan(nominal_dict['NEPTUNE']['mean_Arcsec_FCF']-nominal_dict['NEPTUNE']['SD_Arcsec_FCF'],nominal_dict['NEPTUNE']['mean_Arcsec_FCF']+nominal_dict['NEPTUNE']['SD_Arcsec_FCF'],alpha=0.5,color=colours[i])
                print('NEPTUNE Arcsec_FCF (mean,SD): ',nominal_dict['NEPTUNE']['mean_Arcsec_FCF'],nominal_dict['NEPTUNE']['SD_Arcsec_FCF'],'\n')
                for eachcomparison in comparison_dict.keys():
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_Arcsec_FCF']+comparison_dict[eachcomparison]['NEPTUNE']['SD_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_Arcsec_FCF']-comparison_dict[eachcomparison]['NEPTUNE']['SD_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                    plt.axhspan(comparison_dict[eachcomparison]['NEPTUNE']['mean_Arcsec_FCF']-comparison_dict[eachcomparison]['NEPTUNE']['SD_Arcsec_FCF'],comparison_dict[eachcomparison]['NEPTUNE']['mean_Arcsec_FCF']+comparison_dict[eachcomparison]['NEPTUNE']['SD_Arcsec_FCF'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
                plt.suptitle(wave+' microns - Arcsec FCF, '+source_names[i])
    
            except NameError:
                print('\nNeptune has no previous data with which to compare!\n')
                plt.xticks(np.arange(min(MJD[np.where(SOURCE == 4)]),max(MJD[np.where(SOURCE == 4)]),tick_spacing_mjd),HSTticks_source4,rotation='20')
                plt.suptitle('Arcsec FCF, '+source_names[i])
    
        if source_names[i] == 'Uranus':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['URANUS']['mean_Arcsec_FCF'],color=colours[i],linestyle='solid',label='Uranus '+epoch)
                plt.axhline(y=nominal_dict['URANUS']['mean_Arcsec_FCF']+nominal_dict['URANUS']['SD_Arcsec_FCF'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['URANUS']['mean_Arcsec_FCF']-nominal_dict['URANUS']['SD_Arcsec_FCF'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['URANUS']['mean_Arcsec_FCF']-nominal_dict['URANUS']['SD_Arcsec_FCF'],nominal_dict['URANUS']['mean_Arcsec_FCF']+nominal_dict['URANUS']['SD_Arcsec_FCF'],alpha=0.5,color=colours[i])
            print('URANUS Arcsec_FCF (mean,SD): ',nominal_dict['URANUS']['mean_Arcsec_FCF'],nominal_dict['URANUS']['SD_Arcsec_FCF'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_Arcsec_FCF']+comparison_dict[eachcomparison]['URANUS']['SD_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_Arcsec_FCF']-comparison_dict[eachcomparison]['URANUS']['SD_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['URANUS']['mean_Arcsec_FCF']-comparison_dict[eachcomparison]['URANUS']['SD_Arcsec_FCF'],comparison_dict[eachcomparison]['URANUS']['mean_Arcsec_FCF']+comparison_dict[eachcomparison]['URANUS']['SD_Arcsec_FCF'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            #plt.axhspan(source_5_mean_Arcsec_FCF-source_5_Arcsec_FCF_SD,source_5_mean_Arcsec_FCF+source_5_Arcsec_FCF_SD,alpha=0.5,color=colours[i])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 5)]),max(MJD[np.where(SOURCE == 5)]),tick_spacing_mjd),HSTticks_source5,rotation='20')
            plt.suptitle(wave+' microns - Arcsec FCF, '+source_names[i])

        if source_names[i] == 'Arp220':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['Arp220']['mean_Arcsec_FCF'],color=colours[i],linestyle='solid',label='Arp220 '+epoch)
                plt.axhline(y=nominal_dict['Arp220']['mean_Arcsec_FCF']+nominal_dict['Arp220']['SD_Arcsec_FCF'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['Arp220']['mean_Arcsec_FCF']-nominal_dict['Arp220']['SD_Arcsec_FCF'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['Arp220']['mean_Arcsec_FCF']-nominal_dict['Arp220']['SD_Arcsec_FCF'],nominal_dict['Arp220']['mean_Arcsec_FCF']+nominal_dict['Arp220']['SD_Arcsec_FCF'],alpha=0.5,color=colours[i])
            print('Arp220 Arcsec_FCF (mean,SD): ',nominal_dict['Arp220']['mean_Arcsec_FCF'],nominal_dict['Arp220']['SD_Arcsec_FCF'])
            for eachcomparison in comparison_dict.keys():
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_Arcsec_FCF']+comparison_dict[eachcomparison]['Arp220']['SD_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_Arcsec_FCF']-comparison_dict[eachcomparison]['Arp220']['SD_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
                plt.axhspan(comparison_dict[eachcomparison]['Arp220']['mean_Arcsec_FCF']-comparison_dict[eachcomparison]['Arp220']['SD_Arcsec_FCF'],comparison_dict[eachcomparison]['Arp220']['mean_Arcsec_FCF']+comparison_dict[eachcomparison]['Arp220']['SD_Arcsec_FCF'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
            #plt.axhspan(source_6_mean_Arcsec_FCF-source_6_Arcsec_FCF_SD,source_6_mean_Arcsec_FCF+source_6_Arcsec_FCF_SD,alpha=0.5,color=colours[i])
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 6)]),max(MJD[np.where(SOURCE == 6)]),tick_spacing_mjd),HSTticks_source6,rotation='20')
            plt.suptitle(wave+' microns - Arcsec FCF, '+source_names[i])

    
        if tick_spacing<300:
            plt.xlabel('HST')
        else:
            plt.xlabel('UT Date')
        plt.ylabel('Arcsec FCF')
        plt.legend(loc='upper left')
        if saveplots == False:
            #plt.axvline(x=57693.0,color='k',linestyle='solid',linewidth=2) # 57693 = New Filters (2016/11/01)
            plt.show()
        else:
            #plt.xlim(xmin=57693.0) # 57693 = New Filters (2016/11/01)
            plt.xlim(xmin=plt_date_min,xmax=plot_date_max) # 58484 = January 1st, 2019
            for eachvertline in vert_lines:
                plt.axvline(x=eachvertline,color='k',linestyle='solid',linewidth=2)
            #plt.ylim(ymin=1.5,ymax=30.0)
            #for eachwarmup in RxA_warmups:
            #    plt.axvline(x=eachwarmup,color='k',linestyle='dotted')
            #for eachcooldown in RxA_cooldowns:
            #    plt.axvline(x=eachcooldown,color='k',linestyle='dashdot')
            plt.savefig('ArcsecFCF_vs_HST_'+source_names[i].replace(' ','_')+ending,dpi=300)
        plt.clf()
else:
    plt.scatter(MJD[np.where(np.logical_and(Arcsec_FCF<ARCSEC_UPPER_LIM,Arcsec_FCF>ARCSEC_LOWER_LIM))],Arcsec_FCF[np.where(np.logical_and(Arcsec_FCF<ARCSEC_UPPER_LIM,Arcsec_FCF>ARCSEC_LOWER_LIM))],color='k')
    if plot_nominal_epoch_mean:
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_Arcsec_FCF'],color='k',linestyle='solid',label='All Sources '+epoch)
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_Arcsec_FCF']+nominal_dict['ALLSOURCES']['SD_Arcsec_FCF'],color='k',linestyle='dotted')
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_Arcsec_FCF']-nominal_dict['ALLSOURCES']['SD_Arcsec_FCF'],color='k',linestyle='dotted')
        plt.axhspan(nominal_dict['ALLSOURCES']['mean_Arcsec_FCF']-nominal_dict['ALLSOURCES']['SD_Arcsec_FCF'],nominal_dict['ALLSOURCES']['mean_Arcsec_FCF']+nominal_dict['ALLSOURCES']['SD_Arcsec_FCF'],alpha=0.5,color='k')
    print('ALLSOURCES Arcsec_FCF (mean,SD): ',nominal_dict['ALLSOURCES']['mean_Arcsec_FCF'],nominal_dict['ALLSOURCES']['SD_Arcsec_FCF'],'\n')
    for eachcomparison in comparison_dict.keys():
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_Arcsec_FCF']+comparison_dict[eachcomparison]['ALLSOURCES']['SD_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_Arcsec_FCF']-comparison_dict[eachcomparison]['ALLSOURCES']['SD_Arcsec_FCF'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')
        plt.axhspan(comparison_dict[eachcomparison]['ALLSOURCES']['mean_Arcsec_FCF']-comparison_dict[eachcomparison]['ALLSOURCES']['SD_Arcsec_FCF'],comparison_dict[eachcomparison]['ALLSOURCES']['mean_Arcsec_FCF']+comparison_dict[eachcomparison]['ALLSOURCES']['SD_Arcsec_FCF'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])
    plt.axhspan(allsources_mean_Arcsec_FCF-allsources_Arcsec_FCF_SD,allsources_mean_Arcsec_FCF+allsources_Arcsec_FCF_SD,alpha=0.5,color='k')
    plt.xticks(np.arange(min(MJD),max(MJD),tick_spacing_mjd),HSTticks_allsources,rotation='20')
    plt.suptitle(wave+' microns - Arcsec FCF, All Sources')

    if tick_spacing<300:
        plt.xlabel('HST')
    else:
        plt.xlabel('UT Date')
    plt.ylabel('Arcsec FCF')
    plt.legend(loc='upper left')
    if saveplots == False:
        #plt.axvline(x=58448.0,color='k',linestyle='dashed',linewidth=2)
        plt.show()
    else:
        plt.savefig('ArcsecFCF_vs_HST_ALLSOURCES_'+ending,dpi=300)
    plt.clf()


if individual_source_plots:
    for i in range(len(unique_sources)):
        plt.scatter(MJD[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM))],FWHM1[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM))],color=colours[i])                  
                                                                                                                                                                                                                                     
    # Plot Colouring by weather:                                                                                                                                                                                                     
                                                                                                                                                                                                                                     
        #plt.scatter(MJD[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM),np.logical_and(weather[np.where(SOURCE==unique_sources[i])]>0.12,weather[np.where(SOURCE==unique_sources[i])]<0.2)))],FWHM1[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM),np.logical_and(weather[np.where(SOURCE==unique_sources[i])]>=0.12,weather[np.where(SOURCE==unique_sources[i])]<0.2)))],color='y')                                                                                                                                                                                                                           
        #plt.scatter(MJD[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM),weather[np.where(SOURCE==unique_sources[i])]>=0.2))],FWHM1[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM),weather[np.where(SOURCE==unique_sources[i])]>=0.2))],color='r')                                                                                                                         
                                                                                                                                                                                                                                     
    #                                                                                                                                                                                                                                
                                                                                                                                                                                                                                     
        if source_names[i] == 'CRL 2688':                                                                                                                                                                                            
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['CRL2688']['mean_FWHM1'],color=colours[i],linestyle='solid',label='CRL 2688 '+epoch)                                                                                                     
                plt.axhline(y=nominal_dict['CRL2688']['mean_FWHM1']+nominal_dict['CRL2688']['SD_FWHM1'],color=colours[i],linestyle='dotted')                                                                                   
                plt.axhline(y=nominal_dict['CRL2688']['mean_FWHM1']-nominal_dict['CRL2688']['SD_FWHM1'],color=colours[i],linestyle='dotted')                                                                                   
                plt.axhspan(nominal_dict['CRL2688']['mean_FWHM1']-nominal_dict['CRL2688']['SD_FWHM1'],nominal_dict['CRL2688']['mean_FWHM1']+nominal_dict['CRL2688']['SD_FWHM1'],alpha=0.5,color=colours[i])          
            print('\nCRL2688 FWHM1 (mean,SD): ',nominal_dict['CRL2688']['mean_FWHM1'],nominal_dict['CRL2688']['SD_FWHM1'])                                                                                            
            for eachcomparison in comparison_dict.keys():                                                                                                                                                                            
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)                                 
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_FWHM1']+comparison_dict[eachcomparison]['CRL2688']['SD_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                              
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_FWHM1']-comparison_dict[eachcomparison]['CRL2688']['SD_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                              
                plt.axhspan(comparison_dict[eachcomparison]['CRL2688']['mean_FWHM1']-comparison_dict[eachcomparison]['CRL2688']['SD_FWHM1'],comparison_dict[eachcomparison]['CRL2688']['mean_FWHM1']+comparison_dict[eachcomparison]['CRL2688']['SD_FWHM1'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])                                                                                                               
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 1)]),max(MJD[np.where(SOURCE == 1)]),tick_spacing_mjd),HSTticks_source1,rotation='20')                                                                                   
            plt.suptitle(wave+' microns - FWHM1, '+source_names[i])                                                                                                                                                             
                                                                                                                                                                                                                                     
        if source_names[i] == 'CRL 618':                                                                                                                                                                                             
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['CRL618']['mean_FWHM1'],color=colours[i],linestyle='solid',label='CRL 618 '+epoch)                                                                                                       
                plt.axhline(y=nominal_dict['CRL618']['mean_FWHM1']+nominal_dict['CRL618']['SD_FWHM1'],color=colours[i],linestyle='dotted')                                                                                     
                plt.axhline(y=nominal_dict['CRL618']['mean_FWHM1']-nominal_dict['CRL618']['SD_FWHM1'],color=colours[i],linestyle='dotted')                                                                                     
                plt.axhspan(nominal_dict['CRL618']['mean_FWHM1']-nominal_dict['CRL618']['SD_FWHM1'],nominal_dict['CRL618']['mean_FWHM1']+nominal_dict['CRL618']['SD_FWHM1'],alpha=0.5,color=colours[i])              
            print('CRL618 FWHM1 (mean,SD): ',nominal_dict['CRL618']['mean_FWHM1'],nominal_dict['CRL618']['SD_FWHM1'])                                                                                                 
            for eachcomparison in comparison_dict.keys():                                                                                                                                                                            
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)                                  
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_FWHM1']+comparison_dict[eachcomparison]['CRL618']['SD_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_FWHM1']-comparison_dict[eachcomparison]['CRL618']['SD_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                
                plt.axhspan(comparison_dict[eachcomparison]['CRL618']['mean_FWHM1']-comparison_dict[eachcomparison]['CRL618']['SD_FWHM1'],comparison_dict[eachcomparison]['CRL618']['mean_FWHM1']+comparison_dict[eachcomparison]['CRL618']['SD_FWHM1'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])                                                                                                                   
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 2)]),max(MJD[np.where(SOURCE == 2)]),tick_spacing_mjd),HSTticks_source2,rotation='20')                                                                                   
            plt.suptitle(wave+' microns - FWHM1, '+source_names[i])                                                                                                                                                             
                                                                                                                                                                                                                                     
        if source_names[i] == 'Mars':                                                                                                                                                                                                
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['MARS']['mean_FWHM1'],color=colours[i],linestyle='solid',label='Mars '+epoch)                                                                                                            
                plt.axhline(y=nominal_dict['MARS']['mean_FWHM1']+nominal_dict['MARS']['SD_FWHM1'],color=colours[i],linestyle='dotted')                                                                                         
                plt.axhline(y=nominal_dict['MARS']['mean_FWHM1']-nominal_dict['MARS']['SD_FWHM1'],color=colours[i],linestyle='dotted')                                                                                         
                plt.axhspan(nominal_dict['MARS']['mean_FWHM1']-nominal_dict['MARS']['SD_FWHM1'],nominal_dict['MARS']['mean_FWHM1']+nominal_dict['MARS']['SD_FWHM1'],alpha=0.5,color=colours[i])                      
            print('MARS FWHM1 (mean,SD): ',nominal_dict['MARS']['mean_FWHM1'],nominal_dict['MARS']['SD_FWHM1'])                                                                                                       
            for eachcomparison in comparison_dict.keys():                                                                                                                                                                            
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)                                    
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_FWHM1']+comparison_dict[eachcomparison]['MARS']['SD_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                    
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_FWHM1']-comparison_dict[eachcomparison]['MARS']['SD_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                    
                plt.axhspan(comparison_dict[eachcomparison]['MARS']['mean_FWHM1']-comparison_dict[eachcomparison]['MARS']['SD_FWHM1'],comparison_dict[eachcomparison]['MARS']['mean_FWHM1']+comparison_dict[eachcomparison]['MARS']['SD_FWHM1'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])                                                                                                                           
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 3)]),max(MJD[np.where(SOURCE == 3)]),tick_spacing_mjd),HSTticks_source3,rotation='20')                                                                                   
            plt.suptitle(wave+' microns - FWHM1, '+source_names[i])                                                                                                                                                             
                                                                                                                                                                                                                                     
        if source_names[i] == 'Neptune':                                                                                                                                                                                             
            try:                                                                                                                                                                                                                     
                if plot_nominal_epoch_mean:
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_FWHM1'],color=colours[i],linestyle='solid',label='Neptune '+epoch)                                                                                                  
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_FWHM1']+nominal_dict['NEPTUNE']['SD_FWHM1'],color=colours[i],linestyle='dotted')                                                                               
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_FWHM1']-nominal_dict['NEPTUNE']['SD_FWHM1'],color=colours[i],linestyle='dotted')                                                                               
                    plt.axhspan(nominal_dict['NEPTUNE']['mean_FWHM1']-nominal_dict['NEPTUNE']['SD_FWHM1'],nominal_dict['NEPTUNE']['mean_FWHM1']+nominal_dict['NEPTUNE']['SD_FWHM1'],alpha=0.5,color=colours[i])      
                print('NEPTUNE FWHM1 (mean,SD): ',nominal_dict['NEPTUNE']['mean_FWHM1'],nominal_dict['NEPTUNE']['SD_FWHM1'],'\n')                                                                                     
                for eachcomparison in comparison_dict.keys():                                                                                                                                                                        
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)                             
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_FWHM1']+comparison_dict[eachcomparison]['NEPTUNE']['SD_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                          
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_FWHM1']-comparison_dict[eachcomparison]['NEPTUNE']['SD_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                          
                    plt.axhspan(comparison_dict[eachcomparison]['NEPTUNE']['mean_FWHM1']-comparison_dict[eachcomparison]['NEPTUNE']['SD_FWHM1'],comparison_dict[eachcomparison]['NEPTUNE']['mean_FWHM1']+comparison_dict[eachcomparison]['NEPTUNE']['SD_FWHM1'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])                                                                                                           
                plt.suptitle(wave+' microns - FWHM1, '+source_names[i])                                                                                                                                                         
                                                                                                                                                                                                                                     
            except NameError:                                                                                                                                                                                                        
                print('\nNeptune has no previous data with which to compare!\n')                                                                                                                                                     
                plt.xticks(np.arange(min(MJD[np.where(SOURCE == 4)]),max(MJD[np.where(SOURCE == 4)]),tick_spacing_mjd),HSTticks_source4,rotation='20')                                                                               
                plt.suptitle('FWHM1, '+source_names[i])                                                                                                                                                                         
                                                                                                                                                                                                                                     
        if source_names[i] == 'Uranus':                                                                                                                                                                                              
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['URANUS']['mean_FWHM1'],color=colours[i],linestyle='solid',label='Uranus '+epoch)                                                                                                        
                plt.axhline(y=nominal_dict['URANUS']['mean_FWHM1']+nominal_dict['URANUS']['SD_FWHM1'],color=colours[i],linestyle='dotted')                                                                                     
                plt.axhline(y=nominal_dict['URANUS']['mean_FWHM1']-nominal_dict['URANUS']['SD_FWHM1'],color=colours[i],linestyle='dotted')                                                                                     
                plt.axhspan(nominal_dict['URANUS']['mean_FWHM1']-nominal_dict['URANUS']['SD_FWHM1'],nominal_dict['URANUS']['mean_FWHM1']+nominal_dict['URANUS']['SD_FWHM1'],alpha=0.5,color=colours[i])              
            print('URANUS FWHM1 (mean,SD): ',nominal_dict['URANUS']['mean_FWHM1'],nominal_dict['URANUS']['SD_FWHM1'])                                                                                                 
            for eachcomparison in comparison_dict.keys():                                                                                                                                                                            
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)                                  
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_FWHM1']+comparison_dict[eachcomparison]['URANUS']['SD_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_FWHM1']-comparison_dict[eachcomparison]['URANUS']['SD_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                
                plt.axhspan(comparison_dict[eachcomparison]['URANUS']['mean_FWHM1']-comparison_dict[eachcomparison]['URANUS']['SD_FWHM1'],comparison_dict[eachcomparison]['URANUS']['mean_FWHM1']+comparison_dict[eachcomparison]['URANUS']['SD_FWHM1'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])                                                                                                                   
            plt.axhspan(source_5_mean_FWHM1-source_5_FWHM1_SD,source_5_mean_FWHM1+source_5_FWHM1_SD,alpha=0.5,color=colours[i])                                                                                  
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 5)]),max(MJD[np.where(SOURCE == 5)]),tick_spacing_mjd),HSTticks_source5,rotation='20')                                                                                   
            plt.suptitle(wave+' microns - FWHM1, '+source_names[i])                                                                                                                                                             

        if source_names[i] == 'Arp220':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['Arp220']['mean_FWHM1'],color=colours[i],linestyle='solid',label='Arp220 '+epoch)
                plt.axhline(y=nominal_dict['Arp220']['mean_FWHM1']+nominal_dict['Arp220']['SD_FWHM1'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['Arp220']['mean_FWHM1']-nominal_dict['Arp220']['SD_FWHM1'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['Arp220']['mean_FWHM1']-nominal_dict['Arp220']['SD_FWHM1'],nominal_dict['Arp220']['mean_FWHM1']+nominal_dict['Arp220']['SD_FWHM1'],alpha=0.5,color=colours[i])
            print('Arp220 FWHM1 (mean,SD): ',nominal_dict['Arp220']['mean_FWHM1'],nominal_dict['Arp220']['SD_FWHM1'])                                                                                   
            for eachcomparison in comparison_dict.keys():                                                                                                                                                              
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)                    
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_FWHM1']+comparison_dict[eachcomparison]['Arp220']['SD_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_FWHM1']-comparison_dict[eachcomparison]['Arp220']['SD_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                
                plt.axhspan(comparison_dict[eachcomparison]['Arp220']['mean_FWHM1']-comparison_dict[eachcomparison]['Arp220']['SD_FWHM1'],comparison_dict[eachcomparison]['Arp220']['mean_FWHM1']+comparison_dict[eachcomparison]['Arp220']['SD_FWHM1'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])                                                                                                                   
            plt.axhspan(source_6_mean_FWHM1-source_6_FWHM1_SD,source_6_mean_FWHM1+source_6_FWHM1_SD,alpha=0.5,color=colours[i])                                                                                  
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 6)]),max(MJD[np.where(SOURCE == 6)]),tick_spacing_mjd),HSTticks_source6,rotation='20')                                                                                   
            plt.suptitle(wave+' microns - FWHM1, '+source_names[i])                                                                                                                                                             

    
        if tick_spacing<300:
            plt.xlabel('HST')
        else:                
            plt.xlabel('UT Date')
        plt.ylabel('FWHM1') 
        plt.legend(loc='upper left')
        if saveplots == False:      
            #plt.axvline(x=58448.0,color='k',linestyle='dashed',linewidth=2)
            plt.show()                                                     
        else:                                                              
            plt.savefig('FWHM1_vs_HST_'+source_names[i].replace(' ','_')+ending,dpi=300)
        plt.clf()                                                                           
else:                                                                                       
    plt.scatter(MJD[np.where(np.logical_and(Arcsec_FCF<ARCSEC_UPPER_LIM,Arcsec_FCF>ARCSEC_LOWER_LIM))],FWHM1[np.where(np.logical_and(Arcsec_FCF<ARCSEC_UPPER_LIM,Arcsec_FCF>ARCSEC_LOWER_LIM))],color='k')
    if plot_nominal_epoch_mean:
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_FWHM1'],color='k',linestyle='solid',label='All Sources '+epoch)                                                                                        
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_FWHM1']+nominal_dict['ALLSOURCES']['SD_FWHM1'],color='k',linestyle='dotted')                                                                      
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_FWHM1']-nominal_dict['ALLSOURCES']['SD_FWHM1'],color='k',linestyle='dotted')                                                                      
        plt.axhspan(nominal_dict['ALLSOURCES']['mean_FWHM1']-nominal_dict['ALLSOURCES']['SD_FWHM1'],nominal_dict['ALLSOURCES']['mean_FWHM1']+nominal_dict['ALLSOURCES']['SD_FWHM1'],alpha=0.5,color='k')
    print('ALLSOURCES FWHM1 (mean,SD): ',nominal_dict['ALLSOURCES']['mean_FWHM1'],nominal_dict['ALLSOURCES']['SD_FWHM1'],'\n')                                                                           
    for eachcomparison in comparison_dict.keys():                                                                                                                                                                       
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)                         
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_FWHM1']+comparison_dict[eachcomparison]['ALLSOURCES']['SD_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_FWHM1']-comparison_dict[eachcomparison]['ALLSOURCES']['SD_FWHM1'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                
        plt.axhspan(comparison_dict[eachcomparison]['ALLSOURCES']['mean_FWHM1']-comparison_dict[eachcomparison]['ALLSOURCES']['SD_FWHM1'],comparison_dict[eachcomparison]['ALLSOURCES']['mean_FWHM1']+comparison_dict[eachcomparison]['ALLSOURCES']['SD_FWHM1'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])                                                                                                           
    plt.axhspan(allsources_mean_FWHM1-allsources_FWHM1_SD,allsources_mean_FWHM1+allsources_FWHM1_SD,alpha=0.5,color='k')      
    plt.xticks(np.arange(min(MJD),max(MJD),tick_spacing_mjd),HSTticks_allsources,rotation='20')                                                                                                                                      
    plt.suptitle(wave+' microns - FWHM1, All Sources')                                                                                                                                                                          

    if tick_spacing<300:
        plt.xlabel('HST')
    else:                
        plt.xlabel('UT Date')
    plt.ylabel('FWHM1') 
    plt.legend(loc='upper left')
    if saveplots == False:      
        #plt.axvline(x=58448.0,color='k',linestyle='dashed',linewidth=2)
        plt.show()                                                     
    else:                                                              
        plt.savefig('FWHM1_vs_HST_ALLSOURCES_'+ending,dpi=300)     
    plt.clf()                                                          



###### HERE #######


if individual_source_plots:
    for i in range(len(unique_sources)):
        plt.scatter(MJD[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM))],FWHM2[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM))],color=colours[i])                       
                                                                                                                                                                                                                                     
    # Plot Colouring by weather:                                                                                                                                                                                                     
                                                                                                                                                                                                                                     
        #plt.scatter(MJD[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM),np.logical_and(weather[np.where(SOURCE==unique_sources[i])]>0.12,weather[np.where(SOURCE==unique_sources[i])]<0.2)))],FWHM2[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM),np.logical_and(weather[np.where(SOURCE==unique_sources[i])]>=0.12,weather[np.where(SOURCE==unique_sources[i])]<0.2)))],color='y')                                                                                                                                                                                                                                
        #plt.scatter(MJD[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM),weather[np.where(SOURCE==unique_sources[i])]>=0.2))],FWHM2[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM),weather[np.where(SOURCE==unique_sources[i])]>=0.2))],color='r')                                                                                                                              
                                                                                                                                                                                                                                     
    #                                                                                                                                                                                                                                
                                                                                                                                                                                                                                     
        if source_names[i] == 'CRL 2688':                                                                                                                                                                                            
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['CRL2688']['mean_FWHM2'],color=colours[i],linestyle='solid',label='CRL 2688 '+epoch)                                                                                                          
                plt.axhline(y=nominal_dict['CRL2688']['mean_FWHM2']+nominal_dict['CRL2688']['SD_FWHM2'],color=colours[i],linestyle='dotted')                                                                                             
                plt.axhline(y=nominal_dict['CRL2688']['mean_FWHM2']-nominal_dict['CRL2688']['SD_FWHM2'],color=colours[i],linestyle='dotted')                                                                                             
                plt.axhspan(nominal_dict['CRL2688']['mean_FWHM2']-nominal_dict['CRL2688']['SD_FWHM2'],nominal_dict['CRL2688']['mean_FWHM2']+nominal_dict['CRL2688']['SD_FWHM2'],alpha=0.5,color=colours[i])                              
            print('\nCRL2688 FWHM2 (mean,SD): ',nominal_dict['CRL2688']['mean_FWHM2'],nominal_dict['CRL2688']['SD_FWHM2'])                                                                                                      
            for eachcomparison in comparison_dict.keys():                                                                                                                                                                            
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)                                      
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_FWHM2']+comparison_dict[eachcomparison]['CRL2688']['SD_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                        
                plt.axhline(y=comparison_dict[eachcomparison]['CRL2688']['mean_FWHM2']-comparison_dict[eachcomparison]['CRL2688']['SD_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                        
                plt.axhspan(comparison_dict[eachcomparison]['CRL2688']['mean_FWHM2']-comparison_dict[eachcomparison]['CRL2688']['SD_FWHM2'],comparison_dict[eachcomparison]['CRL2688']['mean_FWHM2']+comparison_dict[eachcomparison]['CRL2688']['SD_FWHM2'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])                                                                                                                                   
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 1)]),max(MJD[np.where(SOURCE == 1)]),tick_spacing_mjd),HSTticks_source1,rotation='20')                                                                                   
            plt.suptitle(wave+' microns - FWHM2, '+source_names[i])                                                                                                                                                                  
                                                                                                                                                                                                                                     
        if source_names[i] == 'CRL 618':                                                                                                                                                                                             
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['CRL618']['mean_FWHM2'],color=colours[i],linestyle='solid',label='CRL 618 '+epoch)                                                                                                            
                plt.axhline(y=nominal_dict['CRL618']['mean_FWHM2']+nominal_dict['CRL618']['SD_FWHM2'],color=colours[i],linestyle='dotted')                                                                                               
                plt.axhline(y=nominal_dict['CRL618']['mean_FWHM2']-nominal_dict['CRL618']['SD_FWHM2'],color=colours[i],linestyle='dotted')                                                                                               
                plt.axhspan(nominal_dict['CRL618']['mean_FWHM2']-nominal_dict['CRL618']['SD_FWHM2'],nominal_dict['CRL618']['mean_FWHM2']+nominal_dict['CRL618']['SD_FWHM2'],alpha=0.5,color=colours[i])                                  
            print('CRL618 FWHM2 (mean,SD): ',nominal_dict['CRL618']['mean_FWHM2'],nominal_dict['CRL618']['SD_FWHM2'])                                                                                                           
            for eachcomparison in comparison_dict.keys():                                                                                                                                                                            
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)                                       
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_FWHM2']+comparison_dict[eachcomparison]['CRL618']['SD_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                          
                plt.axhline(y=comparison_dict[eachcomparison]['CRL618']['mean_FWHM2']-comparison_dict[eachcomparison]['CRL618']['SD_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                          
                plt.axhspan(comparison_dict[eachcomparison]['CRL618']['mean_FWHM2']-comparison_dict[eachcomparison]['CRL618']['SD_FWHM2'],comparison_dict[eachcomparison]['CRL618']['mean_FWHM2']+comparison_dict[eachcomparison]['CRL618']['SD_FWHM2'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])                                                                                                                                       
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 2)]),max(MJD[np.where(SOURCE == 2)]),tick_spacing_mjd),HSTticks_source2,rotation='20')                                                                                   
            plt.suptitle(wave+' microns - FWHM2, '+source_names[i])                                                                                                                                                                  
                                                                                                                                                                                                                                     
        if source_names[i] == 'Mars':                                                                                                                                                                                                
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['MARS']['mean_FWHM2'],color=colours[i],linestyle='solid',label='Mars '+epoch)                                                                                                                 
                plt.axhline(y=nominal_dict['MARS']['mean_FWHM2']+nominal_dict['MARS']['SD_FWHM2'],color=colours[i],linestyle='dotted')                                                                                                   
                plt.axhline(y=nominal_dict['MARS']['mean_FWHM2']-nominal_dict['MARS']['SD_FWHM2'],color=colours[i],linestyle='dotted')                                                                                                   
                plt.axhspan(nominal_dict['MARS']['mean_FWHM2']-nominal_dict['MARS']['SD_FWHM2'],nominal_dict['MARS']['mean_FWHM2']+nominal_dict['MARS']['SD_FWHM2'],alpha=0.5,color=colours[i])                                          
            print('MARS FWHM2 (mean,SD): ',nominal_dict['MARS']['mean_FWHM2'],nominal_dict['MARS']['SD_FWHM2'])                                                                                                                 
            for eachcomparison in comparison_dict.keys():                                                                                                                                                                            
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)                                         
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_FWHM2']+comparison_dict[eachcomparison]['MARS']['SD_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                              
                plt.axhline(y=comparison_dict[eachcomparison]['MARS']['mean_FWHM2']-comparison_dict[eachcomparison]['MARS']['SD_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                              
                plt.axhspan(comparison_dict[eachcomparison]['MARS']['mean_FWHM2']-comparison_dict[eachcomparison]['MARS']['SD_FWHM2'],comparison_dict[eachcomparison]['MARS']['mean_FWHM2']+comparison_dict[eachcomparison]['MARS']['SD_FWHM2'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])                                                                                                                                               
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 3)]),max(MJD[np.where(SOURCE == 3)]),tick_spacing_mjd),HSTticks_source3,rotation='20')                                                                                   
            plt.suptitle(wave+' microns - FWHM2, '+source_names[i])                                                                                                                                                                  
                                                                                                                                                                                                                                     
        if source_names[i] == 'Neptune':                                                                                                                                                                                             
            try:                                                                                                                                                                                                                     
                if plot_nominal_epoch_mean:
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_FWHM2'],color=colours[i],linestyle='solid',label='Neptune '+epoch)                                                                                                       
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_FWHM2']+nominal_dict['NEPTUNE']['SD_FWHM2'],color=colours[i],linestyle='dotted')                                                                                         
                    plt.axhline(y=nominal_dict['NEPTUNE']['mean_FWHM2']-nominal_dict['NEPTUNE']['SD_FWHM2'],color=colours[i],linestyle='dotted')                                                                                         
                    plt.axhspan(nominal_dict['NEPTUNE']['mean_FWHM2']-nominal_dict['NEPTUNE']['SD_FWHM2'],nominal_dict['NEPTUNE']['mean_FWHM2']+nominal_dict['NEPTUNE']['SD_FWHM2'],alpha=0.5,color=colours[i])                          
                print('NEPTUNE FWHM2 (mean,SD): ',nominal_dict['NEPTUNE']['mean_FWHM2'],nominal_dict['NEPTUNE']['SD_FWHM2'],'\n')                                                                                               
                for eachcomparison in comparison_dict.keys():                                                                                                                                                                        
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)                                  
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_FWHM2']+comparison_dict[eachcomparison]['NEPTUNE']['SD_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                    
                    plt.axhline(y=comparison_dict[eachcomparison]['NEPTUNE']['mean_FWHM2']-comparison_dict[eachcomparison]['NEPTUNE']['SD_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                    
                    plt.axhspan(comparison_dict[eachcomparison]['NEPTUNE']['mean_FWHM2']-comparison_dict[eachcomparison]['NEPTUNE']['SD_FWHM2'],comparison_dict[eachcomparison]['NEPTUNE']['mean_FWHM2']+comparison_dict[eachcomparison]['NEPTUNE']['SD_FWHM2'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])                                                                                                                               
                plt.suptitle(wave+' microns - FWHM2, '+source_names[i])                                                                                                                                                              
                                                                                                                                                                                                                                     
            except NameError:                                                                                                                                                                                                        
                print('\nNeptune has no previous data with which to compare!\n')                                                                                                                                                     
                plt.xticks(np.arange(min(MJD[np.where(SOURCE == 4)]),max(MJD[np.where(SOURCE == 4)]),tick_spacing_mjd),HSTticks_source4,rotation='20')                                                                               
                plt.suptitle('FWHM2, '+source_names[i])                                                                                                                                                                              
                                                                                                                                                                                                                                     
        if source_names[i] == 'Uranus':                                                                                                                                                                                              
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['URANUS']['mean_FWHM2'],color=colours[i],linestyle='solid',label='Uranus '+epoch)                                                                                                             
                plt.axhline(y=nominal_dict['URANUS']['mean_FWHM2']+nominal_dict['URANUS']['SD_FWHM2'],color=colours[i],linestyle='dotted')                                                                                               
                plt.axhline(y=nominal_dict['URANUS']['mean_FWHM2']-nominal_dict['URANUS']['SD_FWHM2'],color=colours[i],linestyle='dotted')                                                                                               
                plt.axhspan(nominal_dict['URANUS']['mean_FWHM2']-nominal_dict['URANUS']['SD_FWHM2'],nominal_dict['URANUS']['mean_FWHM2']+nominal_dict['URANUS']['SD_FWHM2'],alpha=0.5,color=colours[i])                                  
            print('URANUS FWHM2 (mean,SD): ',nominal_dict['URANUS']['mean_FWHM2'],nominal_dict['URANUS']['SD_FWHM2'])                                                                                                           
            for eachcomparison in comparison_dict.keys():                                                                                                                                                                            
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)                                       
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_FWHM2']+comparison_dict[eachcomparison]['URANUS']['SD_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                          
                plt.axhline(y=comparison_dict[eachcomparison]['URANUS']['mean_FWHM2']-comparison_dict[eachcomparison]['URANUS']['SD_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                          
                plt.axhspan(comparison_dict[eachcomparison]['URANUS']['mean_FWHM2']-comparison_dict[eachcomparison]['URANUS']['SD_FWHM2'],comparison_dict[eachcomparison]['URANUS']['mean_FWHM2']+comparison_dict[eachcomparison]['URANUS']['SD_FWHM2'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])                                                                                                                                       
            plt.axhspan(source_5_mean_FWHM2-source_5_FWHM2_SD,source_5_mean_FWHM2+source_5_FWHM2_SD,alpha=0.5,color=colours[i])                                                                                            
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 5)]),max(MJD[np.where(SOURCE == 5)]),tick_spacing_mjd),HSTticks_source5,rotation='20')                                                                                   
            plt.suptitle(wave+' microns - FWHM2, '+source_names[i])                                                                                                                                                                  

        if source_names[i] == 'Arp220':
            if plot_nominal_epoch_mean:
                plt.axhline(y=nominal_dict['Arp220']['mean_FWHM2'],color=colours[i],linestyle='solid',label='Arp220 '+epoch)
                plt.axhline(y=nominal_dict['Arp220']['mean_FWHM2']+nominal_dict['Arp220']['SD_FWHM2'],color=colours[i],linestyle='dotted')
                plt.axhline(y=nominal_dict['Arp220']['mean_FWHM2']-nominal_dict['Arp220']['SD_FWHM2'],color=colours[i],linestyle='dotted')
                plt.axhspan(nominal_dict['Arp220']['mean_FWHM2']-nominal_dict['Arp220']['SD_FWHM2'],nominal_dict['Arp220']['mean_FWHM2']+nominal_dict['Arp220']['SD_FWHM2'],alpha=0.5,color=colours[i])
            print('Arp220 FWHM2 (mean,SD): ',nominal_dict['Arp220']['mean_FWHM2'],nominal_dict['Arp220']['SD_FWHM2'])                                                                                   
            for eachcomparison in comparison_dict.keys():                                                                                                                                                              
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)                         
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_FWHM2']+comparison_dict[eachcomparison]['Arp220']['SD_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                          
                plt.axhline(y=comparison_dict[eachcomparison]['Arp220']['mean_FWHM2']-comparison_dict[eachcomparison]['Arp220']['SD_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                          
                plt.axhspan(comparison_dict[eachcomparison]['Arp220']['mean_FWHM2']-comparison_dict[eachcomparison]['Arp220']['SD_FWHM2'],comparison_dict[eachcomparison]['Arp220']['mean_FWHM2']+comparison_dict[eachcomparison]['Arp220']['SD_FWHM2'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])                                                                                                                                       
            plt.axhspan(source_6_mean_FWHM2-source_6_FWHM2_SD,source_6_mean_FWHM2+source_6_FWHM2_SD,alpha=0.5,color=colours[i])                                                                                            
            plt.xticks(np.arange(min(MJD[np.where(SOURCE == 6)]),max(MJD[np.where(SOURCE == 6)]),tick_spacing_mjd),HSTticks_source6,rotation='20')                                                                                   
            plt.suptitle(wave+' microns - FWHM2, '+source_names[i])                                                                                                                                                                  

    
        if tick_spacing<300:
            plt.xlabel('HST')
        else:                
            plt.xlabel('UT Date')
        plt.ylabel('FWHM2')      
        plt.legend(loc='upper left')
        if saveplots == False:      
            #plt.axvline(x=58448.0,color='k',linestyle='dashed',linewidth=2)
            plt.show()                                                     
        else:                                                              
            plt.savefig('FWHM2_vs_HST_'+source_names[i].replace(' ','_')+ending,dpi=300)
        plt.clf()                                                                           
else:                                                                                       
    plt.scatter(MJD[np.where(np.logical_and(Arcsec_FCF<ARCSEC_UPPER_LIM,Arcsec_FCF>ARCSEC_LOWER_LIM))],FWHM2[np.where(np.logical_and(Arcsec_FCF<ARCSEC_UPPER_LIM,Arcsec_FCF>ARCSEC_LOWER_LIM))],color='k')
    if plot_nominal_epoch_mean:
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_FWHM2'],color='k',linestyle='solid',label='All Sources '+epoch)                                                                                             
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_FWHM2']+nominal_dict['ALLSOURCES']['SD_FWHM2'],color='k',linestyle='dotted')                                                                                
        plt.axhline(y=nominal_dict['ALLSOURCES']['mean_FWHM2']-nominal_dict['ALLSOURCES']['SD_FWHM2'],color='k',linestyle='dotted')                                                                                
        plt.axhspan(nominal_dict['ALLSOURCES']['mean_FWHM2']-nominal_dict['ALLSOURCES']['SD_FWHM2'],nominal_dict['ALLSOURCES']['mean_FWHM2']+nominal_dict['ALLSOURCES']['SD_FWHM2'],alpha=0.5,color='k')           
    print('ALLSOURCES FWHM2 (mean,SD): ',nominal_dict['ALLSOURCES']['mean_FWHM2'],nominal_dict['ALLSOURCES']['SD_FWHM2'],'\n')                                                                            
    for eachcomparison in comparison_dict.keys():                                                                                                                                                                       
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='solid',label=eachcomparison)                              
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_FWHM2']+comparison_dict[eachcomparison]['ALLSOURCES']['SD_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                          
        plt.axhline(y=comparison_dict[eachcomparison]['ALLSOURCES']['mean_FWHM2']-comparison_dict[eachcomparison]['ALLSOURCES']['SD_FWHM2'],color=comparison_colours[np.where(epochs==eachcomparison)[0][0]],linestyle='dotted')                                                                                                                                                                                                                                          
        plt.axhspan(comparison_dict[eachcomparison]['ALLSOURCES']['mean_FWHM2']-comparison_dict[eachcomparison]['ALLSOURCES']['SD_FWHM2'],comparison_dict[eachcomparison]['ALLSOURCES']['mean_FWHM2']+comparison_dict[eachcomparison]['ALLSOURCES']['SD_FWHM2'],alpha=0.5,color=comparison_colours[np.where(epochs==eachcomparison)[0][0]])                                                                                                                               
    plt.axhspan(allsources_mean_FWHM2-allsources_FWHM2_SD,allsources_mean_FWHM2+allsources_FWHM2_SD,alpha=0.5,color='k')                                                                                                   
    plt.xticks(np.arange(min(MJD),max(MJD),tick_spacing_mjd),HSTticks_allsources,rotation='20')                                                                                                                                      
    plt.suptitle(wave+' microns - FWHM2, All Sources')                                                                                                                                                                               

    if tick_spacing<300:
        plt.xlabel('HST')
    else:                
        plt.xlabel('UT Date')
    plt.ylabel('FWHM2')      
    plt.legend(loc='upper left')
    if saveplots == False:      
        #plt.axvline(x=58448.0,color='k',linestyle='dashed',linewidth=2)
        plt.show()                                                     
    else:                                                              
        plt.savefig('FWHM2_vs_HST_ALLSOURCES_'+ending,dpi=300)         
    plt.clf()                                                          


if individual_source_plots:
    for i in range(len(unique_sources)):
        plt.scatter(TauTimesAM[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM))],Arcsec_FCF[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<ARCSEC_UPPER_LIM,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>ARCSEC_LOWER_LIM))],color=colours[i])
    #
    
    # Colour by Tau225:
    
        #plt.scatter(TauTimesAM[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<10,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>0),weather[np.where(SOURCE==unique_sources[i])]>=0.2))],Arcsec_FCF[np.where(SOURCE==unique_sources[i])][np.where(np.logical_and(np.logical_and(Arcsec_FCF[np.where(SOURCE==unique_sources[i])]<10,Arcsec_FCF[np.where(SOURCE==unique_sources[i])]>0),weather[np.where(SOURCE==unique_sources[i])]>=0.2))],color='r')
    
    #
    
        if source_names[i] == 'CRL 2688':
            if plot_nominal_epoch_mean:
                plt.axhline(y=source_1_mean_Arcsec_FCF,color=colours[i],linestyle='solid',label='CRL 2688 Nominal')
                plt.axhline(y=source_1_mean_Arcsec_FCF+source_1_Arcsec_FCF_SD,color=colours[i],linestyle='dotted')
                plt.axhline(y=source_1_mean_Arcsec_FCF-source_1_Arcsec_FCF_SD,color=colours[i],linestyle='dotted')
                plt.axhspan(source_1_mean_Arcsec_FCF-source_1_Arcsec_FCF_SD,source_1_mean_Arcsec_FCF+source_1_Arcsec_FCF_SD,alpha=0.5,color=colours[i])
            plt.suptitle('Arcsec FCF, '+source_names[i])
    
        if source_names[i] == 'CRL 618':
            if plot_nominal_epoch_mean:
                plt.axhline(y=source_2_mean_Arcsec_FCF,color=colours[i],linestyle='solid',label='CRL 618 Nominal')
                plt.axhline(y=source_2_mean_Arcsec_FCF+source_2_Arcsec_FCF_SD,color=colours[i],linestyle='dotted')
                plt.axhline(y=source_2_mean_Arcsec_FCF-source_2_Arcsec_FCF_SD,color=colours[i],linestyle='dotted')
                plt.axhspan(source_2_mean_Arcsec_FCF-source_2_Arcsec_FCF_SD,source_2_mean_Arcsec_FCF+source_2_Arcsec_FCF_SD,alpha=0.5,color=colours[i])
            plt.suptitle('Arcsec FCF, '+source_names[i])
    
        if source_names[i] == 'Mars':
            if plot_nominal_epoch_mean:
                plt.axhline(y=source_3_mean_Arcsec_FCF,color=colours[i],linestyle='solid',label='Mars Nominal')
                plt.axhline(y=source_3_mean_Arcsec_FCF+source_3_Arcsec_FCF_SD,color=colours[i],linestyle='dotted')
                plt.axhline(y=source_3_mean_Arcsec_FCF-source_3_Arcsec_FCF_SD,color=colours[i],linestyle='dotted')
                plt.axhspan(source_3_mean_Arcsec_FCF-source_3_Arcsec_FCF_SD,source_3_mean_Arcsec_FCF+source_3_Arcsec_FCF_SD,alpha=0.5,color=colours[i])
            plt.suptitle('Arcsec FCF, '+source_names[i])
    
        if source_names[i] == 'Neptune':
            try:
                if plot_nominal_epoch_mean:
                    plt.axhline(y=source_4_mean_Arcsec_FCF,color=colours[i],linestyle='solid',label='Neptune Nominal')
                    plt.axhline(y=source_4_mean_Arcsec_FCF+source_4_Arcsec_FCF_SD,color=colours[i],linestyle='dotted')
                    plt.axhline(y=source_4_mean_Arcsec_FCF-source_4_Arcsec_FCF_SD,color=colours[i],linestyle='dotted')
                plt.xticks(np.arange(min(MJD[np.where(SOURCE == 4)]),max(MJD[np.where(SOURCE == 4)]),tick_spacing_mjd),HSTticks_source4,rotation='20')
                plt.suptitle('Arcsec FCF, '+source_names[i])
            except NameError:
                print('\nNeptune has no previous data with which to compare!\n')
                plt.suptitle('Arcsec FCF, '+source_names[i])
    
        if source_names[i] == 'Uranus':
            if plot_nominal_epoch_mean:
                plt.axhline(y=source_5_mean_Arcsec_FCF,color=colours[i],linestyle='solid',label='Uranus Nominal')
                plt.axhline(y=source_5_mean_Arcsec_FCF+source_5_Arcsec_FCF_SD,color=colours[i],linestyle='dotted')
                plt.axhline(y=source_5_mean_Arcsec_FCF-source_5_Arcsec_FCF_SD,color=colours[i],linestyle='dotted')
                plt.axhspan(source_5_mean_Arcsec_FCF-source_5_Arcsec_FCF_SD,source_5_mean_Arcsec_FCF+source_5_Arcsec_FCF_SD,alpha=0.5,color=colours[i])
            plt.suptitle('Arcsec FCF, '+source_names[i])

        if source_names[i] == 'Arp220':
            if plot_nominal_epoch_mean:
                plt.axhline(y=source_6_mean_Arcsec_FCF,color=colours[i],linestyle='solid',label='Arp220 Nominal')
                plt.axhline(y=source_6_mean_Arcsec_FCF+source_6_Arcsec_FCF_SD,color=colours[i],linestyle='dotted')
                plt.axhline(y=source_6_mean_Arcsec_FCF-source_6_Arcsec_FCF_SD,color=colours[i],linestyle='dotted')
                plt.axhspan(source_6_mean_Arcsec_FCF-source_6_Arcsec_FCF_SD,source_6_mean_Arcsec_FCF+source_6_Arcsec_FCF_SD,alpha=0.5,color=colours[i])
            plt.suptitle('Arcsec FCF, '+source_names[i])
    
        if tick_spacing<300:
            plt.xlabel('TauTimesAM')
        else:
            plt.xlabel('TauTimesAM')
        plt.ylabel('Arcsec FCF')
        plt.legend(loc='lower left')
        if saveplots == False:
            plt.show()
        else:
            plt.savefig('ArcsecFCF_vs_TauTimesAM_'+source_names[i].replace(' ','_')+ending,dpi=300)
        plt.clf()
else:
    plt.scatter(TauTimesAM[np.where(np.logical_and(Arcsec_FCF<ARCSEC_UPPER_LIM,Arcsec_FCF>ARCSEC_LOWER_LIM))],Arcsec_FCF[np.where(np.logical_and(Arcsec_FCF<ARCSEC_UPPER_LIM,Arcsec_FCF>ARCSEC_LOWER_LIM))],color='k')
    if plot_nominal_epoch_mean:
        plt.axhline(y=allsources_mean_Arcsec_FCF,color='k',linestyle='solid',label='All Sources Nominal')
        plt.axhline(y=allsources_mean_Arcsec_FCF+allsources_Arcsec_FCF_SD,color='k',linestyle='dotted')
        plt.axhline(y=allsources_mean_Arcsec_FCF-allsources_Arcsec_FCF_SD,color='k',linestyle='dotted')
        plt.axhspan(allsources_mean_Arcsec_FCF-allsources_Arcsec_FCF_SD,allsources_mean_Arcsec_FCF+allsources_Arcsec_FCF_SD,alpha=0.5,color='k')
    plt.suptitle('Arcsec FCF, All Sources')

    if tick_spacing<300:
        plt.xlabel('TauTimesAM')
    else:
        plt.xlabel('TauTimesAM')
    plt.ylabel('Arcsec FCF')
    plt.legend(loc='lower left')
    if saveplots == False:
        plt.show()
    else:
        plt.savefig('ArcsecFCF_vs_TauTimesAM_ALLSOURCES_'+ending,dpi=300)
    plt.clf()


if len(epochs_for_comparison) == 13:
    pickle.dump(comparison_dict,open('../catalogues/FCF_AR_PER_EPOCH_'+wave+'.bin','wb'))

#print('\n\nSource: '+source_names[i]+'\n\n')
#print('Now print AR, Peak FCF, and Arcsec FCF info for data and comparison epochs - use comparison_dict\n\n')
