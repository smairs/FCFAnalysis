# This comes from running run_summarise*py for both constrained and unconstrained data, then copying the printout lists and placing them here

import numpy as np
import matplotlib.pyplot as plt


TauTimesAirmassBins  = [0.030,0.084,0.138,0.192,0.246,0.30]
bincolours           = ['r','y','g','b','violet']
bincolours450        = ['r','y','k','black','black']

TauTimesAirmassBins  = [0.030,0.084,0.138,0.192,0.246,0.30]

x_before       = [eachTauAM + 1*(TauTimesAirmassBins[ind+1]-TauTimesAirmassBins[ind])/7.0 for ind,eachTauAM in enumerate(TauTimesAirmassBins) if ind<len(TauTimesAirmassBins)-1]
x_after        = [eachTauAM + 2*(TauTimesAirmassBins[ind+1]-TauTimesAirmassBins[ind])/7.0 for ind,eachTauAM in enumerate(TauTimesAirmassBins) if ind<len(TauTimesAirmassBins)-1]

x_night_before = [eachTauAM + 3*(TauTimesAirmassBins[ind+1]-TauTimesAirmassBins[ind])/7.0 for ind,eachTauAM in enumerate(TauTimesAirmassBins) if ind<len(TauTimesAirmassBins)-1]
x_night_after  = [eachTauAM + 4*(TauTimesAirmassBins[ind+1]-TauTimesAirmassBins[ind])/7.0 for ind,eachTauAM in enumerate(TauTimesAirmassBins) if ind<len(TauTimesAirmassBins)-1]

x_EO_before    = [eachTauAM + 5*(TauTimesAirmassBins[ind+1]-TauTimesAirmassBins[ind])/7.0 for ind,eachTauAM in enumerate(TauTimesAirmassBins) if ind<len(TauTimesAirmassBins)-1]
x_EO_after     = [eachTauAM + 6*(TauTimesAirmassBins[ind+1]-TauTimesAirmassBins[ind])/7.0 for ind,eachTauAM in enumerate(TauTimesAirmassBins) if ind<len(TauTimesAirmassBins)-1]

#x_before       =  [0.0385, 0.0925, 0.14650000000000002, 0.2005, 0.2545]
#x_after        =  [0.06550000101250003, 0.11950000101250002, 0.17350000101250002, 0.2275000010125, 0.2815000010125]

#FCF Arcsec unconstrained


FCFA_450_before     =  np.array([4.480021514289366  , 4.613283001239128  , 4.3879347767712815 , 4.490567925415016   ,  2.8596141297103572])
FCFA_450_before_err =  np.array([0.5515094152302374 , 0.7636179892681594 , 0.9436932635237643 , 0.8461284207351825  , 0.3933417172749852])
FCFA_450_after      =  np.array([3.6852783892088334 , 3.2698885426025197 , 3.5506903796190414 , 3.3815773197673096  , np.nan])
FCFA_450_after_err  =  np.array([0.6071869453981543 , 0.7154814561889367 , 0.8139860725084659 , 0.48844473268050076 , np.nan])
FCFA_850_before     =  np.array([2.2505009498411246 , 2.301967584079901  , 2.314188715066051  , 2.335518424207868   , 2.3134430913879])
FCFA_850_before_err =  np.array([0.1035392886369965 , 0.10193853490445844, 0.11029437015532519, 0.11347750440369536 , 0.09502587648927598])
FCFA_850_after      =  np.array([2.072459502547625  , 2.101146540916701  , 2.1797017037619533 , 2.1731020392562286  , 2.135903919286526])
FCFA_850_after_err  =  np.array([0.08193012287588393, 0.10224306348412762, 0.1112920784098462 , 0.11040409038477567 , 0.18519508568985837])


#FCF Peak unconstrained


FCFP_450_before     =  np.array([551.4496985220447, 548.7064880012377, 413.0079321929186, 519.3488823083178, 336.43925653690223])
FCFP_450_before_err =  np.array([112.75247843116892, 112.30529993254616, 124.16059514554098, 87.15493534441369, 82.69316671793557])
FCFP_450_after      =  np.array([502.60916757197435, 449.4276131524905, 416.6510953786635, 449.297634377957, np.nan])
FCFP_450_after_err  =  np.array([103.7872962480516, 118.80251902316198, 82.63051063905644, 68.00821839167898, np.nan])
FCFP_850_before     =  np.array([530.4599388602943, 545.1259806782343, 553.7137365121889, 558.5551404751509, 567.2029418612623])
FCFP_850_before_err =  np.array([38.273527180059695, 40.34576205059649, 47.39991762293759, 44.176122958784994, 54.527590631965026])
FCFP_850_after      =  np.array([511.12547204734454, 521.7016420512806, 532.931743409866, 545.5471222618037, 533.0239795394134])
FCFP_850_after_err  =  np.array([29.44915614118726, 39.92118802916071, 45.58876515772759, 49.434354824728395, 51.87552803266527])


#FCF Arcsec AM<1.9, FWHM<20, Time between 21:00 and 06:00


FCFA_night_450_before     =  np.array([4.39107652369028, 4.671719831108589, 4.091533393091509, 4.611314673372182, np.nan])
FCFA_night_450_before_err =  np.array([0.4615610774593497, 0.7423663478625906, 0.8617624343452099, 1.0879951104664451, np.nan])
FCFA_night_450_after      =  np.array([3.699074328480581, 3.3808388997124026, 3.624271974317801, 3.3208847013793195, np.nan])
FCFA_night_450_after_err  =  np.array([0.4748320561510639, 0.45125203574861217, 0.804650818523864, 0.34983643389734004, np.nan])
FCFA_night_850_before     =  np.array([2.248199310387599, 2.296300266877474, 2.3000947930522084, 2.341233152749907, 2.2432885905535453])
FCFA_night_850_before_err =  np.array([0.09088824666922223, 0.08353019837904832, 0.0939860857192397, 0.1088834371376563, 0.08885096683888319])
FCFA_night_850_after      =  np.array([2.0613077256288057, 2.1078980948710386, 2.1567244919014223, 2.1657614017316424, 2.0013166452612032])
FCFA_night_850_after_err  =  np.array([0.07075341316960714, 0.07269984549377435, 0.10725179410846326, 0.07578682620695737, 0.2291024442121897])


#FCF Peak AM<1.9, FWHM<20, Time between 21:00 and 06:00

FCFP_night_450_before     =  np.array([523.7893722384534, 534.6015852819148, 433.2723725839703, 508.1512886277077, np.nan])
FCFP_night_450_before_err =  np.array([86.16903438767687, 94.77456612865193, 79.80340661141638, 75.82081866758907, np.nan])
FCFP_night_450_after      =  np.array([504.7546674798686, 457.4397293260542, 352.06788235877855, 446.838587202517, np.nan])
FCFP_night_450_after_err  =  np.array([76.71474537258351, 75.48383831299995, 60.61298606784046, 47.74771732270221, np.nan])
FCFP_night_850_before     =  np.array([522.1823114294436, 538.4238292839461, 546.74509588377, 554.0854744033287, 556.8528576543304])
FCFP_night_850_before_err =  np.array([26.778216492455424, 28.614799646391383, 36.0359723937655, 37.34352018318975, 37.1987830236141])
FCFP_night_850_after      =  np.array([508.2035451906274, 523.5278590016412, 516.0748034792449, 538.9131890080275, 529.0395391419103])
FCFP_night_850_after_err  =  np.array([24.59560358059414, 26.864636589369056, 26.12884239242361, 33.41365989028923, 46.5490134802785])

#FCF Arcsec AM<1.9, FWHM<20, Time between 06:00 and 12:00


FCFA_EO_450_before     =  np.array([4.431592792047677, 4.776974515107194, 4.675015176518045, 4.280551965796654, 2.8596141297103572])
FCFA_EO_450_before_err =  np.array([0.5180659210873274, 0.7751491746081577, 0.8360514693516783, 0.7074598771984784, 0.3933417172749852])
FCFA_EO_450_after      =  np.array([np.nan, 3.917750128052787, 3.168478321789731, np.nan, np.nan])
FCFA_EO_450_after_err  =  np.array([np.nan, 0.4267594338627099, 0.4641627866302031, np.nan, np.nan])
FCFA_EO_850_before     =  np.array([2.2468826663862385, 2.314590654760696, 2.328172862285431, 2.297337329930778, 2.21642045517959])
FCFA_EO_850_before_err =  np.array([0.09690055869393857, 0.09924602360130944, 0.1060584549451298, 0.11240322722146429, 0.09308181828252467])
FCFA_EO_850_after      =  np.array([2.0666499141052093, 2.1133998272622736, 2.1794402894192246, 2.1942105580922204, np.nan])
FCFA_EO_850_after_err  =  np.array([0.1094682294523308, 0.08044477644827853, 0.10382563404243252, 0.10555046793905928, np.nan])


#FCF Peak AM<1.9, FWHM<20, Time between 06:00 and 12:00


FCFP_EO_450_before     =  np.array([554.6214701436161, 566.4773266852328, 508.9591806536469, 607.8191635540161, 336.43925653690223])
FCFP_EO_450_before_err =  np.array([116.07845593245194, 117.16909414458149, 102.48439170325464, 107.41194868364002, 82.69316671793557])
FCFP_EO_450_after      =  np.array([np.nan, 446.37681824003226, 468.3789116823514, np.nan, np.nan])
FCFP_EO_450_after_err  =  np.array([np.nan, 84.96106446893612, 99.86330999021796, np.nan, np.nan])
FCFP_EO_850_before     =  np.array([530.2308959741391, 541.8366487951794, 562.1941670675862, 543.7648867112514, 560.5008449884236])
FCFP_EO_850_before_err =  np.array([42.09964315276942, 41.016018922138116, 48.376690159278176, 52.93826524951387, 46.60307224947925])
FCFP_EO_850_after      =  np.array([512.7618732991328, 509.0443820205881, 535.6232301947732, 578.539218872893, np.nan])
FCFP_EO_850_after_err  =  np.array([43.96379548727547, 32.82428411296843, 52.24962254921132, 26.34758527548278, np.nan])

#############################
#############################
##### 450 FCF PEAK PLOT #####
#############################
#############################

TauTimesAirmassDummy = 0
for i in range(len(TauTimesAirmassBins)-1):
    plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours450[TauTimesAirmassDummy])
    TauTimesAirmassDummy = TauTimesAirmassDummy+1

plt.errorbar(x_before,FCFP_450_before,yerr=FCFP_450_before_err,linestyle='none',marker='o',color='gray',label='Before SMU: Unconstr.',markersize=10)
plt.errorbar(x_after,FCFP_450_after,yerr=FCFP_450_after_err,linestyle='none',marker='o',color='k',label='After SMU: Unconstr',markersize=10)
plt.errorbar(x_night_before,FCFP_night_450_before,yerr=FCFP_night_450_before_err,linestyle='none',marker='s',color='gray',label='21-06',markersize=10)
plt.errorbar(x_night_after,FCFP_night_450_after,yerr=FCFP_night_450_after_err,linestyle='none',marker='s',color='k',markersize=10)
plt.errorbar(x_EO_before,FCFP_EO_450_before,yerr=FCFP_EO_450_before_err,linestyle='none',marker='x',color='gray',label='EO',markersize=10)
plt.errorbar(x_EO_after,FCFP_EO_450_after,yerr=FCFP_EO_450_after_err,linestyle='none',marker='x',color='k',markersize=10)
plt.xticks(TauTimesAirmassBins)
plt.legend(loc='upper right')
plt.ylabel('Peak FCF, 450 microns (Jy/beam/pW)')
plt.xlabel('Tau225 x Airmass')
#plt.show()
plt.savefig('FCFP_by_TauxAM_450_TimeOfNight.png',format='png',dpi=300)
plt.clf()

###############################
###############################
##### 450 FCF Arcsec PLOT #####
###############################
###############################

TauTimesAirmassDummy = 0
for i in range(len(TauTimesAirmassBins)-1):
    plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours450[TauTimesAirmassDummy])
    TauTimesAirmassDummy = TauTimesAirmassDummy+1

plt.errorbar(x_before,FCFA_450_before,yerr=FCFA_450_before_err,linestyle='none',marker='o',color='gray',label='Before SMU: Unconstr.',markersize=10)
plt.errorbar(x_after,FCFA_450_after,yerr=FCFA_450_after_err,linestyle='none',marker='o',color='k',label='After SMU: Unconstr',markersize=10)
plt.errorbar(x_night_before,FCFA_night_450_before,yerr=FCFA_night_450_before_err,linestyle='none',marker='s',color='gray',label='21-06',markersize=10)
plt.errorbar(x_night_after,FCFA_night_450_after,yerr=FCFA_night_450_after_err,linestyle='none',marker='s',color='k',markersize=10)
plt.errorbar(x_EO_before,FCFA_EO_450_before,yerr=FCFA_EO_450_before_err,linestyle='none',marker='x',color='gray',label='EO',markersize=10)
plt.errorbar(x_EO_after,FCFA_EO_450_after,yerr=FCFA_EO_450_after_err,linestyle='none',marker='x',color='k',markersize=10)
plt.xticks(TauTimesAirmassBins)
plt.legend(loc='upper right')
plt.ylabel('Arcsec FCF, 450 microns (Jy/arcsec^2/pW)')
plt.xlabel('Tau225 x Airmass')
#plt.show()
plt.savefig('FCFA_by_TauxAM_450_TimeOfNight.png',format='png',dpi=300)
plt.clf()

#############################
#############################
##### 850 FCF PEAK PLOT #####
#############################
#############################

TauTimesAirmassDummy = 0
for i in range(len(TauTimesAirmassBins)-1):
    plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours[TauTimesAirmassDummy])
    TauTimesAirmassDummy = TauTimesAirmassDummy+1

plt.errorbar(x_before,FCFP_850_before,yerr=FCFP_850_before_err,linestyle='none',marker='o',color='gray',label='Before NF: Unconstr.',markersize=10)
plt.errorbar(x_after,FCFP_850_after,yerr=FCFP_850_after_err,linestyle='none',marker='o',color='k',label='After NF: Unconstr',markersize=10)
plt.errorbar(x_night_before,FCFP_night_850_before,yerr=FCFP_night_850_before_err,linestyle='none',marker='s',color='gray',label='21-06',markersize=10)
plt.errorbar(x_night_after,FCFP_night_850_after,yerr=FCFP_night_850_after_err,linestyle='none',marker='s',color='k',markersize=10)
plt.errorbar(x_EO_before,FCFP_EO_850_before,yerr=FCFP_EO_850_before_err,linestyle='none',marker='x',color='gray',label='EO',markersize=10)
plt.errorbar(x_EO_after,FCFP_EO_850_after,yerr=FCFP_EO_850_after_err,linestyle='none',marker='x',color='k',markersize=10)
plt.xticks(TauTimesAirmassBins)
plt.legend(loc='upper right')
plt.ylabel('Peak FCF, 850 microns (Jy/beam/pW)')
plt.xlabel('Tau225 x Airmass')
#plt.show()
plt.savefig('FCFP_by_TauxAM_850_TimeOfNight.png',format='png',dpi=300)
plt.clf()

###############################
###############################
##### 850 FCF Arcsec PLOT #####
###############################
###############################

TauTimesAirmassDummy = 0
for i in range(len(TauTimesAirmassBins)-1):
    plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours[TauTimesAirmassDummy])
    TauTimesAirmassDummy = TauTimesAirmassDummy+1

plt.errorbar(x_before,FCFA_850_before,yerr=FCFA_850_before_err,linestyle='none',marker='o',color='gray',label='Before NF: Unconstr.',markersize=10)
plt.errorbar(x_after,FCFA_850_after,yerr=FCFA_850_after_err,linestyle='none',marker='o',color='k',label='After NF: Unconstr',markersize=10)
plt.errorbar(x_night_before,FCFA_night_850_before,yerr=FCFA_night_850_before_err,linestyle='none',marker='s',color='gray',label='21-06',markersize=10)
plt.errorbar(x_night_after,FCFA_night_850_after,yerr=FCFA_night_850_after_err,linestyle='none',marker='s',color='k',markersize=10)
plt.errorbar(x_EO_before,FCFA_EO_850_before,yerr=FCFA_EO_850_before_err,linestyle='none',marker='x',color='gray',label='EO',markersize=10)
plt.errorbar(x_EO_after,FCFA_EO_850_after,yerr=FCFA_EO_850_after_err,linestyle='none',marker='x',color='k',markersize=10)
plt.xticks(TauTimesAirmassBins)
plt.legend(loc='lower left')
plt.ylabel('Arcsec FCF, 850 microns (Jy/arcsec^2/pW)')
plt.xlabel('Tau225 x Airmass')
#plt.show()
plt.savefig('FCFA_by_TauxAM_850_TimeOfNight.png',format='png',dpi=300)
plt.clf()


#################################
#################################
######### UNCERTAINTIES #########
#################################
#################################

# 450 FCF Peak

TauTimesAirmassDummy = 0
for i in range(len(TauTimesAirmassBins)-1):
    plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours450[TauTimesAirmassDummy])
    TauTimesAirmassDummy = TauTimesAirmassDummy+1

plt.scatter(x_before,100*FCFP_450_before_err/FCFP_450_before,marker='o',color='gray',label='Before SMU: Unconstr.',s=40)
plt.scatter(x_after,100*FCFP_450_after_err/FCFP_450_after,marker='o',color='k',label='After SMU: Unconstr.',s=40)
plt.scatter(x_night_before,100*FCFP_night_450_before_err/FCFP_night_450_before,marker='s',color='gray',label='21-06',s=40)
plt.scatter(x_night_after,100*FCFP_night_450_after_err/FCFP_night_450_after,marker='s',color='k',s=40)
plt.scatter(x_EO_before,100*FCFP_EO_450_before_err/FCFP_EO_450_before,marker='x',color='gray',label='EO',s=40)
plt.scatter(x_EO_after,100*FCFP_EO_450_after_err/FCFP_EO_450_after,marker='x',color='k',s=40)
plt.ylabel('450 Micron FCF Peak Uncertainty (%)')
plt.xlabel('Tau225 x Airmass')
plt.legend(loc = 'upper left')
#plt.show()
plt.savefig('FCFP_err_by_TauxAM_450_TimeOfNight.png',format='png',dpi=300)
plt.clf()

# 450 FCF Arcsec

TauTimesAirmassDummy = 0
for i in range(len(TauTimesAirmassBins)-1):
    plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours450[TauTimesAirmassDummy])
    TauTimesAirmassDummy = TauTimesAirmassDummy+1

plt.scatter(x_before,100*FCFA_450_before_err/FCFA_450_before,marker='o',color='gray',label='Before SMU: Unconstr.',s=40)
plt.scatter(x_after,100*FCFA_450_after_err/FCFA_450_after,marker='o',color='k',label='After SMU: Unconstr.',s=40)
plt.scatter(x_night_before,100*FCFA_night_450_before_err/FCFA_night_450_before,marker='s',color='gray',label='21-06',s=40)
plt.scatter(x_night_after,100*FCFA_night_450_after_err/FCFA_night_450_after,marker='s',color='k',s=40)
plt.scatter(x_EO_before,100*FCFA_EO_450_before_err/FCFA_EO_450_before,marker='x',color='gray',label='EO',s=40)
plt.scatter(x_EO_after,100*FCFA_EO_450_after_err/FCFA_EO_450_after,marker='x',color='k',s=40)
plt.ylabel('450 Micron FCF Arcsec Uncertainty (%)')
plt.xlabel('Tau225 x Airmass')
plt.legend(loc = 'upper left')
#plt.show()
plt.savefig('FCFA_err_by_TauxAM_450_TimeOfNight.png',format='png',dpi=300)
plt.clf()

# 850 FCF Peak

TauTimesAirmassDummy = 0
for i in range(len(TauTimesAirmassBins)-1):
    plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours[TauTimesAirmassDummy])
    TauTimesAirmassDummy = TauTimesAirmassDummy+1

plt.scatter(x_before,100*FCFP_850_before_err/FCFP_850_before,marker='o',color='gray',label='Before NF: Unconstr.',s=40)
plt.scatter(x_after,100*FCFP_850_after_err/FCFP_850_after,marker='o',color='k',label='After NF: Unconstr.',s=40)
plt.scatter(x_night_before,100*FCFP_night_850_before_err/FCFP_night_850_before,marker='s',color='gray',label='21-06',s=40)
plt.scatter(x_night_after,100*FCFP_night_850_after_err/FCFP_night_850_after,marker='s',color='k',s=40)
plt.scatter(x_EO_before,100*FCFP_EO_850_before_err/FCFP_EO_850_before,marker='x',color='gray',label='EO',s=40)
plt.scatter(x_EO_after,100*FCFP_EO_850_after_err/FCFP_EO_850_after,marker='x',color='k',s=40)
plt.ylabel('850 Micron FCF Peak Uncertainty (%)')
plt.xlabel('Tau225 x Airmass')
plt.legend(loc = 'upper left')
#plt.show()
plt.savefig('FCFP_err_by_TauxAM_850_TimeOfNight.png',format='png',dpi=300)
plt.clf()

# 850 FCF Arcsec

TauTimesAirmassDummy = 0
for i in range(len(TauTimesAirmassBins)-1):
    plt.axvspan(TauTimesAirmassBins[i],TauTimesAirmassBins[i+1],alpha=0.3,color=bincolours[TauTimesAirmassDummy])
    TauTimesAirmassDummy = TauTimesAirmassDummy+1

plt.scatter(x_before,100*FCFA_850_before_err/FCFA_850_before,marker='o',color='gray',label='Before NF: Unconstr.',s=40)
plt.scatter(x_after,100*FCFA_850_after_err/FCFA_850_after,marker='o',color='k',label='After NF: Unconstr.',s=40)
plt.scatter(x_night_before,100*FCFA_night_850_before_err/FCFA_night_850_before,marker='s',color='gray',label='21-06',s=40)
plt.scatter(x_night_after,100*FCFA_night_850_after_err/FCFA_night_850_after,marker='s',color='k',s=40)
plt.scatter(x_EO_before,100*FCFA_EO_850_before_err/FCFA_EO_850_before,marker='x',color='gray',label='EO',s=40)
plt.scatter(x_EO_after,100*FCFA_EO_850_after_err/FCFA_EO_850_after,marker='x',color='k',s=40)
plt.ylabel('850 Micron FCF Arcsec Uncertainty (%)')
plt.xlabel('Tau225 x Airmass')
plt.legend(loc='upper left')
#plt.show()
plt.savefig('FCFA_err_by_TauxAM_850_TimeOfNight.png',format='png',dpi=300)
plt.clf()


########################################
########################################
# Print Out Tables (Visual, Not Latex) #
########################################
########################################

###########
###########
# 850
###########
###########

titles     = ['Tau225xAM     ','Pre-NF, ALL   ','Pre-NF, 21-06 ','Pre-NF, 06-12 ','Post-NF, ALL  ','Post-NF, 21-06','Post-NF, 06-12']
tautimesAM = ['0.030 -> 0.084','0.084 -> 0.138','0.138 -> 0.192','0.192 -> 0.246','0.246 -> 0.300']

#### 850 Peak ####

data       = [titles] + list(zip(tautimesAM,["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFP_850_before_err[i]/x) for i,x in enumerate(FCFP_850_before)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFP_night_850_before_err[i]/x) for i,x in enumerate(FCFP_night_850_before)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFP_EO_850_before_err[i]/x) for i,x in enumerate(FCFP_EO_850_before)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFP_850_after_err[i]/x) for i,x in enumerate(FCFP_850_after)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFP_night_850_after_err[i]/x) for i,x in enumerate(FCFP_night_850_after)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFP_EO_850_after_err[i]/x) for i,x in enumerate(FCFP_EO_850_after)]))

print('\n\n850 microns, FCF Peak\n\n')

for i, d in enumerate(data):
    line = ' | '.join(str(x).ljust(12) for x in d)
    print(line)
    if i == 0:
        print('-' * len(line))
print('\n\n')

#### 850 Arcsec ####

data = [titles] + list(zip(tautimesAM,["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFA_850_before_err[i]/x) for i,x in enumerate(FCFA_850_before)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFA_night_850_before_err[i]/x) for i,x in enumerate(FCFA_night_850_before)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFA_EO_850_before_err[i]/x) for i,x in enumerate(FCFA_EO_850_before)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFA_850_after_err[i]/x) for i,x in enumerate(FCFA_850_after)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFA_night_850_after_err[i]/x) for i,x in enumerate(FCFA_night_850_after)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFA_EO_850_after_err[i]/x) for i,x in enumerate(FCFA_EO_850_after)]))

print('\n\n850 microns, FCF Arcsec\n\n')

for i, d in enumerate(data):
    line = ' | '.join(str(x).ljust(12) for x in d)
    print(line)
    if i == 0:
        print('-' * len(line))
print('\n\n')

#############
#############
# 450
#############
#############

titles     = ['Tau225xAM     ','Pre-SM, ALL   ','Pre-SM, 21-06 ','Pre-SM, 06-12 ','Post-SM, ALL  ','Post-SM, 21-06','Post-SM, 06-12']
tautimesAM = ['0.030 -> 0.084','0.084 -> 0.138','0.138 -> 0.192','0.192 -> 0.246','0.246 -> 0.300']

#### 450 Peak ####

data       = [titles] + list(zip(tautimesAM,["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFP_450_before_err[i]/x) for i,x in enumerate(FCFP_450_before)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFP_night_450_before_err[i]/x) for i,x in enumerate(FCFP_night_450_before)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFP_EO_450_before_err[i]/x) for i,x in enumerate(FCFP_EO_450_before)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFP_450_after_err[i]/x) for i,x in enumerate(FCFP_450_after)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFP_night_450_after_err[i]/x) for i,x in enumerate(FCFP_night_450_after)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFP_EO_450_after_err[i]/x) for i,x in enumerate(FCFP_EO_450_after)]))

print('\n\n450 microns, FCF Peak\n\n')

for i, d in enumerate(data):
    line = ' | '.join(str(x).ljust(12) for x in d)
    print(line)
    if i == 0:
        print('-' * len(line))
print('\n\n')

#### 450 Arcsec ####

data = [titles] + list(zip(tautimesAM,["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFA_450_before_err[i]/x) for i,x in enumerate(FCFA_450_before)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFA_night_450_before_err[i]/x) for i,x in enumerate(FCFA_night_450_before)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFA_EO_450_before_err[i]/x) for i,x in enumerate(FCFA_EO_450_before)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFA_450_after_err[i]/x) for i,x in enumerate(FCFA_450_after)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFA_night_450_after_err[i]/x) for i,x in enumerate(FCFA_night_450_after)],["{:6.2f}".format(x)+' +/-'+"{:4.1f}".format(100*FCFA_EO_450_after_err[i]/x) for i,x in enumerate(FCFA_EO_450_after)]))

print('\n\n450 microns, FCF Arcsec\n\n')

for i, d in enumerate(data):
    line = ' | '.join(str(x).ljust(12) for x in d)
    print(line)
    if i == 0:
        print('-' * len(line))
print('\n\n')

