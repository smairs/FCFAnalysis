def comparetau(a1,b1,a2,b2):
    """
    compare_tau_relations.py
    
    This program allows you to compare how two
    different opacity relations will affect the
    applied extinction correction as a function
    of tau*airmass. 
    
    It will specifically indicate typical Grade 1, 
    Grade 2, and Grade 3 weather and show how much '
    the Extinction Correction/FCF differs.
    
    The extinction correction is calculated by
    a simple plane parallel atmospheric correction:
    
    Im = I0 exp(-tau*A)
    
    Where:
    
    Im is the measured intensity
    I0 is the intensity at the top of the atmosphere
    A is the airmass, and
    
    tau is given by:
    
    tau_wavelength = a*(tau225 - b)
    
    where a and b are coefficients derived
    for both 450 and 850 microns. 
    
    At 450 microns, a = 26 , b = 0.012
    At 850 microns, a = 4.6, b = 0.0043
    
    Tau at 225 GHz is used for historical reasons
    because that is what the CSO tipper used to measure.
    
    These are the nominal values, however.
    For more information, see Dempsey et al. (2013).
    
    
    ###################
    
    Input: 
    
    a1 = the "a" coefficient in the first opacity relation
    
    b1 = the "b" coefficient in the first opacity relation
    
    a2 = the "a" coefficient in the second opacity relation
    
    b2 = the "b" coefficient in the second opacity relation
    
    ##################
    
    Output:
    
    2 plots will be displayed: 
    
    The first plots the transmission as a function of tau225*airmass.
    
    The Second plots the first extinction correction divided by the
    second as a function of tau225*airmass, highlighting typical
    Grade 1, Grade 2, and Grade 3 weather.
    
    The extinction correction 1 / extinction correction 2 factor
    is the same as the FCF 1 / FCF 2.
    
    """

    import numpy as np
    import matplotlib.pyplot as plt

    airmassmodel = np.arange(1.0,2.5,0.0006)
    tau225model  = np.arange(0.0,0.25,0.0001)

    transmissionmodel_1 = np.exp(-1.0*a1*(tau225model-b1)*airmassmodel)
    transmissionmodel_2 = np.exp(-1.0*a2*(tau225model-b2)*airmassmodel)

    plt.scatter(tau225model*airmassmodel,transmissionmodel_1,color='b',label='Ext_cor_1')
    plt.scatter(tau225model*airmassmodel,transmissionmodel_2,color='r',label='Ext_cor_2')
    plt.xlabel('Tau225*Airmass')
    plt.ylabel('Transmission')
    plt.suptitle('Ext_cor_1: a = '+str(a1)+', b = '+str(b1)+' -- Ext_cor_2: a = '+str(a2)+', b = '+str(b2))
    plt.legend(loc = 'upper right')
#    plt.show()
    plt.clf()


    plt.scatter(tau225model*airmassmodel,transmissionmodel_1/transmissionmodel_2,color='k')
    plt.axvline(x=0.040*1.2,label='Typical Band 1, Airmass = 1.2',color='b',linestyle='dashed',linewidth=2)
    plt.axvline(x=0.065*1.2,label='Typical Band 2, Airmass = 1.2',color='g',linestyle='dotted',linewidth=2)
    plt.axvline(x=0.100*1.2,label='Typical Band 3, Airmass = 1.2',color='darkgoldenrod',linestyle='dashdot',linewidth=2)
    plt.axvline(x=0.160*1.2,label='Typical Band 4, Airmass = 1.2',color='red',linestyle='solid',linewidth=2)
    plt.axvline(x=0.260*1.2,label='Typical Band 5, Airmass = 1.2',color='magenta',linestyle='dashed',linewidth=2)
    plt.axhline(y=1.0,linestyle = 'dashed', color = 'k', linewidth=2.0)
    plt.xlabel('Tau225*Airmass')
    plt.ylabel('Ext_cor_1/Ext_cor_2')
    plt.suptitle('Ext_cor_1: a = '+str(a1)+', b = '+str(b1)+' -- Ext_cor_2: a = '+str(a2)+', b = '+str(b2))
    plt.legend(loc='upper right')
    plt.show()
    plt.clf()

