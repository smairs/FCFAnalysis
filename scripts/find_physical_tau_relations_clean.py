def phystau(a1,a2,astep,b1,b2,bstep,wave,physthresh,param_num):                                                                                                                                                  
    '''
    You supply a range of "a" and "b" values to convert the tau225 values into 
    tau345 (850 microns) and tau666 (450 microns):

    tauX = a(tau225 - b)
    
    *step are the grid steps in "ab parameterspace". 
    
    Wave is given as an integer (450 or 850).
    
    We can look at any number of tau*airmass combinations. param_num refers to the number of tau
    and airmass measurements individually, but since they are multiplied
    param_num = 100 -> 10,000 iterations
    
    physthresh is the threshold of expected - measured transmission values for what we are considering physical. 
    physthresh = 0.02 -> 0.05
    
    Then, the code genrates a transmission for every combination of a and b value, 
    tau225 value between 0.03 and 0.20 (you choose the step size), and 
    airmass between 1.0 and 2.0 (you choose the step size). The code then compares this 
    calculated transmission to the expected transmission for all (tau x airmass) values and 
    all combinations of a and b. We can then determine which a and b combinations give us physical 
    transmissions across all the weather bands we are interested in. If less than 2% of the data points
    lie outside the physical threshold - I mark them as "good".
    '''

    import numpy as np
    import matplotlib.pyplot as plt
    import os                      

    # This Tau range covers good grade 1 to the grade 4/5 boundary
    taustart  = 0.03
    tauend    = 0.20

    # These are the most typical airmasses at which we observe- too much over 2.0
    # and the observations become much less certain
    AMstart   = 1.0
    AMend     = 2.0

    # We can look at any number of tau*airmass combinations, but since they are multiplied
    # 100 -> 10,000 iterations
    param_num = 100.0

    tau225s   = np.arange(taustart,tauend+((tauend-taustart)/param_num),(tauend-taustart)/param_num)
    airmasses = np.arange(AMstart,AMend+((AMend-AMstart)/param_num),(AMend-AMstart)/param_num)

    # Import the codes that interpolate over the Transmission versus PWV curve
    # and then give you a transmission value based on an input tau225
    from TauRelAnalysis_20171215 import FitCSOTransvsPWV,CSOtrans

    # Construct the Transmission verusus PWV interpolation
    TransvsPWV = FitCSOTransvsPWV(int(wave))

    # Define a list to "catch" physical solutions
    phys_solutions = []

    os.system('mkdir Figures/simulated_physplots')
    dummy=-1
    for a in np.arange(a1,a2+astep/2.0,astep):
        for b in np.arange(b1,b2+bstep/2.0,bstep):
            dummy = dummy+1
            print('\n Working on Opacity Relation #'+str(dummy+1)+' out of '+str(len(np.arange(a1,a2+astep/2.0,astep))*len(np.arange(b1,b2+bstep/2.0,bstep)))+'\n')
            # Measured Transmission, expected transmission (from CSO), tau*AM combination
            measured_trans = []
            expected_trans = []
            taus_times_am  = []
            for eachtau225 in range(len(tau225s)):
                for eachairmass in range(len(airmasses)):
                    expected_trans.append(CSOtrans(int(wave),tau225s[eachtau225],TransvsPWV)**(airmasses[eachairmass])) #CSOtrans gives you the ZENITH transmission - we want exp(-tau*A) though, so raise Trans_zen to exponent of A
                    tau = a*(tau225s[eachtau225]-b)
                    measured_trans.append(np.exp(-1.0*tau*airmasses[eachairmass]))
                    taus_times_am.append(tau225s[eachtau225]*airmasses[eachairmass])

            expected_trans = np.array(expected_trans)
            measured_trans = np.array(measured_trans)

            # The difference between the expected and the measured transmission (for percent, multiply by 100)
            delta_trans    = expected_trans - measured_trans

            # If less than 2 percent of the points are outside of the physical threshold, mark it as physical
            if len(delta_trans[np.where(abs(delta_trans)>physthresh)])/float(len(delta_trans))<0.02:
                phys_solutions.append((a,b))
                plt.scatter(taus_times_am,delta_trans)
                plt.xlabel('Tau225'+' x Airmass')
                plt.ylabel('Expected Trans - Measured Trans')
                plt.suptitle('a = '+str(a)+' , '+' b = '+str(b))
                plt.savefig('physplot_'+str(wave)+'_'+str(a)+'_'+str(b)+'.pdf',format='pdf')
                plt.clf()
                os.system('mv physplot_'+str(wave)+'_'+str(a)+'_'+str(b)+'.pdf Figures/simulated_physplots/')
    
    os.system('rm -rf tmpADA*') if len(list(glob.glob('tmpADA*'))) > 0
            
    # The code returns (a,b) combinations which passed the "physical" test
    return(phys_solutions)
