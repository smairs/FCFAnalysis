def summarise_epochs(infofilelist, Tau225TimesAirmass_bins, badsourcelists):
    '''
    infofilelist should be of the form:

    ['FCF_AR_PER_EPOCH_450.bin','FCF_AR_PER_EPOCH_850.bin']

    where the lists are one of:

    'FCF_AR_PER_EPOCH_450_EE.bin'
    'FCF_AR_PER_EPOCH_850_EE.bin'
   
    'FCF_AR_PER_EPOCH_450_SS.bin'
    'FCF_AR_PER_EPOCH_850_SS.bin'
   
    'FCF_AR_PER_EPOCH_450_EO.bin'
    'FCF_AR_PER_EPOCH_850_EO.bin'
   
    'FCF_AR_PER_EPOCH_450.bin'
    'FCF_AR_PER_EPOCH_850.bin'

    These are different files containing
    FCF and metadata information on a set
    of observations. 

    EE = Early Evening
    EO = Extended Observing
    SS = Sweet Spot (21-03)

    Just .bin = all data.
 
    Tau225TimesAirmass_bins = [0.030,0.084,0.138,0.192,0.246,0.30] 
    (need at least a left and right edge of bin)

    bad source lists - try deriving FCFs by including or excluding certain sources!

    badsourcelists = [['Arp220','NEPTUNE','MARS'],['Arp220','NEPTUNE','MARS', 'CRL618']]
    '''


    from summarise_epochs import summarise_epochs
    from plot_by_tauxam import plot_by_tauxairmass_2lists
    import numpy as np
    
    for infofile in ['FCF_AR_PER_EPOCH_450.bin','FCF_AR_PER_EPOCH_850.bin']:
        baddummy=-1
        if infofile.split('_')[-1]=='450.bin':
            FCFP_450_before_SMUtrouble1      = []
            FCFP_err_450_before_SMUtrouble1  = []
            FCFP_450_after_SMUtrouble1       = []
            FCFP_err_450_after_SMUtrouble1   = []
            FCFA_450_before_SMUtrouble1      = []
            FCFA_err_450_before_SMUtrouble1  = []
            FCFA_450_after_SMUtrouble1       = []
            FCFA_err_450_after_SMUtrouble1   = []
            FCFP_450_before_SMUtrouble2      = []
            FCFP_err_450_before_SMUtrouble2  = []
            FCFP_450_after_SMUtrouble2       = []
            FCFP_err_450_after_SMUtrouble2   = []
            FCFA_450_before_SMUtrouble2      = []
            FCFA_err_450_before_SMUtrouble2  = []
            FCFA_450_after_SMUtrouble2       = []
            FCFA_err_450_after_SMUtrouble2   = []
        else:    
            FCFP_850_before_newfilters1     = []
            FCFP_err_850_before_newfilters1 = []
            FCFP_850_after_newfilters1      = []
            FCFP_err_850_after_newfilters1  = []
            FCFA_850_before_newfilters1     = []
            FCFA_err_850_before_newfilters1 = []
            FCFA_850_after_newfilters1      = []
            FCFA_err_850_after_newfilters1  = []
            FCFP_850_before_newfilters2     = []
            FCFP_err_850_before_newfilters2 = []
            FCFP_850_after_newfilters2      = []
            FCFP_err_850_after_newfilters2  = []
            FCFA_850_before_newfilters2     = []
            FCFA_err_850_before_newfilters2 = []
            FCFA_850_after_newfilters2      = []
            FCFA_err_850_after_newfilters2  = []
        for badsourcelist in badsourcelists:
            baddummy=baddummy+1
            for i in range(len(Tau225TimesAirmass_bins)-1):
                FCFP_before,FCFP_err_before,FCFP_after,FCFP_err_after,FCFA_before,FCFA_err_before,FCFA_after,FCFA_err_after = summarise_epochs(infofile,[Tau225TimesAirmass_bins[i],Tau225TimesAirmass_bins[i+1]],badsourcelist)
                if infofile.split('_')[-1]=='450.bin':
                    if baddummy == 0:
                        FCFP_450_before_SMUtrouble1.append(FCFP_before)
                        FCFP_err_450_before_SMUtrouble1.append(FCFP_err_before)
                        FCFP_450_after_SMUtrouble1.append(FCFP_after)
                        FCFP_err_450_after_SMUtrouble1.append(FCFP_err_after)
                        FCFA_450_before_SMUtrouble1.append(FCFA_before)
                        FCFA_err_450_before_SMUtrouble1.append(FCFA_err_before)
                        FCFA_450_after_SMUtrouble1.append(FCFA_after)
                        FCFA_err_450_after_SMUtrouble1.append(FCFA_err_after)
                    else:
                        FCFP_450_before_SMUtrouble2.append(FCFP_before)
                        FCFP_err_450_before_SMUtrouble2.append(FCFP_err_before)
                        FCFP_450_after_SMUtrouble2.append(FCFP_after)
                        FCFP_err_450_after_SMUtrouble2.append(FCFP_err_after)
                        FCFA_450_before_SMUtrouble2.append(FCFA_before)
                        FCFA_err_450_before_SMUtrouble2.append(FCFA_err_before)
                        FCFA_450_after_SMUtrouble2.append(FCFA_after)
                        FCFA_err_450_after_SMUtrouble2.append(FCFA_err_after)
                else:
                    if baddummy == 0:
                        FCFP_850_before_newfilters1.append(FCFP_before)
                        FCFP_err_850_before_newfilters1.append(FCFP_err_before)
                        FCFP_850_after_newfilters1.append(FCFP_after)
                        FCFP_err_850_after_newfilters1.append(FCFP_err_after)
                        FCFA_850_before_newfilters1.append(FCFA_before)
                        FCFA_err_850_before_newfilters1.append(FCFA_err_before)
                        FCFA_850_after_newfilters1.append(FCFA_after)
                        FCFA_err_850_after_newfilters1.append(FCFA_err_after)
                    else:
                        FCFP_850_before_newfilters2.append(FCFP_before)
                        FCFP_err_850_before_newfilters2.append(FCFP_err_before)
                        FCFP_850_after_newfilters2.append(FCFP_after)
                        FCFP_err_850_after_newfilters2.append(FCFP_err_after)
                        FCFA_850_before_newfilters2.append(FCFA_before)
                        FCFA_err_850_before_newfilters2.append(FCFA_err_before)
                        FCFA_850_after_newfilters2.append(FCFA_after)
                        FCFA_err_850_after_newfilters2.append(FCFA_err_after) 
    
    FCFP_450_before_SMUtrouble1      = np.array(FCFP_450_before_SMUtrouble1)
    FCFP_err_450_before_SMUtrouble1  = np.array(FCFP_err_450_before_SMUtrouble1)
    FCFP_450_after_SMUtrouble1       = np.array(FCFP_450_after_SMUtrouble1)
    FCFP_err_450_after_SMUtrouble1   = np.array(FCFP_err_450_after_SMUtrouble1)
    FCFA_450_before_SMUtrouble1      = np.array(FCFA_450_before_SMUtrouble1)
    FCFA_err_450_before_SMUtrouble1  = np.array(FCFA_err_450_before_SMUtrouble1)
    FCFA_450_after_SMUtrouble1       = np.array(FCFA_450_after_SMUtrouble1)
    FCFA_err_450_after_SMUtrouble1   = np.array(FCFA_err_450_after_SMUtrouble1)
    FCFP_850_before_newfilters1      = np.array(FCFP_850_before_newfilters1)
    FCFP_err_850_before_newfilters1  = np.array(FCFP_err_850_before_newfilters1)
    FCFP_850_after_newfilters1       = np.array(FCFP_850_after_newfilters1)
    FCFP_err_850_after_newfilters1   = np.array(FCFP_err_850_after_newfilters1)
    FCFA_850_before_newfilters1      = np.array(FCFA_850_before_newfilters1)
    FCFA_err_850_before_newfilters1  = np.array(FCFA_err_850_before_newfilters1)
    FCFA_850_after_newfilters1       = np.array(FCFA_850_after_newfilters1)
    FCFA_err_850_after_newfilters1   = np.array(FCFA_err_850_after_newfilters1)
    
    FCFP_450_before_SMUtrouble2      = np.array(FCFP_450_before_SMUtrouble2)
    FCFP_err_450_before_SMUtrouble2  = np.array(FCFP_err_450_before_SMUtrouble2)
    FCFP_450_after_SMUtrouble2       = np.array(FCFP_450_after_SMUtrouble2)
    FCFP_err_450_after_SMUtrouble2   = np.array(FCFP_err_450_after_SMUtrouble2)
    FCFA_450_before_SMUtrouble2      = np.array(FCFA_450_before_SMUtrouble2)
    FCFA_err_450_before_SMUtrouble2  = np.array(FCFA_err_450_before_SMUtrouble2)
    FCFA_450_after_SMUtrouble2       = np.array(FCFA_450_after_SMUtrouble2)
    FCFA_err_450_after_SMUtrouble2   = np.array(FCFA_err_450_after_SMUtrouble2)
    FCFP_850_before_newfilters2      = np.array(FCFP_850_before_newfilters2)
    FCFP_err_850_before_newfilters2  = np.array(FCFP_err_850_before_newfilters2)
    FCFP_850_after_newfilters2       = np.array(FCFP_850_after_newfilters2)
    FCFP_err_850_after_newfilters2   = np.array(FCFP_err_850_after_newfilters2)
    FCFA_850_before_newfilters2      = np.array(FCFA_850_before_newfilters2)
    FCFA_err_850_before_newfilters2  = np.array(FCFA_err_850_before_newfilters2)
    FCFA_850_after_newfilters2       = np.array(FCFA_850_after_newfilters2)
    FCFA_err_850_after_newfilters2   = np.array(FCFA_err_850_after_newfilters2)
    
    plot_by_tauxairmass_2lists(FCFP_850_before_newfilters1,FCFP_err_850_before_newfilters1,FCFP_850_after_newfilters1,FCFP_err_850_after_newfilters1,FCFA_850_before_newfilters1,FCFA_err_850_before_newfilters1,FCFA_850_after_newfilters1,FCFA_err_850_after_newfilters1,FCFP_450_before_SMUtrouble1,FCFP_err_450_before_SMUtrouble1,FCFP_450_after_SMUtrouble1,FCFP_err_450_after_SMUtrouble1,FCFA_450_before_SMUtrouble1,FCFA_err_450_before_SMUtrouble1,FCFA_450_after_SMUtrouble1,FCFA_err_450_after_SMUtrouble1,FCFP_850_before_newfilters2,FCFP_err_850_before_newfilters2,FCFP_850_after_newfilters2,FCFP_err_850_after_newfilters2,FCFA_850_before_newfilters2,FCFA_err_850_before_newfilters2,FCFA_850_after_newfilters2,FCFA_err_850_after_newfilters2,FCFP_450_before_SMUtrouble2,FCFP_err_450_before_SMUtrouble2,FCFP_450_after_SMUtrouble2,FCFP_err_450_after_SMUtrouble2,FCFA_450_before_SMUtrouble2,FCFA_err_450_before_SMUtrouble2,FCFA_450_after_SMUtrouble2,FCFA_err_450_after_SMUtrouble2,Tau225TimesAirmass_bins,constrained='no')
    
    # Make a nice Latex Table to Summarise the Results
    
    print('\n\n')
    
    print('\\begin{deluxetable*}{c|c|c|c|c}')
    print('\\tablecaption{850 $\mu$m FCFs pre- and post- November 2016 Filter Change organised by $\\tau\\times\mathrm{AM}$.}')
    print('\label{tab:FCFtxAM850}')
    print('\\tablecolumns{5}')
    print('%\\tablenum{1}')
    print('\\tablewidth{0pt}')
    print('\\tablehead{')
    print('\colhead{$\\tau\\times\mathrm{AM}$} & ')
    print('\colhead{FCF$_{\mathrm{peak}}$}  & ')
    print('\colhead{FCF$_{\mathrm{arcsec}}$} & ')
    print('\colhead{FCF$_{\mathrm{peak}}$} & ') 
    print('\colhead{FCF$_{\mathrm{arcsec}}$} \\\\')
    print('&')
    print('\colhead{Before Filter Change} &') 
    print('\colhead{Before Filter Change} &') 
    print('\colhead{After Filter Change} &') 
    print('\colhead{After Filter Change}')
    print('}')
    print('\startdata')
    print('0.030$\\rightarrow$0.084 & '+str(round(FCFP_850_before_newfilters1[0],2))+' $\pm$ '+str(round(100*FCFP_err_850_before_newfilters1[0],1))+'\\% & '+str(round(FCFA_850_before_newfilters1[0],2))+' $\pm$ '+str(round(100*FCFA_err_850_before_newfilters1[0],1))+'\\% & '+str(round(FCFP_850_after_newfilters1[0],2))+' $\pm$ '+str(round(100*FCFP_err_850_after_newfilters1[0],1))+'\\% & '+str(round(FCFA_850_after_newfilters1[0],2))+' $\pm$ '+str(round(100*FCFA_err_850_after_newfilters1[0],1))+'\\% \\\\')
    print('0.084$\\rightarrow$0.138 & '+str(round(FCFP_850_before_newfilters1[1],2))+' $\pm$ '+str(round(100*FCFP_err_850_before_newfilters1[1],1))+'\\% & '+str(round(FCFA_850_before_newfilters1[1],2))+' $\pm$ '+str(round(100*FCFA_err_850_before_newfilters1[1],1))+'\\% & '+str(round(FCFP_850_after_newfilters1[1],2))+' $\pm$ '+str(round(100*FCFP_err_850_after_newfilters1[1],1))+'\\% & '+str(round(FCFA_850_after_newfilters1[1],2))+' $\pm$ '+str(round(100*FCFA_err_850_after_newfilters1[1],1))+'\\% \\\\')
    print('0.138$\\rightarrow$0.192 & '+str(round(FCFP_850_before_newfilters1[2],2))+' $\pm$ '+str(round(100*FCFP_err_850_before_newfilters1[2],1))+'\\% & '+str(round(FCFA_850_before_newfilters1[2],2))+' $\pm$ '+str(round(100*FCFA_err_850_before_newfilters1[2],1))+'\\% & '+str(round(FCFP_850_after_newfilters1[2],2))+' $\pm$ '+str(round(100*FCFP_err_850_after_newfilters1[2],1))+'\\% & '+str(round(FCFA_850_after_newfilters1[2],2))+' $\pm$ '+str(round(100*FCFA_err_850_after_newfilters1[2],1))+'\\% \\\\')
    print('0.192$\\rightarrow$0.246 & '+str(round(FCFP_850_before_newfilters1[3],2))+' $\pm$ '+str(round(100*FCFP_err_850_before_newfilters1[3],1))+'\\% & '+str(round(FCFA_850_before_newfilters1[3],2))+' $\pm$ '+str(round(100*FCFA_err_850_before_newfilters1[3],1))+'\\% & '+str(round(FCFP_850_after_newfilters1[3],2))+' $\pm$ '+str(round(100*FCFP_err_850_after_newfilters1[3],1))+'\\% & '+str(round(FCFA_850_after_newfilters1[3],2))+' $\pm$ '+str(round(100*FCFA_err_850_after_newfilters1[3],1))+'\\% \\\\')
    print('0.246$\\rightarrow$0.320 & '+str(round(FCFP_850_before_newfilters1[4],2))+' $\pm$ '+str(round(100*FCFP_err_850_before_newfilters1[4],1))+'\\% & '+str(round(FCFA_850_before_newfilters1[4],2))+' $\pm$ '+str(round(100*FCFA_err_850_before_newfilters1[4],1))+'\\% & '+str(round(FCFP_850_after_newfilters1[4],2))+' $\pm$ '+str(round(100*FCFP_err_850_after_newfilters1[4],1))+'\\% & '+str(round(FCFA_850_after_newfilters1[4],2))+' $\pm$ '+str(round(100*FCFA_err_850_after_newfilters1[4],1))+'\\%')
    print('\enddata')
    print('%\\tablenotetext{a}{This is a comment}')
    print('%\\tablenotetext{b}{This is a comment}')
    print('%\\tablenotetext{c}{This is a comment}')
    print('%\\tablenotetext{d}{This is a comment}')
    print('\\tablecomments{These are general comments (salutes).}')
    print('\end{deluxetable*}')
    
    print('\n\n')
    
    print('\n\n')
    
    print('\\begin{deluxetable*}{c|c|c|c|c}')
    print('\\tablecaption{450 $\mu$m FCFs pre- and post- 2018 SMU maintenance organised by $\\tau\\times\mathrm{AM}$.}')
    print('\label{tab:tab:FCFtxAM450}')
    print('\\tablecolumns{5}')
    print('%\\tablenum{1}')
    print('\\tablewidth{0pt}')
    print('\\tablehead{')
    print('\colhead{$\\tau\\times\mathrm{AM}$} & ')
    print('\colhead{FCF$_{\mathrm{peak}}$}  & ')
    print('\colhead{FCF$_{\mathrm{arcsec}}$} & ')
    print('\colhead{FCF$_{\mathrm{peak}}$} & ')
    print('\colhead{FCF$_{\mathrm{arcsec}}$} \\\\')
    print('&')
    print('\colhead{Before SMU Fix} &')
    print('\colhead{Before SMU Fix} &')
    print('\colhead{After SMU Fix} &')
    print('\colhead{After SMU Fix}')
    print('}')
    print('\startdata')
    print('0.030$\\rightarrow$0.084 & '+str(round(FCFP_450_before_SMUtrouble1[0],2))+' $\pm$ '+str(round(100*FCFP_err_450_before_SMUtrouble1[0],1))+'\\% & '+str(round(FCFA_450_before_SMUtrouble1[0],2))+' $\pm$ '+str(round(100*FCFA_err_450_before_SMUtrouble1[0],1))+'\\% & '+str(round(FCFP_450_after_SMUtrouble1[0],2))+' $\pm$ '+str(round(100*FCFP_err_450_after_SMUtrouble1[0],1))+'\\% & '+str(round(FCFA_450_after_SMUtrouble1[0],2))+' $\pm$ '+str(round(100*FCFA_err_450_after_SMUtrouble1[0],1))+'\\% \\\\')
    print('0.084$\\rightarrow$0.138 & '+str(round(FCFP_450_before_SMUtrouble1[1],2))+' $\pm$ '+str(round(100*FCFP_err_450_before_SMUtrouble1[1],1))+'\\% & '+str(round(FCFA_450_before_SMUtrouble1[1],2))+' $\pm$ '+str(round(100*FCFA_err_450_before_SMUtrouble1[1],1))+'\\% & '+str(round(FCFP_450_after_SMUtrouble1[1],2))+' $\pm$ '+str(round(100*FCFP_err_450_after_SMUtrouble1[1],1))+'\\% & '+str(round(FCFA_450_after_SMUtrouble1[1],2))+' $\pm$ '+str(round(100*FCFA_err_450_after_SMUtrouble1[1],1))+'\\% \\\\')
    print('0.138$\\rightarrow$0.192 & '+str(round(FCFP_450_before_SMUtrouble1[2],2))+' $\pm$ '+str(round(100*FCFP_err_450_before_SMUtrouble1[2],1))+'\\% & '+str(round(FCFA_450_before_SMUtrouble1[2],2))+' $\pm$ '+str(round(100*FCFA_err_450_before_SMUtrouble1[2],1))+'\\% & '+str(round(FCFP_450_after_SMUtrouble1[2],2))+' $\pm$ '+str(round(100*FCFP_err_450_after_SMUtrouble1[2],1))+'\\% & '+str(round(FCFA_450_after_SMUtrouble1[2],2))+' $\pm$ '+str(round(100*FCFA_err_450_after_SMUtrouble1[2],1))+'\\% \\\\')
    print('0.192$\\rightarrow$0.246 & '+str(round(FCFP_450_before_SMUtrouble1[3],2))+' $\pm$ '+str(round(100*FCFP_err_450_before_SMUtrouble1[3],1))+'\\% & '+str(round(FCFA_450_before_SMUtrouble1[3],2))+' $\pm$ '+str(round(100*FCFA_err_450_before_SMUtrouble1[3],1))+'\\% & '+str(round(FCFP_450_after_SMUtrouble1[3],2))+' $\pm$ '+str(round(100*FCFP_err_450_after_SMUtrouble1[3],1))+'\\% & '+str(round(FCFA_450_after_SMUtrouble1[3],2))+' $\pm$ '+str(round(100*FCFA_err_450_after_SMUtrouble1[3],1))+'\\% \\\\')
    print('0.246$\\rightarrow$0.320 & '+str(round(FCFP_450_before_SMUtrouble1[4],2))+' $\pm$ '+str(round(100*FCFP_err_450_before_SMUtrouble1[4],1))+'\\% & '+str(round(FCFA_450_before_SMUtrouble1[4],2))+' $\pm$ '+str(round(100*FCFA_err_450_before_SMUtrouble1[4],1))+'\\% & '+str(round(FCFP_450_after_SMUtrouble1[4],2))+' $\pm$ '+str(round(100*FCFP_err_450_after_SMUtrouble1[4],1))+'\\% & '+str(round(FCFA_450_after_SMUtrouble1[4],2))+' $\pm$ '+str(round(100*FCFA_err_450_after_SMUtrouble1[4],1))+'\\%')
    print('\enddata')
    print('%\\tablenotetext{a}{This is a comment}')
    print('%\\tablenotetext{b}{This is a comment}')
    print('%\\tablenotetext{c}{This is a comment}')
    print('%\\tablenotetext{d}{This is a comment}')
    print('\\tablecomments{These are general comments (salutes).}')
    print('\end{deluxetable*}')
    
    print('\n\n')