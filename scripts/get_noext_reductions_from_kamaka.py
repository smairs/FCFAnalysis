from jsa_proc.config import get_database
from jsa_proc.state import JSAProcState
from jsa_proc.admin.directories import get_output_dir
import os
from starlink import kappa

db = get_database()


x = db.find_jobs(state=JSAProcState.COMPLETE, task="cal-s2-noext", outputs="%.sdf")

dummy = 0

for eachjob in range(len(x)):

    jobid = x[eachjob].id

    #sdffiles = []
    #for eachfile in x[eachjob].outputs:
    #    sdffiles.append(str(eachfile))

    outputdir = get_output_dir(jobid)
    print outputdir

#    for eachfile in sdffiles:
#        os.system('cp '+outputdir+'/'+eachfile+' .')

