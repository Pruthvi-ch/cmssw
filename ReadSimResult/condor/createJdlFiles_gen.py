import itertools
import os
import sys
import subprocess
import time

#IMPORT MODULES FROM OTHER DIR

Geom_1 = ["Extended2026D86"]
#D86 = ["Extended2026D83"]


if not os.path.exists("tmpSub/log"):
    os.makedirs("tmpSub/log")
condorLogDir = "log"
tarFile = "tmpSub/generator.tar.gz"
if os.path.exists(tarFile):
    os.system("rm %s"%tarFile)
os.system("tar -zcvf %s ../../Configuration ../../ReadSimResult ../../Geometry --exclude condor"%tarFile)
os.system("cp rungen.sh tmpSub/")
common_command = \
'Universe   = vanilla\n\
should_transfer_files = YES\n\
when_to_transfer_output = ON_EXIT\n\
Transfer_Input_Files = generator.tar.gz, rungen.sh\n\
x509userproxy = $ENV(X509_USER_PROXY)\n\
use_x509userproxy = true\n\
RequestCpus = 4\n\
+BenchmarkJob = True\n\
#+JobFlavour = "testmatch"\n\
+MaxRuntime = 259200\n\
Output = %s/log_$(cluster)_$(process).stdout\n\
Error  = %s/log_$(cluster)_$(process).stderr\n\
Log    = %s/log_$(cluster)_$(process).condor\n\n'%(condorLogDir, condorLogDir, condorLogDir)

#----------------------------------------
#Create jdl files
#----------------------------------------
geom_par = 1
sampleList = eval("Geom_%i"%(geom_par))
jdlName = 'submitJobs_%s.jdl'%(geom_par)
jdlFile = open('tmpSub/%s'%jdlName,'w')
jdlFile.write('Executable =  rungen.sh \n')
jdlFile.write(common_command)
jdlFile.write("X=$(step)\n")
for sample in sampleList:
    condorOutDir1="/eos/user/p/psuryade/SimOut/DeltaPt"
    os.system("eos root://eosuser.cern.ch mkdir -p %s/%s"%(condorOutDir1, sample + "_r2"))
    #condorOutDir="/cms/store/user/idas/SimOut/DeltaPt"
    #os.system("xrdfs root://se01.indiacms.res.in/ mkdir -p %s/%s"%(condorOutDir, sample))
    run_command =  'Arguments  = %s $INT(X) \nQueue 10\n\n' %(sample)
    jdlFile.write(run_command)
    #print "condor_submit jdl/%s"%jdlFile
jdlFile.close() 


# subFile = open('tmpSub/condorSubmit.sh','w')
# for year in [2016,2017,2018]:
#     sampleList = eval("samples_%i"%year)
#     jdlName = 'submitJobs_%s.jdl'%(year)
#     jdlFile = open('tmpSub/%s'%jdlName,'w')
#     jdlFile.write('Executable =  rungen.sh \n')
#     jdlFile.write(common_command)
#     condorOutDir1="/eos/user/i/idas/SimOut/DeltaPt"
#     os.system("eos root://eosuser.cern.ch mkdir -p %s/%s"%(condorOutDir1, year))
#     condorOutDir="/cms/store/user/idas/SimOut/DeltaPt"
#     os.system("xrdfs root://se01.indiacms.res.in/ mkdir -p %s/%s"%(condorOutDir, year))
#     jdlFile.write("X=$(step)\n")
    
#     for sample in sampleList:
#         noflines = subprocess.Popen('wc -l ../input/eos/%i/%s_%i.txt | awk \'{print $1}\''%(year,sample,year),shell=True,stdout=subprocess.PIPE).communicate()[0].split('\n')[0]
#         nJob = int(noflines)
#         print "%s %s"%(sample,nJob)
#         if nJob==1:
#             run_command =  'Arguments  = %s %s input/eos/%i/%s_%i.txt 0 base \nQueue 1\n\n' %(year, sample, year, sample, year)
#         else:
#             run_command =  'Arguments  = %s %s input/eos/%i/%s_%i.txt $INT(X) base \nQueue %i\n\n' %(year, sample, year, sample, year, nJob)
#         jdlFile.write(run_command)
# 	#print "condor_submit jdl/%s"%jdlFile
#     subFile.write("condor_submit %s\n"%jdlName)
#     jdlFile.close() 
# subFile.close()
