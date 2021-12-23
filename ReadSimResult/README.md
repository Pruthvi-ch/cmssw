# ReadSimResult
Read the CMS HGCAL simulation results
# Condor job submission
To submit the condor job,
1. One needs to modify the input and output paths in both the .py and .sh files inside the codor directory
2. Then first execute the .py file, this will create folder with name tmp*.
3. Next you have to go inside the tmp* and submit the jdl with condor_submit $filename.jdl
