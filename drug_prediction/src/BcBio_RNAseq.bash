#Have your fastq files in this location

/groups/"groupname"/"yourname"/"projectname"
	example for me this is /groups/springer/sarah/Isog

#Have your yaml file in this directory as well.
	example yaml_singleEnd.yaml

#Make a csv file describing your samples
	example sample_description.csv

	
#Run this line of code

bcbio_nextgen.py -w template yaml_singleEnd.yaml sample_description.csv *.fq (or *fastq depending on how files named)

#Running this sets up the data for bcbio creates project directory and some sub directories.

	NOTE - FILES CAN NOT be named _1 or _anynumber.fastq or will cause FATAL ERROR


#Move into the "work" directory
#Run the lsf from here 


bsub < submit_bcbio.lsf


vim submit_bcbio.lsf
##This is what submission looks like

#!/bin/bash

#BSUB -q priority
#BSUB -J bcbio_mov10
#BSUB -n 1
#BSUB -W 100:0
#BSUB -R "rusage[mem=10000]"
#BSUB -e mov10_project.err
#BSUB -o mov10_project.out

bcbio_nextgen.py ../config/mov10_project.yaml -n 32 -t ipython -s lsf -q parallel -r mincores=2 -r minconcores=2 '-rW=72:00' --retries 3 --timeout 380



## sometimes times out or memory issue - can just re-submit  with different settings - will look at where left off and start from there.
#alternate submission settings:


bcbio_nextgen.py ../config/mov10_project.yaml -n 32 -t ipython -s lsf -q parallel -r mincores=2 -r minconcores=2 '-rW=72:00' --retries 3 --timeout 380


##note hg19 - genome build 37

http://bcbio-nextgen.readthedocs.io/en/latest/index.html
https://github.com/chapmanb/bcbio-nextgen/tree/master/config/templates
