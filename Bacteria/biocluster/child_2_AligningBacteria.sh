#!/bin/bash
#SBATCH -A pipeline
#SBATCH --job-name=child_2_AligningBacteria.sh # write a meaningful name
#SBATCH --time=7-00:00:00 # 100 hours as an example
#SBATCH --output=/home/bioinf/maja488d/logs/%x-%j.log # just write your zih name instead of sxxxxx, and keep the rest as they are 
#SBATCH --error=/home/bioinf/maja488d/logs/%x-%j.err  # just write your zih name instead of sxxxxx,  and keep the rest as they are 
#SBATCH --mem-per-cpu=5000  # set requested memory (in MB) 
#SBATCH --ntasks=1 # leave this for now
#SBATCH --nodes=1 # leave this for n
#SBATCH --cpus-per-task=1 # number of the request cpu if you have parallel processing
#SBATCH --mail-user madina.japakhova@mailbox.tu-dresden.de	### tell the batch system your email address to get updates about the jobs status
#SBATCH --mail-type=FAIL ### specify for what type of events you want to get a mail; valid options beside ALL are: BEGIN, END, FAIL, REQUEUE 


#### you must load the neede modules in advance
#### you must write module name exactly as it should be
# to see the available modules type this in terminal: module avail

## modules
module purge # to prevent conflicts 
module load apps/bwa/0.7.17
module load apps/samtools/1.12

analysis_dir=$1 
fastQ_files_folder=$2                                           
paired_end=$3 
filename_extention=$4 
samplename=$5 
readname=$6 

echo 'Started'

./scripts/2_AligningBacteria.sh $analysis_dir $fastQ_files_folder $paired_end $filename_extention $samplename $readname



echo "Finished at `date`"
echo "__DONE__"
