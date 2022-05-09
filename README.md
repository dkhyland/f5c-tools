# f5c-tools

This repository contains f5c forked from Hasindu Gamaarachchi's repo https://github.com/hasindu2008/f5c and a tool that may be used to tune the GPU parameters on a given machine.

## Usage
1. Clone the repository
2. Compile f5c 
3. Build f5c following the instructions from https://github.com/hasindu2008/f5c 
4. Download some test data sets and run minimap2 (to align the sequences to the genome), samtools sort (to create a .bam file) and then run samtools index on the generated .bam file. Make sure each data set is contained in a folder which has the same name as the .bam files, .fastq files and that all fast5 files are in a subfolder named fast5. The directory tree should have the following structure:
<pre><code>
f5c-tools/
|--data/\
   |--[reference genome].fa\
   |--[reference genome].fa.fai\
   |--[dataset_name]/\
      |--[fast5]/\
          |--\*.fast5\
      |--[dataset_name].bam
      |--[dataset_name].bam.bai
      |--[dataset_name].fastq
      |--[dataset_name].fastq.index./*
|--f5c/
   |--[other directories and files]
   |--scripts/
      |--[other scripts]
      |--param_test.sh
      |--remove_test.sh
      |--tune_parameters.sh
|--results/
   |--[dataset_name]/
   |--process_results.py
|--tools/
</code></pre>
5. Move the data sets to the data folder contained in this repository.
6. cd into the f5c folder.
7. Run scripts/tune_parameters.sh -d {dataset} -n {numruns} where {dataset} is the name of the folder containing the test data in the data folder, and {numruns} is the number of runs to perform for each set of parameters.
8. After the script has finished running (which can take hours-days), the recommended profile to use for your machine can be found in results/{dataset}/best.profile.
9. To run methylation-calling with this profile, run f5c with -x path/to/best.profile as a flag.

## Other information
- Results for individual tests are located in the results/{dataset} folder.
- More information can be found from the scripts by running them with the -h flag
