WORKSHOP = Anolis apletophallus 10x data Illumina (short reads) + Nanopore (long reads) 

 											
 										# DAY 1: CUTADAPT

RUNNER: Renata (pirani) + Kristin (charleskl)
FUNCTION: we are cleaning the DNA cutting the bad sequencing and checking for any contamination
PROGRAM WEBSITE: https://cutadapt.readthedocs.io/en/stable/guide.html)

JOB FILE: /scratch/genomics/piranir/Cutadapt/
		  /scratch/genomics/charleskl/CutAdapt

-> 1 Job: cutadapt_26.job

# /bin/sh                                                                                                                                   
# ----------------Parameters---------------------- #                                                                                        
#$ -S /bin/sh                                                                                                                               
#$ -q mThC.q                                                                                                                                
#$ -pe mthread 1 -l mres=6G,h_data=6G,h_vmem=6G                                                                                             
#$ -cwd                                                                                                                                     
#$ -N cutadapt                                                                                                                      
#$ -o cutadapt.log                                                                                                                  
#$ -m bea                                                                                                                                   
#$ -M piranir@si.edu                                                                                                                         
#                                                                                                                                           
# ----------------Modules------------------------- #                                                                                        
module load bioinformatics/cutadapt/2.4                                                                                                     
#                                                                                                                                           
# ----------------Your Commands------------------- #                                                                                        
#                                                                                                                                           
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                               
echo + NSLOTS = $NSLOTS                                                                                                                     
#                                                                                                                                           
#                                                                                                                                           
cutadapt -u 26 -o E28_26t_val_1.fq.gz /scratch/stri_ap/ariasc_data/anolis_10x/E28_MPS12345004_G06_9489_S3_L004_R1_001_val_1.fq.gz 
#                                                                                                                                           
echo = `date` job $JOB_NAME done


RESULTS
=== Summary ===                                                                                                                   
                                                                                                                                  
Total reads processed:             358,209,584                                                                                    
Reads written (passing filters):   358,209,584 (100.0%)                                                                           
                                                                                                                                  
                                                                                                                                  
Total basepairs processed: 52,891,743,310 bp                                                                                      
Total written (filtered):  43,579,357,893 bp (82.4%)
 											
 											
 											
 											DAY 2: Jellyfish
 											

RUNNER: Renata (piranir) + Kristin (charleskl)
FUNCTION: Run Genomescope, first you need to generate a Jellyfish histogram
		  Genomescope can be used to estimate genome size from short read data
PROGRAM WEBSITE: http://qb.cshl.edu/genomescope/

JOB FILE: /scratch/genomics/piranir/Jellyfish 
		  /scratch/genomics/charleskl/jellyfish

-> 1 job file: jellyfish.job 

# /bin/sh                                                                                                                                                                       
# ----------------Parameters---------------------- #                                                                                                                            
#$ -S /bin/sh                                                                                                                                                                   
#$ -q lThM.q                                                                                                                                                                    
#$ -pe mthread 30 -l mres=300G,h_data=10G,h_vmem=10G,himem                                                                                                                       
#$ -cwd                                                                                                                                                                         
#$ -j y                                                                                                                                                                         
#$ -N jellyfish1                                                                                                                                                                
#$ -o jellyfish1.log                                                                                                                                                            
#$ -m bea                                                                                                                                                                       
#$ -M piranir@si.edu                                                                                                                                                            
#                                                                                                                                                                               
# ----------------Modules------------------------- #                                                                                                                            
module load bioinformatics/jellyfish/2.3.0                                                                                                                                      
#                                                                                                                                                                               
# ----------------Your Commands------------------- #                                                                                                                            
#                                                                                                                                                                               
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                                   
echo + NSLOTS = $NSLOTS                                                                                                                                                         
#                                                                                                                                                                               
#                                                                                                                                                                               
gzip -dc E28_26t_val_1.fq.gz E28_10t_R2_001_val_2.fq.gz | jellyfish count -C -m 21 -s 800000000 -t $NSLOTS -o reads.jf /dev/fd/0                                                
#                                                                                                                                                                               
echo = `date` job $JOB_NAME done



ERROR: changed the permits / typing "chmod 777" , also those files are zip files. So now I need to unzip first and than to run the jellyfish code. 
 
 gzip -dc E28_26t_val_1.fq.gz E28_10t_R2_001_val_2.fq.gz = is to unzip the files. 
 
TIME: 10h to complete



------------------------------
-> 2 job: hist_jellyfish.job

# /bin/sh                                                                                                                    
# ----------------Parameters---------------------- #                                                                         
#$ -S /bin/sh                                                                                                                
#$ -pe mthread 20 -l mres=40G,h_data=2G,h_vmem=2G                                                                            
#$ -q mThC.q                                                                                                                 
#$ -cwd                                                                                                                      
#$ -j y                                                                                                                      
#$ -N jellyfish_histo_ill                                                                                                    
#$ -o jellyfish_histo_ill.log                                                                                                
#$ -m bea                                                                                                                    
#$ -M piranir@si.edu                                                                                                          
#                                                                                                                            
# ----------------Modules------------------------- #                                                                         
module load bioinformatics/jellyfish                                                                                                                                                                                                                                                                                                                                                                                                                               

# ----------------Your Commands------------------- #                                                                         
#                                                                                                                            
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                
echo + NSLOTS = $NSLOTS                                                                                                      
#                                                                                                                            
jellyfish histo -t 20 reads_EF_illumina.jf > reads_EF_ill.histo                                                             
#                                                                                                                            
echo = `date` job $JOB_NAME done   



TIME: 6 min to complete
 
NEXT STEPS: download from the hydra to your computer do:
- qrsh in your folder
- now navegate to the folder where the file is from here.
- type: module load bioinformatics/ffsend
- type: ffsend upload reads.histo and go to the website
- go to the genome website and upload your file (reads.histo)




											DAY 3: Wtdbg2/Redbeans


RUNNER: Renata (piranir) + Kristin (charleskl)
FUNCTION: Wtdbg2 is a sequence assembler for long noisy reads produced by either PacBio or Oxford Nanopore Technologies.
		 The program uses two steps. The first step is the assembler and the second step is the consenser. 
PROGRAM WEBSITE: https://github.com/ruanjue/wtdbg2/blob/master/README-ori.md
														
JOB FILE: /scratch/genomics/piranir/Wtdbg2 
		  /scratch/genomics/charleskl/redbean

	- To check the parameters go to "wtdbg2 --help" 
	- Renata runs with the flag -AS 1 
	- Kristin runs with the flag -A -S 2  


-> 1 Job: readbeans.job

# /bin/sh                                                                                                                                                               
#----------------Parameters---------------------- #                                                                                                                     
#$ -S /bin/sh                                                                                                                                                           
#$ -pe mthread 30                                                                                                                                                       
#$ -q mThM.q                                                                                                                                                            
#$ -l mres=300G,h_data=10G,h_vmem=10G,himem                                                                                                                             
#$ -cwd                                                                                                                                                                 
#$ -j y                                                                                                                                                                 
#$ -N readbeans                                                                                                                                                         
#$ -o readbeans.log                                                                                                                                                     
#$ -m bea                                                                                                                                                               
#$ -M piranir@si.edu                                                                                                                                                    
#                                                                                                                                                                       
# ----------------Modules------------------------- #                                                                                                                    
module load bioinformatics/wtdbg2                                                                                                                                       
#                                                                                                                                                                       
# ----------------Your Commands------------------- #                                                                                                                    
#                                                                                                                                                                       
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                           
echo + NSLOTS = $NSLOTS                                                                                                                                                 
#                                                                                                                                                                       
wtdbg2 -g 1.8g -t $NSLOTS -p 19 -A -S 1 -k 0 -s 0.05 -edge-min 2 --rescue-low-cov-edges -L 1000 
-i /scratch/stri_ap/ariasc_data/anolis_nanopore/anolis_28_filt.fa.gz -fo anolis_28_S1_assemble_wt
#                                                                                                                                                                       
echo = `date` job $JOB_NAME done

  

-> 2 Job: readbeans_2con.job (only Kristin run it, better results than Renata on the previous job)
						folder: /scratch/genomics/charleskl/jobs


# /bin/sh                                                                                                                                                               
#----------------Parameters---------------------- #                                                                                                                     
#$ -S /bin/sh                                                                                                                                                           
#$ -pe mthread 30                                                                                                                                                       
#$ -q mThM.q                                                                                                                                                            
#$ -l mres=300G,h_data=10G,h_vmem=10G,himem                                                                                                                             
#$ -cwd                                                                                                                                                                 
#$ -j y                                                                                                                                                                 
#$ -N readbeans_2con                                                                                                                                                         
#$ -o readbeans_2con.log                                                                                                                                                     
#$ -m bea                                                                                                                                                               
#$ -M charleskl@si.edu                                                                                                                                                    
#                                                                                                                                                                       
# ----------------Modules------------------------- #                                                                                                                    
module load bioinformatics/wtdbg2                                                                                                                                       
#                                                                                                                                                                       
# ----------------Your Commands------------------- #                                                                                                                    
#                                                                                                                                                                       
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                           
echo + NSLOTS = $NSLOTS                                                                                                                                                 
#                                                                                                                                                                       
wtpoa-cns -t $NSLOTS -i anolis_28_S2_assembly.ctg.lay.gz -fo anolis_S2_genomedraft_raw.fa
#                                                                                                                                                                       
echo = `date` job $JOB_NAME done


	
	P.S. to check queues = Go to confluence.si.edu and to "Submitting Jobs" and "available queues". 

	P.S. Renata jobs didn't run because of memory. We will use the assembly of Kristin 

	END: to check how is your assembly data type: 
	module load bioinformatics/assembly_stats
	assembly_stats anolis_S2_genomedraft_raw.fa
	

RESULTS
[piranir@hydra-login02 redbean]$ assembly_stats anolis_S2_genomedraft_raw.fa                                                                 
                                                                                                                                             
{                                                                                                                                            
  "Contig Stats": {                                                                                                                          
    "L10": 146,                                                                                                                              
    "L30": 696,                                                                                                                              
                                                                                                                                             
    "L40": 1095,                                                                                                                             
    "L50": 1597,                                                                                                                             
    "N10": 1252139,                                                                                                                          
    "N30": 676598,                                                                                                                           
                                                                                                                                             
    "N40": 542964,                                                                                                                           
    "N50": 428046,                                                                                                                           
    "gc_content": 43.53442576722115,                                                                                                         
    "mean": 119098.24670711854,                                                                                                              
                                                                                                                                             
    "median": 19088.0,                                                                                                                       
    "sequence_count": 20271,                                                                                                                 
    "shortest": 952,                                                                                                                         
  },"total_bps": 2414240559                                                                                                                  
                                                                                                                                             
  "Scaffold Stats": {                                                                                                                        
    "L10": 146,                                                                                                                              
    "L20": 380,                                                                                                                              
    "L40": 1095,                                                                                                                             
                                                                                                                                             
    "L50": 1597,                                                                                                                             
    "N10": 1252139,                                                                                                                          
    "N20": 874008,                                                                                                                           
    "N40": 542964,                                                                                                                           
                                                                                                                                             
    "N50": 428046,                                                                                                                           
    "gc_content": 43.53442576722115,                                                                                                         
    "longest": 4452698,                                                                                                                      
    "median": 19088.0,70711854,                                                                                                              
                                                                                                                                             
    "sequence_count": 20271,                                                                                                                 
    "shortest": 952,                                                                                                                         
    "total_bps": 2414240559                                                                                                                  
} } 



												
												DAY 4: Minimap/bwa



RUNNER: Kristin (charleskl)
FUNCTION: Polish the Genome 
		  minimap2:  is a versatile sequence alignment program that alings DNA or mRNA sequences against a large reference database. 
PROGRAM WEBSITE: https://github.com/lh3/minimap2
				 https://github.com/isovic/racon
				 http://www.htslib.org/doc/samtools-view.html
				 https://github.com/lh3/bwa
											
JOB FILE: /scratch/genomics/charleskl/jobs

	- We will used long reads data and short reads data to polish your draft genome by using these scripts:
	- For this step we are going to run 4 jobs. The first 2 jobs are for the long reads. 
	- For long reads, you can used minimap2 to map your long reads to your draft genome, after 
		that you will need to sort and filtered by reads that have the higher match to your draft genome, 
		finally you will run de concensus command.


-> 1 job: polish_anolis_S2.job

# /bin/sh                                                                                                                                        
# ----------------Parameters---------------------- #                                                                                             
#$ -S /bin/sh                                                                                                                     
#$ -pe mthread 40                                                                                                                 
#$ -q sThM.q                                                                                                                      
#$ -l mres=240G,h_data=6G,h_vmem=6G,himem                                                                                         
#$ -cwd                                                                                                                           
#$ -j y                                                                                                                           
#$ -N polish_anolis_S2                                                                                                            
#$ -o polish_anolis_S2.log                                                                                                        
#$ -m bea                                                                                                                         
#$ -M charleskl@si.edu                                                                                                                               
#  
# ----------------Modules------------------------- #                                                                                                                                                  
source /home/ariasc/.bashrc                                                                                                                
conda activate minimap2
module load bioinformatics/samtools                                                                                            
#                                                                                                                                                
# ----------------Your Commands------------------- #                                                                                             
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                    
echo + NSLOTS = $NSLOTS                                                                                                                          
#                                                                                                                                                  
minimap2 -t40 -ax map-ont /scratch/genomics/charleskl/redbean/anolis_S2_genomedraft_raw.fa 
/scratch/stri_ap/ariasc_data/anolis_nanopore/anolis_28_filt.fa.gz | samtools sort -@4 >anolis_S2_dbg.bam
#                                                     
echo = `date` job $JOB_NAME done 



	P.S. wait for the results of the Job 1 for you to start the job 2. 

	

-> 2 job: polish_anolis_S2_wtpoa.job
	polishing with long reads, filetering mapping reads, and consensus with wtpoa-cns module in redbean


# ----------------Parameters---------------------- #                                                                     
#$ -S /bin/sh                                                                                                            
#$ -pe mthread 40                                                                                                         
#$ -q sThM.q                                                                                                           
#$ -l mres=240G,h_data=6G,h_vmem=6G,himem                                                                                 
#$ -cwd                                                                                                                  
#$ -j y                                                                                                                   
#$ -N polish_anolis_S2_wtpoa                                                                                                             
#$ -o polish_anolis_S2_wtpoa.log                                                                                                         
#$ -m bea                                                                                                                
#$ -M charleskl@si.edu                                                                                                                               
#  
# ----------------Modules------------------------- #                                                                                                                                                  
module load bioinformatics/samtools                                                                                           
#                                                                                                                                                
# ----------------Your Commands------------------- #                                                                                             
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                    
echo + NSLOTS = $NSLOTS                                                                                                                          
#                                                                                                                                                  
samtools view -F0x900 /scratch/genomics/charleskl/jobs/anolis_S2_dbg.bam | /home/ariasc/programs/wtdbg2/./wtpoa-cns -t 40 -d /scra
tch/genomics/charleskl/redbean/anolis_S2_genomedraft_raw.fa -i - -fo Anolis_dbg_cns.fa
#                                                     
echo = `date` job $JOB_NAME done 

	
	
	SHORT READS
		-Now for short reads, you can used bwa to map your short reads to your draft genome, after 
	  	that you will sort the mapped file and finanlly run the concesus command.
		- We are going to use the BWA program. BWA is a software package for mapping DNA sequences against a large reference genome.
		- BWA first needs to construct the FM-index for the reference genome (the index command)


-> 3 job: bwa_index.job

####Indexing draft genome in bwa

# /bin/sh                                                                                                                
# ----------------Parameters---------------------- #                                                                     
#$ -S /bin/sh                                                                                                            
#$ -q sThM.q                                                                                                             
#$ -l mres=6G,h_data=6G,h_vmem=6G,himem                                                                                        
#$ -cwd                                                                                                                  
#$ -j y                                                                                                                   
#$ -N anolis_bwa_index                                                                                                             
#$ -o anolis_bwa_index.log                                                                                                          
#$ -m bea                                                                                                                
#$ -M charleskl@si.edu                                                                                                                               
# 
# ----------------Modules------------------------- #                                                                                                                                              
module load bioinformatics/bwa/0.7.17                                                                                             
#                                                                                                                                               
# ----------------Your Commands------------------- #                                                                                             
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                    
echo + NSLOTS = $NSLOTS                                                                                                                          
#                                                                                                                                               
bwa index /scratch/genomics/charleskl/jobs/Anolis_dbg_cns.fa  
#                                                    
echo = `date` job $JOB_NAME done 



-> 4 job: polish_anolis_small.job

####bwa polishing with wtpoa-cns module from redbean


# /bin/sh                                                                                                                
# ----------------Parameters---------------------- #                                                                     
#$ -S /bin/sh                                                                                                            
#$ -pe mthread 40                                                                                                        
#$ -q sThM.q                                                                                                            
#$ -l mres=240G,h_data=6G,h_vmem=6G,himem                                                                                
#$ -cwd                                                                                                                   
#$ -j y                                                                                                                  
#$ -N polish_anolis_small                                                                                                 
#$ -o polish_anolis_small.log                                                                                           
#$ -m bea                                                                                                                 
#$ -M charleskl@si.edu                                                                                                                             
# 
# ----------------Modules------------------------- #                                                                                                                                              
module load bioinformatics/bwa/0.7.17                                                                                     
module load bioinformatics/samtools                                                                                             
#                                                                                                                                                
# ----------------Your Commands------------------- #                                                                                             
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                    
echo + NSLOTS = $NSLOTS                                                                                                                          
#                                                                                                                                               
bwa mem -t 40 /scratch/genomics/charleskl/jobs/Anolis_dbg_cns.fa /scratch/genomics/charleskl/E28_26t_val_1.fq.gz /scratch/genomics
/charleskl/E28_10t_R2_001_val_2.fq.gz | samtools sort -O SAM | /home/ariasc/programs/wtdbg2/./wtpoa-cns -t 40 -x sam-sr -d /scratc
h/genomics/charleskl/jobs/Anolis_dbg_cns.fa -i - -fo anolis_dbg.srp.fa
#                                                    
echo = `date` job $JOB_NAME done 


	END: to check how is your assembly data type: 
		module load bioinformatics/assembly_stats
		assembly_stats anolis_dbg.srp.fa
	

RESULTS

	assembly_stats anolis_dbg.srp.fa                                                                    
{                                                                                                                                 
  "Contig Stats": {                                                                                                               
    "L10": 144,                                                                                                                   
    "L20": 377,                                                                                                                   
    "L30": 689,                                                                                                                   
    "L40": 1084,                                                                                                                  
    "L50": 1581,                                                                                                                  
    "N10": 1234042,                                                                                                               
    "N20": 861625,                                                                                                                
    "N30": 666723,                                                                                                                
    "N40": 535005,                                                                                                                
    "N50": 421509,                                                                                                                
    "gc_content": 43.707308313328454,                                                                                             
    "longest": 4407091,                                                                                                           
    "mean": 116353.47664150757,                                                                                                   
    "median": 18271.0,                                                                                                            
    "sequence_count": 20271,                                                                                                      
    "shortest": 249,                                                                                                              
    "total_bps": 2358601325                                                                                                       
  },                                                                                                                              
  "Scaffold Stats": {                                                                                                             
    "L10": 144,                                                                                                                   
    "L20": 377,                                                                                                                   
    "L30": 689,                                                                                                                   
    "L40": 1084,                                                                                                                  
    "L50": 1581,                                                                                                                  
    "N10": 1234042,                                                                                                               
    "N20": 861625,                                                                                                                
    "N30": 666723,                                                                                                                
    "N40": 535005,                                                                                                                
    "N50": 421509,                                                                                                                
    "gc_content": 43.707308313328454,                                                                                             
    "longest": 4407091,                                                                                                           
    "mean": 116353.47664150757,                                                                                                   
    "median": 18271.0,                                                                                                            
    "sequence_count": 20271,                                                                                                      
    "shortest": 249,                                                                                                              
    "total_bps": 2358601325                                                                                                       
  }                                                                                                                               
} 
							
														DAY 5: Scaff10x					
											


RUNNER: Renata (piranir) + Kristin (charleskl)
FUNCTION: Pipeline for scaffolding and breaking a genome assembly using 10x genomics linked-reads
		  make sure that the version is gcc CC= /software/gcc-4.9.2/bin/gcc in the makefile 

PROGRAM WEBSITE: https://github.com/wtsi-hpag/Scaff10X.git

INSTALL: $ git clone  https://github.com/wtsi-hpag/Scaff10X.git 
		 $ cd Scaff10X
		 $ ./install.sh			 

JOB FILES: /scratch/genomics/charleskl/jobs/Anolis_dbg_cns.fa
		   /scratch/stri_ap/ariasc_data/anolis_10x


#input2.dat									

q1=/scratch/stri_ap/ariasc_data/anolis_10x/E28_MPS12345004_G06_9489_S3_L004_R1_001_val_1.fq.gz                          
q2=/scratch/stri_ap/ariasc_data/anolis_10x/E28_MPS12345004_G06_9489_S3_L004_R2_001_val_2.fq.gz

P.S. Kristin will run with -block 100000
	 Renata will run with -block 15000


-> 1 job: /scratch/genomics/piranir/Scaff10/scaff10xAnolis.job
		  /scratch/genomics/charleskl/Scaff10x/scaff10x.job     

####scaffolding
 
# /bin/sh                                                                                                                
# ----------------Parameters---------------------- #                                                                     
#$ -S /bin/sh                                                                                                            
#$ -pe mthread 30                                                                                                        
#$ -q mThM.q                                                                                                              
#$ -l mres=240G,h_data=6G,h_vmem=6G,himem                                                                                      
#$ -cwd                                                                                                                   
#$ -j y                                                                                                                  
#$ -N scaff10x_Anolis_20k                                                                                                          
#$ -o scaff10x_Anolis_20k.log                                                                                                     
#$ -m bea                                                                                                                 
#$ -M charleskl@si.edu                                                                                                             
#                                                                                                                             
# ----------------Modules------------------------- #                                                                          
#                                                                                                                             
# ----------------Your Commands------------------- #                                                                          
#                                                                                                                             
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                 
echo + NSLOTS = $NSLOTS                                                                                                       
#                                                                                                                             
/home/ariasc/programs/Scaff10X/src/./scaff10x -nodes 30 -align bwa -score 10 -read-s1 5 -read-s2 5 -longread 1 -gap 100 -block 20000
 -data input2.dat Anolis_dbg_cns.fa Anolis_20k_scaf1.fa 
#                                                                                                                              
echo = `date` job $JOB_NAME done


- Check your assembly with assembly-stats
- module load bioinformatics/assembly_stats
- assembly_stats Anolis_20k_scaf1.fa


[piranir@hydra-login02 Racon_shortreads]$ assembly_stats Anolis_20k_scaf1.fa
{                                                                                                                                 
  "Contig Stats": {                                                                                                               
    "L10": 145,                                                                                                                   
    "L20": 379,                                                                                                                   
    "L30": 694,                                                                                                                   
    "L40": 1091,                                                                                                                  
    "L50": 1591,                                                                                                                  
    "N10": 1250105,                                                                                                               
    "N20": 873697,                                                                                                                
    "N30": 675791,                                                                                                                
    "N40": 541901,                                                                                                                
    "N50": 428229,                                                                                                                
    "gc_content": 43.51504132095798,                                                                                              
    "longest": 4432229,                                                                                                           
    "mean": 118502.74470919048,                                                                                                   
    "median": 18808.0,                                                                                                            
    "sequence_count": 20271,                                                                                                      
    "shortest": 249,                                                                                                              
    "total_bps": 2402169138                                                                                                       
  },                                                                                                                              
  "Scaffold Stats": {                                                                                                             
    "L10": 60,                                                                                                                    
    "L20": 157,                                                                                                                   
    "L30": 286,                                                                                                                   
    "L40": 451,                                                                                                                   
    "L50": 663,                                                                                                                   
    "N10": 2935956,                                                                                                               
    "N20": 2136871,                                                                                                               
    "N30": 1648216,                                                                                                               
    "N40": 1279402,                                                                                                               
    "N50": 997812,                                                                                                                
    "gc_content": 43.51504132095798,                                                                                              
    "longest": 7218738,                                                                                                           
    "mean": 141473.2091626428,                                                                                                    
    "median": 14259.0,                                                                                                            
    "sequence_count": 16982,                                                                                                      
    "shortest": 249,                                                                                                              
    "total_bps": 2402498038                                                                                                       
  }                                                                                                                               
}                           


	P.S. - Kristin also run got a better N50 = 997812 (file: Anolis_20k_scaf1.fa), so we will use her data to polish.
		 - Renata :  N50 = 974899 (file Anolis_dbg8_RMP_scaf1.fa)
	AFTER: polish with Racun for long read data, we use minimap

	- Python no adapters and no barcodes - short read - map with 	



														DAY 6: minimap2/racon long reads


RUNNER: Kristin (charleskl)
FUNCTION: Polish the Genome 
		  minimap2:  is a versatile sequence alignment program that alings DNA or mRNA sequences against a large reference database. 
		  Polishing with Racon -- Type of data: Nanopore or Pacbio

PROGRAM WEBSITE: https://github.com/lh3/minimap2
				
JOB FILE: /scratch/genomics/charleskl/Racon/Anolis_20k_scaf1.fa

P.S. Kristin data was the best so we are going to use this file: Anolis_20k_scaf1.fa



--> 1 job: anolis_minimap_racon.job 
	First map your clean_raw reads to your assembly with minimap.


# /bin/sh
# ----------------Parameters---------------------- #
#$  -S /bin/sh
#$ -pe mthread 32
#$ -q mThM.q
#$ -l mres=4G,h_data=4G,h_vmem=4G,himem
#$ -cwd
#$ -j y
#$ -N anolis_minimap_racon
#$ -o anolis_minimap_racon.log
#$ -m bea
#$ -M charleskl@si.edu
#
# ----------------Modules------------------------- #
source /home/ariasc/.bashrc                                                                                                                
conda activate minimap2      
# ----------------Your Commands------------------- #
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
echo + NSLOTS = $NSLOTS
#
minimap2 -ax map-ont -t40 /scratch/genomics/charleskl/Racon/Anolis_20k_scaf1.fa /scratch/stri_ap/ariasc_data/anolis_nanopore/anoli
s_28_filterd.fastq > anolis_mapped_data.sam
#
echo = `date` job $JOB_NAME done



--> 2 Job: racon.job


# /bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -pe mthread 30
#$ -q mThM.q
#$ -l mres=240G,h_data=8G,h_vmem=8G,himem
#$ -cwd
#$ -j y
#$ -N anolis_racon
#$ -o anolis_racon.log
#$ -m bea
#$ -M charleskl@si.edu
#
# ----------------Modules------------------------- #
source /home/ariasc/.bashrc 
conda activate racon
# ----------------Your Commands------------------- #
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
echo + NSLOTS = $NSLOTS
#
racon -u -t 30 /scratch/stri_ap/ariasc_data/anolis_nanopore/anolis_28_filterd.fastq anolis_mapped_data.sam /scratch/genomics/charl
eskl/Racon/Anolis_20k_scaf1.fa > anolisgenome_racon_polished.fasta
#
echo = `date` job $JOB_NAME done


RESULTS: 
[piranir@hydra-login02 Racon]$ assembly_stats anolisgenome_racon_polished.fasta                                                   
{                                                                                                                                 
  "Contig Stats": {                                                                                                               
    "L10": 93,                                                                                                                    
    "L20": 240,                                                                                                                   
    "L30": 437,                                                                                                                   
    "L40": 689,                                                                                                                   
    "L50": 1012,                                                                                                                  
    "N10": 1991459,                                                                                                               
    "N20": 1413682,                                                                                                               
    "N30": 1095008,                                                                                                               
    "N40": 851544,                                                                                                                
    "N50": 659956,                                                                                                                
    "gc_content": 43.61210181870258,                                                                                              
    "longest": 4504359,                                                                                                           
    "mean": 133416.50532732866,                                                                                                   
    "median": 15858.0,                                                                                                            
    "sequence_count": 18208,                                                                                                      
    "shortest": 249,                                                                                                              
    "total_bps": 2429247729                                                                                                       
  },                                                                                                                              
  "Scaffold Stats": {                                                                                                             
    "L10": 60,                                                                                                                    
    "L20": 158,                                                                                                                   
    "L30": 286,                                                                                                                   
    "L40": 451,                                                                                                                   
    "L50": 663,                                                                                                                   
    "N10": 2964747,                                                                                                               
    "N20": 2159052,                                                                                                               
    "N30": 1663944,                                                                                                               
    "N40": 1296520,                                                                                                               
    "N50": 1007391,                                                                                                               
    "gc_content": 43.61210181870258,                                                                                              
    "longest": 7302435,                                                                                                           
    "mean": 143055.20710163703,                                                                                                   
    "median": 14419.0,                                                                                                            
    "sequence_count": 16982,                                                                                                      
    "shortest": 249,                                                                                                              
    "total_bps": 2429363527                                                                                                       
  }                                                                                                                               
} 


ERROR: Kristin got warnings about 2 contigs that were possibly chimeric
		We will run Break10x to fix the error. 
		

#########Break10x###################

--> 3 job: /scratch/genomics/piranir/Scaff10/break10xAnolis.job
	User Time = 10:02:18:44

# /bin/sh                                                                                                                
# ----------------Parameters---------------------- #                                                                     
#$ -S /bin/sh                                                                                                            
#$ -pe mthread 30                                                                                                        
#$ -q mThM.q                                                                                                              
#$ -l mres=240G,h_data=8G,h_vmem=8G,himem                                                                                      
#$ -cwd                                                                                                                   
#$ -j y                                                                                                                  
#$ -N break10x                                                                                                          
#$ -o break10x.log                                                                                                     
#$ -m bea                                                                                                                 
#$ -M piranir@si.edu                                                                                                             
#                                                                                                                             
# ----------------Modules------------------------- #                                                                          
#                                                                                                                             
# ----------------Your Commands------------------- #                                                                          
#                                                                                                                             
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                 
echo + NSLOTS = $NSLOTS                                                                                                       
#                                                                                                                             
/home/ariasc/programs/Scaff10X/src/./break10x -nodes 30 -score 10 -gap 100 -data input2.dat 
anolisgenome_racon_polished.fasta scaffolds-break.fasta scaffolds-break.name
#                                                                                                                              
echo = `date` job $JOB_NAME done


RESULTS: 
- Check your assembly with assembly-stats
- module load bioinformatics/assembly_stats
[piranir@hydra-login02 Scaff10x]$ assembly_stats scaffolds-break.fasta
 "Contig Stats": {                                                                                               
    "L10": 99,                                                                                                    
    "L20": 254,                                                                                                   
    "L30": 459,                                                                                                   
    "L40": 721,                                                                                                   
    "L50": 1055,                                                                                                  
    "N10": 1887405,                                                                                               
    "N20": 1333134,                                                                                               
    "N30": 1049363,                                                                                               
    "N40": 821276,                                                                                                
    "N50": 639826,                                                                                                
    "gc_content": 43.61210181870258,                                                                              
    "longest": 4504359,                                                                                           
    "mean": 132369.6452157803,                                                                                    
    "median": 16033.5,                                                                                            
    "sequence_count": 18352,                                                                                      
    "shortest": 249,                                                                                              
    "total_bps": 2429247729                                                                                       
  },                                                                                                              
  "Scaffold Stats": {                                                                                             
    "L10": 66,                                                                                                    
    "L20": 171,                                                                                                   
    "L30": 308,                                                                                                   
    "L40": 485,                                                                                                   
    "L50": 712,                                                                                                   
    "N10": 2762282,                                                                                               
    "N20": 2022627,                                                                                               
    "N30": 1554218,                                                                                               
    "N40": 1218296,                                                                                               
    "N50": 944294,                                                                                                
    "gc_content": 43.61210181870258,                                                                              
    "longest": 7302435,                                                                                           
    "mean": 141761.24333313882,                                                                                   
    "median": 14608.0,                                                                                            
    "sequence_count": 17137,                                                                                      
    "shortest": 249,                                                                                              
    "total_bps": 2429362427                                                                                       
  }                                                                                                               
} 


														DAY 7: BUSCO3 from short reads before Polishing



RUNNER: Carlos (ariasc) + Renata (piranir)
FUNCTION: assess the completeness of genomes, gene sets, and transcriptomes, 
		  using their gene content as a complementary method to common technical metrics. 
		 
PROGRAM WEBSITE: https://busco.ezlab.org/busco_userguide.html

JOB FILE: /scratch/genomics/piranir/Busco/

BEFORE RUNNING: 
	- Running last version of Busco that I have installed
	- To run busco4 you need to first to cp config file from augustus. Please do this:
	- make new directory in your path 
	- cd /scratch/genomics/username/
	- mkdir busco
	- cd busco
	- and Assuming you are in the folder `/scratch/genomics/username/busco`:
	- cp -r /share/apps/bioinformatics/augustus/conda/3.3.2/config/ .

--> 1 Job: busco_before.job

# /bin/sh                                                                                                                       
# ----------------Parameters---------------------- #                                                                            
#$ -S /bin/sh                                                                                                                   
#$ -pe mthread 20 
#$ -l mres=100G,h_data=5G,h_vmem=5G,himem                                                                      
#$ -q  mThM.q                                                                                                                 
#$ -cwd                                                                                                                         
#$ -j y                                                                                                                         
#$ -N busco3_before                                                                                                                
#$ -o busco3_before.log                                                                                                            
#$ -m bea                                                                                                                       
#$ -M piranir@si.edu                                                                                                             
#                                                                                                                               
# ----------------Modules------------------------- #                                                                            
module load bioinformatics/busco/3.0.2                                                                                                           
#                                                                                                                               
# ----------------Your Commands------------------- #                                                                            
#                                                                                                                               
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                   
echo + NSLOTS = $NSLOTS
# 
export AUGUSTUS_CONFIG_PATH="/scratch/genomics/piranir/Busco/augustus/config"                                                              
run_BUSCO.py -m genome -i scaffolds-break.fasta -o Anolis_before -l tetrapoda_odb9 -c $NSLOTS 
#                                                                                                                           
echo = `date` job $JOB_NAME done




RESULTS:
User Time = 4:16:06:52
Summarized benchmarking in BUSCO notation for file scaffolds-break.fasta                                                                    
# BUSCO was run in mode: genome                                                                                                               
                                                                                                                                              
        C:67.5%[S:66.7%,D:0.8%],F:16.5%,M:16.0%,n:3950                                                                                        
                                                                                                                                              
        2664    Complete BUSCOs (C)                                                                                                           
        2633    Complete and single-copy BUSCOs (S)                                                                                           
        31      Complete and duplicated BUSCOs (D)                                                                                            
        652     Fragmented BUSCOs (F)                                                                                                         
        634     Missing BUSCOs (M)                                                                                                            
        3950    Total BUSCO groups searched 



														DAY 8: Minimap2/Racon short reads
														
														
RUNNER: Renata (piranir)
FUNCTION: polishing with racon Illumina short reads
 	  
PROGRAM WEBSITE: https://github.com/isovic/racon
				 https://github.com/lh3/minimap2

JOB FOLDER: /scratch/genomics/piranir/Racon_shortreads/				

FIRST: 
	- rename reads 1 a 2 to used them with racon.
	- for questions go to /scratch/stri_ap/ariasc_data/anolis_10x and the files are there.
	- for the python file use the path /scratch/genomics/ariasc/Efish/
	- copy those files to your folder and unzip them using - gunzip *.gz 

python fix_reads_racon.py E28_10t_R2_001_val_2.fq E28_26t_R1_001_val_1.fq > rename_x_racon.fastq

	- second, map reads to draft genome with minimap2


--> 1 Job: minimap.job

# /bin/sh                                                                                                                               
# ----------------Parameters---------------------- #                                                                                    
#$ -S /bin/sh                                                                                                                           
#$ -pe mthread 30                                                                                                                       
#$ -q sThM.q                                                                                                                            
#$ -l mres=150G,h_data=5G,h_vmem=5G,himem                                                                                               
#$ -cwd                                                                                                                                 
#$ -j y                                                                                                                                 
#$ -N minimap                                                                                                                           
#$ -o minimap.log                                                                                                                       
#$ -m bea                                                                                                                               
#$ -M piranir@si.edu                                                                                                                     
#                                                                                                                                       
# ----------------Modules------------------------- #                                                                                    
source /home/ariasc/.bashrc                                                                                                                        
conda activate minimap2                                                                                                                 
#                                                                                                                                       
# ----------------Your Commands------------------- #                                                                                    
#                                                                                                                                       
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                           
echo + NSLOTS = $NSLOTS                                                                                                                 
#                                                                                                                                       
minimap2 -ax sr -t $NSLOTS /scratch/genomics/piranir/Racon_shortreads/scaffolds-break.fasta rename_x_racon.fastq > Anolis_notC_scaff.sam 
#                                                                                                                                       
echo = `date` job $JOB_NAME done                                                                                                        



--> 2 job: racon_short.job
	run racon with the *.sam file and corrected names short reads files

	
# /bin/sh                                                                                                                              
# ----------------Parameters---------------------- #                                                                                   
#$ -S /bin/sh                                                                                                                           
#$ -pe mthread 30                                                                                                                      
#$ -q uTxlM.rq
#$ -l mres=900G,h_data=30G,h_vmem=30G,himem                                                                                            
#$ -cwd                                                                                                                                 
#$ -j y                                                                                                                                
#$ -N racon                                                                                                                             
#$ -o racon.log                                                                                                                        
#$ -m bea                                                                                                                               
#$ -M piranir@si.edu                                                                                                                    
#                                                                                                                                       
# ----------------Modules------------------------- #                                                                                   
source ~/.bashrc                                                                                                                       
conda activate racon                                                                                                                    
#                                                                                                                                      
# ----------------Your Commands------------------- #                                                                                    
#                                                                                                                                      
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                           
echo + NSLOTS = $NSLOTS                                                                                                                
#                                                                                                                                       
racon -t 30 rename_x_racon.fastq Anolis_notC_scaff.sam scaffolds-break.fasta > Anolis_racon_short.fasta
#                                                                                                                                      
echo = `date` job $JOB_NAME done 



RESULTS: 
User Time = 20:49:10
- Check your assembly with assembly-stats
- module load bioinformatics/assembly_stats

{piranir@hydra-login02 Racon_shortreads]$ assembly_stats Anolis_racon_short.fasta                                                             
  "Contig Stats": {                                                                                                                           
                                                                                                                                              
    "L10": 96,                                                                                                                                
    "L20": 251,                                                                                                                               
    "L40": 716,                                                                                                                               
    "L50": 1048,                                                                                                                              
                                                                                                                                              
    "N10": 1903798,                                                                                                                           
    "N20": 1335717,                                                                                                                           
    "N40": 825896,,                                                                                                                           
    "N50": 642306,                                                                                                                            
                                                                                                                                              
    "gc_content": 43.814300913156785,                                                                                                         
    "longest": 4513562,                                                                                                                       
    "median": 16518.0,4353537,                                                                                                                
    "sequence_count": 17967,                                                                                                                  
                                                                                                                                              
    "shortest": 1,                                                                                                                            
    "total_bps": 2428508377                                                                                                                   
  "Scaffold Stats": {                                                                                                                         
    "L10": 66,                                                                                                                                
                                                                                                                                              
    "L20": 170,                                                                                                                               
    "L30": 308,                                                                                                                               
    "L50": 710,                                                                                                                               
    "N10": 2767450,                                                                                                                           
                                                                                                                                              
    "N20": 2026324,                                                                                                                           
    "N30": 1554505,                                                                                                                           
    "N50": 947285,,                                                                                                                           
    "gc_content": 43.814300913156785,                                                                                                         
                                                                                                                                              
    "longest": 7315638,                                                                                                                       
    "mean": 145748.69081197862,                                                                                                               
    "sequence_count": 16663,                                                                                                                  
    "shortest": 391,                                                                                                                          
                                                                                                                                              
    "total_bps": 2428610435                                                                                                                   
} } 


														DAY 7: BUSCO 3 from short reads AFTER Polishing
														
	
RUNNER: Renata (piranir)
FUNCTION: assess the completeness of genomes, gene sets, and transcriptomes, 
		  using their gene content as a complementary method to common technical metrics. 
		 
PROGRAM WEBSITE: https://busco.ezlab.org/busco_userguide.html

JOB FILE: /scratch/genomics/piranir/Busco/

BEFORE RUNNING: 
	- Running last version of Busco that I have installed
	- To run busco4 you need to first to cp config file from augustus. Please do this:
	- make new directory in your path 
	- cd /scratch/genomics/username/
	- mkdir busco
	- cd busco
	- and Assuming you are in the folder `/scratch/genomics/username/busco`:
	- cp -r /share/apps/bioinformatics/augustus/conda/3.3.2/config/ .
 
 SENT FILES TO HYDRA
	- to send files to hydra go: https://confluence.si.edu/display/HPC/Using+rclone
	- or you can use wget "website link"



--> 1 Job: busco3_after.job

# /bin/sh                                                                                                                       
# ----------------Parameters---------------------- #                                                                            
#$ -S /bin/sh                                                                                                                   
#$ -pe mthread 20 
#$ -l mres=100G,h_data=5G,h_vmem=5G,himem                                                                      
#$ -q  mThM.q                                                                                                                 
#$ -cwd                                                                                                                         
#$ -j y                                                                                                                         
#$ -N busco3_after                                                                                                                
#$ -o busco3_after.log                                                                                                            
#$ -m bea                                                                                                                       
#$ -M piranir@si.edu                                                                                                             
#                                                                                                                               
# ----------------Modules------------------------- #                                                                            
module load bioinformatics/busco/3.0.2                                                                                                           
#                                                                                                                               
# ----------------Your Commands------------------- #                                                                            
#                                                                                                                               
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                     
echo + NSLOTS = $NSLOTS                                                                                           
#                                                                                                                 
export AUGUSTUS_CONFIG_PATH="/scratch/genomics/piranir/Busco/augustus/config"                                                              
run_BUSCO.py -m genome -i Anolis_racon_short.fasta -o Anolis_after -l tetrapoda_odb9 -c $NSLOTS 
#                                                                                                                 
echo = `date` job $JOB_NAME done														
														
														
RESULTS:
User Time = 5:12:39:42	

BUSCO was run in mode: genome                                                                                                               
                                                                                                                                              
        C:88.1%[S:87.1%,D:1.0%],F:7.0%,M:4.9%,n:3950                                                                                          
                                                                                                                                              
        3482    Complete BUSCOs (C)                                                                                                           
        3441    Complete and single-copy BUSCOs (S)                                                                                           
        41      Complete and duplicated BUSCOs (D)                                                                                            
        277     Fragmented BUSCOs (F)                                                                                                         
        191     Missing BUSCOs (M)                                                                                                            
        3950    Total BUSCO groups searched  												
														
		
		
NEXT STEPS:
	- After talking to Carlos, we decided to run polishing again for the Anolis_racon_short.fasta file		
	- So Carlos will run pilon and Renata the minimap and Racon again 
	- Kristin will download the old Anolis genome and compare to our results. 

FILE: Anolis_racon_short.fasta
												
												
												
															DAY 9: Minimap2/Racon to improve the BUSCO RESULTS
														
														
RUNNER: Renata (piranir)

JOB FOLDER: /scratch/genomics/piranir/Racon_shortreads/				

FILE: Anolis_racon_short.fasta


--> 1 Job: minimap2.job

# /bin/sh                                                                                                                               
# ----------------Parameters---------------------- #                                                                                    
#$ -S /bin/sh                                                                                                                           
#$ -pe mthread 30                                                                                                                       
#$ -q sThM.q                                                                                                                            
#$ -l mres=150G,h_data=5G,h_vmem=5G,himem                                                                                               
#$ -cwd                                                                                                                                 
#$ -j y                                                                                                                                 
#$ -N minimap2                                                                                                                           
#$ -o minimap2.log                                                                                                                       
#$ -m bea                                                                                                                               
#$ -M piranir@si.edu                                                                                                                     
#                                                                                                                                       
# ----------------Modules------------------------- #                                                                                    
source /home/ariasc/.bashrc                                                                                                                        
conda activate minimap2                                                                                                                 
#                                                                                                                                       
# ----------------Your Commands------------------- #                                                                                    
#                                                                                                                                       
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                           
echo + NSLOTS = $NSLOTS                                                                                                                 
#                                                                                                                                       
minimap2 -ax sr -t $NSLOTS /scratch/genomics/piranir/Racon_shortreads/Anolis_racon_short.fasta rename_x_racon.fastq > Anolis_notC_scaff2.sam 
#                                                                                                                                       
echo = `date` job $JOB_NAME done                                                                                                        



--> 2 job: racon2.job
	run racon with the *.sam file and corrected names short reads files

	
# /bin/sh                                                                                                                              
# ----------------Parameters---------------------- #                                                                                   
#$ -S /bin/sh                                                                                                                           
#$ -pe mthread 30                                                                                                                      
#$ -q uTxlM.rq
#$ -l mres=900G,h_data=30G,h_vmem=30G,himem                                                                                            
#$ -cwd                                                                                                                                 
#$ -j y                                                                                                                                
#$ -N racon2                                                                                                                             
#$ -o racon2.log                                                                                                                        
#$ -m bea                                                                                                                               
#$ -M piranir@si.edu                                                                                                                    
#                                                                                                                                       
# ----------------Modules------------------------- #                                                                                   
source ~/.bashrc                                                                                                                       
conda activate racon                                                                                                                    
#                                                                                                                                      
# ----------------Your Commands------------------- #                                                                                    
#                                                                                                                                      
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                           
echo + NSLOTS = $NSLOTS                                                                                                                
#                                                                                                                                       
racon -t 30 rename_x_racon.fastq Anolis_notC_scaff2.sam Anolis_racon_short.fasta > Anolis_after_polishing2.fasta
#                                                                                                                                      
echo = `date` job $JOB_NAME done 



RESULTS: 

- Check your assembly with assembly-stats
- module load bioinformatics/assembly_stats											

[piranir@hydra-login02 Racon_shortreads]$ assembly_stats Anolis_after_polishing2.fasta                            
  "Contig Stats": {                                                                                               
                                                                                                                  
    "L10": 95,                                                                                                    
    "L20": 246,                                                                                                   
    "L30": 446,                                                                                                   
    "L50": 1029,                                                                                                  
                                                                                                             
    "N10": 1954367,                                                                                               
    "N20": 1377123,                                                                                               
    "N30": 1062227,                                                                                               
    "N50": 657933,                                                                                                                                                          
    "gc_content": 43.815054602762736,                                                                             
    "longest": 4506118,                                                                                           
    "mean": 135760.43425915556,                                                                                   
    "sequence_count": 17858,                                                                                                                    
    "shortest": 1,                                                                                                
    "total_bps": 2424409835                                                                                       
  },                                                                                                              
    "L10": 66,ats": {                                                                                                                          
    "L20": 170,                                                                                                   
    "L30": 307,                                                                                                   
    "L40": 484,                                                                                                   
    "N10": 2763557,                                                                                                                          
    "N20": 2024021,                                                                                               
    "N30": 1556319,                                                                                               
    "N40": 1219413,                                                                                               
    "gc_content": 43.815054602762736,                                                                                                                                        
    "longest": 7305303,                                                                                           
    "mean": 146098.48418198252,                                                                                   
    "sequence_count": 16595,                                                                                      
    "shortest": 391,                                                                                                                                                  
    "total_bps": 2424504345                                                                                       
  }                  												
	
	
	
	
	
	
												
														DAY 10: BUSCO 3 AFTER2
														
	
RUNNER: Renata (piranir)
FUNCTION: assess the completeness of genomes, gene sets, and transcriptomes, 
		  using their gene content as a complementary method to common technical metrics. 
		 
PROGRAM WEBSITE: https://busco.ezlab.org/busco_userguide.html

JOB FILE: /scratch/genomics/piranir/Busco/Anolis_after_polishing2.fasta


--> 1 Job: busco3_after2.job

# /bin/sh                                                                                                                       
# ----------------Parameters---------------------- #                                                                            
#$ -S /bin/sh                                                                                                                   
#$ -pe mthread 20 
#$ -l mres=100G,h_data=5G,h_vmem=5G,himem                                                                      
#$ -q  mThM.q                                                                                                                 
#$ -cwd                                                                                                                         
#$ -j y                                                                                                                         
#$ -N busco3_after2                                                                                                                
#$ -o busco3_after2.log                                                                                                            
#$ -m bea                                                                                                                       
#$ -M piranir@si.edu                                                                                                             
#                                                                                                                               
# ----------------Modules------------------------- #                                                                            
module load bioinformatics/busco/3.0.2                                                                                                           
#                                                                                                                               
# ----------------Your Commands------------------- #                                                                            
#                                                                                                                               
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                     
echo + NSLOTS = $NSLOTS                                                                                           
#                                                                                                                 
export AUGUSTUS_CONFIG_PATH="/scratch/genomics/piranir/Busco/augustus/config"                                                              
run_BUSCO.py -m genome -i Anolis_after_polishing2.fasta -o Anolis_after2 -l tetrapoda_odb9 -c $NSLOTS 
#                                                                                                                 
echo = `date` job $JOB_NAME done														
														
														
RESULTS:										
	C:88.1%[S:87.1%,D:1.0%],F:7.1%,M:4.8%,n:3950                                                              
                                                                                                                  
        3439    Complete and single-copy BUSCOs (S)                                                               
        40      Complete and duplicated BUSCOs (D)                                                                
                                                                                                                  
        279     Fragmented BUSCOs (F)                                                                             
        192     Missing BUSCOs (M)											
														


######## WE WILL WAIT FOR THE PILON RESULTS



															
															DAY 11: PILON
															
RUNNER: Carlos (ariasc)
FUNCTION: 
		   
		 
PROGRAM WEBSITE: 

JOB FILE: /scratch/genomics/ariasc/anolis/anolis_d1.fasta
														
														
#$ -S /bin/sh                                                                  
# ----------------Parameters---------------------- #                           
#$ -q uTxlM.rq                                                                 
#$ -pe mthread 24 -l mres=1100G,h_data=45G,h_vmem=45G,himem                    
#$ -cwd                                                                        
#$ -j y                                                                        
#$ -N ariasc_pilon_d1                                                          
#$ -o Ariasc_pilon_d1.log                                                      
#$ -m bea                                                                      
#$ -M ariasc@si.edu                                                            
#                                                                              
# ----------------Modules------------------------- #                           
module load bioinformatics/pilon/1.23                                          
#                                                                              
# ----------------Your Commands------------------- #                           
#                                                                              
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME  
echo + NSLOTS = $NSLOTS                                                        
#                                                                              
#                                                                              
PILON_HEAP_SIZE=1000g                                                                                                             
runpilon --genome /scratch/genomics/ariasc/anolis/Anolis_racon_short.fasta --bam /scratch/genomics/ariasc/anolis/Anolis_racon_short_sorted.bam  
--output  anolis_d1 --outdir /scratch/genomics/ariasc/anolis/ --fix bases --changes --vcf --threads $NSLOTS > /scratch/genomics/ariasc/anolis/anolis_d1_out.log                    
#
echo = `date` job $JOB_NAME done 												
														

OUTPUT FILE: anolis_d1.fasta														
															
Assembly stats before kraken
{                                                                                                                                
  "Contig Stats": {                                                                                                               
    "L10": 96,                                                                                                                   
    "L20": 251,                                                                                                                   
    "L30": 455,                                                                                                                  
    "L40": 716,                                                                                                                   
    "L50": 1048,                                                                                                                 
    "N10": 1903798,                                                                                                               
    "N20": 1335717,                                                                                                              
    "N30": 1053586,                                                                                                               
    "N40": 825896,                                                                                                               
    "N50": 642306,                                                                                                                
    "gc_content": 43.814300913156785,                                                                                            
    "longest": 4513562,                                                                                                           
    "mean": 135164.9344353537,                                                                                                   
    "median": 16518.0,                                                                                                            
    "sequence_count": 17967,                                                                                                     
    "shortest": 1,                                                                                                                
    "total_bps": 2428508377                                                                                                      
  },                                                                                                                              
  "Scaffold Stats": {                                                                                                            
    "L10": 66,                                                                                                                    
    "L20": 170,                                                                                                                  
    "L30": 308,                                                                                                                   
    "L40": 484,                                                                                                                  
    "L50": 710,                                                                                                                   
    "N10": 2767450,                                                                                                              
    "N20": 2026324,                                                                                                               
    "N30": 1554505,                                                                                                              
    "N40": 1220439,                                                                                                              
    "N50": 947285,                                                                                                                
    "gc_content": 43.814300913156785,                                                                                            
    "longest": 7315638,                                                                                                           
    "mean": 145748.69081197862,                                                                                                  
    "median": 15073.0,                                                                                                            
    "sequence_count": 16663,                                                                                                     
    "shortest": 391,                                                                                                              
    "total_bps": 2428610435                                                                                                      
  }                                                                                                                               
}	
	
															
																DAY 12: KRAKEN2 / CONTAMINATION

Runner: Renata
Look for contaminations : Kraken2 program. 	
WEBSITE: https://github.com/SmithsonianWorkshops/2020_4_20_STRI_genomics/blob/master/day4-genome_annotation_part1/Additional_Notes01.md
		 Check also: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
	
Program: kraken is a taxonomic sequence classifier that assigns taxonomic labels to DNA sequences. 
		- To check for contamination, how much of contamination we have.
		--confidence = allows the user to specify the threshold score in the interval [0,1]
		
		HELP: https://github.com/DerrickWood/kraken2/issues/399


--> 1 Job: Anolis_kraken.job

# /bin/sh                                                                                                                       
# ----------------Parameters---------------------- #                                                                            
#$ -S /bin/sh                                                                                                                   
#$ -pe mthread 10 
#$ -l mres=40G,h_data=4G,h_vmem=4G                                                                   
#$ -q sThC.q                                                                                                           
#$ -cwd                                                                                                                         
#$ -j y                                                                                                                         
#$ -N Anolis_kraken                                                                                                                
#$ -o Anolis_kraken.log                                                                                                            
#$ -m bea                                                                                                                       
#$ -M piranir@si.edu                                                                                                             
#                                                                                                                               
# ----------------Modules------------------------- #                                                                            
module load bioinformatics/kraken                                                                                                           
#                                                                                                                               
# ----------------Your Commands------------------- #                                                                            
#                                                                                                                               
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                     
echo + NSLOTS = $NSLOTS                                                                                           
#                                                                                                                 
kraken2 --db /data/genomics/db/Kraken/kraken2_db/ --report kraken_report --use-names --confidence 0 --threads $NSLOTS anolis_d1.fasta
#                                                                                                                 
echo = `date` job $JOB_NAME done	


RESULTS: Anolis_kraken.log                                                   
  15176 sequences root (91.08%)                                                                                                 
  1487 sequences unclassified (8.92%)                                                                                            

---------------------------------
Now you can filter those contigs 

# in the interactive queue you can run this.
qrsh -pe mthread 5

#first we will get only the list of contigs
grep "tarseq" kraken4.log > kraken_results_anolis4.txt

#Get list of human and unclassified contigs: By doing this, we are removing anything classified as bacteria, plasmids or viruses. 
We keep the unclassified because they might contain real contigs
grep "unclassified\|Homo sapiens" kraken_results_anolis4.txt | cut -f2 > anolis_r4_d1.list


#Using samtools, we will extract from the assembly only the contigs identified as Human or unclassified. Bacteria, 
viruses and plasmids will be excluded.
module load bioinformatics/samtools
xargs samtools faidx anolis_d1.fasta  < anolis_r4_d1.list > anolis_d1_r4_not_cont.fa


## Anolis genome after Kraken
{                                                                                                                                                
  "Contig Stats": {                                                                                                                              
    "L10": 96,                                                                                                                                    
    "L20": 250,                                                                                                                                  
    "L30": 454,                                                                                                                                   
    "L40": 713,                                                                                                                                  
    "L50": 1044,                                                                                                                                  
    "N10": 1903355,                                                                                                                               
    "N20": 1339109,                                                                                                                              
    "N30": 1054153,                                                                                                                               
    "N40": 827102,                                                                                                                               
    "N50": 644108,                                                                                                                                
    "gc_content": 43.80470629357358,                                                                                                             
    "longest": 4512383,                                                                                                                           
    "mean": 139597.28140848316,                                                                                                                   
    "median": 17245.0,                                                                                                                           
    "sequence_count": 17352,                                                                                                                      
    "shortest": 1,                                                                                                                               
    "total_bps": 2422292027                                                                                                                       
  },                                                                                                                                             
  "Scaffold Stats": {                                                                                                                             
    "L10": 66,                                                                                                                                    
    "L20": 170,                                                                                                                                  
    "L30": 306,                                                                                                                                   
    "L40": 482,                                                                                                                                  
    "L50": 707,                                                                                                                                   
    "N10": 2767002,                                                                                                                              
    "N20": 2026203,                                                                                                                               
    "N30": 1558597,                                                                                                                               
    "N40": 1222272,                                                                                                                              
    "N50": 948263,                                                                                                                                
    "gc_content": 43.80470629357358,                                                                                                             
    "longest": 7314449,                                                                                                                           
    "mean": 150927.9801246106,                                                                                                                   
    "median": 15611.5,                                                                                                                            
    "sequence_count": 16050,                                                                                                                      
    "shortest": 391,                                                                                                                             
    "total_bps": 2422394081  


RESULTS: Also here are the assembly stats before and after kraken.
We just remove 0.25% of the total length of the genome. This was include in 613 contigs that I have deleted from the assembly.



																		DAY 13: BUSCO AFTER KRAKEN2

Runner: piranir
Folder: /scratch/genomics/piranir/Kraken2/
--> 2 Job: busco3_afterKraken.job

# /bin/sh                                                                                                                       
# ----------------Parameters---------------------- #                                                                            
#$ -S /bin/sh                                                                                                                   
#$ -pe mthread 20 
#$ -l mres=100G,h_data=5G,h_vmem=5G,himem                                                                      
#$ -q  mThM.q                                                                                                                 
#$ -cwd                                                                                                                         
#$ -j y                                                                                                                         
#$ -N busco3_afterK                                                                                                                
#$ -o busco3_afterK.log                                                                                                            
#$ -m bea                                                                                                                       
#$ -M piranir@si.edu                                                                                                             
#                                                                                                                               
# ----------------Modules------------------------- #                                                                            
module load bioinformatics/busco/3.0.2                                                                                                           
#                                                                                                                               
# ----------------Your Commands------------------- #                                                                            
#                                                                                                                               
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                     
echo + NSLOTS = $NSLOTS                                                                                           
#                                                                                                                 
export AUGUSTUS_CONFIG_PATH="/scratch/genomics/piranir/Busco/augustus/config"                                                              
run_BUSCO.py -m genome -i anolis_d1_r4_not_cont.fa -o Anolis_afterK -l tetrapoda_odb9 -c $NSLOTS 
#                                                                                                                 
echo = `date` job $JOB_NAME done

--------------------------------

BUSCO RESULTS: 
 C:88.3%[S:87.3%,D:1.0%],F:6.9%,M:4.8%,n:3950                                                                                                  
        3487    Complete BUSCOs (C)                                                                                                                   
        3448    Complete and single-copy BUSCOs (S)                                                                                                                                                                                                                                                 
        39      Complete and duplicated BUSCOs (D)                                                                                                    
        274     Fragmented BUSCOs (F)                                                                                                                 
        3950    Total BUSCO groups searched 


RESULTS: 

- Check your assembly with assembly-stats
- module load bioinformatics/assembly_stats											

[piranir@hydra-login02 Kraken2]$ assembly_stats anolis_d1_r4_not_cont.fa                         

THE SAME AS ABOVE (AFTER KRAKEN2)                                                                                                                             

    

	
	FAZER UM DOCUMENTO NO WORD COM TODOS OS RESULTADOS DO ANTES E DEPOIS DO GENOMA QUE TEMOS AGORA. 
	USANDO O GRAFICO .HISTO QUE FIZEMOS LÁ NO INÍCIO. 
	
	
																DAY 13: ANNOTATION 01
	
	repeatmasker (start annotation) -species (use something close to your orgAnism)
	
	AFTER USE BLAT WITH THE TRANSCRIPTOME
	
	
	- salsa before the annotation. 
	-polishing and annotation
	
	
	
	
	
	
	
	
	
	
 
	
	
	
	
	
	
	
	
	
	
													