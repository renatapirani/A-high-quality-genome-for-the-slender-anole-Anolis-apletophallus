### New high-quality genome assembly of the slender anole (Anolis apletophallus)

#### 10x data Illumina (short reads) + Nanopore (long reads) 

#### Developed by Dr. Renata Pirani (renatapirani) and Dr. Carlos Arias (solracarias)
 
# DAY 1: CUTADAPT

* FUNCTION: we are cleaning the DNA  (short reads) cutting the bad sequencing and checking for any contamination
* PROGRAM WEBSITE: https://cutadapt.readthedocs.io/en/stable/guide.html

* JOB FILE: /scratch/genomics/piranir/Cutadapt/cutadapt_26.job
 		/scratch/genomics/charleskl/CutAdapt/cutadapt_10.job

	+ **module**: ```module load bioinformatics/cutadapt/2.4```
		
	+ **command**: ```cutadapt -u 26 -o E28_26t_val_1.fq.gz /scratch/stri_ap/ariasc_data/anolis_10x/E28_MPS12345004_G06_9489_S3_L004_R1_001_val_1.fq.gz```	
	
	+ **command**: ```cutadapt -u 26 -o E28_10t_R2_001_val_2.fq.gz /scratch/stri_ap/ariasc_data/anolis_10x/E28_MPS12345004_G06_9489_S3_L004_R1_001_val_1.fq.gz```	


#### RESULTS
=== Summary ===                                                                                                                   
                                                                                                                                  
Total reads processed:             358,209,584                                                                                    
Reads written (passing filters):   358,209,584 (100.0%)                                                                           
                                                                                                                                  
                                                                                                                                  
Total basepairs processed: 52,891,743,310 bp                                                                                      
Total written (filtered):  43,579,357,893 bp (82.4%)
 											
 											
 											
# DAY 2: Jellyfish
 											
* FUNCTION: Run Genomescope, first you need to generate a Jellyfish histogram
		  Genomescope can be used to estimate genome size from short read data
* PROGRAM WEBSITE: http://qb.cshl.edu/genomescope/

* JOB FILE: /scratch/genomics/charleskl/jellyfish/jellyfish3.job
		  
	+ **module**: ```module load bioinformatics/jellyfish/2.3.0```  
		                                                                                                                                    
	+ **command**: ```gzip -dc E28_26t_val_1.fq.gz E28_10t_R2_001_val_2.fq.gz | jellyfish count -C -m 21 -s 800000000 -t $NSLOTS -o reads.jf /dev/fd/0```  



* JOB FILE: /scratch/genomics/piranir/Jellyfish/hist_jellyfish.job

	+ **module**: ```module load bioinformatics/jellyfish ```
		
	+ **command**: ```jellyfish histo -t 20 reads_EF_illumina.jf > reads_EF_ill.histo``` 



NEXT STEPS: download the files from the hydra to your computer:
- go to the genome website and upload your file (reads.histo)


# DAY 3: Wtdbg2/Redbeans


* FUNCTION: Wtdbg2 is a sequence assembler for long noisy reads produced by either PacBio or Oxford Nanopore Technologies.
		 The program uses two steps. The first step is the assembler and the second step is the consenser. 
* PROGRAM WEBSITE: https://github.com/ruanjue/wtdbg2/blob/master/README-ori.md
														
* 1 JOB FILE: /scratch/genomics/piranir/Wtdbg2/readbeans.job
	- To check the parameters go to "wtdbg2 --help" 
 

	+ **module**: ```module load bioinformatics/wtdbg2```
		
	+ **command**: ```wtpoa-cns -t $NSLOTS -i anolis_28_S2_assembly.ctg.lay.gz -fo anolis_S2_genomedraft_raw.fa ```
		
			
	P.S. to check queues = Go to confluence.si.edu and to "Submitting Jobs" and "available queues". 
	
	END: to check how is your assembly data type: 
	module load bioinformatics/assembly_stats
	assembly_stats anolis_S2_genomedraft_raw.fa
	

#### RESULTS
assembly_stats anolis_S2_genomedraft_raw.fa                                                                 
                                                                                                                                             
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
    "mean": 119098.24670711854                                                                                                                            
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
    "median": 19088.0,70711854                                                                                                                             
    "sequence_count": 20271,                                                                                                                 
    "shortest": 952,                                                                                                                         
    "total_bps": 2414240559                                                                                                                  
} } 


										
# DAY 4: Minimap/bwa


* FUNCTION: Polish the Genome 
		  - minimap2:  is a versatile sequence alignment program that alings DNA or mRNA sequences against a large reference database. 

* PROGRAM WEBSITE: https://github.com/lh3/minimap2
				 https://github.com/isovic/racon
				 http://www.htslib.org/doc/samtools-view.html
				 https://github.com/lh3/bwa
											
* 1 JOB FILE: /scratch/genomics/charleskl/jobs/polish_anolis_S2.job

	- We will used long reads data and short reads data to polish your draft genome by using these scripts:
	- For this step we are going to run 4 jobs. The first 2 jobs are for the long reads. 
	- For long reads, you can used minimap2 to map your long reads to your draft genome, after 
		that you will need to sort and filtered by reads that have the higher match to your draft genome, 
		finally you will run de concensus command.


	+ **module**: ```source /home/ariasc/.bashrc```                                                                                                   
					 ```conda activate minimap2```
					 ```module load bioinformatics/samtools```
					
	+ **command**: ```minimap2 -t40 -ax map-ont /scratch/genomics/charleskl/redbean/anolis_S2_genomedraft_raw.fa ```
					   ```/scratch/stri_ap/ariasc_data/anolis_nanopore/anolis_28_filt.fa.gz | samtools sort -@4 >anolis_S2_dbg.bam```
                                                                                                                        

P.S. wait for the results of the Job 1 to start the job 2. 

	
* 2 JOB: polish_anolis_S2_wtpoa.job
	polishing with long reads, filetering mapping reads, and consensus with wtpoa-cns module in redbean

	+ **module**: ```module load bioinformatics/samtools```
		
	+ **command**: ```samtools view -F0x900 /scratch/genomics/charleskl/jobs/anolis_S2_dbg.bam | /home/ariasc/programs/wtdbg2/./wtpoa-cns -t 40 -d /scra```
					   ```tch/genomics/charleskl/redbean/anolis_S2_genomedraft_raw.fa -i - -fo Anolis_dbg_cns.fa```                                                                                                                        
                                                                                                                         
* 3 JOB: bwa_index.job
* Indexing draft genome in bwa

	+ **module**: ```module load bioinformatics/bwa/0.7.17 ```
		
	+ **command**: ```bwa index /scratch/genomics/charleskl/jobs/Anolis_dbg_cns.fa```


* SHORT READS
		- Now for short reads, you can used bwa to map your short reads to your draft genome, after 
	  		that you will sort the mapped file and finally run the concesus command.
		- We are going to use the BWA program. BWA is a software package for mapping DNA sequences against a large reference genome.
		- BWA first needs to construct the FM-index for the reference genome (the index command)


* 4 JOB: polish_anolis_small.job
* bwa polishing with wtpoa-cns module from redbean

	+ **module**: ```module load bioinformatics/bwa/0.7.17```
 					 ``` module load bioinformatics/bwa/0.7.17``` 
 					                                                                                                                        
 	+ **command**: ```bwa mem -t 40 /scratch/genomics/charleskl/jobs/Anolis_dbg_cns.fa /scratch/genomics/charleskl/E28_26t_val_1.fq.gz /scratch/genomics```
					 ```/charleskl/E28_10t_R2_001_val_2.fq.gz | samtools sort -O SAM | /home/ariasc/programs/wtdbg2/./wtpoa-cns -t 40 -x sam-sr -d /scratc```
					 ```h/genomics/charleskl/jobs/Anolis_dbg_cns.fa -i - -fo anolis_dbg.srp.fa```


END: to check how is your assembly data type: 
		module load bioinformatics/assembly_stats
		assembly_stats anolis_dbg.srp.fa
	

#### RESULTS

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

							
# DAY 5: Scaff10x					
											

* FUNCTION: Pipeline for scaffolding and breaking a genome assembly using 10x genomics linked-reads
		  make sure that the version is gcc CC= /software/gcc-4.9.2/bin/gcc in the makefile 

* PROGRAM WEBSITE: https://github.com/wtsi-hpag/Scaff10X.git

* INSTALL: $ git clone  https://github.com/wtsi-hpag/Scaff10X.git 
		 $ cd Scaff10X
		 $ ./install.sh			 

* JOB FILES: /scratch/genomics/charleskl/jobs/Anolis_dbg_cns.fa
		   /scratch/stri_ap/ariasc_data/anolis_10x


* input2.dat									
q1=/scratch/stri_ap/ariasc_data/anolis_10x/E28_MPS12345004_G06_9489_S3_L004_R1_001_val_1.fq.gz                          
q2=/scratch/stri_ap/ariasc_data/anolis_10x/E28_MPS12345004_G06_9489_S3_L004_R2_001_val_2.fq.gz

P.S. Best results -block 100000 



* 1 JOB: /scratch/genomics/piranir/Scaff10/scaff10xAnolis.job
		  /scratch/genomics/charleskl/Scaff10x/scaff10x.job     


	+ **command**: ```/home/ariasc/programs/Scaff10X/src/./scaff10x -nodes 30 -align bwa -score 10 -read-s1 5 -read-s2 5 -longread 1 -gap 100 -block 20000```
 					   ```-data input2.dat Anolis_dbg_cns.fa Anolis_20k_scaf1.fa```


- Check your assembly with assembly-stats
- module load bioinformatics/assembly_stats
- assembly_stats Anolis_20k_scaf1.fa


#### RESULTS
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


P.S. - Results N50 = 997812 (file: Anolis_20k_scaf1.fa), so we will use this data to polish.
	AFTER: polish with Racun for long read data, we use minimap




# DAY 6: minimap2/racon long reads


* FUNCTION: Polish the Genome 
		  minimap2: is a versatile sequence alignment program that alings DNA or mRNA sequences against a large reference database. 
		  Polishing with Racon -- Type of data: Nanopore or Pacbio

* PROGRAM WEBSITE: https://github.com/lh3/minimap2
				
* JOB FILE: /scratch/genomics/charleskl/Racon/Anolis_20k_scaf1.fa

	P.S. file used: Anolis_20k_scaf1.fa



--> 1 job: anolis_minimap_racon.job 
	First map your clean_raw reads to your assembly with minimap.
	
+ **module**: ```source /home/ariasc/.bashrc```
		       ```conda activate minimap2```	
		       
	+ **command**: ```minimap2 -ax map-ont -t40 /scratch/genomics/charleskl/Racon/Anolis_20k_scaf1.fa /scratch/stri_ap/ariasc_data/anolis_nanopore/anolis_28_filterd.fastq > anolis_mapped_data.sam```

     

* 2 JOB: racon.job


	+ **module**: ```source /home/ariasc/.bashrc```
					  ```conda activate racon```
					  
	+ **command**: ```racon -u -t 30 /scratch/stri_ap/ariasc_data/anolis_nanopore/anolis_28_filterd.fastq anolis_mapped_data.sam /scratch/genomics/charleskl/Racon/Anolis_20k_scaf1.fa > anolisgenome_racon_polished.fasta```


#### RESULTS: 

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


* ERROR: we got warnings about 2 contigs that were possibly chimeric.
		We will run Break10x to fix the error. 
		

#### Break10x 

* 3 job: /scratch/genomics/piranir/Scaff10/break10xAnolis.job
	
	+ **command**: ```/home/ariasc/programs/Scaff10X/src/./break10x -nodes 30 -score 10 -gap 100 -data input2.dat ```
					```anolisgenome_racon_polished.fasta scaffolds-break.fasta scaffolds-break.name```                                                                              
                                                                         


#### RESULTS: 
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



# DAY 7: BUSCO3 from short reads before Polishing


* FUNCTION: assess the completeness of genomes, gene sets, and transcriptomes, 
		  using their gene content as a complementary method to common technical metrics. 
		 
* PROGRAM WEBSITE: https://busco.ezlab.org/busco_userguide.html

* JOB FILE: /scratch/genomics/piranir/Busco/ busco_before.job

* BEFORE RUNNING: 
	- Running last version of Busco that I have installed
	- To run busco4 you need to first to cp config file from augustus. Please do this:
	- make new directory in your path 
	- cd /scratch/genomics/username/
	- mkdir busco
	- cd busco
	- and Assuming you are in the folder `/scratch/genomics/username/busco`:
	- cp -r /share/apps/bioinformatics/augustus/conda/3.3.2/config/ .


  + **module**: ```module load bioinformatics/busco/3.0.2```	
  		
 + **command**: ```export AUGUSTUS_CONFIG_PATH="/scratch/genomics/piranir/Busco/augustus/config"```                                                        
		 ```run_BUSCO.py -m genome -i scaffolds-break.fasta -o Anolis_before -l tetrapoda_odb9 -c $NSLOTS```


#### RESULTS:
Summarized benchmarking in BUSCO notation for file scaffolds-break.fasta                                                                    
BUSCO was run in mode: genome                                                                                                               
                                                                                                                                              
  C:67.5%[S:66.7%,D:0.8%],F:16.5%,M:16.0%,n:3950                                                                                     
        2664    Complete BUSCOs (C)                                                                                                           
        2633    Complete and single-copy BUSCOs (S)                                                                                           
        31      Complete and duplicated BUSCOs (D)                                                                                            
        652     Fragmented BUSCOs (F)                                                                                                         
        634     Missing BUSCOs (M)                                                                                                            
        3950    Total BUSCO groups searched 



# DAY 8: Minimap2/Racon short reads
														
														
* FUNCTION: polishing with racon Illumina short reads
 	  
* PROGRAM WEBSITE: https://github.com/isovic/racon
				 https://github.com/lh3/minimap2

* 1 JOB FOLDER: /scratch/genomics/piranir/Racon_shortreads/minimap.job			

* FIRST: 
	- rename reads 1 a 2 to used them with racon.
	- for questions go to /scratch/stri_ap/ariasc_data/anolis_10x and the files are there.
	- for the python file use the path /scratch/genomics/ariasc/Efish/
	- copy those files to your folder and unzip them using - gunzip *.gz 



 	+ **command**: ```python fix_reads_racon.py E28_10t_R2_001_val_2.fq E28_26t_R1_001_val_1.fq > rename_x_racon.fastq```



	- second, map reads to draft genome with minimap2



 	+ **module**: ```source /home/ariasc/.bashrc ```
  	+ **module**: ```conda activate minimap2```
  	+ **command**: ```minimap2 -ax sr -t $NSLOTS /scratch/genomics/piranir/Racon_shortreads/scaffolds-break.fasta rename_x_racon.fastq > Anolis_notC_scaff.sam```
  		
                                                                                                    
* 2 JOB FOLDER: /scratch/genomics/piranir/Racon_shortreads/racon_short.job
	run racon with the *.sam file and corrected names short reads files

  
 	+ **module**: ```source ~/.bashrc ``` 
 	+ **module**: ```conda activate racon```               

 	+ **command**: ```racon -t 30 rename_x_racon.fastq Anolis_notC_scaff.sam scaffolds-break.fasta > Anolis_racon_short.fasta```

                                                                                                                    
                                                                                                              
#### RESULTS: 
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


# DAY 9: BUSCO 3 from short reads AFTER Polishing
														
* FUNCTION: assess the completeness of genomes, gene sets, and transcriptomes, 
		  using their gene content as a complementary method to common technical metrics. 
		 
* PROGRAM WEBSITE: https://busco.ezlab.org/busco_userguide.html

* 1 JOB FILE: /scratch/genomics/piranir/Busco/busco3_after.job

* BEFORE RUNNING: 
	- Running last version of Busco that I have installed
	- To run busco4 you need to first to cp config file from augustus. Please do this:
	- make new directory in your path 
	- cd /scratch/genomics/username/
	- mkdir busco
	- cd busco
	- and Assuming you are in the folder `/scratch/genomics/username/busco`:
	- cp -r /share/apps/bioinformatics/augustus/conda/3.3.2/config/ .
 
 * SENT FILES TO HYDRA
	- to send files to hydra go: https://confluence.si.edu/display/HPC/Using+rclone
	- or you can use wget "website link"


	+ **module**: ```module load bioinformatics/busco/3.0.2```
                                                                                                       
  	+ **command**: ```export AUGUSTUS_CONFIG_PATH="/scratch/genomics/piranir/Busco/augustus/config"```                                                          
			```run_BUSCO.py -m genome -i Anolis_racon_short.fasta -o Anolis_after -l tetrapoda_odb9 -c $NSLOTS ```
                                                                                                         					
														
#### RESULTS:

BUSCO was run in mode: genome                                                                                                               
                                                                                                                                              
  C:88.1%[S:87.1%,D:1.0%],F:7.0%,M:4.9%,n:3950     
        3482    Complete BUSCOs (C)                                                                                                           
        3441    Complete and single-copy BUSCOs (S)                                                                                           
        41      Complete and duplicated BUSCOs (D)                                                                                            
        277     Fragmented BUSCOs (F)                                                                                                         
        191     Missing BUSCOs (M)                                                                                                            
        3950    Total BUSCO groups searched  												
														
		
NEXT STEPS:
	- We decided to run polishing again for the Anolis_racon_short.fasta file		
	- We will run pilon, minimap and Racon again 
	- After download Anolis carolinensis and A. sagrei genomes and compare to our results. 


* FILE: Anolis_racon_short.fasta
												

												
# DAY 10: Minimap2/Racon to improve the BUSCO RESULTS
														
														
* 1 JOB FOLDER: /scratch/genomics/piranir/Racon_shortreads/minimap2.job			

* FILE: Anolis_racon_short.fasta

 	+ **module**: ```source /home/ariasc/.bashrc```
           		```conda activate minimap2 ```
                                                                                                                  
 	+ **command**: ```minimap2 -ax sr -t $NSLOTS /scratch/genomics/piranir/Racon_shortreads/Anolis_racon_short.fasta rename_x_racon.fastq > Anolis_notC_scaff2.sam ``` 
  		
                                                                                                  

* 2 JOB: racon2.job
	run racon with the *.sam file and corrected names short reads files

 	+ **module**: ```source ~/.bashrc```
           			```conda activate racon ```	
                                                                                                                  
 	+ **command**: ```racon -t 30 rename_x_racon.fastq Anolis_notC_scaff2.sam Anolis_racon_short.fasta > Anolis_after_polishing2.fasta```                                                                                                       
                                                                                   

#### RESULTS: 

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
	

												
# DAY 11: BUSCO 3 AFTER2
														
	
* FUNCTION: assess the completeness of genomes, gene sets, and transcriptomes, 
		  using their gene content as a complementary method to common technical metrics. 
		 
* PROGRAM WEBSITE: https://busco.ezlab.org/busco_userguide.html

* 1 JOB FILE: /scratch/genomics/piranir/Busco/Anolis_after_polishing2.fasta
* input: busco3_after2.job


 	+ **module**: ```module load bioinformatics/busco/3.0.2```
                                                                                                          

 	+ **command**: ```export AUGUSTUS_CONFIG_PATH="/scratch/genomics/piranir/Busco/augustus/config"```                                                           
					```run_BUSCO.py -m genome -i Anolis_after_polishing2.fasta -o Anolis_after2 -l tetrapoda_odb9 -c $NSLOTS ```

													
																												
#### RESULTS:	
									
C:88.1%[S:87.1%,D:1.0%],F:7.1%,M:4.8%,n:3950                                                                                                         
        3439    Complete and single-copy BUSCOs (S)                                                               
        40      Complete and duplicated BUSCOs (D)                                                                                                          
        279     Fragmented BUSCOs (F)                                                                             
        192     Missing BUSCOs (M)											
														

															
# DAY 12: PILON
															
* FUNCTION: Automatically improve draft assemblies
		   
		 
* JOB FILE: /scratch/genomics/ariasc/anolis/anolis_d1.fasta
														

	+ **module**: ```module load bioinformatics/pilon/1.23```
		
   + **command**: ```PILON_HEAP_SIZE=1000g```                                                                                                           
		```runpilon --genome /scratch/genomics/ariasc/anolis/Anolis_racon_short.fasta --bam /scratch/genomics/ariasc/anolis/Anolis_racon_short_sorted.bam ```
		```--output  anolis_d1 --outdir /scratch/genomics/ariasc/anolis/ --fix bases --changes --vcf --threads $NSLOTS > /scratch/genomics/ariasc/anolis/anolis_d1_out.log ```                 

											
														

#### RESULTS: anolis_d1.fasta														
															
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
	
															
# DAY 13: KRAKEN2 / CONTAMINATION

* Look for contaminations : Kraken2 program. 	
* WEBSITE: https://github.com/SmithsonianWorkshops/2020_4_20_STRI_genomics/blob/master/day4-genome_annotation_part1/Additional_Notes01.md
		 Check also: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
	
* Program: kraken is a taxonomic sequence classifier that assigns taxonomic labels to DNA sequences. 
		-- To check for contamination, how much of contamination we have.
		--confidence = allows the user to specify the threshold score in the interval [0,1]
		
* HELP: https://github.com/DerrickWood/kraken2/issues/399


* 1 JOB: /scratch/genomics/piranir/Kraken2/Anolis_kraken.job

     
	+ **module**: ```module load bioinformatics/kraken``` 
                                                                                                           
 	+ **command**: ```kraken2 --db /data/genomics/db/Kraken/kraken2_db/ --report kraken_report --use-names --confidence 0 --threads $NSLOTS anolis_d1.fasta``` 
                                                                                                       

#### RESULTS: Anolis_kraken.log                                                   
  15176 sequences root (91.08%)                                                                                                 
  1487 sequences unclassified (8.92%)                                                                                            

---------------------------------
Now you can filter those contigs 

* in the interactive queue you can run this.
	qrsh -pe mthread 5

* first we will get only the list of contigs
	grep "tarseq" kraken4.log > kraken_results_anolis4.txt

* Get list of human and unclassified contigs: By doing this, we are removing anything classified as bacteria, plasmids or viruses. 
	We keep the unclassified because they might contain real contigs
	grep "unclassified\|Homo sapiens" kraken_results_anolis4.txt | cut -f2 > anolis_r4_d1.list


* Using samtools, we will extract from the assembly only the contigs identified as Human or unclassified. Bacteria, 
	viruses and plasmids will be excluded.
	module load bioinformatics/samtools
	xargs samtools faidx anolis_d1.fasta  < anolis_r4_d1.list > anolis_d1_r4_not_cont.fa


#### RESULTS: Anolis genome after Kraken
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



# DAY 14: BUSCO AFTER KRAKEN2

* Folder: /scratch/genomics/piranir/Kraken2/
* 2 JOB: busco3_afterKraken.job

 	+ **module**: ```module load bioinformatics/busco/3.0.2 ```

 	+ **command**: ```export AUGUSTUS_CONFIG_PATH="/scratch/genomics/piranir/Busco/augustus/config" ```                                                             
			```run_BUSCO.py -m genome -i anolis_d1_r4_not_cont.fa -o Anolis_afterK -l tetrapoda_odb9 -c $NSLOTS ```                                                                                               
                                                                                          


--------------------------------

#### FINAL BUSCO RESULTS: 
 C:88.3%[S:87.3%,D:1.0%],F:6.9%,M:4.8%,n:3950                                                                                                  
        3487    Complete BUSCOs (C)                                                                                                                   
        3448    Complete and single-copy BUSCOs (S)                                                                                                                                                                                                                                                 
        39      Complete and duplicated BUSCOs (D)                                                                                                    
        274     Fragmented BUSCOs (F)                                                                                                                 
        3950    Total BUSCO groups searched 


RESULTS: 

- Check your assembly with assembly-stats
- module load bioinformatics/assembly_stats											

 
													
