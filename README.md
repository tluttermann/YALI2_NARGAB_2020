# Hybrid genome assembly and genome analyses of the yeast *Yarrowia lipolytica* DSM 3286
The following text contains information on the command lines that were used for the hybrid genome assembly as well as for the genome analyses of *Y. lipolytica* DSM 3286. 

## Hybrid genome assembly

1. Assembly of ONT reads using the CANU long read assembler v.1.6

        canu -p Yli genomeSize=21M useGrid=1 -nanopore-raw ONT*.fastq 'gridEngineThreadsOption=-pe multislot THREADS' 'gridEngineMemoryOption=-l m_mem_total=MEMORY' 'gridEngineSubmitCommand=qsub -P canu_assembly -l idle=1' 'maxThreads=32' 'maxMemory=1024'

2. Polishing with NANOPOLISH

        minimap2 --secondary=no -a -x map-ont Yli.contigs.fasta ONT_*.fastq | samtools sort > ONT_vs_contigs.bam
        
        samtools index ONT_vs_contigs.bam
        
        multi_to_single_fast5 –recursive -i ONT_Seqdata/ -t 32 -s ONT_SingleFast5/
        
        nanopolish extract -d ONT_SingleFast5 -o ONT_Reads.fasta
        
        nanopolish index -d ONT_SingleFast5 ONT_Reads.fasta
        
        nanopolish variants --consensus -r ONT_Reads.fasta -b ONT_vs_contigs.bam -g Yli.contigs.fasta -o Yli.contigs.vcf -t 32 -p 1
        
        nanopolish vcf2fasta -g Yli.contigs.fasta Yli.contigs.vcf > Yli.nanopolished.fasta

3. Polishing with PILON using Illumina data until no further changes are made

        for i in $(seq -f "%02g" 1 20); do
        
        j=$(printf "%02g" $(expr $i + 1));
        
          if [ -e $1/Round$j.fasta ]; then
          
            continue;
            
           else
           
            j=`expr $i + 1`;

			      bwa index Round$i.fasta Round$i.fasta;

			      bwa mem -O1 -E1 -t16 Round$i.fasta TSPf_R1.fastq.gz \
			      TSPf_R2.fastq.gz | /vol/biotools/bin/samtools sort \
			      --threads 16 -o WGS.Round$i.sorted.bam;

			      samtools index WGS.Round$i.sorted.bam;
			      bwa mem -O1 -E1 -t16 Round$i.fasta MP_R1.fastq.gz \
			      MP_R2.fastq.gz | /vol/biotools/bin/samtools sort \
			      --threads 16 -o MP.Round$i.sorted.bam;

			      samtools index MP.Round$i.sorted.bam;

			      java -Xmx80G -jar pilon-1.22.jar --genome Round$i.fasta \
			      --fix all --changes --frags WGS.Round$i.sorted.bam \
			      --jumps MP.Round$i.sorted.bam --threads 16 –output \
			      Round$j | tee Round$i.pilon;

		      fi;

		      sed -i "s/_pilon//g" $1/Round$j.fasta;

		      if ! [ -s $1/Round$j.changes ]; then
			      break
		      fi;
          
	      done
        
	      for i in $(seq -f "%02g" $j $(expr $j + 20)); do
        
		      j=$(printf "%02g" $(expr $i + 1));

		      if [ -e $1/Round$j.fasta ]; then

			      continue;

		      else

			      j=`expr $i + 1`;

			      bowtie2-build --threads 16 Round$i.fasta Round$i > \ 	
            /dev/null;

		       	bowtie2 -X 1000 -x Round$i -1 TSPf_R1.fastq.gz -2 \ 	
            TSPf_R2.fastq.gz --threads 16 2> Round$i.bowtie | \
		      	samtools sort --threads 16 -o WGS.Round$i.sorted.bam;

			      samtools index WGS.Round$i.sorted.bam;
  
			      bowtie2 -X 10000 -x Round$i -1 MP_R1.fastq.gz -2 \ 	
            MP_R2.fastq.gz --threads 16 2>> Round$i.bowtie | \
			      samtools sort --threads 16 -o MP.Round$i.sorted.bam;

			      samtools index MP.Round$i.sorted.bam;

			      rm Round$i*.bt2;

			      java -Xmx80G -jar pilon-1.22.jar --genome Round$i.fasta \
			      --fix all --changes --frags WGS.Round$i.sorted.bam \
			      --jumps MP.Round$i.sorted.bam --threads 16 –output \
			      Round$j | tee Round$i.pilon;

		      fi;

		      sed -i "s/_pilon//g" $1/Round$j.fasta;

		      if ! [ -s $1/Round$j.changes ]; then
			      break
		      fi;
	      done

## Analyses of the rDNA clusters and the mitochondrial genome
The following command line sequence was used for parts of the analysis of the rDNA clusters as well as of the mitochondrial genome of Y. lipolytica DSM 3286. The sequence was used to identify reads that harbor specific features.

1. Conversion of the fastq files obtained from the ONT sequencing runs to fasta files:

        sed -n '1~4s/^@/>/p;2~4p' ONT_run_x.fastq > ONT_run.fasta
  
2. Creation of a BLAST database of the sequencing reads:

        makeblastdb -in ONT_run_x.fasta -dbtype nucl -parse_seqids
  
3. Search for reads with specific features with suitable query sequences (such as the cox1 gene of YALI2 for filtering for mitochondrial reads and the 18S sequence for filtering for rDNA reads; bioproject PRJNA641784):

        blastn -db ONT_run_x.fasta -query query_y_sequence.fasta -outfmt '6 qseqid sseqid pident qstart qend length qlen' > query_y_reads_from_ONT_run_x.txt
  
4. Extraction of the identified reads from the fasta file of the specific ONT sequencing run:

        grep -A1  -f  query_y_reads_from_ONT_run_x.txt ONT_run_x.fasta > query_y_reads.fasta
  

For subsequent analyses both the text files from step 3 and the fasta files from step 4 have been used. Results were visualized with OriginLab software. For detailed analyses, the feature annotation tool of Geneious Prime® 2020.1.2 and the MUSCLE and MAUVE plugins were used.
