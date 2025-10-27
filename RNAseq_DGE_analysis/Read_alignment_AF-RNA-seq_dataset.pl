#!/usr/bin/perl -w
# Program to generate a jobfiles for read alignment of a batch of fastqfiles
# and send the jobfiles to the QMUL apocrita cluster


# Would need to be executed once for each matrix database
use strict;
use warnings;

# [PATHS]
my $homedir="/path/to/lab/directory/";
my $jobs_dir=$homedir."out/";
my $input_dir=$homedir."repository/";


# [PATHS]
my $homedir="/data2/Blizard-VillarLab/";
my $jobs_dir=$homedir."out/";
my $input_dir=$homedir."repository/";

# [PARAMETERS]
my $sqt = "pe";  # Paired-end data
#my $sqt = "se";  # Single-end data


# Read input data from file
my @tmp; my @jobs; my @tmp2;
my ($lane, $i, $line, $jobfile, $string, $jobname, $jobid, $filename); 
open (IN, "GSE128188_AF-Thomas.txt") || die 
"Could 
not 
open input 
file!\n";
  while (<IN>) {
     chomp ($_);
     if (/# /){next;}
     if (/LibraryID/){next;}
     @tmp=split (/\t/,$_); # $tmp[1] LibraryID $tmp[2] SubmissionID $tmp[4] ExpID $tmp[5] Genome $tmp[6] Strain $tmp[7] Tissue $tmp[9] Individual $tmp[10] Sex $tmp[11] Factor
     print "$tmp[0] $tmp[1] $tmp[11]\n";
     # Generate name for output files based on info provided in Inventory
     $filename = $tmp[1]."_".$tmp[11]."_".$tmp[7]."_".$tmp[5]."_".$tmp[10]."_".$tmp[4]."_".$tmp[9]."_".$tmp[2];
     #LibraryID_Factor_Tissue_Genome_Sex_Experiment_AssayType
     print "$filename\n";
     #exit;
     # Rename Novogene fastq files for each library
     if ($sqt eq "pe"){
     $string = "mv $input_dir"."$tmp[1]/$tmp[1]*1.fastq.gz $input_dir"."$tmp[1]/".$filename."_1.fastq.gz";
     print "$string\n";
     system ($string);
     $string = "mv $input_dir"."$tmp[1]/$tmp[1]*2.fastq.gz $input_dir"."$tmp[1]/".$filename."_2.fastq.gz";
      print "$string\n";
      system ($string);
      }
     else{# Single-end data - path to fastq file
      $string = "$input_dir"."$tmp[1]/*.fastq.gz";
      print($string);
     }
     # Generate and submit job files for alignment of each library
     #$name = $jobs_dir.$tmp[0]."_jobfile";
     $jobname = $tmp[1]."_jobfile";
     open (JOB, ">".$jobname."_1") || die "Could not open output job file!\n";
     
     # Step 1. Salmon pseudo-alignment and quantification
     print JOB &WriteHeader($tmp[1], "1", "1"); # Print job file header via subroutine - specify library, step and number of cores
     #print JOB ("#\$ -l h_vmem=16G","\n"); # 16GB RAM for all jobs
     print JOB ("module load salmon\n"); #Load module into cluster
     print JOB ("mkdir $tmp[1]"."_STAR\n"); #Make directory for results
      # Runs quant function, specifying library type, location of input reads, name of output file and default mapping options
     if($sqt eq "pe"){
      print JOB ("salmon quant -i $homedir"."genomes/ensembl/indexed/$tmp[5].cdna.all.index \\
      -l ISF \\
      -1 $input_dir"."$tmp[1]/$filename"."_1.fastq.gz -2 $input_dir"."$tmp[1]/$filename"."_2.fastq.gz \\
      -o $jobs_dir"."$tmp[1]/$tmp[1]\.salmon --seqBias --useVBOpt --validateMappings\n");
                 }
     close JOB;
     # Submit to cluster
     $string="qsub ".$jobname."_1";
     print $string."\n";
     system ($string);
     # Store job id
     system ("qstat > job_ids.txt");
     
     # Step 2. Run STAR alignment
     open (JOB, ">".$jobname."_2") || die "Could not open output job file!\n";
     print JOB &WriteHeader($tmp[1], "2", "6"); # Print job file header via subroutine - specify library, step and number of cores
     print JOB ("module load star/2.7.0f\n");
     if($sqt eq "pe"){
      print JOB ("gunzip $input_dir"."$tmp[1]/$filename"."_1.fastq.gz\n"); #Uncompress fastqfiles for STAR alignment
      print JOB ("gunzip $input_dir"."$tmp[1]/$filename"."_2.fastq.gz\n");
      #Runs STAR alignment using indexed human genome, input files, output folder and default output options for a bam file sorted by coordinates
      print JOB ("STAR --genomeDir $homedir"."genomes/ensembl/indexed/HG38_STAR/ \\
      --runThreadN 6 \\
      --readFilesIn $input_dir"."$tmp[1]/$filename"."_1.fastq $input_dir"."$tmp[1]/$filename"."_2.fastq \\
      --outFileNamePrefix $jobs_dir"."$tmp[1]/ \\
      --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard\n");
      print JOB ("gzip -f $input_dir"."$tmp[1]/$filename"."_1.fastq\n"); #Re-compress fastq files
      print JOB ("gzip -f $input_dir"."$tmp[1]/$filename"."_2.fastq\n");
      }
     close JOB;
     # Submit to cluster as dependent on the previous job - will only start after salmon has finished
     $jobid = &FetchJobID();
     $string="qsub -hold_jid $jobid $jobname"."_2";
     print $string."\n";
     system ($string);
     # Store job id
     system ("qstat > job_ids.txt");


    # Step 3. QC STAR bam file
     open (JOB, ">".$jobname."_3") || die "Could not open output job file!\n";
     print JOB &WriteHeader($tmp[1], "3", "1"); # Print job file header via subroutine - specify library, step and number of cores
     print JOB ("module load samtools\n");
     print JOB ("samtools flagstat $jobs_dir"."$tmp[1]/Aligned.sortedByCoord.out.bam > $jobs_dir"."$tmp[1]/$tmp[1]_STARmappingStats.txt\n"); #Basic QC of STAR bam file saved to a text file
     print JOB ("samtools sort -n -O bam -o $jobs_dir"."$tmp[1]/$filename"."_sorted.bam -T temp $jobs_dir"."$tmp[1]/Aligned.sortedByCoord.out.bam \n"); #Create bam file sorted by read name for featureCounts
     print JOB ("samtools index $jobs_dir"."$tmp[1]/$filename"."_sorted.bam\n");
     print JOB ("mv $jobs_dir"."$tmp[1]/$filename"."_sorted.bam $jobs_dir"."$tmp[1]/$filename.bam\n"); #Index bam file and save
     print JOB ("samtools flagstat $jobs_dir"."$tmp[1]/$filename.bam > $jobs_dir"."$tmp[1]/$tmp[1]_SORTmappingStats.txt\n"); #Confirm QC on sorted bam file
     close JOB;
     
     # Submit to cluster as dependent on the previous job
     $jobid = &FetchJobID();
     $string="qsub -hold_jid $jobid $jobname"."_3";
     print $string."\n";
     system ($string);
     sleep 5;
    }
close (IN);
exit;


sub WriteHeader{
  my $File = $_[0];
  my $Step = $_[1];
  my $Cores= $_[2];
  my $header="#!/bin/bash\n# Job $File $Step\n".
  "#\$ -pe smp $Cores\n".
  "#\$ -l h_vmem=16G\n
#\$ -l h_rt=1:0:0\n
#\$ -l node_type=nxv
#\$ -cwd
#\$ -j y
datapath=\"/path/to/lab/directory/genomes/ensembl\"\n
jobpath=\"/path/to/lab/directory/out\"\n
repopath=\"/path/to/lab/directory/repository\"\n";
  return $header;
}

sub FetchJobID{
my @jobids=(); my @tmp2;
open (DATA, "job_ids.txt") || die "Could not open job ids file!\n";
while(<DATA>){
      chomp ($_);
      #print $_, "\n";
      if (/job-ID/){next;}
      if (/----/){next;}
      @tmp2=split(" ",$_);
      push @jobids, $tmp2[0];
     }
close DATA;
return $jobids[scalar(@jobids)-1]; #Return last job id
}
   