#!/usr/bin/perl -w
# Program to generate jobfiles for read alignment of a batch of fastqfiles
# and send the jobfiles to the QMUL apocrita cluster


# Would need to be executed once for each matrix database
use strict;
use warnings;

# [PATHS]
my $homedir="/path/to/lab/directory/";
my $jobs_dir=$homedir."out/";
my $input_dir=$homedir."repository/";

# [PARAMETERS]
my $sqt = "pe";  # Paired-end data
#my $sqt = "se";  # Single-end data


# Read input data from file
my @tmp; my @jobs; my @tmp2;
my ($lane, $i, $line, $jobfile, $string, $jobname, $jobid, $filename);
open (IN, "AF_samples_ChIPSeq.txt") || die "Could not open input file!\n";
  while (<IN>) {
     chomp ($_);
     if (/# /){next;}
     if (/LibraryID/){next;}
     @tmp=split (/\t/,$_); # $tmp[1] LibraryID $tmp[3] SubmissionID $tmp[4] ExperimentID $tmp[5] Genome $tmp[6] Strain $tmp[7] Tissue $tmp[9] Individual $tmp[10] Sex $tmp[11] Antibody
     # Generate name for output files based on info provided in Inventory
     $filename = $tmp[1]."_".$tmp[11]."_".$tmp[7]."_".$tmp[5]."_".$tmp[4]."_".$tmp[3];
     print "$filename\n";
     # Concatenate fastq files for each library
     if ($sqt eq "pe"){
      $string = "cat $input_dir"."$tmp[1]/*R1*.fastq.gz > $input_dir"."$tmp[1]/".$filename."_R1.fastq.gz"; 
      print "$string\n";
      system ($string);
      $string = "cat $input_dir"."$tmp[1]/*R2*.fastq.gz > $input_dir"."$tmp[1]/".$filename."_R2.fastq.gz";
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
     print JOB & WriteHeader($tmp[1], "1"); # Print job file header via subroutine - same for all jobs


     # 1. BWA alignment
     print JOB ("#\$ -l h_vmem=10G","\n"); # 10GB RAM (maximum amount of memory you can specify in a job)
     print JOB ("module load bwa\n");
     if($sqt eq "se"){print JOB ("bwa aln $homedir"."genomes/ensembl/indexed/$tmp[5].bwaidx $input_dir"."$tmp[1]/$tmp[1]*.fq.gz > $tmp[1].bwa\n");}
     if($sqt eq "pe"){ # Run both alignments on the same job to avoid issue with next job being dependent on two job ids - can improve at a later point
      print JOB ("bwa aln $homedir"."genomes/ensembl/indexed/$tmp[5].bwaidx $input_dir"."$tmp[1]/".$filename."_R1.fastq.gz > $filename"."_R1.sai\n");
      print JOB ("bwa aln $homedir"."genomes/ensembl/indexed/$tmp[5].bwaidx $input_dir"."$tmp[1]/".$filename."_R2.fastq.gz > $filename"."_R2.sai\n");
     }
     close JOB;
     # Submit to cluster
     $string="qsub ".$jobname."_1";
     print $string."\n";  #This is printed to the screen so that you can see what the script is doing
     system ($string); #Sends the first job to the cluster
     # Store job id
     system ("qstat > job_ids.txt"); #Check what the Job ID is and store that in a text file

     # 2. Create sam file
     open (JOB, ">".$jobname."_2") || die "Could not open output job file!\n";
     print JOB & WriteHeader($tmp[1], "2"); # Print job file header via subroutine - same for all jobs
     print JOB ("module load bwa\n");
     if($sqt eq "se"){
       print JOB ("#\$ -l h_vmem=4G","\n"); # 4GB RAM
       print JOB ("bwa samse $homedir"."genomes/ensembl/indexed/$tmp[5].bwaidx $tmp[1].bwa $input_dir"."$tmp[1]/$tmp[1]*.fq.gz > $filename.sam\n");}
     if($sqt eq "pe"){
        print JOB ("#\$ -l h_vmem=10G","\n"); # 10GB RAM - increase memory for paired-end data, seems to fail on and off otherwise
        print JOB ("bwa sampe $homedir"."genomes/ensembl/indexed/$tmp[5].bwaidx $filename"."_R1.sai $filename"."_R2.sai $input_dir"."$tmp[1]/$filename"."_R1.fastq.gz $input_dir"."$tmp[1]/$filename"."_R2.fastq.gz > $filename.sam\n");}
     close JOB;
     # Submit to cluster as dependent on the previous job
     $jobid = &FetchJobID();
     $string="qsub -hold_jid $jobid $jobname"."_2"; #Submits the job to the cluster as a job that is dependent to the previous job (it is sent to the cluster but won't be run until the previous job has finished). It allows you to send the three steps of the alignment in one single time but making sure you don't start running a job that is dependent on a previous job.
     print $string."\n";
     system ($string);
     # Store job id
     system ("qstat > job_ids.txt");

     # 3. Create bam file, sort and index, then save mapping statistics to file
     open (JOB, ">".$jobname."_3") || die "Could not open output job file!\n";
     print JOB &WriteHeader($tmp[1], "3"); # Print job file header via subroutine - same for all jobs
     if($sqt eq "se"){print JOB ("#\$ -l h_vmem=4G","\n");} # 4GB RAM
     if($sqt eq "pe"){print JOB ("#\$ -l h_vmem=10G","\n");} # Increase memory usage for paired-end data; seems to fail on and off with 4G
     print JOB ("module load samtools\n");
     print JOB ("samtools view -bS $filename.sam > $filename.bam\n");
     print JOB ("samtools sort -O bam -o $filename"."_sorted.bam -T temp $filename.bam\n");
     print JOB ("samtools index $filename"."_sorted.bam\n");
     print JOB ("mv $filename"."_sorted.bam $filename.bam\n");
     print JOB ("samtools flagstat $filename.bam > $tmp[1]_mappingStats.txt\n");
     close JOB;

     # Submit to cluster as dependent on the previous job
     $jobid = &FetchJobID(); #JobID is the last job id in the queue
     $string="qsub -hold_jid $jobid $jobname"."_3"; #Submits the job to the cluster as a job that is dependent to the previous job (it is sent to the cluster but won't be run until the previous job has finished). It allows you to send the three steps of the alignment in one single time but making sure you don't start running a job that is dependent on a previous job.
     print $string."\n";
     system ($string);
     sleep 5;
   }
close (IN);
exit;

#WriteHeader is a small subroutine to write the header of each of the jobs (lines that must be specified according to the submission system, cluster documentation!)
sub WriteHeader{
  my $File = $_[0];
  my $Step = $_[1];
  my $header="#!/bin/sh\n# Job $File $Step\n".
  "#\$ -cwd\n".
  "#\$ -pe smp 5\n".
  "#\$ -l h_rt=24:0:0\n
  datapath=\"/path/to/lab/directory/genomes/ensembl\"\n
  jobpath=\"/path/to/lab/directory/out\"\n
  repopath=\"/path/to/lab/directory/repository\"\n";
  return $header;
}

#You can make some of the jobs wait for previous jobs
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
