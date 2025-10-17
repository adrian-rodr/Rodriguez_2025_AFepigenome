#!/usr/bin/perl -w
# Script to pre-process BAM files before running MACS2 peak calling
# and send the jobfiles to the QMUL apocrita cluster

use strict;
use warnings;

# [PATHS]
my $homedir="/path/to/lab/directory/";
my $jobs_dir=$homedir."out/";
my $input_dir=$homedir."repository/";
my $results_dir=$jobs_dir."peaks/";

# Read input data from file
my @tmp; my @jobs; my @tmp2; my @all=();
my ($lane, $i, $line, $jobfile, $string, $jobname, $jobid, $filename, $done);
$i=0;#Line counter
open (IN, "AF_samples_for_PeakCalling.txt") || die "Could not open input file!\n";
  while (<IN>) {
     chomp ($_);
     $i++;
     $done = 0;
     if (/GenomeSize/){next;}
     @tmp=split (/\t/,$_); # $tmp[1] Genome $tmp[2] GenomeSize $tmp[3] ChIPLibID $tmp[4] ChIPExpID $tmp[5] Factor $tmp[6] InputLibID $tmp[7] InputExpID $tmp[8] ExperimentID
     # Generate name for output files based on info provided in Inventory
     $filename = $tmp[3]."_".$tmp[6]."_".$tmp[1]."_".$tmp[5]."_".$tmp[8];
     # $filename = SampleID_InputLibID_Genome_Factor_ExperimentID
     # $filename = vl0211_vl0227_Homo_sapiens.GRCh38_H3K27ac_RightAtrium_47
     print "$filename\n";
     
     # Generate and submit job files for peak calling tasks
     
     # 1. First cluster job - Remove multi-mapping reads from [ChIP Library] and [InputLib] (if new)
     $jobname = $tmp[3]."_".$tmp[6]."_jobfile";
     
     open (JOB, ">".$jobname."_1") || die "Could not open output job file!\n";
     print JOB &WriteHeader($tmp[3], "1", "1"); # Print job file header via subroutine - specifies ChIP file, step and number of cores to use
     print JOB ("#\$ -l h_vmem=4G","\n"); # 4GB RAM memory
     print JOB ("module load samtools\n");
     print JOB ("samtools view -bq 1 ../repository/".$tmp[4]."/".$tmp[3]."*".$tmp[1]."*.bam > $tmp[3].bam\n"); #Remove multi-mapping reads
     print JOB ("samtools flagstat ".$tmp[3].".bam > ".$tmp[3]."_reads.txt\n");#Store number of reads
     print JOB ("samtools view -bq 1 ../repository/".$tmp[7]."/".$tmp[6]."*".$tmp[1]."*.bam > $tmp[6].bam\n"); #Remove multi-mapping reads
     print JOB ("samtools flagstat ".$tmp[6].".bam > ".$tmp[6]."_reads.txt\n");#Store number of reads 
     close JOB;
     # Submit to cluster
     $string="qsub ".$jobname."_1";
     print $string."\n";
     system ($string);
     # Store job id
     system ("qstat > job_ids.txt");
     
     # 2. Subsample to 20 million reads [ChIP Library] and [InputLib]
     open (JOB, ">".$jobname."_2") || die "Could not open output job file!\n";
     print JOB &WriteHeader($tmp[3], "2","1"); # Print job file header via subroutine - specifies ChIP file, step and number of cores to use
     print JOB ("#\$ -l h_vmem=4G","\n"); # 4GB RAM
     print JOB ("ReadsChIP=\$(sed -n 's/ + 0 mapped.*\$//p' $tmp[3]_reads.txt)\n"); #Store uniquely mapped reads in ChIP bam
     print JOB ("ReadsInput=\$(sed -n 's/ + 0 mapped.*\$//p' $tmp[6]_reads.txt)\n"); #Store uniquely mapped reads in input bam
     
     # Calculate fraction of reads to use in subsampling [within job script]
     print JOB ("if [ \$ReadsChIP -lt 20000000 ]\n");
     print JOB ("then\n");
     print JOB (" fChIP=0.99999\n echo \$fChIP\n"); #Use all reads if ner of uniquely mapped is lower than 20 million
     print JOB ("else\n");
     print JOB (" fChIP=\$(printf '%.5f\\n' \$(echo \"20000000/\${ReadsChIP}\" \| bc -l))\n echo \$fChIP\n");
     print JOB ("fi\n");
     
     print JOB ("if [ \$ReadsInput -lt 20000000 ]\n");
     print JOB ("then\n");
     print JOB (" fInput=0.99999\n echo \$fInput\n"); #Use all reads if ner of uniquely mapped is lower than 20 million
     print JOB ("else\n");
     print JOB (" fInput=\$(printf '%.5f\\n' \$(echo \"20000000/\${ReadsInput}\" \| bc -l))\n echo \$fInput\n");
     print JOB ("fi\n");
     
     # Subsample reads to 20 million
     print JOB ("module load samtools\n");
     print JOB ("samtools view -s \$fChIP $tmp[3].bam -b -o $tmp[3]_20M.bam\n");
     print JOB ("samtools view -s \$fInput $tmp[6].bam -b -o $tmp[6]_20M.bam\n");
     #print JOB ("rm $tmp[3].bam\n"); # Erase previous bam without sub-sampling
     #print JOB ("rm $tmp[5].bam\n");
     close JOB;
     
     # Submit to cluster as dependent on the previous job
     $jobid = &FetchJobID();
     $string="qsub -hold_jid $jobid $jobname"."_2";
     print $string."\n";
     system ($string);
     # Store job id
     system ("qstat > job_ids.txt");

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
  my $Cores =$_[2];
  my $header="#!/bin/sh\n# Job $File $Step\n".
  "#\$ -cwd\n".
  "#\$ -pe smp $Cores\n".
  "#\$ -l h_rt=24:0:0\n
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
     
