#!/bin/bash
#$ -cwd           # Set the working directory for the job to the current directory
#$ -j y           # Join stdout and stderr
#$ -pe smp 6       # Request 1 CPU core
#$ -l h_rt=24:0:0  # Request 1 hour runtime
#$ -l h_vmem=20G   # Request 1GB RAM / core, i.e. 1GB total


cd /path/to/project/peak_calling_for_pairwise_comparisons


#load python and macs2 (via virtual environment)
#module load python/3.6.3
module load python
#virtualenv --include-lib macs2
virtualenv macs2
source macs2/bin/activate
pip install numpy
pip install macs2


# The effective sequencing depth is always the smaller number of treatment and control - take note manually

egrep "fragments after filtering in treatment|fragments after filtering in control" 25_AF_RA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 2231_AF_RA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 31_AF_RA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 313_AF_RA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 18_AF_RA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 28_AF_RA_peaks.xls

egrep "fragments after filtering in treatment|fragments after filtering in control" 12b_SR_RA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 27_SR_RA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 33_SR_RA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 20_SR_RA_peaks.xls

egrep "fragments after filtering in treatment|fragments after filtering in control" 25_AF_LA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 28_AF_LA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 47_AF_LA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 31_AF_LA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 2231_AF_LA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 313_AF_LA_peaks.xls

egrep "fragments after filtering in treatment|fragments after filtering in control" 12b_SR_LA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 45_SR_LA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 20_SR_LA_peaks.xls
egrep "fragments after filtering in treatment|fragments after filtering in control" 27_SR_LA_peaks.xls

#25_AF_RA_peaks treatment 4998342
#2231_AF_RA_peaks treatment 9650792
#31_AF_RA_peaks treatment 7011242
#313_AF_RA_peaks control 5125019
#18_AF_RA_peaks treatment 6696145
#28_AF_RA_peaks treatment 5210965

#12b_SR_RA_peaks control 6600510
#27_SR_RA_peaks treatment 5817565
#33_SR_RA_peaks treatment 7708321 
#20_SR_RA_peaks treatment 7641014

#25_AF_LA_peaks treatment 5845266
#28_AF_LA_peaks treatment 6876286
#47_AF_LA_peaks control 9122754
#31_AF_LA_peaks treatment 7403502
#2231_AF_LA_peaks control 8786777
#313_AF_LA_peaks control 6061167

#12b_SR_LA_peaks treatment 3349426
#45_SR_LA_peaks treatment 8782008
#20_SR_LA_peaks control 5002508
#27_SR_LA_peaks treatment 8655863


#####

declare -A arr

arr["2231_AF_RA"]=9650792

arr+=( ["31_AF_RA"]=7011242 ["313_AF_RA"]=5125019 ["18_AF_RA"]=6696145 ["28_AF_RA"]=5210965 ["12b_SR_RA"]=6600510 ["27_SR_RA"]=5817565 ["33_SR_RA"]=7708321 ["20_SR_RA"]=7641014 ["25_AF_LA"]=5845266 ["28_AF_LA"]=6876286 ["47_AF_LA"]=9122754 ["31_AF_LA"]=7403502 ["2231_AF_LA"]=8786777 ["313_AF_LA"]=6061167 ["18_AF_LA"]=3131291 ["12b_SR_LA"]=3349426 ["45_SR_LA"]=8782008 ["20_SR_LA"]=5002508 ["27_SR_LA"]=8655863 )



######

#25_AF_RA
for i in 2231_AF_RA 31_AF_RA 313_AF_RA 18_AF_RA 28_AF_RA 12b_SR_RA 27_SR_RA 33_SR_RA 20_SR_RA 25_AF_LA 28_AF_LA 47_AF_LA 31_AF_LA 2231_AF_LA 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 25_AF_RA_treat_pileup.bdg --c1 25_AF_RA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 4998342 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_25_AF_RA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done


#2231_AF_RA
for i in 31_AF_RA 313_AF_RA 18_AF_RA 28_AF_RA 12b_SR_RA 27_SR_RA 33_SR_RA 20_SR_RA 25_AF_LA 28_AF_LA 47_AF_LA 31_AF_LA 2231_AF_LA 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 2231_AF_RA_treat_pileup.bdg --c1 2231_AF_RA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 9650792 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_2231_AF_RA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done


#31_AF_RA
for i in  313_AF_RA 18_AF_RA 28_AF_RA 12b_SR_RA 27_SR_RA 33_SR_RA 20_SR_RA 25_AF_LA 28_AF_LA 47_AF_LA 31_AF_LA 2231_AF_LA 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 31_AF_RA_treat_pileup.bdg --c1 31_AF_RA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 7011242 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_31_AF_RA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done


#313_AF_RA
for i in  18_AF_RA 28_AF_RA 12b_SR_RA 27_SR_RA 33_SR_RA 20_SR_RA 25_AF_LA 28_AF_LA 47_AF_LA 31_AF_LA 2231_AF_LA 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 313_AF_RA_treat_pileup.bdg --c1 313_AF_RA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 5125019 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_313_AF_RA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done


#18_AF_RA
for i in 28_AF_RA 12b_SR_RA 27_SR_RA 33_SR_RA 20_SR_RA 25_AF_LA 28_AF_LA 47_AF_LA 31_AF_LA 2231_AF_LA 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 18_AF_RA_treat_pileup.bdg --c1 18_AF_RA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 6696145 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_18_AF_RA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done



#28_AF_RA
for i in 12b_SR_RA 27_SR_RA 33_SR_RA 20_SR_RA 25_AF_LA 28_AF_LA 47_AF_LA 31_AF_LA 2231_AF_LA 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 28_AF_RA_treat_pileup.bdg --c1 28_AF_RA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 5210965 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_28_AF_RA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done



#12b_SR_RA
for i in 27_SR_RA 33_SR_RA 20_SR_RA 25_AF_LA 28_AF_LA 47_AF_LA 31_AF_LA 2231_AF_LA 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 12b_SR_RA_treat_pileup.bdg --c1 12b_SR_RA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 6600510 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_12b_SR_RA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done


#27_SR_RA
for i in 33_SR_RA 20_SR_RA 25_AF_LA 28_AF_LA 47_AF_LA 31_AF_LA 2231_AF_LA 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 27_SR_RA_treat_pileup.bdg --c1 27_SR_RA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 5817565 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_27_SR_RA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done



#33_SR_RA
for i in 20_SR_RA 25_AF_LA 28_AF_LA 47_AF_LA 31_AF_LA 2231_AF_LA 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 33_SR_RA_treat_pileup.bdg --c1 33_SR_RA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 7708321 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_33_SR_RA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done



#20_SR_RA
for i in 25_AF_LA 28_AF_LA 47_AF_LA 31_AF_LA 2231_AF_LA 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 20_SR_RA_treat_pileup.bdg --c1 20_SR_RA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 7641014 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_20_SR_RA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done



#25_AF_LA
for i in 28_AF_LA 47_AF_LA 31_AF_LA 2231_AF_LA 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 25_AF_LA_treat_pileup.bdg --c1 25_AF_LA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 5845266 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_25_AF_LA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done



#28_AF_LA
for i in 47_AF_LA 31_AF_LA 2231_AF_LA 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 28_AF_LA_treat_pileup.bdg --c1 28_AF_LA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 6876286 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_28_AF_LA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done



#47_AF_LA
for i in 31_AF_LA 2231_AF_LA 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 47_AF_LA_treat_pileup.bdg --c1 47_AF_LA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 9122754 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_47_AF_LA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done


#31_AF_LA
for i in 2231_AF_LA 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 31_AF_LA_treat_pileup.bdg --c1 31_AF_LA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 7403502 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_31_AF_LA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done



#2231_AF_LA
for i in 313_AF_LA 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 2231_AF_LA_treat_pileup.bdg --c1 2231_AF_LA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 8786777 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_2231_AF_LA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done



#313_AF_LA
for i in 18_AF_LA 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 313_AF_LA_treat_pileup.bdg --c1 313_AF_LA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 6061167 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_313_AF_LA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done



#18_AF_LA
for i in 12b_SR_LA 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 18_AF_LA_treat_pileup.bdg --c1 18_AF_LA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 3131291 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_18_AF_LA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done



#12b_SR_LA
for i in 45_SR_LA 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 12b_SR_LA_treat_pileup.bdg --c1 12b_SR_LA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 3349426 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_12b_SR_LA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done



#45_SR_LA
for i in 20_SR_LA 27_SR_LA
    do
    macs2 bdgdiff --t1 45_SR_LA_treat_pileup.bdg --c1 45_SR_LA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 8782008 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_45_SR_LA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done



#20_SR_LA
for i in 27_SR_LA
    do
    macs2 bdgdiff --t1 20_SR_LA_treat_pileup.bdg --c1 20_SR_LA_control_lambda.bdg -C 2 --t2 "${i}"_treat_pileup.bdg --c2 "${i}"_control_lambda.bdg --d1 5002508 --d2 ${arr[${i}]} -g 100 -l 300 --o-prefix diff_20_SR_LA_vs_"${i}" --outdir bdgdiff_pairwise_analysis_output
    done


## Subset peaksets in subgroups; move to corresponding directory


    #AF:
      ls *AF_*_vs_*_cond1* | wc -l #--> 138
      ls *_vs_*_AF_*_cond2* | wc -l #--> 90
      ls *AF_*_vs_*_AF_*_cond1* | wc -l #--> 66
      ls *AF_*_vs_*_AF_*_cond2* | wc -l #--> 66

      cp *AF_*_vs_*_cond1* AF_enriched_peaks
      cp *_vs_*_AF_*_cond2* AF_enriched_peaks
      cd AF_enriched_peaks
      rm *AF_*_vs_*_AF_*_cond1*
      rm *AF_*_vs_*_AF_*_cond2*

      # 228-66-66 = 96 files


    #SR:
      ls *SR_*_vs_*_cond1* | wc -l #--> 52
      ls *_vs_*_SR_*_cond2* | wc -l #--> 100
      ls *SR_*_vs_*_SR_*_cond1* | wc -l #--> 28
      ls *SR_*_vs_*_SR_*_cond2* | wc -l #--> 28

      cp *SR_*_vs_*_cond1* SR_enriched_peaks
      cp *_vs_*_SR_*_cond2* SR_enriched_peaks
      cd SR_enriched_peaks
      rm *SR_*_vs_*_SR_*_cond1*
      rm *SR_*_vs_*_SR_*_cond2*

      # 152-28-28 = 96 files


  
    #AF-RA: 
      ls *AF_RA_vs_*_cond1* | wc -l #--> 99
      ls *_vs_*_AF_RA_*cond2* | wc -l #--> 15 
      ls *AF_RA_vs_*_AF_RA_*_cond1* | wc -l #--> 15
      ls *AF_RA_vs_*_AF_RA_*_cond2* | wc -l #--> 15

      cp *AF_RA_vs_*_cond1* AF_RA_enriched_peaks
      cp *_vs_*_AF_RA_*cond2* AF_RA_enriched_peaks
      cd AF_RA_enriched_peaks
      rm *AF_RA_vs_*_AF_RA_*_cond1*
      rm *AF_RA_vs_*_AF_RA_*_cond2*

    # 99 + 15 - 15 -15  = 84 files

    #AF-LA: 
      ls *AF_LA_vs_*_cond1* | wc -l #--> 39
      ls *_vs_*_AF_LA_*cond2* | wc -l #--> 75
      ls *AF_LA_vs_*_AF_LA_*_cond1* | wc -l #-> 15
      ls *AF_LA_vs_*_AF_LA_*_cond2* | wc -l #--> 15

      cp *AF_LA_vs_*_cond1* AF_LA_enriched_peaks
      cp *_vs_*_AF_LA_*cond2* AF_LA_enriched_peaks
      cd AF_LA_enriched_peaks
      rm *AF_LA_vs_*_AF_LA_*_cond1*
      rm *AF_LA_vs_*_AF_LA_*_cond2*

    # 39 + 75 - 15 -15 = 84 files

    #SR-RA: 
      ls *_vs_*_SR_RA_*cond2* | wc -l #--> 30
      ls *SR_RA_vs_*_SR_RA_*_cond1* | wc -l #--> 6
      ls *SR_RA_vs_*_SR_RA_*_cond2* | wc -l #--> 6

      cp *SR_RA_vs_*_cond1* SR_RA_enriched_peaks
      cp *_vs_*_SR_RA_*cond2* SR_RA_enriched_peaks
      cd SR_RA_enriched_peaks
      rm *SR_RA_vs_*_SR_RA_*_cond1*
      rm *SR_RA_vs_*_SR_RA_*_cond2*

    # 46 + 30 - 6 - 6 = 64 files

    #SR-LA: 
      ls *SR_LA_vs_*_cond1* | wc -l #--> 6
      ls *_vs_*_SR_LA_*cond2* | wc -l #--> 70
      ls *SR_LA_vs_*_SR_LA_*_cond1* | wc -l #--> 6
      ls *SR_LA_vs_*_SR_LA_*_cond2* | wc -l #--> 6

      cp *SR_LA_vs_*_cond1* SR_LA_enriched_peaks
      cp *_vs_*_SR_LA_*cond2* SR_LA_enriched_peaks
      cd SR_LA_enriched_peaks
      rm *SR_LA_vs_*_SR_LA_*_cond1*
      rm *SR_LA_vs_*_SR_LA_*_cond2*

    # 6 + 70 - 6 - 6 = 64 files




