###  Introduction

This document presents the analyses made to obtain the results that are presented hereafter. For repeatability or different lists the paths are to be changed and adapted.

The versions of the software used are the following:

Here only the modules used in the following scripts and not all modules in the

more module

```bash
#%Module1.0###############################################################

### RAJOUTER PLINK
### RAJOUTER ADMIXTURE
module load bioinfo/seqtk-1.2
module load bioinfo/bwa-0.7.15
module load bioinfo/samtools-1.8
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/bcftools-1.6
module load bioinfo/picard-2.18.2
module load bioinfo/gatk-4.1.2.0
```

### Reference genome

The reference genome used is Amel_HAv3.1 from Genebank: GCA_003254395.2

### fastq file

The individual population fastq files are obtained by Illumina sequencing.

### Preparation of the mapping

For each run we have R1_fastq.gz and R2_fastq.gz with the forward and reverse read

### Recovering the paths of fastq

to get the list :

```bash
ls /genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune* | awk 'BEGIN{FS=OFS="_"}{print $1,$2,$3}' | sort | uniq > corse_yellow.list
```

```bash
head -3 corse_yellow.list

/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune11_TCTCTACT-GAACCGCG-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune12_CTCTCGTC-AGGTTATA-AHJJN2DSX2_L003
```

corse_yellow.list divide into 7 to be able to launch the jobs little by little

### Mapping

Mapping scripts do the mapping with BWA, marking reads with Picard, GATK BQSR (Base Quality Score Replication)

Map_seqapipop_HAV3_1.bash

This mapping script calls 3 scripts in total:

- mappingAV_2019_Dec.sh
- bootstrapingAV_2019_Dec.sh
- callingAV_2019_Dec.sh

Editing the paths:

SAMPLE_FILE OUT

Call the variants with GATK HaplotypeCaller, which give the individual .gvcf files

```bash
#!/bin/bash

# Usage:	Map_seqapipop.bash
# Example:	Map_seqapipop.bash

N=1
PLOIDY=2

#------------------------------------------------------------------
#EDIT FOR ACCESS TO SEQUENCE FILES
# First download links from NG6 to genotoul and update the path here
SAMPLE_FILE=/home/agirardon/work/seqapipopOnHAV3_1/corse_yellow.list2
OUT=/home/agirardon/work/seqapipopOnHAV3_1/Map_2
REF=/home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/GCF_003254395.2_Amel_HAv3.1_genomic.fna
CHROMOSOME=/home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/HAv3_1_Chromosomes.list
UNKNOWN=/home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/HAv3_1_Unknown.list
SCRIPT=/genphyse/cytogen/seqapipop/ScriptsHAV3Mapping/mapping_calling
#------------------------------------------------------------------
mkdir -p ${OUT}/logs/mapping	# /work/ktabetaoul/Project_seqapipop/RESULTATS/logs/mapping
cd ${OUT}

#EDIT FOR mapping.sh
# This is the loop to send each sample to the cluster to be mapped
while read line
do
	ID=`basename ${line}`	# give the name of the sample-population to map
	IN=`dirname ${line}`	# give the path where find the sample fastq
	
	mkdir -p ${OUT}/${ID}
	
	sbatch --cpus-per-task=1 --mem-per-cpu=6G \
		-J ${ID}_mapping -o ${OUT}/logs/${ID}_mapping.o -e ${OUT}/logs/${ID}_mapping.e \
		${SCRIPT}/mappingAV_2019_Dec.sh -s ${ID} -i ${IN} -o ${OUT}/${ID} -p ${PLOIDY} -n ${N} -e ${OUT}/logs -R ${REF} -C ${CHROMOSOME} -U ${UNKNOWN}
		
done < ${SAMPLE_FILE}


# end of file
```



The PLOIDY is set to 2, although normally individuals are haploid, tests have been done with PLOIDY=1 which led to a false genotype in the repeated sequences, and gave false positive SNPs

N=1 as it can be used for a pool of individuals

More detailed explanations in the pdf : Mapping_Explanation_PLOIDY=2.pdf

This script takes a list of paths to sample names:

```bash
head -3 corse_yellow.list

/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune11_TCTCTACT-GAACCGCG-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune12_CTCTCGTC-AGGTTATA-AHJJN2DSX2_L003
```

This will result in several files:

The .bam by calling mappingAV_2019_Dec.sh it performs the mapping and detection of duplicate reads:

- BWA
- Picard MarkDuplicates

The BQSR, by calling bootstrapAV_2019_Dec.sh does a recalibration of the base quality score on each chromosome

- CORyellow10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003.bam
- CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003_bam.list
- CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003.bam.bai

and calls the calling, callingAV_2019_Dec.sh which uses the .bam to detect SNPs through the HaplotypeCaller

- CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003.g.vcf.gz
- CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003.g.vcf.gz.tbi

### Control gvcf are complete

A check is performed to verify that the jobs have been completed correctly (possible interruption, etc...)

```bash
#!/bin/bash

#Usage: controlVcfs.bash

module load -f  /home/gencel/vignal/save/000_ProgramModules/program_module

for i in `ls /home/agirardon/work/seqapipopOnHAV3_1/Map_2/CORjaune*/calling/*.gz`
do
       ID=`basename ${i}`      # give the name of the sample-population to map
       IN=`dirname ${i}`
       sbatch -J test --mem=1G \
       --wrap="module load bioinfo/bcftools-1.6; \
       bcftools query -f '%CHROM\n' ${i} | \ # bcftools query extrait un champ des VCF -- -f format de sortie -- %CHROM on veut la colone CHROM soit print -- \n nouvelle ligne 
       grep ^NC | uniq -c > \ # grep : recherche  un expression regulière commencant par NC -- uniq -c : Compter les occurrences des lignes en double
       /home/agirardon/work/seqapipopOnHAV3_1/controlling/controlVcfs/${ID}.count"
done
```

Output file: 

```bash
more CORjaune11_TCTCTACT-GAACCGCG-AHJJN2DSX2_L003.g.vcf.gz.count 
```

```
7168469 NC_037638.1
4189669 NC_037639.1
3582817 NC_037640.1
3473645 NC_037641.1
3630275 NC_037642.1
4552806 NC_037643.1
3643084 NC_037644.1
3332294 NC_037645.1
3218673 NC_037646.1
3021415 NC_037647.1
4146684 NC_037648.1
3007471 NC_037649.1
2969779 NC_037650.1
2795936 NC_037651.1
2487485 NC_037652.1
1898382 NC_037653.1
   1447 NC_001566.1
```

with the command

In fact, we check that we have reached the mitochondrial chromosome which is the last to be launched, as it is not useful for the analysis even if it is not complete, this is not a handicap for the rest of the analysis.

```bash
grep NC_001566.1 *.count | wc -l
```

### Recuperation list

recovery of paths for the next:

```bash
#!/bin/bash

ls /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/*/calling/*vcf.gz > chemin.all

for i in `awk '{print$1}' /work/genphyse/cytogen/Alain/seqapipopOnHAV3_1/seqApiPopVcfFilteredSonia/plinkAnalyses/WindowSNPs/RFMix/Pure95/in/IndsPopReference.list`
do
	grep ${i}_ chemin.all | awk '{print "    --variant",$1,"\\"}'  >> ref.list 
done 

rm chemin.all
```

this will give :

```bash
head -3 ref.list
```

```
 --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/Ab-PacBio_TCCGCGAA-CCTATCCT-BHCFFJDSXX_L004/calling/Ab-PacBio_TCCGCGAA-CCTATCCT-BHCFFJDSXX_L004.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/BER10_GTCCGC_L003/calling/BER10_GTCCGC_L003.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/BER11_GTGAAA_L003/calling/BER11_GTGAAA_L003.g.vcf.gz \
```

### Combine gvcf

In 3 parts:

- Head:

```bash
#!/bin/bash

#combineGVCFsHAV3Called_slurm.bash

module load -f /home/gencel/vignal/save/000_ProgramModules/program_module

while getopts":c:" opt; do
 case $opt in
c) CHR=${OPTARG};;
 esac
done

gatk --java-options "-Xmx80g" CombineGVCFs \
-R /home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
-L ${CHR} \
```

- List
- Tail:

```bash
-O /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/CORjaune/LesVCF/MetaGenotypes${CHR}.g.vcf.gz
echo "Finnished: "`date`
```

We thus obtain the Called script

This script allows to launch the different chromosomes in parallel, thanks to the passage of the variable ${i} with the option -c.

which will be launched by the combineGVCFsHAV3_1_Lance_slurm.bash script:

```bash
#!/bin/bash

#combineGVCFsHAV3_1_Lance_slurm.bash

for i in `cat /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/genotyping/HAv3_1_Chromosomes.list | cut -f 1`

do

mkdir -p /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/logs/

sbatch --cpus-per-task=1 --mem-per-cpu=100G \
        -J ${i}_genotAll \
        -o /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/logs${i}_combine.o \
        -e /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/logs/${i}_combine.e \
        /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/combineGVCFsHAV3_1_Called_slurm.bash -c ${i}
done
```

One file should be obtained per chromosome (17: 16 autosomes 1 mitochondrial).

### Check Combine

Checks that the jobs are completed (risk of interrupted jobs for various reasons..)

```bash
#!/bin/bash

#combineCheck.bash

module load -f /home/gencel/vignal/save/000_ProgramModules/program_module

for i in `ls /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/*.gz`
do
ID=`basename ${i}`
OUT=${ID/gz}check
echo ${ID}
echo ${ID/gz}check
echo ${i}
sbatch --cpus-per-task=1 --mem-per-cpu=1G \
       --wrap="module load -f /home/gencel/vignal/save/000_ProgramModules/program_module; \
       bcftools query -f '%CHROM\t%POS\n' ${i} | tail > \
       /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/checkCombine/${OUT}"
done
```

The output is the last 10 positions for each chromosome in MetaGenotypesNC_001566.1.g.vcf.check

Positions of the last variants of each chromosome:

```bash
 for i in `ls *.check`; do tail -1 ${i}; done
```

```
NC_001566.1     16343
NC_037638.1     27752957
NC_037639.1     16088822
NC_037640.1     13616014
NC_037641.1     13404415
NC_037642.1     13896578
NC_037643.1     17789067
NC_037644.1     14198576
NC_037645.1     12716919
NC_037646.1     12354637
NC_037647.1     12359965
NC_037648.1     16352576
NC_037649.1     11514127
NC_037650.1     11279598
NC_037651.1     10669577
NC_037652.1     9533787
NC_037653.1     7238523
```

We compare the length of the chromosomes and check that they are close: 

```bash
cat /home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/HAv3_1_Chromosomes.list
```

```
NC_037638.1     27754200
NC_037639.1     16089512
NC_037640.1     13619445
NC_037641.1     13404451
NC_037642.1     13896941
NC_037643.1     17789102
NC_037644.1     14198698
NC_037645.1     12717210
NC_037646.1     12354651
NC_037647.1     12360052
NC_037648.1     16352600
NC_037649.1     11514234
NC_037650.1     11279722
NC_037651.1     10670842
NC_037652.1     9534514
NC_037653.1     7238532
NC_001566.1     16343
```

### Genotyping

Genotyping script :

Basically the same as the previous script

```bash
#!/bin/bash

#genotypeGVCFsHAV3_1_Called_slurm.bash

module load -f /home/gencel/vignal/save/000_ProgramModules/program_module


while getopts ":c:" opt; do
  case $opt in
    c) CHR=${OPTARG};;
  esac
done

gatk --java-options "-Xmx80g" GenotypeGVCFs \
     -R /home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
     -V /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/MetaGenotypes${CHR}.g.vcf.gz \
     --use-new-qual-calculator \
     -O /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/MetaGenotypesCalled${CHR}.vcf.gz

echo "Finnished: "`date`

#end of file
```

Called by the script :

```bash
#!/bin/bash

#genotypeGVCFsHAV3_1_Lance_slurm.bash

for i in `cat /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/genotyping/HAv3_1_Chromosomes.list | cut -f 1`
do
mkdir -p /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/genotyping/logsCalled
sbatch --cpus-per-task=1 --mem-per-cpu=100G \
       -J ${i}_genotAll \
       -o /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/genotyping/logsCalled/${i}_genotype.o \
       -e /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/genotyping/logsCalled/${i}_genotype.e \
       /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/genotyping/genotypeGVCFsHAV3_1_Called_slurm.bash -c ${i}
done
```

### Checking

```bash
grep complete *_genotype.e
```

This line only appears if everything went well, so you should get one line per chromosome



### Concatenate vcf file

Here, a simple concatenation of the VCF file into 1 single VCF comprising the 17 chromosomes is performed

```bash
#!/bin/bash

#concatenateVCFs.bash

module load -f /home/gencel/vignal/save/000_ProgramModules/program_module

bcftools concat -f /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/vcf.list -o MetaGenotypesCalled870.vcf.gz -O z

tabix MetaGenotypesCalled870.vcf.gz
```

less -S MetaGenotypesCalled870.vcf.gz for better readability of the file :



### Count variant per chromosome

```bash
#!/bin/bash

#statsVcfsRaw.bash

zcat MetaGenotypesCalled870.vcf.gz | grep -v '#' | cut -f 1 | sort | uniq -c | \

awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > countVcfSumRaw
```

```
more countVcfSumRaw
```

```
NC_001566.1     1251
NC_037638.1     1475979
NC_037639.1     813108
NC_037640.1     756448
NC_037641.1     657904
NC_037642.1     739419
NC_037643.1     947867
NC_037644.1     733665
NC_037645.1     683901
NC_037646.1     618395
NC_037647.1     583767
NC_037648.1     753642
NC_037649.1     640188
NC_037650.1     608006
NC_037651.1     545183
NC_037652.1     495114
NC_037653.1     416511
Sum     11470348
```

About 11 million variants



### Keeping track of SNPs

```bash
#!/bin/bash

#retainSNPs.bash

module load -f /home/gencel/vignal/save/000_ProgramModules/program_module

OUT=/home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/

gatk --java-options "-Xmx64g" SelectVariants \
     -R /home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
     -V ${OUT}/MetaGenotypesCalled870.vcf.gz \
     --select-type-to-include SNP \
     -O ${OUT}/MetaGenotypesCalled870_raw_snps.vcf.gz

#Count the SNPs
zcat MetaGenotypesCalled870_raw_snps.vcf.gz | grep -v '#' | cut -f 1 | sort | uniq -c | \
awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > countVcfSumRawSNPs
```

```
more countVcfSumRawSNPs
```

```
NC_001566.1	906
NC_037638.1	1036622
NC_037639.1	574994
NC_037640.1	524172
NC_037641.1	463357
NC_037642.1	519761
NC_037643.1	668226
NC_037644.1	512649
NC_037645.1	475086
NC_037646.1	435269
NC_037647.1	415238
NC_037648.1	528650
NC_037649.1	451015
NC_037650.1	434113
NC_037651.1	379640
NC_037652.1	345615
NC_037653.1	288022
Sum	8053335
```

About 8 million SNPs 

### Hold the INDELs

```bash
#!/bin/bash

#retainSNPs.bash

module load -f /home/gencel/vignal/save/000_ProgramModules/program_module

OUT=/home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/Filtrage/retainIndels

gatk --java-options "-Xmx64g" SelectVariants \
     -R /home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
     -V ${OUT}/MetaGenotypesCalled870.vcf.gz \
     --select-type-to-include INDEL \
     -O ${OUT}/MetaGenotypesCalled870_raw_indels.vcf.gz

#Count the SNPs
zcat MetaGenotypesCalled870_raw_snps.vcf.gz | grep -v '#' | cut -f 1 | sort | uniq -c | \
awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > countVcfSumRawINDELs
```

```
more countVcfSumRawINDELs
```

```
NC_001566.1	247
NC_037638.1	374573
NC_037639.1	204149
NC_037640.1	195395
NC_037641.1	166756
NC_037642.1	186348
NC_037643.1	237894
NC_037644.1	187471
NC_037645.1	173642
NC_037646.1	155994
NC_037647.1	142659
NC_037648.1	190957
NC_037649.1	161150
NC_037650.1	147239
NC_037651.1	140678
NC_037652.1	126219
NC_037653.1	107365
Sum	2898736
```

About 2,900,000 Indels 


### SNP list retrieval 

- 2 VCFs with 2 different SNP lists have to be prepared

1 after all filters, with 7 million SNPs 

 - Retrieve the list of 7 million SNPs 

```bash
#!/bin/sh

module load -f /home/gencel/vignal/save/000_ProgramModules/program_module

bcftools query -f '%CHROM %POS\n' /work/genphyse/cytogen/Alain/seqapipopOnHAV3_1/seqApiPopVcfFilteredSonia/vcf_cleanup/MetaGenotypesCalled870_raw_snps_allfilter.vcf.gz -o /home/agirardon/work//seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/list7m.list
```

the file looks like : 

```bash
head -3 recup_list.bash
```

```
NC_037638.1 5671
NC_037638.1 5698
NC_037638.1 6621
```

- As for the list of 600,000 in the form of plink 

```bash
#!/bin/sh



awk '{print $2, $4}' /work/genphyse/cytogen/Alain/seqapipopOnHAV3_1/seqApiPopVcfFilteredSonia/plinkAnalyses/WindowSNPs/SeqApiPop_629_maf001_LD03_pruned.bim > /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/snplist600k.txt

```

the file looks like : 

```bash
head -3 snplist600k.txt
```

```
1:7577[HAV3.1]AG
1:10360[HAV3.1]CT
1:11791[HAV3.1]AG
```

# Filtering 

To perform the filter on our 8,053,335 SNPs, we will retrieve SNPs that have already been filtered by Sonia E. Eynard, located in :

/work/genphyse/cytogen/Alain/seqapipopOnHAV3_1/seqApiPopVcfFilteredSonia/vcf_cleanup/MetaGenotypesCalled870_raw_snps_allfilter.vcf.gz


In fact here we have submitted our VCF file to several quality filters to select the right SNPS:

The first set of filters addresses technical issues related to the sequencing and alignment steps and was therefore used for the total dataset of 870 samples, in order to benefit from its larger size for SNP detection and validation (Supplementary Figure 2). These filters included (i) strand bias and mapping quality metrics (SOR ≥ 3; FS ≤ 60 and MQ ≥ 40), (ii) genotyping quality metrics (QUAL > 200 and DQ < 20) and (iii) individual SNP genotyping metrics (heterozygous calls < 1%; missing genotypes < 5%, number of alleles < 4, and genotypes with individual HQ < 10 < 20%). Distribution and ECDF plots of the values for all filters used on the dataset were used to select the thresholds and are presented in the supplementary file SeqApiPop_2_VcfCleanup.pdf.

The figures can be found in the article (should I add them to my git hub): https://www.biorxiv.org/content/10.1101/2021.09.20.460798v2.full -> Quality filters on SNPs

It is obtained about 7 million good quality SNPs, on which we will base our SNPs filtering. We will match these SNPs, which we are sure of their quality, to the SNPs we have obtained thanks to the filterisec.sh script, and retain only the SNPs in common, this avoids redoing all the quality filtering step: 

/!\ ATTENTION TO THE DIRECTION FOR ISEC the first file will be the one "of support" for the intersection (to have only 403 individual and not 870 if one changes direction)

```bash
#!/bin/sh

module load -f  /home/gencel/vignal/save/000_ProgramModules/program_module


bcftools isec -c none -p /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/outisecautresens -n=2 -w1  /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/Filtrage/retainSNP/MetaGenotypesCalled870_raw_snps.vcf.gz /work/genphyse/cytogen/Alain/seqapipopOnHAV3_1/seqApiPopVcfFilteredSonia/vcf_cleanup/MetaGenotypesCalled870_raw_snps_allfilter.vcf.gz
```

-p creates the output directory with the vcf file and sites 

-c none is the default, so you don't really need to include it in the command, but it tells bcftools to consider two variants as identical only if their chromosome, pos, ref, and alt are all identical. 
Note: this means that A>G and A>G,C are NOT identical.

-n=2 tells bcftools isec to output the variants present in exactly two files.

-w1 Without this option, bcftools will print two files, one with the results of A overlapping with B, and another with the results of B overlapping with A. As these would be redundant, we use this option to print only one file.




Another method with the list of 7 million SNPs retrieved above 

ATTENTION YOU MUST SEPARATE IN THE LIST, WHEN RECEIVING THE SNPs THEY ARE SPACES !!!!

Here we do the same filtering but with bcftools view from the list of 7 million SNPs, thanks to the script filtreview.sh :

```bash
#!/bin/sh

module load -f  /home/gencel/vignal/save/000_ProgramModules/program_module

bcftools view -R list7mtab.list MetaGenotypesCalled870_raw_snps.vcf.gz > MetaGenotypesCalled870_raw_snps_filtreview.vcf.gz
```

(SEE WHICH METHOD IS BEST) -> isec seems to be faster in terms of calculation and handling (no list recuperation)

So we get the file MetaGenotypesCalled870_raw_snps_filtreisec.vcf.gz, where we have all our filtered SNPs 

And we count the number of SNPs we have after "filtering" thanks to a similar script when detecting SNPs: statsVcf_filtred.bash: 

```bash
#!/bin/bash

cat MetaGenotypesCalled870_raw_snps_filtreisec.vcf.gz | grep -v '#' | cut -f 1 | sort | uniq -c | \

awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > countVcfSumRaw_isec

/!\ Ici cat et pas zcat car pas compréssé, cela signifique qu'il faudra le compressé pour le prochain filtre /!\
```

This gives us 5,104,090 SNPs

```
more countVcfSumRaw_isec
```

```
NC_001566.1 218
NC_037638.1 671995
NC_037639.1 386653
NC_037640.1 323763
NC_037641.1 307878
NC_037642.1 336904
NC_037643.1 432453
NC_037644.1 319024
NC_037645.1 285602
NC_037646.1 276177
NC_037647.1 255285
NC_037648.1 326554
NC_037649.1 288652
NC_037650.1 269178
NC_037651.1 242986
NC_037652.1 214940
NC_037653.1 165828
Sum 5104090
```



# Now we will apply the LD filter:

There are several steps to do this: 

- Editing the plink file : 

``` bash
head -3 snplist_plink_600k.txt
```

```
1:7577[HAV3.1]AG
1:10360[HAV3.1]CT
1:11791[HAV3.1]AG
```



But we want to have a CHROM POS file 

- 1) With nedit simply "search" -> "replace" -> character to modify -> to "".

for A, C, G, T, [HAV3.1] 

The file is now in the form : 

```
1:7577
```

a simple 

```bash
sed -i -e 's/:/\t/g' file
```

will produce a file of the form 

``` 
1 7577
1 10360
1 11791
```

This just makes it easier to cut because by default cut uses TAB as a separator

- 2) We now need to change the chromosome numbers to their identifiers:

```
ID Number
NC_001566.1 - mito (I assume delete in plink) 
NC_037638.1 - 1
NC_037639.1 - 2
NC_037640.1 - 3 
NC_037641.1 - 4 
NC_037642.1 - 5
NC_037643.1 - 6
NC_037644.1 - 7
NC_037645.1 - 8
NC_037646.1 - 9
NC_037647.1 - 10
NC_037648.1 - 11
NC_037649.1 - 12
NC_037650.1 - 13 
NC_037651.1 - 14
NC_037652.1 - 15 
NC_037653.1 - 16 
```

only the first 2 characters are retrieved because the chromosome numbers have 2 characters max and we delete 

```bash
cut -c2 snplist_plink_600k_modif2.txt > snplist_plink_600k_modif_col.txt
sed -i -e 's/\t//g' snplist_plink_600k_modif_col.txt
```

- 3) So we just have the chromosome numbers from 1 to 16 and with the following sed commands we replace the number with the identifier

```bash
sed -i -e 's/^1$/NC_037638.1/g' snplist_plink_600k_modif_col.txt.test
sed -i -e 's/^2$/NC_037639.1/g' snplist_plink_600k_modif_col.txt.test
sed -i -e 's/^3$/NC_037640.1/g' snplist_plink_600k_modif_col.txt.test
sed -i -e 's/^4$/NC_037641.1/g' snplist_plink_600k_modif_col.txt.test
```

- 4) We keep the second column in another file ( POS )

```bash
cut -f 2 snplist_plink_600k_modif2.txt > snplist_plink_600k_modifoui.txt
```


- 5) Then simply assemble the two files: 

```bash
paste snplist_plink_600k_modif_col.txt.test snplist_plink_600k_modifoui.txt > snp_list_plink_600k_fini.txt
```

So we get our plink in the desired format:

``` bash
head -3 snp_list_plink_600k_fini.txt
```

```
NC_037638.1 7577
NC_037638.1 10360
NC_037638.1 11791
```



BE IN OUTISEC OTHER WAY AND NOT OUTISEC

```bash
bcftools view -I /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/outisec/0000.vcf -O z -o MetaGenotypesCalled870_raw_snps_filter_isec.vcf.gz 
```

- We can now apply this list as a filter based on the same idea as the LD filters that were applied earlier, then with the filter_list_plink script:

``` bash
#!/bin/sh

module load -f /home/gencel/vignal/save/000_ProgramModules/program_module


bcftools view -I /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/outisec/0000.vcf -O z -o MetaGenotypesCalled403_raw_snps_filter_isec.vcf.gz 
bcftools index /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/outisec/MetaGenotypesCalled403_raw_snps_filter_isec.vcf.gz

bcftools view -R /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCFs/Concatenate/snplist_plink_600k_fini.txt /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCFs/Concatenate/outisec/MetaGenotypesCalled403_raw_snps_filtre_isec. vcf.gz -O z -o /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCFs/Concatenate/outisec/MetaGenotypesCalled403_raw_snps_isec_filter.vcf.gz


```

We now count the number of good quality SNPs after the LD filter with the script statsVcf_isec_sens1.bash :

``` bash
#!/bin/bash

#statsVcfsRaw.bash

zcat MetaGenotypesCalled403_raw_snps_isec_filter_plink.vcf.gz| grep -v '#' | cut -f 1 | sort | uniq -c | \

awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > countVcfSumRaw_isec_sens
```


This gives us 590,724 good quality SNPs

```bash
more countVcfSumRaw_isec_sens
```

```
NC_037638.1 81264
NC_037639.1 46384
NC_037640.1 40374
NC_037641.1 36824
NC_037642.1 38745
NC_037643.1 47025
NC_037644.1 34211
NC_037645.1 31633
NC_037646.1 32432
NC_037647.1 30742
NC_037648.1 33379
NC_037649.1 34858
NC_037650.1 31729
NC_037651.1 28807
NC_037652.1 24934
NC_037653.1 17383
Sum 590724
```

### PCA with Plink

The PCA will be performed with Plink, and the following script:

```bash
#!/bin/bash

module load -f /home/agirardon/work/seqapipopOnHAV3_1/program_module/module


# perform linkage pruning 

VCF=/home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCFs/Concatenate/outisecautresens/MetaGenotypesCalled403_raw_snps_filtre_isec_plink.vcf.gz

plink --vcf ${VCF} --allow-extra-chr \ #allow-extra-chr allows to read the file despite chromosome names and not chromosome numbers
      --indep-pairwise 1749 175 0.3 --out SeqApiPop_403_LD03


# prune and create pca

plink --vcf ${VCF} --allow-extra-chr
      --extract SeqApiPop_403_LD03.prune.in
      --make-bed --pca --out PCA_SeqApiPop_403_LD03
```

In the output we have several files, including the one we will use to plot the PCA: PCA_SeqApiPop_403_LD03.eigenvec

```bash
head -3 PCA_SeqApiPop_403_LD03.eigenvec
```

```
AOC10 AOC10 0.0214301 0.0073294 0.00899092 0.0481961 0.0326964 0.0345816 -0.000382077 -0.0110416 -0.0438832 -0. 000245611 -8.02539e-05 0.00158182 -0.00117293 -0.00122113 0.00110382 0.00133751 -0.000410789 0.00212435 0.0024354 0.00162384
AOC11 AOC11 0.0212416 0.00792505 0.0114335 0.0454902 0.0292201 0.0295928 -0.00482164 -0.00768436 -0.0317105 -0.0011311 0. 000944703 0.00290791 0.00036738 0.00201229 0.000597761 -8.23973e-06 0.00115585 -0.00245061 0.000437679 0.000697422
AOC12 AOC12 0.0203242 0.0084029 0.0106885 0.0467117 0.0326449 0.0326058 0.00100988 -0.00928108 -0.0455404 -0.000493906 -0. 00158167 -0.000253945 -0.000226169 0.000236364 0.000728524 -0.000162989 0.000967244 -0.00127729 0.000300538 -0.000624591
```



However for the rest of the analysis it is important to sort it, because we will have to add to this file, the name of the species to plot the PCA, indeed our list looks more like this :

``` 
ID Sp
Ab.PacBio Mellifera							
AOC10 AncientCorsica
AOC11 Ancient Corsica
AOC12 Ancient Corsica
AOC14 Ancient Corsica
AOC15 Ancient Corsica
AOC16 CorseAnciennes
AOC17 CorseAnciennes
AOC18 CorseAnciennes
AOC19 CorseAnciennes
AOC20 CorseAnciennes
```

Here indeed our first individual is Ab.PacBio, while in the .eigenvec output the first individual is AOC10

So with a simple : 

```bash
sort PCA_SeqApiPop_403_LD03.eigenvec > sort_PCA_SeqApiPop_403_LD03.eigenvec
```

So we have the sort file for the rest of the analysis in the .rdm file PCA_seqapipop

### Additional Analysis :

# MAF Filter 0.01

We perform a MAF 0.01 filter with plink on our filtered data with the 7 millions thanks to the filter: FiltreMaf001.sh

First, we need to prepare a list with our Unique403_test.list individuals:

```bash
head -3 Unique403_test.list 
```

```
Ab-PacBio Ab-PacBio	
BER10 BER10	
BER11 BER11	
```



We can now perform a Maf 0.01 filter on our vcf file with Plink using the following script: 

(here we are doing a filter with maf 0.05 but I won't use it just in case I need it later) 

```bash
#!/bin/bash

module load -f /home/agirardon/work/seqapipopOnHAV3_1/program_module/module


NAME=SeqApiPop_403

VCFin=/home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/outisecautresens/MetaGenotypesCalled403_raw_snps_filtre_isec.vcf.gz
VCFout=~/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/Filtering/Maf001/${NAME}



plink --vcf ${VCFin} \
      --keep-allele-order
      --keep ~/work/seqapipopOnHAV3_1/combineGVCFs/LesVCFs/Concatenate/Filtering/Maf001/Unique403_test.list
      --allow-no-sex \
      --allow-extra-chr \
      --set-missing-var-ids @:#[HAV3.1]\$1\$2 \
      --mind 0.1
      --geno 0.1 \
      --out ${VCFout} \
      --make-bed
      --missing
      
plink --bfile ${VCFout} \
      --allow-extra-chr
      --maf 0.01
      --out ${VCFout}_maf001
      --make-bed
      
plink --bfile ${VCFout} \
      --allow-extra-chr \
      --out ${VCFout}_maf005
      --maf 0.05 \
      --make-bed

```

And we get the .bim file: 

```bash
head -3 SeqApiPop_403_maf001.bim
```

```
NC_037638.1 NC_037638.1:7034[HAV3.1]AG 0 7034 A G
NC_037638.1 NC_037638.1:7092[HAV3.1]AG 0 7092 G A
NC_037638.1 NC_037638.1:7299[HAV3.1]CT 0 7299 T C
```



To get the file in the desired format to be able to filter on my VCF, I will make this file a simple list containing the chromosome and the SNP position: 

with simple :

``` bash
cut -f 1 SeqApiPop_403_maf001.bim > list.Maf001_col1
cut -f 4 SeqApiPop_403_maf001.bim > list.Maf001_col4

paste list.Maf001_col1 list.Maf001_col4 > snp_list_filtremaf001.txt
```

I get the desired format

```bash
head -3 snp_list_filtremaf001.txt
```

```
NC_037638.1 7034
NC_037638.1 7092
NC_037638.1 7299
```

Here we have the individuals with mitochondria that are not useful for the analysis so I will simply delete them with nedit and have the desired list in snp_list_filtremaf001_nomito.txt