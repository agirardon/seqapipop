###  1. Introduction

Ce document présente les analyses faites pour obtenir les résultats qui sont présentés par la suite. Pour la répétabilité ou les différentes listes les chemins sont à changer et à adapter.

Les version des software utilisés sont les suivantes:

``` 
more module
#%Module1.0###############################################################

module load bioinfo/EIG-7.2.1
module load bioinfo/structure_v2.3.4_console
module load bioinfo/sNMF_CL_v1.2
module load bioinfo/CLUMPAK-v1.1
module load bioinfo/admixture_linux-1.3.0
module load bioinfo/blast-2.2.26
module load bioinfo/ncbi-blast-2.7.1+
module load bioinfo/last-956
module load bioinfo/ncbi_tools-v6.1
module load bioinfo/blatSuite.36
module load bioinfo/seqtk-1.2
module load bioinfo/bwa-0.7.15
module load bioinfo/samtools-1.8
module load bioinfo/bedtools-2.27.1
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/bcftools-1.6
module load bioinfo/picard-2.18.2
module load bioinfo/gatk-4.1.2.0
module load bioinfo/plink-v1.90b5.3
module load system/Python-3.6.3
module load system/Java8
module load system/R-3.6.2
#####For installation of some R packages
#module load compiler/gcc-7.2.0 libraries/gdal-2.4.0_gcc-7.2.0 libraries/proj-5.2.0_gcc-7.2.0 libraries/geos-3.4.2_gcc-7.2.0
#####
module load bioinfo/trf-v4.09
module load bioinfo/mafft-7.313
module load bioinfo/popoolation2_1201
```



### Genome de référence

Le génome de référence utiliser est Amel_HAv3.1 de Genebank : GCA_003254395.2



### Recuperation des chemins de fastq 

more corse_yellow.list diviser en 7 pour pouvoir lancer les jobs petits à petit 

### Mapping 



Les script de Mapping effectue le mapping avec BWA, marquage des lectures avec Picard, GATK BQSR( Base Quality Score Replication) 

On appelle les variants avec GATK HaplotypeCaller, qui donnent les fichiers .gvcf individuel 



```
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



La PLOIDY  est fixé à 2, bien que normalement les individus sont haploïdes, des tests ont été fait avec PLOIDY=1 cela conduisait a un faux genotype dans les sequences répétées, et donnaient des SNP faussement positifs 

N=1  car il peut être utilisé pour un pool d'individus



Ce script prend une liste de chemin vers des noms d'echantillons : 



```
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune11_TCTCTACT-GAACCGCG-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune12_CTCTCGTC-AGGTTATA-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune13_CCAAGTCT-TCATCCTT-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune14_TTGGACTC-CTGCTTCC-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune15_GGCTTAAG-GGTCACGA-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune16_AATCCGGA-AACTGTAG-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune17_TAATACAG-GTGAATAT-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune18_CGGCGTGA-ACAGGCGC-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune19_ATGTAAGT-CATAGAGT-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune1_CCGCGGTT-CTAGCGCT-AHJJN2DSX2_L003
```



On obtiendra donc plusieurs fichier :



Les .bam en appelant mappingAV_2019_Dec.sh il effectue le mapping et la detection des lectures dupliquées :



Le bootstraping en appelant bootstrapingAV_2019_Dec.sh fait un recalibrage du score de qualité de base sur chaque chromosome 

- CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003.bam

- CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003_bam.list
- CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003.bam.bai



et appelle le calling, callingAV_2019_Dec.sh qui fusionne les bams des chromosome dans un genome vcf et effectue le genotypage 

- CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003.g.vcf.gz
- CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003.g.vcf.gz.tbi



### Control gvcf are complete 



```
#!/bin/bash

#Usage: controlVcfs.bash

module load -f  /home/gencel/vignal/save/000_ProgramModules/program_module

for i in `ls /home/agirardon/work/seqapipopOnHAV3_1/Map_2/CORjaune*/calling/*.gz`
do
       ID=`basename ${i}`      # give the name of the sample-population to map
       IN=`dirname ${i}`
       sbatch -J test --mem=1G \
       --wrap="module load bioinfo/bcftools-1.6; \
       bcftools query -f '%CHROM\n' ${i} | \
       grep ^NC | uniq -c > \
       /home/agirardon/work/seqapipopOnHAV3_1/controlling/controlVcfs/${ID}.count"
done
```

En effet on verifie qu'on est bien arrivé jusqu'au chromosome mitochondrial qui est le dernier à être lancé 



avec la commande 

``` grep NC_001566.1 *.count | wc -l
 grep NC_001566.1 *.count | wc -l
```



### Recuperation list 



récupération des chemins pour la suite: 

```
#!/bin/bash

ls /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/*/calling/*vcf.gz > chemin.all

for i in `awk '{print$1}' /work/genphyse/cytogen/Alain/seqapipopOnHAV3_1/seqApiPopVcfFilteredSonia/plinkAnalyses/WindowSNPs/RFMix/Pure95/in/IndsPopReference.list`
do
	grep ${i}_ chemin.all | awk '{print "    --variant",$1,"\\"}'  >> ref.list #ATTENTION a supprimer si on relancer sinon ca ajoute tout
done 

rm chemin.all
```



ca donnera : 



``` 
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC10_ATGTCA_L002/calling/AOC10_ATGTCA_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC11_CCGTCC_L002/calling/AOC11_CCGTCC_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC12_GTCCGC_L002/calling/AOC12_GTCCGC_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC14_GTGAAA_L002/calling/AOC14_GTGAAA_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC15_GTGGCC_L002/calling/AOC15_GTGGCC_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC16_GTTTCG_L002/calling/AOC16_GTTTCG_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC17_CGTACG_L002/calling/AOC17_CGTACG_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC18_GAGTGG_L002/calling/AOC18_GAGTGG_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC19_ACTGAT_L002/calling/AOC19_ACTGAT_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC2_CAGATC_L002/calling/AOC2_CAGATC_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC20_ATTCCT_L002/calling/AOC20_ATTCCT_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC21_ATCACG_L002/calling/AOC21_ATCACG_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC22_CGATGT_L002/calling/AOC22_CGATGT_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC25_TTAGGC_L002/calling/AOC25_TTAGGC_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC26_TGACCA_L002/calling/AOC26_TGACCA_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC27_ACAGTG_L002/calling/AOC27_ACAGTG_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC28_GCCAAT_L002/calling/AOC28_GCCAAT_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC29_CAGATC_L002/calling/AOC29_CAGATC_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC3_ACTTGA_L002/calling/AOC3_ACTTGA_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC31_ACTTGA_L002/calling/AOC31_ACTTGA_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC32_GATCAG_L003/calling/AOC32_GATCAG_L003.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC33_TAGCTT_L003/calling/AOC33_TAGCTT_L003.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC34_GGCTAC_L003/calling/AOC34_GGCTAC_L003.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC35_ATCACG_L008/calling/AOC35_ATCACG_L008.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC36_CGATGT_L008/calling/AOC36_CGATGT_L008.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC37_TTAGGC_L008/calling/AOC37_TTAGGC_L008.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC38_TGACCA_L008/calling/AOC38_TGACCA_L008.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC39_ACAGTG_L008/calling/AOC39_ACAGTG_L008.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC4_GATCAG_L002/calling/AOC4_GATCAG_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC40_GCCAAT_L008/calling/AOC40_GCCAAT_L008.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC41_CAGATC_L008/calling/AOC41_CAGATC_L008.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC42_ACTTGA_L008/calling/AOC42_ACTTGA_L008.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC5_TAGCTT_L002/calling/AOC5_TAGCTT_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC6_GGCTAC_L002/calling/AOC6_GGCTAC_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC7_CTTGTA_L002/calling/AOC7_CTTGTA_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC8_AGTCAA_L002/calling/AOC8_AGTCAA_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/AOC9_AGTTCC_L002/calling/AOC9_AGTTCC_L002.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/ITSAP104_GTTTCG_L005/calling/ITSAP104_GTTTCG_L005.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/ITSAP37_TAGCTT_L003/calling/ITSAP37_TAGCTT_L003.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/ITSAP53_ATTCCT_L003/calling/ITSAP53_ATTCCT_L003.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/ITSAP73_GATCAG_L004/calling/ITSAP73_GATCAG_L004.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/ITSAP83_ATGTCA_L004/calling/ITSAP83_ATGTCA_L004.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/ITSAP84_CCGTCC_L004/calling/ITSAP84_CCGTCC_L004.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/TES19_ACTGAT_L005/calling/TES19_ACTGAT_L005.g.vcf.gz \
~

```



### Combine gvcf

En 3 parties:



Head : 



```
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



La list



Tail: 

``` -O /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/CORjaune/LesVCF/MetaGenotypes${CHR}.g.vcf.gz
echo "Finnished: "`date`
```



On obtient donc le script Called 

qui sera lancer par le script 



```
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



On doit obtenir un fichier par chromosomes (17 : 16 autosome 1 mitochondrial)



### Check Combine

```
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



En sorti on a les 10 dernières positions pour chaques chromosome dans MetaGenotypesNC_001566.1.g.vcf.check

Position des derniers variants de chaque chromosomes: 

```
 for i in `ls *.check`; do tail -1 ${i}; done

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



On compare a la longueur des chromosome : 



```
cat /home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/HAv3_1_Chromosomes.list
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

Script de genotypage :



```
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



Appelé par le script :



```
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

```
grep complete *_genotype.e
```



Une ligne par chromosome



### Concatenate vcf file

``` 
#!/bin/bash

#concatenateVCFs.bash

module load -f /home/gencel/vignal/save/000_ProgramModules/program_module

bcftools concat -f /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/vcf.list -o MetaGenotypesCalled870.vcf.gz -O z

tabix MetaGenotypesCalled870.vcf.gz
```



### Count variant par chromosomes



```
#!/bin/bash

#statsVcfsRaw.bash

zcat MetaGenotypesCalled870.vcf.gz | grep -v '#' | cut -f 1 | sort | uniq -c | \

awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > countVcfSumRaw
```



```
more countVcfSumRaw

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

Environ 11 millions de variants



### Retiens seulement les SNPs



```
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


### Retiens les INDELs

```
#!/bin/bash

#retainSNPs.bash

module load -f /home/gencel/vignal/save/000_ProgramModules/program_module

OUT=/home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/

gatk --java-options "-Xmx64g" SelectVariants \
     -R /home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
     -V ${OUT}/MetaGenotypesCalled870.vcf.gz \
     --select-type-to-include INDEL \
     -O ${OUT}/MetaGenotypesCalled870_raw_snps.vcf.gz

#Count the SNPs
zcat MetaGenotypesCalled870_raw_snps.vcf.gz | grep -v '#' | cut -f 1 | sort | uniq -c | \
awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > countVcfSumRawINDELs
```


