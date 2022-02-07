###  1. Introduction

Ce document présente les analyses faites pour obtenir les résultats qui sont présentés par la suite. Pour la répétabilité ou les différentes listes les chemins sont à changer et à adapter.

Les version des software utilisés sont les suivantes:

Ici seulement les modules utilisé dans les scripts suivants et non tous les modules dans le fichiers 

more module 

``` 
 
#%Module1.0###############################################################


module load bioinfo/seqtk-1.2
module load bioinfo/bwa-0.7.15
module load bioinfo/samtools-1.8
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/bcftools-1.6
module load bioinfo/picard-2.18.2
module load bioinfo/gatk-4.1.2.0

```



### Genome de référence

Le génome de référence utilisé est Amel_HAv3.1 de Genebank : GCA_003254395.2

### fastq file 

Les fichiers fastq individuel des population sont obtenus par sequencage Illumina. 

### Préparation du mapping

Pour chaque run on a R1_fastq.gz et R2_fastq.gz avec le forward et le reverse read 

### Recuperation des chemins de fastq 

pour recupurer la liste : 

```
ls /genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune* | awk 'BEGIN{FS=OFS="_"}{print $1,$2,$3}' | sort | uniq > corse_yellow.list
```

head -3 corse_yellow.list

/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune11_TCTCTACT-GAACCGCG-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune12_CTCTCGTC-AGGTTATA-AHJJN2DSX2_L003

corse_yellow.list diviser en 7 pour pouvoir lancer les jobs petits à petit 



### Mapping 



Les script de Mapping effectue le mapping avec BWA, marquage des lectures avec Picard, GATK BQSR( Base Quality Score Replication) 


Map_seqapipop_HAV3_1.bash 

Ce script de mapping appelle 3 scripts en tout : 

- mappingAV_2019_Dec.sh
- bootstrapingAV_2019_Dec.sh
- callingAV_2019_Dec.sh



Editer les chemins:

SAMPLE_FILE
OUT


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

A REVOIR : 

La PLOIDY  est fixé à 2, bien que normalement les individus sont haploïdes, des tests ont été fait avec PLOIDY=1 cela conduisait a un faux genotype dans les sequences répétées, et donnaient des SNP faussement positifs 

N=1  car il peut être utilisé pour un pool d'individus

Explications plus précises dans le pdf : Mapping_Explication_PLOIDY=2.pdf

Ce script prend une liste de chemin vers des noms d'echantillons : 



head -3 corse_yellow.list

```
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune11_TCTCTACT-GAACCGCG-AHJJN2DSX2_L003
/genphyse/cytogen/seqapipop/FastqFromNG6/CORjaune12_CTCTCGTC-AGGTTATA-AHJJN2DSX2_L003
```


On obtiendra donc plusieurs fichier :



Les .bam en appelant mappingAV_2019_Dec.sh il effectue le mapping et la detection des lectures dupliquées :

- BWA
- Picard MarkDuplicates 

Le BQSR, en appelant bootstrapingAV_2019_Dec.sh fait un recalibrage du score de qualité de base sur chaque chromosome 

- CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003.bam
- CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003_bam.list
- CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003.bam.bai



et appelle le calling, callingAV_2019_Dec.sh qui utilise les .bam pour détecter les SNPs grâce au HaplotypeCaller

- CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003.g.vcf.gz
- CORjaune10_GACCTGAA-CTCACCAA-AHJJN2DSX2_L003.g.vcf.gz.tbi



### Control gvcf are complete 

Un contrôle est effectué pour vérifié que les jobs se sont terminés correctement ( interruption possible, etc...)

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
       bcftools query -f '%CHROM\n' ${i} | \ # bcftools query extrait un champ des VCF -- -f format de sortie -- %CHROM on veut la colone CHROM soit print -- \n nouvelle ligne 
       grep ^NC | uniq -c > \ # grep : recherche  un expression regulière commencant par NC -- uniq -c : Compter les occurrences des lignes en double
       /home/agirardon/work/seqapipopOnHAV3_1/controlling/controlVcfs/${ID}.count"
done
```


Fichier de sortie : 

```
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

avec la commande 

En effet on verifie qu'on est bien arrivé jusqu'au chromosome mitochondrial qui est le dernier à être lancé 


``` 
grep NC_001566.1 *.count | wc -l
```



### Recuperation list 



récupération des chemins pour la suite: 

```
#!/bin/bash

ls /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/*/calling/*vcf.gz > chemin.all

for i in `awk '{print$1}' /work/genphyse/cytogen/Alain/seqapipopOnHAV3_1/seqApiPopVcfFilteredSonia/plinkAnalyses/WindowSNPs/RFMix/Pure95/in/IndsPopReference.list`
do
	grep ${i}_ chemin.all | awk '{print "    --variant",$1,"\\"}'  >> ref.list 
done 

rm chemin.all
```



ca donnera : 

head -3 ref.list

``` 

    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/Ab-PacBio_TCCGCGAA-CCTATCCT-BHCFFJDSXX_L004/calling/Ab-PacBio_TCCGCGAA-CCTATCCT-BHCFFJDSXX_L004.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/BER10_GTCCGC_L003/calling/BER10_GTCCGC_L003.g.vcf.gz \
    --variant /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/Haploid/BER11_GTGAAA_L003/calling/BER11_GTGAAA_L003.g.vcf.gz \
  

```



### Combine gvcf

En 3 parties:



- Head : 

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



- La liste



- Tail: 

``` 
-O /genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3_1/CORjaune/LesVCF/MetaGenotypes${CHR}.g.vcf.gz
echo "Finnished: "`date`
```



On obtient donc le script Called 

Ce script permet de lancer les différents chromosomes en parallèle, grâce au passage de la variable ${i} à l'aide de l'option -c.

qui sera lancé par le script 



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

Vérifie que l'on va bien au bout des jobs(risque de jobs interrompus possibles pour diverses raisons..)



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



En sortie on a les 10 dernières positions pour chaques chromosome dans MetaGenotypesNC_001566.1.g.vcf.check

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



On compare a la longueur des chromosome et verifier qu'ils sont proches: 



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

Globalement le même fonctionnement que le script précédent

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


Cette ligne apparait seulement si tout s'est bien passé, donc on doit obtenir une ligne par chromosome



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

```
more countVcfSumRawSNPs
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

Environ 8 millions de SNPs 

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
     -O ${OUT}/MetaGenotypesCalled870_raw_indels.vcf.gz

#Count the SNPs
zcat MetaGenotypesCalled870_raw_snps.vcf.gz | grep -v '#' | cut -f 1 | sort | uniq -c | \
awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > countVcfSumRawINDELs
```
``` 
more countVcfSumRawINDELs
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

Environ 2 900 000 Indels 


### Récuperation liste SNP 

- Il faut préparer 2 VCF avec 2 listes de SNP différentes

1 apres tous les filtres, avec 7 millions de SNP 

 Récupération de la liste des 7 millions de SNP 

```#!/bin/sh



module load -f /home/gencel/vignal/save/000_ProgramModules/program_module

bcftools query -f '%CHROM %POS\n' /work/genphyse/cytogen/Alain/seqapipopOnHAV3_1/seqApiPopVcfFilteredSonia/vcf_cleanup/MetaGenotypesCalled870_raw_snps_allfilter.vcf.gz -o /home/agirardon/work//seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/list7m.list
```

le fichier ressemble à : 


``` 
head -3 recup_list.bash

NC_037638.1 5671
NC_037638.1 5698
NC_037638.1 6621
``` 

Pour ce qui est de la liste des 600 000 sous forme de plink 

```
#!/bin/sh



awk '{print $2, $4}' /work/genphyse/cytogen/Alain/seqapipopOnHAV3_1/seqApiPopVcfFilteredSonia/plinkAnalyses/WindowSNPs/SeqApiPop_629_maf001_LD03_pruned.bim > /home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/snplist600k.txt

```

le fichier ressemble à : 

```
head -3 snplist600k.txt
1:7577[HAV3.1]AG
1:10360[HAV3.1]CT
1:11791[HAV3.1]AG

```


