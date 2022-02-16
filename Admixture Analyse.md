### Admixture Analyse 

Suite au filtrage LD avec les 601 945 SNPs, il est obtenu 590 724 SNPs qui ont donc subit un filtrage LD indirectement, en effet puisqu'ils sont matchés avec des SNPs déjà filtrés par LD=0.3 on considère que ceux la aussi aurait aussi "passer" le filtre LD=0.3 



On va donc utiliser Plink pour obtenir le .bed necessaire à l'Admixture

Ensuite on doit changer les noms de chromosomes en chiffres car Admixture accepte seulement les noms de chromosomes humains. Pour cela on passe les noms de chromosome en chiffres de 1 à 16 

Enfin on effectue l'Admixture



Pour cela, on lance le script Admixture.sh :



```bash
#!/bin/sh 

#Chargement des modules

module load -f /home/agirardon/work/seqapipopOnHAV3_1/program_module/module

#Variables

VCF=/home/agirardon/work/seqapipopOnHAV3_1/combineGVCFs/LesVCF/Concatenate/outisecautresens/MetaGenotypesCalled403_raw_snps_filtre_isec_plink.vcf.gz 
OUT=SeqApiPop403

# Generer les input file dans le format plink

plink --vcf $VCF --allow-extra-chr --make-bed --out $OUT \

# ADMIXTURE accept pas les noms de chromosome qui ne sont pas humains donc on change les noms de chromosomes par chiffre de 1 à 16 

sed -i -e 's/NC_037638.1/1/g' $OUT.bim 
sed -i -e 's/NC_037639.1/2/g' $OUT.bim
sed -i -e 's/NC_037640.1/3/g' $OUT.bim
sed -i -e 's/NC_037641.1/4/g' $OUT.bim
sed -i -e 's/NC_037642.1/5/g' $OUT.bim
sed -i -e 's/NC_037643.1/6/g' $OUT.bim
sed -i -e 's/NC_037644.1/7/g' $OUT.bim
sed -i -e 's/NC_037645.1/8/g' $OUT.bim
sed -i -e 's/NC_037646.1/9/g' $OUT.bim
sed -i -e 's/NC_037647.1/10/g' $OUT.bim
sed -i -e 's/NC_037648.1/11/g' $OUT.bim
sed -i -e 's/NC_037649.1/12/g' $OUT.bim
sed -i -e 's/NC_037650.1/13/g' $OUT.bim
sed -i -e 's/NC_037651.1/14/g' $OUT.bim
sed -i -e 's/NC_037652.1/15/g' $OUT.bim
sed -i -e 's/NC_037653.1/16/g' $OUT.bim




#ADMIXTURE

for K in {2..16}
do
	admixture --cv  $OUT.bed $K > log${K}.out
	
done

```





