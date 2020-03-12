#!/bin/bash
# Author: Yassine Souilmi <yassine.souilmi@adelaide.edu.au>

set -u

# module load BCFtools

bed=$1 # bed file name
ref = $2 # link to the reference fasta file

# Create a polarised VCF for each population:
awk '{print $1}' ${bed}.fam | sort -u | while read pop
do
	echo "${pop}" > pop.txt
	mkdir -pv ${bed}_popSplit

	plink --bfile ${bed} \
	--keep-fam pop.txt \
	--recode vcf-iid bgz \
	--keep-allele-order \
	--out ${bed}_popSplit/${pop}_${bed}

	tabix -p vcf ${bed}_popSplit/${pop}_${bed}.vcf.gz

	bcftools norm \
	-t $(seq 1 22 | tr -s '\n' ',' | sed 's/22,/22/') \
	-O z -o ${bed}_popSplit/${pop}_${bed}_fixedRefAllele.vcf.gz \
	-c s \
	-f ~/fastdir/Refs/gatk-bundle/b37/human_g1k_v37_decoy.fasta \
	${bed}_popSplit/${pop}_${bed}.vcf.gz

	tabix -p vcf ${bed}_popSplit/${pop}_${bed}_fixedRefAllele.vcf.gz

done

# Create SFS files
mkdir -pv ${bed}_popSplit_SF2_inputs

for chr in {1..22}
do
	python ~/bin/vcf2SF.py \
		-g \
		-v ${bed}_popSplit/${pop}_${bed}_fixedRefAllele.vcf.gz \
		-c ${chr} \
		-o ${bed}_popSplit_SF2_inputs/${pop}_${bed}_${chr}.sfs
done


# Running SweepFinder2
wd=$PWD

awk '{print $1}' ${bed}.fam | sort -u | while read pop
	for chr in {1..22}
	do
    	mkdir -pv ${pop}_sf2_wd/${chr}
    	cd ${pop}_sf2_wd/${chr}/ #Making sure all jobs run in a dedicated folder

    	SweepFinder2 -sg 1000 ${wd}/${bed}_popSplit_SF2_inputs/${pop}_${bed}_${chr}.sfs ./${pop}_${chr}.SF2out

	done # for loop (chromosomes)
done # while loop (populations)
