#!/bin/csh -f

# Script to run sNMF.
# to use as 'csh cshfile wgsfile popfile sNMFvalueK species'

set wgsfile = $argv[1]

set popfile = $argv[2]

set MAXK = $argv[3]

set species = $argv[4]


if (-e TEMPFILES) then
    echo "TEMPFILES directory already exist"
else
    echo "Creating TEMPFILES directory"
    mkdir TEMPFILES
endif

if (-e RESULTS) then
    echo "RESULTS directory already exist"
else
    echo "Creating RESULTS directory"
    mkdir RESULTS
endif

if (-e RESULTS/$species) then
    echo "RESULTS/$species directory already exist"
else
    echo "Creating RESULTS/$species directory"
    mkdir RESULTS/$species
endif

if (-e RESULTS/$species/sNMF3) then
    echo "RESULTS/$species/sNMF3 directory already exist"
else
    echo "Creating RESULTS/$species/sNMF3 directory"
    mkdir RESULTS/$species/sNMF3
endif

grep '^#CHROM' $wgsfile | awk '{for (i=10; i<=NF; i++) {print $i}}' > $wgsfile.idlist

cp $wgsfile.idlist TEMPFILES/idlist.list

cp $popfile TEMPFILES/totpop.txt

R --no-save < go-orderedpop2031204.R

cp TEMPFILES/orderedpopfile.txt TEMPFILES/pop.txt

echo "Cross-entropy sNMF" > RESULTS/$species/sNMF3/crossentropy.log
foreach inputvcf ($wgsfile)
	grep -v '^#' $inputvcf | wc -l > TEMPFILES/snpnber.txt
	set snpnber = `cat TEMPFILES/snpnber.txt`
	grep '^#CHROM' $inputvcf | awk '{print NF-9}' | uniq > TEMPFILES/idnber.txt
	set idnber = `cat TEMPFILES/idnber.txt`
	vcf2geno $inputvcf $snpnber $idnber
	foreach j (`awk -v MAX=$MAXK 'BEGIN{for(i=5;i<=MAX;i++) print i;}'`)
		sNMF -x genotype.geno -K $j -p 3 -a 16 | tee RESULTS/$species/sNMF3/$inputvcf.$j.out
		cp genotype.$j.Q RESULTS/$species/sNMF3/$inputvcf.$j.Q		
# 		cp genotype.$j.G RESULTS/$species/sNMF3/$inputvcf.$j.G
		foreach i (1 2 3 4 5)
		    createDataSet -x genotype.geno
		    sNMF -x genotype_I.geno -K $j -p 3 -a 16
		    echo "crossentropy.$inputvcf.$j.$i" >> RESULTS/$species/sNMF3/crossentropy.log
		    crossEntropy -x genotype.geno -K $j | grep Cross-Entropy >> RESULTS/$species/sNMF3/crossentropy.log
		    cp genotype_I.$j.Q RESULTS/$species/sNMF3/$inputvcf.$j.$i.Q		
# 		    cp genotype_I.$j.G RESULTS/$species/sNMF3/$inputvcf.$j.$i.G
                    # just rename files for the R script                                                                                   
		    cp RESULTS/$species/sNMF3/$inputvcf.$j.$i.Q TEMPFILES/snmf.txt
                    printf "Ancestry estimates results for %s with %d classes and rep %d by sNMF" $inputvcf:r:r $j $i > TEMPFILES/title.txt
                    R --no-save < drawsNMF.R | tee TEMPFILES/logRdrawsnmf.$inputvcf.$j.$i.log
                    cp TEMPFILES/output.pdf RESULTS/$species/sNMF3/$inputvcf.$j.$i.pdf			
		end
	end
end
bgzip $wgsfile
rm *.G