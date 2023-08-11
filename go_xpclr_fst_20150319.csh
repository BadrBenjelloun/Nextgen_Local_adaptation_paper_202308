#!/bin/csh -f


# script to prepare and run xpclr with the assumption of 1cM = 1M bp It requires chromosome number ; the name of a gz vcf file (without the extention ".vcf.gz") the name of pop1 with individuals of this pop are listed in a file named $pop1.txt . It's the same for pop2

if (-e TEMPFILES) then
    echo "TEMPFILES directory already exist"
else
    echo "Creating TEMPFILES directory"
    mkdir TEMPFILES
endif

if (-e RESULTS3) then
    echo "RESULTS3 directory already exist"
else
    echo "Creating RESULTS3 directory"
    mkdir RESULTS3
endif

if (-e POP) then
    echo "POP directory already exist"
else
    echo "Creating POP directory"
    mkdir POP
endif


set chrnb = $argv[1]
set vcffile = $argv[2]
set pop1 = $argv[3]
set pop2 = $argv[4]

cp POP/$pop1.txt POP/totpop.txt
cat POP/$pop2.txt >> POP/totpop.txt 
 
foreach chr (`awk -v a=$chrnb 'BEGIN{for(k=1;k<=a;k++) print k;}'`)
    vcftools --gzvcf ../bio7_analysis/$vcffile.vcf.gz --keep POP/totpop.txt --chr $chr --geno 1 --out POP/$vcffile.$pop1.$pop2.chr.$chr --recode
end

foreach chr (`awk -v a=$chrnb 'BEGIN{for(k=1;k<=a;k++) print k;}'`)
	foreach pop ($pop1 $pop2)
		vcftools --vcf POP/$vcffile.$pop1.$pop2.chr.$chr.recode.vcf --remove-indels --max-alleles 2 --min-alleles 2 --keep POP/$pop.txt --chr $chr --from-bp 0 --to-bp 26000000 --out POP/$vcffile.$pop.chr.$chr.0_25M --recode
		vcftools --vcf POP/$vcffile.$pop1.$pop2.chr.$chr.recode.vcf --remove-indels --max-alleles 2 --min-alleles 2 --keep POP/$pop.txt --chr $chr --from-bp 24000001 --to-bp 51000000 --out POP/$vcffile.$pop.chr.$chr.25M_50M --recode
		vcftools --vcf POP/$vcffile.$pop1.$pop2.chr.$chr.recode.vcf --remove-indels --max-alleles 2 --min-alleles 2 --keep POP/$pop.txt --chr $chr --from-bp 49000001 --to-bp 76000000 --out POP/$vcffile.$pop.chr.$chr.50M_75M --recode		
		vcftools --vcf POP/$vcffile.$pop1.$pop2.chr.$chr.recode.vcf --remove-indels --max-alleles 2 --min-alleles 2 --keep POP/$pop.txt --chr $chr --from-bp 74000001 --to-bp 101000000 --out POP/$vcffile.$pop.chr.$chr.75M_100M --recode
		vcftools --vcf POP/$vcffile.$pop1.$pop2.chr.$chr.recode.vcf --remove-indels --max-alleles 2 --min-alleles 2 --keep POP/$pop.txt --chr $chr --from-bp 99000001 --to-bp 126000000 --out POP/$vcffile.$pop.chr.$chr.100M_125M --recode
		vcftools --vcf POP/$vcffile.$pop1.$pop2.chr.$chr.recode.vcf --remove-indels --max-alleles 2 --min-alleles 2 --keep POP/$pop.txt --chr $chr --from-bp 124000001 --to-bp 151000000 --out POP/$vcffile.$pop.chr.$chr.125M_150M --recode
		vcftools --vcf POP/$vcffile.$pop1.$pop2.chr.$chr.recode.vcf --remove-indels --max-alleles 2 --min-alleles 2 --keep POP/$pop.txt --chr $chr --from-bp 149000001 --to-bp 176000000 --out POP/$vcffile.$pop.chr.$chr.150M_175M --recode
		vcftools --vcf POP/$vcffile.$pop1.$pop2.chr.$chr.recode.vcf --remove-indels --max-alleles 2 --min-alleles 2 --keep POP/$pop.txt --chr $chr --from-bp 174000001 --to-bp 201000000 --out POP/$vcffile.$pop.chr.$chr.175M_200M --recode
		vcftools --vcf POP/$vcffile.$pop1.$pop2.chr.$chr.recode.vcf --remove-indels --max-alleles 2 --min-alleles 2 --keep POP/$pop.txt --chr $chr --from-bp 199000001 --to-bp 226000000 --out POP/$vcffile.$pop.chr.$chr.200M_225M --recode
		vcftools --vcf POP/$vcffile.$pop1.$pop2.chr.$chr.recode.vcf --remove-indels --max-alleles 2 --min-alleles 2 --keep POP/$pop.txt --chr $chr --from-bp 224000001 --to-bp 251000000 --out POP/$vcffile.$pop.chr.$chr.225M_250M --recode
		vcftools --vcf POP/$vcffile.$pop1.$pop2.chr.$chr.recode.vcf --remove-indels --max-alleles 2 --min-alleles 2 --keep POP/$pop.txt --chr $chr --from-bp 249000001 --to-bp 276000000 --out POP/$vcffile.$pop.chr.$chr.250M_275M --recode
		vcftools --vcf POP/$vcffile.$pop1.$pop2.chr.$chr.recode.vcf --remove-indels --max-alleles 2 --min-alleles 2 --keep POP/$pop.txt --chr $chr --from-bp 274000001 --to-bp 286000000 --out POP/$vcffile.$pop.chr.$chr.275M_285M --recode
		foreach pos (0_25M 25M_50M 50M_75M 75M_100M 100M_125M 125M_150M 150M_175M 175M_200M 200M_225M 225M_250M 250M_275M 275M_285M)
			vcftools --vcf POP/$vcffile.$pop.chr.$chr.$pos.recode.vcf --remove-indels --IMPUTE --out POP/$vcffile.$pop.chr.$chr.$pos
		end
#		awk '{for (i=1; i<=NF; i++) a[i]=a[i](NR!=1?FS:"")$i} END {for (i=1; i in a; i++) print a[i]}' POP/$vcffile.$pop.chr.$chr.impute.hap > POP/$vcffile.$pop.chr.$chr.ihshap
	end
end

foreach chr (`awk -v a=$chrnb 'BEGIN{for(k=1;k<=a;k++) print k;}'`)
	foreach pos (0_25M 25M_50M 50M_75M 75M_100M 100M_125M 125M_150M 150M_175M 175M_200M 200M_225M 225M_250M 250M_275M 275M_285M)
		echo "chr $chr pos $pos"
		grep -v '^#' POP/$vcffile.$pop1.chr.$chr.$pos.recode.vcf | gawk '{print "rs"$1"0"$2,$2}' > $vcffile.part1.map
#		awk 'FNR > 1 {if (FNR==2) {print "0"} else {print $2}}' RESULTS3/ldhat/$vcffile.$k.chr.$chr.$inipos.$lastpos.res.txt > $vcffile.$k.part2.map
		grep -v '^#' POP/$vcffile.$pop1.chr.$chr.$pos.recode.vcf | awk '{print $4,$5}' > $vcffile.part3.map
		paste $vcffile.part1.map $vcffile.part3.map | awk '{OFS=" "} {if (FNR==1) {print $1,$2,$2/100000000,$3,$4} else {print $1,$2,"0",$3,$4}}' > $vcffile.chr.$chr.$pos.ihsmap
		rm $vcffile.part1.map
#		rm $vcffile.$k.part2.map
		rm $vcffile.part3.map
# 		awk -v a=$chr '{print $1,a,($2/100000000)*1.06556,$2,$4,$5}' $vcffile.$k.chr.$chr.$inipos.$lastpos.ihsmap > $vcffile.$k.chr.$chr.$inipos.$lastpos.xpclrmap
		gawk -v a=$chr '{print $1,a,$2/100000000,$2,$4,$5}' $vcffile.chr.$chr.$pos.ihsmap > POP/$vcffile.chr.$chr.$pos.xpclrmap
		echo "run xpclr for $chr in pos $pos"
		mv POP/$vcffile.$pop1.chr.$chr.$pos.impute.hap POP/$vcffile.$pop1.chr.$chr.$pos.xpclr
		mv POP/$vcffile.$pop2.chr.$chr.$pos.impute.hap POP/$vcffile.$pop2.chr.$chr.$pos.xpclr
		XPCLR -xpclr POP/$vcffile.$pop1.chr.$chr.$pos.xpclr POP/$vcffile.$pop2.chr.$chr.$pos.xpclr POP/$vcffile.chr.$chr.$pos.xpclrmap xpclr.$vcffile.chr.$chr.$pos -w1 0.0005 250 2500 $chr -p1 0.95
	end
	awk '{if($4<=25000000 && $4>0){print}}' xpclr.$vcffile.chr.$chr.0_25M.wtclr.txt > RESULTS3/$vcffile.chr.$chr.$pop1.$pop2.xpclr.out
	rm xpclr.$vcffile.chr.$chr.0_25M.wtclr.txt
	awk '{if($4<=50000000 && $4>25000000){print}}' xpclr.$vcffile.chr.$chr.25M_50M.wtclr.txt >> RESULTS3/$vcffile.chr.$chr.$pop1.$pop2.xpclr.out
	rm xpclr.$vcffile.chr.$chr.25M_50M.wtclr.txt
	awk '{if($4<=75000000 && $4>50000000){print}}' xpclr.$vcffile.chr.$chr.50M_75M.wtclr.txt >> RESULTS3/$vcffile.chr.$chr.$pop1.$pop2.xpclr.out
	rm xpclr.$vcffile.chr.$chr.50M_75M.wtclr.txt
	awk '{if($4<=100000000 && $4>75000000){print}}'  xpclr.$vcffile.chr.$chr.75M_100M.wtclr.txt >> RESULTS3/$vcffile.chr.$chr.$pop1.$pop2.xpclr.out
	rm xpclr.$vcffile.chr.$chr.75M_100M.wtclr.txt
	awk '{if($4<=125000000 && $4>100000000){print}}'  xpclr.$vcffile.chr.$chr.100M_125M.wtclr.txt >> RESULTS3/$vcffile.chr.$chr.$pop1.$pop2.xpclr.out
	rm xpclr.$vcffile.chr.$chr.100M_125M.wtclr.txt
	awk '{if($4<=150000000 && $4>125000000){print}}'  xpclr.$vcffile.chr.$chr.125M_150M.wtclr.txt >> RESULTS3/$vcffile.chr.$chr.$pop1.$pop2.xpclr.out
	rm xpclr.$vcffile.chr.$chr.125M_150M.wtclr.txt
	awk '{if($4<=175000000 && $4>150000000){print}}'  xpclr.$vcffile.chr.$chr.150M_175M.wtclr.txt >> RESULTS3/$vcffile.chr.$chr.$pop1.$pop2.xpclr.out
	rm xpclr.$vcffile.chr.$chr.150M_175M.wtclr.txt
	awk '{if($4<=200000000 && $4>175000000){print}}'  xpclr.$vcffile.chr.$chr.175M_200M.wtclr.txt >> RESULTS3/$vcffile.chr.$chr.$pop1.$pop2.xpclr.out
	rm xpclr.$vcffile.chr.$chr.175M_200M.wtclr.txt
	awk '{if($4<=225000000 && $4>200000000){print}}'  xpclr.$vcffile.chr.$chr.200M_225M.wtclr.txt >> RESULTS3/$vcffile.chr.$chr.$pop1.$pop2.xpclr.out
	rm xpclr.$vcffile.chr.$chr.200M_225M.wtclr.txt
	awk '{if($4<=250000000 && $4>225000000){print}}'  xpclr.$vcffile.chr.$chr.225M_250M.wtclr.txt >> RESULTS3/$vcffile.chr.$chr.$pop1.$pop2.xpclr.out
	rm xpclr.$vcffile.chr.$chr.225M_250M.wtclr.txt
	awk '{if($4<=275000000 && $4>250000000){print}}'  xpclr.$vcffile.chr.$chr.250M_275M.wtclr.txt >> RESULTS3/$vcffile.chr.$chr.$pop1.$pop2.xpclr.out
	rm xpclr.$vcffile.chr.$chr.250M_275M.wtclr.txt
	awk '{if($4<=285000000 && $4>275000000){print}}'  xpclr.$vcffile.chr.$chr.275M_285M.wtclr.txt >> RESULTS3/$vcffile.chr.$chr.$pop1.$pop2.xpclr.out
	rm xpclr.$vcffile.chr.$chr.275M_285M.wtclr.txt
	rm *.ihsmap
	cp RESULTS3/$vcffile.chr.$chr.$pop1.$pop2.xpclr.out TEMPFILES/xpclr.txt
	printf "xpclr score between %s and %s pop for %s in chr %d" $pop1 $pop2 $vcffile $chr > TEMPFILES/xpclrtitle.txt
	R --no-save < drawxpclr.R | tee TEMPFILES/logRdrawxpclr.$vcffile.chr.$chr.log    	
	cp TEMPFILES/outputxpclr.pdf RESULTS3/xpclr.$vcffile.$pop1.$pop2.chr.$chr.pdf
end

rm POP/$vcffile.*.chr.*.log
rm POP/$vcffile.*.chr.*.recode.vcf*
rm POP/$vcffile.*.xpclr
rm POP/$vcffile.*.xpclrmap
rm POP/$vcffile.*.impute.legend
rm POP/$vcffile.*.impute.hap.indv
