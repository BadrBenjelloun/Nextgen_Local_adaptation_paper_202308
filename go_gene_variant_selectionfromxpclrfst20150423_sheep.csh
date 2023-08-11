#!/bin/csh -f                                                                                                                                             


# script to group xpclr scores and to select up to a treshold and to check the genes under selection and find the genes detected in 1 comparison using Fst information also; It requires chromosome number ; the name of a gz vcf file (without the extention ".vcf.gz") the name of pop1 ,pop2, the variant effect predictor file and requires the percenage to keep, the window size used for xpclr analyses and the Fst treshold                                                                                                                                             
# Results are in a directory named RESULTS2/

if (-e TEMPFILES) then
    echo "TEMPFILES directory already exist"
else
    echo "Creating TEMPFILES directory"
    mkdir TEMPFILES
endif

if (-e RESULTS2) then
    echo "RESULTS2 directory already exist"
else
    echo "Creating RESULTS2 directory"
    mkdir RESULTS2
endif


set chrnb = $argv[1]
set vcffile = $argv[2]
set pop1 = $argv[3]
set pop2 = $argv[4]
set percentagekept = $argv[5]
set winxpclr = $argv[6]
set percentagekeptfst = $argv[7]
set vep = $argv[8]
set folder = $argv[9]
set corsp = $argv[10] # for sheep only
set para = $argv[11]

echo "first comparison $pop1 $pop2"
echo "" > RESULTS2/$vcffile.allchr.$pop1.$pop2.xpclr.out
foreach chr (`awk -v a=$chrnb 'BEGIN{for(k=1;k<=a;k++) print k;}'`)
	cat RESULTS2/$vcffile.chr.$chr.$pop1.$pop2.xpclr.out >> RESULTS2/$vcffile.allchr.$pop1.$pop2.xpclr.out
end
awk '(FNR>1){print}' RESULTS2/$vcffile.allchr.$pop1.$pop2.xpclr.out > tmp.txt
mv tmp.txt RESULTS2/$vcffile.allchr.$pop1.$pop2.xpclr.out
set nbvariant = `cat RESULTS2/$vcffile.allchr.$pop1.$pop2.xpclr.out | wc -l`

set tresh = `echo "$nbvariant*$percentagekept" | bc`

awk '{print $6}' RESULTS2/$vcffile.allchr.$pop1.$pop2.xpclr.out | sort -n | awk -v a=$tresh '(FNR>a) {print}' | awk '(FNR==1){print}' > TEMPFILES/treshscore.txt
set treshscore = `cat TEMPFILES/treshscore.txt`
echo "treshold of xpclr score is $treshscore"
awk -v a=$treshscore '($6>=a){print}' RESULTS2/$vcffile.allchr.$pop1.$pop2.xpclr.out > RESULTS2/$vcffile.allchr.$pop1.$pop2.highxpclrscore.$percentagekept.txt

awk '{print $6"\t"$1"\t"$4}' RESULTS2/$vcffile.allchr.$pop1.$pop2.highxpclrscore.$percentagekept.txt | sort -r | awk -v a=$winxpclr '{print $2"\t"$3-(a/2)"\t"$3+(a/2)}' > RESULTS2/$vcffile.allchr.$pop1.$pop2.highxpclrscore.$percentagekept.bed

echo "chr start end score" > RESULTS2/allchr.$pop1.$pop2.highxpclrscore.$percentagekept.score.txt
awk '{print $1"\t"$4"\t"$6}' RESULTS2/$vcffile.allchr.$pop1.$pop2.highxpclrscore.$percentagekept.txt | awk -v a=$winxpclr '{print $1"\t"$2-(a/2)"\t"$2+(a/2)"\t"$3}' >> RESULTS2/allchr.$pop1.$pop2.highxpclrscore.$percentagekept.score.txt
awk -f ../overlapingwin.awk RESULTS2/allchr.$pop1.$pop2.highxpclrscore.$percentagekept.score.txt RESULTS2/allchr.$pop1.$pop2.highxpclrscore.$percentagekept.score.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto.txt
awk -f ../overlapingwin.awk toto.txt toto.txt | awk 'FNR>1{printf ("%s %012d %012d %012s\n",$1,$5,$6,$4);}' | sort > TEMPFILES/scores.txt
awk -f ../overlapingwin.awk toto.txt toto.txt | awk 'FNR>1{printf ("%s %012d %012d\n",$1,$5,$6);}' | sort | uniq -c | awk '{printf ("%s %012d %012d "1" %s\n", $2,$3,$4,$1);}' > TEMPFILES/gridpoints.txt
cat TEMPFILES/scores.txt TEMPFILES/gridpoints.txt | sort | awk 'NF==4{chr=$1;start=$2;end=$3;score=$4;next;}NF==5{if($1==chr&&$2==start&&$3==end){print chr,start,end,score,$5}next;}' | awk '{printf("%012s %s %s %s %s\n", $4,$1,$2,$3,$5);}' | sort -r | awk '{print $2,$3,$4,$1,$5,NR}' > RESULTS2/allchr.$pop1.$pop2.highxpclrscore.$percentagekept.grpwindows.txt


echo "second comparison $pop2 $pop1"
echo "" > RESULTS2/$vcffile.allchr.$pop2.$pop1.xpclr.out
foreach chr (`awk -v a=$chrnb 'BEGIN{for(k=1;k<=a;k++) print k;}'`)
	cat RESULTS2/$vcffile.chr.$chr.$pop2.$pop1.xpclr.out >> RESULTS2/$vcffile.allchr.$pop2.$pop1.xpclr.out
end
awk '(FNR>1){print}' RESULTS2/$vcffile.allchr.$pop2.$pop1.xpclr.out > tmp.txt
mv tmp.txt RESULTS2/$vcffile.allchr.$pop2.$pop1.xpclr.out
set nbvariant = `cat RESULTS2/$vcffile.allchr.$pop2.$pop1.xpclr.out | wc -l`

set tresh = `echo "$nbvariant*$percentagekept" | bc`

awk '{print $6}' RESULTS2/$vcffile.allchr.$pop2.$pop1.xpclr.out | sort -n | awk -v a=$tresh '(FNR>a) {print}' | awk '(FNR==1){print}' > TEMPFILES/treshscore.txt
set treshscore = `cat TEMPFILES/treshscore.txt`
echo "treshold of xpclr score is $treshscore"
awk -v a=$treshscore '($6>=a){print}' RESULTS2/$vcffile.allchr.$pop2.$pop1.xpclr.out > RESULTS2/$vcffile.allchr.$pop2.$pop1.highxpclrscore.$percentagekept.txt

awk '{print $6"\t"$1"\t"$4}' RESULTS2/$vcffile.allchr.$pop2.$pop1.highxpclrscore.$percentagekept.txt | sort -r | awk -v a=$winxpclr '{print $2"\t"$3-(a/2)"\t"$3+(a/2)}' > RESULTS2/$vcffile.allchr.$pop2.$pop1.highxpclrscore.$percentagekept.bed

echo "chr start end score" > RESULTS2/allchr.$pop2.$pop1.highxpclrscore.$percentagekept.score.txt
awk '{print $1"\t"$4"\t"$6}' RESULTS2/$vcffile.allchr.$pop2.$pop1.highxpclrscore.$percentagekept.txt | awk -v a=$winxpclr '{print $1"\t"$2-(a/2)"\t"$2+(a/2)"\t"$3}' >> RESULTS2/allchr.$pop2.$pop1.highxpclrscore.$percentagekept.score.txt
awk -f ../overlapingwin.awk RESULTS2/allchr.$pop2.$pop1.highxpclrscore.$percentagekept.score.txt RESULTS2/allchr.$pop2.$pop1.highxpclrscore.$percentagekept.score.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto.txt
awk -f ../overlapingwin.awk toto.txt toto.txt | awk 'FNR>1{printf ("%s %012d %012d %012s\n",$1,$5,$6,$4);}' | sort > TEMPFILES/scores.txt
awk -f ../overlapingwin.awk toto.txt toto.txt | awk 'FNR>1{printf ("%s %012d %012d\n",$1,$5,$6);}' | sort | uniq -c | awk '{printf ("%s %012d %012d "1" %s\n", $2,$3,$4,$1);}' > TEMPFILES/gridpoints.txt
cat TEMPFILES/scores.txt TEMPFILES/gridpoints.txt | sort | awk 'NF==4{chr=$1;start=$2;end=$3;score=$4;next;}NF==5{if($1==chr&&$2==start&&$3==end){print chr,start,end,score,$5}next;}' | awk '{printf("%012s %s %s %s %s\n", $4,$1,$2,$3,$5);}' | sort -r | awk '{print $2,$3,$4,$1,$5,NR}' > RESULTS2/allchr.$pop2.$pop1.highxpclrscore.$percentagekept.grpwindows.txt

if (! (-e RESULTS2/$vcffile.$pop1.$pop2.nonan.weir.fst) ) then
echo "Fst estimates between $pop1 $pop2"
vcftools --gzvcf ../$folder/$vcffile.vcf.gz --weir-fst-pop POP/$pop1.txt --weir-fst-pop POP/$pop2.txt --out RESULTS2/$vcffile.$pop1.$pop2
grep -v "nan" RESULTS2/$vcffile.$pop1.$pop2.weir.fst | awk '(FNR>1){print}' > RESULTS2/$vcffile.$pop1.$pop2.nonan.weir.fst
rm RESULTS2/$vcffile.$pop1.$pop2.weir.fst   #new
endif

set nbvariants = `cat RESULTS2/$vcffile.$pop1.$pop2.nonan.weir.fst | wc -l`

set treshfst = `echo "$nbvariants*$percentagekeptfst" | bc`
awk '{printf("%012f\n", $3);}' RESULTS2/$vcffile.$pop1.$pop2.nonan.weir.fst | sort -n | awk -v a=$treshfst '(FNR>a) {print}' | awk '(FNR==1){print}' > TEMPFILES/treshfstscore.txt
set treshfstscore = `cat TEMPFILES/treshfstscore.txt`
echo "treshold of fst is $treshfstscore"
awk -v a=$treshfstscore '($3>=a){print}' RESULTS2/$vcffile.$pop1.$pop2.nonan.weir.fst > RESULTS2/$vcffile.$pop1.$pop2.fst.$percentagekeptfst.txt

echo "searching variants detected by xpclr using Fst $pop1 $pop2"
if (! (-e RESULTS2/$vcffile.allchr.$pop1.$pop2.highxpclrscore.$percentagekept.weir.fst) ) then
vcftools --gzvcf ../$folder/$vcffile.vcf.gz --bed RESULTS2/$vcffile.allchr.$pop1.$pop2.highxpclrscore.$percentagekept.bed --out TEMPFILES/positions --recode
vcftools --vcf TEMPFILES/positions.recode.vcf --weir-fst-pop POP/$pop1.txt --weir-fst-pop POP/$pop2.txt --out RESULTS2/$vcffile.allchr.$pop1.$pop2.highxpclrscore.$percentagekept
endif

echo "chr	position	fst" > RESULTS2/$pop1.$pop2.selectionposition.$percentagekeptfst.txt
grep -v "nan" RESULTS2/$vcffile.allchr.$pop1.$pop2.highxpclrscore.$percentagekept.weir.fst | awk -v a=$treshfstscore 'FNR>1 && $3>=a {print $1"\t"$2"\t"$3}' >> RESULTS2/$pop1.$pop2.selectionposition.fst.$percentagekeptfst.txt
awk 'FNR>1{printf ("%s %012d "1" %s\n", $1,$2,$3);}' RESULTS2/$pop1.$pop2.selectionposition.fst.$percentagekeptfst.txt > TEMPFILES/variantsS.txt    # new
awk '{printf ("%s %s "0" %s %s %s %s\n", $1,$2,$3,$4,$5,$6);}' RESULTS2/allchr.$pop1.$pop2.highxpclrscore.$percentagekept.grpwindows.txt > TEMPFILES/scoresS.txt    #new
cat TEMPFILES/scoresS.txt TEMPFILES/variantsS.txt | sort | awk 'NF==7{chr=$1;start=$2;end=$4;xpclr=$5;points=$6;rank=$7;next;}NF==4{if($1==chr && $2>=start && $2<=end){print chr"\t"$2"\t"xpclr"\t"points"\t"rank"\t"$4}next;}' | sort | uniq > RESULTS2/$pop1.$pop2.selectionposition.fst.$percentagekeptfst.xpclr.$percentagekept.txt   #new


awk -F" " '{printf("%s %012d 1 %s\n", $1,$2,$3);}' RESULTS2/$pop1.$pop2.selectionposition.fst.$percentagekeptfst.txt > TEMPFILES/positionSSselection.txt


if (! (-e ../vep.txt) ) then
gzip -d -c $vep | grep -v "#" | awk '{print $2"\t"$4"\t"$7}' | awk -F":" '{print $1"\t"$2}' | awk '{printf("%s %012d 0 %s %s\n",$1,$2,$3,$4);}' > TEMPFILES/vep.txt
awk '{printf("%s "1" %s %s %s %s\n", $4,$1,$2,$3,$5);}' TEMPFILES/vep.txt > TEMPFILES/vep2.txt #for sheep only
awk '{printf ("%s "0" %s\n",$1,$2);}' $corsp > TEMPFILES/corsp.txt #for sheep only
cat TEMPFILES/vep2.txt TEMPFILES/corsp.txt | sort | awk 'NF==3{id=$1;name=$3;next}NF==6{if($1==id && id!="-" && name!="ensembl"){print $3,$4,$5,name,$6}else{print $3,$4,$5,$1,$6}next;}' > ../vep.txt #for sheep only
rm TEMPFILES/vep*txt
endif

awk -F" " '{printf("%s %s 1 %s %s %s %s\n", $1,$2,$3,$4,$5,$6);}' RESULTS2/$pop1.$pop2.selectionposition.fst.$percentagekeptfst.xpclr.$percentagekept.txt > TEMPFILES/positionSSselection2.txt #new
cat ../vep.txt TEMPFILES/positionSSselection2.txt | sort | awk 'NF==5{gene=$4;chr=$1;variant=$2;type=$5;next;}NF==7{if ($1==chr && $2==variant) {print $1,$2,$4,$5,gene,$6,type,$7}next;}' > RESULTS2/$pop1.$pop2.xpclr.$percentagekept.fst.$percentagekeptfst.selectionpositiongenesscores.txt # new
awk '{printf("%04d %s %s %s %s %s\n",$6,$5,$7,$3,$4,$1);}' RESULTS2/$pop1.$pop2.xpclr.$percentagekept.fst.$percentagekeptfst.selectionpositiongenesscores.txt | sort | uniq -c | awk '{print $3"\t"$5"\t"$6"\t"$7"\t"$4"\t"$2"\t"$1}' > RESULTS2/$pop1.$pop2.xpclr.$percentagekept.fst.$percentagekeptfst.selectionpositiongenes.txt #new
awk '{print $7}' RESULTS2/$pop1.$pop2.xpclr.$percentagekept.fst.$percentagekeptfst.selectionpositiongenesscores.txt | sort | grep -v "/" | uniq -c > toto.txt
set all = `awk 'FNR > 1 {OFS="\t"; num+=$1; print num}' toto.txt | tail -1`
awk -v a=$all '{print $1"\t"($1/a)*100"\t"$2"\t"a}' toto.txt > RESULTS2/$pop1.$pop2.xpclr.$percentagekept.fst.$percentagekeptfst.variantKind.txt
rm toto.txt
awk '{print $1}' RESULTS2/$pop1.$pop2.xpclr.$percentagekept.fst.$percentagekeptfst.selectionpositiongenes.txt | grep -v "-" | uniq > RESULTS2/$pop1.$pop2.xpclr.$percentagekept.fst.$percentagekeptfst.genes.txt

# cat ../vep.txt TEMPFILES/positionSSselection.txt | sort | awk 'NF==5{gene=$4;chr=$1;variant=$2;type=$5;next;}NF==4{if ($1==chr && $2==variant) {print $1,$2,$4, gene, type}next;}' > RESULTS2/$pop1.$pop2.selectionposition.fst.$percentagekeptfst.variantunderselection.txt
# awk '{print $4,$5}' RESULTS2/$pop1.$pop2.selectionposition.fst.$percentagekeptfst.variantunderselection.txt | uniq -c | awk '{print $2,$3,$1}' > RESULTS2/$pop1.$pop2.selectionposition.fst.$percentagekeptfst.genesvarianttypesunderselection.txt
# awk '{print $4}' RESULTS2/$pop1.$pop2.selectionposition.fst.$percentagekeptfst.variantunderselection.txt | uniq -c | awk '{print $2}' > RESULTS2/$pop1.$pop2.selectionposition.fst.$percentagekeptfst.genesunderselection.txt

echo "searching variants detected by xpclr using Fst $pop2 $pop1"
if (! (-e RESULTS2/$vcffile.allchr.$pop2.$pop1.highxpclrscore.$percentagekept.weir.fst) ) then
vcftools --gzvcf ../$folder/$vcffile.vcf.gz --bed RESULTS2/$vcffile.allchr.$pop2.$pop1.highxpclrscore.$percentagekept.bed --out TEMPFILES/positions --recode
vcftools --vcf TEMPFILES/positions.recode.vcf --weir-fst-pop POP/$pop2.txt --weir-fst-pop POP/$pop1.txt --out RESULTS2/$vcffile.allchr.$pop2.$pop1.highxpclrscore.$percentagekept
endif

echo "chr	position	fst" > RESULTS2/$pop2.$pop1.selectionposition.$percentagekeptfst.txt
grep -v "nan" RESULTS2/$vcffile.allchr.$pop2.$pop1.highxpclrscore.$percentagekept.weir.fst | awk -v a=$treshfstscore 'FNR>1 && $3>=a {print $1"\t"$2"\t"$3}' >> RESULTS2/$pop2.$pop1.selectionposition.fst.$percentagekeptfst.txt
awk 'FNR>1{printf ("%s %012d "1" %s\n", $1,$2,$3);}' RESULTS2/$pop2.$pop1.selectionposition.fst.$percentagekeptfst.txt > TEMPFILES/variantsS.txt    # new
awk '{printf ("%s %s "0" %s %s %s %s\n", $1,$2,$3,$4,$5,$6);}' RESULTS2/allchr.$pop2.$pop1.highxpclrscore.$percentagekept.grpwindows.txt > TEMPFILES/scoresS.txt    #new
cat TEMPFILES/scoresS.txt TEMPFILES/variantsS.txt | sort | awk 'NF==7{chr=$1;start=$2;end=$4;xpclr=$5;points=$6;rank=$7;next;}NF==4{if($1==chr && $2>=start && $2<=end){print chr"\t"$2"\t"xpclr"\t"points"\t"rank"\t"$4}next;}' | sort | uniq > RESULTS2/$pop2.$pop1.selectionposition.fst.$percentagekeptfst.xpclr.$percentagekept.txt   #new


awk -F" " '{printf("%s %012d 1 %s\n", $1,$2,$3);}' RESULTS2/$pop2.$pop1.selectionposition.fst.$percentagekeptfst.txt > TEMPFILES/positionSSselection.txt

awk -F" " '{printf("%s %s 1 %s %s %s %s\n", $1,$2,$3,$4,$5,$6);}' RESULTS2/$pop2.$pop1.selectionposition.fst.$percentagekeptfst.xpclr.$percentagekept.txt > TEMPFILES/positionSSselection2.txt #new
cat ../vep.txt TEMPFILES/positionSSselection2.txt | sort | awk 'NF==5{gene=$4;chr=$1;variant=$2;type=$5;next;}NF==7{if ($1==chr && $2==variant) {print $1,$2,$4,$5,gene,$6,type,$7}next;}' > RESULTS2/$pop2.$pop1.xpclr.$percentagekept.fst.$percentagekeptfst.selectionpositiongenesscores.txt # new
awk '{printf("%04d %s %s %s %s %s\n",$6,$5,$7,$3,$4,$1);}' RESULTS2/$pop2.$pop1.xpclr.$percentagekept.fst.$percentagekeptfst.selectionpositiongenesscores.txt | sort | uniq -c | awk '{print $3"\t"$5"\t"$6"\t"$7"\t"$4"\t"$2"\t"$1}' > RESULTS2/$pop2.$pop1.xpclr.$percentagekept.fst.$percentagekeptfst.selectionpositiongenes.txt #new
awk '{print $7}' RESULTS2/$pop2.$pop1.xpclr.$percentagekept.fst.$percentagekeptfst.selectionpositiongenesscores.txt | sort | grep -v "/" | uniq -c > toto.txt
set all = `awk 'FNR > 1 {OFS="\t"; num+=$1; print num}' toto.txt | tail -1`
awk -v a=$all '{print $1"\t"($1/a)*100"\t"$2"\t"a}' toto.txt > RESULTS2/$pop2.$pop1.xpclr.$percentagekept.fst.$percentagekeptfst.variantKind.txt
rm toto.txt

awk '{print $1}' RESULTS2/$pop2.$pop1.xpclr.$percentagekept.fst.$percentagekeptfst.selectionpositiongenes.txt | grep -v "-" | uniq > RESULTS2/$pop2.$pop1.xpclr.$percentagekept.fst.$percentagekeptfst.genes.txt

cp RESULTS2/$pop1.$pop2.xpclr.$percentagekept.fst.$percentagekeptfst.genes.txt RESULTS2/$para.xpclr.$percentagekept.fst.$percentagekeptfst.genes.txt
cat RESULTS2/$pop2.$pop1.xpclr.$percentagekept.fst.$percentagekeptfst.genes.txt >>  RESULTS2/$para.xpclr.$percentagekept.fst.$percentagekeptfst.genes.txt

#gzip -d -c $vep | grep -v "#" | awk '{print $2"\t"$4"\t"$7}' | awk -F":" '{print $1"\t"$2}' | awk '{printf("%s %012d 0 %s %s\n",$1,$2,$3,$4);}' > TEMPFILES/vep.txt
# cat ../vep.txt TEMPFILES/positionSSselection.txt | sort | awk 'NF==5{gene=$4;chr=$1;variant=$2;type=$5;next;}NF==4{if ($1==chr && $2==variant) {print $1,$2,$4, gene, type}next;}' > RESULTS2/$pop2.$pop1.selectionposition.fst.$percentagekeptfst.variantunderselection.txt
# awk '{print $4,$5}' RESULTS2/$pop2.$pop1.selectionposition.fst.$percentagekeptfst.variantunderselection.txt | uniq -c | awk '{print $2,$3,$1}' > RESULTS2/$pop2.$pop1.selectionposition.fst.$percentagekeptfst.genesvarianttypesunderselection.txt
# awk '{print $4}' RESULTS2/$pop2.$pop1.selectionposition.fst.$percentagekeptfst.variantunderselection.txt | uniq -c | awk '{print $2}' > RESULTS2/$pop2.$pop1.selectionposition.fst.$percentagekeptfst.genesunderselection.txt

rm TEMPFILES/vep*txt
rm TEMPFILES/positionSSselection*txt
rm TEMPFILES/score*txt
rm TEMPFILES/variant*txt
rm toto.txt

if (! (-e RESULTS2/$vcffile.allchr.$pop1.$pop2.fst.manhattan.pdf) ) then
echo "manhattan Fst"
echo "SNP	CHR	BP	P" > RESULTS2/$vcffile.allchr.$pop1.$pop2.fst.manhattanformat.txt
awk '$1>0 && $1<30 && $3>0.1{print "rs"$1"_"$2"\t"$1"\t"$2"\t"$3}' RESULTS2/$vcffile.$pop1.$pop2.nonan.weir.fst >> RESULTS2/$vcffile.allchr.$pop1.$pop2.fst.manhattanformat.txt
cp RESULTS2/$vcffile.allchr.$pop1.$pop2.fst.manhattanformat.txt TEMPFILES/manhattan.txt
R --no-save < RscriptManhattan.R
cp TEMPFILES/manhattan.pdf RESULTS2/$vcffile.allchr.$pop1.$pop2.fst.manhattan.pdf
endif

if (! (-e RESULTS2/$vcffile.allchr.$pop1.$pop2.xpclr.manhattan.pdf) ) then
echo "manhattan xpclr $pop1"
echo "SNP	CHR	BP	P" > RESULTS2/$vcffile.allchr.$pop1.$pop2.xpclr.manhattanformat.txt
awk '{print "rs"$1"_"$4"\t"$1"\t"$4"\t"$6}' RESULTS2/$vcffile.allchr.$pop1.$pop2.xpclr.out | sed 's/.000000//g' >> RESULTS2/$vcffile.allchr.$pop1.$pop2.xpclr.manhattanformat.txt
cp RESULTS2/$vcffile.allchr.$pop1.$pop2.xpclr.manhattanformat.txt TEMPFILES/manhattan.txt
R --no-save < RscriptManhattan.R
cp TEMPFILES/manhattan.pdf RESULTS2/$vcffile.allchr.$pop1.$pop2.xpclr.manhattan.pdf
endif

if (! (-e RESULTS2/$vcffile.allchr.$pop2.$pop1.xpclr.manhattan.pdf) ) then
echo "manhattan xpclr $pop2"
echo "SNP	CHR	BP	P" > RESULTS2/$vcffile.allchr.$pop2.$pop1.xpclr.manhattanformat.txt
awk '{print "rs"$1"_"$4"\t"$1"\t"$4"\t"$6}' RESULTS2/$vcffile.allchr.$pop2.$pop1.xpclr.out | sed 's/.000000//g' >> RESULTS2/$vcffile.allchr.$pop2.$pop1.xpclr.manhattanformat.txt
cp RESULTS2/$vcffile.allchr.$pop2.$pop1.xpclr.manhattanformat.txt TEMPFILES/manhattan.txt
R --no-save < RscriptManhattan.R
cp TEMPFILES/manhattan.pdf RESULTS2/$vcffile.allchr.$pop2.$pop1.xpclr.manhattan.pdf
endif


set topscoreNb = `wc-l RESULTS/$vcffile.allchr.$pop1.$pop2.highxpclrscore.$percentagekept.txt | awk '{print $1}'`
echo "chr	position	fst" > RESULTS/$pop1.$pop2.selectionposition.$fsttresh.txt
foreach j (`awk -v MAX=$topscoreNb 'BEGIN{for(i=1;i<=MAX;i++) print i;}'`)
	echo "chr	begin	end" > TEMPFILES/position.bed
	awk -v i=$j -v a=$winxpclr 'FNR==i{print $1"\t"$4-a"\t"$4+a}' RESULTS/$vcffile.allchr.$pop1.$pop2.highxpclrscore.$percentagekept.txt > TEMPFILES/position.bed
	vcftools --gzvcf $vcffile.vcf.gz --bed --out TEMPFILES/positions --recode
	vcftools --vcf TEMPFILES/positions.recode.vcf --weir-fst-pop $pop.1.txt --weir-fst-pop $pop.2.txt --out TEMPFILES/fst
	grep -v "nan" TEMPFILES/fst.weir.fst | awk -v a=$fsttresh 'FNR>1 && $3>=a {print $1"\t"$2"\t"$3}' >> RESULTS/$pop1.$pop2.selectionposition.fst.$fsttresh.txt
end

rm TEMPFILES/positions*
rm TEMPFILES/position.bed


awk -F" " '{printf("%s %012d 1 %s\n", $1,$4,$6);}' RESULTS/$vcffile.allchr.$pop1.$pop2.highxpclrscore.$percentagekept.txt > TEMPFILES/positionSSselection.txt
awk -F" " '{printf("%s %012d 1 %s\n", $1,$4,$6);}' RESULTS/$vcffile.allchr.$pop1.$pop2.highxpclrscore.$percentagekept.txt > TEMPFILES/positionSSselection.txt

awk -F" " '{printf("%s %012d 0 %012d %s\n",$1,$2-1500,$3+1500,$4);}' $annotation1 > TEMPFILES/annotation1.txt
awk -F" " '{printf("%s %012d 0 %012d %s\n",$1,$2-1500,$3+1500,$4);}' $annotation2 > TEMPFILES/annotation2.txt
cat TEMPFILES/annotation1.txt TEMPFILES/positionSSselection.txt | sort | awk 'NF==5{gene=$5;chr=$1;begin=$2;end=$4;next;}NF==4{if ($1==chr && $2<=end) {print $0, gene, end-begin-3000}next;}' > RESULTS/$vcffile.allchr.$pop1.$pop2.xpclrscore.$percentagekept.genesunderselection_plus.txt
cat TEMPFILES/annotation2.txt TEMPFILES/positionSSselection.txt | sort | awk 'NF==5{gene=$5;chr=$1;begin=$2;end=$4;next;}NF==4{if ($1==chr && $2<=end) {print $0, gene, end-begin-3000}next;}' > RESULTS/$vcffile.allchr.$pop1.$pop2.xpclrscore.$percentagekept.genesunderselection_minus.txt
mv RESULTS/$vcffile.allchr.$pop1.$pop2.xpclrscore.$percentagekept.genesunderselection_plus.txt RESULTS/$vcffile.allchr.$pop1.$pop2.xpclrscore.$percentagekept.genesunderselection_all.txt
cat RESULTS/$vcffile.allchr.$pop1.$pop2.xpclrscore.$percentagekept.genesunderselection_minus.txt >> RESULTS/$vcffile.allchr.$pop1.$pop2.xpclrscore.$percentagekept.genesunderselection_all.txt
rm RESULTS/$vcffile.allchr.$pop1.$pop2.xpclrscore.$percentagekept.genesunderselection_minus.txt
awk '{print $5,$6,$1,$4}' RESULTS/$vcffile.allchr.$pop1.$pop2.xpclrscore.$percentagekept.genesunderselection_all.txt | sort | uniq > RESULTS/$vcffile.xpclr.$pop1.$pop2.$percentagekept.genesunderselectionscores.txt
awk '{print $5,$6,$1}' RESULTS/$vcffile.allchr.$pop1.$pop2.xpclrscore.$percentagekept.genesunderselection_all.txt | sort | uniq -c > RESULTS/$vcffile.xpclr.$pop1.$pop2.$percentagekept.onlygenesunderselection_count.txt
awk '{print $5,$6,$1}' RESULTS/$vcffile.allchr.$pop1.$pop2.xpclrscore.$percentagekept.genesunderselection_all.txt | sort | uniq > RESULTS/$vcffile.xpclr.$pop1.$pop2.$percentagekept.onlygenesunderselection.txt	



cat RESULTS/$vcffile.xpclr.$pop1.$pop2.$percentagekept.onlygenesunderselection_count.txt | awk '{printf("%s %s %s 0 %s\n", $2,$3,$4,$1);}' | sort >  RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.onlygenesunderselection_geneoccurence.txt
cat RESULTS/$vcffile.xpclr.$pop1.$pop2.$percentagekept.genesunderselectionscores.txt RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.onlygenesunderselection_geneoccurence.txt | sort | awk 'NF==5{gene=$1;leng=$2;chr=$3;zero=$4;occur=$5;next;}NF==4{if ($1==gene) {print $0, occur}next;}' > RESULTS/$vcffile.xpclr.$pop1.$pop2.$percentagekept.genesunderselectionscores_occurence.txt
awk '{print $4,$1,$2,$3,$5}' RESULTS/$vcffile.xpclr.$pop1.$pop2.$percentagekept.genesunderselectionscores_occurence.txt | sort -n -r | awk '{printf ("%s %06d %s %s %s %s\n", $2,FNR,$1,$3,$4,$5);}' > RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.genesunderselection_scores_sorted.txt
awk '{print $1}' RESULTS/$vcffile.xpclr.$pop1.$pop2.$percentagekept.onlygenesunderselection.txt > RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.onlygenesunderselection_genenames.txt
cat RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.genesunderselection_scores_sorted.txt RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.onlygenesunderselection_genenames.txt | sort -n -r | awk 'NF==6{gene=$1;order=$2;score=$3;leng=$4;chr=$5;occur=$6;next;}NF==1{if ($1==gene) {print $0, order, score, leng,chr,occur}next;}' | awk '{print $2,$1,$3,$4,$5,$6}' | sort | awk '{print $2,$5,$6,$4,$4/$6,$3,$1}' > RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.onlygenesunderselection_scoresorted.txt
awk '{print $5,$1,$2,$3,$4,$6,$7}' RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.onlygenesunderselection_scoresorted.txt | sort -n | awk '{print $2,$3,$4,$5,$1,$6,$7}' > RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.onlygenesunderselection_occurencesorted.txt
awk '{print $1}' RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.onlygenesunderselection_scoresorted.txt > RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.onlygenesunderselection_scoresorted_onlygenenames.txt

awk -v a=$bp '$5<=a{print}' RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.onlygenesunderselection_scoresorted.txt > RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.onlygenesunderselection_less.$bp.scoresorted.txt
awk -v a=$bp '$5<=a{print}' RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.onlygenesunderselection_scoresorted.txt | awk '{print $1}' > RESULTS/$vcffile.xpclr.consensuspops.$percentagekept.onlygenesunderselection_less.$bp.scoresorted_onlygenenames.txt

