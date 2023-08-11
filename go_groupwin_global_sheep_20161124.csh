#!/bin/csh -f                                                                                                                                                           

set chrnb = $argv[1]
set vcffile = $argv[2]
#set pop1 = $argv[3]
#set pop2 = $argv[4]
set percentagekept = $argv[3]
set winxpclr = $argv[4]
set percentagekeptfst = $argv[5]
#set var = $argv[8]

if (-e ../groupwindows) then
    echo "../groupwindows directory already exist"
else
    echo "Creating ../groupwindows directory"
    mkdir ../groupwindows
endif


foreach var (bio7 alt temp7)
    set folder=`printf "%s_analysis" $var`  
    set pop1=`printf "low%s_20" $var`
    set pop2=`printf "high%s_20" $var`
	gawk '{printf("%s %d %d %s\n",$1,$2,$3,$4);}' $folder/RESULTS2/allchr.$pop2.$pop1.highxpclrscore.$percentagekept.grpwindows.txt > TEMPFILES/tmp.txt
	gawk '{printf("%s %d %d %s\n",$1,$2,$3,$4);}' $folder/RESULTS2/allchr.$pop1.$pop2.highxpclrscore.$percentagekept.grpwindows.txt >> TEMPFILES/tmp.txt
	sort TEMPFILES/tmp.txt > TEMPFILES/tmp1.txt
	mv TEMPFILES/tmp1.txt TEMPFILES/tmp.txt
	awk -f overlapingwin.awk TEMPFILES/tmp.txt TEMPFILES/tmp.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto.txt
	awk -f overlapingwin.awk toto.txt toto.txt | awk 'FNR>1{printf ("%s %012d %012d %012s\n",$1,$5,$6,$4);}' | sort > TEMPFILES/scores.txt
	awk -f overlapingwin.awk toto.txt toto.txt | awk 'FNR>1{printf ("%s %012d %012d\n",$1,$5,$6);}' | sort | uniq -c | awk '{printf ("%s %012d %012d "1" %s\n", $2,$3,$4,$1);}' > TEMPFILES/gridpoints.txt
	cat TEMPFILES/scores.txt TEMPFILES/gridpoints.txt | sort | awk 'NF==4{chr=$1;start=$2;end=$3;score=$4;next;}NF==5{if($1==chr&&$2==start&&$3==end){print chr,start,end,score,$5}next;}' | awk '{printf("%012s %s %s %s %s\n", $4,$1,$2,$3,$5);}' | sort -r | awk '{print $2,$3,$4,$1,$5,NR}' > $folder/RESULTS2/grpwindows.$percentagekept.txt
	cat $folder/RESULTS2/$pop1.$pop2.selectionposition.fst.$percentagekeptfst.txt > TEMPFILES/tmp1.txt
	cat $folder/RESULTS2/$pop2.$pop1.selectionposition.fst.$percentagekeptfst.txt >> TEMPFILES/tmp1.txt
	sort TEMPFILES/tmp1.txt | uniq > TEMPFILES/tmp2.txt
	awk 'FNR>1{printf ("%s %012d "1" %s\n", $1,$2,$3);}' TEMPFILES/tmp2.txt > TEMPFILES/variantsS.txt    # new                                                              
	awk '{printf ("%s %s "0" %s %s %s %s\n", $1,$2,$3,$4,$5,$6);}' $folder/RESULTS2/grpwindows.$percentagekept.txt > TEMPFILES/scoresS.txt    #new                                           
	cat TEMPFILES/scoresS.txt TEMPFILES/variantsS.txt | sort | awk 'NF==7{chr=$1;start=$2;end=$4;xpclr=$5;points=$6;rank=$7;next;}NF==4{if($1==chr && $2>=start && $2<=end) {print chr"\t"$2"\t"xpclr"\t"points"\t"rank"\t"$4}next;}' | sort | uniq > $folder/RESULTS2/grpWin.fst.$percentagekeptfst.xpclr.$percentagekept.txt                      
                                                                                                                                                                        

	echo "" > $folder/RESULTS2/$var.grpRegions.$percentagekeptfst.txt
	foreach rank (`awk '{print $5}' $folder/RESULTS2/grpWin.fst.$percentagekeptfst.xpclr.$percentagekept.txt | sort | uniq | awk '{printf $1" "}'`)
    	awk -v a=$rank '$6==a{printf("%s %d %d 1\n",$1,$2,$3);}' $folder/RESULTS2/grpwindows.$percentagekept.txt >> $folder/RESULTS2/$var.grpRegions.$percentagekeptfst.txt
	end
    echo "chr	start	end	score" > TEMPFILES/toto.txt
    grep "$var"  ../sambada/MOOA_QvalSign0.1.txt | gawk '{printf ("%s %d %d %1.10Lf\n",$2,$3-50000,$3+50000,$7);}' >> TEMPFILES/toto.txt
	awk -f overlapingwin.awk TEMPFILES/toto.txt TEMPFILES/toto.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto.txt
	awk -f overlapingwin.awk toto.txt toto.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto1.txt
	awk -f overlapingwin.awk toto1.txt toto1.txt | awk 'FNR>1{printf ("%s %d %d 1\n",$1,$5,$6);}' | sort -r | uniq > TEMPFILES/scores.txt
	cat $folder/RESULTS2/$var.grpRegions.$percentagekeptfst.txt >> TEMPFILES/scores.txt
    awk -f overlapingwin.awk TEMPFILES/scores.txt TEMPFILES/scores.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto1.txt
    awk -f overlapingwin.awk toto1.txt toto1.txt | awk 'FNR>1{printf ("%s %d %d 1\n",$1,$5,$6);}' | sort -r | uniq > ../groupwindows/$var.grpwinOutliers.allmeth.$percentagekeptfst.txt
end

foreach var (slope)
    set folder=`printf "%s_analysis" $var`  
    set pop1=`printf "low%s20" $var`
    set pop2=`printf "high%s20" $var`
	gawk '{printf("%s %d %d %s\n",$1,$2,$3,$4);}' $folder/RESULTS2/allchr.$pop2.$pop1.highxpclrscore.$percentagekept.grpwindows.txt > TEMPFILES/tmp.txt
	gawk '{printf("%s %d %d %s\n",$1,$2,$3,$4);}' $folder/RESULTS2/allchr.$pop1.$pop2.highxpclrscore.$percentagekept.grpwindows.txt >> TEMPFILES/tmp.txt
	sort TEMPFILES/tmp.txt > TEMPFILES/tmp1.txt
	mv TEMPFILES/tmp1.txt TEMPFILES/tmp.txt
	awk -f overlapingwin.awk TEMPFILES/tmp.txt TEMPFILES/tmp.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto.txt
	awk -f overlapingwin.awk toto.txt toto.txt | awk 'FNR>1{printf ("%s %012d %012d %012s\n",$1,$5,$6,$4);}' | sort > TEMPFILES/scores.txt
	awk -f overlapingwin.awk toto.txt toto.txt | awk 'FNR>1{printf ("%s %012d %012d\n",$1,$5,$6);}' | sort | uniq -c | awk '{printf ("%s %012d %012d "1" %s\n", $2,$3,$4,$1);}' > TEMPFILES/gridpoints.txt
	cat TEMPFILES/scores.txt TEMPFILES/gridpoints.txt | sort | awk 'NF==4{chr=$1;start=$2;end=$3;score=$4;next;}NF==5{if($1==chr&&$2==start&&$3==end){print chr,start,end,score,$5}next;}' | awk '{printf("%012s %s %s %s %s\n", $4,$1,$2,$3,$5);}' | sort -r | awk '{print $2,$3,$4,$1,$5,NR}' > $folder/RESULTS2/grpwindows.$percentagekept.txt
	cat $folder/RESULTS2/$pop1.$pop2.selectionposition.fst.$percentagekeptfst.txt > TEMPFILES/tmp1.txt
	cat $folder/RESULTS2/$pop2.$pop1.selectionposition.fst.$percentagekeptfst.txt >> TEMPFILES/tmp1.txt
	sort TEMPFILES/tmp1.txt | uniq > TEMPFILES/tmp2.txt
	awk 'FNR>1{printf ("%s %012d "1" %s\n", $1,$2,$3);}' TEMPFILES/tmp2.txt > TEMPFILES/variantsS.txt    # new                                                              
	awk '{printf ("%s %s "0" %s %s %s %s\n", $1,$2,$3,$4,$5,$6);}' $folder/RESULTS2/grpwindows.$percentagekept.txt > TEMPFILES/scoresS.txt    #new                                           
	cat TEMPFILES/scoresS.txt TEMPFILES/variantsS.txt | sort | awk 'NF==7{chr=$1;start=$2;end=$4;xpclr=$5;points=$6;rank=$7;next;}NF==4{if($1==chr && $2>=start && $2<=end) {print chr"\t"$2"\t"xpclr"\t"points"\t"rank"\t"$4}next;}' | sort | uniq > $folder/RESULTS2/grpWin.fst.$percentagekeptfst.xpclr.$percentagekept.txt                      
                                                                                                                                                                        

	echo "" > $folder/RESULTS2/$var.grpRegions.$percentagekeptfst.txt
	foreach rank (`awk '{print $5}' $folder/RESULTS2/grpWin.fst.$percentagekeptfst.xpclr.$percentagekept.txt | sort | uniq | awk '{printf $1" "}'`)
    	awk -v a=$rank '$6==a{printf("%s %d %d 1\n",$1,$2,$3);}' $folder/RESULTS2/grpwindows.$percentagekept.txt >> $folder/RESULTS2/$var.grpRegions.$percentagekeptfst.txt
	end
    echo "chr	start	end	score" > TEMPFILES/toto.txt
    grep "$var"  ../sambada/MOOA_QvalSign0.1.txt | gawk '{printf ("%s %d %d %1.10Lf\n",$2,$3-50000,$3+50000,$7);}' >> TEMPFILES/toto.txt
	awk -f overlapingwin.awk TEMPFILES/toto.txt TEMPFILES/toto.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto.txt
	awk -f overlapingwin.awk toto.txt toto.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto1.txt
	awk -f overlapingwin.awk toto1.txt toto1.txt | awk 'FNR>1{printf ("%s %d %d 1\n",$1,$5,$6);}' | sort -r | uniq > TEMPFILES/scores.txt
	cat $folder/RESULTS2/$var.grpRegions.$percentagekeptfst.txt >> TEMPFILES/scores.txt
    awk -f overlapingwin.awk TEMPFILES/scores.txt TEMPFILES/scores.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto1.txt
    awk -f overlapingwin.awk toto1.txt toto1.txt | awk 'FNR>1{printf ("%s %d %d 1\n",$1,$5,$6);}' | sort -r | uniq > ../groupwindows/$var.grpwinOutliers.allmeth.$percentagekeptfst.txt
end

foreach var (prec4 tmax4 bio15 bio3 ti216 ti2112 bio8 tmin8 tmean2 prec8)
    set folder=`printf "%s_analysis" $var`
    set pop1=`printf "low%s" $var`
    set pop2=`printf "high%s" $var`
	gawk '{printf("%s %d %d %s\n",$1,$2,$3,$4);}' $folder/RESULTS2/allchr.$pop2.$pop1.highxpclrscore.$percentagekept.grpwindows.txt > TEMPFILES/tmp.txt
	gawk '{printf("%s %d %d %s\n",$1,$2,$3,$4);}' $folder/RESULTS2/allchr.$pop1.$pop2.highxpclrscore.$percentagekept.grpwindows.txt >> TEMPFILES/tmp.txt
	sort TEMPFILES/tmp.txt > TEMPFILES/tmp1.txt
	mv TEMPFILES/tmp1.txt TEMPFILES/tmp.txt
	awk -f overlapingwin.awk TEMPFILES/tmp.txt TEMPFILES/tmp.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto.txt
	awk -f overlapingwin.awk toto.txt toto.txt | awk 'FNR>1{printf ("%s %012d %012d %012s\n",$1,$5,$6,$4);}' | sort > TEMPFILES/scores.txt
	awk -f overlapingwin.awk toto.txt toto.txt | awk 'FNR>1{printf ("%s %012d %012d\n",$1,$5,$6);}' | sort | uniq -c | awk '{printf ("%s %012d %012d "1" %s\n", $2,$3,$4,$1);}' > TEMPFILES/gridpoints.txt
	cat TEMPFILES/scores.txt TEMPFILES/gridpoints.txt | sort | awk 'NF==4{chr=$1;start=$2;end=$3;score=$4;next;}NF==5{if($1==chr&&$2==start&&$3==end){print chr,start,end,score,$5}next;}' | awk '{printf("%012s %s %s %s %s\n", $4,$1,$2,$3,$5);}' | sort -r | awk '{print $2,$3,$4,$1,$5,NR}' > $folder/RESULTS2/grpwindows.$percentagekept.txt
	cat $folder/RESULTS2/$pop1.$pop2.selectionposition.fst.$percentagekeptfst.txt > TEMPFILES/tmp1.txt
	cat $folder/RESULTS2/$pop2.$pop1.selectionposition.fst.$percentagekeptfst.txt >> TEMPFILES/tmp1.txt
	sort TEMPFILES/tmp1.txt | uniq > TEMPFILES/tmp2.txt
	awk 'FNR>1{printf ("%s %012d "1" %s\n", $1,$2,$3);}' TEMPFILES/tmp2.txt > TEMPFILES/variantsS.txt    # new                                                              
	awk '{printf ("%s %s "0" %s %s %s %s\n", $1,$2,$3,$4,$5,$6);}' $folder/RESULTS2/grpwindows.$percentagekept.txt > TEMPFILES/scoresS.txt    #new                                           
	cat TEMPFILES/scoresS.txt TEMPFILES/variantsS.txt | sort | awk 'NF==7{chr=$1;start=$2;end=$4;xpclr=$5;points=$6;rank=$7;next;}NF==4{if($1==chr && $2>=start && $2<=end) {print chr"\t"$2"\t"xpclr"\t"points"\t"rank"\t"$4}next;}' | sort | uniq > $folder/RESULTS2/grpWin.fst.$percentagekeptfst.xpclr.$percentagekept.txt                      
                                                                                                                                                                        

	echo "" > $folder/RESULTS2/$var.grpRegions.$percentagekeptfst.txt
	foreach rank (`awk '{print $5}' $folder/RESULTS2/grpWin.fst.$percentagekeptfst.xpclr.$percentagekept.txt | sort | uniq | awk '{printf $1" "}'`)
    	awk -v a=$rank '$6==a{printf("%s %d %d 1\n",$1,$2,$3);}' $folder/RESULTS2/grpwindows.$percentagekept.txt >> $folder/RESULTS2/$var.grpRegions.$percentagekeptfst.txt
	end
    echo "chr	start	end	score" > TEMPFILES/toto.txt
    grep "$var"  ../sambada/MOOA_QvalSign0.1.txt | gawk '{printf ("%s %d %d %1.10Lf\n",$2,$3-50000,$3+50000,$7);}' >> TEMPFILES/toto.txt
	awk -f overlapingwin.awk TEMPFILES/toto.txt TEMPFILES/toto.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto.txt
	awk -f overlapingwin.awk toto.txt toto.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto1.txt
	awk -f overlapingwin.awk toto1.txt toto1.txt | awk 'FNR>1{printf ("%s %d %d 1\n",$1,$5,$6);}' | sort -r | uniq > TEMPFILES/scores.txt
	cat $folder/RESULTS2/$var.grpRegions.$percentagekeptfst.txt >> TEMPFILES/scores.txt
    awk -f overlapingwin.awk TEMPFILES/scores.txt TEMPFILES/scores.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto1.txt
    awk -f overlapingwin.awk toto1.txt toto1.txt | awk 'FNR>1{printf ("%s %d %d 1\n",$1,$5,$6);}' | sort -r | uniq > ../groupwindows/$var.grpwinOutliers.allmeth.$percentagekeptfst.txt
end

echo "" > TEMPFILES/toto.txt
foreach var (alt slope bio7 temp7 prec4 tmax4 bio15 bio3 ti216 ti2112 bio8 tmin8 tmean2 prec8)
	cat ../groupwindows/$var.grpwinOutliers.allmeth.$percentagekeptfst.txt >> TEMPFILES/toto.txt
end
awk -f overlapingwin.awk TEMPFILES/toto.txt TEMPFILES/toto.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto.txt
awk -f overlapingwin.awk toto.txt toto.txt | awk '{print $1"\t"$5"\t"$6"\t"$4}' > toto1.txt
awk -f overlapingwin.awk toto1.txt toto1.txt | awk 'FNR>1{printf ("%s %d %d 1\n",$1,$5,$6);}' | sort | uniq > ../groupwindows/grpwinOutliers.allmeth.allvar.$percentagekeptfst.txt
