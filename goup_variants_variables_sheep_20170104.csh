#!/bin/csh -f
# script to group outliers varaiants in sheep from population based and sambada results. the output file has also the variable associated. The script has been made on 04/01/2017.

foreach var (bio7 alt temp7)
    set folder=`printf "%s_analysis" $var`  
    set pop1=`printf "low%s_20" $var`
    set pop2=`printf "high%s_20" $var`
    if (-e $folder/RESULTS2/$var.fst.xpclr.0.9995.txt) then
    	echo "$folder/RESULTS2/$var.fst.xpclr.0.9995.txt file already exist"
	else
    	echo "Creating $folder/RESULTS2/$var.fst.xpclr.0.9995.txt file"
	    awk '{printf("%s %012d\n",$1,$2);}' $folder/RESULTS2/$pop1.$pop2.xpclr.0.9999.fst.0.9995.selectionpositiongenesscores.txt > toto
    	awk '{printf("%s %012d\n",$1,$2);}' $folder/RESULTS2/$pop2.$pop1.xpclr.0.9999.fst.0.9995.selectionpositiongenesscores.txt >> toto
    	sort toto | uniq -c | awk '{print $2,$3}' > $folder/RESULTS2/$var.fst.xpclr.0.9995.txt
	endif
end

foreach var (slope)
    set folder=`printf "%s_analysis" $var`  
    set pop1=`printf "low%s20" $var`
    set pop2=`printf "high%s20" $var`
    if (-e $folder/RESULTS2/$var.fst.xpclr.0.9995.txt) then
    	echo "$folder/RESULTS2/$var.fst.xpclr.0.9995.txt file already exist"
	else
    	echo "Creating $folder/RESULTS2/$var.fst.xpclr.0.9995.txt file"
	    awk '{printf("%s %012d\n",$1,$2);}' $folder/RESULTS2/$pop1.$pop2.xpclr.0.9999.fst.0.9995.selectionpositiongenesscores.txt > toto
    	awk '{printf("%s %012d\n",$1,$2);}' $folder/RESULTS2/$pop2.$pop1.xpclr.0.9999.fst.0.9995.selectionpositiongenesscores.txt >> toto
    	sort toto | uniq -c | awk '{print $2,$3}' > $folder/RESULTS2/$var.fst.xpclr.0.9995.txt
	endif
end

foreach var (prec4 tmax4 bio15 bio3 ti216 ti2112 bio8 tmin8 tmean2 prec8)
    set folder=`printf "%s_analysis" $var`  
    set pop1=`printf "low%s" $var`
    set pop2=`printf "high%s" $var`
    if (-e $folder/RESULTS2/$var.fst.xpclr.0.9995.txt) then
    	echo "$folder/RESULTS2/$var.fst.xpclr.0.9995.txt file already exist"
	else
    	echo "Creating $folder/RESULTS2/$var.fst.xpclr.0.9995.txt file"
	    awk '{printf("%s %012d\n",$1,$2);}' $folder/RESULTS2/$pop1.$pop2.xpclr.0.9999.fst.0.9995.selectionpositiongenesscores.txt > toto
    	awk '{printf("%s %012d\n",$1,$2);}' $folder/RESULTS2/$pop2.$pop1.xpclr.0.9999.fst.0.9995.selectionpositiongenesscores.txt >> toto
    	sort toto | uniq -c | awk '{print $2,$3}' > $folder/RESULTS2/$var.fst.xpclr.0.9995.txt
	endif
end


echo "" > ../sheepalloutliervariant_variables.txt
foreach var (alt slope bio7 bio3 bio8 temp7 tmax4 tmin8 tmean2 ti216 ti2112 bio15 prec4 prec8)
    set folder=`printf "%s_analysis" $var`  
    awk -v a=$var '{print $1,$2,a}' $folder/RESULTS2/$var.fst.xpclr.0.9995.txt > toto
    grep "$var"  ../sambada/MOOA_QvalSign0.1.txt | awk -v a=$var '{printf("%s %012d %s\n",$2,$3,a);}' >> toto
	sort toto | uniq >> ../sheepalloutliervariant_variables.txt
end

sed 's/bio7/temp/g' ../sheepalloutliervariant_variables.txt | sed 's/bio3/temp/g' | sed 's/bio8/temp/g' | sed 's/temp7/temp/g' | sed 's/tmax4/temp/g' | sed 's/tmin8/temp/g' | sed 's/tmean2/temp/g' | sed 's/ti216/ti/g' | sed 's/ti2112/ti/g' | sed 's/bio15/prec/g' | sed 's/prec4/prec/g' | sed 's/prec8/prec/g' > ../sheepalloutliervariant_variables_groupedVariables.txt