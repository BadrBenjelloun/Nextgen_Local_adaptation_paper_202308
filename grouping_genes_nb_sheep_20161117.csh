#!/bin/csh -f

# script to group popbased results and to make summaries # it needs arguments *ptgkept for xpclr and $ptgkeptfst for fst
set ptgkept = $argv[1]
set ptgkeptfst = $argv[2]

if (-e PopbRes.$ptgkeptfst) then
    echo "PopbRes.$ptgkeptfst directory already exist"
else
    echo "Creating PopbRes.$ptgkeptfst directory"
    mkdir PopbRes.$ptgkeptfst
endif

echo "" > PopbRes.$ptgkeptfst/All_outliergenes_popbased.$ptgkeptfst.txt
echo "" > PopbRes.$ptgkeptfst/summary.genenb_var.$ptgkeptfst.txt
echo "" > PopbRes.$ptgkeptfst/All_outliervariants_popbased.$ptgkeptfst.txt
echo "" > PopbRes.$ptgkeptfst/summary.variantnb_var.$ptgkeptfst.txt
foreach var (bio7 alt temp7)
    set folder=`printf "%s_analysis" $var`  
    set pop1=`printf "low%s_20" $var`
    set pop2=`printf "high%s_20" $var`
    if (-e PopbRes.$ptgkeptfst/$folder) then
    	echo "PopbRes.$ptgkeptfst/$folder directory already exist"
	else
    	echo "Creating PopbRes.$ptgkeptfst/$folder directory"
    	mkdir PopbRes.$ptgkeptfst/$folder
	endif
    echo "" > PopbRes.$ptgkeptfst/$folder/$var.genes.$ptgkeptfst.txt
    echo "" > PopbRes.$ptgkeptfst/$folder/$var.variants.$ptgkeptfst.txt    
    awk '{print $5}' $folder/RESULTS2/$pop1.$pop2.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt | sort | uniq | grep -v ^- > tmp.txt
    awk '{print $5}' $folder/RESULTS2/$pop2.$pop1.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt | sort | uniq | grep -v ^- >> tmp.txt
    cp $folder/RESULTS2/$pop1.$pop2.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt PopbRes.$ptgkeptfst/$folder/$pop1.$pop2.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt
    cp $folder/RESULTS2/$pop2.$pop1.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt PopbRes.$ptgkeptfst/$folder/$pop2.$pop1.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt
    sort tmp.txt | uniq >> PopbRes.$ptgkeptfst/All_outliergenes_popbased.$ptgkeptfst.txt
    sort tmp.txt | uniq >> PopbRes.$ptgkeptfst/$folder/$var.genes.$ptgkeptfst.txt
    echo "$var" >> PopbRes.$ptgkeptfst/summary.genenb_var.$ptgkeptfst.txt
    sort tmp.txt | uniq | wc -l >> PopbRes.$ptgkeptfst/summary.genenb_var.$ptgkeptfst.txt
#    rm tmp.txt
    awk '{print $1,$2}' $folder/RESULTS2/$pop1.$pop2.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt | sort | uniq > tmp.txt
    awk '{print $1,$2}' $folder/RESULTS2/$pop2.$pop1.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt | sort | uniq >> tmp.txt
    sort tmp.txt | uniq >> PopbRes.$ptgkeptfst/All_outliervariants_popbased.$ptgkeptfst.txt
    sort tmp.txt | uniq >> PopbRes.$ptgkeptfst/$folder/$var.variants.$ptgkeptfst.txt
    echo "$var" >> PopbRes.$ptgkeptfst/summary.variantnb_var.$ptgkeptfst.txt
    sort tmp.txt | uniq | wc -l >> PopbRes.$ptgkeptfst/summary.variantnb_var.$ptgkeptfst.txt
end
foreach var (slope)
    set folder=`printf "%s_analysis" $var`  
    set pop1=`printf "low%s20" $var`
    set pop2=`printf "high%s20" $var`
    if (-e PopbRes.$ptgkeptfst/$folder) then
    	echo "PopbRes.$ptgkeptfst/$folder directory already exist"
	else
    	echo "Creating PopbRes.$ptgkeptfst/$folder directory"
    	mkdir PopbRes.$ptgkeptfst/$folder
	endif    
    echo "" > PopbRes.$ptgkeptfst/$folder/$var.genes.$ptgkeptfst.txt
    echo "" > PopbRes.$ptgkeptfst/$folder/$var.variants.$ptgkeptfst.txt    
    awk '{print $5}' $folder/RESULTS2/$pop1.$pop2.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt | sort | uniq | grep -v ^- > tmp.txt
    awk '{print $5}' $folder/RESULTS2/$pop2.$pop1.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt | sort | uniq | grep -v ^- >> tmp.txt
    cp $folder/RESULTS2/$pop1.$pop2.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt PopbRes.$ptgkeptfst/$folder/$pop1.$pop2.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt
    cp $folder/RESULTS2/$pop2.$pop1.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt PopbRes.$ptgkeptfst/$folder/$pop2.$pop1.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt
    sort tmp.txt | uniq >> PopbRes.$ptgkeptfst/All_outliergenes_popbased.$ptgkeptfst.txt
    sort tmp.txt | uniq >> PopbRes.$ptgkeptfst/$folder/$var.genes.$ptgkeptfst.txt
    echo "$var" >> PopbRes.$ptgkeptfst/summary.genenb_var.$ptgkeptfst.txt
    sort tmp.txt | uniq | wc -l >> PopbRes.$ptgkeptfst/summary.genenb_var.$ptgkeptfst.txt
#   rm tmp.txt
    awk '{print $1,$2}' $folder/RESULTS2/$pop1.$pop2.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt | sort | uniq > tmp.txt
    awk '{print $1,$2}' $folder/RESULTS2/$pop2.$pop1.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt | sort | uniq >> tmp.txt
    sort tmp.txt | uniq >> PopbRes.$ptgkeptfst/All_outliervariants_popbased.$ptgkeptfst.txt
    sort tmp.txt | uniq >> PopbRes.$ptgkeptfst/$folder/$var.variants.$ptgkeptfst.txt
    echo "$var" >> PopbRes.$ptgkeptfst/summary.variantnb_var.$ptgkeptfst.txt
    sort tmp.txt | uniq | wc -l >> PopbRes.$ptgkeptfst/summary.variantnb_var.$ptgkeptfst.txt
end
foreach var (prec4 tmax4 bio15 bio3 ti216 ti2112 bio8 tmin8 tmean2 prec8)
    set folder=`printf "%s_analysis" $var`
    set pop1=`printf "low%s" $var`
    set pop2=`printf "high%s" $var`
    if (-e PopbRes.$ptgkeptfst/$folder) then
    	echo "PopbRes.$ptgkeptfst/$folder directory already exist"
	else
    	echo "Creating PopbRes.$ptgkeptfst/$folder directory"
    	mkdir PopbRes.$ptgkeptfst/$folder
	endif
    echo "" > PopbRes.$ptgkeptfst/$folder/$var.genes.$ptgkeptfst.txt
    echo "" > PopbRes.$ptgkeptfst/$folder/$var.variants.$ptgkeptfst.txt    
    awk '{print $5}' $folder/RESULTS2/$pop1.$pop2.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt | sort | uniq | grep -v ^- > tmp.txt
    awk '{print $5}' $folder/RESULTS2/$pop2.$pop1.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt | sort | uniq | grep -v ^- >> tmp.txt
    cp $folder/RESULTS2/$pop1.$pop2.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt PopbRes.$ptgkeptfst/$folder/$pop1.$pop2.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt
    cp $folder/RESULTS2/$pop2.$pop1.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt PopbRes.$ptgkeptfst/$folder/$pop2.$pop1.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt
    sort tmp.txt | uniq >> PopbRes.$ptgkeptfst/All_outliergenes_popbased.$ptgkeptfst.txt
    sort tmp.txt | uniq >> PopbRes.$ptgkeptfst/$folder/$var.genes.$ptgkeptfst.txt
    echo "$var" >> PopbRes.$ptgkeptfst/summary.genenb_var.$ptgkeptfst.txt
    sort tmp.txt | uniq | wc -l >> PopbRes.$ptgkeptfst/summary.genenb_var.$ptgkeptfst.txt
#    rm tmp.txt
    awk '{print $1,$2}' $folder/RESULTS2/$pop1.$pop2.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt | sort | uniq > tmp.txt
    awk '{print $1,$2}' $folder/RESULTS2/$pop2.$pop1.xpclr.$ptgkept.fst.$ptgkeptfst.selectionpositiongenesscores.txt | sort | uniq >> tmp.txt
    sort tmp.txt | uniq >> PopbRes.$ptgkeptfst/All_outliervariants_popbased.$ptgkeptfst.txt
    sort tmp.txt | uniq >> PopbRes.$ptgkeptfst/$folder/$var.variants.$ptgkeptfst.txt
    echo "$var" >> PopbRes.$ptgkeptfst/summary.variantnb_var.$ptgkeptfst.txt
    sort tmp.txt | uniq | wc -l >> PopbRes.$ptgkeptfst/summary.variantnb_var.$ptgkeptfst.txt
end


