
# Capstone 3 project for Computational Biology at Auburn University

# initialize a file to store structural variant count data 
echo "Starting Part A"
echo "Doing the counts"
echo -e "Variant Type \t Count in File \t Filtered Count" > sv-counts.txt

count_INS=$(awk '$0!~/^#/' BD/bdmax9.ctx | awk -F "\t" '($7 ~ /INS/) {print $7}' | wc -l) 
count_DEL=$(awk '$0!~/^#/' BD/bdmax9.ctx | awk -F "\t" '($7 ~ /DEL/) {print $7}' | wc -l) 
count_INV=$(awk '$0!~/^#/' BD/bdmax9.ctx | awk -F "\t" '($7 ~ /INV/) {print $7}' | wc -l) 
count_ITX=$(awk '$0!~/^#/' BD/bdmax9.ctx | awk -F "\t" '($7 ~ /ITX/) {print $7}' | wc -l) 
count_CTX=$(awk '$0!~/^#/' BD/bdmax9.ctx | awk -F "\t" '($7 ~ /CTX/) {print $7}' | wc -l) 
# echo $count_INS
# echo $count_DEL
# echo $count_INV
# echo $count_ITX
# echo $count_CTX

# filtered counts

filtered_INS=$(awk '$0!~/^#/' BD/bdmax9.ctx | awk -F "\t" '($7 ~ /INS/ && $9>95 && $10>10) {print $7, $9}' | wc -l)
filtered_DEL=$(awk '$0!~/^#/' BD/bdmax9.ctx | awk -F "\t" '($7 ~ /DEL/ && $9>95 && $10>10) {print $7, $9}' | wc -l)
filtered_INV=$(awk '$0!~/^#/' BD/bdmax9.ctx | awk -F "\t" '($7 ~ /INV/ && $9>95 && $10>10) {print $7, $9}' | wc -l)
filtered_ITX=$(awk '$0!~/^#/' BD/bdmax9.ctx | awk -F "\t" '($7 ~ /ITX/ && $9>95 && $10>10) {print $7, $9}' | wc -l)
filtered_CTX=$(awk '$0!~/^#/' BD/bdmax9.ctx | awk -F "\t" '($7 ~ /CTX/ && $9>95 && $10>10) {print $7, $9}' | wc -l)
# echo $filtered_INS
# echo $filtered_DEL
# echo $filtered_INV
# echo $filtered_ITX
# echo $filtered_CTX

echo -e "Insertions \t $count_INS \t $filtered_INS" >> sv-counts.txt
echo -e "Deletions \t $count_DEL \t $filtered_DEL" >> sv-counts.txt
echo -e "Inversions \t $count_INV \t $filtered_INV" >> sv-counts.txt
echo -e "Translocations (ITX) \t $count_ITX \t $filtered_ITX" >> sv-counts.txt
echo -e "Translocations (CTX) \t $count_CTX \t $filtered_CTX" >> sv-counts.txt

echo "Finished with Counts"
#####################################################################################

echo "Calculating scores for Part A"
# Calculate the scores for table 2

scores_INS=$(awk '$0!~/^#/' BD/bdmax9.ctx | grep 'INS' | awk -F "\t" '($7 ~ /INS/) {sum+=$5-$2} ($7 ~ /INS/) {avgscore+=$9} ($7 ~ /INS/) {avgreads+=$10} END {print Insertions, avgscore/NR, sum/NR, avgreads/NR}')
scores_DEL=$(awk '$0!~/^#/' BD/bdmax9.ctx | grep 'DEL' |awk -F "\t" ' {sum+=$5-$2} {avgscore+=$9} {avgreads+=$10} END {print Deletions, avgscore/NR, sum/NR, avgreads/NR}')
scores_INV=$(awk '$0!~/^#/' BD/bdmax9.ctx | grep 'INV' |awk -F "\t" ' {sum+=$5-$2} {avgscore+=$9} {avgreads+=$10} END {print Inversions, avgscore/NR, sum/NR, avgreads/NR}')
scores_ITX=$(awk '$0!~/^#/' BD/bdmax9.ctx | grep 'ITX' |awk -F "\t" ' {sum+=$5-$2} {avgscore+=$9} {avgreads+=$10} END {print Translocations (ITX), avgscore/NR, sum/NR, avgreads/NR}')
scores_CTX=$(awk '$0!~/^#/' BD/bdmax9.ctx | grep 'CTX' |awk -F "\t" ' {sum+=$5-$2} {avgscore+=$9} {avgreads+=$10} END {print Translocations (CTX), avgscore/NR, sum/NR, avgreads/NR}')


# creating the second table 
echo -e "Variant Type \t Mean Score \t Mean Size \t Mean Number of Reads" > sv-scores.txt
echo -e  "Insertions \t $scores_INS" >> sv-scores.txt
echo -e  "Deletions \t $scores_DEL" >> sv-scores.txt
echo -e  "Inversions \t $scores_INV" >> sv-scores.txt
echo -e  "Translocations (ITX) \t $scores_ITX" >> sv-scores.txt
echo -e  "Translocations (CTX) \t $scores_CTX" >> sv-scores.txt

echo "Finished calculating scores. Part A complete."
##################################
## Part B - Restriction Sites GAATTC
##################################

echo "Starting Part B"

# cd capstone-projects/capstone-03/Buehler_Project3/

# awk '$0!~/>/' Enzyme/chrXV.fa

# Initialize a table to store restriction enzyme data in
echo -e "Restriction Enzyme \t Cleavage Site \t Number of Fragments \t Average Fragment Length" > enzymes.txt

chr=$(awk '$0!~/>/' Enzyme/chrXV.fa) 
# echo ${#chr}

ressites=(GAATTC GRCGYC GACNNNNNGTC)
# echo ${ressites[0]}
ressitenames=(EcoR1 Acyl Ahdl)
# echo ${ressitenames[@]}

echo $chr | sed 's/GAATTC/&\n/g' > ${ressitenames[0]}.txt
echo $chr | sed -E 's/G[AG]CG[CT]C/&\n/g' > ${ressitenames[1]}.txt
echo $chr | sed -E 's/GAC[ATCG]{5}GTC/&\n/g' > ${ressitenames[2]}.txt


# calculate the average length of an EcoR1 fragment

ecor1_length=($(awk '{ print length }' EcoR1.txt))
# echo ${ecor1_length[@]}
# echo ${#ecor1_length[@]}

eco_sum=0
for i in ${ecor1_length[@]}; do 
 let eco_sum+=$i
done
# echo $eco_sum

ecor1_meanlength=$(expr $eco_sum / ${#ecor1_length[@]})
# echo $ecor1_meanlength

# store the EcoR1 data into the file we initialized earlier
echo -e "${ressitenames[0]} \t ${ressites[0]} \t  ${#ecor1_length[@]} \t $ecor1_meanlength" >> enzymes.txt


# calculate the average lenght of an Acyl fragment

acyl_length=($(awk '{ print length }' Acyl.txt))
# echo ${acyl_length[@]}
# echo ${#acyl_length[@]}

acyl_sum=0
for i in ${acyl_length[@]}; do 
 let acyl_sum+=$i
done
# echo $acyl_sum

acyl_meanlength=$(expr $acyl_sum / ${#acyl_length[@]})
# echo $acyl_meanlength

echo -e "${ressitenames[1]} \t ${ressites[1]} \t  ${#acyl_length[@]} \t $acyl_meanlength" >> enzymes.txt


# calculate the average lenght of an ahdl fragment

ahdl_length=($(awk '{ print length }' Ahdl.txt))
# echo ${ahdl_length[@]}
# echo ${#ahdl_length[@]}

ahdl_sum=0
for i in ${ahdl_length[@]}; do 
 let ahdl_sum+=$i
done
# echo $ahdl_sum

ahdl_meanlength=$(expr $ahdl_sum / ${#ahdl_length[@]})
# echo $ahdl_meanlength

echo -e "${ressitenames[2]} \t ${ressites[2]} \t  ${#ahdl_length[@]} \t $ahdl_meanlength" >> enzymes.txt

echo "Finished Part B"

rm EcoR1.txt
rm Acyl.txt
rm Ahdl.txt