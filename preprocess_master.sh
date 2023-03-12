antibodys=(...)
for antibody in ${antibodys[*]}
do
    if [ ! -d $antibody ]; then mkdir $antibody; fi 
    bash preprocess_day.sh ${antibody} &
    
done
