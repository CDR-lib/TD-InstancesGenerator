#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
echo "n \t prop \t size \t nc \t maxc \n" >> output_tests
for n in 10 15 20 25 30 35 40 45 50 75 100 
do
  for prop in  0.05 0.10 0.15 0.20 0.25
  do
    for size in 125 150 175 200 225 250 275 300   
    do
      maxConf="$n*($n-1)/2"
      nc=$(bc -l <<< "scale=2; $prop*$maxConf")
      nc=$( printf "%.0f" $nc)
      min_maxc=$(bc -l <<< "scale=2; 4*($nc/$n)")
      min_maxc=$( printf "%.0f" $min_maxc)
      for step in 1 2 3 4 5
      do
        let maxc=$min_maxc+$step
        printf "\n--------------------------------------------------------------------------------------------------------------------\n" >> output_tests
        printf "$n\t$prop\t$size\t$nc\t$maxc" >> output_tests
        ./generator -mode 6 -n $n -nc $nc -maxc $maxc -w $size -h $size >> output_tests
      done
    done
  done
done
