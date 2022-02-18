#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
echo -e "n \t prop \t size \t nc \t maxc \t pc \t n'c \t distmin \t duration" >> output_tests
for n in 10 15 20 25 30 35 40 45 50 75 100 
do
  for prop in  0.05 0.10 0.15 0.20 0.25
  do
    for size in 50 60 70 80 90 100 125 150
    do
      maxConf="$n*($n-1)/2"
      nc=$(bc -l <<< "scale=2; $prop*$maxConf")
      nc=$( printf "%.0f" $nc)
      min_maxc=$(bc -l <<< "scale=2; 4*($nc/$n)")
      min_maxc=$( printf "%.0f" $min_maxc)
      for step in 1 2 3 4 5
      do
        let maxc=$min_maxc+$step
        printf "$n\t$prop\t$size\t$nc\t$maxc" >> output_tests
        ./generator -mode 14 -n $n -nc $nc -maxc $maxc -w $size -h $size -a $size >> output_tests
      done
    done
  done
done
