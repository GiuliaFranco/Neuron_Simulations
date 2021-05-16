#!/bin/bash

### All values of frequency from 1 to 50 Hz ###

#for i in 1000 200 100 66.6 50 40 33.3 28.7 25 22.2 20
#folder=WBSG
#mkdir ${folder}

for i in 33.3  
do
   ./a.out ${i} 7.96  &>  Gates_${i}.txt
   mv Latency.txt Latency_${i}.txt
   echo ${i}
done
