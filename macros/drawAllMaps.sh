startCore=${1}
stopCore=${2}
for core in `seq ${startCore} ${stopCore}`; do
    root -b -q drawAvgMaps.C\(${core}\) &> ${core}.log
done
