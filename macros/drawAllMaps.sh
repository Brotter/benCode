for core in `seq 10 256`; do
    root drawAvgMaps.C\(${core}\) &
done