########
#
#  In case you want to parse the final event rate from the log files to plot
#
#######


for file in `ls`; do
    core=`echo ${file} | cut -d"." -f1`
    rateStr=`tail -3 ${file} | head -1 | cut -d" " -f4`
    rate=${rateStr:0:3}

    echo ${core} ${rate} >> parseRate.txt
done