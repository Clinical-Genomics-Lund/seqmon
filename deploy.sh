# Generate and copy plots to seqmon
cp illuminaseqstat.R /data/bnf/proj/seqmon
cp plotLoad.R /data/bnf/proj/seqmon
Rscript /data/bnf/proj/seqmon/illuminaseqstat.R
Rscript /data/bnf/proj/seqmon/plotLoad.R

# Deploy web page
scp www/seqmon.html pi@10.0.224.47:/home/pi/seqmon/
scp www/animation.css pi@10.0.224.47:/home/pi/seqmon/
