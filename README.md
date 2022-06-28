# QTL


git clone https://github.com/PBGLMichaelHall/QTL.git

mamba env create --file env/QTLseqr.yaml 

conda activate QTLseqr

Rscript QTL.R -i ../../Desktop/freebayes_D2.filtered.vcf -H D2_F2_tt -l D2_F2_TT -c Chr.txt -r 0.20 -m 100 -d 400 -D 40 -G 99 -w 5000000 -f 0.1 -P F2 -B bulksize.txt -q 0.01 -a 0.01
