module load gaussian/16-c.01
mkdir -p /scratch/$USER
g16 -scrdir=/scratch/$USER < qc.in > qc.out
