#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH--mem=64G
#SBATCH -c 32
#SBATCH --time=12-00:00:00
##SBATCH --output=my.stdout
#SBATCH --job-name="FastStructure_SP"
#SBATCH -p batch

module load miniconda3
source activate faststructure

for i in 1 2 3 4 5 6 7 8 9 10
do
  python ~/.conda/envs/faststructure/bin/structure.py \
  -K $i \
  --input=Pc_only136.p2.f1SV.vcft_dp4_Q20_NoSingleton_noMD_NoRepeats.583k.STRUCTURE.v3 \
  --format=str --output=FS_SP_583k.K
done;

python ~/.conda/envs/faststructure/bin/chooseK.py \
--input=FS_SP_583k.K > FS_SP_583k_ChooseK.txt
