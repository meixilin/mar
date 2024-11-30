#!/bin/bash
#
#SBATCH --account=ac_moilab
#SBATCH --partition=savio3_htc
#SBATCH --time=23:59:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

set -eo pipefail

HOMEDIR="/global/scratch/projects/fc_moilab/meixilin/mar"
COMMITID=$(git --git-dir="${HOMEDIR}/.git" --work-tree="${HOMEDIR}/" rev-parse master)

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOBID}; GIT commit id ${COMMITID}; workdir = $(pwd)"

# download imputed matrix and check md5sum
wget -nv https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/1001_SNP_MATRIX.tar.gz
tar -xzvf 1001_SNP_MATRIX.tar.gz
cd 1001_SNP_MATRIX
md5sum --check md5sum.txt
cd ..

# download accessions
wget -O 1001_accessions.csv https://tools.1001genomes.org/api/accessions.csv?query=SELECT%20*%20FROM%20tg_accessions%20ORDER%20BY%20id

# run custom conversion scripts
module load python
python gm1001g_format.py

# gzip the genotype matrix
gzip 1001g_genotypes.txt

# create RDA file (use RDA)
module load r-spatial
Rscript --vanilla gm1001g_rda.R

# # clean up downloaded files
# rm -r 1001_SNP_MATRIX/
# rm 1001_accessions.csv
# rm 1001_SNP_MATRIX.tar.gz

# # rebuild the package
# R CMD INSTALL --preclean --no-multiarch --with-keep.source mar

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOBID} Done"
