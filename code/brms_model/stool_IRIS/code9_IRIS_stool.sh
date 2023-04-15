#!/bin/bash

#of sbatch options.
#SBATCH --job-name=brms_stool_cytokine_9
#SBATCH --nodes=1
#SBATCH --mem=20000
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --partition=nih_s10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xzhou7@stanford.edu
#SBATCH --account=mpsnyder
#SBATCH --time=167:00:00
 
#dir是过滤之后的数据存放路径
module load R
Rscript --vanilla code9_IRIS_stool.R
