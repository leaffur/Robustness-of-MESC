chr=$1
simul=$2
h_g=$3
N=$4
Ne=$5
dir=$6
h_med_fix=$7

module load PLINK/1.9b_6.21-x86_64; module load R/4.2.0-foss-2020b;


mkdir -p $dir
echo "sim_gcta.R"
if test -f "$dir/temp/sumstat_hm0.2_t100.sumstats.gz"
	then 
		echo "Step02 has done"
	else
		 Rscript /gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/s1_samplesize/sim_gcta.R $chr $simul $h_g $N $Ne $dir $h_med_fix > $dir/sim_gcta.Rout	
	fi
		echo "Step02 done!"

echo "step03_scores.R"




if test -f "$dir/expr_mesc.$chr.expscore.gz"
	then 
		echo "Step03 has done"
	else
		Rscript /gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/s1_samplesize/step03_scores.R $chr $dir 2>&1 >/dev/null

		module purge;
		module load miniconda; conda activate ldsc; sh $dir/ldscore.sh
		echo "ldscore done"

		sh $dir/exprscore_ldsc.sh
		echo "exprscore_ldsc done"

		conda deactivate; conda activate mesc
		sh $dir/exprscore_mesc.sh
		echo "exprscore_mesc done"
	fi
		echo "Step03 done!"



module purge;
module load R/4.2.0-foss-2020b;
echo "run_ldsc_mesc.R"
Rscript /gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/s1_samplesize/run_ldsc_mesc.R $chr $simul $h_g $dir $h_med_fix 2>&1 >/dev/null

module purge; module load miniconda; conda activate mesc
sh $dir/mesc.sh
echo "mesc done"

