while read p
do
n=$(echo $p | xargs -n 1 basename | awk -F'_L002' '{print $1}')
sbatch -J ${n}_lane2 --mem=50G --time=1-24:00 --cpus-per-task 12 --mail-user=rohit.satyam@kaust.edu.sa --mail-type=FAIL --partition=batch -o ${n}.out -e ${n}.err --wrap="nextflow run main.nf --input '$p' --outdir lane1_results --ref /ibex/scratch/projects/c2077/rohit/rahuls_rnaseq_march2024/RNAgrinder/resources/ToxoDB-67_TgondiiME49_Genome.fasta --gtf /ibex/scratch/projects/c2077/rohit/rahuls_rnaseq_march2024/RNAgrinder/resources/ToxoDB-67_TgondiiME49.gtf --mode PE --rrnaUse ribodetector --index_dir /ibex/scratch/projects/c2077/rohit/rahuls_rnaseq_march2024/RNAgrinder/resources/00_index --skipKraken=true --cpus 12 -w ${n}_work"
done < file
