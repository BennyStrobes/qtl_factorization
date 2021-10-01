#SBATCH
#SBATCH --time=30:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


sldsc_input_data_dir="$1"
surge_interaction_eqtl_dir="$2"
sldsc_processed_data_dir="$3"
ldsc_source_code_dir="$4"

if false; then
python generate_sldsc_annotation_file_based_on_surge_eqtls.py $sldsc_input_data_dir $surge_interaction_eqtl_dir $sldsc_processed_data_dir

python generate_sldsc_annotation_file_based_on_surge_eqtls_and_base.py $sldsc_input_data_dir $surge_interaction_eqtl_dir $sldsc_processed_data_dir
fi
if false; then
python generate_sldsc_annotation_file_based_on_real_valued_surge_eqtls.py $sldsc_input_data_dir $surge_interaction_eqtl_dir $sldsc_processed_data_dir
fi

if false; then
for chrom_num in $(seq 1 $((22))); do 
	echo $chrom_num
	python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_data_dir}"surge_egenes_real_valued_05."${chrom_num}".annot" --out ${sldsc_processed_data_dir}"surge_egenes_real_valued_05."${chrom_num} --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"
done 

for chrom_num in $(seq 1 $((22))); do 
	echo $chrom_num
	python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_data_dir}"joint_surge_egenes_real_valued_05."${chrom_num}".annot" --out ${sldsc_processed_data_dir}"joint_surge_egenes_real_valued_05."${chrom_num} --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"
done 
fi






if false; then
for chrom_num in $(seq 1 $((22))); do 
	echo $chrom_num
	python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_data_dir}"joint_surge_egenes_w_base_05."${chrom_num}".annot" --out ${sldsc_processed_data_dir}"joint_surge_egenes_w_base_05."${chrom_num} --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"
done 


for chrom_num in $(seq 1 $((22))); do 
	echo $chrom_num
	python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_data_dir}"joint_surge_egenes_05."${chrom_num}".annot" --out ${sldsc_processed_data_dir}"joint_surge_egenes_05."${chrom_num} --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"
done 
fi



if false; then
for chrom_num in $(seq 1 $((22))); do 
	echo $chrom_num
	python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_data_dir}"surge_egenes_w_base_2."${chrom_num}".annot" --out ${sldsc_processed_data_dir}"surge_egenes_w_base_2."${chrom_num} --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"
done 
fi

if false; then
for chrom_num in $(seq 1 $((22))); do 
	echo $chrom_num
	python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_data_dir}"surge_egenes_2."${chrom_num}".annot" --out ${sldsc_processed_data_dir}"surge_egenes_2."${chrom_num} --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"
done 
fi


if false; then
for chrom_num in $(seq 1 $((22))); do 
	echo $chrom_num
	cp ${sldsc_input_data_dir}"baseline_v1.2/baseline."$chrom_num".annot.gz" ${sldsc_processed_data_dir}
	python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_data_dir}"baseline."$chrom_num".annot.gz" --out ${sldsc_processed_data_dir}"baseline."$chrom_num --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"
done 
fi
