#job name
#$ -N Can_megalondon_call_GPU_5hmC
#name of logfile 
#$ -o Can_megalondon_call_GPU_5hmC.log
# Set the working directory 
#$ -wd /SAN/mottlab/ArabidopsisAssembly/Methylation_analysis/Ziming_analysis/Nanopore_analysis_all
#memory and runtime options 
#$ -l tmem=30G
# -l h_vmem=50G
#$ -l h_rt=150:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l gpu=true
#$ -pe gpu 2
#$ -R y
# -t 1-15 
# -tc 2 
# -l tscratch=10G

PATH=/share/apps/mottlab/Guppy/ont-guppy/bin/:$PATH
PATH=/share/apps/mottlab/miniconda3/bin/:$PATH
PATH=/share/apps/cuda-11.4/bin/:$PATH

conda activate megalodon

echo $(date)

#megalodon /SAN/mottlab/ArabidopsisAssembly/Methylation_analysis/Nanopore_data_all/Col/ --guppy-config dna_r9.4.1_450bps_sup_prom.cfg --remora-modified-bases dna_r9.4.1_e8 sup 0.0.0 5mc CG 0 --outputs mod_mappings mods --output-directory Col_megalodon --reference ./reference/Col_final_23022022_new_hapog_patch_hapog_reorder.fasta --devices 0 1 --processes 4 --guppy-params "--num_callers 3 --ipc_threads 5" --guppy-server-path /share/apps/mottlab/Guppy/ont-guppy/bin/guppy_basecall_server --overwrite --guppy-timeout 300 

#megalodon /SAN/mottlab/ArabidopsisAssembly/Methylation_analysis/Nanopore_data_all/Can/ --guppy-config dna_r9.4.1_450bps_sup_prom.cfg --remora-modified-bases dna_r9.4.1_e8 sup 0.0.0 5mc CG 0 --outputs mod_mappings mods --output-directory Can_megalodon --reference ./reference/Can_final_23022022_new_hapog.fasta --devices 0 1 --processes 4 --guppy-params "--num_callers 3 --ipc_threads 5" --guppy-server-path /share/apps/mottlab/Guppy/ont-guppy/bin/guppy_basecall_server --overwrite --guppy-timeout 300 

#megalodon /SAN/mottlab/ArabidopsisAssembly/Methylation_analysis/Nanopore_data_all/Can/ --guppy-config dna_r9.4.1_450bps_sup_prom.cfg --remora-modified-bases dna_r9.4.1_e8 sup 0.0.0 5mc CG 0 --outputs mod_mappings mods --output-directory Can_megalodon_Col_reference --reference ./reference/Col_final_23022022_new_hapog_patch_hapog_reorder.fasta --devices 0 1 --processes 4 --guppy-params "--num_callers 3 --ipc_threads 5" --guppy-server-path /share/apps/mottlab/Guppy/ont-guppy/bin/guppy_basecall_server --overwrite --guppy-timeout 300 


megalodon /SAN/mottlab/ArabidopsisAssembly/Methylation_analysis/Nanopore_data_all/Can/ --guppy-config dna_r9.4.1_450bps_sup_prom.cfg --remora-modified-bases dna_r9.4.1_e8 sup 0.0.0 5hmc_5mc CG 0 --outputs mod_mappings mods --output-directory Can_megalodon_5hmC_5mC --reference ./reference/Can_final_23022022_new_hapog.fasta --devices 0 1 --processes 4 --guppy-params "--num_callers 3 --ipc_threads 5" --guppy-server-path /share/apps/mottlab/Guppy/ont-guppy/bin/guppy_basecall_server --overwrite --guppy-timeout 300

#megalodon /SAN/mottlab/ArabidopsisAssembly/Methylation_analysis/Nanopore_data_all/Col/ --guppy-config dna_r9.4.1_450bps_sup_prom.cfg --remora-modified-bases dna_r9.4.1_e8 sup 0.0.0 5hmc_5mc CG 0 --outputs mod_mappings mods --output-directory Col_megalodon_5hmC_5mC --reference ./reference/Col_final_23022022_new_hapog_patch_hapog_reorder.fasta --devices 0 1 --processes 4 --guppy-params "--num_callers 3 --ipc_threads 3" --guppy-server-path /share/apps/mottlab/Guppy/ont-guppy/bin/guppy_basecall_server --overwrite --guppy-timeout 300

#megalodon /SAN/mottlab/ArabidopsisAssembly/Methylation_analysis/Nanopore_data_all/Can/ --guppy-config dna_r9.4.1_450bps_sup_prom.cfg --remora-modified-bases dna_r9.4.1_e8 sup 0.0.0 5hmc_5mc CG 0 --outputs mod_mappings mods --output-directory Can_megalodon_Col_reference_5hmC_5mC --reference ./reference/Col_final_23022022_new_hapog_patch_hapog_reorder.fasta --devices 0 1 --processes 4 --guppy-params "--num_callers 3 --ipc_threads 5" --guppy-server-path /share/apps/mottlab/Guppy/ont-guppy/bin/guppy_basecall_server --overwrite --guppy-timeout 300

