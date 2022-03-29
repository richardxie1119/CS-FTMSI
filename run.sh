#!/bin/bash
export PATH=/mnt/c/Users/Richard/Anaconda3/:%PATH:/mnt/d/Work/projects/MSI/subspaceMSI_py/CS-FTMSI
out_dir="coronal_40perc"
path_files="./file_dir_coronal_sparse.json"
ROI="R00"
basis_ROI="R00"
get_data_basis="False"
estimate_basis="False"
basis_path="./output/coronal_40perc/basis_transients.npy"
tissue_autoxcute_path="./Imaging_20210419_Coronal_sparse0.4.xml"
sample_autoxcute_path="./sampled_Imaging_20210419_Coronal_sparse0.4.xml"
n_basis_data=4000
lamb="0.1"
mu="1e-4"
max_iter=50
r=100
mz_range="400 1000"

python.exe create_data.py --out_dir $out_dir --path_file $path_files --tissue_autoxcute_path $tissue_autoxcute_path --sample_autoxcute_path $sample_autoxcute_path --ROI $ROI --get_data_basis $get_data_basis --basis_ROI $basis_ROI
sparse_subspace.exe --out_dir $out_dir --path_file $path_files --estimate_basis $estimate_basis --n_basis_data $n_basis_data --basis_path $basis_path --r $r --lamb $lamb --mu $mu --max_iter $max_iter
python.exe process_data.py --out_dir $out_dir --path_file $path_files  --r $r --mz_range $mz_range