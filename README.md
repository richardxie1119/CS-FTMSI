# Joint Compressed Sensing and Subspace Modeling of the FT Mass Spectrometry Imaging Data

<p align="center">
  <img src="TOC_git.png" /width="700"> 
</p>

## What's included
This is the code repository containing the processing and algorithmic implementation as described in our paper **Enhancing the Throughput of FT Mass Spectrometry Imaging Using Joint Compressed Sensing and Subspace Modeling. Xie, YR., Castro D.C., Rubakhin, SS., Sweedler J.V., Lam, F., Anal. Chem., doi:10.1021/acs.analchem.1c05279** [[paper]](https://pubs.acs.org/doi/full/10.1021/acs.analchem.1c05279)

The repository contains:
- Preprocessing dependencies to prepare data sets from Bruker FT-ICR acquisitions and flexImaging generated files.
- A Pyhon executable binary for the joint CS and subspace reconstruction.
- Post analysis and processing pipeline to obtain hyperspectral data from reconstructed transients.


## Run the program
Install all dependencies (numpy >=1.19.5, scipy>=1.16.3, h5py). The Python binary executable currently only supports running on Windows system.
In the bash script `run.sh`, the following arguments need to be defined.
```
export PATH=/mnt/c/Users/Richard/Anaconda3/:%PATH:/mnt/d/Work/projects/MSI/subspaceMSI_py/CS-FTMSI
out_dir **where the output files will be stored in the output folder**
path_files="./file_dir_coronal_sparse.json" #file that indicates the required data paths
ROI="R00" #ROI header for the data imaging file.
basis_ROI="R00" #ROI header for the basis imaging file.
get_data_basis="False" #if to prepare data for the basis estimation. set to True will require .ser file to sample the long transients for basis estimation.
estimate_basis="False" #set to True to estimate the basis transients if no predefined basis file is available.

#if both get_data_basis and estimate_basis are False, then the program will look for basis transients file in the output folder.
basis_path="./output/coronal_40perc/basis_transients.npy" 
tissue_autoxcute_path="./Imaging_20210419_Coronal_sparse0.4.xml" #autoxcute file generated by Bruker.
sample_autoxcute_path="./sampled_Imaging_20210419_Coronal_sparse0.4.xml" #sampled autoxcute file for acquisition.
n_basis_data=4000 #number of transients loaded for the basis estimation. only if get_data_basis and estimate_basis set to True.
lamb="0.1" #lambda parameter for ADMM.
mu="1e-4" #mu parameter for ADMM.
max_iter=50 #maximum number of iteration for the optimization algorithm.
r=100 #number of basis transients used for reconstruction.
mz_range="400 1000" #m/z range for the data processing.
```


## How to run the demo
1. Download the repo. Be sure to manually download the `sparse_subspace.exe` from the repo (it won't automatically download due to the file size) and place it under the same project folder.

2. To download the demo data, follow the link: 

3. You should find `20210419_Coronal_HR.d`, `20210419_Coronal_sparse_40perc.d`, and `basis_transients.npy` within the zip file. Place the `.d` folders under the main folder, and `basis_transients.npy` under the `output/coronal_40perc` folder. `.d` folders are data folders generated by Bruker Solarix and flexImaging, which contain the experimental parameters and raw transients. 




