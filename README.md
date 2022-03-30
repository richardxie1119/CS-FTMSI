# Joint Compressed Sensing and Subspace Modeling of the FT Mass Spectrometry Imaging Data

<p align="center">
  <img src="https://github.com/richardxie1119/CS-FTMSI/blob/master/TOC_git.png" /width="600"> 
</p>

## What's included
This is the code repository containing the processing and algorithmic implementation as described in the paper:
- [ ]**Enhancing the Throughput of FT Mass Spectrometry Imaging Using Joint Compressed Sensing and Subspace Modeling. Xie, YR., Castro D.C., Rubakhin, SS., Sweedler J.V., Lam, F., Anal. Chem. doi:10.1021/acs.analchem.1c05279** [[paper]](https://pubs.acs.org/doi/full/10.1021/acs.analchem.1c05279)


The repository contains:
- Preprocessing dependencies to prepare data sets from Bruker FT-ICR acquisitions and flexImaging generated files.
- A Pyhon executable binary for the joint CS and subspace reconstruction.
- Post analysis and processing pipeline to obtain hyperspectral data from reconstructed transients.

<p align="center">
  <img src="https://github.com/richardxie1119/CS-FTMSI/blob/master/workflow.png" /width="600"> 
</p>

## Run the program
Install all dependencies (numpy >=1.19.5, scipy>=1.16.3, h5py). The Python binary executable currently only supports running on Windows system.
Within the bash script `run.sh`, the following arguments need to be defined.

`export PATH=PATH_TO_PYTHON:%PATH:PATH_TO_PROJECT`: define the paths to Python and project folder\
`out_dir` :where the output files will be stored in the output folder\
`path_files` :file that indicates the required data paths\
`ROI` :ROI header for the data imaging file.\
`basis_ROI` :ROI header for the basis imaging file.\
`get_data_basis` :if to prepare data for the basis estimation. set to True will require .ser file to sample the long transients for basis estimation.\
`estimate_basis` :set to True to estimate the basis transients if no predefined basis file is available.\
`basis_path` :if both get_data_basis and estimate_basis are False, then the program will look for basis transients file in the output folder.\
`tissue_autoxcute_path` :autoxcute file generated by Bruker. within the file, relative coordinates for each transient acquisition can be found and parsed.\
`sample_autoxcute_path` :sampled autoxcute file for acquisition.\
`n_basis_data` :number of transients loaded for the basis estimation. only if get_data_basis and estimate_basis set to True.\
`lamb` :lambda parameter for ADMM.\
`mu` :mu parameter for ADMM.\
`max_iter` :maximum number of iteration for the optimization algorithm.\
`r` :number of basis transients used for reconstruction.\
`mz_range` :m/z range for the data processing.\

In the bash script, three commands are run for the whole pipeline.\
Create data sets:
```
python.exe create_data.py --out_dir $out_dir --path_file $path_files --tissue_autoxcute_path $tissue_autoxcute_path --sample_autoxcute_path $sample_autoxcute_path --ROI $ROI --get_data_basis $get_data_basis --basis_ROI $basis_ROI
```
Run binary for CS and subspace reconstruction on the prepared data:
```
sparse_subspace.exe --out_dir $out_dir --path_file $path_files --estimate_basis $estimate_basis --n_basis_data $n_basis_data --basis_path $basis_path --r $r --lamb $lamb --mu $mu --max_iter $max_iter
```
Process the reconstructed data:
```
python.exe process_data.py --out_dir $out_dir --path_file $path_files  --r $r --mz_range $mz_range
```

## How to run the demo
1. Download the repo. Be sure to manually download the `sparse_subspace.exe` from the repo (it won't automatically download due to the file size) and place it under the same project folder.
2. To download the demo data, follow the link [[demo data]](https://uofi.box.com/s/dkip85acls48owqbn4oxkhd8tnymjzbw)
3. You should find `20210419_Coronal_HR.d`, `20210419_Coronal_sparse_40perc.d`, and `basis_transients.npy` within the zip file. Place the `.d` folders under the main folder, and place `basis_transients.npy` under the `output/coronal_40perc` folder. `.d` folders are data folders generated by Bruker Solarix and flexImaging, which contain the experimental parameters and raw transients. Due to the large size of data for basis estimation, basis transients were precomputed and provided for the reconstruction.
4. In command prompt, run ```bash run.sh```

## Analyze the data
Some post analysis can be found in the [[notebook]](https://github.com/richardxie1119/CS-FTMSI/blob/master/demo/demo_40perc.ipynb)

Visualize single-pixel reconstructed spectra:
```
idx = 62310
fid_recon = spatial_coef[idx].dot(V_hat)
mz_recon,sp_recon = fid2spec(fid_recon,params_basis['m'],mz_range=(700,900))

mz_filter = (mz_recon>865)&(mz_recon<875)
plt.figure(figsize=(10,5))
plt.plot(mz_recon[mz_filter], sp_recon[mz_filter],c='k')
plt.xlabel('m/z')
plt.ylabel('intensity')
```
Visualize reconstructed ion images:
```
mz_idx = [341,344]
fig,axes = plt.subplots(1,2,figsize=(8,4))
ax = axes.ravel()
for i,idx in enumerate(mz_idx):
    ion_img = IonImg(peak_data_recon['intens_mtx'][:,idx]/peak_data_recon['intens_mtx'].sum(1),tissue_coords,True,False)
    ax[i].imshow(ion_img,'hot')
    ax[i].set_title('Ion image for m/z {}'.format(np.round(peak_data_recon['mz'][idx],5)))
plt.show()
```



