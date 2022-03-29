# Joint Compressed Sensing and Subspace Modeling of the FT Mass Spectrometry Imaging Data

<p align="center">
  <img src="TOC_git.png" /width="700"> 
</p>

## What's included
This is the code repository containing the processing and algorithmic implementation as described in our paper **Enhancing the Throughput of FT Mass Spectrometry Imaging Using Joint Compressed Sensing and Subspace Modeling (Xie, YR., Castro D.C., Rubakhin, SS., Sweedler J.V., Lam, F.)** [[paper]](https://pubs.acs.org/doi/full/10.1021/acs.analchem.1c05279)

The repository contains:
- Preprocessing dependencies to prepare data sets from Bruker FT-ICR acquisitions and flexImaging generated files.
- A python executable binary for the joint CS and subspace reconstruction.
- Post analysis and processing pipeline to obtain hyperspectral data from reconstructed transients.

## How to run the demo
1. Download the repo. Be sure to manually download the `sparse_subspace.exe` from the repo (it won't automatically download due to the file size) and place it under the same project folder.

2. To download the demo data, follow the link: 

3. You should find `20210419_Coronal_HR.d`, `20210419_Coronal_sparse_40perc.d`, and `basis_transients.npy` within the zip file. Place the `.d` folders under the main folder, and `basis_transients.npy` under the `output/coronal_40perc` folder.




