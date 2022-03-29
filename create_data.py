import numpy as np
import os
from os import path
from glob import glob
import sys
from utils import *
from processing import *
import pickle
import argparse
import json
import h5py

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Create and prepare data for reconstruction.'
        )

    parser.add_argument(
        '--out_dir', required = True, type = str,
        help='An output directory to store the processed results.'
        )
    parser.add_argument(
        '--path_file', required = True, type = str,
        help='Path to the .json file that specifies the raw .ser and imaginginfo files.'
        )
    parser.add_argument(
        '--tissue_autoxcute_path', type = str,
        help='tissue_autoxcute_path.'
        )
    parser.add_argument(
        '--sample_autoxcute_path', type = str,
        help='sample_autoxcute_path.'
        )
    parser.add_argument(
        '--ROI', nargs = '?', const = 'False', default = 'R00', type = str)
    parser.add_argument(
        '--get_data_basis', nargs = '?', const = 'False', default = 'False', type = str, choices = ['True','False'],
        help='If process the raw data.'
        )
    parser.add_argument(
        '--basis_ROI', nargs = '?', const = 'False', default = 'R00', type = str)
    parser.add_argument(
        '--n_basis_data', type = int, default = 4000,
        help='n_basis_data.'
        )

    args = parser.parse_args()

    print('-'*40, 'script parameters','-'*40)
    print(args)

    if not os.path.isdir('output/'+args.out_dir):
        os.makedirs('output/'+args.out_dir)

    with open(args.path_file, 'r') as fp:
        path_dict = json.load(fp)

    params = getParams(path_dict['parameter_data_path'])
    params_basis = getParams(path_dict['parameter_basis_path'])

    sampled_coords = parse_coords(args.sample_autoxcute_path)
    tissue_coords = parse_coords(args.tissue_autoxcute_path)

    np.save('./output/'+args.out_dir+'/sampled_coords.npy',sampled_coords)
    np.save('./output/'+args.out_dir+'/tissue_coords.npy',tissue_coords)

    data_imaging_info = parseImagingInfo(path_dict['data_imaging_info'])
    basis_imaging_info = parseImagingInfo(path_dict['basis_imaging_info'])

    if args.get_data_basis == 'True':
        print('loading transients for basis estimation...')
        
        fid_idx_sampled = np.array(random.sample(list(basis_imaging_info[args.basis_ROI]['scan_index']),args.n_basis_data))
        fid_idx_sampled = fid_idx_sampled[np.argsort(fid_idx_sampled)]

        with h5py.File('./output/'+args.out_dir+"/basis_data.h5", "w") as data_file:
            for i in tqdm(range(fid_idx_sampled.shape[0])):
                fid = loadBrukerFIDs(path_dict['basis_data_path'], params_basis['parameters']['TD'],
                    params_basis['parameters']['TD'], fid_idx_sampled[i])
                group = data_file.create_group(str(i))
                group.create_dataset('transient',data=fid)


    print('loading available transients...')

    with h5py.File('./output/'+args.out_dir+"/data.h5", "w") as data_file:
        for i in tqdm(range(sampled_coords.shape[0])):
            fid = loadBrukerFIDs(path_dict['data_path'], params['parameters']['TD'],
                params['parameters']['TD'], data_imaging_info[args.ROI]['scan_index'][i])
            group = data_file.create_group(str(i))
            group.create_dataset('transient',data=fid)


