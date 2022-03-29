import numpy as np
import pandas as pd
import os
from os import path
from glob import glob
import sys
from utils import *
from processing import *
import pickle
import random
import argparse
import re
from scipy.stats import median_abs_deviation as mad
import json


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Reconstruct and process the reconstructed data.'
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
        '--r', type = str,
        help='r.'
        )
    parser.add_argument(
        '--mz_range', required = True, type = int, nargs = '+',
        help='Two values that defines the lower and upper bound of m/z range for analyzing spectra.'
        )

    args = parser.parse_args()

    with open(args.path_file, 'r') as fp:
        path_dict = json.load(fp)

    params = getParams(path_dict['parameter_data_path'])
    params_basis = getParams(path_dict['parameter_basis_path'])

    V_hat = np.load('./output/'+args.out_dir+'/basis_transients.npy')

    Ui = np.load('output/{}/updatedU.npy'.format(args.out_dir))
    mz = params_basis['m'][(params_basis['m']>args.mz_range[0])&(params_basis['m']<args.mz_range[1])]
    average_sp = np.zeros(mz.shape)

    r = int(args.r)

    print('calculating average spectra...')
    for i in tqdm(range(0, Ui.shape[0],100)):
        fid = Ui[i:i+100].dot(V_hat[:r])
        mz,sp = fid2spec(fid, params_basis['m'], args.mz_range)
        average_sp += np.sum(sp,0)
    
    average_sp = average_sp/Ui.shape[0]

    np.save('output/{}/avgspec_recon.npy'.format(args.out_dir),average_sp)

    mzs_use = peak_detection(mz,average_sp,5*mad(average_sp), 5*mad(average_sp))['mz']

    print('propogating intensity values...')

    peak_data = []
    for i in tqdm(range(0, Ui.shape[0], 100)):
        fid = Ui[i:i+100].dot(V_hat[:r])
        mz,sp = fid2spec(fid, params_basis['m'], args.mz_range)
        peak_data.append(sp[:,np.in1d(mz, mzs_use)])
    peak_data = np.concatenate(peak_data)

    with open('output/{}/peakdata_recon.pkl'.format(args.out_dir), 'wb') as f:
        pickle.dump({'intens_mtx':peak_data,'mz':mzs_use}, f, protocol=pickle.HIGHEST_PROTOCOL)



