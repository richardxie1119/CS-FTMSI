import os
from os import path
import numpy as np
import xml.etree.ElementTree as ET
import re
from tqdm import tqdm
import pickle
from processing import *

def loadBrukerFIDs(file_path, fid_length, read_length, fid_idx, verbose = False):
    """

    """
    fids = []
    if path.exists(file_path):

        f = open(file_path,'r')

        if type(fid_idx) == list or type(fid_idx) == np.ndarray:

            #print('loading {} FID from file...'.format(len(fid_idx)))

            for i in range(len(fid_idx)):

                if verbose:
                    print('loading {} FID from file...'.format(i))

                f.seek(4*(fid_idx[i]-1)*fid_length) #seek FID locations within the .ser file

                if read_length == 'all':
                    fid = np.fromfile(f, count = fid_length, dtype = 'int32')
                    fids.append(fid)
                else:
                    fid = np.fromfile(f, count = read_length, dtype = 'int32')
                    fids.append(fid)
            
        else:

            f.seek(4*(fid_idx-1)*fid_length) #seek FID locations within the .ser file

            if read_length == 'all':
                fid = np.fromfile(f, count = fid_length, dtype = 'int32')
                fids.append(fid)
            else:
                fid = np.fromfile(f, count = read_length, dtype = 'int32')
                fids.append(fid)

        f.close()

    else:
        raise Exception('ser file does not exist in the provided file path. please double check.')

    return np.array(fids,dtype='float64')


def loadBrukerMethod(file_path):

    """TODO"""

    return 'A'


def parseImagingInfo(file_path):

    """parses the ImagingInfo.xml file in the imaging .d folder, and returns the relative
    coordinates for each imaged regions, starting with RXX. The parsed dictionary contains
    the arrays of relative coordinates Xs and Ys under keys named as the regions (RXX).

    """

    tree = ET.parse(file_path)
    root = tree.getroot()

    parsed_spots = {}
    spotNames = []
    scan = []
    TIC = []
    for child in root:
        spotNames.append(child.find('spotName').text)
        scan.append(child.find('count').text)
        TIC.append(child.find('tic').text)

    ROI = set([spot[:3] for spot in spotNames])

    for roi in ROI:

        coord = []
        scan_idx = []
        tic = []

        for i in range(len(spotNames)):
            spot = spotNames[i]

            if roi in spot:
                coord.append([int(re.search('X(.*)Y' ,spot).group(1)),
                int(re.search('Y(.*)' ,spot).group(1))])
                scan_idx.append(scan[i])
                tic.append(TIC[i])
                
        scan_idx = np.array(scan_idx, dtype='int64')
        tic = np.array(tic, dtype='float')

        coord = np.array(coord)
        coord[:,0] -= coord[:,0].min()
        coord[:,1] -= coord[:,1].min()

        parsed_spots[roi] = {'coordinates':coord, 'scan_index':scan_idx, 'tic':tic}

    return parsed_spots


def parseBrukerMethod(file_path):

    tree = ET.parse(file_path)
    root = tree.getroot()

    for type_tag in root.findall('paramlist')[0]:
        name = type_tag.get('name')
        if name == 'SW_h':
            SW_h = float(type_tag.findall('value')[0].text)
        if name == 'TD':
            TD = int(type_tag.findall('value')[0].text)
        if name == 'ML1':
            ML1 = float(type_tag.findall('value')[0].text)
        if name == 'ML2':
            ML2 = float(type_tag.findall('value')[0].text)
        if name == 'ML3':
            ML3 = float(type_tag.findall('value')[0].text)

    return {'SW_h':SW_h,'TD':TD,'ML1':ML1,'ML2':ML2,'ML3':ML3}


def getParams(method_file_path):

    """
    """
    parameters = parseBrukerMethod(method_file_path)

    T = 1/(parameters['SW_h']*2)
    t = np.arange(0, parameters['TD'])*T
    f = parameters['SW_h']*2*np.arange(0, parameters['TD']/2+1)/parameters['TD']
    m = fticr_mass_axis(f, [parameters['ML1'],parameters['ML2'],parameters['ML3']])

    return {'parameters':parameters,'T':T,'t':t,'f':f,'m':m}


def parse_coords(file_dir):

    tree = ET.parse(file_dir)
    root = tree.getroot()

    coords = []
    for child in root:
        point = child.attrib['Pos_on_Scout']
        coords.append([int(re.search('X(.*)Y' ,point).group(1)),int(re.search('Y(.*)' ,point).group(1))])

    coords = np.array(coords)
    coords[:,0] -= coords[:,0].min()
    coords[:,1] -= coords[:,1].min()

    return coords
