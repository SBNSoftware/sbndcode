#!/usr/bin/env python3

# setup fhiclpy v4_01_04 -q e20:prof

import fhicl



# raise Exception('Ahhhhhhhh')

name_reco = 'reco_sbnd.fcl'
name_g4 = 'standard_g4_sbnd.fcl'



pset = fhicl.make_pset(name_g4)
vis_g4 = pset['physics']['producers']['pdfastsim']['VISHits']
vuv_g4 = pset['physics']['producers']['pdfastsim']['VUVHits']

pset = fhicl.make_pset(name_reco)
vis_reco = pset['physics']['producers']['opt0finder']['VIVHits']
vuv_reco = pset['physics']['producers']['opt0finder']['VUVHits']

if (vis_g4 != vis_reco):
    raise Exception('VISHits')

if (vuv_g4 != vuv_reco):
    raise Exception('VISHits')

# import os, glob

# def find_fhich(name):
#     file_path = None

#     paths = os.environ['FHICL_FILE_PATH'].split(os.pathsep)
#     for p in paths:
#         file_path = os.path.join(p, name)
#         if(os.path.exists(file_path)):
#             print('Found:', file_path)
#             return file_path

# n = find_fhich(name)




