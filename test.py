# -*- coding: utf-8 -*-

import os

from tgeconet import TGECONET

###
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: tGeCoNet: a framework for constructing (t)emporal (Ge)ne (Co)-expression (Net)works
# @url: https://github.com/pietrocinaglia/tgeconet
###

WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"
OUTPUT_PATH = WORKSPACE + 'results/use_cases/MAFLD/'

genes_of_interest = ['PNPLA3', 'TM6SF2', 'GCKR', 'MBOAT7', 'HSD17B13', 'TOR1B', 'COBLL1', 'GRB14', 'INSR', 'SREBF1']
tissues_of_interest = ['Liver','Adipose_Subcutaneous', 'Adipose_Visceral_Omentum','Kidney_Cortex','Kidney_Medulla','Small_Intestine_Terminal_Ileum']

print("################################################")
print("> TESTING <")
print("- Genes: " + str(genes_of_interest))
print("- Tissues: " + str(tissues_of_interest))
print("################################################")
print()

# Instantiate
tgeconet = TGECONET(
    genes_of_interest=genes_of_interest,
    tissues_of_interest=tissues_of_interest,
    threshold=0.05,
    verbose=True
)

# Build temporal network
temporal_network = tgeconet.construct_temporal_network()

# Statistics concerning the temporal network
stats = tgeconet.analyze_temporal_network(temporal_network, output_path=OUTPUT_PATH)

# Extract top (n) genes
# (if 'output_path' is defined, then each snapshot will be stored as image; it is not mandatory)
top_genes = tgeconet.extract_high_degree_genes(temporal_network, top_n=10, output_path=OUTPUT_PATH)

# Export temporal network as adjacency_matrices with pvalues (layer by layer)
tgeconet.export_adjacency_matrices(temporal_network, output_path=OUTPUT_PATH)

# Export temporal network by using snapshot-based representation as set of edgelist
tgeconet.save_as_files(temporal_network,output_path=OUTPUT_PATH)

# Plotting temporal network and storing plots as image
# (if 'output_path' is defined, then each snapshot will be stored as image; it is not mandatory)
tgeconet.plot(temporal_network, with_labels=True, output_path=OUTPUT_PATH)

print(">>> DONE <<<")
