import numpy as np
import pandas as pd
import os 
from util import read_json_from

PATH_ROOT = os.getcwd().replace("\\","/").replace("/notebooks","")

proteomeHD_path = f"{PATH_ROOT}/data_sources/ProteomeHD/ProteomeHD_v1_1.csv"
proteomeHD_df = pd.read_csv(proteomeHD_path)
proteomeHD_simplified_protein_ids = proteomeHD_df["Simplified_protein_ID"].to_numpy()
proteomeHD_majority_protein_ids = proteomeHD_df["Majority_protein_IDs"].to_numpy()
proteomeHD_feature_matrix = proteomeHD_df.iloc[:,4:].fillna(0).to_numpy()
proteomeHD_feature_matrix_with_na = proteomeHD_df.iloc[:,4:].to_numpy()
major_simplified_idx_lookup_path = f"{PATH_ROOT}/data_sources/ProteomeHD/major_simplified_to_idx_lookup.json"
major_simplified_idx_lookup = read_json_from(major_simplified_idx_lookup_path)

pQTL_protein_path = f"{PATH_ROOT}/data_sources/pQTL/pQTL_protein_converted.csv"
pQTL_protein_df = pd.read_csv(pQTL_protein_path)
pQTL_protein_ids = pQTL_protein_df['uniprotswissprot'].to_numpy()
pQTL_protein_feature_matrix = pQTL_protein_df.iloc[:,2:].fillna(0).to_numpy()
pQTL_protein_feature_matrix_with_na = pQTL_protein_df.iloc[:,2:].to_numpy()
pQTL_protein_idx_lookup_path = f"{PATH_ROOT}/data_sources/pQTL/pQTL_protein_converted_idx_lookup.json"
pQTL_protein_idx_lookup = read_json_from(pQTL_protein_idx_lookup_path)

protein_consensus_path = f"{PATH_ROOT}/data_sources/Tissue/protein_consensus_converted.csv"
protein_consensus_df = pd.read_csv(protein_consensus_path)
protein_consensus_ids = protein_consensus_df['uniprotswissprot'].to_numpy()
protein_consensus_feature_matrix_with_na = protein_consensus_df.iloc[:,2:].to_numpy()

pQTL_rna_seq_path = f"{PATH_ROOT}/data_sources/pQTL/rna_seq_shared_with_protein.csv"
pQTL_rna_seq_df = pd.read_csv(pQTL_rna_seq_path)
pQTL_rna_seq_ids = pQTL_rna_seq_df['Uniprot_Id'].to_numpy()
pQTL_rna_seq_feature_matrix_with_na = pQTL_rna_seq_df.iloc[:,2:].to_numpy()
# pQTL_rna_seq_idx_lookup_path = f"{PATH_ROOT}/data_sources/pQTL/pQTL_protein_converted_idx_lookup.json"
# pQTL_rna_seq_idx_lookup = read_json_from(pQTL_protein_idx_lookup_path)

nikolai_protein_path = f"{PATH_ROOT}/data_sources/Nikolai/Proteins-processed.csv"
nikolai_protein_df = pd.read_csv(nikolai_protein_path)
nikolai_protein_ids = nikolai_protein_df['uniprot_id'].to_numpy()
nikolai_protein_feature_matrix = nikolai_protein_df.iloc[:,1:].fillna(0).to_numpy()
nikolai_protein_idx_lookup_path = f"{PATH_ROOT}/data_sources/Nikolai/protein_processed_lookup.json"
nikolai_protein_idx_lookup = read_json_from(nikolai_protein_idx_lookup_path)