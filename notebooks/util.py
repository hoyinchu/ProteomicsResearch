import numpy as np
import pandas as pd
from functools import reduce
import json
from csv import writer
import requests, sys
import time
import os

PATH_ROOT = os.getcwd().replace("\\","/").replace("/notebooks","")

# Given two vectors, return a subsets where neither are 0s
def find_common_observations(vector1,vector2):
  vector1_bool = np.where(vector1 != 0, 1, 0)
  vector2_bool = np.where(vector2 != 0, 1, 0)
  take_indices = np.logical_and(vector1_bool,vector2_bool)
  take_indices = take_indices.nonzero()[0]
  x1 = np.take(vector1,take_indices)
  x2 = np.take(vector2,take_indices)
  return x1,x2

# Appends list of tuples to csv
def append_to_csv(file_name,column_names,cor_tuples):
  with open(file_name, 'a+', newline='') as write_obj:
    csv_writer = writer(write_obj)
    if column_names is not None:
      csv_writer.writerow(column_names)
    for cor_tuple in cor_tuples:
      csv_writer.writerow(cor_tuple)

# Writes the given dict as json to path
def write_json_to(json_dict,path):
  json_to_write = json.dumps(json_dict)
  write_file = open(path,"w")
  write_file.write(json_to_write)
  write_file.close()

# Returns a json object from given path
def read_json_from(path):
  with open(path, "r") as read_file:
    init_dict = json.load(read_file)
    if isinstance(init_dict,str):
        return eval(init_dict)
    return init_dict

# Returns a json that contains all uniprot ids in the given protein vector
def request_protein_info(protein_vec):
    json_to_build = {}
    base_url = "https://www.ebi.ac.uk/proteins/api/proteins"
    # 100 requests per second
    for protein in protein_vec:
        time.sleep(0.01)
        request_url = f"{base_url}/{protein}"
        r = requests.get(request_url, headers={ "Accept" : "application/json"})
        if not r.ok:
          r.raise_for_status()
          sys.exit()
        responseBody = json.loads(r.text)
        json_to_build[protein] = responseBody
    return json_to_build

coexpression_lookup_path = f"{PATH_ROOT}/data_sources/StringDB/human/medium_confidence_coexpression_relation_lookup.json"
cooccurence_lookup_path = f"{PATH_ROOT}/data_sources/StringDB/human/medium_confidence_cooccurence_relation_lookup.json"
experiments_lookup_path = f"{PATH_ROOT}/data_sources/StringDB/human/medium_confidence_experiments_relation_lookup.json"
fusion_lookup_path = f"{PATH_ROOT}/data_sources/StringDB/human/medium_confidence_fusion_relation_lookup.json"
homology_lookup_path = f"{PATH_ROOT}/data_sources/StringDB/human/medium_confidence_homology_relation_lookup.json"
cocomplex_lookup_path = f"{PATH_ROOT}/data_sources/Corum/all_corum_complex_pairs_size_only.json"
database_lookup_path = f"{PATH_ROOT}/data_sources/StringDB/human/medium_confidence_database_relation_lookup.json"

coexpression_lookup = read_json_from(coexpression_lookup_path)
cooccurence_lookup = read_json_from(cooccurence_lookup_path)
experiments_lookup = read_json_from(experiments_lookup_path)
fusion_lookup = read_json_from(fusion_lookup_path)
homology_lookup = read_json_from(homology_lookup_path)
cocomplex_lookup = read_json_from(cocomplex_lookup_path)
database_lookup = read_json_from(database_lookup_path)

lookup_list = [coexpression_lookup,cooccurence_lookup,experiments_lookup,fusion_lookup,homology_lookup,cocomplex_lookup,database_lookup]


# returns report json containing relations between two proteins
def eval_relation(protein1,protein2):
    report_json = {}
    for lookup in lookup_list:
        lookup_type = lookup['relation_type']
        relation_score = float('NaN')
        try:
            relation_score = lookup[protein1][protein2]
        except KeyError:
            pass
        report_json[lookup_type] = relation_score
    return report_json

# def generate_link_candidates(protein,name_vector,embedding,lookup,candidate_size=5,neighbor_size=30,chunk_size=5,dist=euclidean_dist,return_all=False):
#     protein_vec = embedding[lookup[protein]]
#     diff_embeddeding = embedding - protein_vec
#     chunked_embeddings = chunk_features(diff_embeddeding,chunk_size)
#     candidates = np.zeros((len(chunked_embeddings),neighbor_size))
#     for idx,chunk in enumerate(chunked_embeddings):
#         norm_vec = np.linalg.norm(chunk, axis=1)
#         top_indices = take_top_k(norm_vec,neighbor_size+1)[1:]
#         candidates[idx] = top_indices
#     interaction_freq_indices = calc_appearance_freq(candidates)
#     top_interacting_indices = take_top_k(interaction_freq_indices[:,1],candidate_size,desc=True)
#     interacting_counts = interaction_freq_indices[top_interacting_indices][:,1]
#     top_interacting = name_vector[interaction_freq_indices[top_interacting_indices][:,0].astype(int)]
#     return top_interacting,interacting_counts

corum_complexes_path = f"{PATH_ROOT}/data_sources/Corum/allComplexes.txt"
all_corum_complex_pairs_json_path = f"{PATH_ROOT}/data_sources/Corum/all_corum_complex_pairs.json"

corum_complexes_dataframe = pd.read_csv(corum_complexes_path,sep='\t')
all_corum_subunits = corum_complexes_dataframe["subunits(UniProt IDs)"].to_numpy()
all_corum_subunits_list = list(map(lambda x: x.split(";"), all_corum_subunits))
all_corum_complexes_id = corum_complexes_dataframe["ComplexID"].to_numpy()

corum_complex_pairs_lookup_dict = read_json_from(all_corum_complex_pairs_json_path)
all_unique_corum_proteins = list(corum_complex_pairs_lookup_dict.keys())

protein_to_included_complex_id_path = f"{PATH_ROOT}/data_sources/Corum/protein_to_included_complex_id.json"
protein_to_included_complex_id_json = read_json_from(protein_to_included_complex_id_path)

proteomeHD_uniprot_to_idx_lookup_json_path = f"{PATH_ROOT}/data_sources/ProteomeHD/uniprot_id_idx_lookup.json"
proteomeHD_uniprot_to_idx_lookup_json = read_json_from(proteomeHD_uniprot_to_idx_lookup_json_path)

proteomeHD_verifiable_corum_complexes_amount_json_path = f"{PATH_ROOT}/data_sources/Corum/proteomeHD_verifiable_complexes_amount.json"
proteomeHD_verifiable_corum_complexes_amount_json = read_json_from(proteomeHD_verifiable_corum_complexes_amount_json_path)

def is_verifiable(to_verify,verify_list,min_verify_count):
  count = 0
  for node in to_verify:
    if node in verify_list:
      count += 1
  return count >= min_verify_count

def is_corum_verifiable(to_verify):
    if len(to_verify) < 2: return False
    else: return is_verifiable(to_verify,corum_complex_pairs_lookup_dict.keys(),len(to_verify))

def is_proteomeHD_verifiable(subunits,min_count=2):
  return is_verifiable(subunits.split(";"),proteomeHD_uniprot_to_idx_lookup_json.keys(),min_count)

all_proteomeHD_verifiable_corum_subunits = list(filter(lambda x: is_proteomeHD_verifiable(x), all_corum_subunits))


# Checks if a given clique (tuple or list) is a hit
# When the proteins in the clique is a subset of a protein subunits in a protein
# complex in the Corum dataset
def is_hit(clique,return_complexes=False):
  shared_complexes = []
  for protein in clique:
    shared_complexes.append(protein_to_included_complex_id_json[protein])
  overlapped_complexes = reduce(np.intersect1d,shared_complexes)
  overlapped_complexes = list(map(lambda x: int(x),overlapped_complexes))
  return overlapped_complexes if return_complexes else len(overlapped_complexes) > 0
      
# Calculates the precision of the given list of cliques
# precison = number of hitting cliques / total number of cliques
def calc_precision(cliques,loud=False):
  hit_num = 0
  for clique in cliques:
    if(is_hit(clique)):
      hit_num += 1
  if loud:
    print(f"Total Number of Cliques: {len(cliques)}")
    print(f"Hitting Cliques: {hit_num}")
    print(f"Precision: {hit_num / len(cliques)}")
  return hit_num / len(cliques)

# Calculates the recall of the given cliques
# recall = Number of complexes hit by cliques / total number of candidate complexes
# candidate complexes are complexes that contain at least 1 protein listed in cliques
def calc_recall(cliques,loud=False,return_raw=False,recall_denominator="all_corum_complexes"):
  denominator = len(all_proteomeHD_verifiable_corum_subunits)
  flattened_protein_list = [item for sublist in cliques for item in sublist]
  # All the unique proteins in the cliques
  all_clique_proteins = np.unique(flattened_protein_list)

  # Find all complexes such that it contains at least of the proteins listed in cliques
  all_candidate_complexes_id = []
  for idx in range(len(corum_complexes_dataframe)):
    if len(np.intersect1d(all_clique_proteins,all_corum_subunits_list[idx])) > 0:
      all_candidate_complexes_id.append(all_corum_complexes_id[idx])
  if (recall_denominator == "complex_candidate"):
    denominator = len(all_candidate_complexes_id)
  elif (recall_denominator == "all_corum_complexes"):
    denominator = len(all_corum_subunits)
  # Find all complexes hit by cliques
  all_complexes_hit = []
  for clique in cliques:
    all_complexes_hit.append(is_hit(clique,return_complexes=True))
  flat_list = [item for sublist in all_complexes_hit for item in sublist]
  all_complexes_hit_unique = np.unique(flat_list)
  if loud:
    print(f"Total Number of Complex Hit: {len(all_complexes_hit_unique)}")
    print(f"Total Number of Complex Candidates: {len(all_candidate_complexes_id)}")
    print(f"Recall: {len(all_complexes_hit_unique) / len(all_candidate_complexes_id)}")
  if return_raw:
    all_complexes_hit_unique = list(map(lambda x: int(x),all_complexes_hit_unique))
    all_candidate_complexes_id = list(map(lambda x: int(x),all_candidate_complexes_id))
    return all_complexes_hit_unique,all_candidate_complexes_id
  return len(all_complexes_hit_unique) / denominator

# Give cliques, return all cliques that are hitting, non hitting, and the complexes
# they hit for all hitting cliques
def partition_by_hitting_status(cliques):
  hitting_cliques = []
  corresponding_complexes = []
  non_hitting_cliques = []
  for clique in cliques:
    hitting_complex = is_hit(clique,return_complexes=True)
    if (len(hitting_complex) > 0):
      hitting_cliques.append(clique)
      corresponding_complexes.append(hitting_complex)
    else:
      non_hitting_cliques.append(clique)
  return hitting_cliques,non_hitting_cliques,corresponding_complexes

def compute_cliques_info_dict(cliques,recall_denominator="all_verifiable_corum_complexes_size_adjusted"):
  info_dict = {}
  hitting,non_hitting,corr_complexes = partition_by_hitting_status(cliques)
  complex_hit,complex_candidates = calc_recall(cliques,return_raw=True)
  info_dict["total_cliques"] = len(cliques)
  info_dict['hitting_cliques'] = hitting
  info_dict['non_hitting_cliques'] = non_hitting
  info_dict['corresponding_complexes'] = list(map(list,corr_complexes))
  info_dict['precision'] = len(hitting) / len(cliques)
  info_dict['complex_hit'] = list(complex_hit)
  info_dict['complex_candidates_id'] = complex_candidates
  info_dict['recall'] = {}
  info_dict['recall']["complex_candidate"] = len(complex_hit) / len(complex_candidates)
  info_dict['recall']["all_corum_complexes"] = len(complex_hit) / len(all_corum_subunits)
  info_dict['recall']["all_verifiable_corum_complexes"] = len(complex_hit) / len(all_proteomeHD_verifiable_corum_subunits)
  try:
    denom = proteomeHD_verifiable_corum_complexes_amount_json[str(len(cliques[0]))]
    if denom > 0:
        info_dict['recall']["all_verifiable_corum_complexes_size_adjusted"] = len(complex_hit) / denom
    else:
        info_dict['recall']["all_verifiable_corum_complexes_size_adjusted"] = float("NaN")
  except KeyError:
    info_dict['recall']["all_verifiable_corum_complexes_size_adjusted"] = float("NaN")
  return info_dict

def full_report_calculation_by_size(cliques):
  clique_partition = partition_cliques_by_size(cliques)
  report_dict = {}
  for size in clique_partition:
    info_dict = compute_cliques_info_dict(clique_partition[size])
    report_dict[size] = info_dict
  return report_dict

def partition_cliques_by_size(cliques):
  counter = {}
  for clique in cliques:
    try:
      counter[len(clique)].append(clique)
    except KeyError:
      counter[len(clique)] = [clique]
  return counter

def simplifed_report_df(report_dict):
    columns = ['clique_size','hitting_cliques','total_cliques','precision','complex_hit']
    row_values = []
    for size in report_dict:
        row_values.append((size,len(report_dict[size]['hitting_cliques']),report_dict[size]['total_cliques'],report_dict[size]['precision'],len(report_dict[size]['complex_hit'])))
    df = pd.DataFrame(data=row_values,columns=columns)
    return df


# def calc_combined_counts(proteins,candidates_1,candidates_2):
#     assert len(proteins) == len(candidates_1) and len(candidates_1) == len(candidates_2)
#     count_list = []
#     for i in range(len(proteins)):
#         short_list = []
#         intersection = np.intersect1d(candidates_1[i],candidates_2[i])
#         jaccard = len(intersection) / (len(candidates_1[i]) + len(candidates_2[i]) - len(intersection))
#         # Number of interacting proteins
#         short_list.append(calc_interaction_count(proteins[i],intersection))
#         # Number of overlapping proteins
#         short_list.append(len(intersection))
#         # Jaccard Similarity
#         short_list.append(jaccard)
#         count_list.append(short_list)
#     return count_list

# def combined_count_statistics(combined_count):
#     combined_count = np.array(combined_count)
#     # ALL
#     combined_filtered_all = np.array(list(filter(lambda x: x[1] >= 1, combined_count)))
#     combined_filtered_all_precision = np.sum(combined_filtered_all[:,0]) / np.sum(combined_filtered_all[:,1])
#     # J >= 0.2
#     combined_filtered_J = np.array(list(filter(lambda x: x[2] >= 0.2, combined_count)))
#     combined_filtered_J_precision = 0
#     if (len(combined_filtered_J) > 0):
#         combined_filtered_J_precision = np.sum(combined_filtered_J[:,0]) / np.sum(combined_filtered_J[:,1])
#     # Known >= 1
#     combined_filtered_known = np.array(list(filter(lambda x: x[0] >= 1 and x[1] >= 2, combined_count)))
#     combined_filtered_known_precision = 0
#     if (len(combined_filtered_known) > 0):
#         combined_filtered_known_precision = np.sum(combined_filtered_known[:,0]) / np.sum(combined_filtered_known[:,1])
#     return combined_filtered_all_precision,combined_filtered_J_precision,combined_filtered_known_precision,np.sum(combined_filtered_all[:,1]),np.sum(combined_filtered_J[:,1]),np.sum(combined_filtered_known[:,1])