import numpy as np
import pandas as pd
import os
from util import read_json_from
from functools import reduce
import time


PATH_ROOT = os.getcwd().replace("\\","/").replace("/notebooks","")
core_complex_id_to_subunits = read_json_from(f"{PATH_ROOT}/data_sources/Corum/core_complexes_id_to_subunits.json")
all_complex_id_to_subunits = read_json_from(f"{PATH_ROOT}/data_sources/Corum/all_complexes_id_to_subunits.json")

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

def validate_relations_in_df(df,timed=False):
    start_time = time.time()
    binary_interaction = []
    cocomplex = []
    #subcellular_location = []
    for idx,row in df.iterrows():
        p1 = row['protein_1']
        p2 = row['protein_2']
        relation_report = eval_relation(p1,p2)
        if relation_report['experiments'] > 0 or relation_report['database'] > 0:
            binary_interaction.append(1)
        else:
            binary_interaction.append(0)
        if relation_report['cocomplex'] > 0:
            cocomplex.append(1)
        else:
            cocomplex.append(0)
        if timed and idx % 10000 == 1:
            time_used = time.time() - start_time
            rate = idx / time_used
            eta = (len(df) - idx) / rate
            print(f"Current Index: {idx}, ETA: {eta}")
        # p1_loc = get_protein_location(p1)
        # p2_loc = get_protein_location(p2)
        # if len(np.intersect1d(p1_loc,p2_loc)) > 0:
        #     subcellular_location.append(1)
        # else:
        #     subcellular_location.append(0)
    to_return_df = df.copy()
    to_return_df['binary_interaction'] = binary_interaction
    to_return_df['cocomplex'] = cocomplex
    #to_return_df['subcellular_location'] = subcellular_location
    return to_return_df

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

protein_to_included_core_complex_id_path = f'{PATH_ROOT}/data_sources/Corum/protein_to_included_core_complex_id.json'
protein_to_included_core_complex_id_json = read_json_from(protein_to_included_core_complex_id_path)

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
def is_hit(clique,return_complexes=False,exact_hit=False):
    shared_complexes = []
    for protein in clique:
        shared_complexes.append(protein_to_included_complex_id_json[protein])
    overlapped_complexes = reduce(np.intersect1d,shared_complexes)
    overlapped_complexes = list(map(lambda x: int(x),overlapped_complexes))
    if exact_hit:
        for com in overlapped_complexes:
          try:
            if len(all_complex_id_to_subunits[str(com)]) == len(clique):
                return [com]
          except KeyError:
            pass
        return [] if return_complexes else False
    return overlapped_complexes if return_complexes else len(overlapped_complexes) > 0
# def is_hit(clique,return_complexes=False):
#   shared_complexes = []
#   for protein in clique:
#     shared_complexes.append(protein_to_included_complex_id_json[protein])
#   overlapped_complexes = reduce(np.intersect1d,shared_complexes)
#   overlapped_complexes = list(map(lambda x: int(x),overlapped_complexes))
#   return overlapped_complexes if return_complexes else len(overlapped_complexes) > 0
      
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
  exact_hitting_cliques = []
  corresponding_exact_complexes = []
  hitting_cliques = []
  corresponding_complexes = []
  non_hitting_cliques = []
  for clique in cliques:
    hitting_complex = is_hit(clique,return_complexes=True)
    if (len(hitting_complex) > 0):
      exact_hit = is_hit(clique,return_complexes=True,exact_hit=True)
      if (len(exact_hit) > 0):
        exact_hitting_cliques.append(clique)
        corresponding_exact_complexes.append(exact_hit)
      hitting_cliques.append(clique)
      corresponding_complexes.append(hitting_complex)
    else:
      non_hitting_cliques.append(clique)
  return exact_hitting_cliques,hitting_cliques,non_hitting_cliques,corresponding_exact_complexes,corresponding_complexes

def compute_cliques_info_dict(cliques,recall_denominator="all_verifiable_corum_complexes_size_adjusted"):
  info_dict = {}
  exact_hitting,hitting,non_hitting,corr_exact_complexes,corr_complexes = partition_by_hitting_status(cliques)
  complex_hit,complex_candidates = calc_recall(cliques,return_raw=True)
  info_dict["total_cliques"] = len(cliques)
  info_dict['exact_hitting_cliques'] = exact_hitting
  info_dict['hitting_cliques'] = hitting
  info_dict['non_hitting_cliques'] = non_hitting
  info_dict['exact_complexes'] = list(map(list,corr_exact_complexes))
  info_dict['corresponding_complexes'] = list(map(list,corr_complexes))
  info_dict['exact_precision'] = len(exact_hitting) / len(cliques)
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
    columns = ['clique_size','exact_hitting_cliques','hitting_cliques','total_cliques','exact_precision','precision','complex_hit']
    row_values = []
    for size in report_dict:
        row_values.append((size,
        len(report_dict[size]['exact_hitting_cliques']),
        len(report_dict[size]['hitting_cliques']),
        report_dict[size]['total_cliques'],
        report_dict[size]['exact_precision'],
        report_dict[size]['precision'],
        len(report_dict[size]['complex_hit'])))
    df = pd.DataFrame(data=row_values,columns=columns)
    return df

# validation
# def calc_interaction_count(protein,candidates,interaction='physical'):
#     count = 0
#     for candidate in candidates:
#         relation_report = eval_relation(protein,candidate)
#         if interaction == 'physical':
#             if relation_report['experiments'] > 0 or relation_report['database'] > 0:
#                 count += 1
#         elif interaction == 'cocomplex':
#             if relation_report['cocomplex'] > 0:
#                 count += 1
#     return count

# def overall_interaction_counts(proteins,interation_candidates,interaction='physical'):
#     assert len(proteins) == len(interation_candidates)
#     count_list = []
#     for i in range(len(proteins)):
#         short_list = []
#         # Number of interacting proteins
#         short_list.append(calc_interaction_count(proteins[i],interation_candidates[i],interaction))
#         # Number of overlapping proteins
#         short_list.append(len(interation_candidates[i]))
#         count_list.append(short_list)
#     return count_list