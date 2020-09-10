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
def request_protein_info(protein_vec,timed=True):
    if timed:
      starttime = time.time()
    json_to_build = {}
    base_url = "https://www.ebi.ac.uk/proteins/api/proteins"
    # 100 requests per second
    for idx,protein in enumerate(protein_vec):
        time.sleep(0.01)
        request_url = f"{base_url}/{protein}"
        r = requests.get(request_url, headers={ "Accept" : "application/json"})
        if not r.ok:
          print(f"{protein} not found, skipped")
          #r.raise_for_status()
          continue
        responseBody = json.loads(r.text)
        json_to_build[protein] = responseBody
        if idx%100==1 and timed:
          calc_eta(starttime,idx,len(protein_vec))
    return json_to_build

def df_to_json(df,col1='protein_1',col2='protein_2',val='appearance'):
    to_return = {}
    for idx,row in df.iterrows():
        p1 = row[col1]
        p2 = row[col2]
        if not p1 in to_return:
            to_return[p1] = {}
        if not p2 in to_return:
            to_return[p2] = {}
        to_return[p1][p2] = row[val]
        to_return[p2][p1] = row[val]
    return to_return

def generate_pair_distances(p1_col,p2_col,dist_function,feature_matrix,lookup_json):
    assert(len(p1_col) == len(p2_col))
    distances = []
    for i in range(len(p1_col)):
      try:
        protein1_vec = feature_matrix[lookup_json[p1_col[i]]]
        protein2_vec = feature_matrix[lookup_json[p2_col[i]]]
        dist = dist_function(protein1_vec,protein2_vec)
        distances.append(dist)
      except KeyError:
        distances.append(float('NaN'))
    return distances

def construct_complex_df(complex_arr):
    col_names = ['ComplexID','ComplexName','subunit_1','subunit_2']
    com_data = []
    for com in complex_arr:
        subunits = com['subunits(UniProt IDs)'].split(";")
        cur_tuple = [com['ComplexID'],com['ComplexName'],subunits[0],subunits[1]]
        com_data.append(cur_tuple)
    to_return_df = pd.DataFrame(columns=col_names,data=com_data)
    return to_return_df

def get_new_sampled_dataset(pos_df,pos_amount,neg_df,neg_amount):
    pos_rand = np.random.choice(len(pos_df), pos_amount, replace=False)
    neg_rand = np.random.choice(len(neg_df), neg_amount, replace=False)
    pos_df = pos_df.iloc[pos_rand,:]
    neg_df = neg_df.iloc[neg_rand,:]
    new_df = pd.concat([pos_df,neg_df])
    return new_df

def calc_eta(start_time,cur_idx,total_idx):
  rate = cur_idx / (time.time() - start_time)
  eta = (total_idx-cur_idx) / rate
  print(f"ETA: {eta} seconds")

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