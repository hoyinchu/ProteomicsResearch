import numpy as np
from numpy import nan
import pandas as pd

# kNN
# def get_top_k_nearest_neighbors(vector,k,candidates,dist_function):
#     k += 1
#     dist_to_all = np.array(list(map(lambda x: dist_function(vector,x),candidates)))
#     neighbor_indices = np.argpartition(dist_to_all, k)[0:k]
#     neighbor_indices = neighbor_indices[np.argsort(dist_to_all[neighbor_indices])]
#     neighbor_indices = neighbor_indices[1:]
#     top_k_neighbors_dist = dist_to_all[neighbor_indices]
#     return neighbor_indices,top_k_neighbors_dist

# def get_closest_neighbors_in_embeddings(protein,k,embeddings,lookup,dist_func):
#     neighbor_lists = []
#     for embedding in embeddings:
#         protein_vec = embedding[lookup[protein]]
#         neighbor_indices,_ = get_top_k_nearest_neighbors(protein_vec,k,embedding,dist_func)
#         neighbor_lists.append(neighbor_indices)
#     return neighbor_lists

# def get_common_closest_neighbors_in_embeddings(protein,k,embeddings,lookup,dist_func):
#     neighbor_lists = get_closest_neighbors_in_embeddings(protein,k,embeddings,lookup,dist_func)
#     shared_proteins_indices = reduce(np.intersect1d, neighbor_lists)
#     name_list = list(map(lambda x: shared_proteins[x], shared_proteins_indices))
#     return name_list


def convert_prediction_to_df(proteins,interation_candidates,count_matrix,feature_col='appearance'):
    assert len(proteins) == len(interation_candidates) and len(proteins) == len(count_matrix)
    predictions = []
    for i in range(len(proteins)):
        for j in range(len(interation_candidates[i])):
            predictions.append((proteins[i],interation_candidates[i][j],count_matrix[i][j]))
    df = pd.DataFrame(predictions, columns =['protein_1', 'protein_2', feature_col]) 
    return df

def euclidean_dist(vec1,vec2):
    return np.linalg.norm(vec1-vec2)
    
# Cut an embedding into chunks with size n
def chunk_features(embedding,n):
    columns_to_cut = [ i for i in range(n,embedding.shape[1],n)]
    return np.hsplit(embedding,columns_to_cut)

def calc_appearance_freq(neighbors):
    unique, counts = np.unique(neighbors, return_counts=True)
    return np.asarray((unique, counts)).T

# Indices of the top k element in the given array
def take_top_k(arr,k,desc=False):
    if desc: arr = -arr
    if len(arr) <= k:
        return np.argsort(arr)
    part = np.argpartition(arr,k)[:k]
    correct_order = np.argsort(arr[part])
    return part[correct_order]

def sample_embedding(embedding,chunk_size,sample_size):
    return np.array([embedding[:,np.random.choice(embedding.shape[1], chunk_size, replace=False)] for i in range(sample_size)])

# Keep chunk_size small but neighbor size big
def generate_link_candidates(protein,name_vector,embedding,lookup,neighbor_size=5,chunk_size=5,sample_size=100,candidate_size=5,dist=euclidean_dist,return_all=False):
    protein_vec = embedding[lookup[protein]]
    if dist.__name__ == 'euclidean_dist':
        diff_embeddeding = embedding - protein_vec
        chunked_embeddings = sample_embedding(diff_embeddeding,chunk_size,sample_size)
    else:
        chunked_embeddings = sample_embedding(embedding,chunk_size,sample_size)
    #chunked_embeddings = chunk_features(diff_embeddeding,chunk_size)
    candidates = np.zeros((len(chunked_embeddings),neighbor_size))
    for idx,chunk in enumerate(chunked_embeddings):
        if dist.__name__ == 'euclidean_dist':
            norm_vec = np.linalg.norm(chunk, axis=1)
            top_indices = take_top_k(norm_vec,neighbor_size+1)[1:]
        elif dist.__name__ == 'pearson_corr':
            r_vec = np.corrcoef(chunk)[lookup[protein]]
            top_indices = take_top_k(r_vec,neighbor_size+1,desc=True)[1:]
        candidates[idx] = top_indices
    interaction_freq_indices = calc_appearance_freq(candidates)
    if return_all:
        top_interacting = name_vector[interaction_freq_indices[:,0].astype(int)]
        interacting_counts = interaction_freq_indices[:,1]
    else:
        top_interacting_indices = take_top_k(interaction_freq_indices[:,1],candidate_size,desc=True)
        interacting_counts = interaction_freq_indices[top_interacting_indices][:,1]
        top_interacting = name_vector[interaction_freq_indices[top_interacting_indices][:,0].astype(int)]
    return top_interacting,interacting_counts


def generate_knn_interaction_score(protein1,protein2,name_vector,embedding,lookup,neighbor_size=5,chunk_size=5,sample_size=100,candidate_size=5,dist=euclidean_dist,return_all=False):
    protein1_neighbors,protein1_neighbor_counts = generate_link_candidates(protein1,name_vector,embedding,lookup,neighbor_size,chunk_size,sample_size,candidate_size,dist,return_all)
    protein2_neighbors,protein2_neighbor_counts = generate_link_candidates(protein2,name_vector,embedding,lookup,neighbor_size,chunk_size,sample_size,candidate_size,dist,return_all)
    score = 0
    # print(protein1_neighbors)
    # print(protein2_neighbors)
    if protein2 in protein1_neighbors:
        score += protein1_neighbor_counts[np.where(protein1_neighbors == protein2)][0]
    if protein1 in protein2_neighbors:
        score += protein2_neighbor_counts[np.where(protein2_neighbors == protein1)][0]
    return score