import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import PrecisionRecallDisplay
import matplotlib.pyplot as plt
from sklearn.metrics import auc
import copy

# Visualization

def calc_roc_auc(vector_list,label_vector):
    return list(map(lambda x: roc_auc_score(label_vector,x),vector_list))

def calc_prc_auc(vector_list,label_vector):
    prc_tuple_list = []
    for vector in vector_list:
        prec, recall, prc_threshold = precision_recall_curve(label_vector,vector)
        prc_tuple_list.append((prec,recall,prc_threshold))
    return prc_tuple_list

def calc_roc_fpr_tpr(vector_list,label_vector):
    tuples_to_return = []
    for vector in vector_list:
        fpr,tpr,_ = roc_curve(label_vector, vector)
        tuples_to_return.append((fpr,tpr))
    return tuples_to_return

def draw_roc_curve(vector_list,label_vector,vector_names,title,auc_only=False):
    fpr_rpr_list = calc_roc_fpr_tpr(vector_list,label_vector)
    auc_list = calc_roc_auc(vector_list,label_vector)
    if auc_only:
        return auc_list
    for i in range(len(vector_list)):
        plt.plot(fpr_rpr_list[i][0], fpr_rpr_list[i][1], linestyle='--', label=f'{vector_names[i]} AUC: {int(10000*auc_list[i])/10000}')
    # axis labels
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    # show the legend
    plt.legend()
    # show the plot
    plt.show()
    return auc_list
    
def df_roc_analysis(df,vector_col_names,label_col_name,title,auc_only=False):
    vectors_to_analyze = []
    for name in vector_col_names:
        vectors_to_analyze.append(df[name].to_numpy())
    #print(vectors_to_analyze)
    vectors_to_analyze.append(np.zeros(len(df)))
    #print(vectors_to_analyze)
    label_vector = df[label_col_name].to_numpy()
    vector_col_names_with_no_skill = copy.copy(vector_col_names)
    vector_col_names_with_no_skill.append("No Skill")
    return draw_roc_curve(vectors_to_analyze,label_vector,vector_col_names_with_no_skill,title,auc_only)


def df_precision_recall_analysis(df,vector_col_names,label_col_name,title,ylim=None,auc_only=False,legend_placement='best'):
    if ylim:
        plt.ylim(ylim)
    ground_truth = df[label_col_name].to_numpy()
    threshold_vectors = [df[name].to_numpy() for name in vector_col_names]
    tuple_list = calc_prc_auc(threshold_vectors,ground_truth)
    model_aucs = []
    for tup in tuple_list:
        model_aucs.append(auc(tup[1],tup[0]))
    vector_col_names_with_baseline = copy.copy(vector_col_names)
    vector_col_names_with_baseline.append('Baseline')
    no_skill_prec,no_skill_recall, no_skill_thresholds = precision_recall_curve(ground_truth,np.zeros(len(ground_truth)))
    no_skill_prec = [np.sum(ground_truth) / len(ground_truth) for i in range(len(no_skill_recall))]
    no_skill_auc = auc(no_skill_recall,no_skill_prec)
    tuple_list.append((no_skill_prec,no_skill_recall,no_skill_thresholds))
    model_aucs.append(no_skill_auc)
    if auc_only:
        return model_aucs
    for i in range(len(vector_col_names_with_baseline)):
        plt.plot(tuple_list[i][1], tuple_list[i][0], linestyle='--', label=f'{vector_col_names_with_baseline[i]} AUC: {int(10000*model_aucs[i])/10000}')
    # plt.plot(recall, prec, marker='.', label=f'{model_label}, AUC: {model_auc}')
    # plt.plot(recall, no_skill_prec, marker='.', label=f'Baseline, AUC: {no_skill_auc}')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(title)
    plt.legend(loc=legend_placement)
    plt.show()
    return model_aucs