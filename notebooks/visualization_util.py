import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import PrecisionRecallDisplay
from matplotlib import pyplot
from sklearn.metrics import auc

# Visualization

def calc_roc_auc(vector_list,label_vector):
    return list(map(lambda x: roc_auc_score(label_vector,x),vector_list))

def calc_roc_fpr_tpr(vector_list,label_vector):
    tuples_to_return = []
    for vector in vector_list:
        fpr,tpr,_ = roc_curve(label_vector, vector)
        tuples_to_return.append((fpr,tpr))
    return tuples_to_return

def draw_roc_curve(vector_list,label_vector,vector_names,sample_source,validation_source):
    fpr_rpr_list = calc_roc_fpr_tpr(vector_list,label_vector)
    auc_list = calc_roc_auc(vector_list,label_vector)
    for i in range(len(vector_list)):
        pyplot.plot(fpr_rpr_list[i][0], fpr_rpr_list[i][1], linestyle='--', label=f'{vector_names[i]} AUC: {auc_list[i]}')
    # axis labels
    pyplot.xlabel('False Positive Rate')
    pyplot.ylabel('True Positive Rate')
    pyplot.title(f'ROC Curve of {sample_source}, validated against {validation_source}')
    # show the legend
    pyplot.legend()
    # show the plot
    pyplot.show()
    
def df_roc_analysis(df,vector_col_names,label_col_name,df_sample_source,df_validation_source,include_no_skill=True):
    vectors_to_analyze = []
    for name in vector_col_names:
        vectors_to_analyze.append(df[name].to_numpy())
    label_vector = df[label_col_name].to_numpy()
    if include_no_skill:
        vectors_to_analyze.append(np.zeros(len(label_vector)))
        vector_col_names.append("No Skill")
#     normalized_manhattan_vector = 1 - convert_nan_to_one(df['normalized_manhattan_distance'].to_numpy())
    draw_roc_curve(vectors_to_analyze,label_vector,vector_col_names,df_sample_source,df_validation_source)


def df_precision_recall_analysis(ground_truth,thresholds,model_label,title):
    prec, recall, prc_threhsold = precision_recall_curve(ground_truth,thresholds)
    no_skill_prec,no_skill_recall, _ = precision_recall_curve(ground_truth,np.zeros(len(ground_truth)))
    model_auc = auc(recall,prec)
    no_skill_auc = auc(no_skill_recall,no_skill_prec)
    pyplot.plot(recall, prec, marker='.', label=f'{model_label}, AUC: {model_auc}')
    pyplot.plot(no_skill_recall, no_skill_prec, marker='.', label=f'No Skill, AUC: {no_skill_auc}')
    pyplot.xlabel('Recall')
    pyplot.ylabel('Precision')
    pyplot.title(title)
    pyplot.legend()
    pyplot.show()