## load libriaries

import os
import pandas as pd
import numpy as np
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
import scipy 
import scipy.stats as ss
import statsmodels
from statsmodels import stats
from statsmodels.stats import multitest
from gsva_py import gsva_py
#### User input: User can select functional moduels and input the gene expression matrix 

ROOT_DIR = os.path.abspath("../")

#### Load functional modules
def load_function_modules(FM):
    Functional_modules = pd.read_csv(FM,header= None, names = ("name","member"))
    Functional_modules['description'] = Functional_modules['name']
    Functional_modules.columns = ["name",'member','description']
    dic_module = {}
    module_list = list(set(Functional_modules['name']))#
    for module in module_list:
        dic_module[module] = list(set(Functional_modules[Functional_modules['name'] == module]['member']))
    return(dic_module, Functional_modules)


#### Define functional module BIQs. Here we defined gene_level_Up, gene_level_Down, pathway_level_ssGSEA, and TFs_network_level_TF as shown below. User can also define their own BIQs. 

def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3 

#Define BIQ 1: gene_level_Up
def ratio_up(x):
    result = len(x[x > 1.96])/len(x)  #1.6
    return(result)

def gene_level_Up(data_matrix, FM_50,Functional_modules):
    matrix_gene_level_Up = pd.DataFrame()
    for i in range(0,len(FM_50)):
        matrix_module = data_matrix[intersection(data_matrix.columns, list(Functional_modules[Functional_modules['name'] == FM_50[i]]['member']))]
        
        module_up_result = []
        #for sample in range(0,matrix_module.shape[0]):
        for sample in matrix_module.index.values:
            x = matrix_module.loc[sample]
            module_up_result.append(ratio_up(x))
        matrix_gene_level_Up[FM_50[i]] = module_up_result
    #matrix_gene_level_Up.index.names  = data_matrix.index
    
    matrix_gene_level_Up.columns = matrix_gene_level_Up.columns.values + '_Up'
    matrix_gene_level_Up.index = list(data_matrix.index.values)
    return(matrix_gene_level_Up)

#Define BIQ 2: gene_level_Down

def ratio_down(x):
    result = -1 * len(x[x < -1.96])/len(x)  #-1.6
    return(result)

def gene_level_Down(data_matrix, FM_50,Functional_modules):
    matrix_gene_level_Down = pd.DataFrame()
    for i in range(0,len(FM_50)):
        matrix_module = data_matrix[intersection(data_matrix.columns,list(Functional_modules[Functional_modules['name'] == FM_50[i]]['member']))]
        #matrix_module.index = data_matrix['Unnamed: 0']
        
        module_down_result = []
        #for sample in range(0,matrix_module.shape[0]):
        for sample in matrix_module.index.values:
            x = matrix_module.loc[sample]
            module_down_result.append(ratio_down(x))
        
        matrix_gene_level_Down[FM_50[i]] = module_down_result

    matrix_gene_level_Down.columns = matrix_gene_level_Down.columns.values+'_Down'
    matrix_gene_level_Down.index = list(data_matrix.index.values)
    #print(matrix_gene_level_Down.index.values)
    return(matrix_gene_level_Down)

#Define BIQ 3: pathway_level_ssGSEA

def pathway_level_ssGSEA(data_matrix, FM_50,Functional_modules):
    
    tempdir = ROOT_DIR+"/temp"

    if os.path.exists(tempdir) == False:
        try:
            os.makedirs(tempdir)
        except OSError:
            print ("Creation of the directory %s failed" % tempdir)
        else:
            print ("Successfully created the directory %s " % tempdir)
    else:
        print ("INfO:  %s already exists!" % tempdir)

    #from GSVA import gsva, gmt_to_dataframe
    import math
    limit = 1000
    n = math.ceil(data_matrix.shape[0] / limit)
    
    module_selected_gmt = Functional_modules.loc[Functional_modules['name'].isin( FM_50) ]
    matrix_BIQ_gsva = pd.DataFrame()
    result_list = []

    for i in range(0,n):
        print(i)
        data = data_matrix.iloc[i*1000: i*1000+1000]
        pathways_df = gsva_py(data.T, geneset_df=module_selected_gmt,
                        method='ssgsea',
                        kcdf='Gaussian',
                        abs_ranking=False,
                        min_sz=1,
                        max_sz=None,
                        parallel_sz=4,
                        mx_diff=True,
                        ssgsea_norm=True,
                        verbose=False,
                        tempdir= tempdir)

        result = pathways_df.T
        result.columns = result.columns.values + '_ssGSEA'
        result.index = list(data.index.values)
        result_list.append(result)
    matrix_BIQ_gsva = pd.concat(result_list, axis = 0)
    return(matrix_BIQ_gsva)

#Define BIQ 4: TFs_network_level_TF

def TFs_network_level_TF(data_matrix,FM_50,TFs): 
    matrix_TFs_network_level_TF = pd.DataFrame()
    for module in range(0,len(FM_50)):
        expr_tf = data_matrix[intersection(list(set(TFs.loc[TFs['Pathway'] == FM_50[module]]['TFs'])),data_matrix.columns.values )]
        tf_module = intersection(list(set(TFs.loc[TFs['Pathway'] == FM_50[module]]['TFs'])),expr_tf.columns.values )
        TFs_g = TFs.loc[TFs['TFs'].isin(tf_module)]
        TFs_g_m = TFs_g.loc[TFs_g['Pathway'] == FM_50[module]]
        expr_tf_order = expr_tf[list(TFs_g_m['TFs'])]
        
        f = TFs_g_m['Score']/sum(TFs_g_m['Score'])
        
        tf_list = list()

        for sample in range(0,expr_tf_order.shape[0]):
            sum_tf = 0
            for gene in range(0,len(expr_tf_order.iloc[1])):
                sum_tf = sum_tf + expr_tf_order.iloc[sample][gene] * f.iloc[gene]
            tf_list.append(sum_tf)
        
        matrix_TFs_network_level_TF[FM_50[module]] = tf_list

    matrix_TFs_network_level_TF.columns = matrix_TFs_network_level_TF.columns.values + '_TF'
    
    matrix_TFs_network_level_TF.index = data_matrix.index
    return(matrix_TFs_network_level_TF)


def generate_BIQ(data_matrix, FM_50, dic_module, TFs, UP = False, DOWN = False, ssGSEA = False, TF = False): 
    
    BIQ_matrix_list = list()
    matrix_gene_level_Up = pd.DataFrame()
    matrix_gene_level_Down = pd.DataFrame()
    matrix_pathway_level_ssGSEA = pd.DataFrame()
    matrix_TFs_network_level_TF = pd.DataFrame()
    if UP == True:
        matrix_gene_level_Up = gene_level_Up(data_matrix, FM_50, dic_module)
        BIQ_matrix_list.append(matrix_gene_level_Up)
        
    if DOWN == True:
        matrix_gene_level_Down = gene_level_Down(data_matrix, FM_50, dic_module)
        BIQ_matrix_list.append(matrix_gene_level_Down)
        
    if ssGSEA == True:
        matrix_pathway_level_ssGSEA = pathway_level_ssGSEA(data_matrix,  FM_50, dic_module)
        BIQ_matrix_list.append(matrix_pathway_level_ssGSEA)
        
    if TF == True:
        matrix_TFs_network_level_TF = TFs_network_level_TF(data_matrix, FM_50, TFs)
        BIQ_matrix_list.append(matrix_TFs_network_level_TF)
        

    if len(BIQ_matrix_list) == 0:
        print("Error!! No BIQ is selected!")

    matrix_BIQ = pd.concat(BIQ_matrix_list, axis = 1)
    return(matrix_BIQ)
