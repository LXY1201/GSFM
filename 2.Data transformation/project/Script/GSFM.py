## load libriaries

import sys, os
import time
import subprocess
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
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE
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

#gsva_py

def gsva_py(expression_df,geneset_df=None,
         method='ssgsea',
         kcdf='Gaussian',
         abs_ranking=False,
         min_sz=1,
         max_sz=None,
         parallel_sz=4,
         #parallel_type="SOCK",
         mx_diff=True,
         tau=None,
         ssgsea_norm=True,
         verbose=False,
         tempdir= None
         ):
    df = expression_df
    gmt_df = geneset_df

    if not tempdir:
        tempdir =  mkdtemp(prefix="weirathe.",dir=gettempdir().rstrip('/'))
    if verbose:
        sys.stderr.write("Caching to "+tempdir+"\n")
    # Remove genes from the genesets that do not occur in the dataset
    members = gmt_df['member'].unique()
    missing = set(members)-set(df.index)
    original = df.index
    if len(missing) > 0:
        if verbose: sys.stderr.write("WARNING removing "+str(len(missing))+\
          " genes from gene sets that don't exist in the data\n"+\
          ",".join(sorted(list(missing)))+"\n")
    gmt_df = gmt_df[~gmt_df['member'].isin(list(missing))]
    # Write our gene sets
    gmt_df = gmt_df.groupby(['name']).\
        apply(lambda x: "\t".join(sorted(list(x['member'])))).reset_index().rename(columns={0:'members'})
    of = open(os.path.join(tempdir,"gs.gmt"),'w')
    for row in gmt_df.itertuples():
        name = row.name
        description = 'description'
        fields = row.members
        of.write(name+"\t"+description+"\t"+fields+"\n")
    of.close()
    df.to_csv(os.path.join(tempdir,"expr.csv"))
    cur = os.path.dirname(os.path.realpath(__file__))
    rscript = os.path.join(cur,"gsva.r")
    cmd = ["Rscript",rscript]+[str(x) for x in [method,kcdf,abs_ranking,min_sz,max_sz,parallel_sz,mx_diff,ssgsea_norm,verbose,tempdir]]
    if verbose: sys.stderr.write(" ".join(cmd)+"\n")
    destination = PIPE
    if verbose: destination = sys.stderr
    sp = Popen(cmd,stdout=PIPE,stderr=destination)
    print(sp)
    sp.communicate()
    output = pd.read_csv(os.path.join(tempdir,"pathways.csv"),index_col=0)
    output.index = output.index.astype(str)
    output.columns = output.columns.astype(str)
    output.index.name = 'name'
    return output

#TFs

def sele_pairs(evidence):
    ## evidence should be a list, among the list, only 'chipSeq','TFbindingMotif','curated' are effective selection.
    TF_pair_all = pd.read_csv(os.path.join(ROOT_DIR, 'Dataset/database.csv')) #Data from Luz Garcia-Alonso et al., Genome biology, 2019. 
    if 'chipSeq' in evidence:
        TF_pair_chipSeq = TF_pair_all.loc[TF_pair_all['is_evidence_chipSeq'] == True]
    else:
        TF_pair_chipSeq = pd.DataFrame()

    if 'TFbindingMotif' in evidence: 
        TF_pair_TFbindingMotif = TF_pair_all.loc[TF_pair_all['is_evidence_TFbindingMotif'] == True]
    else:
        TF_pair_chipSeq = pd.DataFrame()

    if 'curated' in evidence:
        TF_pair_curated = TF_pair_all.loc[TF_pair_all['is_evidence_curateddatabase'] == True]
    else:
        TF_pair_curated = {}

    result = pd.concat([TF_pair_curated, TF_pair_chipSeq,TF_pair_TFbindingMotif])
    return(result)

#### 2. select transcription factor and target gene pairs with co-expression evidences
def pair_cor_test(test, data_matrix):
    from scipy.stats import pearsonr
    TF_remain = []
    Target_remain = []
    corr_list = []
    dic_TF_target_pair = {}
    for i in range(0,test.shape[0]):
        dic_TF_target_pair[test.iloc[i]['TF'] + ':' + test.iloc[i]['target']] = '' 
    #test = x_A.loc[x_A['target'].isin(dic_module[list(set(KEGG_level2['name']))[i]])]
    
    for TF in list(set(test['TF'])):
        for Target in list(set(test['target'])):
            if (TF+':'+Target) in dic_TF_target_pair and (TF in data_matrix.columns.tolist()) and (Target in data_matrix.columns.tolist()) :
                corr, p = pearsonr(data_matrix[TF],data_matrix[Target])
                
                if p < 0.05 and (abs(corr)) > 0.2:  ##abs(corr) > 0.2 for both activator and supressor or corr > 0.2 for only activator
                
                    TF_remain.append(TF)
                    Target_remain.append(Target)
                    corr_list.append(corr)

    pair_cor = pd.DataFrame({"TF":TF_remain,"Target":Target_remain,"corr":corr_list})
    return(pair_cor)

#### 3. select transcription factor and target gene pairs with both co-expression evidence and high confident evidences from the orignal literatures.
def get_tf_target_pair(data_matrix,pathway_select,dic_module):
    pair_cor_df = pd.DataFrame({})
    for i in range(0,len(pathway_select)):
        x = sele_pairs(['chipSeq','TFbindingMotif','curated']) 
        x_AB = x.loc[x['score'].isin(['A','B'])]
        x_AB_module = x_AB.loc[x_AB['target'].isin(dic_module[pathway_select[i]])]
        pair_cor = pair_cor_test(x_AB_module, data_matrix)
        pair_cor_df = pd.concat([pair_cor_df, pair_cor])
        
    
    Resulting_pairs = pair_cor_df.drop_duplicates()
    return(Resulting_pairs)

#### 4. Estimate which transcription factors have the key regulation role for specific module.
def get_tfpairs_for_select_pathways(data_matrix,pathway_select,dic_module): ## remove parameter active = True
    
    Resulting_pairs = get_tf_target_pair(data_matrix,pathway_select,dic_module)
    #if active == True:
    #    Resulting_pairs = Resulting_pairs.loc[Resulting_pairs['corr'] > 0]
    #elif active == False:
    #    Resulting_pairs = Resulting_pairs.loc[Resulting_pairs['corr'] < 0]
        
    result_df = pd.DataFrame()
    
    for i in range(0,len(pathway_select)):
        pathway = pathway_select[i]
        test = Resulting_pairs.loc[Resulting_pairs['Target'].isin(dic_module[pathway])]
        pathway_list = []
        TF_list = []
        Occupancy_inPathway_list = []
        Occupancy_inTF_list = []
        Score_list = []

        for gene in list(set(test['TF'])):
            A = len(test.loc[test['TF'] == gene])
            B = test.shape[0] - A
            C = Resulting_pairs[Resulting_pairs['TF'] == gene].shape[0] - A
            D = Resulting_pairs.shape[0] - B - C + A
            oddsratio, pvalue = scipy.stats.fisher_exact([[A, B], [C, D]], alternative= 'greater')
            
            if pvalue < 0.05 and oddsratio > 0:
                Occupancy_inPathway = test.loc[test['TF'] == gene].shape[0]/test.shape[0]
                #Occupancy_inPathway estimate the percent of genes regulated by this transcription factor. The higher values means more enriched in this pathway.
                Occupancy_inTF = test.loc[test['TF'] == gene].shape[0]/Resulting_pairs.loc[Resulting_pairs['TF'] == gene].shape[0]
                #Occupancy_inTF estimate the percent of target genes for this transcription factor in this pathway. The higher values means more enriched in this pathway.
                Score = Occupancy_inPathway * Occupancy_inTF  
                pathway_list.append(pathway)
                TF_list.append(gene)
                Occupancy_inPathway_list.append(Occupancy_inPathway)
                Occupancy_inTF_list.append(Occupancy_inTF)
                Score_list.append(Score)

        result = pd.DataFrame({"Pathway":pathway_list, "TFs":TF_list, "Occupancy_inPathway":Occupancy_inPathway_list, "Occupancy_inTF":Occupancy_inTF_list, "Score":Score_list})
        result_df = pd.concat([result_df, result])
    return(result_df)
