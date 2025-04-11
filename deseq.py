# -*- coding: utf-8 -*-

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

#choose how to change expression: new exp = exp * fraction from fractions
fractions = [1.2,1.5,2]

#open and format input files
counts = pd.read_csv('counts.csv', index_col=0).transpose()
for col in counts.columns:
    counts[col]=counts[col].astype(int)
counts.index.name='sample'
coldata = pd.read_csv('coldata.csv', index_col=0)
with open('genes.txt') as f:
    genes = f.read().splitlines()

#do deseq on normal counts
inference = DefaultInference(n_cpus=6)
dds = DeseqDataSet(counts=counts, metadata=coldata,design="~condition",
    refit_cooks=True,inference=inference,)
dds.deseq2()
ds = DeseqStats(dds, contrast=["condition", "disease", "control"], 
                inference=inference)        
ds.summary()
results_df = ds.results_df
results_df.to_csv('normal_deseq_noshrink.csv')

#change gene expression and recalculate lfc and pvalue
#get disease cols id
disease_ids=[]
for row in range(len(dds.obs['condition'])):
    if dds.obs.iloc[row,0] == 'disease':
        disease_ids.append(row)
#get gene id
for gene in genes:
    result_df1 = pd.DataFrame(columns=['baseMean','log2FoldChange','lfcSE','stat',
                                      'pvalue','padj'])
    print(f'Changing {gene} expression...')
    gene_id = dds.var.index.to_list().index(gene)
    #change counts and normed counts for gene in disease
    for fraction in fractions:
        dds1=dds.copy()
        for row in disease_ids:
            dds1.X[row,gene_id]=dds.X[row,gene_id]*fraction
            dds1.layers['normed_counts'][row,gene_id]=dds1.layers['normed_counts'][row,gene_id]*fraction
        #recalculate mean, dispersion, lfc, pvalue
        dds1=DeseqDataSet(adata=dds1)
        dds1.varm['_normed_means']=dds1.layers["normed_counts"].mean(0)
        dds1.fit_genewise_dispersions()
        dds1.fit_MAP_dispersions()
        dds1.fit_LFC()
        ds1 = DeseqStats(dds1, contrast=["condition", "disease", "control"], 
                        inference=inference)        
        ds1.summary()
        result_df1 = pd.concat([result_df1,ds1.results_df.loc[[gene]]])
        #result_df1.to_csv(f'{gene}x{fraction}.csv')
    result_df1.to_csv(f'{gene}.csv')

