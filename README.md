this is code and input files to analyze if deseq will find the change to raw expression of one gene significant and consider changed gene a deg
it uses pydeseq2 (https://github.com/owkin/PyDESeq2/tree/main) for easier coding, the result is consistent with deseq2 R package (with slight value variation, tested on this input data)

input files include:
1) counts.csv = count matrix with rnaseq data
2) coldata.csv = metadata for DEG experiment. Column "condition" denotes random "disease" and "control" samples (3x3 genotypes)
3) genes.txt = list of genes to examine. gene ids must be from gencode.v36.annotation.gtf

deseq.py will:
1) run standard deseq experiment on counts.csv using metadata from coldata.csv. saves the result_df as a file 'normal_deseq_noshrink.csv' in the current directory. Here and later lgc shrink is NOT performed on the final lgfcs.
2) change raw counts of each gene from genes.txt by multipling it by 'fraction' (you can change fractions list in code)
3) recalculate mean, dispersion, lfc and pvalue for the changed gene using the same size coefs and trend coefs as before. Saves results to {gene}.csv in the current directory
