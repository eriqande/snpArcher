import sys

import pandas as pd


vcf_table = pd.read_csv(config["vcf_list"]).set_index(["genome", "final_prefix"])
vcf_table_noindex = pd.read_csv(config["vcf_list"])




refgenome_list = vcf_table_noindex['genome'].tolist()
final_prefix_list = vcf_table_noindex['final_prefix'].tolist()

def get_vcf_path(wildcards):
    return vcf_table.loc[ (wildcards.refGenome, wildcards.prefix), "vcf_path" ]

def get_file_type(wildcards, input):
    if input.vcf.endswith(".bcf"):
        return "--bcf"
    elif input.vcf.endswith(".vcf.gz"):
        return "--gzvcf"
    elif input.vcf.endswith(".vcf"):
        return "--vcf"
    else:
        raise ValueError("Unsupported file type. The file must end with either '.bcf' or '.vcf.gz'.")

def get_sample_info_path(wildcards):
    return vcf_table.loc[ (wildcards.refGenome, wildcards.prefix), "sample_info_path" ]


def get_excluded_scaffolds(wildcards):
    return vcf_table.loc[ (wildcards.refGenome, wildcards.prefix), "excluded_scaffolds_path" ]

def get_coords_if_available(wildcards):
    si_path = get_sample_info_path(wildcards)
    samples = pd.read_csv(si_path)
    if 'lat' in samples.columns and 'long' in samples.columns:
        return "results/{refGenome}/QC/{prefix}.coords.txt"
    return []

def check_contig_names(fai, touch_file):
    dffai = pd.read_table(fai, sep='\t', header = None)
    fai_result=pd.to_numeric(dffai[0], errors='coerce').notnull().all()
    if fai_result==True:
        print("QC plots not generated because contig names are numeric and plink does not accept numeric contig names")
    elif fai_result==False:
        with open(touch_file, "w") as writer:
            writer.write("contigs are strings")
