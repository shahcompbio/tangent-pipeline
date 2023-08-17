import os
import subprocess

if not os.path.exists(config['log_dir']): subprocess.run(f'mkdir -p {config["log_dir"]}', shell=True)
if not os.path.exists(config['tmp_dir']): subprocess.run(f'mkdir -p {config["tmp_dir"]}', shell=True)
if config['genome_version'] == 'hg19':
    ref_fasta = "juno/work/shah/isabl_data_lake/assemblies/GRCh37/GRCh37-lite/GRCh37-lite.fa"

PATIENTS = ['

def _fetch_bam_metadata(patient):
    meta = pd.read_table(config['metadata'])
    meta = meta[meta['isabl_patient_id']==patient]
    meta = meta[meta['result_type']=='bam']
    return meta

def _get_normal_bam(wildcards):
    meta = _fetch_bam_metadata(wildcards.patient)
    meta = meta[meta['sample_category']=='NORMAL']
    paths = meta['result_filepath'].values
    assert len(paths) == 1, f'paths:{paths} for {wildcards.patient}'
    path = paths[0]
    return path

def _get_tumor_bam(wildcards):
    meta = _fetch_bam_metadata(wildcards.patient)
    meta = meta[meta['sample_category']=='TUMOR']
    paths = meta['result_filepath'].values
    assert len(paths) == 1, f'paths:{paths} for {wildcards.patient}'
    path = paths[0]
    return path

rule all:
    input:
        expand('results/bam/{patient}.NORMAL.bam.bai', patient=PATIENTS),
        expand('results/bam/{patient}.TUMOR.bam.bai', patient=PATIENTS),

rule symbolic_link:
    input:
        normal_bam = _get_normal_bam,
        tumor_bam = _get_tumor_bam,
    output:
        normal_bam = 'results/bam/{patient}.NORMAL.bam',
        normal_bai = 'results/bam/{patient}.NORMAL.bam.bai',
        tumor_bam = 'results/bam/{patient}.TUMOR.bam',
        tumor_bai = 'results/bam/{patient}.TUMOR.bam.bai',
    shell:
        'ln -s $(realpath {input.normal_bam}) $(realpath {output.normal_bam} && '
        'ln -s $(realpath {input.normal_bai}.bai) $(realpath {output.normal_bai} && '
        'ln -s $(realpath {input.tumor_bam}) $(realpath {output.tumor_bam} && '
        'ln -s $(realpath {input.tumor_bai}.bai) $(realpath {output.tumor_bai}'
        
