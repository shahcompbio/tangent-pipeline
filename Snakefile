import os
import subprocess
import pandas as pd

# settings
wildcard_constraints:
    type = 'TUMOR|NORMAL',

PATIENTS = [s.strip() for s in open(config['patientlist']).readlines()]
PATIENTS = [s for s in PATIENTS if not s.startswith('DO')]
# PATIENTS = PATIENTS[:1] ##@##
# PATIENTS = PATIENTS[1:] ##@##
# PATIENTS.remove('DO46325') # bam not present
# PATIENTS.remove('DO46327') # bam not present
# PATIENTS.remove('DO46328') # bam not present
# PATIENTS.remove('DO46329') # bam not present
# PATIENTS.remove('DO46330') # bam not present
TYPES = ['TUMOR', 'NORMAL']

if not os.path.exists(config['log_dir']): subprocess.run(f'mkdir -p {config["log_dir"]}', shell=True)
if not os.path.exists(config['tmp_dir']): subprocess.run(f'mkdir -p {config["tmp_dir"]}', shell=True)
if config['genome_version'] == 'hg19':
    ref_fasta = "/juno/work/shah/isabl_data_lake/assemblies/GRCh37/GRCh37-lite/GRCh37-lite.fa"

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


# rules
rule all:
    input:
        expand('results/doc/{patient}.{type}.bam_summary.sample_summary', patient=PATIENTS, type=TYPES,),
        'results/intervals/genome.intervals',
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
        'ln -s $(realpath {input.normal_bam}) $(realpath {output.normal_bam}) && '
        'ln -s $(realpath {input.normal_bam}.bai) $(realpath {output.normal_bai}) && '
        'ln -s $(realpath {input.tumor_bam}) $(realpath {output.tumor_bam}) && '
        'ln -s $(realpath {input.tumor_bam}.bai) $(realpath {output.tumor_bai})'

rule create_intervals:
    input:
        ref_fasta,
    output:
        intervals = 'results/intervals/genome.intervals',
    run:
        import wgs_analysis.refgenome as refgenome
        refgenome.set_genome_version(config['genome_version'])
        bin_size = 500000
        assert refgenome.info.chromosome_lengths.shape[0] > 0
        with open(output.intervals, 'w') as out:
            for chrom, length in refgenome.info.chromosome_lengths.items():
                for i in range(0, length, bin_size):
                    start = i + 1
                    end = min(length, i + bin_size)
                    out.write(f'{chrom}:{start}-{end}\n')

rule depth_of_coverage:
    input:
        bam = 'results/bam/{patient}.{type}.bam',
        intervals = 'results/intervals/genome.intervals',
    output:
        doc = 'results/doc/{patient}.{type}.bam_summary.sample_summary',
        cov_counts = 'results/doc/{patient}.{type}.bam_summary.sample_cumulative_coverage_counts',
        int_summary = 'results/doc/{patient}.{type}.bam_summary.sample_interval_summary',
        int_stats = 'results/doc/{patient}.{type}.bam_summary.sample_interval_statistics',
        sample_stats = 'results/doc/{patient}.{type}.bam_summary.sample_statistics',
        cov_proportions = 'results/doc/{patient}.{type}.bam_summary.sample_cumulative_coverage_proportions',
    params:
        prefix = lambda w: 'results/doc/{w.patient}.{w.type}.bam',
        ref = ref_fasta,
    singularity:
        '/juno/work/shah/users/chois7/singularity/sif/gatk4.sif',
    shell:
        'gatk DepthOfCoverage '
        '-R {params.ref} -L {input.intervals} '
        '-I {input.bam} -O {params.prefix} '
        '--omit-depth-output-at-each-base true '
