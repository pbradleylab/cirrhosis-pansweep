### Rules to run MIDAS2 on cirrhosis data

import tempfile
import os
from pathlib import Path


## SET AS APPROPRIATE FOR YOUR ENVIRONMENT
TMPDIR="/fs/ess/scratch/PAS2276/"
SCRATCHDIR="/fs/ess/scratch/PAS2276/"
(SMPS, PAIRS) = glob_wildcards("data/processed/passed_qc/{sample}_{pair}.fastq.gz")
SAMPLES = set(list(SMPS))
# print(list(SAMPLES)[0:5])

(PREFIXES, SPECIES, FILENAMES) = glob_wildcards("data/raw/uhgg_catalogue/{prefix}/{species}/{filename}")
SPECIES_UNIQ = set([s for s in SPECIES if s.startswith("MGYG")])

rule midas_coverage:
  conda: "midas2"
  threads: 6
  resources:
    runtime=600
  input:
    r1="data/processed/passed_qc/{sample}_1.fastq.gz",
    r2="data/processed/passed_qc/{sample}_2.fastq.gz"
  output: "data/processed/midas_coverage/{sample}/species/species_profile.tsv"
  shell: """
      midas2 run_species \
             --sample_name {wildcards.sample} \
             -1 {input.r1} \
             -2 {input.r2} \
             --midasdb_name uhgg \
             --midasdb_dir data/raw/midasdb_uhgg \
             --num_cores {threads} \
             data/processed/midas_coverage
    """

rule all_coverage:
  input: expand("data/processed/midas_coverage/{sample}/species/species_profile.tsv", sample=SAMPLES)
  output: "data/processed/midas_coverage/.coverage_done.txt"
  run: shell("touch {output}")

rule sample_manifest:
  output: "data/processed/sample_manifest.tsv"
  run:
    with(open(output[0], 'w')) as fh:
      fh.write("sample_name\tmidas_outdir\n")
      for smp in SAMPLES:
        #fh.write(f"{smp}\tdata/processed/midas_coverage/{smp}\n")
        fh.write(f"{smp}\tdata/processed/midas_coverage/\n")

rule midas_combine:
  conda: "midas2"
  input:
    cd="data/processed/midas_coverage/.coverage_done.txt",
    sm="data/processed/sample_manifest.tsv"
  output: directory("data/processed/midas_merged")
  shell: "midas2 merge_species --samples_list {input.sm} --min_cov 2 {output}"

rule copy_db_to_scratch:
  input: "data/raw/midasdb_uhgg"
  output: directory(SCRATCHDIR+"midasdb_uhgg")
  shell: "rsync -avz {input} {output}"

# Todo: add retries programmatically instead of manually
rule midas_snv:
  conda: "midas2"
  input:
    r1="data/processed/passed_qc/{sample}_1.fastq.gz",
    r2="data/processed/passed_qc/{sample}_2.fastq.gz",
    db=SCRATCHDIR+"midasdb_uhgg",
    sp_prof="data/processed/midas_coverage/{sample}/species/species_profile.tsv"
  output: "data/processed/midas_coverage/{sample}/genes/genes_summary.tsv"
  #resources: mem_mb=64000, runtime=900 # try this first
  resources: mem_mb=256000, runtime=4320, partition="hugemem"
  threads: 12
  shell: """
    mkdir -p /fs/ess/scratch/PAS2276/midas_snv/ &&
    rsync -avz data/processed/midas_coverage/{wildcards.sample} /fs/ess/scratch/PAS2276/midas_snv/ &&
    midas2 run_genes \
      --sample_name {wildcards.sample} \
      -1 {input.r1} \
      -2 {input.r2} \
      --midasdb_name uhgg \
      --midasdb_dir {input.db} \
      --select_by median_marker_coverage,unique_fraction_covered \
      --select_threshold=2,0.5 \
      --num_cores {threads} \
      /fs/ess/scratch/PAS2276/midas_snv/ &&
    rsync -avz /fs/ess/scratch/PAS2276/midas_snv/{wildcards.sample} data/processed/midas_coverage/
    """

# don't re-run unless asked
rule all_snv:
  input: expand("data/processed/midas_coverage/{sample}/genes/genes_summary.tsv", sample=SAMPLES)
  output: "data/processed/midas_coverage/.snv_done.txt"
  run: shell("touch {output}")

rule midas_snv_merge:
  conda: "midas2"
  input:
    md="data/processed/midas_coverage/.snv_done.txt",
    sm="data/processed/sample_manifest.tsv"
  output: directory("data/processed/midas_snv_merge/")
  threads: 4
  resources: mem_mb=8000 
  shell:
    """ \
    midas2 merge_genes \
    --samples_list {input.sm} \
    --midasdb_name uhgg \
    --midasdb_dir data/raw/midasdb_uhgg \
    --num_cores {threads} \
    {output}
  """

