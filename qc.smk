#command run iterate snakemake: snakemake --profile d-slurm all_fastqc all -s snakefile_IF_pair -c 20 -j 

MIN_READS=1e6
[SMP,PAIR]=glob_wildcards("data/raw/fastq/{samplename}/{samplename}_{pair}.fastq.gz")
SMP=list(set(SMP))
PAIR=list(set(PAIR))

rule all_fastqc:
  input: expand("data/processed/fastqc/{samplename}_{pair}_fastqc.html", samplename=SMP, pair=PAIR)
    
rule fastqc:
   input:
    in1="data/raw/fastq/{samplename}/{samplename}_1.fastq.gz",
    in2="data/raw/fastq/{samplename}/{samplename}_2.fastq.gz"
   output:
    "data/processed/fastqc/{samplename}_1_fastqc.html",
    "data/processed/fastqc/{samplename}_2_fastqc.html",
   shell:
    "fastqc -o data/processed/fastqc/ {input.in1} {input.in2}"

rule multiqc:
  input: rules.all_fastqc.output
  output:
    "data/processed/multiqc/multiqc_report.html"
  shell:
    "multiqc -o data/processed/multiqc/ data/processed/fastqc/"

def bbduk_runtime_retries(wildcard, attempt):
  if (attempt == 1): return(10)
  elif (attempt == 2): return(30)
  elif (attempt == 3): return(60)
  else: return(300)

# hardcode qin at 33 since it seems to be getting confused
# also, these are mostly very fast so we can use the debug queue; however, we won't do that here since we can only run 5 jobs at a time then
rule bbduk:
  resources:
    runtime=bbduk_runtime_retries,
    mem_mb=8000
  input:
   in1="data/raw/fastq/{samplename}/{samplename}_1.fastq.gz",
   in2="data/raw/fastq/{samplename}/{samplename}_2.fastq.gz"
  output:
   out1="data/processed/bbduk/{samplename}_1.fastq.gz",
   out2="data/processed/bbduk/{samplename}_2.fastq.gz"
  shell:
    "bbduk.sh qin=33 in1={input.in1} in2={input.in2} out1={output.out1} out2={output.out2} qtrim=r trimq=20 minlength=60"

rule bbduk_all:
   input: expand("data/processed/bbduk/{samplename}_{pair}.fastq.gz", samplename=SMP, pair=[1,2])

def filt_read_counts_mem_retries(wildcard, attempt):
#  if (attempt == 1): return(8000)
  if (attempt == 1): return(16000)
  elif (attempt == 2): return(24000)
  elif (attempt == 3): return(64000)
  else: return(112000)

rule filt_read_counts:
   resources:
      mem_mb=filt_read_counts_mem_retries
   input: "data/processed/bbduk/{samplename}_{pair}.fastq.gz"
   output: "data/processed/filt_read_counts/{samplename}_{pair}.txt"
   shell: "scripts/count_fastq.sh {input} > {output}" 

rule filt_read_counts_all:
   input: expand("data/processed/filt_read_counts/{samplename}_{pair}.txt", samplename=SMP, pair=[1,2])
   output: "data/processed/filt_read_counts/overall.txt"
   shell: "cat {input} | sort -gk2 > {output}"

rule passed_qc:
   input: "data/processed/filt_read_counts/overall.txt"
   output: directory("data/processed/passed_qc/")
   run:
      samples=dict()
      with open(input[0]) as fh:
        for line in fh:
           row=line[:-1].split(' ')
           samples[row[0]]=int(row[1])
      passed=[s for s in samples.keys() if samples[s] >= MIN_READS]
      for p in passed:
         in_p = f"../bbduk/{p}.fastq.gz"
         shell("export LASTCD=$PWD; mkdir -p {output}; cd {output}; ln -s {in_p}; cd $LASTCD") 
      print(f"{len(passed)} / {len(samples.keys())} had >{MIN_READS} after filtering")
