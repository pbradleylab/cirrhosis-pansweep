
srr_ids = []
with open("SRR_cirrhosis.txt") as fh:
    for line in fh:
        srr_ids.append(line[:-1])
print(len(srr_ids))

def get_all(wc):
    return(expand("fastq/{run}/{run}_{pair}.fastq.gz", run=srr_ids, pair=[1, 2]))

print(len(get_all(None)))

rule fasterq_dump:
    resources:
        load=25,
        runtime=60
    threads: 6 
    input: "sra/{run}/{run}.sra"
    output: r1="fastq/{run}/{run}_1.fastq", r2="fastq/{run}/{run}_2.fastq"
    shell: "cd sra; fasterq-dump --outdir=\"../fastq/{wildcards.run}/\" {wildcards.run} -e {threads}"

rule prefetch:
    output: "sra/{run}/{run}.sra"
    shell: "prefetch {wildcards.run} -O ./sra/"

rule gzip:
    input: "fastq/{run}/{run}_{pair}.fastq"
    output: "fastq/{run}/{run}_{pair}.fastq.gz"
    shell: "gzip {input}"

rule all_data:
    input: get_all
   
