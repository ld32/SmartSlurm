[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)

# Using ssbatch with other pipeline managers

Because `ssbatch` accepts standard sbatch syntax, you can point Snakemake, Nextflow, or Cromwell at it as their submit command and get smart sizing for free.

## Snakemake

```bash
git clone https://github.com/ld32/SmartSlurm.git $HOME/SmartSlurm

mkdir -p SmartSlurmTest && cd SmartSlurmTest
source $HOME/SmartSlurm/bin/smartslurm.sh   # load the smartslurm command in this shell
smartslurm activate                          # ssbatch on PATH (run 'bash smartSlurmSetup.sh' once to make it permanent — see Installation and setup)

cp $HOME/SmartSlurm/bin/Snakefile .
cp $HOME/SmartSlurm/config/config.yaml .

module load conda/miniforge3/24.11.3-0
mamba env create --name snakemakeEnv --file $HOME/SmartSlurm/config/snakemakeEnv.yaml
conda activate snakemakeEnv

# Use ssbatch as the cluster submit command:
snakemake -p -j 999 --latency-wait=80 --cluster "ssbatch -t 100 --mem 1G -p short"

# With a specific Slurm account:
snakemake -p -j 999 --latency-wait=80 --cluster "ssbatch -A mySlurmAccount -t 100 --mem 1G"

checkRun
```

## Nextflow

```bash
git clone https://github.com/ld32/SmartSlurm.git $HOME/SmartSlurm

mkdir -p SmartSlurmTest && cd SmartSlurmTest
module load conda/miniforge3/24.11.3-0
mamba create -n nextflowEnv -c bioconda -y nextflow

source $HOME/SmartSlurm/bin/smartslurm.sh   # load the smartslurm command in this shell
smartslurm activate                          # also aliases sbatch -> ssbatch, which Nextflow calls directly
conda activate nextflowEnv

cp $HOME/SmartSlurm/bin/nextflow.nf .
cp $HOME/SmartSlurm/config/nextflow.config .

# For a specific Slurm account, edit nextflow.config and uncomment/set:
#   process.clusterOptions = '--account=mySlurmAcc'

nextflow run nextflow.nf -profile slurm
checkRun

# When finished, restore the system sbatch:
smartslurm deactivate
```

## Cromwell

> [!NOTE]
> Cromwell support is planned and not yet documented.

---


---

[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)
