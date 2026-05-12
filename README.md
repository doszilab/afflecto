# AFFlecto-Nextflow

Automated Nextflow workflow for generating conformational ensembles of intrinsically disordered proteins from AlphaFold structures.

The workflow supports:

- Automatic disorder region classification from AlphaFold structures (based on pLDDT)
- Optional user-provided region classification TSV files
- Containerized ensemble generation with FReSa
- Parallel execution with Nextflow
- Local execution or HPC/Slurm deployment
- Resume/restart-safe execution

---

# Overview

AFFlecto-Nextflow processes one or multiple AlphaFold protein structures (`.pdb`, `.cif`, `.mmcif`) and generates conformational ensembles using the AFFlecto/FReSa framework.

# Sampling strategy

AFFlecto uses two complementary conformational sampling strategies to generate realistic ensembles of intrinsically disordered proteins. The single-residue-based sampling (SRS) strategy generates conformations independently for each residue, while the three-residue-based sampling (TRS) strategy incorporates local sequence and structural context, making it particularly suitable for partially structured or conditionally folded regions. In AlphaFold structures, such partially structured regions are identified as secondary structural elements with limited tertiary contacts, and are automatically sampled using the TRS strategy to preserve their local structural preferences during ensemble generation. The TRS sampling strategy can be disabled using --use_TRS false, in which case all flexible regions, including partially structured regions, are sampled exclusively using the SRS strategy.


Two execution modes are supported:

## 1. Automatic region classification

The workflow automatically:

- reads AlphaFold structures
- detects flexible regions from pLDDT
- runs DSSP to identify and filter out falsely predicted secondary structural elements from AlphaFold structures
- classifies the identified disordered regions as tails, loops, or linkers
- launches ensemble generation

## 2. User-provided region classification

The user provides `.tsv` region files.

In this case:

- DSSP and automatic classification are skipped
- provided TSV annotations are validated

If TRS regions are manually defined in the TSV file, all TRS regions must be fully contained within non-rigid (disordered) regions such as tail, loop, linker, or flexible_prot.

---

# Installation

## Requirements

Install:

- Java 17+
- Nextflow
- Apptainer

---

## Install Java

```bash
sudo apt update
sudo apt install -y openjdk-17-jdk
```

Check:

```bash
java -version
```

---

## Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
```

Check:

```bash
nextflow -version
```

---

## Install Apptainer

Ubuntu:

```bash
sudo apt install -y apptainer
```

Check:

```bash
apptainer --version
```

---

# Running the workflow

---

# 1. Standard mode: automatic region classification

In the standard mode, the user only provides a directory containing input structure files. AFFlecto automatically identifies and classifies disordered regions, and launches ensemble generation.

Example command:

```bash
nextflow run doszilab/afflecto \
  -r v0.1.1 \
  -profile apptainer,standard \
  --input_dir /path/to/input_dir \
  --outdir results_standard \
  --n_conformers 1000 \
  --fresa_threads 2 \
  --fresa_timeout_min 20
```

This mode:

- computes disordered regions automatically with classification
- generates region TSV files
- runs ensemble generation


# 2. TSV mode: user-provided region classification

In TSV mode, the user provides precomputed region classification files. In this case, AFFlecto skips the automatic region classification step and directly launches ensemble generation.

Example command:

```bash
nextflow run doszilab/afflecto \
  -r v0.1.1 \
  -profile apptainer,tsv \
  --input_dir /path/to/input_dir_tsv \
  --outdir results_tsv \
  --n_conformers 1000 \
  --fresa_threads 2 \
  --fresa_timeout_min 20
```

This mode:

- uses user-provided region classification TSV files
- runs ensemble generation

---

## TSV region files

Each input structure file must have a corresponding TSV file with the same protein identifier.

Example:

```text
AF-P35222.pdb
AF-P35222.tsv
```

Example TSV content:

```text
AF-P35222	1	88	tail
AF-P35222	89	553	rigid
AF-P35222	554	559	loop
AF-P35222	560	612	rigid
AF-P35222	4	28	TRS
```

The protein identifier in the TSV file must match the structure filename.

Supported region types:

- rigid
- tail
- loop
- linker
- flexible_prot
- TRS



---

# Resume interrupted runs

One of the main advantages of Nextflow is restart-safe execution.

If the workflow stops or crashes, resume with:

```bash
  -resume
```

Only unfinished tasks are rerun.

---

# Input files

---

## Structure inputs

Supported formats:

- `.pdb`
- `.cif`
- `.mmcif`

Requirements:

- single-chain structures
- AlphaFold-style pLDDT values stored in B-factors

---

# Main parameters

| Parameter | Description |
|---|---|
| `-profile` | Execution profile(s), e.g. `apptainer,standard` or `apptainer,tsv` || `-resume` | Resume a previously interrupted or completed workflow execution |
| `--input_dir` | Directory containing input structure files (`.pdb`, `.cif`, `.mmcif`) |
| `--outdir` | Output directory |
| `--n_conformers` | Number of conformations generated per protein |
| `--fresa_threads` | Number of CPU threads used internally by FReSa |
| `--fresa_timeout_min` | Maximum allowed runtime (in minutes) for one FReSa job |
| `--use_TRS` | Enable or disable TRS-based sampling (`true` or `false`). When disabled, all regions are sampled using the SRS strategy (default: `true`) |

# Output structure

Example:

```text
ensembles/AF-P35222/
├── AF-P35222.fresa_meta.json
├── AF-P35222.params_meta.json
├── AF-P35222.pdb
├── AF-P35222.region_meta.json
├── AF-P35222.tsv
├── fresa_stdout.log
├── fresa_stderr.log
└── PDBs/
```

---

# Nextflow working directories

During execution, Nextflow creates:

```text
work/
.nextflow/
.nextflow.log
```

## Important

### `work/`

Contains task cache and intermediate files.

Required for:

```bash
-resume
```

Do not delete if you want resume support.

---

### `.nextflow.log`

Workflow execution logs.

Safe to delete.

---

# Cleaning workflow cache

Preview cleanup:

```bash
nextflow clean -n
```

Clean cache:

```bash
nextflow clean -f
```

---

# Parallelization model

AFFlecto-Nextflow uses two levels of parallelization:

## Protein-level parallelization

Handled by Nextflow.

Multiple proteins can run simultaneously.

---

## Internal FReSa threading

Handled by FReSa itself:

```text
nb_threads
```


---

# Development status

Current development stage:

- basic Nextflow structure
- container integration
- automatic classification
- TSV mode
- local execution
- Slurm integration (in progress)

---

# Citation

If you use AFFlecto-Nextflow in scientific work, please cite:

## AFFlecto

Pajkos M., Clerc I., Zanon C., Bernadó P., Cortés J.  
**AFflecto: A web server to generate conformational ensembles of flexible proteins from AlphaFold models**  
*Journal of Molecular Biology* (2025), 437(15):169003.  

---

## FReSa

Estaña A., Sibille N., Delaforge E., Vaisset M., Cortés J., Bernadó P.  
**Realistic ensemble models of intrinsically disordered proteins using a structure-encoding coil database**  
*Structure* (2019), 27(2), 381–391.e2.

