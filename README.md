# Zambian Elephant Genomics: Inbreeding and Population Ancestry

A collaborative bioinformatics project examining **runs of homozygosity (ROH)** and **forest–savanna ancestry** in African savanna elephants (*Loxodonta africana*) sampled in Zambia.

**Collaborators:** Camille Morales and Xiaoyue Tang  
**Project PI:** Matthew Hahn  
**Institution:** Indiana University  
**Project status:** In progress; analytical methods and storage locations are being finalized.

## Research questions

1. **Genomic inbreeding:** What do the number, length, and genomic proportion of ROH reveal about homozygosity in the focal Zambian elephants? How do these estimates compare *descriptively* with published results from other African elephant populations?
2. **Population ancestry:** How do the focal Zambian elephants relate to forest and savanna elephant reference populations? Is there detectable forest elephant-associated ancestry in any focal individuals?

ROH and ancestry are complementary analyses, not interchangeable measures. Ancestry components inferred by clustering methods are not, on their own, proof of recent hybridization.

## Samples and scope

### Focal cohort

- **Planned cohort:** 10 Zambian savanna elephants, sampled for whole-genome sequencing (WGS).
- **Expected sequencing depth:** approximately 5–8× per individual.
- **Source data:** paired-end Illumina FASTQ files. The sample manifest should record which of the planned samples have complete data available.
- **ROH estimates will be calculated for the focal Zambian elephants only.**

### Public reference panel

Ten publicly available elephants were selected from the sequencing data associated with **Pečnerová et al. (2026)**, ENA project [PRJEB77639](https://www.ebi.ac.uk/ena/browser/view/PRJEB77639). These individuals will be included in **population structure and ancestry analyses**, not in new ROH calculations.

| Sample ID | Population / location | Intended role |
|---|---|---|
| KR0027 | Kruger | Savanna reference |
| KR0068 | Kruger | Savanna reference |
| KR0116 | Kruger | Savanna reference |
| KR0138 | Kruger | Savanna reference |
| LO3503 | Lopé | Forest reference |
| LO3505 | Lopé | Forest reference |
| DS1506 | Dzanga Sangha | Forest reference |
| DS1507 | Dzanga Sangha | Forest reference |
| QE0870 | Queen Elizabeth NP, Mweya | Admixed comparison |
| QE2798 | Queen Elizabeth NP, Ishasha | Admixed comparison |

The two admixed elephants are included as **comparison individuals**, not as separate ancestral source populations. A sample manifest should link each biological sample ID to **all associated ENA sequencing runs**, locations, published ancestry estimates, and QC status; one sequencing run does not necessarily represent a complete individual's data.

## Planned analysis workflow

Methods listed below are **proposed** and may be revised following initial data quality assessment.

| Stage | Planned approach | Main outputs |
|---|---|---|
| 1. Inventory and QC | Verify samples and FASTQs; FastQC and MultiQC | Sample manifest, read-quality reports |
| 2. Read alignment | BWA-MEM against LoxAfr4; SAMtools for sorting, duplicate marking/handling, and indexing | Aligned BAM/CRAM files and mapping statistics |
| 3. Coverage assessment | Examine depth, mapping quality, missingness, and callable genomic regions | Per-sample QC summaries |
| 4. Population structure | Evaluate genotype-likelihood-based methods (e.g., ANGSD and PCAngsd) across focal and public genomes | PCA and population-structure plots |
| 5. Ancestry estimation | Evaluate NGSadmix; explore K = 2 initially, with additional settings/replicates as warranted | Individual ancestry estimates and bar plots |
| 6. ROH inference | Analyze focal genomes only; select and document an approach suitable for 5–8× depth (e.g., BCFtools roh or PLINK) | ROH segments, ROH length distributions, and FROH |
| 7. Contextual comparison | Compare Zambian ROH summaries with published estimates in Pečnerová et al. (2026) | Descriptive tables and interpretation of comparability limits |

**Important analysis distinctions:**

- Prepare **appropriate, potentially different filtering datasets** for ancestry inference and ROH detection. Do not automatically reuse an LD-pruned ancestry SNP panel for ROH calling.
- When calculating **FROH**, record both the total accepted ROH length and the autosomal genomic length considered callable/eligible in the analysis.
- Comparisons against published ROH results are **descriptive**, because sequencing depth, reference genome, imputation, filtering, and ROH-calling methods may differ.
- The modest reference-panel size and possible reference-mapping bias (from using a savanna-derived assembly for both species) may limit the detection and interpretation of subtle forest ancestry.

## Where files belong

**GitHub (this repository):** scripts, small sample manifests, configurations, documented decisions, run summaries, figures, and manuscript/report drafts as appropriate.

**IU storage:** raw FASTQs, large alignment files, variant files, intermediate datasets, and large output files. The precise shared **Slate-Project** allocation and path should be confirmed before adding paths to scripts or documentation.

**IU Quartz:** computational jobs using the project's allocation, when appropriate. The currently identified SLURM account is `c02431`:

```bash
#SBATCH -A c02431
```

This is a **compute account**, not a filesystem path or confirmation of Slate-Project storage access. Follow IU policies for large data transfers and jobs; avoid intensive transfers on login nodes if a transfer/compute node is required.

**Do not commit** FASTQ, BAM, CRAM, VCF/BCF, large temporary files, credentials, or unpublished sensitive metadata to GitHub. Keep the repository private unless public release is approved by the project PI and collaborators.

## Suggested repository structure

```text
zambia-elephant-genomics/
├── README.md
├── docs/
│   ├── analysis_log.md       # Dated summaries of runs and results
│   └── decisions.md          # Methodological decisions and rationale
├── metadata/
│   ├── focal_samples.csv     # Sample IDs, read pairs, QC status
│   └── public_references.csv # Reference IDs, ENA runs, population labels
├── scripts/                 # Shell, Python, R, and SLURM scripts
├── config/                  # Parameter files and environment records
├── notebooks/               # Exploratory analyses and plotting
├── figures/                 # Small, report-ready figures
└── .gitignore               # Excludes large data and secrets
```

*This is a suggested layout; files and directories will be created as the workflow develops.*

## Documenting each analysis

For every meaningful computational run, record the following in `docs/analysis_log.md` (or link to a tracked job record):

- **Date, analyst, and purpose** of the run
- **Input sample IDs and file paths** (prefer stable project-relative paths over personal paths)
- **Command/script and Git commit** used
- **Software versions**, reference genome assembly, and parameter/filter thresholds
- **SLURM job ID** and relevant CPU, memory, and time settings
- **Output locations**, QC observations, errors, and whether the run should be rerun

Record important methodological choices (such as the ROH caller, reference panel, depth filters, and ancestry K values) in `docs/decisions.md`, including alternatives considered and reasons for the final choice.

## Collaboration workflow

1. Use **GitHub Issues** to assign tasks, flag problems, and discuss open methodological questions.
2. Keep scripts under version control; use separate branches for changes likely to conflict and review each other's updates before merging when practical.
3. Avoid editing the same notebook simultaneously. For reusable pipeline steps, prefer standalone scripts over notebook-only code.
4. Update the analysis log when a result is generated, not only when preparing the final report.
5. Keep data paths and parameters configurable rather than hard-coding personal directories.

## Expected deliverables

- A documented, quality-checked sample and reference manifest
- Reproducible scripts and recorded run parameters
- PCA and ancestry-proportion visualizations for focal and selected public elephants
- Individual ROH metrics and length distributions for the focal Zambian elephants
- A cautious comparison of Zambian ROH patterns with published African elephant results
- A final report describing findings, limitations, and possible follow-up analyses

## Key references

- Pečnerová, P., et al. (2026). *The genomic impact of population connectivity and decline in Africa's elephants*. **Nature Communications**. https://doi.org/10.1038/s41467-026-71262-w
- Skotte, L., Korneliussen, T. S., & Albrechtsen, A. (2013). *Estimating individual admixture proportions from next generation sequencing data*. **Genetics**. https://doi.org/10.1534/genetics.113.154138
- Meisner, J., & Albrechtsen, A. (2018). *Inferring population structure and admixture proportions in low-depth NGS data*. **Genetics**. https://doi.org/10.1534/genetics.118.301336

---

**Status note:** This README describes the intended study design. Exact input manifests, computing procedures, quality thresholds, and final software choices should be updated as they are confirmed.
