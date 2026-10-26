
# Zambia Elephant Genomics — Analysis Log

## Date: October 21, 2026 — FASTQ Transfer

**Researcher:** Camille  
**Task:** Transfer focal Zambian elephant FASTQs to shared Slate-Project storage.

### Objective
Copy eight paired-end elephant sequencing samples (16 compressed FASTQ files) from personal Slate storage to the shared project directory so both collaborators can access the data.

### File locations
- Source: `/N/slate/cgmorale/elephant_fastqs/`
- Destination: `/N/project/ZambianElephants/elephant_fastqs/`

### Script
[`transfer_fastqs.sbatch`](../scripts/transfer_fastqs.sbatch)

### Commands executed
```bash
sbatch transfer_fastqs.sbatch
```

### SLURM job information
- Job ID: `8680666`
- Memory requested: 2 GB
- CPUs: 1
- Time limit: 12 hours
- Log file: `transfer_fastqs_8680666.log`

### Results
- Status: Submitted; completion not yet verified.
- Files expected: 16 (8 paired-end samples)
- Transfer method: `rsync`

### Next steps
- Confirm successful job completion.
- Verify that all FASTQ files transferred correctly.
- Confirm Xiaoyue can access the shared directory.
- Begin FastQC analysis.

## Date: 10/26/2026

**Researcher:**   Xiaoyue
**Task:** qc and multiqc for selected 10 samples(84 runs) from PRJEB77639, 
          qc and multiqc for 8 Zambian elephant samples(8 runs), 
          test: processing fastp on all Zambian Elephant samples
          test: processing fastp on all selected PRJEB77639 samples.

### Objective
    for evaluating raw-read quality and technical heterogeneity 
    -- including :
    read length, sequencing yield, GC content, duplication, adapter content, and other FastQC metrics.
    
    for testing whether standardized common read preprocessing should be applied before alignment.

### Script
- `qc_public.sbatch`
- `zambia_samples_qc.sbatch`
- `multiqc_public_selected.sh`
- `zambia_samples_multiqc.sbatch`
- `zambia_fastp_pilot.sbatch`
- `zambia_trimmed_fastqc.sbatch`
- 
### Commands executed
```
```
### SLURM job information

### Results
Raw-read FastQC and MultiQC were completed for 84 FASTQ files representing 10 selected public elephants and 16 FASTQ files representing eight focal Zambian elephants.

| Feature | Selected public dataset | Focal Zambian dataset | Interpretation |
|---|---|---|---|
| Dataset structure | 10 individuals; 84 FASTQs from multiple sequencing runs/libraries | 8 individuals; 16 FASTQs, one R1/R2 pair per individual | Public data have substantially greater run-level complexity |
| Read length before preprocessing | Predominantly 150–151 bp; several archived FASTQs averaged ~147–150 bp | Uniform 151 bp | Public data show greater archival/preprocessing heterogeneity |
| GC content | ~40–43% | ~41–43% | Highly comparable between datasets |
| Raw-read duplication | ~12–29% | ~26–30% | Focal samples cluster near the upper end of the public range; post-alignment duplicate rates will be more informative |
| Sequencing yield | Highly variable among individuals and runs | ~138–213 million read pairs per individual; relatively uniform | Sequencing effort is substantially more heterogeneous in the public dataset |
| R2-specific poly-G pattern | Present in a subset of libraries, notably QE0870 | Present in all 8 focal R2 libraries | Stronger cohort-wide sequencing-batch signature in the focal data |
| Focal poly-G abundance | — | ~0.39–1.03% of R2 reads before preprocessing | Zambia06 exceeded the FastQC failure threshold |
| Dominant focal adapter-associated signal | — | Nextera transposase-associated sequence in Zambia01, Zambia04, and Zambia07 | Signal reached ~5–6% by read position 140 in flagged libraries |

| QC outcome | Zambia | Public | Overall interpretation |
|---|---:|---:|---|
| Adapter WARN before | 6 | 4 | Present in both datasets, more widespread across focal individuals |
| Adapter WARN after | **0** | **0** | **Uniform fastp strategy successfully removed detectable adapter signal** |
| Overrepresented before | 8 R2; 1 FAIL | ~24 WARN/FAIL | Technical sequences present in both datasets |
| Overrepresented after | **3 R2 WARN; 0 FAIL** | **8 R2 WARN; 0 FAIL** | Strong reduction in both datasets |
| Per-base content | unchanged | unchanged | Not corrected by trimming |
| Per-tile quality | essentially unchanged | essentially unchanged | Not corrected by trimming |
| Length-distribution WARN | increased | **84/84** | Expected consequence of selective trimming |

| Metric | Selected public (42 pairs) | Zambian focal (8 pairs) | Interpretation |
|---|---:|---:|---|
| Read pairs/units analyzed | 42 | 8 | 84 vs 16 FASTQs |
| Reads retained | **99.188–100.000%** | **99.017–99.645%** | Very little read loss in both datasets |
| Median reads retained | **99.948%** | ~99.4% | Public data were minimally affected overall |
| Bases retained | **97.016–99.847%** | **97.125–97.691%** | Both datasets retained the great majority of sequence |
| Median bases retained | **99.672%** | ~97.25% | More trimming was required in the focal data |
| Adapter-trimmed reads / pair | 109,506–23,498,817 | ~16.1–31.8 M | Adapter-associated sequence was present in both datasets |
| Adapter FastQC WARN after | **0** | **0** | Detectable adapter warnings were eliminated in both |
| Overrepresented-sequence FAIL after | **0** | **0** | Strong technical signals were reduced |

| QC signal | Response to conservative fastp preprocessing | Interpretation / decision |
|---|---|---|
| Adapter content | Strongly improved | Consistent with removable library-associated technical sequence |
| R2 poly-G overrepresentation | Strongly reduced | Consistent with a sequencing-associated artifact targeted by poly-G trimming |
| Per-base sequence content | Largely unchanged | Not primarily caused by adapter/poly-G contamination; no justification for more aggressive trimming |
| Per-tile sequence quality | Largely unchanged | Consistent with flow-cell/tile-level sequencing effects that trimming cannot correct |
| Sequence-length distribution | More warnings after fastp | Expected consequence of selective trimming; not considered evidence of poorer read quality |
| Read/base retention | High in both datasets | Conservative preprocessing removed targeted technical sequence without extensive read loss |
| Sequencing-depth heterogeneity | Not addressed by fastp | Must be evaluated after mapping; potential public-data downsampling should be considered separately |


### Methodological decisions

The raw-read QC indicates technical heterogeneity between and within the public and focal datasets. Public sequencing runs vary substantially in total sequencing yield and show evidence of heterogeneous archival or preprocessing histories, whereas the focal Zambian data are more uniform in read length and sequencing yield. Both datasets nevertheless show broadly comparable GC composition and contain sequencing- or library-associated FastQC signals.

The focal fastp pilot demonstrates that conservative adapter and poly-G trimming can remove identifiable technical sequence while retaining more than 99% of reads. A final decision is still required on whether this preprocessing should be applied uniformly to the complete focal and public cohort. To evaluate this consistently, the same fastp parameters will next be tested across all 42 paired sequencing units (84 FASTQs) from the selected public dataset, followed by repeat FastQC/MultiQC assessment.

Sequencing depth will be treated separately from read preprocessing. Standardized trimming can reduce preprocessing heterogeneity but cannot eliminate depth-dependent bias in ROH detection. After alignment, actual mapped depth will therefore be estimated for all focal and public individuals before determining whether public genomes should be computationally downsampled.

For ROH sensitivity analysis, selected public elephants may be analyzed both at their original mapped depth and after downsampling to coverage comparable to the focal elephants. ROH estimates from these depth-matched datasets can then be compared with estimates from the same individuals at original depth to isolate the effect of sequencing depth while holding the analysis pipeline constant. Published ROH estimates from Pečnerová et al. will provide an additional external benchmark, but differences between published and newly generated ROH estimates cannot be attributed solely to sequencing depth because the published study used its own genotype filtering, imputation, and ROH-calling framework.

### Next steps
- Add Namibian elephant samples into the integrated framework and do comparison analysis.
- Compare pre- and post-fastp QC metrics for the public dataset and make a final cohort-wide preprocessing decision.
- Prepare the LoxAfr4 reference(✅) and construct a mapping manifest with sequencing-unit-specific read groups and biological-individual sample IDs.（✅）
- Align the eight focal and 10 selected public elephants using the same mapping and post-alignment QC workflow.
- Estimate actual mapped depth, coverage breadth, mapping rate, duplicate rate, and related post-alignment metrics for every individual.
- Determine depth-matching targets only after the observed focal coverage distribution is known.
- For selected public individuals, compare ROH estimates at original depth and at downsampled depth to quantify coverage-dependent ROH bias.
- Compare newly generated ROH estimates with published African-elephant ROH summaries while explicitly accounting for methodological differences.
