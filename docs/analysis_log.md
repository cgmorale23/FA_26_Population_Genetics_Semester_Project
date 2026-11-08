
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

## Date: 10/27/2026

**Researcher:**   Xiaoyue
**Task:**  qc and multiqc for 4 Namibian elephant samples(4 runs), 
          test: processing fastp on all Namibian Elephant samples and re-check qualification

### Objective
    for evaluating raw-read quality and (technical heterogeneity) 
    -- including :
    read length, sequencing yield, GC content, duplication, adapter content, and other FastQC metrics.
    
    for testing whether standardized common read preprocessing should be applied before alignment.

### Script
- 'namibia_raw_fastqc.sbatch'
- 'namibia_raw_multiqc.sbatch'
- 'namibia_fastp.sbatch'
- 'namibia_post_fastp_fastqc.sbatch'


### Results
| FastQC metric | Raw | Post-fastp | Change | Interpretation |
|---|---:|---:|---|---|
| Per-base sequence content | **8 FAIL** | **8 FAIL** | No change | Systematic sequence-composition pattern; not caused by adapter/poly-G |
| Per-sequence GC content | **3 WARN** | **3 WARN** | No change | Mild GC-distribution signal persists |
| Overrepresented sequences | **4 WARN** | **0** | **Resolved** | R2-specific technical sequences were successfully removed |
| Adapter Content | 0 | 0 | No issue detected | No FastQC-level adapter warning in either condition |
| Per-tile sequence quality | 0 | 0 | No issue detected | No evident tile-quality problem |
| Sequence Length Distribution | **8 WARN** | **8 WARN** | No change in status | Already variable in archived raw FASTQs; therefore not solely introduced by fastp |


| Metric | Namibia |
|---|---:|
| Individuals | 4 |
| PE units | 4 |
| FASTQs | 8 |
| Reads retained | **99.903–99.932%** |
| Median reads retained | **99.918%** |
| Bases retained | **99.692–99.735%** |
| Median bases retained | **99.701%** |
| Adapter-trimmed reads / pair | 262,752–1,158,816 |
| Raw R2 poly-G | **4/4 WARN, 0.107–0.157%** |
| Post-fastp overrepresented sequences | **0 WARN / 0 FAIL** |


 
## Date: 11/2/2026

**Researcher:**   Xiaoyue
**Task:**  Sensitivity analysis on QC on three branches

### Objective
While reviewing our initial QC pipeline against the workflow used by Pečnerová et al.(STAMPEDE workflow), It is noticed that our initial preprocessing used fastp rather than AdapterRemoval, and it included poly-G trimming but did not merge overlapping paired-end reads. In contrast, the published PALEOMIX workflow used AdapterRemoval with overlapping-read collapsing and did not explicitly include poly-G trimming.

To evaluate whether these differences affect downstream results, we set up three parallel preprocessing sensitivity experiments:

Current fastp (merge OFF): Adapter and poly-G trimming, retaining reads as paired-end.

fastp (merge ON): The same adapter and poly-G trimming, with overlapping read pairs merged into single reads. Unmerged pairs are retained.

AdapterRemoval (Pečnerová-style): Adapter trimming and overlapping-read collapsing using the settings specified in their published PALEOMIX configuration, without an additional poly-G trimming step.

The goal is to assess whether the preprocessing choices influence mapping statistics, effective coverage, variant calling, and ultimately ROH inference. This will help us determine which approach is most appropriate while evaluating sensitivity to the published methodology.

### Script
fastp_merge_on.sbatch
fastp_merge_on_fastqc.sbatch
adapterremoval_pecnerova.sbatch
adapterremoval_post_fastqc.sbatch


### Results
QC and Sensitivity Comparison Tables
qc/fastp_merge/fastp_merge_summary.tsv
qc/adapterremoval_pecnerova/adapterremoval_summary.tsv
qc/preprocessing_merge_comparison.tsv

Metric	Current fastp	Fastp merge ON	AdapterRemoval
Sequencing units	54	54	54
Adapter trimming	Yes	Yes	Yes
Poly-G trimming	Yes	Yes	No explicit trimming
Overlapping-read merging	OFF	ON	ON (collapse)
Public samples: mean merge/collapse rate	—	13.91%	18.93%
Namibia: mean merge/collapse rate	—	43.36%	53.37%
Zambia: mean merge/collapse rate	—	23.21%	28.55%
R2 files with FastQC poly-G overrepresentation signals	Not yet consolidated	17/54	40/54


Note: AdapterRemoval consistently produced higher read-collapsing rates than fastp. This is consistent with differences in their overlap-detection criteria: AdapterRemoval permits shorter overlaps (11 bp by default) and uses a more permissive mismatch threshold (--mm 3), whereas fastp requires at least 30 bp of overlap and applies stricter mismatch limits by default. Therefore, the higher collapse rate should not be interpreted as evidence of superior preprocessing quality without examining downstream mapping and genotype-level metrics.

## Date: 11/2/2026

**Researcher:**   Xiaoyue
**Task:**  BWA-MEM alignment and PALEOMIX cleanup across three preprocessing branches

### Script
mapping_align_cleanup.sbatch
qc/mapping/master_mapping_manifest.tsv
qc/mapping/current_fastp_mapping.tsv
qc/mapping/fastp_merge_mapping.tsv
qc/mapping/adapterremoval_mapping.tsv

### Results
Complete and validate BWA-MEM alignment and PLAEOMIX cleanup across three preprocessing sensitivity branches (11/05)

Successfully extracted all 270 expected mapping records from 305 SLURM output logs.

- `current_fastp`: 54/54 BAMs
- `fastp_merge`: 108/108 BAMs
- `adapterremoval`: 108/108 BAMs

All BAMs passed `samtools quickcheck` and had corresponding index files.

Branch	Read type	Units	Primary mapped	Properly paired
Current fastp	Paired	54	98.7151%	94.7862%
Fastp merge ON	Paired	54	98.6253%	94.4504%
Fastp merge ON	Merged	54	98.4097%	N/A
AdapterRemoval	Paired	54	98.6756%	94.4728%
AdapterRemoval	Collapsed	54	98.4892%	N/A


Initial post-alignment QC has been completed for all 270 BAMs across the three preprocessing branches. Primary mapping rates were consistently high (98.41–98.72%), and properly paired rates ranged from 94.45% to 94.79% for paired-end alignments. 

These preliminary results suggest broadly comparable overall alignment performance among the preprocessing strategies. 

However, read-type-specific mapping rates are not directly equivalent because merging/collapsing changes the composition of the remaining paired-read subsets. Further sample-level MAPQ, mismatch, and coverage analyses are needed to evaluate alignment quality and potential reference-mapping bias.

**2. Singleton-rate analysis**
Run-to-run variation in singleton rates was observed across multiple public elephant samples.

Under `current_fastp`, the largest within-sample ranges were observed in DS1507 (0.574 percentage points) and DS1506 (0.516 percentage points).

Comparing `fastp_merge` and `adapterremoval`, AdapterRemoval showed a smaller within-sample singleton-rate range in 7 of the 9 samples with multiple sequencing units. However, the remaining paired-read populations differ between preprocessing branches, so this does not establish superior alignment quality.


### Next step
post-alignment processing, mapping quality assessment and preprocessing sensitivity evaluation


## Date: 11/05/2026 — 64-window pilot mapping QC

**Task:** Pilot mapping QC and soft-clipping assessment using reads sampled from 64 genomic windows.

### Script

`extract_mapping_flagstats.py`
`unit_mapping_flagstat.tsv`, 
`mapping_qc_pilot.py`
`mapping_qc_worker.py` 
`mapping_qc_pilot.sbatch`

### Results

| Branch | Read type | Sampled reads | MAPQ ≥20 | MAPQ ≥30 | Soft-clipped reads | Soft-clipped bases |
|---|---|---:|---:|---:|---:|---:|
| AdapterRemoval | Paired | 2,078,595 | 97.67% | 97.41% | 5.49% | 1.46% |
| Current fastp | Paired | 2,364,844 | 97.67% | 97.41% | 5.25% | 1.36% |
| fastp merge ON | Paired | 2,157,082 | 97.67% | 97.41% | 5.40% | 1.38% |
| AdapterRemoval | Collapsed | 145,223 | 97.14% | 96.58% | 13.46% | 6.80% |
| fastp merge ON | Merged | 104,557 | 97.12% | 96.58% | 3.74% | 0.88% |

**Interpretation:** This was a window-based sample, not a whole-BAM assessment. These values are preliminary and should not be used as the final mapping QC results.

## Date: 11/06/2026 onward — Full-BAM mapping QC and soft-clipping screening

**Task:** Assess mapping quality and clipping across the complete BAMs for all 54 sequencing units in each preprocessing branch.

### Script

- `mapping_qc_full.py`
- `mapping_qc_full_scan.c`
- `mapping_qc_full.sbatch`
- `mapping_qc_length_diagnostic.sbatch`

### Results: paired BAMs

| Branch | Units | MAPQ ≥20 | MAPQ ≥30 | MAPQ ≥40 | Soft-clipped reads | Soft-clipped bases |
|---|---:|---:|---:|---:|---:|---:|
| AdapterRemoval | 54 | 92.60% | 90.97% | 90.29% | 6.82% | 2.01% |
| Current fastp | 54 | 92.09% | 90.36% | 89.63% | 6.96% | 2.01% |
| fastp merge ON | 54 | 92.44% | 90.72% | 89.99% | 6.96% | 2.00% |

Across paired BAMs, the three branches had similar MAPQ distributions and soft-clipped-base percentages. The full-BAM results differ from the 64-window pilot: the pilot reported higher MAPQ ≥20 percentages (about 97.7%) and lower paired-read soft-clipping (about 5.3–5.5%). The full-BAM screening is the appropriate basis for final conclusions; the difference also indicates that the window sample did not represent the full-BAM distribution.

Soft-clipped reads varied among cohorts. In paired BAMs, public samples had approximately 6.98–7.40% soft-clipped reads, Namibia 5.20–5.54%, and focal Zambia 3.81–5.57% across branches.

### Results: collapsed and merged single-end BAMs

| Branch / read type | Units | MAPQ ≥20 | MAPQ ≥30 | Soft-clipped reads | Soft-clipped bases |
|---|---:|---:|---:|---:|---:|
| AdapterRemoval collapsed | 54 | 83.03% | 80.93% | 19.47% | 8.65% |
| fastp merge ON merged | 54 | 84.95% | 83.11% | 6.79% | 1.55% |

These read types show larger differences than the paired BAMs. Their results should be interpreted separately because AdapterRemoval collapsed reads and fastp merged reads are different preprocessing products. AdapterRemoval collapsed-read clipping was especially high in the public cohort (26.65% of reads; 12.04% of bases), compared with Namibia (5.84%; 1.69%) and focal Zambia (7.84%; 3.22%).

### Read-pair diagnostic

A targeted diagnostic was available for one public unit, KR0116 / ERR14017898, with 50 sampled reads per stratum. In the 250–269 bp low-clipping stratum, 94% of sampled reads were merged by fastp; in the matching high-clipping stratum, 8% were merged, while 56% mapped as proper pairs, 26% mapped to different references, and 10% had one or both reads unmapped. This is a one-unit, stratified diagnostic and should not be generalized to the full cohort.

`missing_flagstat_logs.tsv` contains only its header, indicating no missing flagstat entries in this summary.

### Output Files

- `qc/mapping/post_alignment/full_bam_qc/single/full_bam_qc_all.tsv`
- `qc/mapping/post_alignment/full_bam_qc/single/full_bam_qc_by_readtype.tsv`
- `qc/mapping/post_alignment/full_bam_qc/single/full_bam_qc_length_bins.tsv`
- `qc/mapping/post_alignment/length_diagnostic/`

### Next steps

- Use the full-BAM screening as the final basis for mapping QC; retain the 64-window results as a clearly labeled pilot.
- Review overall mapping-rate and depth/breadth summaries separately; those measures are not included in the attached tables.
- Expand the read-pair/clipping diagnostic beyond KR0116 before making cohort-wide claims about read-level outcomes.
- **Library-group note:** The Pečnerová replication track will follow their cleanup approach for comparability. Their handling may combine independently prepared libraries when marking duplicates. Keep the library-aware track separate and mark duplicates within verified physical libraries; resolve the Uganda sample library assignments before treating those groups as confirmed.
