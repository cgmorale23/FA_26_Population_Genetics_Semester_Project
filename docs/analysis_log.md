
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

## Date: 

**Researcher:**   
**Task:** 

### Objective

### Script

### Commands executed
```
```
### SLURM job information

### Results

### Methodological decisions

### Next steps

