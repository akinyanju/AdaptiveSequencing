# Plasmid QC Pipeline

## Overview
This pipeline performs automated quality control (QC) for Oxford Nanopore plasmid sequencing runs using Nextflow workflows from EPI2ME.

## Deliverables
- HTML QC report  
- Final FASTA file(s)  
- Organized output by project  

---

## Requirements

### Environment
- SLURM scheduler
- nextflow/23.04
- singularity (or Apptainer)

### Directory Access
- /gt/research_development  
- /gt/data/seqdma  
- /gt/seqdata  

---

## Input Requirements

### 1. ONT Run Directory (--ontdir)
Must contain:
- fastq_pass/ directory  
- report*  
- final_summary*  
- pore_activity*  
- throughput*.csv  

---

### 2. Signal File (Required)
plasmid.txt must exist somewhere in the ONT directory to trigger processing.

---

### 3. Sample Sheet (CSV)

Required columns:
- alias  
- barcode  
- approx_size  
- project_id  
- delivery_folder  
- sample_type  

Example:
alias,barcode,approx_size,project_ID,delivery_folder,sample_type  
Sample1,barcode01,6631,ProjectA,DeliveryFolderA,plasmid  

---

## How to Run

sbatch plasmidQC.sh --ontdir <ONT_RUN_PATH> --qcdir <OUTPUT_PATH>

---

## Parameters

--ontdir   Path to ONT run directory (REQUIRED)  
--qcdir    Output directory (default: /gt/data/seqdma/plasmid_epi2me)  
-h         Show help  

---

## Workflow Logic

1. Detects plasmid runs using plasmid.txt  
2. Locates correct run folder (auto-detection)  
3. Identifies valid SampleSheet with required headers  
4. Determines workflow type:
   - plasmid/bac → wf-clone-validation  
   - amplicon → wf-amplicon  
5. Executes Nextflow pipeline  
6. Organizes output by project_id  
7. Sends email report with attachments  

---

## Output Structure

<QCdir>/<RunFolder>/  
├── wf-clone-validation-report.html  
├── sample_status.txt  
├── feature_table.txt  
├── *.fasta  
├── work/  
└── <project_id>/  

---

## Email Notifications

Automatically sends completion emails with:
- HTML report  
- SampleSheet  
- FASTA files (if <20MB)  

If files exceed size limit, user is notified to transfer manually.

---

## Auto Cleanup

- Removes Nextflow work directories older than 14 days  
- Deletes SLURM logs older than 1 day  

---

## Common Issues

Missing plasmid.txt  
→ Run will be skipped  

Missing SampleSheet  
→ Email notification sent, run skipped  

Invalid SampleSheet headers  
→ Error + email notification  

Missing fastq_pass  
→ Warning email sent  

---

## Notes

- Pipeline skips already completed runs  
- SampleSheet is strictly validated  
- Supports both barcoded and unbarcoded runs  
- Fully automated email reporting  
