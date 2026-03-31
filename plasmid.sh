#!/bin/bash
#SBATCH --job-name=plasmidQC
#SBATCH --mem=64G
#SBATCH -n 10
#SBATCH --time=120:00:00
#SBATCH --output=/gt/data/seqdma/plasmid_epi2me/slurmlog/wf-clone-validation_%j.log
#SBATCH -p gt_compute 

####################################################################
#Deliverables
#html report, 
#final.fasta file, 
####################################################################
# Make sure that a file named plasmid.txt is present in the run directory to be able to determine that it is a plasmid sequencing

# SampleSheet content looks like: 
#alias,barcode,approx_size,project_ID,delivery_folder,sample_type
#GT25-03693_H2H2,barcode01,6631,GTBH25-MurrayS-82,TechnologyEvaluationDevelopment_Group_BH,plasmid
#GT25-03694_T1H2,barcode02,6799,GTBH25-MurrayS-82,TechnologyEvaluationDevelopment_Group_BH,plasmid
#GT25-03695_T12H2,barcode03,8299,GTBH25-MurrayS-82,TechnologyEvaluationDevelopment_Group_BH,plasmid
#GT25-03696_T4H2,barcode04,8191,GTBH25-MurrayS-82,TechnologyEvaluationDevelopment_Group_BH,plasmid
#GT25-03697_H2,barcode05,6124,GTBH25-MurrayS-82,TechnologyEvaluationDevelopment_Group_BH,plasmid
#GT25-03698_Gli3,barcode06,11048,GTBH25-MurrayS-82,TechnologyEvaluationDevelopment_Group_BH,plasmid
####################################################################
####################################################################
####################################################################

###############################################
### Argument Parsing
###############################################

# Defaults
ONTdir_default="/gt/seqdata/Runs/ME_OxfordNanopore_GridIONX5_GXB03074"
QCdir_default="/gt/data/seqdma/plasmid_epi2me"

ONTdir="$ONTdir_default"
QCdir="$QCdir_default"

show_help () {
    echo ""
    echo "Usage: sbatch plasmidQC.sh [options]"
    echo ""
    echo "Options:"
    echo "  --ontdir PATH      Path to ONT run directory (RECOMMENDED: must change for each run)"
    echo "                     If not provided, script will use: $ONTdir_default"
    echo "                     WARNING: This default may be outdated."
    echo ""
    echo "  --qcdir PATH       Output QC directory (default: $QCdir_default)"
    echo ""
    echo "  -h, --help         Show this help message and exit"
    echo ""
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --ontdir)
            ONTdir="$2"
            shift 2
            ;;
        --qcdir)
            QCdir="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown argument: $1"
            show_help
            exit 1
            ;;
    esac
done

# Warn if ONTdir not provided
if [[ "$ONTdir" == "$ONTdir_default" ]]; then
    echo "WARNING: No --ontdir argument provided."
    echo "         Using default ONTdir:"
    echo "           $ONTdir_default"
    echo "         This path may NOT be correct for the new sequencing run!"
    echo "         Please pass --ontdir <path> to avoid data mismatch."
fi

# Create directories
mkdir -p "$QCdir"
mkdir -p "$QCdir/tmp"

# ---- Cleanup old Nextflow work directories (older than 14 days) ----
find "$QCdir" -type d -name "work" -mtime +14 -exec rm -rf {} +
### ── Environment Setup ───────────────────────────
module use --append /gt/research_development/qifa/elion/modulefiles
module load nextflow/23.04
module load singularity

# # Paths and defaults
# ONTdir="/gt/seqdata/Runs/ME_OxfordNanopore_GridIONX5_GXB03074/251205_GTBH25-MurrayS-146_GXB03074_005"
# QCdir="/gt/data/seqdma/plasmid_epi2me"
# mkdir -p "$QCdir" #create the directory if missing
# mkdir -p "$QCdir/tmp"

#export slumrlog="$QCdir/slurmlog/wf-clone-validation_%j.log"
export PERMISSION_DENIED_FILE="$QCdir/tmp/plasmid_qc_skip_dirs.txt"
touch "$PERMISSION_DENIED_FILE"

export NXF_SINGULARITY_CACHEDIR="$QCdir/epi2me/singularityCache"
export NXF_HOME="$QCdir/epi2me/nextflow"
export SINGULARITY_BIND="/gt/research_development"
export APPTAINER_BIND="/gt/research_development"
export SIGNAL="plasmid.txt"

senderDefault="GTdrylab@jax.org"

recipientCopy=(
  "Raman.Lawal@jax.org"
  "Harianto.Tjong@jax.org"
  "dave.john.harrison@jax.org"
)

recipients=(
  "Chrystal.Snow@jax.org"
  "Qingchang.Meng@jax.org"
)


#recipients="akinyanju.lawal@jax.org raman.lawal@jax.org"
#recipients="Chrystal.Snow@jax.org Qingchang.Meng@jax.org"

#delete logs older than 1 day
find /gt/data/seqdma/plasmid_epi2me/slurmlog/ \
    -type f -name "wf-clone-validation_*.log" \
    -mtime +0 -exec rm -f {} \;

### ── Main QC Functions ───────────────────────────
run_barcoded_qc_wf_clone_validation() {
  local log_file="$QCdir/$RunFolder/qc_run.log"
  local stderr_file="$QCdir/$RunFolder/qc_run.err"

  nextflow run epi2me-labs/wf-clone-validation \
    -r 'v1.6.0' \
    -c '/gt/research_development/epi2me/elion_config.cfg' \
    -profile 'elion2_singularity' \
    --sample_sheet "$SampleSheet" \
    --fastq "$ProjDir/fastq_pass" \
    --out_dir "$QCdir/$RunFolder" \
    -w "$QCdir/$RunFolder/work" \
    >"$log_file" 2>"$stderr_file"

  # Parse SampleSheet Header for Column Indices
  local header col_alias=-1 col_project_id=-1 col_delivery_folder=-1
  header=$(head -n1 "$SampleSheet" | tr -d '\r')
  IFS=',' read -r -a columns <<< "$header"

  for i in "${!columns[@]}"; do
    col_name=$(echo "${columns[$i]}" | tr '[:upper:]' '[:lower:]' | tr -d '\r')
    case "$col_name" in
      alias) col_alias=$i ;;
      project_id) col_project_id=$i ;;
      delivery_folder) col_delivery_folder=$i ;;
    esac
  done

  if [[ $col_alias -lt 0 ]]; then
    echo "[WARN] Alias column not found. Skipping file organization."
    return 0
  fi

  # Organize Files by project_ID 
  declare -A project_to_delivery

  tail -n +2 "$SampleSheet" | while IFS=',' read -r -a fields; do
    local alias project_ID delivery_folder

    alias=$(echo "${fields[$col_alias]}" | xargs)
    [[ -z "$alias" ]] && continue

    [[ $col_project_id -ge 0 ]] && project_ID=$(echo "${fields[$col_project_id]}" | xargs)
    [[ $col_delivery_folder -ge 0 ]] && delivery_folder=$(echo "${fields[$col_delivery_folder]}" | xargs)

    if [[ -n "$project_ID" ]]; then
      mkdir -p "$QCdir/$RunFolder/$project_ID"

      # Move matching alias files into the project_ID directory
      mapfile -t matched_files < <(find "$QCdir/$RunFolder" -type f -iname "*$alias*")
      for f in "${matched_files[@]}"; do
        mv -v "$f" "$QCdir/$RunFolder/$project_ID/"
      done

      # Track delivery destination
      project_to_delivery["$project_ID"]="$delivery_folder"
    fi
  done

  # Print rsync Suggestions ─
  header=$(head -n1 "$SampleSheet" | tr -d '\r')
  IFS=',' read -ra columns <<< "$header"

  col_project_id=-1
  col_delivery_folder=-1

  for i in "${!columns[@]}"; do
    normalized=$(echo "${columns[$i]}" | tr '[:upper:]' '[:lower:]' | tr -d '\r')
    case "$normalized" in
      project_id) col_project_id=$i ;;
      delivery_folder) col_delivery_folder=$i ;;
    esac
  done

  rsync_suggestions=""
  tail -n +2 "$SampleSheet" | while IFS=',' read -ra fields; do
    project_ID="${fields[$col_project_id]}"
    delivery_folder="${fields[$col_delivery_folder]}"
    rsync_suggestions+="rsync -vahP \"$QCdir/$RunFolder/$pid/\" \"/gt/gt_delivery/gt_secondary_analysis/$dpath/$pid/\"\n"
  done
}

run_barcoded_qc_wf-amplicon() {
  local log_file="$QCdir/$RunFolder/qc_run_amplicon.log"
  local stderr_file="$QCdir/$RunFolder/qc_run_amplicon.err"

  nextflow run epi2me-labs/wf-amplicon \
    -r 'v1.2.1' \
    -c '/gt/research_development/epi2me/elion_config.cfg' \
    -profile 'elion2_singularity' \
    --sample_sheet "$SampleSheet" \
    --fastq "$ProjDir/fastq_pass" \
    --out_dir "$QCdir/$RunFolder" \
    -w "$QCdir/$RunFolder/work" \
    >"$log_file" 2>"$stderr_file"

      # Parse SampleSheet Header for Column Indices
  local header col_alias=-1 col_project_id=-1 col_delivery_folder=-1
  header=$(head -n1 "$SampleSheet" | tr -d '\r')
  IFS=',' read -r -a columns <<< "$header"

  for i in "${!columns[@]}"; do
    col_name=$(echo "${columns[$i]}" | tr '[:upper:]' '[:lower:]' | tr -d '\r')
    case "$col_name" in
      alias) col_alias=$i ;;
      project_id) col_project_id=$i ;;
      delivery_folder) col_delivery_folder=$i ;;
    esac
  done

  if [[ $col_alias -lt 0 ]]; then
    echo "[WARN] Alias column not found. Skipping file organization."
    return 0
  fi

  # Organize Files by project_ID 
  declare -A project_to_delivery

  tail -n +2 "$SampleSheet" | while IFS=',' read -r -a fields; do
    local alias project_ID delivery_folder

    alias=$(echo "${fields[$col_alias]}" | xargs)
    [[ -z "$alias" ]] && continue

    [[ $col_project_id -ge 0 ]] && project_ID=$(echo "${fields[$col_project_id]}" | xargs)
    [[ $col_delivery_folder -ge 0 ]] && delivery_folder=$(echo "${fields[$col_delivery_folder]}" | xargs)

    if [[ -n "$project_ID" ]]; then
      mkdir -p "$QCdir/$RunFolder/$project_ID"

      # Move matching alias files into the project_ID directory
      mapfile -t matched_files < <(find "$QCdir/$RunFolder" -type f -iname "*$alias*")
      for f in "${matched_files[@]}"; do
        mv -v "$f" "$QCdir/$RunFolder/$project_ID/"
      done

      # Track delivery destination
      project_to_delivery["$project_ID"]="$delivery_folder"
    fi
  done

  # Print rsync Suggestions ─
  header=$(head -n1 "$SampleSheet" | tr -d '\r')
  IFS=',' read -ra columns <<< "$header"

  col_project_id=-1
  col_delivery_folder=-1

  for i in "${!columns[@]}"; do
    normalized=$(echo "${columns[$i]}" | tr '[:upper:]' '[:lower:]' | tr -d '\r')
    case "$normalized" in
      project_id) col_project_id=$i ;;
      delivery_folder) col_delivery_folder=$i ;;
    esac
  done

  rsync_suggestions=""
  tail -n +2 "$SampleSheet" | while IFS=',' read -ra fields; do
    project_ID="${fields[$col_project_id]}"
    delivery_folder="${fields[$col_delivery_folder]}"
    rsync_suggestions+="rsync -vahP \"$QCdir/$RunFolder/$pid/\" \"/gt/gt_delivery/gt_secondary_analysis/$dpath/$pid/\"\n"
  done
}


function run_unbarcoded_qc_wf_clone_validation() {
  SAMPLEID=$(grep "sample_id" "$checksummarysheet" | sed 's+=+\t+g' | cut -f2)
  nextflow run epi2me-labs/wf-clone-validation \
    -r 'v1.6.0' \
    -c '/gt/research_development/epi2me/elion_config.cfg' \
    -profile 'elion2_singularity' \
    --sample "$SAMPLEID" \
    --fastq "$ProjDir/fastq_pass" \
    -w "$QCdir/$RunFolder/work" \
    --out_dir "$QCdir/$RunFolder"
}

workflow_type() {
  local header sample_type_col=-1
  header=$(head -n 1 "$SampleSheet" | tr -d '\r')
  IFS=',' read -ra columns <<< "$header"

  # Find column index of sample_type (case-insensitive)
  for i in "${!columns[@]}"; do
    col_lower=$(echo "${columns[$i]}" | tr '[:upper:]' '[:lower:]' | xargs)
    if [[ "$col_lower" == "sample_type" ]]; then
      sample_type_col=$i
      break
    fi
  done

  if [[ $sample_type_col -lt 0 ]]; then
    echo "[WARN] sample_type column not found in $SampleSheet."
    # Build array of CC arguments
    cc_args=()
    for cc in "${recipientCopy[@]}"; do
        cc_args+=( -c "$cc" )
    done

    echo -e "SampleSheet: $SampleSheet\n\nMissing 'sample_type' column (case-insensitive).\nUnable to determine correct workflow." | \
      mailx -s "SampleSheet Missing sample_type Column" \
      -r "$senderDefault" \
      "${cc_args[@]}" \
      "${recipients[@]}" \
    return 1
  fi

  # Check values in sample_type column
  local found_bac_or_plasmid=0
  local found_amplicon=0

  while IFS=',' read -ra fields; do
    sample_type_val=$(echo "${fields[$sample_type_col]}" | tr -dc '[:print:]' | tr '[:upper:]' '[:lower:]' | xargs)
    case "$sample_type_val" in
      bac|plasmid) found_bac_or_plasmid=1 ;;
      amplicon) found_amplicon=1 ;;
    esac
  done < <(tail -n +2 "$SampleSheet")

  # Route workflow
  if [[ $found_bac_or_plasmid -eq 1 ]]; then
    run_barcoded_qc_wf_clone_validation && send_email_report_clone_validation
  elif [[ $found_amplicon -eq 1 ]]; then
    run_barcoded_qc_wf-amplicon && send_email_report_amplicon
  fi
}

send_email_report_clone_validation() {

    ############################################
    # Resolve GroupID and email subject text
    ############################################
    GroupID=$(grep "protocol_group_id" "$checksummarysheet" | sed 's+=+\t+g' | cut -f2)
    subject="COMPLETE: Plasmid QC for $GroupID"

    body_header="Find attached plasmid QC output files (HTML, sample sheet, and FASTA if small enough).

$rsync_suggestions

Note: This email is auto-generated. Contact gtdrylab@jax.org if you experience issues."

    mkdir -p "$QCdir/.mail"
    cd "$QCdir/$RunFolder" || return

    ############################################
    # CHECK IF RUN HAS COMPLETED
    # Must have ALL of these:
    #   feature_table.txt
    #   wf-clone-validation-report.html
    #   sample_status.txt
    ############################################
    if [[ -f "feature_table.txt" \
          && -f "wf-clone-validation-report.html" \
          && -f "sample_status.txt" ]]; then
        :
        # OK — proceed to email
    else
        # Not complete yet → exit silently (NO email)
        return 0
    fi

    ############################################
    # PREPARE FASTA FILES (compress and size-check)
    ############################################
    maxsize=$((20 * 1024 * 1024))   # 20 MB max per file
    attachments_fasta=()
    fasta_note=""

    # Recursively find all FASTA files
    fasta_files=$(find . -type f -name "*.fasta")
    if [[ -n "$fasta_files" ]]; then
        while IFS= read -r fasta_file; do
            compressed="${fasta_file}.gz"
            gzip -c "$fasta_file" > "$compressed"

            filesize=$(stat -c%s "$compressed")

            if [[ $filesize -le $maxsize ]]; then
                attachments_fasta+=("-a" "$compressed")
            else
                rm -f "$compressed"
                fasta_note+="\n⚠ FASTA file too large to attach: $fasta_file (>20 MB). Please send manually.\n"
            fi
        done <<< "$fasta_files"
    else
        fasta_note="\n⚠ No FASTA files were found for this run.\n"
    fi

    ############################################
    # BUILD ATTACHMENTS
    ############################################
    attachments=(
        -a "wf-clone-validation-report.html"
        -a "sample_status.txt"
    )

    # Add sample sheet if present
    if [[ -n "$SampleSheet" && -f "$SampleSheet" ]]; then
        attachments+=(-a "$SampleSheet")
    fi

    # Add compressed FASTA files
    if [[ ${#attachments_fasta[@]} -gt 0 ]]; then
        attachments+=("${attachments_fasta[@]}")
    fi

    ############################################
    # COMPOSE FINAL BODY TEXT
    ############################################
    final_body="${body_header}${fasta_note}"

    ############################################
    # SEND SINGLE FINAL EMAIL
    ############################################
    # Build array of CC arguments
    cc_args=()
    for cc in "${recipientCopy[@]}"; do
        cc_args+=( -c "$cc" )
    done

    # Send email
    echo -e "$final_body" | mailx -s "$subject" \
        -r "$senderDefault" \
        "${cc_args[@]}" \
        "${attachments[@]}" \
        "${recipients[@]}" \
        | tee "$QCdir/.mail/mailx_email_log.txt"
}

send_email_report_amplicon() {

    cd "$QCdir/$RunFolder" || return

    # ---- Completion check ----
    if [[ ! -f "wf-amplicon-report.html" \
          || ! -f "all-consensus-seqs.fasta" ]]; then
        return 0
    fi

    subject="COMPLETE: Amplicon QC for $RunFolder"

    body="Find attached amplicon QC results.

Attached:
- wf-amplicon-report.html
- consensus_fasta_bundle.tar.gz
- bam.tar.gz
- SampleSheet

Note: Large files may require manual transfer.

This email is auto-generated."

    maxsize=$((20 * 1024 * 1024))   # 20MB
    attachments=()

    # ----------------------------------------------------
    # Attach HTML
    # ----------------------------------------------------
    attachments+=(-a "wf-amplicon-report.html")

    # ----------------------------------------------------
    # Bundle FASTA + FAI
    # ----------------------------------------------------
    if [[ -f "all-consensus-seqs.fasta.fai" ]]; then
        tar -czf consensus_fasta_bundle.tar.gz \
            all-consensus-seqs.fasta \
            all-consensus-seqs.fasta.fai
    else
        tar -czf consensus_fasta_bundle.tar.gz \
            all-consensus-seqs.fasta
    fi

    fasta_size=$(stat -c%s consensus_fasta_bundle.tar.gz)

    if [[ $fasta_size -le $maxsize ]]; then
        attachments+=(-a "consensus_fasta_bundle.tar.gz")
    else
        body+="\n⚠ FASTA bundle too large to attach. Please transfer manually."
        rm -f consensus_fasta_bundle.tar.gz
    fi

    # ----------------------------------------------------
    # Extract project_ID from SampleSheet (FIRST DATA ROW)
    # ----------------------------------------------------
    project_id=$(awk -F',' 'NR==2 {gsub(/\r/,""); print $4; exit}' "$SampleSheet")

    if [[ -n "$project_id" && -d "$project_id" ]]; then

        tar -czf bam.tar.gz "$project_id"
        proj_size=$(stat -c%s bam.tar.gz)

        if [[ $proj_size -le $maxsize ]]; then
            attachments+=(-a "bam.tar.gz")
        else
            body+="\n⚠ Project folder too large to attach. Please transfer manually."
            rm -f bam.tar.gz
        fi
    else
        body+="\n⚠ Project directory not found: $project_id"
    fi

    # ----------------------------------------------------
    # Attach SampleSheet
    # ----------------------------------------------------
    if [[ -n "$SampleSheet" && -f "$SampleSheet" ]]; then
        attachments+=(-a "$SampleSheet")
    fi

    # ----------------------------------------------------
    # CC handling
    # ----------------------------------------------------
    cc_args=()
    for cc in "${recipientCopy[@]}"; do
        cc_args+=( -c "$cc" )
    done

    # ----------------------------------------------------
    # Send Email
    # ----------------------------------------------------
    echo -e "$body" | mailx -s "$subject" \
        -r "$senderDefault" \
        "${cc_args[@]}" \
        "${attachments[@]}" \
        "${recipients[@]}"

    # ----------------------------------------------------
    # Cleanup temp files
    # ----------------------------------------------------
    rm -f consensus_fasta_bundle.tar.gz
    rm -f bam.tar.gz
}
#Auto-detect the True Run Folder
find_true_runfolder() {
    local start="$1"

    # Search downwards first for fastq_pass
    candidate=$(find "$start" -type d -name "fastq_pass" -exec dirname {} \; | head -n 1)
    if [[ -n "$candidate" ]]; then
        echo "$candidate"
        return 0
    fi

    # Fallback: search upwards if no fastq_pass found in subdirs
    while [[ "$start" != "/" ]]; do
        if [[ -d "$start/fastq_pass" ]]; then
            echo "$start"
            return 0
        fi
        start=$(dirname "$start")
    done

    # Nothing found
    echo ""
    return 1
}


### ── Per-Project Processing ───────────────────────
function process_plasmid_qc_projects() {
  local checkrplasmid
    checkrplasmid=$(find "$ONTdir" -type f -iname "$SIGNAL")

  if [[ -z "$checkrplasmid" ]]; then
    echo "No signal file ($SIGNAL) found in $ONTdir. Nothing to process."
    return 0
  fi

  mapfile -t project_list <<< "$checkrplasmid"
  for plasmid_file in "${project_list[@]}"; do

    #ProjDir="${plasmid_file%/*}"
    #RunFolder=$(basename "$ProjDir")
    # ----------------------------------------
    # SKIP YEAR ARCHIVE FOLDERS
    # ----------------------------------------
    if [[ "$(basename "$raw_plasmid_dir")" =~ ^20[0-4][0-9]$ ]]; then
        echo "[INFO] Skipping archive folder: $(basename "$raw_plasmid_dir")"
        continue
    fi
    raw_plasmid_dir="${plasmid_file%/*}"
    ProjDir=$(find_true_runfolder "$raw_plasmid_dir")
    if [[ -z "$ProjDir" ]]; then
      echo "[ERROR] Could not determine run folder for $raw_plasmid_dir"
      continue
    fi

    RunFolder=$(basename "$ProjDir")
    
    #SampleSheet=$(find "$ProjDir" -type f -iregex '.*/sample[_-]*sheet\.csv' | head -n 1)
    # 1. Look specifically for plasmid_samplesheet.csv
    # ───────────────────────────────────────────────────────────────
    # Smart SampleSheet detection with strict header validation
    # Priority 1 → plasmid_samplesheet.csv
    # Priority 2 → any *sample*sheet*.csv that contains required headers
    # ───────────────────────────────────────────────────────────────
        required_headers=("alias" "barcode" "approx_size" "project_id" "delivery_folder" "sample_type")
        SampleSheet=""
        search_dir="$ProjDir"

        while true; do
            # Find all CSVs in current directory
            candidates=($(find "$search_dir" -maxdepth 1 -type f -iname "*sample*sheet*.csv"))

            for file in "${candidates[@]}"; do
                header=$(head -n 1 "$file" | tr -d '\r' | tr '[:upper:]' '[:lower:]')
                valid=1

                for h in "${required_headers[@]}"; do
                    # Convert required header to regex accepting underscore or space
                    regex_header=$(echo "$h" | sed 's/_/[_ ]/g')
                    if ! grep -qiE "\b$regex_header\b" <<< "$header"; then
                        valid=0
                        echo "[DEBUG] Missing header: $h in $file"
                        break
                    fi
                done

                if [[ $valid -eq 1 ]]; then
                    SampleSheet="$file"
                    break 2  # Exit both loops
                fi
            done

            # Stop if we reached $ONTdir or root
            if [[ "$search_dir" == "$ONTdir" ]] || [[ "$search_dir" == "/" ]]; then
                break
            fi

            search_dir=$(dirname "$search_dir")
        done

        echo "[DEBUG] SampleSheet final: $SampleSheet"

        # Finally, check $ONTdir itself if nothing found
        if [[ -z "$SampleSheet" ]]; then
            candidate=$(find "$ONTdir" -maxdepth 1 -type f -iname "plasmid_samplesheet.csv" | head -n 1)
            if [[ -n "$candidate" ]]; then
                header=$(head -n 1 "$candidate" | tr -d '\r' | tr '[:upper:]' '[:lower:]')
                valid=1
                for h in "${required_headers[@]}"; do
                    regex_header=$(echo "$h" | sed 's/_/[_ ]/g')
                    if ! grep -qiE "\b$regex_header\b" <<< "$header"; then
                        valid=0
                        break
                    fi
                done
                [[ $valid -eq 1 ]] && SampleSheet="$candidate"
            fi
        fi
    # ---- Final check: missing or invalid SampleSheet ----
    if [[ -z "$SampleSheet" ]]; then
        echo "[ERROR] No valid SampleSheet detected for $RunFolder. Skipping."

        message="No valid SampleSheet found in:
$ProjDir

Files checked:
- plasmid_samplesheet.csv
- all *sample*sheet*.csv files

None contained ALL required headers:
alias, barcode, approx_size, project_id, delivery_folder, sample_type

Run skipped: $RunFolder"

    # Build array of CC arguments
    cc_args=()
    for cc in "${recipientCopy[@]}"; do
        cc_args+=( -c "$cc" )
    done
        echo "$message" | mailx \
          -s "Missing or Invalid SampleSheet: $RunFolder" \
          -r "$senderDefault" \
          "${cc_args[@]}" \
          "${recipients[@]}" \

        rm -r "$QCdir/$RunFolder"
        continue
    fi

    export SampleSheet
    echo "[INFO] Selected SampleSheet: $SampleSheet"

	
	checkreport=$(find $ProjDir -type f -name "report*")
	checksummarysheet=$(find $ProjDir -type f -name "final_summary*")
	checkporeactivity=$(find $ProjDir -type f -name "pore_activity*")
	checkthrougput=$(find $ProjDir -type f -name "throughput*.csv")
	#Only proceed if all four of these variables are non-empty (i.e., all required output files were found).
    if [[ -z "$checkreport" || -z "$checksummarysheet" || -z "$checkporeactivity" || -z "$checkthrougput" ]]; then
      echo "[WARN] Required files missing or empty in $ProjDir. Skipping QC..."
      #echo -e "Required file(s) missing or empty in:\n$ProjDir\n\nPlease ensure all of the following exist and are not empty:\n- report*\n- final_summary*\n- pore_activity*\n- throughput*.csv\n\nThis run will be skipped." | \
      #  mailx -s "Skipped: Missing Key Files for $RunFolder" -r "$senderDefault" -c "$recipientCopy" $recipients
      continue
    fi

    # if [[ -d "$QCdir/$RunFolder" ]]; then
    #   if [[ -f "$QCdir/$RunFolder/feature_table.txt" && -f "$QCdir/$RunFolder/wf-clone-validation-report.html" && -f "$QCdir/$RunFolder/sample_status.txt" ]]; then
    #     echo "⚠️  QC previously completed for $RunFolder. Skipping..."
		#     echo "To re-run, delete the folder"
    #     continue
    #   fi
    # fi

    if [[ -d "$QCdir/$RunFolder" ]]; then

      # ---- Clone validation already done ----
      if [[ -f "$QCdir/$RunFolder/feature_table.txt" \
        && -f "$QCdir/$RunFolder/wf-clone-validation-report.html" \
        && -f "$QCdir/$RunFolder/sample_status.txt" ]]; then
          echo "⚠️  Clone-validation already completed for $RunFolder. Skipping..."
          continue
      fi
      # ---- Amplicon already done ----
      if [[ -f "$QCdir/$RunFolder/wf-amplicon-report.html" \
        && -f "$QCdir/$RunFolder/all-consensus-seqs.fasta" ]]; then
          echo "⚠️  Amplicon QC already completed for $RunFolder. Skipping..."
          continue
      fi
    fi

    mkdir -p "$QCdir/$RunFolder"

    if [[ -d "$ProjDir/fastq_pass" ]]; then
      if ls "$ProjDir/fastq_pass" | grep -q "^barcode"; then
        if [[ -n "$SampleSheet" ]]; then
            # ── Header validation block begins here ──
            required_headers=("alias" "barcode" "approx_size" "project_id" "delivery_folder" "sample_type")
            # Read and normalize header
            header_line=$(head -n 1 "$SampleSheet" | tr -d '\r' | tr '[:upper:]' '[:lower:]')
            missing_headers=()

            for h in "${required_headers[@]}"; do
                if ! grep -qw "$h" <<< "$header_line"; then
                    missing_headers+=("$h")
                fi
            done

            if [[ ${#missing_headers[@]} -gt 0 ]]; then
                echo "[ERROR] SampleSheet missing required headers: ${missing_headers[*]}"
                msg="SampleSheet in $ProjDir is missing required column(s): ${missing_headers[*]}.

Required headers:
- alias
- barcode
- approx_size
- project_id
- delivery_folder
- sample_type

Please correct the file and re-run the QC pipeline.

SampleSheet path: $SampleSheet

RunFolder: $RunFolder

Note: This email is auto-generated. Contact gtdrylab@jax.org for help."
            # Build array of CC arguments
            cc_args=()
            for cc in "${recipientCopy[@]}"; do
              cc_args+=( -c "$cc" )
            done
                echo "$msg" | mailx -s "Invalid SampleSheet Headers: $RunFolder" \
                  -r "$senderDefault" \
                  "${cc_args[@]}" \
                  "${recipients[@]}" \
                  exit 1
            fi

          export SampleSheet
          workflow_type
        else
            # Build array of CC arguments
            cc_args=()
            for cc in "${recipientCopy[@]}"; do
              cc_args+=( -c "$cc" )
            done
            echo -e "SampleSheet missing in $ProjDir. QC not run." | mailx -s "Missing SampleSheet in $RunFolder" \
              -r "$senderDefault" \
              "${cc_args[@]}" \
              "${recipients[@]}" \
          rm -r "$QCdir/$RunFolder"
        fi
      else
	  	run_unbarcoded_qc_wf_clone_validation && send_email_report_clone_validation
		echo "unbarcoded"
      fi
    else
      optional_files=(
        "$(find "$ProjDir" -name "report*" -print -quit)"
        "$(find "$ProjDir" -name "final_summary*" -print -quit)"
        "$(find "$ProjDir" -name "pore_activity*" -print -quit)"
        "$(find "$ProjDir" -name "throughput*.csv" -print -quit)"
      )

      if [[ ${#optional_files[@]} -ge 2 ]]; then
        echo -e "Missing fastq_pass directory for project $ProjDir.\nPlease investigate." | \
          mailx -s "⚠️  Missing fastq_pass @ $ProjDir" -r "$senderDefault" -c "$recipientCopy" $recipientCopy
      fi
    fi
  done
}

### ── Pipeline Trigger ─────────────────────────────
process_plasmid_qc_projects
