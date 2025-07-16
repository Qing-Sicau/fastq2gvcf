#!/bin/bash
#================================================================================================
#
#  VCF Quality Control Pipeline with BCFtools (Resumable & Robust Manual Stats)
#
#  Description: This script performs strict QC on a VCF file. It is designed to be
#               resumable and uses a robust manual calculation for average depth
#               to avoid issues with 'bcftools stats'.
#
#================================================================================================

#-----------------------------------------------------------------------------------------------
# 1. Variable Definitions
#-----------------------------------------------------------------------------------------------

# --- Input and Output File Definitions ---
RAW_VCF="raw.vcf.gz"
STEP1_OUT="01_hard_filtered.vcf.gz"
STEP2_OUT="02_genotype_filtered.vcf.gz"
FINAL_VCF="quality_controlled.vcf.gz"

# --- Resource Allocation ---
N_THREADS=120

# --- Filtering Thresholds ---
N_SAMPLES=400
MIN_DEPTH_PER_SAMPLE=20
MIN_GQ=30
MAX_MISSING_RATE=0.1
MIN_MAF=0.01

# --- GATK Hard-Filtering Thresholds ---
QD_THRESH=5.0
FS_THRESH=20.0
MQ_THRESH=50.0
SOR_THRESH=2.0
MQRankSum_THRESH=-5.0
ReadPosRankSum_THRESH=-3.0


#MIN_GQ: 30 (or 40 for ultra-strict scenarios)
#MAX_MISSING_RATE: 0.1 (or 0.05 for ultra-strict)
#MIN_MAF: 0.01 (or 0.05 for population studies)
#QD_THRESH: 5.0 (or 8.0 for ultra-strict)
#FS_THRESH: 20.0 (SNPs), 100.0 (indels)
#MQ_THRESH: 50.0 (or 60.0 for ultra-strict)
#SOR_THRESH: 2.0 (SNPs), 4.0 (indels)
#MQRankSum_THRESH: -5.0 (or -3.0 for ultra-strict)
#ReadPosRankSum_THRESH: -3.0 (or -2.0 for ultra-strict)


#-----------------------------------------------------------------------------------------------
# 2. Quality Control Execution
#-----------------------------------------------------------------------------------------------

echo "Starting VCF quality control pipeline..."
echo "Script will check for existing files to resume from the last completed step."

# --- Step 2.1: Initial Site-Level Filtering ---
if [ ! -f "${STEP1_OUT}" ]; then
    echo "--- Step 1/3: Selecting biallelic SNPs and applying GATK hard filters ---"
    bcftools view \
      --threads ${N_THREADS} \
      -m2 -M2 -v snps \
      ${RAW_VCF} | \
    bcftools filter \
      --threads ${N_THREADS} \
      --exclude "QD < ${QD_THRESH} || FS > ${FS_THRESH} || SOR > ${SOR_THRESH} || MQ < ${MQ_THRESH} || MQRankSum < ${MQRankSum_THRESH} || ReadPosRankSum < ${ReadPosRankSum_THRESH}" \
      -o ${STEP1_OUT} -Oz

    echo "Indexing ${STEP1_OUT}..."
    bcftools index --threads ${N_THREADS} -f ${STEP1_OUT}
    echo "--- Step 1/3 completed. ---"
else
    echo "--- [SKIP] Step 1: Output file '${STEP1_OUT}' already exists. ---"
fi


# --- Step 2.2: Genotype-Level Filtering ---
if [ ! -f "${STEP2_OUT}" ]; then
    echo "--- Step 2/3: Filtering individual genotypes based on depth (DP) and quality (GQ)... ---"
    bcftools filter \
      --threads ${N_THREADS} \
      -S . \
      -e "FMT/DP < ${MIN_DEPTH_PER_SAMPLE} || FMT/GQ < ${MIN_GQ}" \
      -o ${STEP2_OUT} -Oz \
      ${STEP1_OUT}

    echo "Indexing ${STEP2_OUT}..."
    bcftools index --threads ${N_THREADS} -f ${STEP2_OUT}
    echo "--- Step 2/3 completed. ---"
else
    echo "--- [SKIP] Step 2: Output file '${STEP2_OUT}' already exists. ---"
fi


# --- Step 2.3: Final Variant-Level Filtering ---
if [ ! -f "${FINAL_VCF}" ]; then
    echo "--- Step 3/3: Applying final filters on missing rate, MAF, and optionally, depth... ---"
    echo "Manually calculating average depth from FORMAT/DP fields..."

    AVG_DEPTH_PER_SAMPLE=$(bcftools query -f '[%DP\t]\n' ${STEP2_OUT} | awk '
        BEGIN { total_dp = 0; count = 0; }
        {
            for (i=1; i<=NF; i++) {
                if ($i != "." && $i != "0") { # Exclude missing and zero-depth values
                    total_dp += $i;
                    count++;
                }
            }
        }
        END {
            if (count > 0) {
                print total_dp / count;
            } else {
                print 0; # Avoid division by zero if file is empty or all depth is missing
            }
        }
    ')

    # Define the base filter expression that will always be applied.
    FILTER_EXPRESSION="F_MISSING > ${MAX_MISSING_RATE} || MAF[0] < ${MIN_MAF}"

    # Check if a valid, non-zero average depth was calculated.
    if [[ $(echo "$AVG_DEPTH_PER_SAMPLE > 0" | bc) -eq 1 ]]; then
        # If successful, calculate the max depth threshold and add it to the filter expression.
        MAX_TOTAL_DP=$(echo "($AVG_DEPTH_PER_SAMPLE * 3 * $N_SAMPLES)/1" | bc)
        echo "Manual calculation successful. Average depth per sample is ${AVG_DEPTH_PER_SAMPLE}."
        echo "Adding high-depth filter (INFO/DP > ${MAX_TOTAL_DP})."
        FILTER_EXPRESSION="${FILTER_EXPRESSION} || INFO/DP > ${MAX_TOTAL_DP}"
    else
        # If it failed, print a warning and proceed without the depth filter.
        echo "Warning: Manual calculation of average depth resulted in zero or an error."
        echo "Proceeding without the high-depth filter."
    fi

    echo "Applying final filter expression: ${FILTER_EXPRESSION}"
    bcftools filter \
      --threads ${N_THREADS} \
      -e "${FILTER_EXPRESSION}" \
      -o ${FINAL_VCF} -Oz \
      ${STEP2_OUT}

    echo "Indexing final file ${FINAL_VCF}..."
    bcftools index --threads ${N_THREADS} -f ${FINAL_VCF}
    echo "--- Step 3/3 completed. ---"
else
    echo "--- [SKIP] Step 3: Final output file '${FINAL_VCF}' already exists. ---"
fi

echo ""
echo "VCF quality control pipeline has finished. The final high-quality file is: ${FINAL_VCF}"
