#!/bin/bash

####################################################################################################
#          This script is used to run diamond jobs, from a list of species                         #
#                                                                                                  #
#                                     Command example:                                             #
#          screen -dmS diamond_run -L -Logfile ./logs/Diamond/RunDiamond20250402.log bash          #
#       ./scripts/program_runners/diamond.sh ./files/FANTASIA_project/inputs/diamond_input.csv     #
####################################################################################################

# Check if a file argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <species_list_file>"
    exit 1
fi

# The file containing the list of species
SPECIES_FILE=$1

source /data/users/demartini/miniconda3/etc/profile.d/conda.sh
# Loop through each line in the file
while IFS=$',' read -r SPECIES FASTA GFF; do

    echo "Processing: $SPECIES"
    
    WORKDIR="/data/users/demartini/FANTASIA_project/species/$SPECIES/04_FunctionalAnnotation/Diamond/"
    # mkdir -p "$WORKDIR"
    cd "$WORKDIR"
    
    # # Extract the prefix from the species
    PREFIX=$(echo "$SPECIES" | awk -F'_' '{print substr($1,1,2) substr($2,1,3)}')

    # # Search for the FASTA file with the longest isoform of each protein
    INPUT=$(find ../../ -type f -name '*_longest_isof_proteins.fasta') 

    # Check if the input file exists and is not empty
    if [[ ! -s "$INPUT" ]]; then
        echo "LONGEST_ISOF_FASTA file not found. Creating it..."
        
        if [[ -s "$GFF" && -s "$FASTA" ]]; then
            # Obtain only the longest isoform per gene model using AGAT
            agat_sp_keep_longest_isoform.pl -gff "$GFF" -o "${GFF%.gff*}_longest_isof.gff"
            
            # Obtain the protein sequence with gffread
            gffread -y "${GFF%.gff*}_longest_isof_proteins.fasta" -g "$FASTA" "${GFF%.gff*}_longest_isof.gff"

            # Reassign INPUT after creating it
            INPUT="${GFF%.gff*}_longest_isof_proteins.fasta"
        else
            echo "Error: GFF file or genomic FASTA file not found. Cannot generate FASTA."
            continue
        fi
    fi
    # Remove dots from sequences (but not from headers)
    echo "Removing dots from sequences in $INPUT..."
    awk '/^>/ { print; next }{ gsub(/\./, "", $0); print }' "$INPUT" > tmpfile && mv tmpfile "$INPUT"
    # 1. Run DIAMOND with Swiss-Prot database
    if [ ! -s "${PREFIX}_sprot.tsv" ]; then
        echo "Running DIAMOND with Swiss-Prot..."
        diamond blastp -d /data/shared_dbs/swissprot/uniprot_sprot_r2025_01.dmnd -q $INPUT -o "${PREFIX}_sprot.tsv" --outfmt 6 --evalue 1e-5 --max-target-seqs 1 --threads 8 > "${PREFIX}_sprot.log" 2>&1
    else
        echo "File ${PREFIX}_sprot.tsv already exists, skipping DIAMOND with Swiss-Prot."
    fi

    # 2. Run DIAMOND with TrEMBL database
    if [ ! -s "${PREFIX}_trembl.tsv" ]; then
        echo "Running DIAMOND with TrEMBL..."
        diamond blastp -d /data/shared_dbs/swissprot/uniprot_trembl_r2025_01.dmnd -q $INPUT -o "${PREFIX}_trembl.tsv" --outfmt 6 --evalue 1e-5 --max-target-seqs 1 --threads 8 > "${PREFIX}_trembl.log" 2>&1
    else
        echo "File ${PREFIX}_trembl.tsv already exists, skipping DIAMOND with TrEMBL."
    fi

    # 3. Extract the Swiss-Prot and TrEMBL IDs from the DIAMOND results 
    echo "Extracting Swiss-Prot and TrEMBL IDs from DIAMOND outputs..."
    cut -f 1,2 "${PREFIX}_sprot.tsv" > "${PREFIX}_sp_ids.txt"
    cut -f 1,2 "${PREFIX}_trembl.tsv" > "${PREFIX}_tr_ids.txt"

    awk -F'\t' '{split($2,a,"|"); print $1"\t"a[2]}' "${PREFIX}_sp_ids.txt" > "${PREFIX}_sprot_ids.txt"
    awk -F'\t' '{split($2,a,"|"); print $1"\t"a[2]}' "${PREFIX}_tr_ids.txt" > "${PREFIX}_trembl_ids.txt"

    rm "${PREFIX}_sp_ids.txt" "${PREFIX}_tr_ids.txt"

    python3 /data/users/demartini/FANTASIA_project/scripts/uniprot_GO_mapper.py "${PREFIX}_sprot_ids.txt" "${PREFIX}_trembl_ids.txt" "$PREFIX"
    echo "GO terms extraction completed."

    # 4. Use AHRD to extract GOs from uniprot terms

    # Paths
    AHRD_JAR="/data/software/AHRD/dist/ahrd.jar"
    BLACKLIST="/data/software/AHRD/test/resources/blacklist_descline.txt"
    FILTER_SPROT="/data/software/AHRD/test/resources/filter_descline_sprot.txt"
    FILTER_TREMBL="/data/software/AHRD/test/resources/filter_descline_trembl.txt"
    TOKEN_BLACKLIST="/data/software/AHRD/test/resources/blacklist_token.txt"
    GO_GAF="/data/shared_dbs/swissprot/goa_uniprot_all.gaf"
    UNIPROT_SPROT="/data/shared_dbs/swissprot/uniprot_sprot_r2025_01.fasta"
    UNIPROT_TREMBL="/data/shared_dbs/swissprot/uniprot_trembl_r2025_01.fasta"
    OUT_FILE="${WORKDIR}/${PREFIX}.proteins.funct_ahrd.tsv"
    SPROT="${WORKDIR}/${PREFIX}_sprot.tsv"
    TREMBL="${WORKDIR}/${PREFIX}_trembl.tsv"

    # YAML temporal
    YAML_TMP="tmp_ahrd_conf.yml"

    # Check if required files exist
    if [[ ! -f "$INPUT" || ! -f "$SPROT" || ! -f "$TREMBL" ]]; then
        echo "Missing files for $SPECIES skipping."
        continue
    fi

    cat > "$YAML_TMP" <<EOF
proteins_fasta: $INPUT
token_score_bit_score_weight: 0.468
token_score_database_score_weight: 0.2098
token_score_overlap_score_weight: 0.3221
gene_ontology_result: $GO_GAF
reference_go_regex: '^UniProtKB\\t(?<shortAccession>[^\\t]+)\\t[^\\t]+\\t(?!NOT\\|)[^\\t]*\\t(?<goTerm>GO:\\d{7})'
prefer_reference_with_go_annos: true
output: $OUT_FILE
blast_dbs:
  swissprot:
    weight: 653
    description_score_bit_score_weight: 2.717061
    file: $SPROT
    database: $UNIPROT_SPROT
    blacklist: $BLACKLIST
    filter: $FILTER_SPROT
    token_blacklist: $TOKEN_BLACKLIST

  trembl:
    weight: 904
    description_score_bit_score_weight: 2.590211
    file: $TREMBL
    database: $UNIPROT_TREMBL
    blacklist: $BLACKLIST
    filter: $FILTER_TREMBL
    token_blacklist: $TOKEN_BLACKLIST
EOF

    # Run AHRD
    echo "Running AHRD for $SPECIES..."
    java -Xmx2g -jar "$AHRD_JAR" "$YAML_TMP"

done < "$SPECIES_FILE"


# Update file of files (remove species already annotated)
# cd /data/users/demartini/FANTASIA_project/species
# bash ../scripts/program_runners/check_annotation_status.sh
# python3 /data/users/demartini/FANTASIA_project/scripts/program_runners/file_of_files_manager.py '/data/users/demartini/FANTASIA_project/files/FANTASIA_project/inputs/diamond_input.csv'