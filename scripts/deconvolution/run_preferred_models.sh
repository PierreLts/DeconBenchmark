#!/bin/bash
# run_preferred_models.sh
# This script runs a predefined list of recommended deconvolution models
#
# Usage: ./run_preferred_models.sh [INPUT_RDA_FILE]

# Default input file
DEFAULT_INPUT="/work/gr-fe/lorthiois/project2/data/Batch1_Downsampled.rda"
INPUT_RDA_FILE=${1:-$DEFAULT_INPUT}

# Reference-based methods with highest scores according to benchmark
RECOMMENDED_MODELS="AdRoit,ARIC,AutoGeneS,BayCount,BayesCCE,BayesPrism,BayICE,BisqueMarker,BisqueRef,BseqSC,CellDistinguisher,CIBERSORT,CPM,DAISM,debCAM,Deblender,DeCompress,deconf,DeconICA,DeconPeaker,DeconRNASeq,deconvSeq,DecOT,DeMixT,DESeq2,digitalDLSorter,DSA,dtangle,DWLS,EMeth,EPIC,FARDEEP,ImmuCellAI,LinDeconSeq,Linseed,MCPcounter,MethylResolver,MIXTURE,MOMF,MuSic,MySort,NITUMID,PREDE,quanTIseq,ReFACTor,RNA-Sieve,scaden,SCDC,spatialDWLS,TOAST"

# Run the main script
echo "Running recommended models: $RECOMMENDED_MODELS"
echo "Using input file: $INPUT_RDA_FILE"

./run_multiple_models.sh "$INPUT_RDA_FILE" "$RECOMMENDED_MODELS"