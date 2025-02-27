#!/bin/bash
# run_top_methods.sh - Run the top deconvolution methods in parallel
#
# Usage: ./run_top_methods.sh [INPUT_RDA_FILE]

# Default input file
DEFAULT_INPUT="/work/gr-fe/lorthiois/DeconBenchmark/data/Batch1_Downsampled.rda"
INPUT_RDA_FILE=${1:-$DEFAULT_INPUT}

# Top recommended methods based on the review paper
TOP_METHODS="AdRoit,ARIC,AutoGeneS,BayCount,BayesCCE,BayesPrism,BayICE,BisqueMarker,BisqueRef,BseqSC,CellDistinguisher,CIBERSORT,CPM,DAISM,debCAM,Deblender,DeCompress,deconf,DeconICA,DeconPeaker,DeconRNASeq,deconvSeq,DecOT,DeMixT,DESeq2,digitalDLSorter,DSA,dtangle,DWLS,EMeth,EPIC,FARDEEP,ImmuCellAI,LinDeconSeq,Linseed,MCPcounter,MethylResolver,MIXTURE,MOMF,MuSic,MySort,NITUMID,PREDE,quanTIseq,ReFACTor,RNA-Sieve,scaden,SCDC,spatialDWLS,TOAST"

# Run the parallel submission script
./parallel_deconv.sh "$INPUT_RDA_FILE" "$TOP_METHODS"

echo "All jobs submitted. Run './evaluate_results.sh' after completion to process results."

# "AdRoit,ARIC,AutoGeneS,BayCount,BayesCCE,BayesPrism,BayICE,BisqueMarker,BisqueRef,BseqSC,CellDistinguisher,CIBERSORT,CPM,DAISM,debCAM,Deblender,DeCompress,deconf,DeconICA,DeconPeaker,DeconRNASeq,deconvSeq,DecOT,DeMixT,DESeq2,digitalDLSorter,DSA,dtangle,DWLS,EMeth,EPIC,FARDEEP,ImmuCellAI,LinDeconSeq,Linseed,MCPcounter,MethylResolver,MIXTURE,MOMF,MuSic,MySort,NITUMID,PREDE,quanTIseq,ReFACTor,RNA-Sieve,scaden,SCDC,spatialDWLS,TOAST"
