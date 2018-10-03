# download TCGA datasets from UCSC Xena
install.packages("UCSCXenaTools")

p_load(UCSCXenaTools, tidyverse)

xe = XenaHub(hostName = "TCGA")
xe

xe %>% filterXena(filterDatasets = "clinical") -> xe_clinical

# clinical data query and download
xe_clinical.query = XenaQuery(xe_clinical)
xe_clinical.download = XenaDownload(xe_clinical.query, destdir = "G:/biodata/TCGA/UCSC_Xena/TCGA/Clinical")

# pan-cancer data query and download
xe %>% filterXena(filterDatasets = "HiSeqV2_PANCAN$") -> xe_rna_pancan
xe_rna_pancan
xe_rna_pancan.query = XenaQuery(xe_rna_pancan)
xe_rna_pancan.download = XenaDownload(xe_rna_pancan.query, destdir = "G:/biodata/TCGA/UCSC_Xena/TCGA/RNAseq_Pancan")
