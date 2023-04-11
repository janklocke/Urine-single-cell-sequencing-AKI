# Urine-single-cell-sequencing-AKI

This repository contains code and data from the publication "Urinary single-cell sequencing captures intrarenal injury and repair processes in human acute kidney injury", currently uploaded as a preprint to biorxiv: https://doi.org/10.1101/2022.02.15.479234 

----
You'll find several scripts containing the code for analysis and plotting found in the manuscript:

1) "AKI_Urine_Sediment_SO_script" contains the code for generating the whole data Seurat object ("SO_all_urine_cells.rds"/"URINE") and subsetted renal cell Seurat object ("SO_all_kidney_cells.rds"/"RENAL") from raw data
2) "AKI_Urine_Sediment_Demultiplexing_pooled_samples_script" contains the code used to demultiplex barcoded and pooled urine samples as well as the code for Suppl.Fig. 3
3) "AKI_Urine_Sediment_Figure1_script" contains the code for Fig. 1 and Suppl. Fig. 4-5
4) "AKI_Urine_Sediment_Figure2_script" contains the code for Fig. 2 and Suppl. Fig. 6-8, 17
5) "AKI_Urine_Sediment_Figure3_script" contains the code for Fig. 3 and Suppl. Fig. 9, 10, 13
6) "AKI_Urine_Sediment_Figure4_script" contains the code for Fig. 4 
7) "AKI_Urine_Sediment_Figure5_script" contains the code for Fig. 5 and Suppl. Fig. 11-12
8) "AKI_Urine_Sediment_Figure6_script" contains the code for Fig. 6 and Suppl. Fig. 14
9) "AKI_Urine_Sediment_SupplFig15_script" contains the code for Suppl. Fig. 15
10) "AKI_Urine_Sediment_SupplFig16-20_script" contains the code for Suppl. Fig. 16, 18-20
11) "AKI_Urine_Sediment_SupplFig21_script" contains the code for Suppl. Fig. 21

----
Some additional data is provided in form of these documents: 

DEG_all_urine_cells.rds  	= Table with differentially expressed genes from SO_all_urine_cells.rds  
DEG_kidney_urine_cells.rds  	= Table with differentially expressed genes from SO_kidney_urine_cells.rds  

Urine_AKI_patient_data.csv 		= patient information added as metadata in "AKI_Urine_Sediment_SO_script"
Urine_AKI_barcode_list.csv 		= barcode patient information used in "AKI_Urine_Sediment_Demultiplexing_pooled_samples_script"
Urine_AKI_urine_output_crea_data.csv	= longitudinal info on serum creatinine and urine output for patients sampled multiple times (fig. 6)

----
The raw data is uploaded to the NCBI GEO repository GSE199321. Seurat objects used for the analysis and generation of figures can be found here:
Seurat object containing all urine cells post QC: https://figshare.com/articles/dataset/KI_2022_Klocke_et_al_SO_all_urine_cells_rds/22567201
Seurat object containing only kidney parenchymal cells post QC: https://figshare.com/articles/dataset/KI_2022_Klocke_et_al_SO_kidney_urine_cells_rds/22567195


