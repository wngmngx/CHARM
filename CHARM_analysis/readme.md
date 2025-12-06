# CHARM Analysis

This directory contains the main analysis code for the CHARM (Chromatin Architecture Mapping) project, divided into two core components corresponding to different manuscript figures.

## 📁 mesc_related/
Mouse embryonic stem cell (mESC) related data analysis code, primarily corresponding to **Main Figures 1-3** in the manuscript.

### Part 1: Data Preprocessing
- `1.1_generate_pairs_info.ipynb` - Generate pairing information from CHARM data and calculate descriptive statistics like near% and mitotic% from Hi-C pairs
- `1.2_create_seurat_object.ipynb` - Create Seurat objects for integrated single-cell multi-omics analysis
- `1.3_calc_tss_enrichment.ipynb` - Document Signac’s TSS enrichment implementation and adjust it to an ENCODE/ArchR-style metric used throughout the manuscript.

### Part 2: Validation of CHARM data
- `2.1_sox2_scr_distance.ipynb` - Analyze spatial distances and organization of Sox2 locus
- `2.2_validationPython.ipynb` - Python-based code for method validation.
- `2.3_validationR.ipynb` - R-based code for method validation.
- `2.4_bivalent_bins.ipynb` - Analysis of bivalent chromatin domains and their properties.
- `2.5_ordermat_python.ipynb` - Python implementation for ordering and clustering chromatin interaction matrices, related to figure2.
- `2.6_genome_structure_impact_on_me3_restore.ipynb` - Estimate how 3D chromatin structure modulates H3K27me3 restoration by propagating contact-derived exposure scores and correlating them with ΔM across cell-cycle phases.

### Part 3: 3D structure Analysis
- `3.1_reconstruct_structure_validation.ipynb` - Validate 5 kb reconstructed structures through Rabl polarity metrics, radial positioning, and agreement with bulk Hi-C before exporting exemplar conformations.
- `3.2_spatial_clustering_1_cr.ipynb` - Build allele-resolved Cell3D objects with accessibility/H3K27me3/RNA annotations, classify particle types, and compute spatial cross-correlations against RNA and bivalent peaks.
- `3.2_spatial_clustering_2_gr.ipynb` - Perform pair-correlation (G(r)) analysis—including radial position controls and randomized backgrounds—to quantify local enrichments of accessibility, H3K27me3, and RNA particles.
- `3.3_cluster.ipynb` - DBSCAN clustering analysis of accessibility cluster.
- `3.4_local_enrich_region.ipynb` - Identification and analysis of locally enriched chromatin regions(2D genome, adapt from "rose" pipeline originally for calling super enhancer)
- `3.5_co_expression.ipynb` - Test whether co-expressed gene pairs preferentially fall within the same 3D accessibility clusters by aggregating single-cell odds ratios against shuffled controls.
- `my_rose_py3.py` - Python3 re-implementation of ROSE algorithm for super-enhancer identification. 

---

## 📁 brain_related/
Brain dataset related analysis code, primarily corresponding to **Main Figures 4-5** in the manuscript.

### Part 1: Object Creation
- `1.1_create_obj.ipynb` - Create and initialize Seurat objects for brain single-cell multi-omics data

### Part 2: Gene Type Analysis (Extended Data Fig10)
- `2.1_generate_different_types_genes.ipynb` - Classify and generate gene sets based on expression patterns(Cell type-specific or housekeeping).
- `2.2_cluster_enrichment.ipynb` - Perform enrichment analysis across different gene categories.

### Part 3: Correlation Analysis (Figure4)
- `3.1_createobj.ipynb` - Genearate Cell3D and MultiCell3D object for correlation analysis.
- `3.2_createmetacells.ipynb` - Generate metacells to reduce noise(mESC + brian).
- `3.3_calccor.py` - Calculate comprehensive correlation matrices between genes and regulatory elements
- `3.4_cor_post_process_AUC.ipynb` - Post-process correlation results and calculate AUC metrics
- `3_generatedata.py` - Generate regression-ready inputs per gene by exporting contact-normalized accessibility/H3K27me3 features around each TSS across all cells.

### Part 4: Regression Analysis (Figure5)
- `4.1_createmetacells.ipynb` - Create metacells with improved parameters(brain dataset only).
- `4.2_createobj.ipynb` - Similar as 3.1_createobj.ipynb, but for brain dataset.
- `4.3_calccor.py` - Similar as 3.3_calccor.py, but for brain dataset.
- `4.4_generatedata.py` - Similar as 3.5_generatedata.py, but for brain dataset.
- `4.5_cor_post_process.ipynb` - Similar as 3.4_cor_post_process_AUC.ipynb, but for brain dataset.
- `4.6_regression_model.ipynb` - Build regression model for gene regulation prediction.
- `4.7_shap_describe.ipynb` - SHAP (SHapley Additive exPlanations) analysis for model interpretability
- `4.8_split_frag_by_subclass.py` - Split Li2021 fragments into subclass-specific files by mapping CellIDs to subclasses for downstream motif modeling.
- `4.8_run_chrombpnet.sh` - Wrapper script to filter subclass fragments and launch ChromBPNet training with the appropriate genome, peaks, and bias models.
- `4.9_compare_with_scarlink` - Compare CHARM regression model results with SCARlink predictions for validation and benchmarking.   
- `4.10_compare_with_encode_re2g` - Compare CHARM regression model results with ENCODE_re2g model predictions for validation and benchmarking.
