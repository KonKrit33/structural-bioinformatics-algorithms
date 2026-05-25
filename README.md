# Structural Bioinformatics Algorithms

This repository contains MSc coursework projects in Algorithms in Structural Bioinformatics, focusing on sequence comparison, protein structure comparison, molecular visualization, RMSD analysis, distance geometry, and protein backbone reconstruction.

The projects combine algorithmic bioinformatics, structural biology, protein modelling, and computational analysis. They are relevant to computational drug discovery, protein structure analysis, molecular modelling workflows, and AI-ready structural bioinformatics pipelines.

## Project Overview

The repository includes three assignments:

### Assignment 1: Sequence and Protein Structure Comparison

This assignment focuses on both sequence-level and structure-level comparison of proteins.

The sequence comparison part evaluates four classical protein sequence search and alignment methods:

- Needleman-Wunsch global alignment
- Smith-Waterman local alignment
- FASTA heuristic sequence search
- BLAST heuristic sequence search

The analysis includes:

- Generation of a controlled toy protein database
- Comparison of strong homologs, weak homologs, random sequences, fragment queries, and motif queries
- Runtime comparison of exact dynamic programming methods and heuristic methods
- Sensitivity analysis of FASTA and BLAST against dynamic programming methods
- BLAST seed-length sensitivity analysis using different word sizes
- Top-hit and top-5 ranking comparison using Spearman rank correlation

The structural comparison part uses protein structures 6LU7 and 2AMD and includes:

- Extraction of Cα atoms from residues 50-99
- Optimal rigid-body alignment using rotation and translation
- Coordinate RMSD calculation
- Distance RMSD calculation using all pairwise Cα distances
- Approximate d-RMSD estimation using subsets of pairwise distances

### Assignment 2: VanA Protein Structure Visualization and RMSD Analysis

This assignment focuses on molecular structure inspection and visualization of ancient VanA using UCSF ChimeraX.

The analysis includes:

- Identification of the experimental structure determination method and resolution
- Inspection of protein chains, residues, water molecules, and ligands
- Cartoon visualization of the overall VanA structure
- Molecular surface visualization of the protein
- Ball-and-stick visualization of ligands
- Analysis of ligand-protein hydrogen bonding interactions
- Visualization of Mg²⁺ coordination interactions
- Structural superposition of VanA structures 3SE7 and 1E4E
- RMSD analysis using ChimeraX Matchmaker
- Focused comparison of the Mg²⁺ metal-binding site

### Assignment 3: Protein Backbone Reconstruction from Distance Geometry

This assignment focuses on distance-based reconstruction of protein backbone geometry.

The analysis uses the Cα coordinates of residues 102-152 from chain A of PDB structure 6LU7 and includes:

- Construction of a Gram matrix from pairwise distance data
- Singular value decomposition of the Gram matrix
- Rank analysis of exact and noisy distance matrices
- Reconstruction of 3D coordinates from distance data
- Coordinate RMSD calculation after optimal alignment
- Analysis of reconstruction ambiguity under reflection
- Reconstruction under random distance noise
- Reconstruction using different reference atoms
- Reconstruction from sparse partial contact-like distances
- Optimization-based refinement using distance bounds and backbone regularization

## Repository Structure

```text
.
├── assignment_1_sequence_and_structure_comparison/
│   ├── report/
│   │   └── Assignment_1_Sequence_and_Structure_Comparison_Report.pdf
│   ├── sequence_comparison/
│   │   ├── a_code.py
│   │   ├── b_code.py
│   │   ├── c_code.py
│   │   ├── d_code.py
│   │   ├── e_code.py
│   │   ├── db.fasta
│   │   ├── queries.fasta
│   │   ├── reference.fasta
│   │   ├── a_results/
│   │   ├── b_results/
│   │   ├── c_results/
│   │   ├── d_results/
│   │   └── e_results/
│   └── structural_comparison/
│       ├── 2a.py
│       ├── 2b.py
│       ├── 2c.py
│       ├── 6LU7.pdb
│       └── 2AMD.pdb
│
├── assignment_2_vana_structure_visualization/
│   └── Assignment_2_VanA_Structure_Visualization_Report.pdf
│
├── assignment_3_backbone_reconstruction_distance_geometry/
│   └── Assignment_3_Backbone_Reconstruction_Distance_Geometry_Report.pdf
│
├── README.md
└── .gitignore
```

## Included Assignments

### Assignment 1

`Assignment_1_Sequence_and_Structure_Comparison_Report.pdf`

This report documents sequence comparison experiments using Needleman-Wunsch, Smith-Waterman, FASTA, and BLAST, together with structural comparison of protein fragments from 6LU7 and 2AMD. It includes runtime analysis, sensitivity comparison, Spearman rank correlation, Kabsch-style alignment, c-RMSD, d-RMSD, and approximate distance-based comparison.

### Assignment 2

`Assignment_2_VanA_Structure_Visualization_Report.pdf`

This report documents structural analysis of ancient VanA using UCSF ChimeraX. It includes structure metadata, ligand inspection, molecular visualization, hydrogen-bond and metal-coordination analysis, structural superposition of 3SE7 and 1E4E, RMSD analysis, and Mg²⁺ binding-site comparison.

### Assignment 3

`Assignment_3_Backbone_Reconstruction_Distance_Geometry_Report.pdf`

This report documents distance-geometry reconstruction of a protein Cα backbone from pairwise distances. It includes Gram matrix construction, SVD, rank analysis, reconstruction from exact and noisy distances, reference-atom comparison, sparse contact-based reconstruction, and optimization refinement with distance bounds and backbone regularization.

## Tools and Technologies

- Python
- NumPy
- FASTA36
- BLAST+
- Needleman-Wunsch alignment
- Smith-Waterman alignment
- Spearman rank correlation
- PDB structure analysis
- UCSF ChimeraX
- Protein structure visualization
- Kabsch alignment
- c-RMSD
- d-RMSD
- Distance geometry
- Gram matrix reconstruction
- Singular value decomposition
- Optimization-based coordinate refinement
- Structural bioinformatics

## Relevance

This repository demonstrates practical and algorithmic skills in structural bioinformatics, including protein sequence comparison, protein structure alignment, molecular visualization, and distance-based 3D reconstruction.

The work is relevant to:

- Computational drug discovery
- Protein structure analysis
- Molecular modelling
- Structural bioinformatics pipelines
- Protein similarity search
- Structure-based comparison of biomolecules
- Distance geometry and conformational reconstruction
- Scientific computing for biomolecular data

## Notes

The repository contains solution reports, code, selected input files, and selected output files. Large intermediate files or unnecessary temporary files are not included.

System-generated files such as `.DS_Store`, cache folders, and temporary outputs are intentionally excluded.
