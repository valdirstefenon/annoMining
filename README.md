AnnoMining v3.6 User Manual
A Comprehensive Tool for Secondary Metabolism, Pharmaceutical Potential, and Disease Resistance Analysis in Plant Genomes
Version 3.6 | September 2026

Table of Contents
Overview

What AnnoMining Does

New in Version 3.6

Pipeline Workflow

System Requirements

Installation Guide

Input Files

Output Files

Using the Graphical Interface

Command-Line Usage

Understanding the Results

Disease Resistance Analysis

TNJ Gene Detection

Troubleshooting

Frequently Asked Questions

Citation

Overview
AnnoMining is a Python-based pipeline designed to identify and characterize genes involved in secondary metabolism, pharmaceutical potential, and disease resistance from annotated plant genomes. Version 3.6 introduces comprehensive disease resistance analysis with specialized detection of TNJ (TIR-NBS-Jacalin) genes, a Myrtaceae-specific resistance gene family.

The tool combines multiple evidence sources:

KEGG pathway mapping from EggNOG annotations

PFAM domain analysis from InterProScan

HMMER profile searches against curated domain databases

Multi-layer scoring for pharmaceutical and resistance potential

This multi-layered strategy ensures maximum recovery of genes of interest, including those that are species-specific or poorly characterized in non-model plant species.

What AnnoMining Does
Core Functions
Secondary Metabolism Gene Identification

Detects genes associated with plant secondary metabolism pathways

Covers: Terpenoids, Phenylpropanoids, Alkaloids, Glucosinolates, Cannabinoids, Xanthines

Uses KEGG pathway mapping, PFAM domains, and HMMER searches

Pharmaceutical Potential Assessment

Assigns pharmaceutical potential scores (HIGH, MEDIUM, LOW, NONE)

Predicts compound classes based on domain architecture

Identifies genes with nutraceutical and therapeutic applications

Disease Resistance Analysis ✨ NEW in v3.6

Comprehensive detection of NLR (Nucleotide-binding Leucine-rich Repeat) genes

TNJ gene detection (TIR-NBS-Jacalin) - Myrtaceae-specific resistance genes

Resistance pathway mapping (plant-pathogen interaction, MAPK signaling, hormone signaling)

Resistance class distribution and enrichment analysis

Comprehensive Visualization

Heatmaps of gene-pathway associations

Dot plots of enriched PFAM domains

Pathway completeness analysis

Class distribution charts

Functional networks of high-potential genes

Resistance network visualization

Compound predictions based on domain architecture

Comparative Analysis

Multi-genome comparisons for evolutionary studies

Species ranking by pharmaceutical potential

Pathway presence matrix across species

New in Version 3.6
Major Additions
Feature	Description
Disease Resistance Analysis	Complete new module for identifying and characterizing disease resistance genes
TNJ Gene Detection	Specialized detection of TIR-NBS-Jacalin genes (Myrtaceae-specific)
Resistance Scoring System	Multi-factor scoring for resistance potential (HIGH/MEDIUM/LOW)
Resistance Pathway Mapping	Plant-pathogen interaction, MAPK signaling, hormone signaling
NLR Classification	TIR-NLR, CC-NLR, TNL, TNJ, TNJ-like, NLR-Jacalin
Resistance Network	Functional network visualization of resistance genes
PFAM Enrichment for Resistance	Statistical enrichment of resistance-associated PFAM domains
HMMER Integration	Full HMMER support for resistance domain detection
Improvements
Enhanced GFF/HMMER integration across all analysis modes

Improved validation with mutual exclusivity checks for gene groups

Extended pharmaceutical domain dictionary with more PFAMs

Better error handling and logging

Optimized memory usage for large genomes

Pipeline Workflow
text
┌─────────────────────────────────────────────────────────────────────────────┐
│                              INPUT FILES                                    │
├─────────────────┬─────────────────┬─────────────────┬───────────────────────┤
│    EggNOG       │   InterProScan  │  GFF + Genome   │   HMM Profiles        │
│  annotations    │   annotations   │     FASTA       │   Directory           │
└────────┬────────┴────────┬────────┴────────┬────────┴──────────┬────────────┘
         │                 │                 │                   │
         ▼                 ▼                 ▼                   ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         MERGE ANNOTATION DATA                               │
│  • Unify genes from all sources                                            │
│  • Extract KEGG pathways, KOs, EC numbers, GO terms                        │
│  • Collect PFAM domains                                                    │
└─────────────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         HMMER PROFILE SEARCH                                │
│  • Extract proteins from GFF using gffread                                 │
│  • Search against curated HMM profiles                                     │
│  • Identify genes with secondary, pharma, or resistance domains            │
└─────────────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    SECONDARY METABOLISM CLASSIFICATION                      │
│  • Map KEGG pathways to secondary metabolism classes                       │
│  • Classify genes by metabolite class                                      │
│  • Integrate HMMER-detected genes                                          │
└─────────────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      PHARMACEUTICAL SCORING                                │
│  • Score genes based on pathway, KO, and PFAM evidence                    │
│  • Assign HIGH/MEDIUM/LOW/NONE classification                              │
│  • Predict compound classes                                                │
└─────────────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    DISEASE RESISTANCE SCORING  ✨ NEW                      │
│  • Detect NLR, TIR-NLR, CC-NLR, TNJ genes                                 │
│  • Score based on resistance domains and pathways                         │
│  • Classify resistance gene architectures                                  │
└─────────────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    VISUALIZATION & REPORTING                                │
│  • Generate publication-ready figures                                      │
│  • Create summary tables (CSV, Excel)                                      │
│  • Produce statistics reports                                              │
│  • Build functional networks                                               │
└─────────────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                              OUTPUT FILES                                   │
├─────────────────────────────────────────────────────────────────────────────┤
│  • secondary_statistics.txt          • pharma_statistics.txt               │
│  • secondary_genes_table.csv         • pharma_genes_table.csv              │
│  • resistance_statistics.txt  ✨     • resistance_genes_table.csv ✨       │
│  • *_combined_high_medium_pfam_dotplot.pdf                                 │
│  • *_heatmap_scores.pdf              • *_pharma_heatmap_scores.pdf         │
│  • functional_network.pdf            • resistance_network.pdf ✨           │
│  • compound_predictions.csv          • tnj_analysis.txt ✨                 │
│  • pathway_completeness_analysis.csv                                       │
└─────────────────────────────────────────────────────────────────────────────┘
System Requirements
Minimum Hardware
CPU: 2+ cores recommended

RAM: 8 GB minimum, 16 GB+ recommended for large genomes

Storage: 10 GB free space for input files and outputs

Software Dependencies
Operating System: Linux, macOS, or Windows (with WSL)

Python: 3.8 or higher

Conda: For environment management (recommended)

Required Tools
gffread: For extracting proteins from GFF files

HMMER: For profile searches (hmmsearch)

Perl: Sometimes required for gffread dependencies

tkinter: For GUI (install python3-tk)

Installation Guide
Step 1: Create a Conda Environment (Recommended)
bash
# Create the environment with Python 3.10
conda create -n annomining python=3.10 -y

# Activate the environment
conda activate annomining
Step 2: Install Core Dependencies
bash
# Install bioinformatics tools
conda install -c conda-forge -c bioconda gffread hmmer -y

# Install Python scientific libraries
conda install -c conda-forge matplotlib seaborn pandas numpy scipy statsmodels -y

# Install additional file format support
conda install -c conda-forge openpyxl xlsxwriter -y

# Install tkinter (for GUI)
conda install -c conda-forge tk -y
Step 3: Install HMMER Profiles for Secondary Metabolism and Resistance
bash
# Create profiles directory
mkdir -p hmmer_profiles
cd hmmer_profiles

# Download Pfam database (large file ~400 MB compressed, ~2 GB uncompressed)
wget -O Pfam-A.hmm.gz "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"

# Decompress (this may take a few minutes)
gunzip Pfam-A.hmm.gz

# Extract specific profiles using Python
python3 << 'EOF'
import re

# Secondary metabolism profiles
names = {
    'PF01397': 'Terpene_synthase',
    'PF03936': 'Terpene_synthase_N',
    'PF00067': 'Cytochrome_P450',
    'PF00195': 'Chalcone_synthase',
    'PF00891': 'O_methyltransferase',
    'PF01596': 'SAM_methyltransferase',
    'PF00248': 'Cyclase',
    'PF00494': 'Prenyltransferase',
    'PF02797': 'Flavonoid_3_hydroxylase',
    'PF00201': 'Alkaloid_synthase',
    'PF00514': 'Carotenoid_biosynthesis',
}

# Disease resistance profiles (NEW in v3.6)
resistance_names = {
    'PF00931': 'NB_ARC',
    'PF01582': 'TIR',
    'PF00560': 'LRR',
    'PF01419': 'Jacalin',
    'PF05659': 'RPW8',
    'PF13676': 'TIR_2',
    'PF07723': 'LRR_5',
}

all_names = {**names, **resistance_names}

with open('Pfam-A.hmm', 'r') as f:
    content = f.read()

blocks = re.split(r'(?=HMMER3/f)', content)

for pfam, name in all_names.items():
    for block in blocks:
        if re.search(rf'^ACC\s+{pfam}\.\d+', block, re.MULTILINE):
            with open(f"{name}.hmm", 'w') as f:
                f.write(block)
            print(f"✅ {name}.hmm extracted")
            break

EOF

# Return to main directory
cd ..
Step 4: Download AnnoMining Script
Save the script as annoMining3.6.py in your working directory.

Step 5: Verify Installation
bash
# Check if all tools are available
gffread --version
hmmsearch -h | head -3

# Verify Python libraries
python3 -c "import pandas, numpy, matplotlib, seaborn, scipy, statsmodels; print('✅ All libraries installed')"
Input Files
Mandatory Inputs
File	Format	Description	Source
EggNOG annotations	.emapper.annotations	Functional annotation from EggNOG-mapper	EggNOG-mapper v2
InterProScan annotations	.tsv	Domain predictions from InterProScan	InterProScan v5
Optional but Recommended Inputs
File	Format	Description	Source
GFF file	.gff3	Filtered gene predictions	BRAKER3 + hintsFilter
Genome FASTA	.fasta	Genome assembly file	Assembly pipeline
HMM profiles	.hmm	HMMER profile directory	Downloaded from Pfam
Note: The GFF and genome FASTA are required for HMMER-based detection. Without them, AnnoMining will only use EggNOG and InterProScan annotations.

File Preparation
For EggNOG:
Use the --output_type annotations option when running EggNOG-mapper

The file should contain columns with KEGG pathways, KOs, EC numbers, and PFAM domains

For InterProScan:
Use the --output-format tsv option

The file should include PFAM domain predictions

For GFF:
The GFF should contain gene features with IDs

Best if filtered using hintsFilter to ensure only confident predictions

Output Files
Core Reports
File	Description
*_statistics.txt	Summary statistics including total genes, secondary genes, class distribution
*_secondary_genes_table.csv	Complete list of secondary metabolism genes with annotations
*_pharma_genes_table.csv	Filtered list of genes with pharmaceutical potential
*_resistance_statistics.txt	✨ Disease resistance summary statistics
*_resistance_genes_table.csv	✨ Complete list of resistance genes with classifications
Secondary Metabolism Visualizations
File	Description
*_combined_high_medium_pfam_dotplot.pdf	Dot plot showing enriched PFAM domains in HIGH+MEDIUM genes
*_heatmap_scores.pdf	Heatmap of gene-pathway associations with pharma scores
*_pharma_heatmap_scores.pdf	Heatmap focused on pharmaceutical pathways
*_class_distribution.pdf	Pie charts of metabolite class distribution
*_pharma_score_ranking.pdf	Top 20 genes ranked by pharmaceutical score
*_ko_dotplot.pdf	KEGG Ortholog distribution by metabolite class
*_pathway_completeness.pdf	Genes per pathway analysis
*_functional_network.pdf	Network of functionally similar high-potential genes
*_compound_predictions.pdf	Predicted compound classes based on domain architecture
Disease Resistance Visualizations ✨ NEW
File	Description
*_resistance_score_ranking.pdf	Top 20 genes by disease resistance potential
*_resistance_class_distribution.pdf	Distribution of resistance gene classes (NLR, TNJ, etc.)
*_resistance_pathway_completeness.pdf	Resistance pathway gene counts
*_resistance_donut.pdf	Donut chart of resistance potential distribution
*_resistance_network.pdf	Functional network of resistance genes
*_tnj_analysis.txt	✨ Detailed TNJ gene analysis
*_tnj_architectures.pdf	✨ TNJ gene architecture distribution
*_tnj_confidence.pdf	✨ TNJ gene confidence levels
*_resistance_pfam_dotplot.pdf	✨ Enriched PFAM domains in resistance genes
Data Files
File	Description
*_pathway_completeness_analysis.csv	Detailed pathway statistics
*_functional_network.csv	Edge list for functional network
*_compound_predictions.csv	Compound predictions for high-scoring genes
*_combined_pfam_category_summary.csv	Summary of enriched PFAM categories
*_resistance_network.csv	✨ Resistance network edge list
*_resistance_pfam_enrichment.csv	✨ Statistical enrichment of resistance PFAMs
*_tnj_genes_table.csv	✨ Complete table of TNJ genes
Using the Graphical Interface
Launching the GUI
bash
conda activate annomining
python3 annoMining3.6.py
Main Tabs
1. Secondary Metabolism Tab
Identifies secondary metabolism genes using EggNOG, InterProScan, and HMMER.

Step-by-Step Instructions:

Select Input Files

EggNOG: Click "Browse" and select your .emapper.annotations file

InterProScan: Click "Browse" and select your .tsv file

GFF (Optional): Click "Browse" and select your filtered GFF file

Genome FASTA (Optional): Click "Browse" and select your genome assembly

Configure HMMER (Optional)

Specify the directory containing your HMM profiles (e.g., hmmer_profiles/)

Set Output Directory

Default: secondary_metabolism

Run Analysis

Click the "▶ Run Analysis" button

Progress will appear in the log window

2. Pharma Analysis Tab
For detailed pharmaceutical potential analysis of secondary metabolism genes.

Steps:

Select Input Files - Same EggNOG and InterProScan files

Set Parameters

Threshold Score: Default 4 (lower = more MEDIUM genes)

Network Nodes: Number of genes in functional network (default 15)

Network Method: Selection method (enrichment, score, or balanced)

Run Analysis

3. Disease Resistance Tab ✨ NEW
For comprehensive disease resistance gene analysis.

Steps:

Select Input Files

EggNOG and InterProScan files (as above)

GFF and Genome FASTA for HMMER detection

Set Parameters

Threshold Score: Default 5

Network Nodes: Number of genes in resistance network (default 15)

Network Method: Selection method (enrichment, score, or balanced)

Run Analysis - Click "▶ Run Analysis"

Output includes:

Resistance gene classification (NLR, TIR-NLR, CC-NLR, TNJ, etc.)

TNJ gene detection with architecture and confidence levels

Resistance pathway mapping

Functional network visualization

4. eggnog2kegg Tab
Utility to extract KEGG Orthologs from EggNOG files.

Command-Line Usage
Basic Usage
bash
python3 annoMining3.6.py [arguments]
Secondary Metabolism Analysis
bash
python3 annoMining3.6.py --mode secondary \
    --eggnog eggnog.emapper.annotations \
    --interpro interpro.tsv \
    --gff braker_confident_only.gff3 \
    --genome genome.fasta \
    --hmm_dir hmmer_profiles/ \
    --output_dir secondary_metabolism
Pharma Analysis
bash
python3 annoMining3.6.py --mode pharma \
    --eggnog eggnog.emapper.annotations \
    --interpro interpro.tsv \
    --gff braker_confident_only.gff3 \
    --genome genome.fasta \
    --hmm_dir hmmer_profiles/ \
    --output_dir pharma_analysis \
    --threshold 4 \
    --nodes 15 \
    --method enrichment
Disease Resistance Analysis ✨ NEW
bash
python3 annoMining3.6.py --mode resistance \
    --eggnog eggnog.emapper.annotations \
    --interpro interpro.tsv \
    --gff braker_confident_only.gff3 \
    --genome genome.fasta \
    --hmm_dir hmmer_profiles/ \
    --output_dir resistance_analysis \
    --threshold 5 \
    --nodes 15 \
    --method enrichment
Multi-Genome Analysis
bash
python3 annoMining3.6.py --mode multi \
    --eggnog file1.emapper.annotations,file2.emapper.annotations \
    --interpro interpro1.tsv,interpro2.tsv \
    --names Species1,Species2 \
    --output_dir multi_analysis
Complete Argument Reference
Argument	Description	Required
--mode	Analysis mode: secondary, pharma, resistance, multi, eggnog2kegg	Yes
--eggnog	Path to EggNOG annotations file	Yes
--interpro	Path to InterProScan TSV file	For secondary/pharma
--gff	Path to GFF file	No (recommended)
--genome	Path to genome FASTA file	Required if using GFF
--hmm_dir	Directory with HMM profiles	No
--output_dir	Output directory name	No (defaults vary)
--threshold	Score threshold (default: 4 for pharma, 5 for resistance)	No
--nodes	Network nodes (default: 15)	No
--method	Network selection method (default: enrichment)	No
--names	Comma-separated species names (multi mode)	Yes for multi
Understanding the Results
Pharmaceutical Potential Scores
AnnoMining assigns scores based on:

Pathway contribution (+2 points per pharma-related KEGG pathway)

KEGG Ortholog presence (+3 points per pharma-related KO)

PFAM domain presence (+2 points per core pharma domain)

Classification thresholds (default threshold = 4):

HIGH: Score ≥ 10

MEDIUM: Score ≥ 4 and < 10

LOW: Score ≥ 1 and < 4

NONE: Score = 0

Disease Resistance Scores ✨ NEW
AnnoMining assigns resistance scores based on:

Domain architecture (NLR, TIR, NB-ARC, LRR, Jacalin, etc.)

Resistance class detection (TNJ, TNL, CC-NLR, etc.)

Pathway involvement (plant-pathogen interaction, MAPK signaling)

KO presence (NLR genes, MAPK signaling, defense proteins)

Classification thresholds (default threshold = 5):

HIGH: Score ≥ 10

MEDIUM: Score ≥ 5 and < 10

LOW: Score ≥ 1 and < 5

NONE: Score = 0

TNJ Gene Classification ✨ NEW
TNJ (TIR-NBS-Jacalin) genes are classified based on domain architecture:

Architecture	Domains	Confidence
TNJ canônico	TIR + NB-ARC + Jacalin	HIGH
TNJ-like	TIR + NB-ARC + Jacalin + LRR	MEDIUM
NLR-Jacalin	NB-ARC + Jacalin (no TIR)	MEDIUM
TNL clássico	TIR + NB-ARC + LRR	HIGH
Interpreting the Statistics File
Secondary Metabolism Example:
text
SECONDARY METABOLISM ANALYSIS - SUMMARY
Total genes: 29,698
Secondary genes: 746
Percentage: 2.5%

Genes by class:
  Phenylpropanoids: 335
  Terpenoids: 286
  Alkaloids: 58
  Glucosinolates: 31
  Cannabinoids: 14
  Xanthines: 4

Pharmaceutical potential subset:
  Total pharma genes: 725
  Percentage of secondary: 97.2%
  HIGH: 4
  MEDIUM: 205
  LOW: 516
Disease Resistance Example ✨ NEW:
text
DISEASE RESISTANCE GENE ANALYSIS
Total genes: 29,698

Genes with disease resistance potential: 1,247
  Percentage of total: 4.2%
  HIGH: 23
  MEDIUM: 189
  LOW: 1,035

TNJ genes detected: 12
  Architectures:
    TIR-NB-ARC-Jacalin (TNJ canônico): 8
    TIR-NB-ARC-Jacalin-LRR (TNJ-like): 3
    NB-ARC-Jacalin (NLR com Jacalin): 1

Resistance classes:
  NLR: 345
  TIR_NLR: 189
  LRR_repeat: 178
  NB_ARC: 156
  TNJ: 12
Disease Resistance Analysis
The disease resistance module provides comprehensive characterization of plant resistance genes:

Resistance Gene Classes Detected
Class	Domains	Description
NLR	NB-ARC + LRR	Classic NLR resistance genes
TIR_NLR	TIR + NB-ARC + LRR	TIR-domain containing NLRs
TNJ	TIR + NB-ARC + Jacalin	Myrtaceae-specific resistance genes
TNJ_like	TIR + NB-ARC + Jacalin + LRR	Extended TNJ architecture
CC_NLR	CC + NB-ARC + LRR	Coiled-coil NLRs
TNL	TIR + NB-ARC + LRR	TIR-NBS-LRR
RPW8	RPW8 domains	RPW8-type resistance
NB_ARC	NB-ARC	Nucleotide-binding domain
TIR	TIR	Toll/Interleukin-1 receptor
LRR_repeat	LRR	Leucine-rich repeats
Jacalin	Jacalin	Jacalin lectin domain
NLR_Jacalin	NB-ARC + Jacalin	NLR with Jacalin (no TIR)
Resistance Pathways Mapped
Pathway	KEGG ID	Description
Plant-pathogen interaction	map04626	Plant defense signaling
MAPK signaling - plant	map04016	MAPK cascade
Hormone signal transduction	map04075	SA, JA, ET, ABA signaling
Phenylpropanoid biosynthesis	map00940	Lignin and phytoalexins
Flavonoid biosynthesis	map00941	Flavonoid phytoalexins
Terpenoid backbone	map00900	Terpenoid phytoalexins
TNJ Gene Detection
TNJ genes are a Myrtaceae-specific family of resistance genes characterized by the TIR-NBS-Jacalin domain architecture.

Detection Criteria
Mandatory Domains (for canonical TNJ):

PF01582 (TIR domain)

PF00931 (NB-ARC domain)

PF01419 (Jacalin domain)

Optional Domains:

PF00560 (LRR domain) - for TNJ-like

PF13676 (TIR_2) - alternative TIR domain

TNJ Architecture Types
Type	Domain Composition	Confidence
Canonical TNJ	TIR + NB-ARC + Jacalin	HIGH
TNJ-like	TIR + NB-ARC + Jacalin + LRR	MEDIUM
NLR-Jacalin	NB-ARC + Jacalin (no TIR)	MEDIUM
TNJ Analysis Output
The tnj_analysis.txt file provides:

Total number of TNJ genes detected

Architecture distribution

Confidence levels

Resistance potential distribution

Complete gene list with details

Interpreting TNJ Results
Presence of TNJ genes indicates potential for species-specific resistance mechanisms

Canonical TNJ (TIR-NB-ARC-Jacalin) suggests typical TNJ function

TNJ-like genes may represent evolutionary intermediates or specialized variants

No TNJ genes detected is consistent with the variation observed in Myrtaceae; some species naturally lack these genes

Troubleshooting
Common Issues and Solutions
Issue	Solution
"gffread not found"	Ensure gffread is installed: conda install -c bioconda gffread
"HMMER profiles not found"	Download profiles or specify correct directory. See installation guide.
"No secondary genes found"	Check input files. Ensure EggNOG annotations include KEGG pathways.
"No resistance genes detected"	This may be biologically accurate for some species. Check that resistance domains are present in your HMM profiles.
"GUI does not start"	Install tkinter: conda install -c conda-forge tk
"ValueError: could not convert string to float"	Check HMMER output parsing. Ensure using the corrected version 3.6.
"MemoryError"	Reduce number of genes analyzed or increase RAM.
"TNJ detection not working"	Ensure Jacalin domain (PF01419) is in your HMM profiles.
Debug Mode
To run with more verbose output:

bash
python3 annoMining3.6.py 2>&1 | tee output.log
Checking HMMER Profiles
To verify a profile is valid:

bash
hmmsearch --check hmmer_profiles/TIR.hmm
hmmsearch --check hmmer_profiles/Jacalin.hmm
hmmsearch --check hmmer_profiles/NB_ARC.hmm
Frequently Asked Questions
Q: Why do I need both EggNOG and InterProScan?
A: They provide complementary information. EggNOG gives KEGG pathways and KOs, while InterProScan provides PFAM domains. Combining them improves detection.

Q: What does HMMER add to the analysis?
A: HMMER can detect genes with characteristic domains even if they lack functional annotation in databases. This is crucial for non-model species and for detecting resistance genes that may be species-specific.

Q: What are TNJ genes and why are they important?
A: TNJ (TIR-NBS-Jacalin) genes are a Myrtaceae-specific family of disease resistance genes. They combine TIR and NB-ARC domains with a Jacalin lectin domain. Detection of these genes can indicate species-specific resistance mechanisms.

Q: How do I know if HMMER found new genes?
A: The log will show "HMMER identified X genes with domains". The statistics file will show the total secondary/resistance gene count.

Q: Can I run AnnoMining without the GUI?
A: Yes. Use command-line mode with the --mode argument.

Q: What if my genome has no EggNOG or InterProScan annotations?
A: AnnoMining requires at least EggNOG annotations. Consider running EggNOG-mapper on your proteome first.

Q: How long does the analysis take?
A: For a typical plant genome (~30,000 genes), secondary metabolism analysis takes 5-15 minutes. Resistance analysis adds 5-10 minutes. HMMER searches depend on the number of profiles used.

Q: What does the "Other" category represent?
A: The "Other" category includes pathways that don't fit the main classes. Check the specific pathways in the table.

Q: Can I add my own HMM profiles?
A: Yes. Add your profiles to the hmmer_profiles directory and update the appropriate dictionary in the script.

Q: Why didn't my genome detect any TNJ genes?
A: TNJ genes are specific to Myrtaceae and may not be present in all species. Additionally, some Myrtaceae species naturally lack TNJ genes. Check that your genome has TIR and NB-ARC domains (these are more broadly distributed).

Q: How do I interpret the resistance network?
A: The resistance network shows connections between genes based on shared resistance classes, pathways, and PFAM domains. Node color indicates resistance potential (HIGH/MEDIUM/LOW), and 🔬 indicates TNJ genes.

Citation
If you use AnnoMining in your research, please cite:

text
AnnoMining v3.6: A comprehensive pipeline for secondary metabolism,
pharmaceutical potential, and disease resistance discovery in plant genomes.
License
AnnoMining is distributed under the MIT License.

Version History
Version	Date	Changes
v3.6	Sep 2026	✨ Disease resistance module with TNJ detection
v3.5	Aug 2026	Full GFF/HMMER integration
v3.1	Aug 2026	HMMER support added
v3.0	Jul 2026	Initial release
AnnoMining v3.6 | Making secondary metabolism, pharmaceutical, and disease resistance discovery accessible to everyone

