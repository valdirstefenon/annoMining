annoMining
Version 2.2
________________________________________
Table of Contents
1.	Introduction
2.	Installation and Requirements
3.	Overview of Methodology
4.	Input File Formats
5.	Analysis Modules
o	5.1 Secondary Metabolism Analysis
o	5.2 Pharmaceutical & Nutraceutical Analysis
o	5.3 Multi-Genome Comparative Analysis
o	5.4 eggNOG to KEGG Converter
6.	Output Files and Interpretation
7.	Data Sources and Reference Databases
8.	Scoring and Ranking Criteria
9.	Troubleshooting
10.	Citation
________________________________________
 
Introduction
AnnoMining is an integrated computational platform designed for comprehensive analysis of secondary metabolism genes from genomic annotations. The software leverages functional annotation data from eggNOG-mapper and InterProScan to identify, classify, and prioritize genes involved in the biosynthesis of secondary metabolites with pharmaceutical and nutraceutical potential.
The platform is particularly valuable for researchers working with plant genomes, microbial genomes, or any organism where secondary metabolite production is of interest for drug discovery, natural product research, or metabolic engineering.

Key Features

•	Secondary Metabolism Mapping: Identifies genes involved in known secondary metabolite pathways (terpenoids, phenylpropanoids, alkaloids, glucosinolates, cannabinoids, and xanthines)

•	Pharmaceutical Potential Scoring: Ranks genes by their predicted pharmaceutical and nutraceutical relevance using a multi-criteria weighted system

•	Score-Based Heatmaps: Visualizes gene-pathway associations with continuous color gradients reflecting pharmaceutical potential scores

•	Multi-Genome Comparison: Enables comparative analysis across multiple species or strains with pathway score matrices

•	Pathway Completeness Analysis: Evaluates the completeness of each metabolic pathway based on gene counts and pharmaceutical scores

•	Publication-Ready Visualizations: Generates high-resolution figures with legends positioned outside plots to prevent overlap

•	Interactive GUI: User-friendly interface for researchers without command-line expertise

•	Batch Processing: Command-line mode for high-throughput analysis

________________________________________
 
Installation and Requirements

System Requirements

•	Operating System: Linux, macOS, or Windows (with Python 3.8+)

•	Memory: 4GB RAM minimum (8GB+ recommended for large genomes)

•	Disk Space: 1GB for software + space for output files


Dependencies

AnnoMining requires the following Python packages:


python >= 3.8

pandas >= 1.3.0

numpy >= 1.21.0

matplotlib >= 3.4.0

seaborn >= 0.11.0

scipy >= 1.7.0

statsmodels >= 0.12.0

tkinter (typically included with Python)


Installation

Option 1: From GitHub

bash

git clone https://github.com/valdirstefenon/annoMining.git

cd annomining

pip install -r requirements.txt

python annomining.py


Option 2: Direct Download

Download the annomining.py file and run with Python:

bash

python annomining.py



Verifying Installation

To verify that all dependencies are correctly installed, run:

bash

python -c "import pandas, numpy, matplotlib, seaborn, scipy, statsmodels; print('All dependencies OK')"

________________________________________
 
Overview of Methodology
AnnoMining operates through a multi-step workflow that transforms raw functional annotations into prioritized lists of pharmaceutically relevant genes.

Core Workflow
text
Input Files
    ↓
[Parser] - Extracts KOs, ECs, GO, PFAM, Pathways
    ↓
[Integration] - Merges eggNOG and InterProScan data
    ↓
[Pathway Mapping] - Maps genes to secondary metabolic pathways
    ↓
[Scoring] - Calculates pharmaceutical potential scores
    ↓
[Validation] - Validates group exclusivity (HIGH/MEDIUM/LOW)
    ↓
[Visualization] - Generates heatmaps, networks, rankings, completeness analysis
    ↓
Output Files

Conceptual Framework
1.	Secondary Metabolism Classification: Genes are classified into seven major classes:
o	Terpenoids (isoprenoid-based compounds)
o	Phenylpropanoids (aromatic amino acid derivatives)
o	Alkaloids (nitrogen-containing compounds)
o	Glucosinolates (sulfur-containing compounds)
o	Xanthines (purine alkaloids)
o	Cannabinoids (terpenophenolic compounds)
o	Other (miscellaneous secondary metabolites)
2.	Pharmaceutical Potential Assessment: Each gene is scored based on:
o	Presence in known pharmaceutical pathways (2 points per pathway)
o	Association with pharmaceutically relevant KOs (3 points per KO)
o	Co-occurrence with core PFAM domains (2 points per domain)
o	Classification into HIGH (≥10), MEDIUM (5-9), LOW (1-4), or NONE (0)
3.	Comparative Analysis: Multiple genomes are compared using:
o	Average pharmaceutical scores per pathway
o	Compound class diversity metrics
o	Pathway completeness indices
o	Network-based similarity metrics
4.	PFAM Enrichment Analysis: Compares HIGH+MEDIUM genes against LOW+NONE background to identify significantly enriched domains
5.	Pathway Completeness Analysis: Evaluates each pathway's gene count and pharmaceutical score distribution
________________________________________
 
Input File Formats
1. eggNOG-mapper Output
AnnoMining accepts the standard eggNOG-mapper annotation file format (.emapper.annotations).
Expected Format:
text
#query_name    seed_ortholog    e-value    score    eggNOG_OG    max_annot_lvl    COG_category    Description    Preferred_name    GOs    ECs    KEGG_ko    KEGG_Pathway    KEGG_Module    ...
gene1    XXX   1e-150   100    OG5_123456    1    S    hypothetical protein    -    GO:0001234   3.4.11.1    K00001    map00900    M00001
gene2    YYY   1e-140   95    OG5_234567    1    C    kinase    -    GO:0005678   2.7.1.1    K00002    map00901    M00002
Critical Columns:
•	Column 1: Gene/query identifier
•	Column 8: COG category (optional)
•	Column 10: Gene Ontology (GO) terms
•	Column 11: EC numbers
•	Column 12: KEGG Orthologs (KO) numbers
•	Column 13: KEGG pathway identifiers (mapXXXXX)
The parser is flexible and can extract KEGG information from multiple column formats (including ko:KXXXXX or just KXXXXX).

2. InterProScan Output (Optional)
InterProScan TSV format files can be provided for additional domain and GO term information.
Expected Format (InterProScan TSV):
text
#protein_acc    md5    length    analysis    signature_acc    signature_desc    start    end    score    status    date    interpro_acc    interpro_desc    ...
gene1    ABC123    456    Pfam    PF00195    ABC transporter    1    200    45.6    T    2024-01-01    IPR001234    ABC transporter domain
Key Fields Used:
•	Column 1: Gene identifier
•	Column 5: Domain signature (PFAM)
•	Column 12: InterPro domain
Note: InterProScan file is optional but enhances predictions through additional PFAM domain data.
________________________________________
 
Analysis Modules
5.1 Secondary Metabolism Analysis
Purpose: Identify and characterize genes involved in secondary metabolite biosynthesis.
Workflow:
1.	Parse eggNOG and InterProScan files
2.	Extract KEGG pathways, KO numbers, EC numbers, GO terms, and PFAM domains
3.	Map genes to secondary metabolism pathways (see SECONDARY_PATHWAYS dictionary)
4.	Classify genes into metabolite classes
5.	Calculate pharmaceutical potential scores for all secondary metabolism genes
6.	Validate group exclusivity (HIGH, MEDIUM, LOW)
7.	Generate comprehensive visualizations and tables
8.	Perform pathway completeness analysis
9.	Generate PFAM enrichment dotplots (combined HIGH+MEDIUM vs LOW+NONE)
Configuration Parameters:
Parameter	Description	Default
eggNOG file	eggNOG annotations file	Required
InterPro file	InterProScan TSV file (optional)	None
Output folder	Directory for output files	secondary_metabolism



Outputs Generated:
File	Description
secondary_heatmap_scores.pdf/png	Gene-pathway heatmap with pharma scores (white→red gradient)
secondary_pharma_heatmap_scores.pdf/png	Pharmaceutical pathway heatmap with scores
secondary_ko_dotplot.pdf/png	Top KOs by metabolite class
secondary_top_kos_by_class.pdf/png	Most abundant KOs by class
secondary_pathway_completeness.pdf/png	Distribution of genes across pathways
secondary_class_distribution.pdf/png	Pie/donut charts (legend outside plot)
secondary_upset_venn.pdf/png	Gene overlap between metabolite classes
secondary_top_ec_numbers.pdf/png	Most frequent EC numbers
secondary_top_core_pfam_domains.pdf/png	Most frequent PFAM domains
secondary_combined_high_medium_pfam_dotplot.pdf/png	PFAM enrichment: HIGH+MEDIUM vs LOW+NONE
secondary_pathway_completeness_analysis.csv	Detailed pathway completeness data
secondary_secondary_genes_table.csv/xlsx	Detailed gene-level data table
secondary_statistics.txt	Summary statistics including pharmaceutical subset
________________________________________
 
5.2 Pharmaceutical & Nutraceutical Analysis
Purpose: Identify and rank genes with pharmaceutical and nutraceutical potential using a multi-criteria scoring system.
Workflow:
1.	Parse and integrate annotation data
2.	Map to pharmaceutical-relevant pathways (see PHARMA_DETAILS)
3.	Score genes using weighted criteria (pathways ×2, KOs ×3, PFAM domains ×2)
4.	Classify genes as HIGH, MEDIUM, LOW, or NONE potential
5.	Validate group exclusivity with automatic correction
6.	Generate prioritization visualizations and tables
7.	Predict compound classes based on domain composition
8.	Perform functional network analysis
9.	Perform pathway completeness analysis
10.	Generate PFAM enrichment dotplots

Configuration Parameters:
Parameter	Description	Default
eggNOG file	eggNOG annotations file	Required
InterPro file	InterProScan TSV file (optional)	None
Threshold	Minimum score for pharmaceutical potential	5
Network nodes	Number of genes in functional network	15
Network method	Node selection method	enrichment
Output folder	Directory for output files	pharma_analysis
Network Methods:
•	score: Selects top genes by pharmaceutical score
•	enrichment: Selects genes based on combined importance metrics (score×2 + KOs×3 + pathways + PFAMs×0.5)
•	balanced: Samples from multiple compound classes
Outputs Generated:
File	Description
pharma_heatmap_scores.pdf/png	Gene-pathway heatmap with pharma scores (white→red gradient)
pharma_pharma_heatmap_scores.pdf/png	Pharmaceutical pathway heatmap with scores
pharma_score_ranking.pdf/png	Top-ranked genes by pharmaceutical potential
pharma_compound_class_distribution.pdf/png	Distribution of bioactive compound classes
pharma_pharmaceutical_uses.pdf/png	Predicted pharmaceutical applications
pharma_donut_potential.pdf/png	Distribution of potential levels
pharma_combined_high_medium_pfam_dotplot.pdf/png	PFAM enrichment: HIGH+MEDIUM vs LOW+NONE
pharma_functional_network.pdf/png	Functional interaction network
pharma_compound_predictions.pdf/png	Predicted compound-producing genes
pharma_pathway_completeness_analysis.csv	Detailed pathway completeness data
pharma_genes_table.csv/xlsx	Detailed gene-level pharmaceutical data
pharma_statistics.txt	Summary pharmaceutical statistics with consistency checks
pharma_compound_predictions.csv	Compound class predictions
pharma_combined_pfam_enrichment.csv	PFAM enrichment analysis results
pharma_combined_pfam_category_summary.csv	PFAM enrichment summary by category
pharma_functional_network.csv	Network edge data
________________________________________
5.3 Multi-Genome Comparative Analysis
Purpose: Compare pharmaceutical potential across multiple genomes or strains.
Workflow:
1.	Process each genome independently using the pharma analysis pipeline
2.	Validate group exclusivity for each genome
3.	Generate comparative statistics for each species
4.	Create multi-species visualizations with score-based matrices
5.	Build species similarity networks
6.	Rank species by pharmaceutical potential metrics
Configuration Parameters:
Parameter	Description	Default
eggNOG files	Comma-separated eggNOG files	Required
InterPro files	Comma-separated InterPro files (optional)	None
Species names	Comma-separated species labels	Required
Threshold	Pharmaceutical score threshold	5
Output folder	Directory for output files	pharma_multi_analysis
Input Requirements:
•	All files must be in the same order (file 1 corresponds to species 1, etc.)
•	Species names should be provided in the same order
•	InterPro files can be omitted (use empty string: ,,,)
Outputs Generated:
File	Description
multi_comparative_summary.csv/xlsx	Cross-species summary statistics
multi_pathway_score_matrix.pdf/png	Heatmap of pathway scores (white→red gradient)
multi_comparative_bars.pdf/png	Stacked bar chart of potential levels
multi_comparative_radar.pdf/png	Multi-dimensional species comparison
multi_class_score_heatmap.pdf/png	Average scores by compound class (white→red)
multi_species_network.pdf/png	Species similarity network
multi_rankings.pdf/png	Species ranking by various metrics
multi_ranking.csv	Ranking data table
multi_comparative_combined_pfam_dotplot.pdf/png	Comparative PFAM enrichment
multi_comparative_report.txt	Comprehensive comparative report
________________________________________
5.4 eggNOG to KEGG Converter
Purpose: Extract KEGG Ortholog (KO) numbers from eggNOG annotations for use with KEGG Mapper tools.
Workflow:
1.	Parse eggNOG annotation file
2.	Extract KO numbers from the KEGG_ko column
3.	Generate a simple two-column format: gene → KO
Input/Output:
•	Input: eggNOG .emapper.annotations file
•	Output: Tab-delimited file with gene IDs and KO numbers
Output Format:
text
gene1    K00001
gene1    K00002
gene2    K00003
Usage Example: This converted file can be used directly with the KEGG Mapper online tool for pathway enrichment analysis.
________________________________________
Output Files and Interpretation
Score-Based Heatmaps (NEW in v2.2)
Gene-Pathway Association Heatmaps with Scores:
•	Rows represent genes with their pharmaceutical score in parentheses
•	Columns represent pathways
•	White-to-red gradient: white = low score, red = high score
•	Intensity reflects quantitative pharmaceutical potential, not just presence/absence
•	More informative than binary matrices for identifying high-potential gene-pathway associations
•	Legends positioned outside the plot to prevent overlap
PFAM Enrichment Analysis (NEW in v2.2)
Combined HIGH+MEDIUM vs LOW+NONE:
•	Compares all genes with pharmaceutical potential (HIGH+MEDIUM) against those without (LOW+NONE)
•	Eliminates duplicate PFAM entries that appeared in separate HIGH/MEDIUM analyses
•	Identifies domains significantly enriched in pharmaceutically relevant genes
•	Results sorted by fold enrichment
•	Includes p-value adjustment (FDR) for multiple testing
Visualization:
•	Dotplot with fold enrichment on x-axis
•	PFAM domains on y-axis
•	Point size proportional to gene count
•	Color-coded by pathway category (Flavonoids, Terpenoids, Alkaloids, etc.)
•	Statistical significance indicated by p.adjust < 0.05
Pathway Completeness Analysis (NEW in v2.2)
For each secondary metabolism pathway, the analysis provides:
•	Gene Count: Number of genes associated with the pathway
•	Pharma Score Sum: Total pharmaceutical score across all genes in the pathway
•	High/Medium Gene Count: Number of HIGH and MEDIUM potential genes
•	Average Pharma Score: Mean pharmaceutical score per gene
Output Files:
•	*_pathway_completeness_analysis.csv: Complete data table
•	*_pathway_completeness_analysis.pdf/png: Visualization of top pathways by gene count
Score Rankings
Pharmaceutical Score Ranking:
•	Higher scores indicate greater pharmaceutical relevance
•	Color coding:
o	🟢 Green (HIGH): Score ≥ 10
o	🟡 Yellow (MEDIUM): Score 5–9
o	🔴 Red (LOW): Score 1–4
•	Genes with score 0 are excluded from ranking
Functional Networks
Network Visualization:
•	Nodes represent genes
•	Node size reflects pharmaceutical score
•	Edge thickness represents functional similarity
•	Edge weight calculated from shared:
o	Pathways (weight × 2)
o	KOs (weight × 3)
o	GO terms (weight × 1)
Compound Predictions
Prediction Confidence:
•	Based on PFAM domain presence
•	Classes predicted:
o	Flavonoids: PF00195, PF02797, PF02458
o	Terpenoids: PF03936, PF01397, PF00067
o	Alkaloids: PF01596, PF00891, PF00248
o	Phenylpropanoids: PF00195, PF00430
o	Carotenoids: PF00514, PF06444
•	Confidence score = (matched domains / required domains) × class confidence
Comparative Reports
Species Rankings:
•	Percent_Pharma_of_Secondary: Percentage of secondary metabolism genes with pharmaceutical potential (more meaningful than percentage of total)
•	Unique_Compound_Classes: Diversity of predicted compound classes
•	Avg_Pharma_Score: Average score across pharma-relevant genes
•	Best Performer identified based on highest percent of pharma genes relative to secondary metabolism
________________________________________
Data Sources and Reference Databases
KEGG Pathways (Secondary Metabolism)
AnnoMining includes a curated set of KEGG pathways relevant to secondary metabolism:
Pathway	Class	Description
map00900	Terpenoids	Terpenoid backbone biosynthesis
map00902	Terpenoids	Monoterpenoid biosynthesis
map00904	Terpenoids	Diterpenoid biosynthesis
map00909	Terpenoids	Sesquiterpenoid/triterpenoid biosynthesis
map00940	Phenylpropanoids	Phenylpropanoid biosynthesis
map00941	Phenylpropanoids	Flavonoid biosynthesis
map00950	Alkaloids	Isoquinoline alkaloid biosynthesis
map00960	Alkaloids	Tropane/piperidine biosynthesis
map00966	Glucosinolates	Glucosinolate biosynthesis
map00232	Xanthines	Caffeine metabolism
map00903	Cannabinoids	Cannabinoid biosynthesis
Pharmaceutical Pathways
Extended pathway set for pharmaceutical and nutraceutical applications:
Pathway	Class	Bioactive Compounds	Applications
map00950	Alkaloids	Morphine, codeine, berberine	Analgesic, antimicrobial
map00909	Terpenoids	Artemisinin, taxol	Antimalarial, anticancer
map00941	Flavonoids	Quercetin, catechins	Cardioprotective
map00945	Stilbenoids	Resveratrol, curcumin	Chemopreventive
map00966	Glucosinolates	Sulforaphane	Anticancer
map00903	Cannabinoids	THC, CBD	Analgesic, anti-epileptic
map00232	Xanthines	Caffeine	Stimulant
Key KEGG Orthologs
Important KOs for pharmaceutical biosynthesis:
Category	KOs
Bioactive flavonoids	K00660, K05275, K05276, K13065, K13066, K05265, K05266
Anticancer terpenes	K15891, K15892, K15893, K00487, K00507, K00509
Nutraceutical carotenoids	K00514, K00515, K06444, K06445, K06446, K06447, K06448, K06449, K06450, K06451
Therapeutic alkaloids	K01799, K01800, K01900, K01901, K01902, K01903
Core PFAM Domains
Pharmaceutically relevant PFAM domains:
Category	Domains
Flavonoids	PF00195, PF02797, PF05834
Terpenoids	PF01397, PF03936, PF00067, PF00494
Alkaloids	PF01596, PF00891, PF00201, PF00155, PF00141
Phenylpropanoids	PF00195, PF00141
Carotenoids	PF00494, PF00514
P450	PF00067
Methyltransferases	PF01596, PF00891, PF01799
________________________________________
Scoring and Ranking Criteria
Pharmaceutical Potential Score Calculation
The pharmaceutical score for each gene is calculated as:
text
Score = (Pathway Contribution × 2) + (KO Contribution × 3) + Domain Bonus
Pathway Contribution (×2 each):
•	Gene associated with a pathway in PHARMA_DETAILS → +2
•	Gene associated with multiple pharma pathways → cumulative
KO Contribution (×3 each):
•	Gene contains a KO listed in PHARMA_KOS → +3
•	Multiple important KOs → cumulative
Domain Bonus:
•	Gene contains a PFAM domain in CORE_PHARMA_DOMAINS → +2 per domain
•	Only non-generic domains (not in GENERIC_DOMAINS) are considered
Potential Classification
Category	Score Range	Color	Description
HIGH	≥ 10	🟢	Strong pharmaceutical potential
MEDIUM	5–9	🟡	Moderate pharmaceutical potential
LOW	1–4	🔴	Limited pharmaceutical potential
NONE	0	⚪	No detected pharmaceutical potential
Group Exclusivity Validation
AnnoMining automatically validates and corrects group assignments:
•	Priority: HIGH > MEDIUM > LOW
•	Genes found in multiple groups are reassigned to the highest priority group
•	Validation report printed to console showing overlap statistics
Gene Importance Metric (Network Analysis)
For network node selection using the enrichment method:
text
Importance = (Pharma Score × 2) + (KO Count × 3) + (Pathway Count × 1) + (PFAM Count × 0.5)
Species Comparison Metrics
Similarity Score (Species Network):
text
Similarity = 1 / (1 + ΔHIGH × 0.5 + ΔPercent × 0.3 + ΔClasses × 0.2)
Where:
•	ΔHIGH = difference in HIGH-potential gene counts
•	ΔPercent = difference in percentage of pharma genes relative to secondary metabolism
•	ΔClasses = difference in compound class diversity
PFAM Enrichment Statistics
Enrichment analysis uses Fisher's exact test with FDR correction:
•	Fold Enrichment: Frequency in group / Frequency in background
•	p-value: Statistical significance from Fisher's exact test
•	p.adjust: FDR-corrected p-value
•	Significant: p.adjust < 0.05
________________________________________
Troubleshooting
Common Issues
1. Tkinter Not Available
Error: ImportError: No module named tkinter
Solution: Install tkinter:
•	Ubuntu/Debian: sudo apt-get install python3-tk
•	CentOS/RHEL: sudo yum install python3-tkinter
•	macOS: Tkinter is typically included with Python.org installation
•	Windows: Tkinter is included with standard Python installation
2. File Format Issues
Symptom: Parsing produces few or no results
Solution:
•	Verify eggNOG file is in the correct format (tab-delimited)
•	Check that KEGG pathway columns contain mapXXXXX format
•	Ensure file contains at least one headerless line with data
•	Try the command-line converter first to verify extraction
3. Memory Errors
Symptom: Program crashes with memory errors on large genomes
Solution:
•	Close other applications to free RAM
•	Run analysis in segments using command-line options
•	Increase system swap space
4. Missing Visualizations
Symptom: Some plots are not generated
Solution:
•	These are skipped when insufficient data exists (e.g., <2 genes)
•	Check that genes have associated pathways
•	Verify that input files contain the required data
•	Ensure at least 2 genes are present in the group being analyzed
5. InterPro File Issues
Symptom: InterPro data not appearing in outputs
Solution:
•	Verify file is in TSV format
•	Check column 1 contains gene IDs matching eggNOG file
•	Ensure PFAM domains are in PFXXXXX format
•	GO terms should be in GO:XXXXXXX format
6. Text Overlap in Figures
Symptom: Text overlapping in class distribution plots
Solution (automatic in v2.2):
•	Legends are now positioned outside plots using bbox_to_anchor
•	Figure sizes adjust dynamically based on the number of items
•	Labels include counts to reduce crowding
•	If issues persist, increase figure size by modifying figsize parameters
7. Inconsistent Group Counts
Symptom: Pharma genes > Secondary genes in statistics
Solution (automatic in v2.2):
•	The pipeline validates that pharmaceutical potential is a subset of secondary metabolism
•	Automatic correction of group assignments with priority hierarchy
•	Consistency checks printed to console
•	If issues persist, verify that PHARMA_DETAILS is a subset of SECONDARY_PATHWAYS
Command-Line Usage
The GUI interface is the primary mode of operation, but the analysis functions can also be called from the command line:
python
# Secondary metabolism analysis
from annomining import run_secondary_analysis
run_secondary_analysis('eggnog.emapper.annotations', 'interpro.tsv', 'output_dir', 'prefix')

# Pharmaceutical analysis
from annomining import run_pharma_analysis
run_pharma_analysis('eggnog.emapper.annotations', 'interpro.tsv', 'pharma_out', 'pharma', threshold=5)

# Multi-genome analysis
from annomining import process_multiple_genomes
process_multiple_genomes(['sp1.emapper', 'sp2.emapper'], ['sp1.tsv', 'sp2.tsv'], ['Species1', 'Species2'])
Logging and Debugging
All modules generate detailed console output including:
•	Start time with timestamps
•	Progress indicators
•	Validation reports for group exclusivity
•	Warning messages for missing data
•	Error messages with stack traces
•	Completion messages with output locations
________________________________________
Citation
If you use AnnoMining in your research, please cite:
[Author Name] et al. (2024). AnnoMining: Integrated Platform for Secondary Metabolism Analysis. Version 2.2. GitHub Repository. https://github.com/yourusername/annomining
________________________________________
License
AnnoMining is provided under the MIT License.
________________________________________
Contact and Support
•	GitHub Issues: https://github.com/yourusername/annomining/issues
•	Email: valdir.stefenon@ufsc.br
________________________________________

