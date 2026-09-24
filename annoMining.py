#!/usr/bin/env python3
"""
annoMining - Secondary Metabolism, Pharmaceutical Potential and Disease Resistance Analysis
Version: 3.5
"""

import sys
import os
import re
import argparse
import subprocess
import tempfile
from collections import defaultdict, Counter
from datetime import datetime
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch, Circle
from itertools import combinations
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.multitest import multipletests

plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'DejaVu Sans'

try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext
except ImportError:
    print("Erro: tkinter não disponível. Instale python3-tk.")
    sys.exit(1)

# ============================================================================
# DICIONÁRIOS UNIFICADOS
# ============================================================================
SECONDARY_PATHWAYS = {
    'map00900': {'name': 'Terpenoid backbone', 'class': 'Terpenoids', 'color': '#2E86AB', 'pharma': True},
    'map00902': {'name': 'Monoterpenoid', 'class': 'Terpenoids', 'color': '#2E86AB', 'pharma': True},
    'map00904': {'name': 'Diterpenoid', 'class': 'Terpenoids', 'color': '#2E86AB', 'pharma': True},
    'map00906': {'name': 'Carotenoid', 'class': 'Terpenoids', 'color': '#2E86AB', 'pharma': True},
    'map00909': {'name': 'Sesquiterpenoid/Triterpenoid', 'class': 'Terpenoids', 'color': '#2E86AB', 'pharma': True},
    'map00940': {'name': 'Phenylpropanoid', 'class': 'Phenylpropanoids', 'color': '#A23B72', 'pharma': True},
    'map00941': {'name': 'Flavonoid', 'class': 'Phenylpropanoids', 'color': '#A23B72', 'pharma': True},
    'map00942': {'name': 'Anthocyanin', 'class': 'Phenylpropanoids', 'color': '#A23B72', 'pharma': True},
    'map00943': {'name': 'Isoflavonoid', 'class': 'Phenylpropanoids', 'color': '#A23B72', 'pharma': True},
    'map00944': {'name': 'Flavone/Flavonol', 'class': 'Phenylpropanoids', 'color': '#A23B72', 'pharma': True},
    'map00945': {'name': 'Stilbenoid/Curcuminoid', 'class': 'Phenylpropanoids', 'color': '#A23B72', 'pharma': True},
    'map00950': {'name': 'Isoquinoline alkaloid', 'class': 'Alkaloids', 'color': '#F18F01', 'pharma': True},
    'map00960': {'name': 'Tropane/Piperidine', 'class': 'Alkaloids', 'color': '#F18F01', 'pharma': True},
    'map00965': {'name': 'Betalain', 'class': 'Alkaloids', 'color': '#F18F01', 'pharma': True},
    'map00966': {'name': 'Glucosinolate', 'class': 'Glucosinolates', 'color': '#9B5DE5', 'pharma': True},
    'map00903': {'name': 'Cannabinoid biosynthesis', 'class': 'Cannabinoids', 'color': '#FF006E', 'pharma': True},
    'map00232': {'name': 'Caffeine metabolism', 'class': 'Xanthines', 'color': '#F15BB5', 'pharma': True},
}

SECONDARY_DOMAINS_HMM = {
    'Terpene_synthase': {
        'pfam': 'PF01397',
        'hmm': 'Terpene_synthase.hmm',
        'class': 'Terpenoids',
        'evalue': 1e-10,
        'description': 'Terpene synthase family'
    },
    'Terpene_synthase_N': {
        'pfam': 'PF03936',
        'hmm': 'Terpene_synthase_N.hmm',
        'class': 'Terpenoids',
        'evalue': 1e-10,
        'description': 'Terpene synthase N-terminal domain'
    },
    'Cytochrome_P450': {
        'pfam': 'PF00067',
        'hmm': 'Cytochrome_P450.hmm',
        'class': 'Terpenoids',
        'evalue': 1e-10,
        'description': 'Cytochrome P450 family'
    },
    'Chalcone_synthase': {
        'pfam': 'PF00195',
        'hmm': 'Chalcone_synthase.hmm',
        'class': 'Phenylpropanoids',
        'evalue': 1e-10,
        'description': 'Chalcone/stilbene synthase family'
    },
    'O_methyltransferase': {
        'pfam': 'PF00891',
        'hmm': 'O_methyltransferase.hmm',
        'class': 'Alkaloids',
        'evalue': 1e-10,
        'description': 'O-methyltransferase domain'
    },
    'SAM_methyltransferase': {
        'pfam': 'PF01596',
        'hmm': 'SAM_methyltransferase.hmm',
        'class': 'Alkaloids',
        'evalue': 1e-10,
        'description': 'SAM-dependent methyltransferase'
    },
    'Cyclase': {
        'pfam': 'PF00248',
        'hmm': 'Cyclase.hmm',
        'class': 'Terpenoids',
        'evalue': 1e-10,
        'description': 'Aldo/keto reductase family'
    },
    'Prenyltransferase': {
        'pfam': 'PF00494',
        'hmm': 'Prenyltransferase.hmm',
        'class': 'Terpenoids',
        'evalue': 1e-10,
        'description': 'Prenyltransferase family'
    },
    'Flavonoid_3_hydroxylase': {
        'pfam': 'PF02797',
        'hmm': 'Flavonoid_3_hydroxylase.hmm',
        'class': 'Flavonoids',
        'evalue': 1e-10,
        'description': 'Flavonoid 3-hydroxylase'
    },
    'Alkaloid_synthase': {
        'pfam': 'PF00201',
        'hmm': 'Alkaloid_synthase.hmm',
        'class': 'Alkaloids',
        'evalue': 1e-10,
        'description': 'Isoquinoline alkaloid synthase'
    },
}

# Domínios HMM específicos para Pharma (expandido)
PHARMA_DOMAINS_HMM = {
    'Terpene_synthase': {
        'pfam': 'PF01397',
        'hmm': 'Terpene_synthase.hmm',
        'class': 'Terpenoids',
        'evalue': 1e-10,
        'pharma_category': 'Terpenoids'
    },
    'Cytochrome_P450': {
        'pfam': 'PF00067',
        'hmm': 'Cytochrome_P450.hmm',
        'class': 'Terpenoids',
        'evalue': 1e-10,
        'pharma_category': 'Terpenoids'
    },
    'Chalcone_synthase': {
        'pfam': 'PF00195',
        'hmm': 'Chalcone_synthase.hmm',
        'class': 'Phenylpropanoids',
        'evalue': 1e-10,
        'pharma_category': 'Flavonoids'
    },
    'O_methyltransferase': {
        'pfam': 'PF00891',
        'hmm': 'O_methyltransferase.hmm',
        'class': 'Alkaloids',
        'evalue': 1e-10,
        'pharma_category': 'Alkaloids'
    },
    'SAM_methyltransferase': {
        'pfam': 'PF01596',
        'hmm': 'SAM_methyltransferase.hmm',
        'class': 'Alkaloids',
        'evalue': 1e-10,
        'pharma_category': 'Alkaloids'
    },
    'Prenyltransferase': {
        'pfam': 'PF00494',
        'hmm': 'Prenyltransferase.hmm',
        'class': 'Terpenoids',
        'evalue': 1e-10,
        'pharma_category': 'Terpenoids'
    },
    'Flavonoid_3_hydroxylase': {
        'pfam': 'PF02797',
        'hmm': 'Flavonoid_3_hydroxylase.hmm',
        'class': 'Flavonoids',
        'evalue': 1e-10,
        'pharma_category': 'Flavonoids'
    },
    'Alkaloid_synthase': {
        'pfam': 'PF00201',
        'hmm': 'Alkaloid_synthase.hmm',
        'class': 'Alkaloids',
        'evalue': 1e-10,
        'pharma_category': 'Alkaloids'
    },
    'Carotenoid_biosynthesis': {
        'pfam': 'PF00514',
        'hmm': 'Carotenoid_biosynthesis.hmm',
        'class': 'Carotenoids',
        'evalue': 1e-10,
        'pharma_category': 'Carotenoids'
    }
}

# Domínios HMM específicos para Resistance
RESISTANCE_DOMAINS_HMM = {
    'NB_ARC': {
        'pfam': 'PF00931',
        'hmm': 'NB_ARC.hmm',
        'class': 'NLR',
        'evalue': 1e-10,
        'resistance_class': 'NB_ARC'
    },
    'TIR': {
        'pfam': 'PF01582',
        'hmm': 'TIR.hmm',
        'class': 'TIR_NLR',
        'evalue': 1e-10,
        'resistance_class': 'TIR'
    },
    'LRR': {
        'pfam': 'PF00560',
        'hmm': 'LRR.hmm',
        'class': 'NLR',
        'evalue': 1e-10,
        'resistance_class': 'LRR_repeat'
    },
    'Jacalin': {
        'pfam': 'PF01419',
        'hmm': 'Jacalin.hmm',
        'class': 'Jacalin',
        'evalue': 1e-10,
        'resistance_class': 'Jacalin'
    },
    'CC_domain': {
        'pfam': 'PF00560',
        'hmm': 'CC_domain.hmm',
        'class': 'CC_NLR',
        'evalue': 1e-10,
        'resistance_class': 'CC'
    },
    'RPW8': {
        'pfam': 'PF05659',
        'hmm': 'RPW8.hmm',
        'class': 'RPW8',
        'evalue': 1e-10,
        'resistance_class': 'RPW8'
    }
}

PHARMA_DETAILS = {
    'map00950': {'name': 'Isoquinoline alkaloid', 'class': 'Alkaloids',
                 'compounds': 'morphine, codeine, berberine', 'pharma_use': 'Analgesic, antimicrobial',
                 'color': '#E63946'},
    'map00960': {'name': 'Tropane/Piperidine', 'class': 'Alkaloids',
                 'compounds': 'atropine, cocaine, nicotine', 'pharma_use': 'Anticholinergic',
                 'color': '#E63946'},
    'map00965': {'name': 'Betalain', 'class': 'Alkaloids',
                 'compounds': 'betanin, indicaxanthin', 'pharma_use': 'Antioxidant',
                 'color': '#E63946'},
    'map00900': {'name': 'Terpenoid backbone', 'class': 'Terpenoids',
                 'compounds': 'precursors for all terpenes', 'pharma_use': 'Diverse bioactivities',
                 'color': '#2A9D8F'},
    'map00909': {'name': 'Sesquiterpenoid/Triterpenoid', 'class': 'Terpenoids',
                 'compounds': 'artemisinin, taxol, ginsenosides', 'pharma_use': 'Antimalarial, anticancer',
                 'color': '#2A9D8F'},
    'map00904': {'name': 'Diterpenoid', 'class': 'Terpenoids',
                 'compounds': 'tanshinone, forskolin', 'pharma_use': 'Cardiovascular',
                 'color': '#2A9D8F'},
    'map00906': {'name': 'Carotenoid', 'class': 'Terpenoids',
                 'compounds': 'β-carotene, lycopene', 'pharma_use': 'Antioxidant',
                 'color': '#2A9D8F'},
    'map00940': {'name': 'Phenylpropanoid', 'class': 'Phenylpropanoids',
                 'compounds': 'curcumin, resveratrol', 'pharma_use': 'Antioxidant',
                 'color': '#E9C46A'},
    'map00941': {'name': 'Flavonoid', 'class': 'Flavonoids',
                 'compounds': 'quercetin, kaempferol, catechins', 'pharma_use': 'Cardioprotective',
                 'color': '#E9C46A'},
    'map00942': {'name': 'Anthocyanin', 'class': 'Flavonoids',
                 'compounds': 'cyanidin, delphinidin', 'pharma_use': 'Antioxidant',
                 'color': '#E9C46A'},
    'map00944': {'name': 'Flavone/Flavonol', 'class': 'Flavonoids',
                 'compounds': 'apigenin, luteolin', 'pharma_use': 'Anti-inflammatory',
                 'color': '#E9C46A'},
    'map00945': {'name': 'Stilbenoid/Curcuminoid', 'class': 'Stilbenoids',
                 'compounds': 'resveratrol, curcumin', 'pharma_use': 'Chemopreventive',
                 'color': '#E9C46A'},
    'map00966': {'name': 'Glucosinolate', 'class': 'Glucosinolates',
                 'compounds': 'sulforaphane', 'pharma_use': 'Anticancer',
                 'color': '#9B5DE5'},
    'map00232': {'name': 'Caffeine metabolism', 'class': 'Xanthines',
                 'compounds': 'caffeine', 'pharma_use': 'Stimulant',
                 'color': '#F15BB5'},
    'map00903': {'name': 'Cannabinoid biosynthesis', 'class': 'Cannabinoids',
                 'compounds': 'THC, CBD', 'pharma_use': 'Analgesic, anti-epileptic',
                 'color': '#FF006E'},
}

# ============================================================================
# DICIONÁRIOS PARA RESISTÊNCIA A DOENÇAS COM SUPORTE A TNJ
# ============================================================================

DISEASE_RESISTANCE_PFAMS = {
    'NLR': {
        'pfams': ['PF00931', 'PF01582', 'PF00560', 'PF07723'],
        'domains': ['NB-ARC', 'LRR', 'TIR', 'CC'],
        'description': 'NLR resistance genes (NB-ARC + LRR)',
        'color': '#E63946'
    },
    'TIR_NLR': {
        'pfams': ['PF00931', 'PF01582', 'PF13676', 'PF00560'],
        'domains': ['TIR', 'NB-ARC', 'LRR'],
        'description': 'TIR-NLR resistance genes',
        'color': '#F4A261'
    },
    'TNJ': {
        'pfams': ['PF01582', 'PF00931', 'PF01419'],
        'domains': ['TIR', 'NB-ARC', 'Jacalin'],
        'description': 'TIR-NBS-Jacalin (TNJ) resistance genes - Myrtaceae-specific',
        'color': '#2A9D8F',
        'requires': {
            'mandatory': ['PF01582', 'PF00931'],
            'optional': ['PF01419']
        }
    },
    'TNJ_like': {
        'pfams': ['PF01582', 'PF00931', 'PF01419', 'PF00560'],
        'domains': ['TIR', 'NB-ARC', 'Jacalin', 'LRR'],
        'description': 'TNJ-like genes (TIR-NBS-Jacalin-LRR)',
        'color': '#264653'
    },
    'CC_NLR': {
        'pfams': ['PF00931', 'PF00560'],
        'domains': ['CC', 'NB-ARC', 'LRR'],
        'description': 'CC-NLR resistance genes',
        'color': '#E9C46A'
    },
    'RPW8': {
        'pfams': ['PF05659', 'PF05660'],
        'domains': ['RPW8', 'CC'],
        'description': 'RPW8-type resistance genes',
        'color': '#9B5DE5'
    },
    'Jacalin': {
        'pfams': ['PF01419'],
        'domains': ['Jacalin'],
        'description': 'Jacalin lectin domain',
        'color': '#F15BB5'
    },
    'LRR_repeat': {
        'pfams': ['PF00560', 'PF07723', 'PF13855'],
        'domains': ['LRR', 'LRR_repeat'],
        'description': 'Leucine-rich repeat domains',
        'color': '#FF006E'
    },
    'NB_ARC': {
        'pfams': ['PF00931'],
        'domains': ['NB-ARC'],
        'description': 'NB-ARC domain (nucleotide binding)',
        'color': '#00BBF9'
    },
    'TIR': {
        'pfams': ['PF01582', 'PF13676'],
        'domains': ['TIR'],
        'description': 'TIR domain (Toll/Interleukin-1 receptor)',
        'color': '#00F5D4'
    },
    'CC': {
        'pfams': ['PF00560'],
        'domains': ['CC'],
        'description': 'Coiled-coil domain',
        'color': '#FEE440'
    },
    'NLR_Jacalin': {
        'pfams': ['PF00931', 'PF01419'],
        'domains': ['NB-ARC', 'Jacalin'],
        'description': 'NLR with Jacalin domain (no TIR)',
        'color': '#8338EC'
    },
    'TNL': {
        'pfams': ['PF01582', 'PF00931', 'PF00560'],
        'domains': ['TIR', 'NB-ARC', 'LRR'],
        'description': 'TNL resistance genes (TIR-NBS-LRR)',
        'color': '#FB8500'
    }
}

DISEASE_RESISTANCE_PATHWAYS = {
    'plant_pathogen_interaction': {
        'name': 'Plant-pathogen interaction',
        'kegg': 'map04626',
        'description': 'Plant defense signaling pathway',
        'color': '#E63946'
    },
    'MAPK_signaling_plant': {
        'name': 'MAPK signaling pathway - plant',
        'kegg': 'map04016',
        'description': 'Mitogen-activated protein kinase cascade',
        'color': '#F4A261'
    },
    'hormone_signal_transduction': {
        'name': 'Plant hormone signal transduction',
        'kegg': 'map04075',
        'description': 'SA, JA, ET, ABA signaling',
        'color': '#2A9D8F'
    },
    'phenylpropanoid_biosynthesis': {
        'name': 'Phenylpropanoid biosynthesis',
        'kegg': 'map00940',
        'description': 'Lignin and phytoalexin biosynthesis',
        'color': '#E9C46A'
    },
    'flavonoid_biosynthesis': {
        'name': 'Flavonoid biosynthesis',
        'kegg': 'map00941',
        'description': 'Flavonoid phytoalexins',
        'color': '#9B5DE5'
    },
    'terpenoid_backbone': {
        'name': 'Terpenoid backbone biosynthesis',
        'kegg': 'map00900',
        'description': 'Terpenoid phytoalexins',
        'color': '#F15BB5'
    }
}

DISEASE_RESISTANCE_KOS = {
    'NLR_genes': ['K13475', 'K13476', 'K13477', 'K13478', 'K13479'],
    'MAPK_signaling': ['K02540', 'K04345', 'K04436', 'K04437', 'K04438', 'K04439', 'K04440'],
    'Plant_pathogen_interaction': ['K02607', 'K03081', 'K03120', 'K03244', 'K03245', 'K03246', 'K03247'],
    'Hormone_signaling': ['K06445', 'K06446', 'K06447', 'K06448', 'K06449', 'K06450'],
    'Defense_proteins': ['K13480', 'K13481', 'K13482', 'K13483', 'K13484']
}

PHARMA_KOS = {
    'Bioactive_flavonoids': ['K00660','K05275','K05276','K13065','K13066','K05265','K05266','K13083','K13084'],
    'Anticancer_terpenes': ['K15891','K15892','K15893','K00487','K00507','K00509','K00511','K00512'],
    'Nutraceutical_carotenoids': ['K00514','K00515','K06444','K06445','K06446','K06447','K06448','K06449','K06450','K06451'],
    'Resveratrol_biosynthesis': ['K13071','K13072','K00430','K00431'],
    'Therapeutic_alkaloids': ['K01799','K01800','K01900','K01901','K01902','K01903'],
}

CORE_PHARMA_DOMAINS = {
    'Flavonoids': ['PF00195', 'PF02797', 'PF05834'],
    'Terpenoids': ['PF01397', 'PF03936', 'PF00067', 'PF00494'],
    'Alkaloids': ['PF01596', 'PF00891', 'PF00201', 'PF00155', 'PF00141'],
    'Phenylpropanoids': ['PF00195', 'PF00141'],
    'Carotenoids': ['PF00494', 'PF00514'],
    'P450': ['PF00067'],
    'Methyltransferases': ['PF01596', 'PF00891', 'PF01799']
}

GENERIC_DOMAINS = [
    'PF00232', 'PF00221', 'PF02728', 'PF01179', 'PF01154',
    'PF01255', 'PF01315', 'PF01370', 'PF01398', 'PF01494',
    'PF01593', 'PF01715', 'PF01915', 'PF02431', 'PF02670',
    'PF02738', 'PF02798', 'PF02803', 'PF02900', 'PF03171',
    'PF03450', 'PF03754', 'PF05368', 'PF08436', 'PF08502',
    'PF08540', 'PF08544', 'PF09265', 'PF12142', 'PF12143',
    'PF13193', 'PF13243', 'PF13249', 'PF13292', 'PF13460',
    'PF13561', 'PF14226', 'PF14497', 'PF08100', 'PF13450',
    'PF00348', 'PF00368', 'PF00107', 'PF00108', 'PF00111',
    'PF00264', 'PF00282', 'PF00501', 'PF00682', 'PF00933',
    'PF00941', 'PF01315', 'PF01370', 'PF01565', 'PF02727',
    'PF02803', 'PF03055', 'PF03450', 'PF03754', 'PF08491'
]

# ============================================================================
# FUNÇÕES DE PARSING
# ============================================================================

def parse_eggnog_universal(eggnog_file):
    data = {
        'genes': {},
        'kegg_pathways': defaultdict(set),
        'ec_numbers': defaultdict(set),
        'ko_numbers': defaultdict(set),
        'go_terms': defaultdict(set),
        'pfam_domains': defaultdict(set),
        'cog_categories': Counter(),
    }
    patterns = {
        'ko': re.compile(r'(?:ko:)?(K\d+)', re.IGNORECASE),
        'pathway': re.compile(r'(?:ko)?(map\d{5})', re.IGNORECASE),
        'ec': re.compile(r'(\d+\.\d+\.\d+\.\d+)'),
        'go': re.compile(r'(GO:\d+)', re.IGNORECASE),
        'pfam': re.compile(r'(PF\d+)', re.IGNORECASE),
    }
    with open(eggnog_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
            gene = fields[0]
            description = ''
            for field in fields[1:8]:
                if field and field != '-' and not field.startswith('COG') and not field.startswith('KOG'):
                    if len(field) < 100 and not any(x in field for x in ['@','|']):
                        description = field
                        break
            cog_category = '-'
            for field in fields:
                if re.match(r'^[A-Z]+$', field) and len(field) <= 5:
                    cog_category = field
                    break
            data['genes'][gene] = {'description': description if description else '-', 'cog_category': cog_category}
            for c in cog_category:
                if c.isalpha():
                    data['cog_categories'][c] += 1
            for ko in patterns['ko'].findall(line):
                if ko.startswith('K'):
                    data['ko_numbers'][gene].add(ko)
            for path in patterns['pathway'].findall(line):
                if path.startswith('map'):
                    data['kegg_pathways'][gene].add(path)
            for ec in patterns['ec'].findall(line):
                data['ec_numbers'][gene].add(ec)
            for go in patterns['go'].findall(line):
                data['go_terms'][gene].add(go)
            for pfam in patterns['pfam'].findall(line):
                data['pfam_domains'][gene].add(pfam)
    return data

def parse_interpro_universal(interpro_file):
    data = {
        'genes': set(),
        'pfam_domains': defaultdict(set),
        'go_terms': defaultdict(set),
        'interpro_domains': defaultdict(set)
    }
    if not os.path.exists(interpro_file):
        return data
    patterns = {
        'pfam': re.compile(r'(PF\d+)', re.IGNORECASE),
        'go': re.compile(r'(GO:\d+)', re.IGNORECASE),
    }
    with open(interpro_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            gene = fields[0]
            data['genes'].add(gene)
            for field in fields:
                pfam_match = patterns['pfam'].search(field)
                if pfam_match:
                    data['pfam_domains'][gene].add(pfam_match.group(1))
            for go in patterns['go'].findall(line):
                data['go_terms'][gene].add(go)
            if len(fields) > 12 and fields[12] and fields[12] != '-':
                data['interpro_domains'][gene].add(fields[12])
    return data

def parse_gff_pfams(gff_file):
    """Extrai domínios PFAM do GFF (de attributes)."""
    pfam_data = defaultdict(set)
    if not os.path.exists(gff_file):
        return pfam_data
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            feature = fields[2]
            attrs_str = fields[8]
            
            # Procurar por PFAM na linha de atributos
            pfam_matches = re.findall(r'(PF\d+)', attrs_str, re.IGNORECASE)
            
            # Extrair gene ID
            gene_id = None
            for item in attrs_str.split(';'):
                if item.startswith('ID='):
                    gene_id = item[3:]
                    break
                elif item.startswith('gene_id='):
                    gene_id = item[8:]
                    break
            
            if gene_id and pfam_matches:
                for pfam in pfam_matches:
                    pfam_data[gene_id].add(pfam)
    
    return pfam_data

def find_hmm_file(hmm_file, hmm_dir):
    """Encontra arquivo HMM em vários diretórios possíveis."""
    if os.path.exists(hmm_file):
        return hmm_file
    if hmm_dir and os.path.exists(os.path.join(hmm_dir, hmm_file)):
        return os.path.join(hmm_dir, hmm_file)
    
    common_dirs = ['/usr/share/hmmer/profiles/', '/opt/hmmer/profiles/', 
                  os.path.expanduser('~/hmmer_profiles/'), os.path.expanduser('~/hmmer/hmmer_profiles/')]
    for d in common_dirs:
        if os.path.exists(os.path.join(d, hmm_file)):
            return os.path.join(d, hmm_file)
    return None

# ============================================================================
# FUNÇÃO PARA EXTRAIR PROTEÍNAS DO GFF
# ============================================================================

def extract_proteins_from_gff(gff_file, genome_file, output_fasta):
    if not genome_file or not os.path.exists(genome_file):
        print("  ⚠️ Arquivo do genoma não fornecido ou não encontrado.")
        return False
    
    if not gff_file or not os.path.exists(gff_file):
        print("  ⚠️ Arquivo GFF não encontrado.")
        return False
    
    try:
        subprocess.run(['gffread', '--version'], capture_output=True, check=True)
    except:
        print("  ⚠️ gffread não encontrado. Instale: conda install -c bioconda gffread")
        return False
    
    cmd = f"gffread {gff_file} -g {genome_file} -y {output_fasta}"
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ⚠️ Erro no gffread: {result.stderr[:200] if result.stderr else 'desconhecido'}")
            return False
        
        if os.path.exists(output_fasta) and os.path.getsize(output_fasta) > 0:
            print(f"  ✅ Proteínas extraídas: {output_fasta} ({os.path.getsize(output_fasta)} bytes)")
            return True
        else:
            print("  ⚠️ Nenhuma proteína extraída. Verifique se o GFF contém CDS.")
            return False
    except Exception as e:
        print(f"  ⚠️ Erro ao executar gffread: {e}")
        return False

# ============================================================================
# FUNÇÃO PARA EXECUTAR HMMER (GENÉRICA)
# ============================================================================

def run_hmmer_search_generic(protein_fasta, hmm_dir, output_dir, domains_dict, evalue=1e-10, prefix=''):
    """Executa HMMER genérico para qualquer conjunto de domínios."""
    hmmer_results = defaultdict(dict)
    
    try:
        subprocess.run(['hmmsearch', '-h'], capture_output=True, check=True)
    except:
        print("  ⚠️ hmmsearch não encontrado. Instale: conda install -c bioconda hmmer")
        return hmmer_results
    
    for domain_name, domain_info in domains_dict.items():
        hmm_file = domain_info.get('hmm')
        hmm_path = find_hmm_file(hmm_file, hmm_dir)
        
        if not hmm_path:
            print(f"  ⚠️ Arquivo HMM não encontrado para {domain_name}: {hmm_file}")
            continue
        
        output_file = os.path.join(output_dir, f"{prefix}_{domain_name}_hits.tsv")
        
        cmd = [
            'hmmsearch',
            '--domtblout', output_file,
            '--noali',
            '-E', str(evalue),
            hmm_path,
            protein_fasta
        ]
        
        try:
            subprocess.run(cmd, capture_output=True, check=True)
            
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                with open(output_file, 'r') as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        fields = line.strip().split()
                        if len(fields) >= 5:
                            gene = fields[0]
                            domain = fields[3]
                            evalue_hit = float(fields[5])
                            
                            if gene not in hmmer_results:
                                hmmer_results[gene] = {
                                    'domains': [],
                                    'classes': set(),
                                    'pfams': [],
                                    'max_evalue': float('inf')
                                }
                            
                            hmmer_results[gene]['domains'].append({
                                'domain_name': domain_name,
                                'domain': domain,
                                'evalue': evalue_hit
                            })
                            if 'class' in domain_info:
                                hmmer_results[gene]['classes'].add(domain_info['class'])
                            if 'pfam' in domain_info:
                                hmmer_results[gene]['pfams'].append(domain_info['pfam'])
                            hmmer_results[gene]['max_evalue'] = min(hmmer_results[gene]['max_evalue'], evalue_hit)
            
        except subprocess.CalledProcessError as e:
            print(f"  ⚠️ Erro ao executar hmmsearch para {domain_name}: {e}")
    
    return hmmer_results

# ============================================================================
# FUNÇÃO PARA INTEGRAR RESULTADOS HMMER (GENÉRICA)
# ============================================================================

def integrate_hmmer_results_generic(eggnog_data, interpro_data, hmmer_results, gff_genes, 
                                   domain_type='secondary', pharma_category=None, resistance_class=None):
    """
    Integra resultados HMMER aos dados existentes.
    domain_type: 'secondary', 'pharma', 'resistance'
    """
    merged_data = merge_annotation_data(eggnog_data, interpro_data, SECONDARY_PATHWAYS)
    
    genes_added = 0
    genes_updated = 0
    
    for gene, hmmer_info in hmmer_results.items():
        if gene in merged_data:
            # Atualizar gene existente
            if 'hmmer_domains' not in merged_data[gene]:
                merged_data[gene]['hmmer_domains'] = []
                merged_data[gene]['hmmer_classes'] = set()
            merged_data[gene]['hmmer_domains'].extend(hmmer_info['domains'])
            merged_data[gene]['hmmer_classes'].update(hmmer_info['classes'])
            
            # Adicionar PFAMs
            for pfam in hmmer_info['pfams']:
                merged_data[gene]['pfam_domains'].add(pfam)
            
            # Adicionar pathways secundários
            for class_name in hmmer_info['classes']:
                for pathway, p_info in SECONDARY_PATHWAYS.items():
                    if p_info['class'] == class_name:
                        merged_data[gene]['secondary_pathways'].add(pathway)
                        merged_data[gene]['secondary_class'].add(class_name)
                        break
            
            if 'source' not in merged_data[gene]:
                merged_data[gene]['source'] = []
            if 'HMMER' not in merged_data[gene]['source']:
                merged_data[gene]['source'].append('HMMER')
            genes_updated += 1
            
        elif gene in gff_genes or True:  # Adicionar mesmo que não esteja no GFF
            merged_data[gene] = {
                'source': ['HMMER'],
                'description': f'HMMER: {", ".join(hmmer_info["classes"])}',
                'cog_category': '-',
                'go_terms': set(),
                'pfam_domains': set(hmmer_info['pfams']),
                'kegg_pathways': set(),
                'ec_numbers': set(),
                'ko_numbers': set(),
                'secondary_pathways': set(),
                'secondary_class': set(),
                'hmmer_domains': hmmer_info['domains'],
                'hmmer_classes': hmmer_info['classes']
            }
            
            for class_name in hmmer_info['classes']:
                for pathway, p_info in SECONDARY_PATHWAYS.items():
                    if p_info['class'] == class_name:
                        merged_data[gene]['secondary_pathways'].add(pathway)
                        merged_data[gene]['secondary_class'].add(class_name)
                        break
            
            genes_added += 1
    
    if genes_added > 0:
        print(f"  ✅ {genes_added} genes adicionados via HMMER ({domain_type})")
    if genes_updated > 0:
        print(f"  ✅ {genes_updated} genes existentes atualizados via HMMER ({domain_type})")
    
    return merged_data

# ============================================================================
# FUNÇÕES DE MERGE - INDEPENDENTES PARA CADA ANÁLISE
# ============================================================================

def merge_annotation_data(eggnog_data, interpro_data, pathway_dict):
    merged = {}
    all_genes = set(eggnog_data['genes'].keys()) | interpro_data['genes']
    
    for gene in all_genes:
        merged[gene] = {
            'source': [],
            'description': eggnog_data['genes'].get(gene, {}).get('description', '-'),
            'cog_category': eggnog_data['genes'].get(gene, {}).get('cog_category', '-'),
            'go_terms': eggnog_data['go_terms'].get(gene, set()) | interpro_data['go_terms'].get(gene, set()),
            'pfam_domains': eggnog_data['pfam_domains'].get(gene, set()) | interpro_data['pfam_domains'].get(gene, set()),
            'kegg_pathways': eggnog_data['kegg_pathways'].get(gene, set()),
            'ec_numbers': eggnog_data['ec_numbers'].get(gene, set()),
            'ko_numbers': eggnog_data['ko_numbers'].get(gene, set()),
            'secondary_pathways': set(),
            'secondary_class': set(),
            'hmmer_domains': [],
            'hmmer_classes': set()
        }
        
        if gene in eggnog_data['genes']:
            merged[gene]['source'].append('eggNOG')
        if gene in interpro_data['genes']:
            merged[gene]['source'].append('InterProScan')
        
        for path in merged[gene]['kegg_pathways']:
            if path in pathway_dict:
                merged[gene]['secondary_pathways'].add(path)
                merged[gene]['secondary_class'].add(pathway_dict[path]['class'])
    
    return merged

# ============================================================================
# FUNÇÃO DE VALIDAÇÃO DE GRUPOS
# ============================================================================

def validate_pharma_groups(merged_data):
    high = set([g for g, info in merged_data.items() if info.get('pharma_potential') == 'HIGH'])
    medium = set([g for g, info in merged_data.items() if info.get('pharma_potential') == 'MEDIUM'])
    low = set([g for g, info in merged_data.items() if info.get('pharma_potential') == 'LOW'])
    
    print("\n🔍 VALIDAÇÃO DE GRUPOS FARMACOLÓGICOS:")
    print(f"   HIGH: {len(high)} genes")
    print(f"   MEDIUM: {len(medium)} genes")
    print(f"   LOW: {len(low)} genes")
    
    overlap_hm = high & medium
    overlap_hl = high & low
    overlap_ml = medium & low
    
    issues_found = False
    
    if overlap_hm:
        print(f"   ⚠️ SOBREPOSIÇÃO HIGH-MEDIUM: {len(overlap_hm)} genes")
        issues_found = True
    if overlap_hl:
        print(f"   ⚠️ SOBREPOSIÇÃO HIGH-LOW: {len(overlap_hl)} genes")
        issues_found = True
    if overlap_ml:
        print(f"   ⚠️ SOBREPOSIÇÃO MEDIUM-LOW: {len(overlap_ml)} genes")
        issues_found = True
    
    if not issues_found:
        print("   ✅ Grupos são mutuamente exclusivos")
    else:
        print("\n🔧 CORRIGINDO SOBREPOSIÇÕES (prioridade: HIGH > MEDIUM > LOW):")
        for gene in overlap_hm:
            merged_data[gene]['pharma_potential'] = 'HIGH'
        if overlap_hl:
            for gene in overlap_hl:
                merged_data[gene]['pharma_potential'] = 'HIGH'
        if overlap_ml:
            for gene in overlap_ml:
                merged_data[gene]['pharma_potential'] = 'MEDIUM'
    
    return {
        'high': high,
        'medium': medium,
        'low': low,
        'overlap_hm': overlap_hm,
        'overlap_hl': overlap_hl,
        'overlap_ml': overlap_ml,
        'issues_found': issues_found
    }

# ============================================================================
# FUNÇÕES DE ENRIQUECIMENTO
# ============================================================================

def calculate_enrichment(group_genes, background_genes, merged_data, min_count=2, min_fold=1.0, filter_generic=True):
    if not group_genes or not background_genes:
        return pd.DataFrame()
    
    group_pfam_counts = Counter()
    for gene in group_genes:
        if gene in merged_data:
            for pfam in merged_data[gene]['pfam_domains']:
                if filter_generic and pfam in GENERIC_DOMAINS:
                    continue
                group_pfam_counts[pfam] += 1
    
    bg_pfam_counts = Counter()
    for gene in background_genes:
        if gene in merged_data:
            for pfam in merged_data[gene]['pfam_domains']:
                if filter_generic and pfam in GENERIC_DOMAINS:
                    continue
                bg_pfam_counts[pfam] += 1
    
    total_group = len(group_genes)
    total_bg = len(background_genes)
    
    results = []
    all_pfams = set(group_pfam_counts.keys()) | set(bg_pfam_counts.keys())
    
    for pfam in all_pfams:
        count_in_group = group_pfam_counts.get(pfam, 0)
        count_in_bg = bg_pfam_counts.get(pfam, 0)
        
        if count_in_group < min_count:
            continue
        
        freq_group = count_in_group / total_group
        freq_bg = count_in_bg / total_bg if total_bg > 0 else 0
        
        if freq_bg == 0:
            fold_enrich = freq_group * 100
        else:
            fold_enrich = freq_group / freq_bg
        
        if fold_enrich < min_fold:
            continue
        
        table = [
            [count_in_group, total_group - count_in_group],
            [count_in_bg, total_bg - count_in_bg]
        ]
        
        try:
            oddsratio, p_fisher = fisher_exact(table, alternative='two-sided')
        except:
            oddsratio = np.nan
            p_fisher = 1.0
        
        try:
            chi2, p_chi2, dof, expected = chi2_contingency(table)
        except:
            p_chi2 = 1.0
        
        p_combined = min(p_fisher, p_chi2)
        
        pathway_category = 'Other'
        for category, domains in CORE_PHARMA_DOMAINS.items():
            if pfam in domains:
                pathway_category = category
                break
        
        results.append({
            'pfam_domain': pfam,
            'pathway_category': pathway_category,
            'count_in_group': count_in_group,
            'total_group': total_group,
            'freq_group': round(freq_group, 4),
            'count_in_background': count_in_bg,
            'total_background': total_bg,
            'freq_background': round(freq_bg, 4),
            'fold_enrichment': round(fold_enrich, 2),
            'odds_ratio': round(oddsratio, 2) if not np.isnan(oddsratio) else np.nan,
            'p_value': p_combined,
            'p_adjust': np.nan,
            'significant': p_combined < 0.05
        })
    
    if results:
        df = pd.DataFrame(results)
        if len(df) > 0:
            _, p_adjust, _, _ = multipletests(df['p_value'].values, method='fdr_bh')
            df['p_adjust'] = p_adjust
            df['significant'] = df['p_adjust'] < 0.05
        return df
    return pd.DataFrame()

# ============================================================================
# FUNÇÕES DE VISUALIZAÇÃO (mantidas do original)
# ============================================================================

def create_ko_dotplot(merged_data, output_prefix, top_n=15):
    ko_by_class = defaultdict(Counter)
    class_counts = Counter()
    for gene, info in merged_data.items():
        for pathway in info['secondary_pathways']:
            if pathway in SECONDARY_PATHWAYS:
                pc = SECONDARY_PATHWAYS[pathway]['class']
                class_counts[pc] += 1
                for ko in info['ko_numbers']:
                    ko_by_class[pc][ko] += 1
    all_data = []
    for class_name, ko_counter in ko_by_class.items():
        total_in_class = class_counts[class_name]
        for ko, count in ko_counter.most_common(top_n):
            freq_in_class = count / total_in_class if total_in_class > 0 else 0
            all_data.append({
                'class': class_name,
                'ko': ko,
                'count': count,
                'frequency': freq_in_class,
                'label': f"{ko} [{class_name}]"
            })
    if not all_data:
        return
    df = pd.DataFrame(all_data)
    df = df.sort_values('frequency', ascending=False).head(top_n)
    fig_height = max(6, len(df) * 0.4)
    fig, ax = plt.subplots(figsize=(14, fig_height))
    class_colors = {'Terpenoids': '#2E86AB', 'Phenylpropanoids': '#A23B72',
                    'Alkaloids': '#F18F01', 'Other': '#6A4E3A', 
                    'Glucosinolates': '#9B5DE5', 'Xanthines': '#F15BB5',
                    'Cannabinoids': '#FF006E'}
    colors = [class_colors.get(c, '#999999') for c in df['class']]
    unique_classes = df['class'].unique()
    legend_elements = [Patch(facecolor=class_colors.get(c, '#999999'), label=c)
                       for c in unique_classes if c in class_colors]
    ax.scatter(
        df['frequency'],
        df['ko'],
        s=df['count'] * 20 + 30,
        c=colors,
        alpha=0.8,
        edgecolors='black',
        linewidth=0.5
    )
    ax.set_xlabel('Frequency in class', fontsize=12, fontweight='bold')
    ax.set_ylabel('KEGG Ortholog', fontsize=12, fontweight='bold')
    ax.set_title('Top KOs by Secondary Metabolite Class', fontsize=14, fontweight='bold')
    ax.legend(handles=legend_elements,
             loc='center left',
             bbox_to_anchor=(1.02, 0.5),
             fontsize=10,
             title='Class',
             title_fontsize=11)
    plt.subplots_adjust(right=0.78)
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(f"{output_prefix}_ko_dotplot.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_ko_dotplot.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   ✅ KO dotplot salvo: {output_prefix}_ko_dotplot.pdf")

def plot_top_kos_by_category(merged_data, output_prefix):
    ko_by_class = defaultdict(Counter)
    for gene, info in merged_data.items():
        for pathway in info['secondary_pathways']:
            if pathway in SECONDARY_PATHWAYS:
                pc = SECONDARY_PATHWAYS[pathway]['class']
                for ko in info['ko_numbers']:
                    ko_by_class[pc][ko] += 1
    classes = [c for c in ko_by_class.keys() if sum(ko_by_class[c].values()) > 0]
    if not classes:
        return
    
    n_cols = min(2, len(classes))
    n_rows = (len(classes) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
    if n_rows == 1 and n_cols == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    for idx, class_name in enumerate(classes):
        if idx >= len(axes):
            break
        ax = axes[idx]
        top = ko_by_class[class_name].most_common(10)
        if not top:
            continue
        kos, counts = zip(*top)
        bars = ax.barh(range(len(kos)), counts, color='#2E86AB')
        ax.set_yticks(range(len(kos)))
        ax.set_yticklabels(kos, fontsize=9)
        ax.set_xlabel('Number of genes', fontsize=10)
        ax.set_title(f'{class_name} - Top KOs', fontsize=11)
        ax.invert_yaxis()
        for i, (bar, count) in enumerate(zip(bars, counts)):
            ax.text(count, bar.get_y()+bar.get_height()/2, f' {count}', va='center', fontsize=8)
    
    for idx in range(len(classes), len(axes)):
        axes[idx].set_visible(False)
    
    plt.suptitle('Most Abundant KEGG Orthologs by Secondary Metabolite Class', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"{output_prefix}_top_kos_by_class.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_top_kos_by_class.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_pathway_completeness(merged_data, output_prefix):
    pathway_counts = Counter()
    for gene, info in merged_data.items():
        for pathway in info['secondary_pathways']:
            pathway_counts[pathway] += 1
    if not pathway_counts:
        return
    
    pathways = list(pathway_counts.keys())
    counts = [pathway_counts[p] for p in pathways]
    names = [SECONDARY_PATHWAYS[p]['name'] for p in pathways]
    colors = [SECONDARY_PATHWAYS[p]['color'] for p in pathways]
    
    fig, ax = plt.subplots(figsize=(14, max(6, len(pathways)*0.4)))
    bars = ax.barh(range(len(pathways)), counts, color=colors, alpha=0.8)
    ax.set_yticks(range(len(pathways)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel('Number of genes', fontsize=12, fontweight='bold')
    ax.set_title('Genes Associated with Secondary Metabolism Pathways', fontsize=14, fontweight='bold')
    
    for bar, count in zip(bars, counts):
        ax.text(count + 0.5, bar.get_y() + bar.get_height()/2, str(count), 
                va='center', fontsize=9, fontweight='bold')
    
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_pathway_completeness.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_pathway_completeness.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_pie_chart_secondary_classes(merged_data, output_prefix):
    class_counts = Counter()
    for gene, info in merged_data.items():
        for c in info['secondary_class']:
            class_counts[c] += 1
    if not class_counts:
        return
    
    class_colors = {'Terpenoids':'#2E86AB', 'Phenylpropanoids':'#A23B72', 
                    'Alkaloids':'#F18F01', 'Other':'#6A4E3A',
                    'Glucosinolates':'#9B5DE5', 'Xanthines':'#F15BB5',
                    'Cannabinoids':'#FF006E'}
    colors = [class_colors.get(c, '#999999') for c in class_counts.keys()]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
    
    wedges1, texts1, autotexts1 = ax1.pie(
        class_counts.values(), 
        labels=class_counts.keys(), 
        colors=colors, 
        autopct='%1.1f%%', 
        startangle=90, 
        explode=[0.05]*len(class_counts),
        pctdistance=0.75,
        labeldistance=1.1
    )
    ax1.set_title('Distribution by Metabolite Class', fontsize=12, fontweight='bold')
    
    wedges2, texts2, autotexts2 = ax2.pie(
        class_counts.values(), 
        colors=colors, 
        autopct='%1.1f%%', 
        startangle=90, 
        wedgeprops=dict(width=0.3),
        pctdistance=0.7,
        labeldistance=1.1
    )
    ax2.set_title('Donut View', fontsize=12, fontweight='bold')
    ax2.legend(
        wedges2, 
        [f"{c} ({count})" for c, count in zip(class_counts.keys(), class_counts.values())],
        title="Classes", 
        loc="center left", 
        bbox_to_anchor=(1.15, 0.5), 
        fontsize=8,
        title_fontsize=9
    )
    
    plt.suptitle('Secondary Metabolism Gene Distribution', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 0.85, 0.95])
    plt.savefig(f"{output_prefix}_class_distribution.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_class_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_upset_venn(merged_data, output_prefix):
    classes = ['Terpenoids','Phenylpropanoids','Alkaloids','Other','Glucosinolates','Xanthines','Cannabinoids']
    classes = [c for c in classes if any(c in info['secondary_class'] for info in merged_data.values())]
    if len(classes) < 2:
        return
    
    matrix = []
    for gene, info in merged_data.items():
        row = [1 if c in info['secondary_class'] else 0 for c in classes]
        if sum(row) > 0:
            matrix.append(row)
    if len(matrix) < 2:
        return
    
    df_pres = pd.DataFrame(matrix, columns=classes)
    combo_counts = {}
    max_comb = min(4, len(classes))
    for i in range(1, max_comb + 1):
        for combo in combinations(classes, i):
            mask = pd.Series([True]*len(df_pres))
            for c in combo:
                mask &= (df_pres[c]==1)
            count = mask.sum()
            if count > 0:
                combo_counts[combo] = count
    
    sorted_combos = sorted(combo_counts.items(), key=lambda x: len(x[0]), reverse=False)
    combos_names = ['∩'.join(c[0]) for c in sorted_combos]
    counts = [c[1] for c in sorted_combos]
    
    fig, ax = plt.subplots(figsize=(max(10, len(combos_names)*0.6), 6))
    bars = ax.bar(range(len(combos_names)), counts, color='#2E86AB', alpha=0.7)
    ax.set_xticks(range(len(combos_names)))
    ax.set_xticklabels(combos_names, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Number of genes', fontsize=12, fontweight='bold')
    ax.set_title('Gene Overlap Between Metabolite Classes (UpSet-style)', fontsize=14, fontweight='bold')
    
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5, 
                str(count), ha='center', va='bottom', fontsize=9)
    
    ax.grid(True, axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_upset_venn.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_upset_venn.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_top_ec_numbers(merged_data, output_prefix):
    ec_counter = Counter()
    for gene, info in merged_data.items():
        for ec in info['ec_numbers']:
            ec_counter[ec] += 1
    if not ec_counter:
        return
    
    top = ec_counter.most_common(15)
    ec_names, ec_counts = zip(*top)
    
    fig, ax = plt.subplots(figsize=(14, max(6, len(ec_names)*0.4)))
    bars = ax.barh(range(len(ec_names)), ec_counts, color='#A23B72', alpha=0.8)
    ax.set_yticks(range(len(ec_names)))
    ax.set_yticklabels(ec_names, fontsize=8)
    ax.set_xlabel('Number of genes', fontsize=12, fontweight='bold')
    ax.set_title('Most Abundant EC Numbers in Secondary Metabolism', fontsize=14, fontweight='bold')
    
    for bar, count in zip(bars, ec_counts):
        ax.text(count + 0.5, bar.get_y() + bar.get_height()/2, str(count), 
                va='center', fontsize=9, fontweight='bold')
    
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_top_ec_numbers.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_top_ec_numbers.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_top_pfam_domains(merged_data, output_prefix):
    pfam_counter = Counter()
    for gene, info in merged_data.items():
        for pfam in info['pfam_domains']:
            if pfam in GENERIC_DOMAINS:
                continue
            pfam_counter[pfam] += 1
    if not pfam_counter:
        return
    
    top = pfam_counter.most_common(15)
    pfam_names, pfam_counts = zip(*top)
    
    fig, ax = plt.subplots(figsize=(14, max(6, len(pfam_names)*0.4)))
    bars = ax.barh(range(len(pfam_names)), pfam_counts, color='#F18F01', alpha=0.8)
    ax.set_yticks(range(len(pfam_names)))
    ax.set_yticklabels(pfam_names, fontsize=9)
    ax.set_xlabel('Number of genes', fontsize=12, fontweight='bold')
    ax.set_title('Most Abundant Core PFAM Domains in Secondary Metabolism', fontsize=14, fontweight='bold')
    
    for bar, count in zip(bars, pfam_counts):
        ax.text(count + 0.5, bar.get_y() + bar.get_height()/2, str(count), 
                va='center', fontsize=9, fontweight='bold')
    
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_top_core_pfam_domains.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_top_core_pfam_domains.png", dpi=300, bbox_inches='tight')
    plt.close()

def generate_publication_table_secondary(merged_data, output_prefix):
    rows = []
    for gene, info in merged_data.items():
        if info['secondary_pathways']:
            rows.append({
                'Gene_ID': gene,
                'Source': ', '.join(info['source']),
                'Description': info['description'][:100],
                'COG_Category': info['cog_category'],
                'Secondary_Classes': ', '.join(info['secondary_class']),
                'Secondary_Pathways': ', '.join([f"{p} ({SECONDARY_PATHWAYS[p]['name']})" for p in info['secondary_pathways']]),
                'KEGG_Orthologs': ', '.join(list(info['ko_numbers'])[:10]),
                'EC_Numbers': ', '.join(list(info['ec_numbers'])[:10]),
                'PFAM_Domains': ', '.join(list(info['pfam_domains'])[:10]),
                'GO_Terms_Count': len(info['go_terms'])
            })
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(f"{output_prefix}_secondary_genes_table.csv", index=False)
        try:
            df.to_excel(f"{output_prefix}_secondary_genes_table.xlsx", index=False)
        except:
            pass
        with open(f"{output_prefix}_statistics.txt", 'w') as f:
            f.write("SECONDARY METABOLISM ANALYSIS - SUMMARY\n")
            f.write(f"Total genes: {len(merged_data)}\n")
            f.write(f"Secondary genes: {len(rows)}\n")
            f.write(f"Percentage: {len(rows)/len(merged_data)*100:.1f}%\n")
            class_summary = Counter()
            for row in rows:
                for c in row['Secondary_Classes'].split(', '):
                    class_summary[c] += 1
            f.write("\nGenes by class:\n")
            for c, cnt in class_summary.most_common():
                f.write(f"  {c}: {cnt}\n")

# ============================================================================
# FUNÇÕES DE ANÁLISE FARMACOLÓGICA - INDEPENDENTE
# ============================================================================

def calculate_pharma_potential(merged_data, pharma_threshold=4):
    pharma_scores = {}
    
    for gene, info in merged_data.items():
        if not info['secondary_pathways']:
            pharma_scores[gene] = {
                'score': 0,
                'potential': 'NONE',
                'pharma_pathways': [],
                'compound_classes': [],
                'pharma_uses': []
            }
            continue
        
        score = 0
        pathways_found = []
        compound_classes = set()
        pharma_uses = set()
        
        for pathway in info['secondary_pathways']:
            if pathway in PHARMA_DETAILS:
                score += 2
                pathways_found.append(PHARMA_DETAILS[pathway]['name'])
                compound_classes.add(PHARMA_DETAILS[pathway]['class'])
                pharma_uses.add(PHARMA_DETAILS[pathway]['pharma_use'])
        
        if score == 0:
            pharma_scores[gene] = {
                'score': 0,
                'potential': 'NONE',
                'pharma_pathways': [],
                'compound_classes': [],
                'pharma_uses': []
            }
            continue
        
        for ko in info['ko_numbers']:
            for category, ko_list in PHARMA_KOS.items():
                if ko in ko_list:
                    score += 3
                    break
        
        for pfam in info['pfam_domains']:
            if pfam in GENERIC_DOMAINS:
                continue
            for category, domain_list in CORE_PHARMA_DOMAINS.items():
                if pfam in domain_list:
                    score += 2
                    break
        
        if score >= 10:
            potential = "HIGH"
        elif score >= pharma_threshold:
            potential = "MEDIUM"
        else:
            potential = "LOW"
        
        pharma_scores[gene] = {
            'score': score,
            'potential': potential,
            'pharma_pathways': list(set(pathways_found)),
            'compound_classes': list(compound_classes),
            'pharma_uses': list(pharma_uses)
        }
    
    return pharma_scores

def merge_pharma_data(eggnog_data, interpro_data, pharma_threshold=4):
    merged_temp = merge_annotation_data(eggnog_data, interpro_data, SECONDARY_PATHWAYS)
    pharma_scores = calculate_pharma_potential(merged_temp, pharma_threshold)
    
    merged = {}
    for gene in merged_temp:
        merged[gene] = merged_temp[gene].copy()
        merged[gene]['pharma_potential'] = pharma_scores.get(gene, {}).get('potential', 'NONE')
        merged[gene]['pharma_score'] = pharma_scores.get(gene, {}).get('score', 0)
        merged[gene]['pharma_pathways'] = pharma_scores.get(gene, {}).get('pharma_pathways', [])
        merged[gene]['compound_classes'] = pharma_scores.get(gene, {}).get('compound_classes', [])
        merged[gene]['pharma_uses'] = pharma_scores.get(gene, {}).get('pharma_uses', [])
    
    return merged

# ============================================================================
# FUNÇÕES DE ANÁLISE DE RESISTÊNCIA - INDEPENDENTE
# ============================================================================

def calculate_disease_resistance_score(merged_data, resistance_threshold=4):
    scores = {}
    
    for gene, info in merged_data.items():
        score = 0
        resistance_classes = set()
        resistance_pathways = set()
        resistance_pfams = set()
        
        gene_pfams = set(info['pfam_domains'])
        
        has_tir = 'PF01582' in gene_pfams or 'PF13676' in gene_pfams
        has_nb_arc = 'PF00931' in gene_pfams
        has_jacalin = 'PF01419' in gene_pfams
        has_lrr = any(pfam in gene_pfams for pfam in ['PF00560', 'PF07723', 'PF13855'])
        has_cc = 'PF00560' in gene_pfams
        
        if has_tir and has_nb_arc and has_jacalin and not has_lrr:
            resistance_classes.add('TNJ')
            score += 10
            resistance_pfams.update(['PF01582', 'PF00931', 'PF01419'])
            info['tnj_detected'] = True
            info['tnj_architecture'] = 'TIR-NB-ARC-Jacalin (TNJ canônico)'
            info['tnj_confidence'] = 'HIGH'
        
        elif has_tir and has_nb_arc and has_jacalin and has_lrr:
            resistance_classes.add('TNJ_like')
            score += 8
            resistance_pfams.update(['PF01582', 'PF00931', 'PF01419'])
            info['tnj_detected'] = True
            info['tnj_architecture'] = 'TIR-NB-ARC-Jacalin-LRR (TNJ-like)'
            info['tnj_confidence'] = 'MEDIUM'
        
        elif has_nb_arc and has_jacalin and not has_tir:
            resistance_classes.add('NLR_Jacalin')
            score += 6
            resistance_pfams.update(['PF00931', 'PF01419'])
            info['tnj_detected'] = False
            info['tnj_architecture'] = 'NB-ARC-Jacalin (NLR com Jacalin)'
        
        elif has_tir and has_nb_arc and has_lrr and not has_jacalin:
            resistance_classes.add('TNL')
            score += 8
            resistance_pfams.update(['PF01582', 'PF00931', 'PF00560'])
            info['tnj_detected'] = False
            info['tnj_architecture'] = 'TIR-NB-ARC-LRR (TNL clássico)'
        
        elif has_cc and has_nb_arc and has_lrr and not has_tir:
            resistance_classes.add('CC_NLR')
            score += 7
            resistance_pfams.update(['PF00560', 'PF00931'])
        
        for pfam in gene_pfams:
            for class_name, class_info in DISEASE_RESISTANCE_PFAMS.items():
                if pfam in class_info['pfams']:
                    if class_name not in ['TNJ', 'TNJ_like', 'NLR_Jacalin', 'TNL', 'CC_NLR']:
                        resistance_classes.add(class_name)
                        resistance_pfams.add(pfam)
                        if class_name in ['NLR', 'TIR_NLR']:
                            score += 3
                        elif class_name in ['LRR_repeat', 'NB_ARC', 'TIR', 'CC']:
                            score += 2
                        else:
                            score += 1
                        break
        
        for ko in info['ko_numbers']:
            for category, ko_list in DISEASE_RESISTANCE_KOS.items():
                if ko in ko_list:
                    score += 2
                    break
        
        for pathway in info['kegg_pathways']:
            for rp, rp_info in DISEASE_RESISTANCE_PATHWAYS.items():
                if pathway == rp_info.get('kegg', ''):
                    score += 2
                    resistance_pathways.add(rp)
                    break
        
        if 'NLR' in resistance_classes and 'Jacalin' in resistance_classes:
            if 'TNJ' not in resistance_classes and 'NLR_Jacalin' not in resistance_classes:
                score += 3
                resistance_classes.add('NLR_Jacalin_hybrid')
        
        if len(resistance_classes) >= 3:
            score += 2
        
        if score >= 10:
            potential = "HIGH"
        elif score >= resistance_threshold:
            potential = "MEDIUM"
        elif score > 0:
            potential = "LOW"
        else:
            potential = "NONE"
        
        scores[gene] = {
            'score': score,
            'potential': potential,
            'resistance_classes': list(resistance_classes),
            'resistance_pathways': list(resistance_pathways),
            'resistance_pfams': list(resistance_pfams),
            'tnj_detected': info.get('tnj_detected', False),
            'tnj_architecture': info.get('tnj_architecture', ''),
            'tnj_confidence': info.get('tnj_confidence', 'NONE')
        }
    
    return scores

def merge_disease_resistance_data(eggnog_data, interpro_data, resistance_threshold=5):
    merged_temp = merge_annotation_data(eggnog_data, interpro_data, SECONDARY_PATHWAYS)
    resistance_scores = calculate_disease_resistance_score(merged_temp, resistance_threshold)
    
    merged = {}
    for gene in merged_temp:
        merged[gene] = merged_temp[gene].copy()
        merged[gene]['resistance_potential'] = resistance_scores.get(gene, {}).get('potential', 'NONE')
        merged[gene]['resistance_score'] = resistance_scores.get(gene, {}).get('score', 0)
        merged[gene]['resistance_classes'] = resistance_scores.get(gene, {}).get('resistance_classes', [])
        merged[gene]['resistance_pathways'] = resistance_scores.get(gene, {}).get('resistance_pathways', [])
        merged[gene]['resistance_pfams'] = resistance_scores.get(gene, {}).get('resistance_pfams', [])
        merged[gene]['tnj_detected'] = resistance_scores.get(gene, {}).get('tnj_detected', False)
        merged[gene]['tnj_architecture'] = resistance_scores.get(gene, {}).get('tnj_architecture', '')
        merged[gene]['tnj_confidence'] = resistance_scores.get(gene, {}).get('tnj_confidence', 'NONE')
    
    return merged

# ============================================================================
# FUNÇÕES DE REDE PARA RESISTÊNCIA
# ============================================================================

def build_resistance_network(merged_data, output_prefix, min_edge_weight=1, max_nodes=15, selection_method='enrichment'):
    resistance_genes = {g: info for g, info in merged_data.items() if info.get('resistance_score', 0) > 0}
    if len(resistance_genes) < 2:
        print("   ⚠️ Poucos genes de resistência para construir rede")
        return
    
    if selection_method == 'score':
        selected = sorted(resistance_genes.keys(), key=lambda x: resistance_genes[x]['resistance_score'], reverse=True)[:max_nodes]
    elif selection_method == 'enrichment':
        importance = {}
        for gene, info in resistance_genes.items():
            score = info['resistance_score']
            n_classes = len(info.get('resistance_classes', []))
            n_pathways = len(info.get('resistance_pathways', []))
            n_pfams = len(info.get('resistance_pfams', []))
            tnj_bonus = 5 if info.get('tnj_detected', False) else 0
            importance[gene] = (score*2) + (n_classes*3) + n_pathways + (n_pfams*0.5) + tnj_bonus
        selected = sorted(importance.keys(), key=lambda x: importance[x], reverse=True)[:max_nodes]
    else:
        class_genes = defaultdict(list)
        for gene, info in resistance_genes.items():
            for rc in info.get('resistance_classes', []):
                if rc:
                    class_genes[rc].append((gene, info['resistance_score']))
        selected = []
        for rc, genes in class_genes.items():
            selected.extend([g[0] for g in sorted(genes, key=lambda x: x[1], reverse=True)[:3]])
        if len(selected) < max_nodes:
            remaining = [g for g in resistance_genes.keys() if g not in selected]
            selected.extend(sorted(remaining, key=lambda x: resistance_genes[x]['resistance_score'], reverse=True)[:max_nodes-len(selected)])
        selected = selected[:max_nodes]
    
    edges = []
    for i, g1 in enumerate(selected):
        for g2 in selected[i+1:]:
            shared_classes = len(set(merged_data[g1].get('resistance_classes', [])) & set(merged_data[g2].get('resistance_classes', [])))
            shared_pathways = len(set(merged_data[g1].get('resistance_pathways', [])) & set(merged_data[g2].get('resistance_pathways', [])))
            shared_pfams = len(set(merged_data[g1].get('resistance_pfams', [])) & set(merged_data[g2].get('resistance_pfams', [])))
            shared_kos = len(set(merged_data[g1]['ko_numbers']) & set(merged_data[g2]['ko_numbers']))
            
            tnj_bonus = 3 if (merged_data[g1].get('tnj_detected', False) and merged_data[g2].get('tnj_detected', False)) else 0
            
            total_sim = (shared_classes*3) + (shared_pathways*2) + shared_pfams + (shared_kos*2) + tnj_bonus
            if total_sim >= min_edge_weight:
                edges.append({
                    'gene1': g1, 'gene2': g2,
                    'shared_classes': shared_classes,
                    'shared_pathways': shared_pathways,
                    'shared_pfams': shared_pfams,
                    'shared_kos': shared_kos,
                    'similarity_score': total_sim
                })
    
    if edges:
        df_edges = pd.DataFrame(edges)
        df_edges.to_csv(f"{output_prefix}_resistance_network.csv", index=False)
        plot_resistance_network(merged_data, selected, df_edges, output_prefix)
    else:
        print("   ⚠️ Nenhuma aresta com similaridade suficiente para rede")

def plot_resistance_network(merged_data, genes_list, edges_df, output_prefix):
    resistance_colors = {
        'HIGH': '#2ECC71',
        'MEDIUM': '#F39C12', 
        'LOW': '#E74C3C',
        'NONE': '#95A5A6'
    }
    
    tnj_color = '#2A9D8F'
    
    n_nodes = len(genes_list)
    angles = np.linspace(0, 2*np.pi, n_nodes, endpoint=False)
    pos = {gene: (np.cos(angle), np.sin(angle)) for gene, angle in zip(genes_list, angles)}
    
    fig, ax = plt.subplots(figsize=(14, 12))
    
    if len(edges_df) > 0:
        max_sim = edges_df['similarity_score'].max() or 1
        for _, edge in edges_df.iterrows():
            g1, g2 = edge['gene1'], edge['gene2']
            if g1 in pos and g2 in pos:
                x1, y1 = pos[g1]
                x2, y2 = pos[g2]
                weight = edge['similarity_score'] / max_sim
                lw = 0.5 + weight*4
                if merged_data[g1].get('tnj_detected', False) and merged_data[g2].get('tnj_detected', False):
                    color = '#2A9D8F'
                    alpha = 0.8
                else:
                    color = '#3498DB'
                    alpha = 0.5
                ax.plot([x1, x2], [y1, y2], color, alpha=alpha, linewidth=lw, zorder=1)
    
    node_colors = []
    node_sizes = []
    for gene in genes_list:
        potential = merged_data[gene].get('resistance_potential', 'NONE')
        base_color = resistance_colors.get(potential, '#95A5A6')
        
        if merged_data[gene].get('tnj_detected', False):
            node_colors.append(tnj_color)
        else:
            node_colors.append(base_color)
        
        size = 500 + merged_data[gene].get('resistance_score', 0) * 60
        node_sizes.append(min(size, 1800))
    
    x_coords = [pos[gene][0] for gene in genes_list]
    y_coords = [pos[gene][1] for gene in genes_list]
    
    scatter = ax.scatter(x_coords, y_coords, s=node_sizes, c=node_colors, alpha=0.85, 
                         edgecolors='black', linewidth=1.5, zorder=2)
    
    for gene in genes_list:
        x, y = pos[gene]
        label = gene if len(gene) <= 25 else gene[:22] + '...'
        if merged_data[gene].get('tnj_detected', False):
            label = '🔬 ' + label
        ax.annotate(label, (x, y), xytext=(0, 0), textcoords='offset points', 
                    ha='center', va='center', fontsize=9, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85))
    
    ax.set_title('Functional Network of Disease Resistance Genes\n(🔬 = TNJ genes)', fontsize=14, fontweight='bold')
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect('equal')
    ax.axis('off')
    
    legend_elements = [
        Patch(facecolor='#2ECC71', label='HIGH'), 
        Patch(facecolor='#F39C12', label='MEDIUM'), 
        Patch(facecolor='#E74C3C', label='LOW'),
        Patch(facecolor='#2A9D8F', label='🔬 TNJ gene')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_resistance_network.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_resistance_network.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   ✅ Rede de resistência salva: {output_prefix}_resistance_network.pdf")

# ============================================================================
# FUNÇÕES DE VISUALIZAÇÃO PARA RESISTÊNCIA
# ============================================================================

def plot_resistance_score_ranking(merged_data, output_prefix):
    positive = {g: info for g, info in merged_data.items() if info.get('resistance_score', 0) > 0}
    if not positive:
        return
    
    sorted_genes = sorted(positive.items(), key=lambda x: x[1]['resistance_score'], reverse=True)[:20]
    genes = [g[0] for g in sorted_genes]
    scores = [g[1]['resistance_score'] for g in sorted_genes]
    potentials = [g[1]['resistance_potential'] for g in sorted_genes]
    tnj_status = [g[1].get('tnj_detected', False) for g in sorted_genes]
    
    colors = {'HIGH':'#2ECC71', 'MEDIUM':'#F39C12', 'LOW':'#E74C3C'}
    bar_colors = [colors.get(p, '#95A5A6') for p in potentials]
    
    fig, ax = plt.subplots(figsize=(16, max(8, len(genes)*0.4)))
    bars = ax.barh(range(len(genes)), scores, color=bar_colors, alpha=0.8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=9)
    ax.set_xlabel('Disease Resistance Score', fontsize=12, fontweight='bold')
    ax.set_title('Top 20 Genes by Disease Resistance Potential', fontsize=14, fontweight='bold')
    
    for bar, score, potential, tnj in zip(bars, scores, potentials, tnj_status):
        tnj_label = " 🔬TNJ" if tnj else ""
        ax.text(score + 0.5, bar.get_y() + bar.get_height()/2, 
                f'{score} ({potential}){tnj_label}', va='center', fontsize=10, fontweight='bold')
    
    legend_elements = [
        Patch(facecolor='#2ECC71', label='HIGH (≥10)'), 
        Patch(facecolor='#F39C12', label='MEDIUM (5-9)'), 
        Patch(facecolor='#E74C3C', label='LOW (1-4)'),
        Patch(facecolor='#2A9D8F', label='🔬 TNJ gene')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=10)
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_resistance_score_ranking.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_resistance_score_ranking.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_resistance_class_distribution(merged_data, output_prefix):
    class_counts = Counter()
    tnj_count = 0
    tnj_architectures = Counter()
    
    for gene, info in merged_data.items():
        if info.get('tnj_detected', False):
            tnj_count += 1
            arch = info.get('tnj_architecture', 'Unknown')
            tnj_architectures[arch] += 1
            class_counts['TNJ'] += 1
        
        for rc in info.get('resistance_classes', []):
            if info.get('tnj_detected', False) and rc in ['TIR_NLR', 'NLR', 'TNL']:
                continue
            class_counts[rc] += 1
    
    if tnj_count == 0 and 'TNJ' not in class_counts:
        class_counts['TNJ (não detectado)'] = 0
    
    if not class_counts:
        return
    
    class_colors = {
        'NLR': '#E63946',
        'TIR_NLR': '#F4A261',
        'TNJ': '#2A9D8F',
        'TNJ_like': '#264653',
        'TNL': '#FB8500',
        'CC_NLR': '#E9C46A',
        'RPW8': '#9B5DE5',
        'Jacalin': '#F15BB5',
        'LRR_repeat': '#FF006E',
        'NB_ARC': '#00BBF9',
        'TIR': '#00F5D4',
        'CC': '#FEE440',
        'NLR_Jacalin': '#8338EC',
        'NLR_Jacalin_hybrid': '#3A86FF',
        'TNJ (não detectado)': '#CCCCCC'
    }
    
    classes = list(class_counts.keys())
    counts = [class_counts[c] for c in classes]
    colors = [class_colors.get(c, '#999999') for c in classes]
    idx = np.argsort(counts)
    classes = [classes[i] for i in idx]
    counts = [counts[i] for i in idx]
    colors = [colors[i] for i in idx]
    
    fig, ax = plt.subplots(figsize=(14, max(6, len(classes)*0.5)))
    bars = ax.barh(classes, counts, color=colors, alpha=0.8)
    ax.set_xlabel('Number of genes', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of Disease Resistance Gene Classes\n(🔬 TNJ highlighted in green)', fontsize=14, fontweight='bold')
    
    for bar, count in zip(bars, counts):
        ax.text(count + 0.5, bar.get_y() + bar.get_height()/2, 
                str(count), va='center', fontsize=10, fontweight='bold')
    
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_resistance_class_distribution.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_resistance_class_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    with open(f"{output_prefix}_tnj_summary.txt", 'w') as f:
        f.write("TNJ GENE ANALYSIS - SUMMARY\n")
        f.write("=" * 50 + "\n")
        f.write(f"Total genes com TNJ detectados: {tnj_count}\n")
        if tnj_count > 0:
            f.write("\nArquiteturas TNJ encontradas:\n")
            for arch, count in tnj_architectures.most_common():
                f.write(f"  {arch}: {count}\n")
        else:
            f.write("\nNenhum gene TNJ detectado neste genoma.\n")
            f.write("Isso é consistente com a variação observada em Myrtaceae.\n")
    
    return tnj_count, tnj_architectures

def plot_resistance_pathway_completeness(merged_data, output_prefix):
    pathway_counts = Counter()
    for gene, info in merged_data.items():
        for pathway in info.get('resistance_pathways', []):
            pathway_counts[pathway] += 1
    if not pathway_counts:
        return
    
    pathways = list(pathway_counts.keys())
    counts = [pathway_counts[p] for p in pathways]
    names = [DISEASE_RESISTANCE_PATHWAYS.get(p, {}).get('name', p) for p in pathways]
    colors = [DISEASE_RESISTANCE_PATHWAYS.get(p, {}).get('color', '#2E86AB') for p in pathways]
    
    fig, ax = plt.subplots(figsize=(14, max(6, len(pathways)*0.4)))
    bars = ax.barh(range(len(pathways)), counts, color=colors, alpha=0.8)
    ax.set_yticks(range(len(pathways)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel('Number of genes', fontsize=12, fontweight='bold')
    ax.set_title('Genes Associated with Disease Resistance Pathways', fontsize=14, fontweight='bold')
    
    for bar, count in zip(bars, counts):
        ax.text(count + 0.5, bar.get_y() + bar.get_height()/2, str(count), 
                va='center', fontsize=9, fontweight='bold')
    
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_resistance_pathway_completeness.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_resistance_pathway_completeness.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_resistance_donut_chart(merged_data, output_prefix):
    potential_counts = Counter(info['resistance_potential'] for info in merged_data.values() 
                              if info.get('resistance_potential', 'NONE') != 'NONE')
    if not potential_counts:
        return
    
    colors_donut = {'HIGH':'#2ECC71','MEDIUM':'#F39C12','LOW':'#E74C3C'}
    colors = [colors_donut.get(p, '#95A5A6') for p in potential_counts.keys()]
    labels = [f"{p}\n({count} genes)" for p, count in potential_counts.items()]
    
    fig, ax = plt.subplots(figsize=(10, 10))
    wedges, texts, autotexts = ax.pie(potential_counts.values(), labels=labels, colors=colors, 
                                      autopct='%1.1f%%', startangle=90, pctdistance=0.85, 
                                      textprops={'fontsize':11, 'fontweight':'bold'})
    centre_circle = Circle((0, 0), 0.70, fc='white', linewidth=2, edgecolor='black')
    ax.add_artist(centre_circle)
    total = sum(potential_counts.values())
    ax.text(0, 0, f'Total\n{total}\ngenes', ha='center', va='center', 
            fontsize=14, fontweight='bold')
    ax.set_title('Distribution of Disease Resistance Potential', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_resistance_donut.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_resistance_donut.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_tnj_analysis(merged_data, output_prefix):
    tnj_genes = []
    
    for gene, info in merged_data.items():
        if info.get('tnj_detected', False):
            tnj_genes.append({
                'gene': gene,
                'architecture': info.get('tnj_architecture', ''),
                'confidence': info.get('tnj_confidence', 'NONE'),
                'resistance_score': info.get('resistance_score', 0),
                'resistance_potential': info.get('resistance_potential', 'NONE'),
                'pfams': ', '.join(list(info.get('resistance_pfams', []))),
                'description': info.get('description', '')[:100]
            })
    
    if not tnj_genes:
        print("   ⚠️ Nenhum gene TNJ detectado neste genoma")
        with open(f"{output_prefix}_tnj_analysis.txt", 'w') as f:
            f.write("TNJ GENE ANALYSIS\n")
            f.write("=" * 50 + "\n")
            f.write("Nenhum gene TNJ (TIR-NBS-Jacalin) foi detectado neste genoma.\n")
            f.write("Isso é consistente com a variação observada em Myrtaceae.\n")
        return
    
    df = pd.DataFrame(tnj_genes)
    df.to_csv(f"{output_prefix}_tnj_genes_table.csv", index=False)
    
    arch_counts = Counter([g['architecture'] for g in tnj_genes])
    if arch_counts:
        fig, ax = plt.subplots(figsize=(10, 6))
        arch_names = list(arch_counts.keys())
        arch_counts_values = list(arch_counts.values())
        colors = ['#2A9D8F' if 'canônico' in a else '#264653' for a in arch_names]
        
        bars = ax.bar(arch_names, arch_counts_values, color=colors, alpha=0.8)
        ax.set_xlabel('TNJ Architecture', fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of genes', fontsize=12, fontweight='bold')
        ax.set_title('TNJ Gene Architectures Detected', fontsize=14, fontweight='bold')
        
        for bar, count in zip(bars, arch_counts_values):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, 
                    str(count), ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        plt.xticks(rotation=15, ha='right')
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_tnj_architectures.pdf", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_prefix}_tnj_architectures.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    conf_counts = Counter([g['confidence'] for g in tnj_genes])
    if conf_counts:
        fig, ax = plt.subplots(figsize=(8, 6))
        conf_names = list(conf_counts.keys())
        conf_counts_values = list(conf_counts.values())
        colors_conf = {'HIGH': '#2ECC71', 'MEDIUM': '#F39C12', 'LOW': '#E74C3C'}
        bar_colors = [colors_conf.get(c, '#95A5A6') for c in conf_names]
        
        bars = ax.bar(conf_names, conf_counts_values, color=bar_colors, alpha=0.8)
        ax.set_xlabel('Confidence Level', fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of genes', fontsize=12, fontweight='bold')
        ax.set_title('TNJ Gene Confidence Levels', fontsize=14, fontweight='bold')
        
        for bar, count in zip(bars, conf_counts_values):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, 
                    str(count), ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_tnj_confidence.pdf", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_prefix}_tnj_confidence.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    with open(f"{output_prefix}_tnj_analysis.txt", 'w') as f:
        f.write("TNJ GENE ANALYSIS\n")
        f.write("=" * 50 + "\n")
        f.write(f"Total genes TNJ detectados: {len(tnj_genes)}\n")
        f.write(f"\nArquiteturas encontradas:\n")
        for arch, count in arch_counts.most_common():
            f.write(f"  {arch}: {count}\n")
        f.write(f"\nDistribuição por confiança:\n")
        for conf, count in conf_counts.most_common():
            f.write(f"  {conf}: {count}\n")
        f.write(f"\nDistribuição por potencial de resistência:\n")
        potential_counts = Counter([g['resistance_potential'] for g in tnj_genes])
        for pot, count in potential_counts.most_common():
            f.write(f"  {pot}: {count}\n")

def generate_publication_table_resistance(merged_data, output_prefix):
    rows = []
    for gene, info in merged_data.items():
        if info.get('resistance_score', 0) > 0:
            tnj_flag = "✅" if info.get('tnj_detected', False) else ""
            rows.append({
                'Gene_ID': gene,
                'Resistance_Potential': info['resistance_potential'],
                'Resistance_Score': info['resistance_score'],
                'TNJ': tnj_flag,
                'TNJ_Architecture': info.get('tnj_architecture', ''),
                'Description': info['description'][:100],
                'Resistance_Classes': ', '.join(info.get('resistance_classes', [])),
                'Resistance_Pathways': ', '.join(info.get('resistance_pathways', [])[:3]),
                'Resistance_PFAMs': ', '.join(info.get('resistance_pfams', [])[:5]),
                'KEGG_Orthologs': ', '.join(list(info['ko_numbers'])[:10]),
                'EC_Numbers': ', '.join(list(info['ec_numbers'])[:5])
            })
    
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(f"{output_prefix}_resistance_genes_table.csv", index=False)
        try:
            df.to_excel(f"{output_prefix}_resistance_genes_table.xlsx", index=False)
        except:
            pass
        
        with open(f"{output_prefix}_resistance_statistics.txt", 'w') as f:
            f.write("DISEASE RESISTANCE GENE ANALYSIS\n")
            f.write(f"Total genes: {len(merged_data)}\n")
            
            resistance_genes = len(rows)
            f.write(f"\nGenes with disease resistance potential: {resistance_genes}\n")
            f.write(f"  Percentage of total: {resistance_genes/len(merged_data)*100:.1f}%\n")
            
            potential_counts = Counter([r['Resistance_Potential'] for r in rows])
            for level in ['HIGH','MEDIUM','LOW']:
                f.write(f"  {level}: {potential_counts.get(level, 0)} genes\n")
            
            tnj_genes = [r for r in rows if r['TNJ'] == '✅']
            f.write(f"\nTNJ genes detected: {len(tnj_genes)}\n")
            if tnj_genes:
                arch_counts = Counter([r['TNJ_Architecture'] for r in tnj_genes])
                f.write("  Architectures:\n")
                for arch, count in arch_counts.most_common():
                    f.write(f"    {arch}: {count}\n")
            
            class_counter = Counter()
            for row in rows:
                for c in row['Resistance_Classes'].split(', '):
                    if c:
                        class_counter[c] += 1
            f.write("\nResistance classes:\n")
            for c, cnt in class_counter.most_common():
                f.write(f"  {c}: {cnt}\n")

def create_resistance_pfam_dotplot(merged_data, output_prefix, top_n=15):
    high_medium_genes = [g for g, info in merged_data.items() 
                        if info.get('resistance_potential') in ['HIGH', 'MEDIUM']]
    
    background_genes = [g for g, info in merged_data.items() 
                       if info.get('resistance_potential') in ['LOW', 'NONE']]
    
    if len(high_medium_genes) < 2:
        print("   ⚠️ Poucos genes HIGH+MEDIUM para análise")
        return
    
    if len(background_genes) < 2:
        print("   ⚠️ Poucos genes background para análise")
        return
    
    print(f"   Grupo HIGH+MEDIUM: {len(high_medium_genes)} genes")
    print(f"   Background (LOW+NONE): {len(background_genes)} genes")
    
    df_enrich = calculate_enrichment(
        high_medium_genes, background_genes, merged_data,
        min_count=2, min_fold=1.0, filter_generic=False
    )
    
    if len(df_enrich) == 0:
        print("   ⚠️ Nenhum PFAM enriquecido encontrado")
        return
    
    df_sig = df_enrich[df_enrich['significant']]
    
    if len(df_sig) == 0:
        print("   ⚠️ Nenhum PFAM significativamente enriquecido (p.adjust < 0.05)")
        return
    
    df_sig.to_csv(f"{output_prefix}_resistance_pfam_enrichment.csv", index=False)
    print(f"   ✅ {len(df_sig)} PFAMs significativamente enriquecidos")
    
    top_terms = df_sig.sort_values('fold_enrichment', ascending=False).head(top_n)
    
    resistance_colors = {c: DISEASE_RESISTANCE_PFAMS.get(c, {}).get('color', '#CCCCCC') 
                        for c in DISEASE_RESISTANCE_PFAMS.keys()}
    
    pfam_to_class = {}
    for class_name, class_info in DISEASE_RESISTANCE_PFAMS.items():
        for pfam in class_info['pfams']:
            if pfam not in pfam_to_class:
                pfam_to_class[pfam] = class_name
    
    colors = [resistance_colors.get(pfam_to_class.get(pfam, 'Other'), '#CCCCCC') 
              for pfam in top_terms['pfam_domain']]
    
    fig_height = max(6, len(top_terms) * 0.4)
    fig, ax = plt.subplots(figsize=(16, fig_height))
    
    scatter = ax.scatter(
        top_terms['fold_enrichment'],
        top_terms['pfam_domain'],
        s=top_terms['count_in_group'] * 15 + 30,
        c=colors,
        alpha=0.8,
        edgecolors='black',
        linewidth=0.5
    )
    
    ax.set_xlabel('Fold Enrichment', fontsize=12, fontweight='bold')
    ax.set_ylabel('PFAM Domain', fontsize=12, fontweight='bold')
    ax.set_title('Enriched PFAM Domains in HIGH+MEDIUM Disease Resistance Genes\n(Compared to LOW+NONE Background)', 
                 fontsize=14, fontweight='bold')
    ax.axvline(x=1, color='red', linestyle='--', alpha=0.5, label='Fold = 1')
    
    unique_classes = set(pfam_to_class.get(pfam, 'Other') for pfam in top_terms['pfam_domain'])
    legend_elements = [Patch(facecolor=resistance_colors.get(c, '#CCCCCC'), label=c)
                       for c in unique_classes if c in resistance_colors]
    
    if legend_elements:
        ax.legend(handles=legend_elements,
                 loc='center left',
                 bbox_to_anchor=(1.02, 0.5),
                 fontsize=9,
                 title='Resistance Class',
                 title_fontsize=10)
    
    plt.subplots_adjust(right=0.78)
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(f"{output_prefix}_resistance_pfam_dotplot.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_resistance_pfam_dotplot.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   ✅ Dotplot de resistência salvo: {output_prefix}_resistance_pfam_dotplot.pdf")
    
    return df_sig

# ============================================================================
# FUNÇÕES DE PHARMA ESPECÍFICAS - VISUALIZAÇÕES
# ============================================================================

def plot_pharma_score_ranking(merged_data, output_prefix):
    positive = {g: info for g, info in merged_data.items() if info.get('pharma_score', 0) > 0}
    if not positive:
        return
    
    sorted_genes = sorted(positive.items(), key=lambda x: x[1]['pharma_score'], reverse=True)[:20]
    genes = [g[0] for g in sorted_genes]
    scores = [g[1]['pharma_score'] for g in sorted_genes]
    potentials = [g[1]['pharma_potential'] for g in sorted_genes]
    colors = {'HIGH':'#2ECC71', 'MEDIUM':'#F39C12', 'LOW':'#E74C3C'}
    bar_colors = [colors.get(p, '#95A5A6') for p in potentials]
    
    fig, ax = plt.subplots(figsize=(16, max(8, len(genes)*0.4)))
    bars = ax.barh(range(len(genes)), scores, color=bar_colors, alpha=0.8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=9)
    ax.set_xlabel('Pharmaceutical Potential Score', fontsize=12, fontweight='bold')
    ax.set_title('Top 20 Genes by Pharmaceutical & Nutraceutical Potential', fontsize=14, fontweight='bold')
    
    for bar, score, potential in zip(bars, scores, potentials):
        ax.text(score + 0.5, bar.get_y() + bar.get_height()/2, 
                f'{score} ({potential})', va='center', fontsize=10, fontweight='bold')
    
    legend_elements = [Patch(facecolor='#2ECC71', label='HIGH (≥10)'), 
                       Patch(facecolor='#F39C12', label='MEDIUM (4-9)'), 
                       Patch(facecolor='#E74C3C', label='LOW (1-3)')]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=10)
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_pharma_score_ranking.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_pharma_score_ranking.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_compound_class_distribution(merged_data, output_prefix):
    class_counts = defaultdict(int)
    for gene, info in merged_data.items():
        for cc in info.get('compound_classes', []):
            class_counts[cc] += 1
    if not class_counts:
        return
    
    class_colors = {'Alkaloids':'#E63946','Terpenoids':'#2A9D8F','Phenylpropanoids':'#E9C46A',
                    'Flavonoids':'#F4A261','Stilbenoids':'#E76F51','Glucosinolates':'#9B5DE5',
                    'Xanthines':'#F15BB5','Antibiotics':'#00BBF9','Cannabinoids':'#FF006E','Other':'#CCCCCC'}
    
    classes = list(class_counts.keys())
    counts = [class_counts[c] for c in classes]
    colors = [class_colors.get(c, '#999999') for c in classes]
    idx = np.argsort(counts)
    classes = [classes[i] for i in idx]
    counts = [counts[i] for i in idx]
    colors = [colors[i] for i in idx]
    
    fig, ax = plt.subplots(figsize=(14, max(6, len(classes)*0.5)))
    bars = ax.barh(classes, counts, color=colors, alpha=0.8)
    ax.set_xlabel('Number of genes', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of Bioactive Compound Classes', fontsize=14, fontweight='bold')
    
    for bar, count in zip(bars, counts):
        ax.text(count + 0.5, bar.get_y() + bar.get_height()/2, 
                str(count), va='center', fontsize=10, fontweight='bold')
    
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_compound_class_distribution.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_compound_class_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_pharmaceutical_uses(merged_data, output_prefix):
    uses_counter = Counter()
    for gene, info in merged_data.items():
        for use in info.get('pharma_uses', []):
            uses_counter[use] += 1
    if not uses_counter:
        return
    
    top = uses_counter.most_common(15)
    uses_names, uses_counts = zip(*top)
    
    fig, ax = plt.subplots(figsize=(14, max(6, len(uses_names)*0.4)))
    bars = ax.barh(range(len(uses_names)), uses_counts, color='#9B5DE5', alpha=0.8)
    ax.set_yticks(range(len(uses_names)))
    ax.set_yticklabels(uses_names, fontsize=9)
    ax.set_xlabel('Number of associated genes', fontsize=12, fontweight='bold')
    ax.set_title('Predicted Pharmaceutical Applications', fontsize=14, fontweight='bold')
    
    for bar, count in zip(bars, uses_counts):
        ax.text(count + 0.5, bar.get_y() + bar.get_height()/2, 
                str(count), va='center', fontsize=9, fontweight='bold')
    
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_pharmaceutical_uses.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_pharmaceutical_uses.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_sunburst_hierarchical(merged_data, output_prefix):
    pathway_class_counts = defaultdict(lambda: defaultdict(int))
    for gene, info in merged_data.items():
        for pathway in info['secondary_pathways']:
            if pathway in PHARMA_DETAILS:
                pc = PHARMA_DETAILS[pathway]['class']
                pn = PHARMA_DETAILS[pathway]['name']
                pathway_class_counts[pc][pn] += 1
    if not pathway_class_counts:
        return
    
    fig, ax = plt.subplots(figsize=(14, 8))
    classes = list(pathway_class_counts.keys())
    bottom = np.zeros(len(classes))
    color_idx = 0
    total_pathways = sum(len(pathway_class_counts[c]) for c in classes)
    colors_pathways = plt.cm.Set3(np.linspace(0, 1, total_pathways))
    
    for i, class_name in enumerate(classes):
        pathways = list(pathway_class_counts[class_name].keys())
        counts = [pathway_class_counts[class_name][p] for p in pathways]
        bars = ax.bar([class_name], [sum(counts)], bottom=bottom[i], 
                      color=colors_pathways[color_idx:color_idx+len(pathways)], alpha=0.8)
        bottom[i] += sum(counts)
        color_idx += len(pathways)
    
    ax.set_ylabel('Number of genes', fontsize=12, fontweight='bold')
    ax.set_title('Hierarchical Distribution of Bioactive Pathways', fontsize=14, fontweight='bold')
    ax.grid(True, axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_sunburst_pharma.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_sunburst_pharma.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_donut_chart_potential(merged_data, output_prefix):
    potential_counts = Counter(info['pharma_potential'] for info in merged_data.values() 
                              if info.get('pharma_potential', 'NONE') != 'NONE')
    if not potential_counts:
        return
    
    colors_donut = {'HIGH':'#2ECC71','MEDIUM':'#F39C12','LOW':'#E74C3C'}
    colors = [colors_donut.get(p, '#95A5A6') for p in potential_counts.keys()]
    labels = [f"{p}\n({count} genes)" for p, count in potential_counts.items()]
    
    fig, ax = plt.subplots(figsize=(10, 10))
    wedges, texts, autotexts = ax.pie(potential_counts.values(), labels=labels, colors=colors, 
                                      autopct='%1.1f%%', startangle=90, pctdistance=0.85, 
                                      textprops={'fontsize':11, 'fontweight':'bold'})
    centre_circle = Circle((0, 0), 0.70, fc='white', linewidth=2, edgecolor='black')
    ax.add_artist(centre_circle)
    total = sum(potential_counts.values())
    ax.text(0, 0, f'Total\n{total}\ngenes', ha='center', va='center', 
            fontsize=14, fontweight='bold')
    ax.set_title('Distribution of Pharmaceutical Potential', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_donut_potential.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_donut_potential.png", dpi=300, bbox_inches='tight')
    plt.close()

def build_functional_network(merged_data, output_prefix, min_edge_weight=1, max_nodes=15, selection_method='enrichment'):
    pharma_genes = {g: info for g, info in merged_data.items() if info.get('pharma_score', 0) > 0}
    if len(pharma_genes) < 2:
        return
    
    if selection_method == 'score':
        selected = sorted(pharma_genes.keys(), key=lambda x: pharma_genes[x]['pharma_score'], reverse=True)[:max_nodes]
    elif selection_method == 'enrichment':
        importance = {}
        for gene, info in pharma_genes.items():
            score = info['pharma_score']
            n_kos = len(info['ko_numbers'])
            n_pathways = len(info['secondary_pathways'])
            n_pfams = len([p for p in info['pfam_domains'] if p not in GENERIC_DOMAINS])
            importance[gene] = (score*2) + (n_kos*3) + n_pathways + (n_pfams*0.5)
        selected = sorted(importance.keys(), key=lambda x: importance[x], reverse=True)[:max_nodes]
    else:
        class_genes = defaultdict(list)
        for gene, info in pharma_genes.items():
            for cc in info.get('compound_classes', []):
                if cc:
                    class_genes[cc].append((gene, info['pharma_score']))
        selected = []
        for cc, genes in class_genes.items():
            selected.extend([g[0] for g in sorted(genes, key=lambda x: x[1], reverse=True)[:3]])
        if len(selected) < max_nodes:
            remaining = [g for g in pharma_genes.keys() if g not in selected]
            selected.extend(sorted(remaining, key=lambda x: pharma_genes[x]['pharma_score'], reverse=True)[:max_nodes-len(selected)])
        selected = selected[:max_nodes]
    
    edges = []
    for i, g1 in enumerate(selected):
        for g2 in selected[i+1:]:
            shared_pathways = len(set(merged_data[g1]['secondary_pathways']) & set(merged_data[g2]['secondary_pathways']))
            shared_go = len(set(merged_data[g1]['go_terms']) & set(merged_data[g2]['go_terms']))
            shared_kos = len(set(merged_data[g1]['ko_numbers']) & set(merged_data[g2]['ko_numbers']))
            total_sim = (shared_pathways*2) + shared_go + (shared_kos*3)
            if total_sim >= min_edge_weight:
                edges.append({
                    'gene1': g1, 'gene2': g2,
                    'shared_pathways': shared_pathways,
                    'shared_go': shared_go,
                    'shared_kos': shared_kos,
                    'similarity_score': total_sim
                })
    
    if edges:
        df_edges = pd.DataFrame(edges)
        df_edges.to_csv(f"{output_prefix}_functional_network.csv", index=False)
        plot_functional_network(merged_data, selected, df_edges, output_prefix)

def plot_functional_network(merged_data, genes_list, edges_df, output_prefix):
    potential_colors = {'HIGH':'#2ECC71','MEDIUM':'#F39C12','LOW':'#E74C3C','NONE':'#95A5A6'}
    n_nodes = len(genes_list)
    angles = np.linspace(0, 2*np.pi, n_nodes, endpoint=False)
    pos = {gene: (np.cos(angle), np.sin(angle)) for gene, angle in zip(genes_list, angles)}
    
    fig, ax = plt.subplots(figsize=(14, 12))
    
    if len(edges_df) > 0:
        max_sim = edges_df['similarity_score'].max() or 1
        for _, edge in edges_df.iterrows():
            g1, g2 = edge['gene1'], edge['gene2']
            if g1 in pos and g2 in pos:
                x1, y1 = pos[g1]
                x2, y2 = pos[g2]
                weight = edge['similarity_score'] / max_sim
                lw = 0.5 + weight*4
                ax.plot([x1, x2], [y1, y2], '#3498DB', alpha=0.6, linewidth=lw, zorder=1)
    
    node_colors = [potential_colors.get(merged_data[gene]['pharma_potential'], '#95A5A6') for gene in genes_list]
    node_sizes = [500 + merged_data[gene]['pharma_score']*80 for gene in genes_list]
    node_sizes = [min(s, 1800) for s in node_sizes]
    
    x_coords = [pos[gene][0] for gene in genes_list]
    y_coords = [pos[gene][1] for gene in genes_list]
    ax.scatter(x_coords, y_coords, s=node_sizes, c=node_colors, alpha=0.85, 
               edgecolors='black', linewidth=1.5, zorder=2)
    
    for gene in genes_list:
        x, y = pos[gene]
        label = gene if len(gene) <= 25 else gene[:22]+'...'
        ax.annotate(label, (x, y), xytext=(0, 0), textcoords='offset points', 
                    ha='center', va='center', fontsize=9, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85))
    
    ax.set_title('Functional Network of Pharmaceutical Potential Genes', fontsize=14, fontweight='bold')
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect('equal')
    ax.axis('off')
    
    legend_elements = [Patch(facecolor='#2ECC71', label='HIGH'), 
                       Patch(facecolor='#F39C12', label='MEDIUM'), 
                       Patch(facecolor='#E74C3C', label='LOW')]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_functional_network.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_functional_network.png", dpi=300, bbox_inches='tight')
    plt.close()

def predict_novel_compounds(merged_data, output_prefix):
    compound_pfam_map = {
        'Flavonoids': {'pfams':['PF00195','PF02797','PF02458','PF05834'], 'compounds':['naringenin','quercetin','kaempferol','anthocyanins'], 'confidence':0.85},
        'Terpenoids': {'pfams':['PF03936','PF01397','PF00067','PF00494'], 'compounds':['artemisinin','taxol','limonene','ginsenosides'], 'confidence':0.80},
        'Alkaloids': {'pfams':['PF01596','PF00891','PF00248','PF00201'], 'compounds':['morphine','codeine','berberine','caffeine'], 'confidence':0.70},
        'Phenylpropanoids': {'pfams':['PF00195','PF00430','PF00141'], 'compounds':['curcumin','resveratrol','lignin'], 'confidence':0.75},
        'Carotenoids': {'pfams':['PF00514','PF06444','PF00494'], 'compounds':['β-carotene','lycopene','astaxanthin'], 'confidence':0.80},
        'Diterpenoids': {'pfams':['PF00067','PF03936'], 'compounds':['tanshinone','forskolin','gibberellin'], 'confidence':0.75},
        'Flavonoid_advanced': {'pfams':['PF00195','PF02797'], 'compounds':['quercetin','kaempferol','myricetin'], 'confidence':0.82},
        'P450_terpenoid': {'pfams':['PF00067','PF01397'], 'compounds':['artemisinin','taxol','limonene'], 'confidence':0.78}
    }
    
    predictions = []
    for gene, info in merged_data.items():
        if info.get('pharma_score', 0) == 0:
            continue
        gene_pfams = set(info['pfam_domains'])
        predicted = []
        for cc, data in compound_pfam_map.items():
            required = set(data['pfams'])
            matched = required & gene_pfams
            if matched:
                match_score = len(matched)/len(required)
                conf = data['confidence'] * match_score
                predicted.append({'class': cc, 'compounds': data['compounds'][0], 'confidence': conf})
        if predicted:
            max_conf = max([p['confidence'] for p in predicted])
            predictions.append({
                'gene': gene,
                'pharma_score': info['pharma_score'],
                'pharma_potential': info['pharma_potential'],
                'predicted_compound_classes': ', '.join([p['class'] for p in predicted[:3]]),
                'predicted_compounds': ', '.join([p['compounds'] for p in predicted[:3]]),
                'max_confidence': round(max_conf, 3)
            })
    
    if predictions:
        df = pd.DataFrame(predictions)
        df.to_csv(f"{output_prefix}_compound_predictions.csv", index=False)
        
        fig, ax = plt.subplots(figsize=(14, max(8, len(df.head(15))*0.4)))
        top = df.head(15)
        colors = {'HIGH':'#2ECC71','MEDIUM':'#F39C12','LOW':'#E74C3C'}
        bar_colors = [colors.get(p, '#95A5A6') for p in top['pharma_potential']]
        bars = ax.barh(range(len(top)), top['max_confidence'].values, color=bar_colors, alpha=0.8)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top['gene'].values, fontsize=9)
        ax.set_xlabel('Prediction Confidence', fontsize=12, fontweight='bold')
        ax.set_title('Top Predicted Compound-Producing Genes', fontsize=14, fontweight='bold')
        
        for bar, conf, compounds in zip(bars, top['max_confidence'].values, top['predicted_compounds'].values):
            ax.text(conf + 0.02, bar.get_y() + bar.get_height()/2, 
                    f'{conf:.0%} - {compounds[:40]}', va='center', fontsize=8)
        
        ax.set_xlim(0, 1)
        ax.grid(True, axis='x', alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_compound_predictions.pdf", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_prefix}_compound_predictions.png", dpi=300, bbox_inches='tight')
        plt.close()

def generate_publication_table_pharma(merged_data, output_prefix):
    rows = []
    for gene, info in merged_data.items():
        if info.get('pharma_score', 0) > 0:
            rows.append({
                'Gene_ID': gene,
                'Pharma_Potential': info['pharma_potential'],
                'Pharma_Score': info['pharma_score'],
                'Description': info['description'][:100],
                'Compound_Classes': ', '.join(info.get('compound_classes', [])),
                'Pharmaceutical_Uses': ', '.join(info.get('pharma_uses', [])[:3]),
                'Pharma_Pathways': ', '.join(info.get('pharma_pathways', [])[:3]),
                'KEGG_Orthologs': ', '.join(list(info['ko_numbers'])[:10]),
                'EC_Numbers': ', '.join(list(info['ec_numbers'])[:5]),
                'PFAM_Domains': ', '.join([p for p in list(info['pfam_domains'])[:5] if p not in GENERIC_DOMAINS])
            })
    
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(f"{output_prefix}_pharma_genes_table.csv", index=False)
        try:
            df.to_excel(f"{output_prefix}_pharma_genes_table.xlsx", index=False)
        except:
            pass
        
        with open(f"{output_prefix}_pharma_statistics.txt", 'w') as f:
            f.write("PHARMACEUTICAL AND NUTRACEUTICAL POTENTIAL ANALYSIS\n")
            f.write(f"Total genes: {len(merged_data)}\n")
            
            secondary_genes = sum(1 for g in merged_data.values() if g['secondary_pathways'])
            f.write(f"Secondary metabolism genes: {secondary_genes}\n")
            f.write(f"  Percentage of total: {secondary_genes/len(merged_data)*100:.1f}%\n")
            
            pharma_genes = len(rows)
            f.write(f"\nGenes with pharmaceutical potential: {pharma_genes}\n")
            f.write(f"  Percentage of secondary genes: {pharma_genes/secondary_genes*100:.1f}%\n")
            
            potential_counts = Counter([r['Pharma_Potential'] for r in rows])
            for level in ['HIGH','MEDIUM','LOW']:
                f.write(f"  {level}: {potential_counts.get(level, 0)} genes\n")
            
            class_counter = Counter()
            for row in rows:
                for c in row['Compound_Classes'].split(', '):
                    if c:
                        class_counter[c] += 1
            f.write("\nCompound classes:\n")
            for c, cnt in class_counter.most_common():
                f.write(f"  {c}: {cnt}\n")

def analyze_pathway_completeness(merged_data, output_prefix):
    pathway_genes = defaultdict(list)
    for gene, info in merged_data.items():
        for pathway in info['secondary_pathways']:
            if pathway in SECONDARY_PATHWAYS:
                pathway_genes[pathway].append(gene)
    
    results = []
    for pathway, genes in pathway_genes.items():
        pathway_info = SECONDARY_PATHWAYS.get(pathway, {})
        results.append({
            'pathway_id': pathway,
            'pathway_name': pathway_info.get('name', pathway),
            'class': pathway_info.get('class', 'Unknown'),
            'gene_count': len(genes),
            'pharma_score_sum': sum(merged_data[g].get('pharma_score', 0) for g in genes),
            'high_genes': sum(1 for g in genes if merged_data[g].get('pharma_potential') == 'HIGH'),
            'medium_genes': sum(1 for g in genes if merged_data[g].get('pharma_potential') == 'MEDIUM'),
            'avg_pharma_score': round(sum(merged_data[g].get('pharma_score', 0) for g in genes) / len(genes), 2) if genes else 0
        })
    
    if results:
        df = pd.DataFrame(results)
        df = df.sort_values('pharma_score_sum', ascending=False)
        df.to_csv(f"{output_prefix}_pathway_completeness_analysis.csv", index=False)
        
        fig, ax = plt.subplots(figsize=(14, max(6, len(df)*0.3)))
        top_pathways = df.head(20)
        bars = ax.barh(range(len(top_pathways)), top_pathways['gene_count'].values, 
                       color='#2E86AB', alpha=0.7)
        ax.set_yticks(range(len(top_pathways)))
        ax.set_yticklabels(top_pathways['pathway_name'].values, fontsize=9)
        ax.set_xlabel('Number of genes', fontsize=12, fontweight='bold')
        ax.set_title('Top Pathways by Gene Count', fontsize=14, fontweight='bold')
        
        for bar, count in zip(bars, top_pathways['gene_count'].values):
            ax.text(count + 0.5, bar.get_y() + bar.get_height()/2, str(count), 
                    va='center', fontsize=9, fontweight='bold')
        
        ax.invert_yaxis()
        ax.grid(True, axis='x', alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_pathway_completeness_analysis.pdf", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_prefix}_pathway_completeness_analysis.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   ✅ Pathway completeness analysis salvo: {output_prefix}_pathway_completeness_analysis.csv")
    
    return df

# ============================================================================
# FUNÇÕES DE ENRIQUECIMENTO COMBINADO
# ============================================================================

def create_pfam_dotplot_combined(merged_data, output_prefix, top_n=15, group_key='pharma_potential'):
    print(f"\n📊 Gerando dotplot combinado HIGH+MEDIUM vs LOW+NONE...")
    
    high_medium_genes = [g for g, info in merged_data.items() 
                        if info.get(group_key) in ['HIGH', 'MEDIUM']]
    
    background_genes = [g for g, info in merged_data.items() 
                       if info.get(group_key) in ['LOW', 'NONE']]
    
    if len(high_medium_genes) < 2:
        print("   ⚠️ Poucos genes HIGH+MEDIUM para análise")
        return
    
    if len(background_genes) < 2:
        print("   ⚠️ Poucos genes background para análise")
        return
    
    print(f"   Grupo HIGH+MEDIUM: {len(high_medium_genes)} genes")
    print(f"   Background (LOW+NONE): {len(background_genes)} genes")
    
    df_enrich = calculate_enrichment(
        high_medium_genes, background_genes, merged_data,
        min_count=2, min_fold=1.0, filter_generic=True
    )
    
    if len(df_enrich) == 0:
        print("   ⚠️ Nenhum PFAM enriquecido encontrado")
        return
    
    df_sig = df_enrich[df_enrich['significant']]
    
    if len(df_sig) == 0:
        print("   ⚠️ Nenhum PFAM significativamente enriquecido (p.adjust < 0.05)")
        return
    
    df_sig.to_csv(f"{output_prefix}_combined_high_medium_pfam_enrichment.csv", index=False)
    print(f"   ✅ {len(df_sig)} PFAMs significativamente enriquecidos")
    
    top_terms = df_sig.sort_values('fold_enrichment', ascending=False).head(top_n)
    
    pathway_colors = {
        'Flavonoids': '#E76F51',
        'Terpenoids': '#2A9D8F',
        'Alkaloids': '#E63946',
        'Phenylpropanoids': '#E9C46A',
        'Carotenoids': '#F4A261',
        'P450': '#9B5DE5',
        'Methyltransferases': '#F15BB5',
        'Other': '#CCCCCC'
    }
    
    fig_height = max(6, len(top_terms) * 0.4)
    fig, ax = plt.subplots(figsize=(16, fig_height))
    
    colors = [pathway_colors.get(cat, '#CCCCCC') for cat in top_terms['pathway_category']]
    
    scatter = ax.scatter(
        top_terms['fold_enrichment'],
        top_terms['pfam_domain'],
        s=top_terms['count_in_group'] * 15 + 30,
        c=colors,
        alpha=0.8,
        edgecolors='black',
        linewidth=0.5
    )
    
    ax.set_xlabel('Fold Enrichment', fontsize=12, fontweight='bold')
    ax.set_ylabel('PFAM Domain', fontsize=12, fontweight='bold')
    ax.set_title('Enriched PFAM Domains in HIGH+MEDIUM Pharmaceutical Potential Genes\n(Compared to LOW+NONE Background)', 
                 fontsize=14, fontweight='bold')
    ax.axvline(x=1, color='red', linestyle='--', alpha=0.5, label='Fold = 1')
    
    unique_cats = top_terms['pathway_category'].unique()
    legend_elements = [Patch(facecolor=pathway_colors.get(cat, '#CCCCCC'), label=cat)
                       for cat in unique_cats if cat in pathway_colors]
    
    if legend_elements:
        ax.legend(handles=legend_elements,
                 loc='center left',
                 bbox_to_anchor=(1.02, 0.5),
                 fontsize=10,
                 title='PFAM Category',
                 title_fontsize=11)
    
    plt.subplots_adjust(right=0.78)
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(f"{output_prefix}_combined_high_medium_pfam_dotplot.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_combined_high_medium_pfam_dotplot.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   ✅ Dotplot combinado salvo: {output_prefix}_combined_high_medium_pfam_dotplot.pdf")
    
    category_summary = df_sig.groupby('pathway_category').size().reset_index(name='count')
    category_summary = category_summary.sort_values('count', ascending=False)
    category_summary.to_csv(f"{output_prefix}_combined_pfam_category_summary.csv", index=False)
    print("\n   📊 Resumo por categoria:")
    for _, row in category_summary.iterrows():
        print(f"      {row['pathway_category']}: {row['count']} PFAMs")
    
    return df_sig

# ============================================================================
# HEATMAPS
# ============================================================================

def create_heatmap_presence_absence(merged_data, output_prefix):
    secondary_genes = {g: info for g, info in merged_data.items() if info['secondary_pathways']}
    if len(secondary_genes) < 2:
        return
    
    genes_sorted = sorted(secondary_genes.keys(), 
                         key=lambda x: (secondary_genes[x].get('pharma_score', 0), len(secondary_genes[x]['secondary_pathways'])), 
                         reverse=True)[:50]
    
    pathways_list = list(SECONDARY_PATHWAYS.keys())
    
    matrix = []
    for gene in genes_sorted:
        row = []
        for p in pathways_list:
            if p in secondary_genes[gene]['secondary_pathways']:
                score = secondary_genes[gene].get('pharma_score', 0)
                row.append(score if score > 0 else 1)
            else:
                row.append(0)
        matrix.append(row)
    
    df = pd.DataFrame(matrix, 
                      index=[f"{g[:20]} ({secondary_genes[g].get('pharma_score', 0)})" for g in genes_sorted], 
                      columns=[SECONDARY_PATHWAYS[p]['name'][:20] for p in pathways_list])
    
    fig, ax = plt.subplots(figsize=(18, 14))
    cmap_white_red = sns.light_palette("#d62728", as_cmap=True, reverse=False)
    
    sns.heatmap(df, 
                cmap=cmap_white_red,
                cbar_kws={'label': 'Pharma Score / Presence'},
                linewidths=0.5, 
                linecolor='lightgray', 
                ax=ax,
                vmin=0,
                square=False)
    
    ax.set_title('Gene-Pathway Association with Pharma Scores\n(White = low potential, Red = high potential)', 
                 fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(fontsize=7)
    plt.subplots_adjust(right=0.92)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_heatmap_scores.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_heatmap_scores.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   ✅ Heatmap com scores salvo: {output_prefix}_heatmap_scores.pdf")

def create_pharma_heatmap(merged_data, output_prefix):
    pharma_genes = {g: info for g, info in merged_data.items() if info.get('pharma_score', 0) > 0}
    if len(pharma_genes) < 2:
        return
    
    genes_sorted = sorted(pharma_genes.keys(), 
                         key=lambda x: pharma_genes[x]['pharma_score'], 
                         reverse=True)[:50]
    
    pathways_list = list(PHARMA_DETAILS.keys())
    
    matrix = []
    for gene in genes_sorted:
        row = []
        for p in pathways_list:
            if p in merged_data[gene]['secondary_pathways']:
                row.append(merged_data[gene]['pharma_score'])
            else:
                row.append(0)
        matrix.append(row)
    
    df = pd.DataFrame(matrix, 
                      index=[f"{g[:20]} ({pharma_genes[g]['pharma_score']})" for g in genes_sorted], 
                      columns=[f"{PHARMA_DETAILS[p]['name'][:20]}" for p in pathways_list])
    
    fig, ax = plt.subplots(figsize=(20, 16))
    cmap_white_red = sns.light_palette("#d62728", as_cmap=True, reverse=False)
    
    sns.heatmap(df, 
                cmap=cmap_white_red,
                cbar_kws={'label': 'Pharma Score'},
                linewidths=0.5, 
                linecolor='lightgray', 
                ax=ax,
                vmin=0,
                square=False)
    
    ax.set_title('Pharmaceutical & Nutraceutical Pathways - Gene Scores\n(White = low potential, Red = high potential)', 
                 fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(fontsize=7)
    plt.subplots_adjust(right=0.92)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_pharma_heatmap_scores.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_pharma_heatmap_scores.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   ✅ Pharma heatmap com scores salvo: {output_prefix}_pharma_heatmap_scores.pdf")

# ============================================================================
# FUNÇÕES COMPARATIVAS (MULTI-GENOME) - MANTIDAS
# ============================================================================

def create_comparative_dotplot(genomes_data, output_prefix, top_n=15):
    all_enrichment = []
    for species, merged in genomes_data.items():
        high_medium_genes = [g for g, info in merged.items() 
                            if info.get('pharma_potential') in ['HIGH', 'MEDIUM']]
        background_genes = [g for g, info in merged.items() 
                           if info.get('pharma_potential') in ['LOW', 'NONE']]
        if len(high_medium_genes) < 2 or len(background_genes) < 2:
            continue
        df_enrich = calculate_enrichment(
            high_medium_genes, background_genes, merged,
            min_count=2, min_fold=1.0, filter_generic=True
        )
        if len(df_enrich) > 0:
            df_enrich['species'] = species
            all_enrichment.append(df_enrich)
    
    if not all_enrichment:
        return
    
    df_all = pd.concat(all_enrichment, ignore_index=True)
    df_all = df_all[df_all['significant']]
    
    if len(df_all) == 0:
        return
    
    df_all.to_csv(f"{output_prefix}_comparative_combined_pfam_enrichment.csv", index=False)
    
    top_terms = df_all.sort_values('fold_enrichment', ascending=False).head(top_n)
    
    fig_height = max(6, len(top_terms) * 0.4)
    fig, ax = plt.subplots(figsize=(14, fig_height))
    
    species_list = list(genomes_data.keys())
    color_palette = plt.cm.Set3(np.linspace(0, 1, len(species_list)))
    species_colors = {s: color_palette[i] for i, s in enumerate(species_list)}
    colors = [species_colors.get(s, '#999999') for s in top_terms['species']]
    unique_species = top_terms['species'].unique()
    legend_elements = [Patch(facecolor=species_colors.get(s, '#999999'), label=s)
                       for s in unique_species if s in species_colors]
    
    ax.scatter(
        top_terms['fold_enrichment'],
        top_terms['pfam_domain'],
        s=top_terms['count_in_group'] * 15 + 30,
        c=colors,
        alpha=0.8,
        edgecolors='black',
        linewidth=0.5
    )
    
    ax.set_xlabel('Fold Enrichment', fontsize=12, fontweight='bold')
    ax.set_ylabel('PFAM Domain', fontsize=12, fontweight='bold')
    ax.set_title('Enriched PFAM Domains in HIGH+MEDIUM Genes by Species\n(Compared to LOW+NONE Background)', 
                 fontsize=14, fontweight='bold')
    ax.axvline(x=1, color='red', linestyle='--', alpha=0.5, label='Fold = 1')
    ax.legend(handles=legend_elements,
             loc='center left',
             bbox_to_anchor=(1.02, 0.5),
             fontsize=10,
             title='Species',
             title_fontsize=11)
    
    plt.subplots_adjust(right=0.78)
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(f"{output_prefix}_comparative_combined_pfam_dotplot.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_comparative_combined_pfam_dotplot.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   ✅ Comparative combined PFAM dotplot salvo: {output_prefix}_comparative_combined_pfam_dotplot.pdf")

def create_pathway_presence_matrix(genomes_data, all_pathways, output_prefix):
    species_list = list(genomes_data.keys())
    pathways_list = sorted(list(all_pathways), key=lambda x: SECONDARY_PATHWAYS.get(x, {}).get('class','Z'))
    
    matrix = []
    for pathway in pathways_list:
        row = []
        for species in species_list:
            total_score = 0
            count = 0
            for gene, info in genomes_data[species].items():
                if pathway in info['secondary_pathways']:
                    total_score += info.get('pharma_score', 0)
                    count += 1
            avg_score = total_score / count if count > 0 else 0
            row.append(avg_score)
        matrix.append(row)
    
    pathway_names = [SECONDARY_PATHWAYS.get(p, {}).get('name', p)[:35] for p in pathways_list]
    df = pd.DataFrame(matrix, index=pathway_names, columns=species_list)
    
    fig, ax = plt.subplots(figsize=(max(12, len(species_list)*2), max(8, len(pathways_list)*0.4)))
    cmap_white_red = sns.light_palette("#d62728", as_cmap=True, reverse=False)
    
    sns.heatmap(df, annot=True, fmt='.1f', cmap=cmap_white_red, linewidths=0.5, ax=ax, 
                cbar_kws={'label':'Average Pharma Score per Pathway'})
    ax.set_title('Pathway Scores Across Species\n(White = low score, Red = high score)', fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_pathway_score_matrix.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_pathway_score_matrix.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   ✅ Pathway score matrix salvo: {output_prefix}_pathway_score_matrix.pdf")

def create_comparative_bar_chart(genomes_data, output_prefix):
    species = list(genomes_data.keys())
    high = [sum(1 for g in genomes_data[s].values() if g.get('pharma_potential') == 'HIGH') for s in species]
    medium = [sum(1 for g in genomes_data[s].values() if g.get('pharma_potential') == 'MEDIUM') for s in species]
    low = [sum(1 for g in genomes_data[s].values() if g.get('pharma_potential') == 'LOW') for s in species]
    
    fig, ax = plt.subplots(figsize=(14, max(6, len(species)*0.5)))
    y = np.arange(len(species))
    ax.barh(y, high, color='#2ECC71', label='HIGH')
    ax.barh(y, medium, left=high, color='#F39C12', label='MEDIUM')
    ax.barh(y, low, left=[h+m for h,m in zip(high,medium)], color='#E74C3C', label='LOW')
    ax.set_yticks(y)
    ax.set_yticklabels(species, fontsize=11)
    ax.set_xlabel('Number of genes with pharmaceutical potential', fontsize=12, fontweight='bold')
    ax.set_title('Comparative Pharmaceutical Potential Across Species', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_comparative_bars.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_comparative_bars.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_comparative_radar_chart(genomes_data, output_prefix):
    from math import pi
    metrics = ['Percent_Pharma_of_Secondary', 'Unique_Compound_Classes', 'Avg_Pharma_Score']
    rows = []
    for species, merged in genomes_data.items():
        total = len(merged)
        secondary = sum(1 for g in merged.values() if g['secondary_pathways'])
        pharma = sum(1 for g in merged.values() if g.get('pharma_score', 0) > 0)
        avg = sum(g.get('pharma_score', 0) for g in merged.values()) / pharma if pharma > 0 else 0
        classes = set()
        for g in merged.values():
            classes.update(g.get('compound_classes', []))
        rows.append({
            'Species': species,
            'Percent_Pharma_of_Secondary': pharma/secondary*100 if secondary > 0 else 0,
            'Unique_Compound_Classes': len(classes),
            'Avg_Pharma_Score': avg
        })
    df = pd.DataFrame(rows).set_index('Species')
    df_norm = (df - df.min()) / (df.max() - df.min() + 1e-9)
    angles = [n * 2*pi / len(metrics) for n in range(len(metrics))]
    angles += angles[:1]
    
    fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(projection='polar'))
    colors = plt.cm.Set3(np.linspace(0, 1, len(df_norm)))
    for i, (species, row) in enumerate(df_norm.iterrows()):
        values = row[metrics].tolist()
        values += values[:1]
        ax.plot(angles, values, 'o-', linewidth=2, label=species, color=colors[i])
        ax.fill(angles, values, alpha=0.1, color=colors[i])
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(metrics, fontsize=10)
    ax.set_ylim(0, 1)
    ax.set_title('Multi-dimensional Comparison of Pharmaceutical Potential', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0), fontsize=9)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_comparative_radar.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_comparative_radar.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_heatmap_comparative(genomes_data, output_prefix):
    species_list = list(genomes_data.keys())
    all_classes = set()
    class_scores = defaultdict(lambda: defaultdict(list))
    for species, merged in genomes_data.items():
        for gene, info in merged.items():
            for cc in info.get('compound_classes', []):
                if cc:
                    all_classes.add(cc)
                    class_scores[species][cc].append(info.get('pharma_score', 0))
    all_classes = sorted(all_classes)
    matrix = []
    for species in species_list:
        row = []
        for cc in all_classes:
            scores = class_scores[species].get(cc, [0])
            row.append(np.mean(scores))
        matrix.append(row)
    df = pd.DataFrame(matrix, index=species_list, columns=all_classes)
    
    fig, ax = plt.subplots(figsize=(max(10, len(all_classes)*1.2), max(8, len(species_list)*0.8)))
    cmap_white_red = sns.light_palette("#d62728", as_cmap=True, reverse=False)
    
    sns.heatmap(df, annot=True, fmt='.1f', cmap=cmap_white_red, linewidths=0.5, ax=ax, 
                cbar_kws={'label':'Average Pharma Score per Compound Class'})
    ax.set_title('Average Pharmaceutical Score by Compound Class Across Species\n(White = low score, Red = high score)', fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_class_score_heatmap.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_class_score_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_species_similarity_network(genomes_data, output_prefix):
    species_list = list(genomes_data.keys())
    n = len(species_list)
    if n < 2:
        return
    
    rows = []
    for species, merged in genomes_data.items():
        total = len(merged)
        secondary = sum(1 for g in merged.values() if g['secondary_pathways'])
        pharma = sum(1 for g in merged.values() if g.get('pharma_score', 0) > 0)
        high = sum(1 for g in merged.values() if g.get('pharma_potential') == 'HIGH')
        classes = set()
        for g in merged.values():
            classes.update(g.get('compound_classes', []))
        rows.append({'Species': species, 'Total': total, 'Secondary': secondary, 'Pharma': pharma, 'HIGH': high, 'Classes': len(classes)})
    df = pd.DataFrame(rows)
    
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            diff_high = abs(df.iloc[i]['HIGH'] - df.iloc[j]['HIGH'])
            diff_pct = abs(df.iloc[i]['Pharma']/df.iloc[i]['Secondary'] - df.iloc[j]['Pharma']/df.iloc[j]['Secondary']) if df.iloc[i]['Secondary']>0 and df.iloc[j]['Secondary']>0 else 1
            diff_classes = abs(df.iloc[i]['Classes'] - df.iloc[j]['Classes'])
            sim = 1 / (1 + diff_high*0.5 + diff_pct*0.3 + diff_classes*0.2)
            if sim > 0.3:
                edges.append({'species1': species_list[i], 'species2': species_list[j], 'similarity': sim})
    if not edges:
        return
    
    angles = np.linspace(0, 2*np.pi, n, endpoint=False)
    pos = {species: (np.cos(angle), np.sin(angle)) for species, angle in zip(species_list, angles)}
    
    fig, ax = plt.subplots(figsize=(14, 12))
    
    for edge in edges:
        x1, y1 = pos[edge['species1']]
        x2, y2 = pos[edge['species2']]
        lw = 1 + edge['similarity']*5
        ax.plot([x1, x2], [y1, y2], '#3498DB', alpha=0.5, linewidth=lw, zorder=1)
    
    node_sizes = [1000 + df[df['Species']==s]['HIGH'].values[0]*80 for s in species_list]
    node_sizes = [min(s, 3000) for s in node_sizes]
    pcts = [df[df['Species']==s]['Pharma'].values[0]/df[df['Species']==s]['Secondary'].values[0]*100 if df[df['Species']==s]['Secondary'].values[0]>0 else 0 for s in species_list]
    colors = ['#2ECC71' if p>5 else '#F39C12' if p>2 else '#E74C3C' for p in pcts]
    
    x_coords = [pos[s][0] for s in species_list]
    y_coords = [pos[s][1] for s in species_list]
    ax.scatter(x_coords, y_coords, s=node_sizes, c=colors, alpha=0.85, 
               edgecolors='black', linewidth=2, zorder=2)
    
    for species in species_list:
        x, y = pos[species]
        ax.annotate(species, (x, y), xytext=(0, 0), textcoords='offset points', 
                    ha='center', va='center', fontsize=10, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85))
    
    ax.set_title('Species Similarity Network', fontsize=14, fontweight='bold')
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.axis('off')
    
    legend_elements = [Patch(facecolor='#2ECC71', label='High potential (>5%)'), 
                       Patch(facecolor='#F39C12', label='Medium (2-5%)'), 
                       Patch(facecolor='#E74C3C', label='Low (<2%)')]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_species_network.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_species_network.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_ranking_table(genomes_data, output_prefix):
    rows = []
    for species, merged in genomes_data.items():
        total = len(merged)
        secondary = sum(1 for g in merged.values() if g['secondary_pathways'])
        pharma = sum(1 for g in merged.values() if g.get('pharma_score', 0) > 0)
        high = sum(1 for g in merged.values() if g.get('pharma_potential') == 'HIGH')
        avg = sum(g.get('pharma_score', 0) for g in merged.values()) / pharma if pharma > 0 else 0
        classes = set()
        for g in merged.values():
            classes.update(g.get('compound_classes', []))
        rows.append({
            'Species': species,
            'Percent_Pharma_of_Secondary': pharma/secondary*100 if secondary > 0 else 0,
            'HIGH': high,
            'Classes': len(classes),
            'Avg_Score': avg
        })
    df = pd.DataFrame(rows)
    df.to_csv(f"{output_prefix}_ranking.csv", index=False)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes.flatten()
    metrics = ['Percent_Pharma_of_Secondary','HIGH','Classes','Avg_Score']
    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        df_sorted = df.sort_values(metric, ascending=True).tail(10)
        ax.barh(df_sorted['Species'], df_sorted[metric], color='#3498DB', alpha=0.8)
        ax.set_xlabel(metric.replace('_',' '), fontsize=12, fontweight='bold')
        ax.set_title(f'Ranking by {metric.replace("_"," ")}', fontsize=12, fontweight='bold')
        ax.grid(True, axis='x', alpha=0.3)
    plt.suptitle('Species Rankings by Pharmaceutical Potential Metrics', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"{output_prefix}_rankings.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_rankings.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_comparative_summary_table(genomes_data, output_prefix):
    rows = []
    for species, merged in genomes_data.items():
        total = len(merged)
        secondary = sum(1 for g in merged.values() if g['secondary_pathways'])
        pharma = sum(1 for g in merged.values() if g.get('pharma_score', 0) > 0)
        high = sum(1 for g in merged.values() if g.get('pharma_potential') == 'HIGH')
        medium = sum(1 for g in merged.values() if g.get('pharma_potential') == 'MEDIUM')
        low = sum(1 for g in merged.values() if g.get('pharma_potential') == 'LOW')
        total_score = sum(g.get('pharma_score', 0) for g in merged.values())
        avg_score = total_score / pharma if pharma > 0 else 0
        compound_classes = set()
        for g in merged.values():
            compound_classes.update(g.get('compound_classes', []))
        rows.append({
            'Species': species, 
            'Total_Genes': total,
            'Secondary_Genes': secondary,
            'Pharma_Genes': pharma,
            'Percent_Pharma_of_Secondary': round(pharma/secondary*100, 2) if secondary > 0 else 0,
            'Percent_Pharma_of_Total': round(pharma/total*100, 2) if total > 0 else 0,
            'HIGH': high, 
            'MEDIUM': medium, 
            'LOW': low,
            'Avg_Pharma_Score': round(avg_score, 2),
            'Unique_Compound_Classes': len(compound_classes)
        })
    df = pd.DataFrame(rows)
    df.to_csv(f"{output_prefix}_comparative_summary.csv", index=False)
    try:
        df.to_excel(f"{output_prefix}_comparative_summary.xlsx", index=False)
    except:
        pass
    return df

def generate_comparative_report(genomes_data, output_prefix):
    rows = []
    for species, merged in genomes_data.items():
        total = len(merged)
        secondary = sum(1 for g in merged.values() if g['secondary_pathways'])
        pharma = sum(1 for g in merged.values() if g.get('pharma_score', 0) > 0)
        high = sum(1 for g in merged.values() if g.get('pharma_potential') == 'HIGH')
        medium = sum(1 for g in merged.values() if g.get('pharma_potential') == 'MEDIUM')
        low = sum(1 for g in merged.values() if g.get('pharma_potential') == 'LOW')
        classes = set()
        for g in merged.values():
            classes.update(g.get('compound_classes', []))
        rows.append({
            'Species': species,
            'Total': total,
            'Secondary': secondary,
            'Pharma': pharma,
            'Percent_Pharma_of_Secondary': pharma/secondary*100 if secondary > 0 else 0,
            'HIGH': high,
            'MEDIUM': medium,
            'LOW': low,
            'Classes': len(classes)
        })
    df = pd.DataFrame(rows).sort_values('Percent_Pharma_of_Secondary', ascending=False)
    
    with open(f"{output_prefix}_comparative_report.txt", 'w') as f:
        f.write("COMPARATIVE ANALYSIS OF PHARMACEUTICAL POTENTIAL ACROSS GENOMES\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Number of species: {len(df)}\n\n")
        f.write("SPECIES RANKING (by % pharmaceutical genes relative to secondary):\n")
        for i, (_, row) in enumerate(df.iterrows(), 1):
            f.write(f"{i}. {row['Species']}: {row['Percent_Pharma_of_Secondary']:.2f}% (HIGH:{row['HIGH']}, MEDIUM:{row['MEDIUM']}, LOW:{row['LOW']})\n")
        best = df.iloc[0]
        f.write(f"\nBEST PERFORMER: {best['Species']}\n")
        f.write(f"  • {best['Percent_Pharma_of_Secondary']:.2f}% of secondary genes have potential\n")
        f.write(f"  • {best['HIGH']} HIGH potential genes\n")
        f.write(f"  • {best['Classes']} compound classes\n")


# ============================================================================
# FUNÇÃO AUXILIAR PARA PROCESSAR GFF NAS ANÁLISES
# ============================================================================

def process_gff_for_analysis(gff_file, genome_file, hmm_dir, output_dir, prefix,
                            domains_dict, merge_func, *args, **kwargs):
    """
    Função genérica para processar GFF em qualquer análise.
    domains_dict: Dicionário de domínios HMM para busca
    merge_func: Função para integrar resultados
    """
    if not gff_file or not os.path.exists(gff_file):
        print("  ⚠️ GFF não fornecido ou não encontrado. Pulando busca HMMER.")
        return None
    
    if not genome_file or not os.path.exists(genome_file):
        print("  ⚠️ Arquivo do genoma não fornecido ou não encontrado. Pulando busca HMMER.")
        return None
    
    print(f"\n🔬 Processando GFF para {prefix}...")
    
    # Extrair proteínas
    protein_fasta = os.path.join(output_dir, f"{prefix}_proteins.fasta")
    if not extract_proteins_from_gff(gff_file, genome_file, protein_fasta):
        print("  ⚠️ Falha na extração de proteínas. Pulando busca HMMER.")
        return None
    
    # Executar HMMER
    hmmer_dir = hmm_dir or os.path.join(os.path.dirname(__file__), 'hmmer_profiles')
    hmmer_results = run_hmmer_search_generic(protein_fasta, hmmer_dir, output_dir, 
                                            domains_dict, evalue=1e-10, prefix=prefix)
    
    if not hmmer_results:
        print("  ⚠️ Nenhum domínio HMMER encontrado.")
        return None
    
    print(f"  ✅ HMMER identificou {len(hmmer_results)} genes com domínios relevantes")
    
    # Extrair genes do GFF
    gff_genes = set()
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 9 and fields[2] == 'gene':
                attrs = fields[8]
                gene_id = None
                for item in attrs.split(';'):
                    if item.startswith('ID='):
                        gene_id = item[3:]
                        break
                    elif item.startswith('gene_id='):
                        gene_id = item[8:]
                        break
                if gene_id:
                    gff_genes.add(gene_id)
    
    # Integrar resultados
    merged_data = merge_func(hmmer_results, gff_genes, *args, **kwargs)
    return merged_data

# ============================================================================
# FUNÇÕES PRINCIPAIS COM SUPORTE A GFF
# ============================================================================

def run_secondary_analysis(eggnog_file, interpro_file, gff_file=None, genome_file=None, hmm_dir=None, 
                           output_dir='secondary_metabolism', prefix='secondary'):
    print(f"\n{'='*80}")
    print("SECONDARY METABOLISM ANALYSIS v3.5 (with GFF support)")
    print(f"{'='*80}")
    print(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    os.makedirs(output_dir, exist_ok=True)
    output_prefix = os.path.join(output_dir, prefix)
    
    # Parse anotações
    eggnog_data = parse_eggnog_universal(eggnog_file) if os.path.exists(eggnog_file) else None
    interpro_data = parse_interpro_universal(interpro_file) if os.path.exists(interpro_file) else None
    
    if eggnog_data is None:
        print("  ⚠️ eggNOG não encontrado. Usando apenas InterPro/HMMER.")
        eggnog_data = {'genes': {}, 'kegg_pathways': defaultdict(set), 'ko_numbers': defaultdict(set)}
    if interpro_data is None:
        print("  ⚠️ InterPro não encontrado. Usando apenas eggNOG/HMMER.")
        interpro_data = {'genes': set(), 'pfam_domains': defaultdict(set)}
    
    merged_data = merge_annotation_data(eggnog_data, interpro_data, SECONDARY_PATHWAYS)
    print(f"\n📊 Genes com anotação funcional: {len(merged_data):,}")
    
    # Processar GFF se fornecido
    if gff_file and genome_file:
        print("\n🔬 Executando busca HMMER no GFF...")
        hmmer_results = process_gff_for_analysis(
            gff_file, genome_file, hmm_dir, output_dir, prefix,
            SECONDARY_DOMAINS_HMM,
            lambda hmmer_results, gff_genes: integrate_hmmer_results_generic(
                eggnog_data, interpro_data, hmmer_results, gff_genes, 'secondary'
            )
        )
        if hmmer_results:
            merged_data = hmmer_results
            secondary_count = sum(1 for g in merged_data.values() if g['secondary_pathways'])
            print(f"  ✅ Genes com metabolismo secundário (após HMMER): {secondary_count}")
    
    # Gerar visualizações
    create_ko_dotplot(merged_data, output_prefix, top_n=15)
    plot_top_kos_by_category(merged_data, output_prefix)
    plot_pathway_completeness(merged_data, output_prefix)
    plot_pie_chart_secondary_classes(merged_data, output_prefix)
    plot_upset_venn(merged_data, output_prefix)
    plot_top_ec_numbers(merged_data, output_prefix)
    plot_top_pfam_domains(merged_data, output_prefix)
    generate_publication_table_secondary(merged_data, output_prefix)
    analyze_pathway_completeness(merged_data, output_prefix)
    
    secondary_count = sum(1 for g in merged_data.values() if g['secondary_pathways'])
    print(f"\n📊 Estatísticas de metabolismo secundário:")
    print(f"   Genes com metabolismo secundário: {secondary_count}")
    print(f"   Proporção do total: {secondary_count/len(merged_data)*100:.1f}%")
    
    print(f"\n✅ Secondary metabolism analysis completed. Outputs in: {output_dir}")
    return True

def run_pharma_analysis(eggnog_file, interpro_file, gff_file=None, genome_file=None, hmm_dir=None,
                        output_dir='pharma_analysis', prefix='pharma', 
                        threshold=4, network_nodes=15, network_method='enrichment'):
    print(f"\n{'='*80}")
    print("PHARMACEUTICAL SECONDARY METABOLISM ANALYSIS v3.5 (with GFF support)")
    print(f"{'='*80}")
    print(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    os.makedirs(output_dir, exist_ok=True)
    output_prefix = os.path.join(output_dir, prefix)
    
    # Parse anotações
    eggnog_data = parse_eggnog_universal(eggnog_file) if os.path.exists(eggnog_file) else None
    interpro_data = parse_interpro_universal(interpro_file) if os.path.exists(interpro_file) else None
    
    if eggnog_data is None:
        print("  ⚠️ eggNOG não encontrado. Usando apenas InterPro/HMMER.")
        eggnog_data = {'genes': {}, 'kegg_pathways': defaultdict(set), 'ko_numbers': defaultdict(set)}
    if interpro_data is None:
        print("  ⚠️ InterPro não encontrado. Usando apenas eggNOG/HMMER.")
        interpro_data = {'genes': set(), 'pfam_domains': defaultdict(set)}
    
    merged_data = merge_pharma_data(eggnog_data, interpro_data, threshold)
    print(f"\n📊 Genes com anotação funcional: {len(merged_data):,}")
    
    # Processar GFF com domínios farmacológicos
    if gff_file and genome_file:
        print("\n🔬 Executando busca HMMER para domínios farmacológicos...")
        hmmer_results = process_gff_for_analysis(
            gff_file, genome_file, hmm_dir, output_dir, prefix,
            PHARMA_DOMAINS_HMM,
            lambda hmmer_results, gff_genes: integrate_hmmer_results_pharma(
                eggnog_data, interpro_data, hmmer_results, gff_genes, threshold
            )
        )
        if hmmer_results:
            merged_data = hmmer_results
    
    validate_pharma_groups(merged_data)
    
    secondary_count = sum(1 for g in merged_data.values() if g['secondary_pathways'])
    pharma_count = sum(1 for g in merged_data.values() if g.get('pharma_score', 0) > 0)
    
    print(f"\n📊 Verificação de consistência:")
    print(f"   Genes com metabolismo secundário: {secondary_count}")
    print(f"   Genes com potencial farmacológico: {pharma_count}")
    print(f"   Proporção farmacológico/secundário: {pharma_count/secondary_count*100:.1f}%")
    
    if pharma_count > secondary_count:
        print("   ⚠️ ATENÇÃO: Pharma > Secondary - isso NÃO deveria acontecer!")
    else:
        print("   ✅ CONSISTENTE: Pharma é subconjunto do Secondary")
    
    # Gerar visualizações
    create_heatmap_presence_absence(merged_data, output_prefix)
    create_pharma_heatmap(merged_data, output_prefix)
    plot_pharma_score_ranking(merged_data, output_prefix)
    plot_compound_class_distribution(merged_data, output_prefix)
    plot_pharmaceutical_uses(merged_data, output_prefix)
    plot_sunburst_hierarchical(merged_data, output_prefix)
    plot_donut_chart_potential(merged_data, output_prefix)
    build_functional_network(merged_data, output_prefix, min_edge_weight=1, max_nodes=network_nodes, selection_method=network_method)
    predict_novel_compounds(merged_data, output_prefix)
    generate_publication_table_pharma(merged_data, output_prefix)
    analyze_pathway_completeness(merged_data, output_prefix)
    
    print("\n📊 Gerando dotplots de enriquecimento combinado...")
    create_pfam_dotplot_combined(merged_data, output_prefix, top_n=15)
    
    print(f"\n✅ Pharma analysis completed. Outputs in: {output_dir}")
    return True

def integrate_hmmer_results_pharma(eggnog_data, interpro_data, hmmer_results, gff_genes, pharma_threshold=4):
    """Integra resultados HMMER específicos para análise farmacológica."""
    merged_temp = merge_annotation_data(eggnog_data, interpro_data, SECONDARY_PATHWAYS)
    
    genes_added = 0
    genes_updated = 0
    
    for gene, hmmer_info in hmmer_results.items():
        if gene in merged_temp:
            # Atualizar gene existente
            if 'hmmer_domains' not in merged_temp[gene]:
                merged_temp[gene]['hmmer_domains'] = []
                merged_temp[gene]['hmmer_classes'] = set()
            merged_temp[gene]['hmmer_domains'].extend(hmmer_info['domains'])
            merged_temp[gene]['hmmer_classes'].update(hmmer_info['classes'])
            
            for pfam in hmmer_info['pfams']:
                merged_temp[gene]['pfam_domains'].add(pfam)
            
            for class_name in hmmer_info['classes']:
                for pathway, p_info in SECONDARY_PATHWAYS.items():
                    if p_info['class'] == class_name:
                        merged_temp[gene]['secondary_pathways'].add(pathway)
                        merged_temp[gene]['secondary_class'].add(class_name)
                        break
            
            if 'source' not in merged_temp[gene]:
                merged_temp[gene]['source'] = []
            if 'HMMER' not in merged_temp[gene]['source']:
                merged_temp[gene]['source'].append('HMMER')
            genes_updated += 1
            
        elif gene in gff_genes or True:
            merged_temp[gene] = {
                'source': ['HMMER'],
                'description': f'HMMER: {", ".join(hmmer_info["classes"])}',
                'cog_category': '-',
                'go_terms': set(),
                'pfam_domains': set(hmmer_info['pfams']),
                'kegg_pathways': set(),
                'ec_numbers': set(),
                'ko_numbers': set(),
                'secondary_pathways': set(),
                'secondary_class': set(),
                'hmmer_domains': hmmer_info['domains'],
                'hmmer_classes': hmmer_info['classes']
            }
            
            for class_name in hmmer_info['classes']:
                for pathway, p_info in SECONDARY_PATHWAYS.items():
                    if p_info['class'] == class_name:
                        merged_temp[gene]['secondary_pathways'].add(pathway)
                        merged_temp[gene]['secondary_class'].add(class_name)
                        break
            
            genes_added += 1
    
    if genes_added > 0:
        print(f"  ✅ {genes_added} genes adicionados via HMMER (pharma)")
    if genes_updated > 0:
        print(f"  ✅ {genes_updated} genes existentes atualizados via HMMER (pharma)")
    
    # Recalcular pontuações farmacológicas
    pharma_scores = calculate_pharma_potential(merged_temp, pharma_threshold)
    
    merged = {}
    for gene in merged_temp:
        merged[gene] = merged_temp[gene].copy()
        merged[gene]['pharma_potential'] = pharma_scores.get(gene, {}).get('potential', 'NONE')
        merged[gene]['pharma_score'] = pharma_scores.get(gene, {}).get('score', 0)
        merged[gene]['pharma_pathways'] = pharma_scores.get(gene, {}).get('pharma_pathways', [])
        merged[gene]['compound_classes'] = pharma_scores.get(gene, {}).get('compound_classes', [])
        merged[gene]['pharma_uses'] = pharma_scores.get(gene, {}).get('pharma_uses', [])
    
    return merged

def run_disease_resistance_analysis(eggnog_file, interpro_file, gff_file=None, genome_file=None, hmm_dir=None,
                                    output_dir='resistance_analysis', prefix='resistance',
                                    threshold=5, network_nodes=15, network_method='enrichment'):
    print(f"\n{'='*80}")
    print("DISEASE RESISTANCE GENE ANALYSIS v3.5 (with GFF support)")
    print(f"{'='*80}")
    print(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    os.makedirs(output_dir, exist_ok=True)
    output_prefix = os.path.join(output_dir, prefix)
    
    # Parse anotações
    eggnog_data = parse_eggnog_universal(eggnog_file) if os.path.exists(eggnog_file) else None
    interpro_data = parse_interpro_universal(interpro_file) if os.path.exists(interpro_file) else None
    
    if eggnog_data is None:
        print("  ⚠️ eggNOG não encontrado. Usando apenas InterPro/HMMER.")
        eggnog_data = {'genes': {}, 'kegg_pathways': defaultdict(set), 'ko_numbers': defaultdict(set)}
    if interpro_data is None:
        print("  ⚠️ InterPro não encontrado. Usando apenas eggNOG/HMMER.")
        interpro_data = {'genes': set(), 'pfam_domains': defaultdict(set)}
    
    merged_data = merge_disease_resistance_data(eggnog_data, interpro_data, threshold)
    print(f"\n📊 Genes com anotação funcional: {len(merged_data):,}")
    
    # Processar GFF com domínios de resistência
    if gff_file and genome_file:
        print("\n🔬 Executando busca HMMER para domínios de resistência (TNJ, NLR, etc)...")
        hmmer_results = process_gff_for_analysis(
            gff_file, genome_file, hmm_dir, output_dir, prefix,
            RESISTANCE_DOMAINS_HMM,
            lambda hmmer_results, gff_genes: integrate_hmmer_results_resistance(
                eggnog_data, interpro_data, hmmer_results, gff_genes, threshold
            )
        )
        if hmmer_results:
            merged_data = hmmer_results
    
    # Validar grupos
    high = set([g for g, info in merged_data.items() if info.get('resistance_potential') == 'HIGH'])
    medium = set([g for g, info in merged_data.items() if info.get('resistance_potential') == 'MEDIUM'])
    low = set([g for g, info in merged_data.items() if info.get('resistance_potential') == 'LOW'])
    tnj_genes = [g for g, info in merged_data.items() if info.get('tnj_detected', False)]
    
    print("\n🔍 VALIDAÇÃO DE GRUPOS DE RESISTÊNCIA:")
    print(f"   HIGH: {len(high)} genes")
    print(f"   MEDIUM: {len(medium)} genes")
    print(f"   LOW: {len(low)} genes")
    print(f"   Total genes com potencial de resistência: {len(high)+len(medium)+len(low)}")
    print(f"   🔬 TNJ genes detectados: {len(tnj_genes)}")
    
    overlap_hm = high & medium
    overlap_hl = high & low
    overlap_ml = medium & low
    
    if overlap_hm or overlap_hl or overlap_ml:
        print("   ⚠️ Sobreposições encontradas (prioridade: HIGH > MEDIUM > LOW)")
        for gene in overlap_hm:
            merged_data[gene]['resistance_potential'] = 'HIGH'
        for gene in overlap_hl:
            merged_data[gene]['resistance_potential'] = 'HIGH'
        for gene in overlap_ml:
            merged_data[gene]['resistance_potential'] = 'MEDIUM'
    else:
        print("   ✅ Grupos são mutuamente exclusivos")
    
    # Gerar visualizações
    plot_resistance_score_ranking(merged_data, output_prefix)
    plot_resistance_class_distribution(merged_data, output_prefix)
    plot_resistance_pathway_completeness(merged_data, output_prefix)
    plot_resistance_donut_chart(merged_data, output_prefix)
    
    print(f"\n🌐 Construindo rede funcional de genes de resistência (max_nodes={network_nodes})...")
    build_resistance_network(merged_data, output_prefix, min_edge_weight=1, 
                            max_nodes=network_nodes, selection_method=network_method)
    
    print("\n🔬 Analisando genes TNJ (TIR-NBS-Jacalin)...")
    plot_tnj_analysis(merged_data, output_prefix)
    
    generate_publication_table_resistance(merged_data, output_prefix)
    
    print("\n📊 Gerando dotplot de enriquecimento para genes de resistência...")
    create_resistance_pfam_dotplot(merged_data, output_prefix, top_n=15)
    
    resistance_count = sum(1 for g in merged_data.values() if g.get('resistance_score', 0) > 0)
    print(f"\n📊 Estatísticas de resistência a doenças:")
    print(f"   Genes com potencial de resistência: {resistance_count}")
    print(f"   Proporção do total: {resistance_count/len(merged_data)*100:.1f}%")
    print(f"   Distribuição: HIGH={len(high)}, MEDIUM={len(medium)}, LOW={len(low)}")
    print(f"   🔬 TNJ genes: {len(tnj_genes)}")
    
    print(f"\n✅ Disease resistance analysis completed. Outputs in: {output_dir}")
    return True

def integrate_hmmer_results_resistance(eggnog_data, interpro_data, hmmer_results, gff_genes, resistance_threshold=5):
    """Integra resultados HMMER específicos para análise de resistência."""
    merged_temp = merge_annotation_data(eggnog_data, interpro_data, SECONDARY_PATHWAYS)
    
    genes_added = 0
    genes_updated = 0
    
    for gene, hmmer_info in hmmer_results.items():
        if gene in merged_temp:
            # Atualizar gene existente
            if 'hmmer_domains' not in merged_temp[gene]:
                merged_temp[gene]['hmmer_domains'] = []
                merged_temp[gene]['hmmer_classes'] = set()
            merged_temp[gene]['hmmer_domains'].extend(hmmer_info['domains'])
            merged_temp[gene]['hmmer_classes'].update(hmmer_info['classes'])
            
            for pfam in hmmer_info['pfams']:
                merged_temp[gene]['pfam_domains'].add(pfam)
            
            if 'source' not in merged_temp[gene]:
                merged_temp[gene]['source'] = []
            if 'HMMER' not in merged_temp[gene]['source']:
                merged_temp[gene]['source'].append('HMMER')
            genes_updated += 1
            
        elif gene in gff_genes or True:
            merged_temp[gene] = {
                'source': ['HMMER'],
                'description': f'HMMER: {", ".join(hmmer_info["classes"])}',
                'cog_category': '-',
                'go_terms': set(),
                'pfam_domains': set(hmmer_info['pfams']),
                'kegg_pathways': set(),
                'ec_numbers': set(),
                'ko_numbers': set(),
                'secondary_pathways': set(),
                'secondary_class': set(),
                'hmmer_domains': hmmer_info['domains'],
                'hmmer_classes': hmmer_info['classes']
            }
            genes_added += 1
    
    if genes_added > 0:
        print(f"  ✅ {genes_added} genes adicionados via HMMER (resistance)")
    if genes_updated > 0:
        print(f"  ✅ {genes_updated} genes existentes atualizados via HMMER (resistance)")
    
    # Recalcular pontuações de resistência
    resistance_scores = calculate_disease_resistance_score(merged_temp, resistance_threshold)
    
    merged = {}
    for gene in merged_temp:
        merged[gene] = merged_temp[gene].copy()
        merged[gene]['resistance_potential'] = resistance_scores.get(gene, {}).get('potential', 'NONE')
        merged[gene]['resistance_score'] = resistance_scores.get(gene, {}).get('score', 0)
        merged[gene]['resistance_classes'] = resistance_scores.get(gene, {}).get('resistance_classes', [])
        merged[gene]['resistance_pathways'] = resistance_scores.get(gene, {}).get('resistance_pathways', [])
        merged[gene]['resistance_pfams'] = resistance_scores.get(gene, {}).get('resistance_pfams', [])
        merged[gene]['tnj_detected'] = resistance_scores.get(gene, {}).get('tnj_detected', False)
        merged[gene]['tnj_architecture'] = resistance_scores.get(gene, {}).get('tnj_architecture', '')
        merged[gene]['tnj_confidence'] = resistance_scores.get(gene, {}).get('tnj_confidence', 'NONE')
    
    return merged

# ============================================================================
# FUNÇÕES AUXILIARES (EGGNOG2KEGG)
# ============================================================================

def run_eggnog2kegg(eggnog_file, output_file, output_dir='eggnog2kegg_output', prefix='kegg'):
    print(f"\n{'='*80}")
    print("EGGNOG TO KEGG CONVERTER")
    print(f"{'='*80}")
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, output_file)
    parse_eggnog_annotation(eggnog_file, out_path)
    print(f"\n✅ Eggnog to KEGG conversion completed. Output: {out_path}")
    return True

def parse_eggnog_annotation(eggnog_file, output_file):
    results = []
    total = 0
    with_ko = 0
    with open(eggnog_file, 'r') as f:
        for line in f:
            if line.startswith('##') or line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            query = fields[0]
            kegg_ko = fields[11] if len(fields) > 11 else ''
            total += 1
            if kegg_ko and kegg_ko != '-' and kegg_ko != 'ko:' and kegg_ko != '':
                matches = re.findall(r'ko:?(K\d{5})', kegg_ko)
                if matches:
                    for ko in matches:
                        results.append(f"{query}\t{ko}")
                        with_ko += 1
                else:
                    kos = re.findall(r'K\d{5}', kegg_ko)
                    for ko in kos:
                        results.append(f"{query}\t{ko}")
                        with_ko += 1
            else:
                results.append(query)
    with open(output_file, 'w') as f:
        for line in results:
            f.write(line + '\n')
    print(f"Total genes: {total}, Genes with KO: {with_ko} ({with_ko/total*100:.2f}%)")

# ============================================================================
# INTERFACE GRÁFICA
# ============================================================================

class MetabAnalysisApp:
    def __init__(self, root):
        self.root = root
        self.root.title("annoMining - v3.5")
        self.root.geometry("1050x950")
        self.root.resizable(True, True)
        self.bg_color = "#ffffff"
        self.fg_color = "#000000"
        self.accent_color = "#2563eb"
        self.button_bg = "#2563eb"
        self.button_fg = "#ffffff"
        self.entry_bg = "#f3f4f6"
        self.tab_bg = "#ffffff"
        style = ttk.Style()
        style.theme_use('clam')
        style.configure('TNotebook', background='#ffffff', borderwidth=0)
        style.configure('TNotebook.Tab', background='#f3f4f6', foreground='#1f2937',
                        padding=[15, 8], font=('Segoe UI', 11, 'bold'))
        style.map('TNotebook.Tab', background=[('selected', '#2563eb')],
                  foreground=[('selected', 'white')])
        style.configure('TFrame', background='#ffffff')
        style.configure('TLabelframe', background='#ffffff', foreground='#1f2937')
        style.configure('TLabelframe.Label', background='#ffffff', foreground='#1f2937')
        main_frame = ttk.Frame(root, style='TFrame')
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        title = tk.Label(main_frame, text="annoMining v3.5",
                         font=('Segoe UI', 18, 'bold'), bg='#ffffff', fg='#1f2937')
        title.pack(pady=10)
        subtitle = tk.Label(main_frame, text="Secondary Metabolism, Pharmaceutical & Disease Resistance Analysis\n🔬 with Full GFF/HMMER Support",
                            font=('Segoe UI', 10, 'bold'), bg='#ffffff', fg='#334e76')
        subtitle.pack(pady=(0,5))
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        self.create_secondary_tab()
        self.create_pharma_tab()
        self.create_resistance_tab()
        self.create_eggnog_tab()
        footer = tk.Label(main_frame, text="Version 3.5 - Full GFF Integration",
                          font=('Segoe UI', 9), bg='#ffffff', fg='#6b7280')
        footer.pack(side=tk.BOTTOM, pady=5)

    def create_secondary_tab(self):
        tab = tk.Frame(self.notebook, bg='#ffffff')
        self.notebook.add(tab, text="Secondary Metabolism")
        container = tk.Frame(tab, bg='#ffffff', padx=20, pady=20)
        container.pack(fill=tk.BOTH, expand=True)
        tk.Label(container, text="Secondary Metabolism Analysis (with GFF/HMMER support)",
                 font=('Segoe UI', 14, 'bold'), bg='#ffffff', fg='#1f2937').grid(row=0, column=0, columnspan=3, pady=(0, 20), sticky='w')
        tk.Label(container, text="eggNOG annotations file (*.emapper.annotations):",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=1, column=0, sticky='e', padx=5, pady=4)
        self.secondary_eggnog = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                        relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.secondary_eggnog.grid(row=1, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_file(self.secondary_eggnog),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=1, column=2, padx=5)
        tk.Label(container, text="InterProScan file (*.tsv):",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=2, column=0, sticky='e', padx=5, pady=4)
        self.secondary_interpro = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                          relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.secondary_interpro.grid(row=2, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_file(self.secondary_interpro),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=2, column=2, padx=5)
        tk.Label(container, text="Braker GFF file (filtered):", font=('Segoe UI', 10),
                 bg='#ffffff', fg='#1f2937').grid(row=3, column=0, sticky='e', padx=5, pady=4)
        self.secondary_gff = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                      relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.secondary_gff.grid(row=3, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_file(self.secondary_gff),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=3, column=2, padx=5)
        tk.Label(container, text="Genome FASTA file:", font=('Segoe UI', 10),
                 bg='#ffffff', fg='#1f2937').grid(row=4, column=0, sticky='e', padx=5, pady=4)
        self.secondary_genome = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                         relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.secondary_genome.grid(row=4, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_file(self.secondary_genome),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=4, column=2, padx=5)
        tk.Label(container, text="HMM profiles directory (optional):",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=5, column=0, sticky='e', padx=5, pady=4)
        self.secondary_hmm = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                      relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.secondary_hmm.grid(row=5, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_dir(self.secondary_hmm),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=5, column=2, padx=5)
        tk.Label(container, text="Output folder:",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=6, column=0, sticky='e', padx=5, pady=4)
        self.secondary_outdir = tk.Entry(container, width=30, bg='#f3f4f6', fg='#1f2937',
                                        relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.secondary_outdir.insert(0, "secondary_metabolism")
        self.secondary_outdir.grid(row=6, column=1, sticky='w', padx=5, pady=4)
        tk.Button(container, text="▶ Run Analysis", command=self.run_secondary,
                  bg='#2563eb', fg='#ffffff', font=('Segoe UI', 12, 'bold'),
                  relief='flat', padx=30, pady=8).grid(row=7, column=0, columnspan=3, pady=15)
        self.secondary_log = scrolledtext.ScrolledText(container, height=8, bg='#1e293b', fg='#e2e8f0',
                                                       font=('Consolas', 9), relief='flat')
        self.secondary_log.grid(row=8, column=0, columnspan=3, padx=10, pady=10, sticky='we')

    def create_pharma_tab(self):
        tab = tk.Frame(self.notebook, bg='#ffffff')
        self.notebook.add(tab, text="Pharma Analysis")
        container = tk.Frame(tab, bg='#ffffff', padx=20, pady=20)
        container.pack(fill=tk.BOTH, expand=True)
        tk.Label(container, text="Pharmaceutical Analysis (with GFF/HMMER support)",
                 font=('Segoe UI', 14, 'bold'), bg='#ffffff', fg='#1f2937').grid(row=0, column=0, columnspan=3, pady=(0, 20), sticky='w')
        tk.Label(container, text="eggNOG annotations file:",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=1, column=0, sticky='e', padx=5, pady=4)
        self.pharma_eggnog = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                      relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.pharma_eggnog.grid(row=1, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_file(self.pharma_eggnog),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=1, column=2, padx=5)
        tk.Label(container, text="InterProScan file (optional):",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=2, column=0, sticky='e', padx=5, pady=4)
        self.pharma_interpro = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                        relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.pharma_interpro.grid(row=2, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_file(self.pharma_interpro),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=2, column=2, padx=5)
        tk.Label(container, text="Braker GFF file (filtered):", font=('Segoe UI', 10),
                 bg='#ffffff', fg='#1f2937').grid(row=3, column=0, sticky='e', padx=5, pady=4)
        self.pharma_gff = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                   relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.pharma_gff.grid(row=3, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_file(self.pharma_gff),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=3, column=2, padx=5)
        tk.Label(container, text="Genome FASTA file:", font=('Segoe UI', 10),
                 bg='#ffffff', fg='#1f2937').grid(row=4, column=0, sticky='e', padx=5, pady=4)
        self.pharma_genome = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                      relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.pharma_genome.grid(row=4, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_file(self.pharma_genome),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=4, column=2, padx=5)
        tk.Label(container, text="HMM profiles directory (optional):",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=5, column=0, sticky='e', padx=5, pady=4)
        self.pharma_hmm = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                   relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.pharma_hmm.grid(row=5, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_dir(self.pharma_hmm),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=5, column=2, padx=5)
        tk.Label(container, text="Threshold score:",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=6, column=0, sticky='e', padx=5, pady=4)
        self.pharma_threshold = tk.Entry(container, width=10, bg='#f3f4f6', fg='#1f2937',
                                         relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.pharma_threshold.insert(0, "4")
        self.pharma_threshold.grid(row=6, column=1, sticky='w', padx=5, pady=4)
        tk.Label(container, text="Network nodes:",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=7, column=0, sticky='e', padx=5, pady=4)
        self.pharma_nodes = tk.Entry(container, width=10, bg='#f3f4f6', fg='#1f2937',
                                     relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.pharma_nodes.insert(0, "15")
        self.pharma_nodes.grid(row=7, column=1, sticky='w', padx=5, pady=4)
        tk.Label(container, text="Network method:",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=8, column=0, sticky='e', padx=5, pady=4)
        self.pharma_method = ttk.Combobox(container, values=['score', 'enrichment', 'balanced'], width=15)
        self.pharma_method.set('enrichment')
        self.pharma_method.grid(row=8, column=1, sticky='w', padx=5, pady=4)
        tk.Label(container, text="Output folder:",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=9, column=0, sticky='e', padx=5, pady=4)
        self.pharma_outdir = tk.Entry(container, width=30, bg='#f3f4f6', fg='#1f2937',
                                      relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.pharma_outdir.insert(0, "pharma_analysis")
        self.pharma_outdir.grid(row=9, column=1, sticky='w', padx=5, pady=4)
        tk.Button(container, text="▶ Run Analysis", command=self.run_pharma,
                  bg='#2563eb', fg='#ffffff', font=('Segoe UI', 12, 'bold'),
                  relief='flat', padx=30, pady=8).grid(row=10, column=0, columnspan=3, pady=15)
        self.pharma_log = scrolledtext.ScrolledText(container, height=8, bg='#1e293b', fg='#e2e8f0',
                                                    font=('Consolas', 9), relief='flat')
        self.pharma_log.grid(row=11, column=0, columnspan=3, padx=10, pady=10, sticky='we')

    def create_resistance_tab(self):
        tab = tk.Frame(self.notebook, bg='#ffffff')
        self.notebook.add(tab, text="Disease Resistance")
        container = tk.Frame(tab, bg='#ffffff', padx=20, pady=20)
        container.pack(fill=tk.BOTH, expand=True)
        tk.Label(container, text="Disease Resistance Analysis (with GFF/HMMER support)",
                 font=('Segoe UI', 14, 'bold'), bg='#ffffff', fg='#1f2937').grid(row=0, column=0, columnspan=3, pady=(0, 20), sticky='w')
        tk.Label(container, text="eggNOG annotations file:",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=1, column=0, sticky='e', padx=5, pady=4)
        self.resistance_eggnog = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                          relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.resistance_eggnog.grid(row=1, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_file(self.resistance_eggnog),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=1, column=2, padx=5)
        tk.Label(container, text="InterProScan file (optional):",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=2, column=0, sticky='e', padx=5, pady=4)
        self.resistance_interpro = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                            relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.resistance_interpro.grid(row=2, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_file(self.resistance_interpro),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=2, column=2, padx=5)
        tk.Label(container, text="Braker GFF file (filtered):", font=('Segoe UI', 10),
                 bg='#ffffff', fg='#1f2937').grid(row=3, column=0, sticky='e', padx=5, pady=4)
        self.resistance_gff = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                       relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.resistance_gff.grid(row=3, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_file(self.resistance_gff),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=3, column=2, padx=5)
        tk.Label(container, text="Genome FASTA file:", font=('Segoe UI', 10),
                 bg='#ffffff', fg='#1f2937').grid(row=4, column=0, sticky='e', padx=5, pady=4)
        self.resistance_genome = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                          relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.resistance_genome.grid(row=4, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_file(self.resistance_genome),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=4, column=2, padx=5)
        tk.Label(container, text="HMM profiles directory (optional):",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=5, column=0, sticky='e', padx=5, pady=4)
        self.resistance_hmm = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                       relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.resistance_hmm.grid(row=5, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_dir(self.resistance_hmm),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=5, column=2, padx=5)
        tk.Label(container, text="Threshold score:",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=6, column=0, sticky='e', padx=5, pady=4)
        self.resistance_threshold = tk.Entry(container, width=10, bg='#f3f4f6', fg='#1f2937',
                                             relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.resistance_threshold.insert(0, "5")
        self.resistance_threshold.grid(row=6, column=1, sticky='w', padx=5, pady=4)
        tk.Label(container, text="Network nodes:",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=7, column=0, sticky='e', padx=5, pady=4)
        self.resistance_nodes = tk.Entry(container, width=10, bg='#f3f4f6', fg='#1f2937',
                                         relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.resistance_nodes.insert(0, "15")
        self.resistance_nodes.grid(row=7, column=1, sticky='w', padx=5, pady=4)
        tk.Label(container, text="Network method:",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=8, column=0, sticky='e', padx=5, pady=4)
        self.resistance_method = ttk.Combobox(container, values=['score', 'enrichment', 'balanced'], width=15)
        self.resistance_method.set('enrichment')
        self.resistance_method.grid(row=8, column=1, sticky='w', padx=5, pady=4)
        tk.Label(container, text="Output folder:",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=9, column=0, sticky='e', padx=5, pady=4)
        self.resistance_outdir = tk.Entry(container, width=30, bg='#f3f4f6', fg='#1f2937',
                                          relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.resistance_outdir.insert(0, "resistance_analysis")
        self.resistance_outdir.grid(row=9, column=1, sticky='w', padx=5, pady=4)
        tk.Button(container, text="▶ Run Analysis", command=self.run_resistance,
                  bg='#2563eb', fg='#ffffff', font=('Segoe UI', 12, 'bold'),
                  relief='flat', padx=30, pady=8).grid(row=10, column=0, columnspan=3, pady=15)
        self.resistance_log = scrolledtext.ScrolledText(container, height=8, bg='#1e293b', fg='#e2e8f0',
                                                        font=('Consolas', 9), relief='flat')
        self.resistance_log.grid(row=11, column=0, columnspan=3, padx=10, pady=10, sticky='we')

    def create_eggnog_tab(self):
        tab = tk.Frame(self.notebook, bg='#ffffff')
        self.notebook.add(tab, text="eggnog2kegg")
        container = tk.Frame(tab, bg='#ffffff', padx=20, pady=20)
        container.pack(fill=tk.BOTH, expand=True)
        tk.Label(container, text="eggNOG to KEGG Converter",
                 font=('Segoe UI', 14, 'bold'), bg='#ffffff', fg='#1f2937').grid(row=0, column=0, columnspan=3, pady=(0, 20), sticky='w')
        tk.Label(container, text="eggNOG annotations file:",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=1, column=0, sticky='e', padx=5, pady=4)
        self.eggnog2kegg_input = tk.Entry(container, width=50, bg='#f3f4f6', fg='#1f2937',
                                          relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.eggnog2kegg_input.grid(row=1, column=1, padx=5, pady=4)
        tk.Button(container, text="📂 Browse", command=lambda: self.browse_file(self.eggnog2kegg_input),
                  bg='#2563eb', fg='#ffffff', relief='flat', padx=10, pady=4,
                  font=('Segoe UI', 9)).grid(row=1, column=2, padx=5)
        tk.Label(container, text="Output file name:",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=2, column=0, sticky='e', padx=5, pady=4)
        self.eggnog2kegg_output = tk.Entry(container, width=30, bg='#f3f4f6', fg='#1f2937',
                                           relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.eggnog2kegg_output.insert(0, "kegg_mapper.txt")
        self.eggnog2kegg_output.grid(row=2, column=1, sticky='w', padx=5, pady=4)
        tk.Label(container, text="Output folder:",
                 bg='#ffffff', fg='#1f2937', font=('Segoe UI', 10)).grid(row=3, column=0, sticky='e', padx=5, pady=4)
        self.eggnog2kegg_outdir = tk.Entry(container, width=30, bg='#f3f4f6', fg='#1f2937',
                                           relief='flat', highlightthickness=1, highlightcolor='#2563eb')
        self.eggnog2kegg_outdir.insert(0, "eggnog2kegg_output")
        self.eggnog2kegg_outdir.grid(row=3, column=1, sticky='w', padx=5, pady=4)
        tk.Button(container, text="▶ Convert", command=self.run_eggnog2kegg,
                  bg='#2563eb', fg='#ffffff', font=('Segoe UI', 12, 'bold'),
                  relief='flat', padx=30, pady=8).grid(row=4, column=0, columnspan=3, pady=15)
        self.eggnog2kegg_log = scrolledtext.ScrolledText(container, height=8, bg='#1e293b', fg='#e2e8f0',
                                                         font=('Consolas', 9), relief='flat')
        self.eggnog2kegg_log.grid(row=5, column=0, columnspan=3, padx=10, pady=10, sticky='we')

    def browse_file(self, entry_widget):
        filename = filedialog.askopenfilename(title="Select file", filetypes=[("All files", "*.*")])
        if filename:
            entry_widget.delete(0, tk.END)
            entry_widget.insert(0, filename)

    def browse_files(self, entry_widget):
        files = filedialog.askopenfilenames(title="Select files", filetypes=[("All files", "*.*")])
        if files:
            entry_widget.delete(0, tk.END)
            entry_widget.insert(0, ",".join(files))

    def browse_dir(self, entry_widget):
        directory = filedialog.askdirectory(title="Select directory")
        if directory:
            entry_widget.delete(0, tk.END)
            entry_widget.insert(0, directory)

    def log(self, widget, msg):
        widget.insert(tk.END, msg + "\n")
        widget.see(tk.END)

    def run_secondary(self):
        eggnog = self.secondary_eggnog.get().strip()
        interpro = self.secondary_interpro.get().strip()
        gff = self.secondary_gff.get().strip()
        genome = self.secondary_genome.get().strip()
        hmm_dir = self.secondary_hmm.get().strip()
        outdir = self.secondary_outdir.get().strip()
        
        if not eggnog and not interpro and not gff:
            messagebox.showerror("Error", "Please provide at least eggNOG, InterPro, or GFF file.")
            return
        
        if (gff and not genome) or (genome and not gff):
            messagebox.showerror("Error", "Both GFF and Genome files are required together.")
            return
        
        try:
            self.log(self.secondary_log, "Starting secondary metabolism analysis...")
            run_secondary_analysis(eggnog, interpro, gff_file=gff if gff else None, 
                                  genome_file=genome if genome else None, hmm_dir=hmm_dir if hmm_dir else None,
                                  output_dir=outdir, prefix="secondary")
            self.log(self.secondary_log, f"✅ Analysis completed. Outputs in: {outdir}")
            messagebox.showinfo("Success", f"Secondary metabolism analysis completed.\nOutputs in: {outdir}")
        except Exception as e:
            self.log(self.secondary_log, f"❌ Error: {e}")
            import traceback
            self.log(self.secondary_log, traceback.format_exc())
            messagebox.showerror("Error", str(e))

    def run_pharma(self):
        eggnog = self.pharma_eggnog.get().strip()
        interpro = self.pharma_interpro.get().strip()
        gff = self.pharma_gff.get().strip()
        genome = self.pharma_genome.get().strip()
        hmm_dir = self.pharma_hmm.get().strip()
        outdir = self.pharma_outdir.get().strip()
        threshold = int(self.pharma_threshold.get().strip() or 4)
        nodes = int(self.pharma_nodes.get().strip() or 15)
        method = self.pharma_method.get().strip() or 'enrichment'
        
        if not eggnog and not interpro and not gff:
            messagebox.showerror("Error", "Please provide at least eggNOG, InterPro, or GFF file.")
            return
        
        if (gff and not genome) or (genome and not gff):
            messagebox.showerror("Error", "Both GFF and Genome files are required together.")
            return
        
        try:
            self.log(self.pharma_log, "Starting pharmaceutical analysis...")
            run_pharma_analysis(eggnog, interpro, gff_file=gff if gff else None,
                               genome_file=genome if genome else None, hmm_dir=hmm_dir if hmm_dir else None,
                               output_dir=outdir, prefix="pharma",
                               threshold=threshold, network_nodes=nodes, network_method=method)
            self.log(self.pharma_log, f"✅ Analysis completed. Outputs in: {outdir}")
            messagebox.showinfo("Success", f"Pharma analysis completed.\nOutputs in: {outdir}")
        except Exception as e:
            self.log(self.pharma_log, f"❌ Error: {e}")
            messagebox.showerror("Error", str(e))

    def run_resistance(self):
        eggnog = self.resistance_eggnog.get().strip()
        interpro = self.resistance_interpro.get().strip()
        gff = self.resistance_gff.get().strip()
        genome = self.resistance_genome.get().strip()
        hmm_dir = self.resistance_hmm.get().strip()
        outdir = self.resistance_outdir.get().strip()
        threshold = int(self.resistance_threshold.get().strip() or 5)
        nodes = int(self.resistance_nodes.get().strip() or 15)
        method = self.resistance_method.get().strip() or 'enrichment'
        
        if not eggnog and not interpro and not gff:
            messagebox.showerror("Error", "Please provide at least eggNOG, InterPro, or GFF file.")
            return
        
        if (gff and not genome) or (genome and not gff):
            messagebox.showerror("Error", "Both GFF and Genome files are required together.")
            return
        
        try:
            self.log(self.resistance_log, "Starting disease resistance analysis...")
            run_disease_resistance_analysis(eggnog, interpro, gff_file=gff if gff else None,
                                           genome_file=genome if genome else None, hmm_dir=hmm_dir if hmm_dir else None,
                                           output_dir=outdir, prefix="resistance",
                                           threshold=threshold, network_nodes=nodes, network_method=method)
            self.log(self.resistance_log, f"✅ Analysis completed. Outputs in: {outdir}")
            messagebox.showinfo("Success", f"Disease resistance analysis completed.\nOutputs in: {outdir}")
        except Exception as e:
            self.log(self.resistance_log, f"❌ Error: {e}")
            messagebox.showerror("Error", str(e))

    def run_eggnog2kegg(self):
        infile = self.eggnog2kegg_input.get().strip()
        outfile = self.eggnog2kegg_output.get().strip()
        outdir = self.eggnog2kegg_outdir.get().strip()
        if not infile or not os.path.exists(infile):
            messagebox.showerror("Error", "Invalid input file.")
            return
        if not outfile:
            outfile = "kegg_mapper.txt"
        try:
            self.log(self.eggnog2kegg_log, "Converting eggNOG to KEGG format...")
            run_eggnog2kegg(infile, outfile, output_dir=outdir, prefix="kegg")
            self.log(self.eggnog2kegg_log, f"✅ Conversion completed. Output: {os.path.join(outdir, outfile)}")
            messagebox.showinfo("Success", f"Conversion completed.\nOutput: {os.path.join(outdir, outfile)}")
        except Exception as e:
            self.log(self.eggnog2kegg_log, f"❌ Error: {e}")
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    root = tk.Tk()
    app = MetabAnalysisApp(root)
    root.mainloop()
