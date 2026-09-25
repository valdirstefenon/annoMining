#!/usr/bin/env python3
import sys
import os
import re
import json
import gzip
import subprocess
import threading
import traceback
from datetime import datetime
from collections import defaultdict, Counter
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.multitest import multipletests

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch, Circle

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


SECONDARY_PATHWAYS = {
    'map00900': {'name': 'Terpenoid backbone', 'class': 'Terpenoids', 'color': '#2E86AB'},
    'map00902': {'name': 'Monoterpenoid', 'class': 'Terpenoids', 'color': '#2E86AB'},
    'map00904': {'name': 'Diterpenoid', 'class': 'Terpenoids', 'color': '#2E86AB'},
    'map00906': {'name': 'Carotenoid', 'class': 'Terpenoids', 'color': '#2E86AB'},
    'map00909': {'name': 'Sesquiterpenoid/Triterpenoid', 'class': 'Terpenoids', 'color': '#2E86AB'},
    'map00940': {'name': 'Phenylpropanoid', 'class': 'Phenylpropanoids', 'color': '#A23B72'},
    'map00941': {'name': 'Flavonoid', 'class': 'Phenylpropanoids', 'color': '#A23B72'},
    'map00942': {'name': 'Anthocyanin', 'class': 'Phenylpropanoids', 'color': '#A23B72'},
    'map00943': {'name': 'Isoflavonoid', 'class': 'Phenylpropanoids', 'color': '#A23B72'},
    'map00944': {'name': 'Flavone/Flavonol', 'class': 'Phenylpropanoids', 'color': '#A23B72'},
    'map00945': {'name': 'Stilbenoid/Curcuminoid', 'class': 'Phenylpropanoids', 'color': '#A23B72'},
    'map00950': {'name': 'Isoquinoline alkaloid', 'class': 'Alkaloids', 'color': '#F18F01'},
    'map00960': {'name': 'Tropane/Piperidine', 'class': 'Alkaloids', 'color': '#F18F01'},
    'map00965': {'name': 'Betalain', 'class': 'Alkaloids', 'color': '#F18F01'},
    'map00966': {'name': 'Glucosinolate', 'class': 'Glucosinolates', 'color': '#9B5DE5'},
    'map00903': {'name': 'Cannabinoid biosynthesis', 'class': 'Cannabinoids', 'color': '#FF006E'},
    'map00232': {'name': 'Caffeine metabolism', 'class': 'Xanthines', 'color': '#F15BB5'},
}

PHARMA_DETAILS = {
    'map00950': {'name': 'Isoquinoline alkaloid', 'class': 'Alkaloids', 'pharma_use': 'Analgesic, antimicrobial'},
    'map00960': {'name': 'Tropane/Piperidine', 'class': 'Alkaloids', 'pharma_use': 'Anticholinergic'},
    'map00965': {'name': 'Betalain', 'class': 'Alkaloids', 'pharma_use': 'Antioxidant'},
    'map00900': {'name': 'Terpenoid backbone', 'class': 'Terpenoids', 'pharma_use': 'Diverse bioactivities'},
    'map00909': {'name': 'Sesquiterpenoid/Triterpenoid', 'class': 'Terpenoids', 'pharma_use': 'Antimalarial, anticancer'},
    'map00904': {'name': 'Diterpenoid', 'class': 'Terpenoids', 'pharma_use': 'Cardiovascular'},
    'map00906': {'name': 'Carotenoid', 'class': 'Terpenoids', 'pharma_use': 'Antioxidant'},
    'map00940': {'name': 'Phenylpropanoid', 'class': 'Phenylpropanoids', 'pharma_use': 'Antioxidant'},
    'map00941': {'name': 'Flavonoid', 'class': 'Flavonoids', 'pharma_use': 'Cardioprotective'},
    'map00942': {'name': 'Anthocyanin', 'class': 'Flavonoids', 'pharma_use': 'Antioxidant'},
    'map00944': {'name': 'Flavone/Flavonol', 'class': 'Flavonoids', 'pharma_use': 'Anti-inflammatory'},
    'map00945': {'name': 'Stilbenoid/Curcuminoid', 'class': 'Stilbenoids', 'pharma_use': 'Chemopreventive'},
    'map00966': {'name': 'Glucosinolate', 'class': 'Glucosinolates', 'pharma_use': 'Anticancer'},
    'map00232': {'name': 'Caffeine metabolism', 'class': 'Xanthines', 'pharma_use': 'Stimulant'},
    'map00903': {'name': 'Cannabinoid biosynthesis', 'class': 'Cannabinoids', 'pharma_use': 'Analgesic, anti-epileptic'},
}

PHARMA_KOS = {
    'Bioactive_flavonoids': ['K00660', 'K05275', 'K05276', 'K13065', 'K13066', 'K05265', 'K05266', 'K13083', 'K13084'],
    'Anticancer_terpenes': ['K15891', 'K15892', 'K15893', 'K00487', 'K00507', 'K00509', 'K00511', 'K00512'],
    'Nutraceutical_carotenoids': ['K00514', 'K00515', 'K06444', 'K06445', 'K06446', 'K06447', 'K06448', 'K06449', 'K06450', 'K06451'],
    'Resveratrol_biosynthesis': ['K13071', 'K13072', 'K00430', 'K00431'],
    'Therapeutic_alkaloids': ['K01799', 'K01800', 'K01900', 'K01901', 'K01902', 'K01903'],
}

CORE_PHARMA_DOMAINS = {
    'Flavonoids': ['PF00195', 'PF02797', 'PF05834'],
    'Terpenoids': ['PF01397', 'PF03936', 'PF00067', 'PF00494'],
    'Alkaloids': ['PF01596', 'PF00891', 'PF00201', 'PF00155', 'PF00141'],
    'Phenylpropanoids': ['PF00195', 'PF00141'],
    'Carotenoids': ['PF00494', 'PF00514'],
    'P450': ['PF00067'],
    'Methyltransferases': ['PF01596', 'PF00891', 'PF01799'],
}

DISEASE_RESISTANCE_PFAMS = {
    'NLR': {'pfams': ['PF00931', 'PF01582', 'PF00560', 'PF07723'], 'color': '#E63946'},
    'TIR_NLR': {'pfams': ['PF00931', 'PF01582', 'PF13676', 'PF00560'], 'color': '#F4A261'},
    'TNJ': {'pfams': ['PF01582', 'PF00931', 'PF01419'], 'color': '#2A9D8F'},
    'TNJ_like': {'pfams': ['PF01582', 'PF00931', 'PF01419', 'PF00560'], 'color': '#264653'},
    'CC_NLR': {'pfams': ['PF00931', 'PF00560'], 'color': '#E9C46A'},
    'RPW8': {'pfams': ['PF05659', 'PF05660'], 'color': '#9B5DE5'},
    'Jacalin': {'pfams': ['PF01419'], 'color': '#F15BB5'},
    'LRR_repeat': {'pfams': ['PF00560', 'PF07723', 'PF13855'], 'color': '#FF006E'},
    'NB_ARC': {'pfams': ['PF00931'], 'color': '#00BBF9'},
    'TIR': {'pfams': ['PF01582', 'PF13676'], 'color': '#00F5D4'},
    'CC': {'pfams': ['PF00560'], 'color': '#FEE440'},
    'NLR_Jacalin': {'pfams': ['PF00931', 'PF01419'], 'color': '#8338EC'},
    'TNL': {'pfams': ['PF01582', 'PF00931', 'PF00560'], 'color': '#FB8500'},
}

DISEASE_RESISTANCE_PATHWAYS = {
    'plant_pathogen_interaction': {'name': 'Plant-pathogen interaction', 'kegg': 'map04626', 'color': '#E63946'},
    'MAPK_signaling_plant': {'name': 'MAPK signaling - plant', 'kegg': 'map04016', 'color': '#F4A261'},
    'hormone_signal_transduction': {'name': 'Plant hormone signal transduction', 'kegg': 'map04075', 'color': '#2A9D8F'},
    'phenylpropanoid_biosynthesis': {'name': 'Phenylpropanoid biosynthesis', 'kegg': 'map00940', 'color': '#E9C46A'},
    'flavonoid_biosynthesis': {'name': 'Flavonoid biosynthesis', 'kegg': 'map00941', 'color': '#9B5DE5'},
    'terpenoid_backbone': {'name': 'Terpenoid backbone biosynthesis', 'kegg': 'map00900', 'color': '#F15BB5'},
}

DISEASE_RESISTANCE_KOS = {
    'NLR_genes': ['K13475', 'K13476', 'K13477', 'K13478', 'K13479'],
    'MAPK_signaling': ['K02540', 'K04345', 'K04436', 'K04437', 'K04438', 'K04439', 'K04440'],
    'Plant_pathogen_interaction': ['K02607', 'K03081', 'K03120', 'K03244', 'K03245', 'K03246', 'K03247'],
    'Hormone_signaling': ['K06445', 'K06446', 'K06447', 'K06448', 'K06449', 'K06450'],
    'Defense_proteins': ['K13480', 'K13481', 'K13482', 'K13483', 'K13484'],
}

SECONDARY_HMM = {
    'Terpene_synthase': {'pfam': 'PF01397', 'hmm': 'Terpene_synthase.hmm', 'class': 'Terpenoids'},
    'Terpene_synthase_N': {'pfam': 'PF03936', 'hmm': 'Terpene_synthase_N.hmm', 'class': 'Terpenoids'},
    'Cytochrome_P450': {'pfam': 'PF00067', 'hmm': 'Cytochrome_P450.hmm', 'class': 'Terpenoids'},
    'Chalcone_synthase': {'pfam': 'PF00195', 'hmm': 'Chalcone_synthase.hmm', 'class': 'Phenylpropanoids'},
    'O_methyltransferase': {'pfam': 'PF00891', 'hmm': 'O_methyltransferase.hmm', 'class': 'Alkaloids'},
    'SAM_methyltransferase': {'pfam': 'PF01596', 'hmm': 'SAM_methyltransferase.hmm', 'class': 'Alkaloids'},
    'Prenyltransferase': {'pfam': 'PF00494', 'hmm': 'Prenyltransferase.hmm', 'class': 'Terpenoids'},
    'Flavonoid_3_hydroxylase': {'pfam': 'PF02797', 'hmm': 'Flavonoid_3_hydroxylase.hmm', 'class': 'Flavonoids'},
}

RESISTANCE_HMM = {
    'NB_ARC': {'pfam': 'PF00931', 'hmm': 'NB_ARC.hmm', 'class': 'NB_ARC'},
    'TIR': {'pfam': 'PF01582', 'hmm': 'TIR.hmm', 'class': 'TIR'},
    'LRR': {'pfam': 'PF00560', 'hmm': 'LRR.hmm', 'class': 'LRR_repeat'},
    'Jacalin': {'pfam': 'PF01419', 'hmm': 'Jacalin.hmm', 'class': 'Jacalin'},
    'RPW8': {'pfam': 'PF05659', 'hmm': 'RPW8.hmm', 'class': 'RPW8'},
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
    'PF00941', 'PF01565', 'PF02727', 'PF03055', 'PF08491',
]


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def timestamp():
    return datetime.now().strftime('%Y%m%d_%H%M%S')


def find_hmm_file(hmm_file, hmm_dir):
    if os.path.exists(hmm_file):
        return hmm_file
    if hmm_dir and os.path.exists(os.path.join(hmm_dir, hmm_file)):
        return os.path.join(hmm_dir, hmm_file)
    for d in ['/usr/share/hmmer/profiles/', '/opt/hmmer/profiles/',
              os.path.expanduser('~/hmmer_profiles/'),
              os.path.expanduser('~/hmmer/hmmer_profiles/')]:
        c = os.path.join(d, hmm_file)
        if os.path.exists(c):
            return c
    return None


def check_executable(name):
    try:
        subprocess.run([name, '--version'], capture_output=True, check=True)
        return True
    except Exception:
        return False


def write_session_info(output_dir):
    import platform
    versions = {'python': platform.python_version()}
    for pkg in ('pandas', 'numpy', 'matplotlib', 'seaborn', 'scipy', 'statsmodels'):
        try:
            mod = __import__(pkg)
            versions[pkg] = getattr(mod, '__version__', 'unknown')
        except Exception:
            versions[pkg] = 'not installed'
    with open(os.path.join(output_dir, 'session_info.txt'), 'w') as f:
        f.write(f"Analysis date: {datetime.now().isoformat()}\n")
        for k, v in versions.items():
            f.write(f"{k}: {v}\n")


def write_parameters(output_dir, params):
    with open(os.path.join(output_dir, 'parameters.json'), 'w') as f:
        json.dump(params, f, indent=2, default=str)


def parse_eggnog(eggnog_file):
    data = {
        'genes': {},
        'kegg_pathways': defaultdict(set),
        'ec_numbers': defaultdict(set),
        'ko_numbers': defaultdict(set),
        'go_terms': defaultdict(set),
        'pfam_domains': defaultdict(set),
        'cog_categories': Counter(),
    }
    if not eggnog_file or not os.path.exists(eggnog_file):
        return data
    ko_re = re.compile(r'(?:ko:)?(K\d{5})', re.IGNORECASE)
    path_re = re.compile(r'(?:ko)?(map\d{5})', re.IGNORECASE)
    ec_re = re.compile(r'(\d+\.\d+\.\d+\.\d+)')
    go_re = re.compile(r'(GO:\d+)', re.IGNORECASE)
    pfam_re = re.compile(r'(PF\d{5})', re.IGNORECASE)
    with open(eggnog_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 12:
                continue
            gene = fields[0]
            description = fields[7] if len(fields) > 7 and fields[7] != '-' else '-'
            preferred = fields[8] if len(fields) > 8 and fields[8] != '-' else ''
            cog = fields[6] if len(fields) > 6 and fields[6] != '-' else '-'
            data['genes'][gene] = {'description': description, 'preferred_name': preferred, 'cog_category': cog}
            for c in cog:
                if c.isalpha():
                    data['cog_categories'][c] += 1
            for ko in ko_re.findall(line):
                data['ko_numbers'][gene].add(ko)
            for path in path_re.findall(line):
                data['kegg_pathways'][gene].add(path)
            for ec in ec_re.findall(line):
                data['ec_numbers'][gene].add(ec)
            for go in go_re.findall(line):
                data['go_terms'][gene].add(go)
            for pfam in pfam_re.findall(line):
                data['pfam_domains'][gene].add(pfam)
    return data


def parse_interpro(interpro_file):
    data = {'genes': set(), 'pfam_domains': defaultdict(set),
            'go_terms': defaultdict(set), 'interpro_domains': defaultdict(set)}
    if not interpro_file or not os.path.exists(interpro_file):
        return data
    pfam_re = re.compile(r'(PF\d{5})', re.IGNORECASE)
    go_re = re.compile(r'(GO:\d+)', re.IGNORECASE)
    with open(interpro_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 5:
                continue
            gene = fields[0]
            data['genes'].add(gene)
            for field in fields:
                for m in pfam_re.findall(field):
                    data['pfam_domains'][gene].add(m)
            for go in go_re.findall(line):
                data['go_terms'][gene].add(go)
            if len(fields) > 12 and fields[12] and fields[12] != '-':
                data['interpro_domains'][gene].add(fields[12])
    return data


def parse_eggnog_kegg_mapper(eggnog_file, output_file):
    results = []
    total = 0
    with_ko = 0
    ko_re = re.compile(r'(K\d{5})')
    with open(eggnog_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 12:
                continue
            query = fields[0]
            kegg_ko = fields[11] if len(fields) > 11 else ''
            total += 1
            kos = ko_re.findall(kegg_ko) if kegg_ko and kegg_ko != '-' else []
            if kos:
                for ko in kos:
                    results.append(f"{query}\t{ko}")
                    with_ko += 1
            else:
                results.append(query)
    with open(output_file, 'w') as f:
        f.write('\n'.join(results) + '\n')
    return total, with_ko


def extract_proteins_from_gff(gff_file, genome_file, output_fasta):
    if not gff_file or not os.path.exists(gff_file):
        return False
    if not genome_file or not os.path.exists(genome_file):
        return False
    if not check_executable('gffread'):
        return False
    cmd = ['gffread', gff_file, '-g', genome_file, '-y', output_fasta]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True)
        return r.returncode == 0 and os.path.exists(output_fasta) and os.path.getsize(output_fasta) > 0
    except Exception:
        return False


def run_hmmer(protein_fasta, hmm_dir, output_dir, domains_dict, evalue=1e-10, prefix=''):
    results = defaultdict(lambda: {'domains': [], 'classes': set(), 'pfams': [], 'max_evalue': float('inf')})
    if not check_executable('hmmsearch'):
        return results
    if not protein_fasta or not os.path.exists(protein_fasta):
        return results
    for domain_name, info in domains_dict.items():
        hmm_file = info.get('hmm')
        hmm_path = find_hmm_file(hmm_file, hmm_dir)
        if not hmm_path:
            continue
        out_file = os.path.join(output_dir, f"{prefix}_{domain_name}_hits.tsv")
        cmd = ['hmmsearch', '--domtblout', out_file, '--noali', '-E', str(evalue), hmm_path, protein_fasta]
        try:
            subprocess.run(cmd, capture_output=True, check=True)
        except Exception:
            continue
        if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
            continue
        with open(out_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.split()
                if len(fields) < 12:
                    continue
                gene = fields[0]
                domain = fields[3]
                try:
                    ev = float(fields[11])
                except ValueError:
                    continue
                entry = results[gene]
                entry['domains'].append({'domain_name': domain_name, 'domain': domain, 'evalue': ev})
                if 'class' in info:
                    entry['classes'].add(info['class'])
                if 'pfam' in info:
                    entry['pfams'].append(info['pfam'])
                entry['max_evalue'] = min(entry['max_evalue'], ev)
    return results


def merge_annotations(eggnog_data, interpro_data, pathway_dict):
    merged = {}
    all_genes = set(eggnog_data['genes'].keys()) | interpro_data['genes']
    for gene in all_genes:
        merged[gene] = {
            'source': [],
            'description': eggnog_data['genes'].get(gene, {}).get('description', '-'),
            'preferred_name': eggnog_data['genes'].get(gene, {}).get('preferred_name', ''),
            'cog_category': eggnog_data['genes'].get(gene, {}).get('cog_category', '-'),
            'go_terms': eggnog_data['go_terms'].get(gene, set()) | interpro_data['go_terms'].get(gene, set()),
            'pfam_domains': eggnog_data['pfam_domains'].get(gene, set()) | interpro_data['pfam_domains'].get(gene, set()),
            'kegg_pathways': eggnog_data['kegg_pathways'].get(gene, set()),
            'ec_numbers': eggnog_data['ec_numbers'].get(gene, set()),
            'ko_numbers': eggnog_data['ko_numbers'].get(gene, set()),
            'secondary_pathways': set(),
            'secondary_class': set(),
            'hmmer_domains': [],
            'hmmer_classes': set(),
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


def integrate_hmmer(merged_data, hmmer_results, pathway_dict):
    added = updated = 0
    for gene, info in hmmer_results.items():
        if gene in merged_data:
            target = merged_data[gene]
            target['hmmer_domains'].extend(info['domains'])
            target['hmmer_classes'].update(info['classes'])
            for pfam in info['pfams']:
                target['pfam_domains'].add(pfam)
            for class_name in info['classes']:
                for pathway, p_info in pathway_dict.items():
                    if p_info['class'] == class_name:
                        target['secondary_pathways'].add(pathway)
                        target['secondary_class'].add(class_name)
                        break
            if 'HMMER' not in target['source']:
                target['source'].append('HMMER')
            updated += 1
        else:
            merged_data[gene] = {
                'source': ['HMMER'],
                'description': f"HMMER: {', '.join(info['classes'])}",
                'preferred_name': '',
                'cog_category': '-',
                'go_terms': set(),
                'pfam_domains': set(info['pfams']),
                'kegg_pathways': set(),
                'ec_numbers': set(),
                'ko_numbers': set(),
                'secondary_pathways': set(),
                'secondary_class': set(),
                'hmmer_domains': info['domains'],
                'hmmer_classes': info['classes'],
            }
            for class_name in info['classes']:
                for pathway, p_info in pathway_dict.items():
                    if p_info['class'] == class_name:
                        merged_data[gene]['secondary_pathways'].add(pathway)
                        merged_data[gene]['secondary_class'].add(class_name)
                        break
            added += 1
    return added, updated


def calculate_pharma_score(merged_data, threshold=4):
    for gene, info in merged_data.items():
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
            info['pharma_score'] = 0
            info['pharma_potential'] = 'NONE'
            info['pharma_pathways'] = []
            info['compound_classes'] = []
            info['pharma_uses'] = []
            continue
        for ko in info['ko_numbers']:
            for ko_list in PHARMA_KOS.values():
                if ko in ko_list:
                    score += 3
                    break
        for pfam in info['pfam_domains']:
            if pfam in GENERIC_DOMAINS:
                continue
            for domain_list in CORE_PHARMA_DOMAINS.values():
                if pfam in domain_list:
                    score += 2
                    break
        if score >= 10:
            potential = 'HIGH'
        elif score >= threshold:
            potential = 'MEDIUM'
        else:
            potential = 'LOW'
        info['pharma_score'] = score
        info['pharma_potential'] = potential
        info['pharma_pathways'] = list(set(pathways_found))
        info['compound_classes'] = list(compound_classes)
        info['pharma_uses'] = list(pharma_uses)


def calculate_resistance_score(merged_data, threshold=5):
    for gene, info in merged_data.items():
        score = 0
        classes = set()
        pathways = set()
        pfams = set()
        gene_pfams = set(info['pfam_domains'])
        has_tir = 'PF01582' in gene_pfams or 'PF13676' in gene_pfams
        has_nb = 'PF00931' in gene_pfams
        has_jac = 'PF01419' in gene_pfams
        has_lrr = any(p in gene_pfams for p in ['PF00560', 'PF07723', 'PF13855'])
        tnj_detected = False
        tnj_arch = ''
        tnj_conf = 'NONE'
        if has_tir and has_nb and has_jac and not has_lrr:
            classes.add('TNJ')
            score += 10
            pfams.update(['PF01582', 'PF00931', 'PF01419'])
            tnj_detected = True
            tnj_arch = 'TIR-NB-ARC-Jacalin (canonical TNJ)'
            tnj_conf = 'HIGH'
        elif has_tir and has_nb and has_jac and has_lrr:
            classes.add('TNJ_like')
            score += 8
            pfams.update(['PF01582', 'PF00931', 'PF01419'])
            tnj_detected = True
            tnj_arch = 'TIR-NB-ARC-Jacalin-LRR (TNJ-like)'
            tnj_conf = 'MEDIUM'
        elif has_nb and has_jac and not has_tir:
            classes.add('NLR_Jacalin')
            score += 6
            pfams.update(['PF00931', 'PF01419'])
            tnj_arch = 'NB-ARC-Jacalin'
        elif has_tir and has_nb and has_lrr and not has_jac:
            classes.add('TNL')
            score += 8
            pfams.update(['PF01582', 'PF00931', 'PF00560'])
            tnj_arch = 'TIR-NB-ARC-LRR'
        elif has_nb and has_lrr and not has_tir:
            classes.add('CC_NLR')
            score += 7
            pfams.update(['PF00560', 'PF00931'])
        for pfam in gene_pfams:
            for class_name, class_info in DISEASE_RESISTANCE_PFAMS.items():
                if pfam in class_info['pfams']:
                    if class_name not in ['TNJ', 'TNJ_like', 'NLR_Jacalin', 'TNL', 'CC_NLR']:
                        classes.add(class_name)
                        pfams.add(pfam)
                        if class_name in ['NLR', 'TIR_NLR']:
                            score += 3
                        elif class_name in ['LRR_repeat', 'NB_ARC', 'TIR', 'CC']:
                            score += 2
                        else:
                            score += 1
                        break
        for ko in info['ko_numbers']:
            for ko_list in DISEASE_RESISTANCE_KOS.values():
                if ko in ko_list:
                    score += 2
                    break
        for pathway in info['kegg_pathways']:
            for rp, rp_info in DISEASE_RESISTANCE_PATHWAYS.items():
                if pathway == rp_info.get('kegg', ''):
                    score += 2
                    pathways.add(rp)
                    break
        if len(classes) >= 3:
            score += 2
        if score >= 10:
            potential = 'HIGH'
        elif score >= threshold:
            potential = 'MEDIUM'
        elif score > 0:
            potential = 'LOW'
        else:
            potential = 'NONE'
        info['resistance_score'] = score
        info['resistance_potential'] = potential
        info['resistance_classes'] = list(classes)
        info['resistance_pathways'] = list(pathways)
        info['resistance_pfams'] = list(pfams)
        info['tnj_detected'] = tnj_detected
        info['tnj_architecture'] = tnj_arch
        info['tnj_confidence'] = tnj_conf


def calculate_enrichment(group_genes, background_genes, merged_data,
                         min_count=3, min_fold=1.0, filter_generic=True):
    if not group_genes or not background_genes:
        return pd.DataFrame()
    group_counts = Counter()
    for gene in group_genes:
        if gene not in merged_data:
            continue
        for pfam in merged_data[gene]['pfam_domains']:
            if filter_generic and pfam in GENERIC_DOMAINS:
                continue
            group_counts[pfam] += 1
    bg_counts = Counter()
    for gene in background_genes:
        if gene not in merged_data:
            continue
        for pfam in merged_data[gene]['pfam_domains']:
            if filter_generic and pfam in GENERIC_DOMAINS:
                continue
            bg_counts[pfam] += 1
    total_g = len(group_genes)
    total_b = len(background_genes)
    rows = []
    for pfam in set(group_counts) | set(bg_counts):
        cg = group_counts.get(pfam, 0)
        cb = bg_counts.get(pfam, 0)
        if cg < min_count or cb < min_count:
            continue
        fg = cg / total_g
        fb = cb / total_b
        fold = fg / fb if fb > 0 else np.inf
        if fold < min_fold:
            continue
        table = [[cg, total_g - cg], [cb, total_b - cb]]
        try:
            _, p_f = fisher_exact(table, alternative='two-sided')
        except Exception:
            p_f = 1.0
        try:
            _, p_c, _, _ = chi2_contingency(table)
        except Exception:
            p_c = 1.0
        rows.append({
            'pfam_domain': pfam,
            'count_in_group': cg, 'total_group': total_g,
            'freq_group': round(fg, 4),
            'count_in_background': cb, 'total_background': total_b,
            'freq_background': round(fb, 4),
            'fold_enrichment': round(fold, 2),
            'p_value': min(p_f, p_c),
            'p_adjust': np.nan, 'significant': False,
        })
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    _, p_adj, _, _ = multipletests(df['p_value'].values, method='fdr_bh')
    df['p_adjust'] = p_adj
    df['significant'] = df['p_adjust'] < 0.05
    return df


def build_network(merged_data, score_key, class_key, pathway_key, pfam_key,
                  min_edge_weight=1, max_nodes=15, method='enrichment'):
    candidates = {g: i for g, i in merged_data.items() if i.get(score_key, 0) > 0}
    if len(candidates) < 2:
        return pd.DataFrame(), []
    if method == 'score':
        selected = sorted(candidates.keys(), key=lambda x: candidates[x][score_key], reverse=True)[:max_nodes]
    elif method == 'enrichment':
        importance = {}
        for g, i in candidates.items():
            importance[g] = (i[score_key] * 2) + (len(i.get(class_key, [])) * 3) + \
                            len(i.get(pathway_key, [])) + (len(i.get(pfam_key, [])) * 0.5)
        selected = sorted(importance.keys(), key=lambda x: importance[x], reverse=True)[:max_nodes]
    else:
        class_genes = {}
        for g, i in candidates.items():
            for rc in i.get(class_key, []):
                class_genes.setdefault(rc, []).append((g, i[score_key]))
        selected = []
        for genes in class_genes.values():
            selected.extend([x[0] for x in sorted(genes, key=lambda x: x[1], reverse=True)[:3]])
        if len(selected) < max_nodes:
            remaining = [g for g in candidates if g not in selected]
            selected.extend(sorted(remaining, key=lambda x: candidates[x][score_key], reverse=True)[:max_nodes - len(selected)])
        selected = selected[:max_nodes]
    edges = []
    for i, g1 in enumerate(selected):
        for g2 in selected[i + 1:]:
            sc = len(set(merged_data[g1].get(class_key, [])) & set(merged_data[g2].get(class_key, [])))
            sp = len(set(merged_data[g1].get(pathway_key, [])) & set(merged_data[g2].get(pathway_key, [])))
            spf = len(set(merged_data[g1].get(pfam_key, [])) & set(merged_data[g2].get(pfam_key, [])))
            sk = len(set(merged_data[g1]['ko_numbers']) & set(merged_data[g2]['ko_numbers']))
            total = (sc * 3) + (sp * 2) + spf + (sk * 2)
            if total >= min_edge_weight:
                edges.append({'gene1': g1, 'gene2': g2, 'shared_classes': sc,
                              'shared_pathways': sp, 'shared_pfams': spf,
                              'shared_kos': sk, 'similarity_score': total})
    return (pd.DataFrame(edges) if edges else pd.DataFrame()), selected


def configure_style():
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['font.size'] = 10
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42


def save_figure(fig, prefix, suffix):
    for ext in ('pdf', 'png', 'svg'):
        fig.savefig(f"{prefix}_{suffix}.{ext}", bbox_inches='tight')
    plt.close(fig)


def plot_ko_dotplot(merged_data, output_prefix, top_n=15):
    ko_by_class = defaultdict(Counter)
    class_counts = Counter()
    for info in merged_data.values():
        for pathway in info['secondary_pathways']:
            if pathway in SECONDARY_PATHWAYS:
                pc = SECONDARY_PATHWAYS[pathway]['class']
                class_counts[pc] += 1
                for ko in info['ko_numbers']:
                    ko_by_class[pc][ko] += 1
    rows = []
    for cn, counter in ko_by_class.items():
        total = class_counts[cn]
        for ko, count in counter.most_common(top_n):
            rows.append({'class': cn, 'ko': ko, 'count': count,
                         'frequency': count / total if total else 0})
    if not rows:
        return
    df = pd.DataFrame(rows).sort_values('frequency', ascending=False).head(top_n)
    fig, ax = plt.subplots(figsize=(14, max(6, len(df) * 0.4)))
    colors_map = {'Terpenoids': '#2E86AB', 'Phenylpropanoids': '#A23B72',
                  'Alkaloids': '#F18F01', 'Glucosinolates': '#9B5DE5',
                  'Xanthines': '#F15BB5', 'Cannabinoids': '#FF006E'}
    colors = [colors_map.get(c, '#999999') for c in df['class']]
    ax.scatter(df['frequency'], df['ko'], s=df['count'] * 20 + 30, c=colors,
               alpha=0.8, edgecolors='black', linewidth=0.5)
    ax.set_xlabel('Frequency in class')
    ax.set_ylabel('KEGG Ortholog')
    ax.set_title('Top KOs by Secondary Metabolite Class')
    handles = [Patch(facecolor=colors_map[c], label=c) for c in df['class'].unique() if c in colors_map]
    if handles:
        ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.5))
        plt.subplots_adjust(right=0.78)
    plt.tight_layout()
    save_figure(fig, output_prefix, 'ko_dotplot')


def plot_pathway_completeness(merged_data, output_prefix):
    counts = Counter()
    for info in merged_data.values():
        for p in info['secondary_pathways']:
            counts[p] += 1
    if not counts:
        return
    pathways = list(counts.keys())
    fig, ax = plt.subplots(figsize=(14, max(6, len(pathways) * 0.4)))
    names = [SECONDARY_PATHWAYS[p]['name'] for p in pathways]
    colors = [SECONDARY_PATHWAYS[p]['color'] for p in pathways]
    bars = ax.barh(range(len(pathways)), [counts[p] for p in pathways], color=colors, alpha=0.8)
    ax.set_yticks(range(len(pathways)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel('Number of genes')
    ax.set_title('Genes Associated with Secondary Metabolism Pathways')
    for bar, p in zip(bars, pathways):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2, str(counts[p]), va='center')
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    save_figure(fig, output_prefix, 'pathway_completeness')


def plot_class_distribution(merged_data, output_prefix):
    counts = Counter()
    for info in merged_data.values():
        for c in info['secondary_class']:
            counts[c] += 1
    if not counts:
        return
    colors_map = {'Terpenoids': '#2E86AB', 'Phenylpropanoids': '#A23B72',
                  'Alkaloids': '#F18F01', 'Glucosinolates': '#9B5DE5',
                  'Xanthines': '#F15BB5', 'Cannabinoids': '#FF006E'}
    colors = [colors_map.get(c, '#999999') for c in counts.keys()]
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.pie(counts.values(), labels=counts.keys(), colors=colors,
           autopct='%1.1f%%', startangle=90, pctdistance=0.8)
    ax.set_title('Distribution by Metabolite Class')
    plt.tight_layout()
    save_figure(fig, output_prefix, 'class_distribution')


def plot_top_ec(merged_data, output_prefix, top_n=15):
    counts = Counter()
    for info in merged_data.values():
        for ec in info['ec_numbers']:
            counts[ec] += 1
    if not counts:
        return
    top = counts.most_common(top_n)
    names, values = zip(*top)
    fig, ax = plt.subplots(figsize=(14, max(6, len(names) * 0.4)))
    bars = ax.barh(range(len(names)), values, color='#A23B72', alpha=0.8)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel('Number of genes')
    ax.set_title('Most Abundant EC Numbers')
    for bar, v in zip(bars, values):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2, str(v), va='center')
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    save_figure(fig, output_prefix, 'top_ec_numbers')


def plot_top_pfam(merged_data, output_prefix, top_n=15):
    counts = Counter()
    for info in merged_data.values():
        for pfam in info['pfam_domains']:
            if pfam in GENERIC_DOMAINS:
                continue
            counts[pfam] += 1
    if not counts:
        return
    top = counts.most_common(top_n)
    names, values = zip(*top)
    fig, ax = plt.subplots(figsize=(14, max(6, len(names) * 0.4)))
    bars = ax.barh(range(len(names)), values, color='#F18F01', alpha=0.8)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel('Number of genes')
    ax.set_title('Most Abundant Core PFAM Domains')
    for bar, v in zip(bars, values):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2, str(v), va='center')
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    save_figure(fig, output_prefix, 'top_core_pfam_domains')


def plot_score_ranking(merged_data, score_key, potential_key, output_prefix, title):
    positive = {g: i for g, i in merged_data.items() if i.get(score_key, 0) > 0}
    if not positive:
        return
    sorted_genes = sorted(positive.items(), key=lambda x: x[1][score_key], reverse=True)[:20]
    genes = [g[0] for g in sorted_genes]
    scores = [g[1][score_key] for g in sorted_genes]
    potentials = [g[1][potential_key] for g in sorted_genes]
    colors_map = {'HIGH': '#2ECC71', 'MEDIUM': '#F39C12', 'LOW': '#E74C3C'}
    bar_colors = [colors_map.get(p, '#95A5A6') for p in potentials]
    fig, ax = plt.subplots(figsize=(16, max(8, len(genes) * 0.4)))
    bars = ax.barh(range(len(genes)), scores, color=bar_colors, alpha=0.8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=9)
    ax.set_xlabel(title)
    ax.set_title(f'Top 20 Genes by {title}')
    for bar, s, p in zip(bars, scores, potentials):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2, f'{s} ({p})', va='center')
    handles = [Patch(facecolor='#2ECC71', label='HIGH'),
               Patch(facecolor='#F39C12', label='MEDIUM'),
               Patch(facecolor='#E74C3C', label='LOW')]
    ax.legend(handles=handles, loc='lower right')
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    save_figure(fig, output_prefix, 'score_ranking')


def plot_resistance_class_distribution(merged_data, output_prefix):
    counts = Counter()
    for info in merged_data.values():
        if info.get('tnj_detected', False):
            counts['TNJ'] += 1
            continue
        for rc in info.get('resistance_classes', []):
            counts[rc] += 1
    if not counts:
        return
    colors_map = {c: v.get('color', '#999999') for c, v in DISEASE_RESISTANCE_PFAMS.items()}
    colors_map['TNJ'] = '#2A9D8F'
    classes = list(counts.keys())
    values = [counts[c] for c in classes]
    idx = np.argsort(values)
    classes = [classes[i] for i in idx]
    values = [values[i] for i in idx]
    colors = [colors_map.get(c, '#999999') for c in classes]
    fig, ax = plt.subplots(figsize=(14, max(6, len(classes) * 0.5)))
    bars = ax.barh(classes, values, color=colors, alpha=0.8)
    ax.set_xlabel('Number of genes')
    ax.set_title('Distribution of Disease Resistance Gene Classes')
    for bar, v in zip(bars, values):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2, str(v), va='center')
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    save_figure(fig, output_prefix, 'resistance_class_distribution')


def plot_resistance_pathways(merged_data, output_prefix):
    counts = Counter()
    for info in merged_data.values():
        for p in info.get('resistance_pathways', []):
            counts[p] += 1
    if not counts:
        return
    pathways = list(counts.keys())
    fig, ax = plt.subplots(figsize=(14, max(6, len(pathways) * 0.4)))
    names = [DISEASE_RESISTANCE_PATHWAYS.get(p, {}).get('name', p) for p in pathways]
    colors = [DISEASE_RESISTANCE_PATHWAYS.get(p, {}).get('color', '#2E86AB') for p in pathways]
    bars = ax.barh(range(len(pathways)), [counts[p] for p in pathways], color=colors, alpha=0.8)
    ax.set_yticks(range(len(pathways)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel('Number of genes')
    ax.set_title('Genes Associated with Disease Resistance Pathways')
    for bar, p in zip(bars, pathways):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2, str(counts[p]), va='center')
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    save_figure(fig, output_prefix, 'resistance_pathway_completeness')


def plot_donut(merged_data, potential_key, output_prefix, title):
    counts = Counter(i[potential_key] for i in merged_data.values() if i.get(potential_key) not in (None, 'NONE'))
    if not counts:
        return
    colors_map = {'HIGH': '#2ECC71', 'MEDIUM': '#F39C12', 'LOW': '#E74C3C'}
    colors = [colors_map.get(p, '#95A5A6') for p in counts.keys()]
    labels = [f"{p}\n({c} genes)" for p, c in counts.items()]
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.pie(counts.values(), labels=labels, colors=colors, autopct='%1.1f%%',
           startangle=90, pctdistance=0.85)
    ax.add_artist(Circle((0, 0), 0.70, fc='white', edgecolor='black'))
    ax.text(0, 0, f'Total\n{sum(counts.values())}\ngenes', ha='center', va='center',
            fontsize=14, fontweight='bold')
    ax.set_title(title)
    plt.tight_layout()
    save_figure(fig, output_prefix, 'donut')


def plot_network(merged_data, selected, edges_df, score_key, potential_key, output_prefix, title, tnj_key=None):
    if not selected:
        return
    colors_map = {'HIGH': '#2ECC71', 'MEDIUM': '#F39C12', 'LOW': '#E74C3C', 'NONE': '#95A5A6'}
    n = len(selected)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    pos = {g: (np.cos(a), np.sin(a)) for g, a in zip(selected, angles)}
    fig, ax = plt.subplots(figsize=(14, 12))
    if edges_df is not None and len(edges_df) > 0:
        max_sim = edges_df['similarity_score'].max() or 1
        for _, edge in edges_df.iterrows():
            g1, g2 = edge['gene1'], edge['gene2']
            if g1 in pos and g2 in pos:
                x1, y1 = pos[g1]
                x2, y2 = pos[g2]
                lw = 0.5 + (edge['similarity_score'] / max_sim) * 4
                ax.plot([x1, x2], [y1, y2], '#3498DB', alpha=0.5, linewidth=lw, zorder=1)
    node_colors = []
    node_sizes = []
    for g in selected:
        info = merged_data[g]
        if tnj_key and info.get(tnj_key, False):
            node_colors.append('#2A9D8F')
        else:
            node_colors.append(colors_map.get(info.get(potential_key, 'NONE'), '#95A5A6'))
        node_sizes.append(min(500 + info.get(score_key, 0) * 60, 1800))
    xs = [pos[g][0] for g in selected]
    ys = [pos[g][1] for g in selected]
    ax.scatter(xs, ys, s=node_sizes, c=node_colors, alpha=0.85,
               edgecolors='black', linewidth=1.5, zorder=2)
    for g in selected:
        x, y = pos[g]
        label = g if len(g) <= 25 else g[:22] + '...'
        if tnj_key and merged_data[g].get(tnj_key, False):
            label = '[TNJ] ' + label
        ax.annotate(label, (x, y), xytext=(0, 0), textcoords='offset points',
                    ha='center', va='center', fontsize=9, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85))
    ax.set_title(title)
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect('equal')
    ax.axis('off')
    handles = [Patch(facecolor='#2ECC71', label='HIGH'),
               Patch(facecolor='#F39C12', label='MEDIUM'),
               Patch(facecolor='#E74C3C', label='LOW')]
    if tnj_key:
        handles.append(Patch(facecolor='#2A9D8F', label='TNJ'))
    ax.legend(handles=handles, loc='upper right')
    plt.tight_layout()
    save_figure(fig, output_prefix, 'network')


def plot_pharma_uses(merged_data, output_prefix):
    counts = Counter()
    for info in merged_data.values():
        for u in info.get('pharma_uses', []):
            counts[u] += 1
    if not counts:
        return
    top = counts.most_common(15)
    names, values = zip(*top)
    fig, ax = plt.subplots(figsize=(14, max(6, len(names) * 0.4)))
    bars = ax.barh(range(len(names)), values, color='#9B5DE5', alpha=0.8)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel('Number of associated genes')
    ax.set_title('Predicted Pharmaceutical Applications')
    for bar, v in zip(bars, values):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2, str(v), va='center')
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    save_figure(fig, output_prefix, 'pharmaceutical_uses')


def plot_compound_classes(merged_data, output_prefix):
    counts = Counter()
    for info in merged_data.values():
        for cc in info.get('compound_classes', []):
            counts[cc] += 1
    if not counts:
        return
    colors_map = {'Alkaloids': '#E63946', 'Terpenoids': '#2A9D8F',
                  'Phenylpropanoids': '#E9C46A', 'Flavonoids': '#F4A261',
                  'Stilbenoids': '#E76F51', 'Glucosinolates': '#9B5DE5',
                  'Xanthines': '#F15BB5', 'Cannabinoids': '#FF006E'}
    classes = list(counts.keys())
    values = [counts[c] for c in classes]
    idx = np.argsort(values)
    classes = [classes[i] for i in idx]
    values = [values[i] for i in idx]
    colors = [colors_map.get(c, '#999999') for c in classes]
    fig, ax = plt.subplots(figsize=(14, max(6, len(classes) * 0.5)))
    bars = ax.barh(classes, values, color=colors, alpha=0.8)
    ax.set_xlabel('Number of genes')
    ax.set_title('Distribution of Bioactive Compound Classes')
    for bar, v in zip(bars, values):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2, str(v), va='center')
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    save_figure(fig, output_prefix, 'compound_class_distribution')


def plot_enrichment_dotplot(df_enriched, output_prefix, title, category_col, category_colors, top_n=15):
    if df_enriched is None or len(df_enriched) == 0:
        return
    df_sig = df_enriched[df_enriched['significant']]
    if len(df_sig) == 0:
        return
    top = df_sig.sort_values('fold_enrichment', ascending=False).head(top_n)
    fig, ax = plt.subplots(figsize=(16, max(6, len(top) * 0.4)))
    colors = [category_colors.get(c, '#CCCCCC') for c in top[category_col]] if category_col in top.columns else ['#3498DB'] * len(top)
    ax.scatter(top['fold_enrichment'], top['pfam_domain'],
               s=top['count_in_group'] * 15 + 30, c=colors, alpha=0.8,
               edgecolors='black', linewidth=0.5)
    ax.axvline(x=1, color='red', linestyle='--', alpha=0.5, label='Fold = 1')
    ax.set_xlabel('Fold Enrichment')
    ax.set_ylabel('PFAM Domain')
    ax.set_title(title)
    if category_col in top.columns:
        unique = top[category_col].unique()
        handles = [Patch(facecolor=category_colors.get(c, '#CCCCCC'), label=c) for c in unique if c in category_colors]
        if handles:
            ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.5))
            plt.subplots_adjust(right=0.78)
    plt.tight_layout()
    save_figure(fig, output_prefix, 'enrichment_dotplot')


def plot_tnj_architectures(merged_data, output_prefix):
    arch_counts = Counter()
    for info in merged_data.values():
        if info.get('tnj_detected', False):
            arch_counts[info.get('tnj_architecture', 'Unknown')] += 1
    if not arch_counts:
        return
    names = list(arch_counts.keys())
    values = list(arch_counts.values())
    colors = ['#2A9D8F' if 'canonical' in a.lower() else '#264653' for a in names]
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(names, values, color=colors, alpha=0.8)
    ax.set_xlabel('TNJ Architecture')
    ax.set_ylabel('Number of genes')
    ax.set_title('TNJ Gene Architectures Detected')
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5, str(v), ha='center')
    plt.xticks(rotation=15, ha='right')
    plt.tight_layout()
    save_figure(fig, output_prefix, 'tnj_architectures')


def plot_tnj_confidence(merged_data, output_prefix):
    counts = Counter()
    for info in merged_data.values():
        if info.get('tnj_detected', False):
            counts[info.get('tnj_confidence', 'NONE')] += 1
    if not counts:
        return
    colors_map = {'HIGH': '#2ECC71', 'MEDIUM': '#F39C12', 'LOW': '#E74C3C'}
    names = list(counts.keys())
    values = list(counts.values())
    bar_colors = [colors_map.get(c, '#95A5A6') for c in names]
    fig, ax = plt.subplots(figsize=(8, 6))
    bars = ax.bar(names, values, color=bar_colors, alpha=0.8)
    ax.set_xlabel('Confidence Level')
    ax.set_ylabel('Number of genes')
    ax.set_title('TNJ Gene Confidence Levels')
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5, str(v), ha='center')
    plt.tight_layout()
    save_figure(fig, output_prefix, 'tnj_confidence')


def process_hmmer(gff, genome, hmm_dir, outdir, prefix, domains_dict, pathway_dict,
                  eggnog_data, interpro_data):
    merged = merge_annotations(eggnog_data, interpro_data, pathway_dict)
    if not (gff and genome):
        return merged
    fasta = os.path.join(outdir, f"{prefix}_proteins.fasta")
    if not extract_proteins_from_gff(gff, genome, fasta):
        return merged
    results = run_hmmer(fasta, hmm_dir, outdir, domains_dict, prefix=prefix)
    if results:
        added, updated = integrate_hmmer(merged, results, pathway_dict)
        print(f"HMMER: added={added}, updated={updated}")
    return merged


def run_secondary_analysis(eggnog_file, interpro_file, gff_file, genome_file, hmm_dir, output_dir):
    ensure_dir(output_dir)
    write_session_info(output_dir)
    write_parameters(output_dir, {'eggnog': eggnog_file, 'interpro': interpro_file,
                                  'gff': gff_file, 'genome': genome_file, 'hmm_dir': hmm_dir})
    eggnog_data = parse_eggnog(eggnog_file) if eggnog_file else parse_eggnog(None)
    interpro_data = parse_interpro(interpro_file) if interpro_file else parse_interpro(None)
    merged = process_hmmer(gff_file, genome_file, hmm_dir, output_dir, 'secondary',
                           SECONDARY_HMM, SECONDARY_PATHWAYS, eggnog_data, interpro_data)
    prefix = os.path.join(output_dir, 'secondary')
    plot_ko_dotplot(merged, prefix)
    plot_pathway_completeness(merged, prefix)
    plot_class_distribution(merged, prefix)
    plot_top_ec(merged, prefix)
    plot_top_pfam(merged, prefix)
    print(f"Secondary analysis: {len(merged)} genes processed")


def run_pharma_analysis(eggnog_file, interpro_file, gff_file, genome_file, hmm_dir, output_dir,
                        threshold=4, network_nodes=15, network_method='enrichment'):
    ensure_dir(output_dir)
    write_session_info(output_dir)
    write_parameters(output_dir, {'eggnog': eggnog_file, 'interpro': interpro_file,
                                  'gff': gff_file, 'genome': genome_file, 'hmm_dir': hmm_dir,
                                  'threshold': threshold, 'network_nodes': network_nodes,
                                  'network_method': network_method})
    eggnog_data = parse_eggnog(eggnog_file) if eggnog_file else parse_eggnog(None)
    interpro_data = parse_interpro(interpro_file) if interpro_file else parse_interpro(None)
    merged = process_hmmer(gff_file, genome_file, hmm_dir, output_dir, 'pharma',
                           SECONDARY_HMM, SECONDARY_PATHWAYS, eggnog_data, interpro_data)
    calculate_pharma_score(merged, threshold)
    prefix = os.path.join(output_dir, 'pharma')
    plot_score_ranking(merged, 'pharma_score', 'pharma_potential', prefix, 'Pharmaceutical Potential Score')
    plot_compound_classes(merged, prefix)
    plot_pharma_uses(merged, prefix)
    plot_donut(merged, 'pharma_potential', prefix, 'Distribution of Pharmaceutical Potential')
    edges_df, selected = build_network(merged, 'pharma_score', 'compound_classes',
                                       'secondary_pathways', 'pfam_domains',
                                       max_nodes=network_nodes, method=network_method)
    if selected:
        plot_network(merged, selected, edges_df, 'pharma_score', 'pharma_potential',
                     prefix, 'Functional Network of Pharmaceutical Genes')
    high_med = [g for g, i in merged.items() if i.get('pharma_potential') in ('HIGH', 'MEDIUM')]
    low_none = [g for g, i in merged.items() if i.get('pharma_potential') in ('LOW', 'NONE')]
    if len(high_med) >= 2 and len(low_none) >= 2:
        df = calculate_enrichment(high_med, low_none, merged, min_count=3)
        if len(df) > 0:
            def categorize(pfam):
                for cat, pfams in CORE_PHARMA_DOMAINS.items():
                    if pfam in pfams:
                        return cat
                return 'Other'
            df['pathway_category'] = df['pfam_domain'].apply(categorize)
            cat_colors = {'Flavonoids': '#E76F51', 'Terpenoids': '#2A9D8F',
                          'Alkaloids': '#E63946', 'Phenylpropanoids': '#E9C46A',
                          'Carotenoids': '#F4A261', 'P450': '#9B5DE5',
                          'Methyltransferases': '#F15BB5'}
            plot_enrichment_dotplot(df, prefix, 'Enriched PFAMs in HIGH+MEDIUM Pharma Genes',
                                    'pathway_category', cat_colors)
    print(f"Pharma analysis: {len(merged)} genes processed")


def run_resistance_analysis(eggnog_file, interpro_file, gff_file, genome_file, hmm_dir, output_dir,
                            threshold=5, network_nodes=15, network_method='enrichment'):
    ensure_dir(output_dir)
    write_session_info(output_dir)
    write_parameters(output_dir, {'eggnog': eggnog_file, 'interpro': interpro_file,
                                  'gff': gff_file, 'genome': genome_file, 'hmm_dir': hmm_dir,
                                  'threshold': threshold, 'network_nodes': network_nodes,
                                  'network_method': network_method})
    eggnog_data = parse_eggnog(eggnog_file) if eggnog_file else parse_eggnog(None)
    interpro_data = parse_interpro(interpro_file) if interpro_file else parse_interpro(None)
    merged = process_hmmer(gff_file, genome_file, hmm_dir, output_dir, 'resistance',
                           RESISTANCE_HMM, SECONDARY_PATHWAYS, eggnog_data, interpro_data)
    calculate_resistance_score(merged, threshold)
    prefix = os.path.join(output_dir, 'resistance')
    plot_score_ranking(merged, 'resistance_score', 'resistance_potential', prefix, 'Disease Resistance Score')
    plot_resistance_class_distribution(merged, prefix)
    plot_resistance_pathways(merged, prefix)
    plot_donut(merged, 'resistance_potential', prefix, 'Distribution of Disease Resistance Potential')
    plot_tnj_architectures(merged, prefix)
    plot_tnj_confidence(merged, prefix)
    edges_df, selected = build_network(merged, 'resistance_score', 'resistance_classes',
                                       'resistance_pathways', 'resistance_pfams',
                                       max_nodes=network_nodes, method=network_method)
    if selected:
        plot_network(merged, selected, edges_df, 'resistance_score', 'resistance_potential',
                     prefix, 'Functional Network of Disease Resistance Genes', tnj_key='tnj_detected')
    high_med = [g for g, i in merged.items() if i.get('resistance_potential') in ('HIGH', 'MEDIUM')]
    low_none = [g for g, i in merged.items() if i.get('resistance_potential') in ('LOW', 'NONE')]
    if len(high_med) >= 2 and len(low_none) >= 2:
        df = calculate_enrichment(high_med, low_none, merged, min_count=3)
        if len(df) > 0:
            def categorize(pfam):
                for cat, info in DISEASE_RESISTANCE_PFAMS.items():
                    if pfam in info['pfams']:
                        return cat
                return 'Other'
            df['pathway_category'] = df['pfam_domain'].apply(categorize)
            colors = {c: info['color'] for c, info in DISEASE_RESISTANCE_PFAMS.items()}
            plot_enrichment_dotplot(df, prefix, 'Enriched PFAMs in HIGH+MEDIUM Resistance Genes',
                                    'pathway_category', colors)
    print(f"Resistance analysis: {len(merged)} genes processed")


def run_eggnog2kegg(eggnog_file, output_file, output_dir):
    ensure_dir(output_dir)
    path = os.path.join(output_dir, output_file)
    total, with_ko = parse_eggnog_kegg_mapper(eggnog_file, path)
    print(f"Total: {total}, with KO: {with_ko} ({with_ko / total * 100:.2f}%)")


class App:
    def __init__(self, root):
        self.root = root
        self.root.title("annoMining v4.0")
        self.root.geometry("1050x950")
        style = ttk.Style()
        style.theme_use('clam')
        style.configure('TNotebook.Tab', padding=[15, 8], font=('Segoe UI', 11, 'bold'))
        style.map('TNotebook.Tab', background=[('selected', '#2563eb')],
                  foreground=[('selected', 'white')])
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        self.build_secondary_tab()
        self.build_pharma_tab()
        self.build_resistance_tab()
        self.build_eggnog_tab()

    def _row(self, parent, row, label, attr, browse='file'):
        tk.Label(parent, text=label).grid(row=row, column=0, sticky='e', padx=5, pady=4)
        entry = tk.Entry(parent, width=50)
        entry.grid(row=row, column=1, padx=5, pady=4)
        setattr(self, attr, entry)
        cmd = (lambda: self.browse_dir(entry)) if browse == 'dir' else (lambda: self.browse_file(entry))
        tk.Button(parent, text="Browse", command=cmd).grid(row=row, column=2, padx=5)

    def browse_file(self, entry):
        p = filedialog.askopenfilename()
        if p:
            entry.delete(0, tk.END)
            entry.insert(0, p)

    def browse_dir(self, entry):
        p = filedialog.askdirectory()
        if p:
            entry.delete(0, tk.END)
            entry.insert(0, p)

    def log(self, widget, msg):
        widget.insert(tk.END, msg + "\n")
        widget.see(tk.END)

    def run_thread(self, target, log_widget, success_msg):
        def wrapper():
            try:
                self.log(log_widget, "Starting analysis...")
                target()
                self.log(log_widget, f"Completed: {success_msg}")
                messagebox.showinfo("Success", success_msg)
            except Exception as e:
                self.log(log_widget, f"Error: {e}")
                self.log(log_widget, traceback.format_exc())
                messagebox.showerror("Error", str(e))
        threading.Thread(target=wrapper, daemon=True).start()

    def build_secondary_tab(self):
        tab = tk.Frame(self.notebook)
        self.notebook.add(tab, text="Secondary Metabolism")
        c = tk.Frame(tab, padx=20, pady=20)
        c.pack(fill=tk.BOTH, expand=True)
        tk.Label(c, text="Secondary Metabolism Analysis", font=('Segoe UI', 14, 'bold')).grid(
            row=0, column=0, columnspan=3, pady=(0, 20), sticky='w')
        self._row(c, 1, "eggNOG file:", 'sec_eggnog')
        self._row(c, 2, "InterProScan file:", 'sec_interpro')
        self._row(c, 3, "GFF file:", 'sec_gff')
        self._row(c, 4, "Genome FASTA:", 'sec_genome')
        self._row(c, 5, "HMM profiles dir:", 'sec_hmm', browse='dir')
        tk.Label(c, text="Output folder:").grid(row=6, column=0, sticky='e', padx=5, pady=4)
        self.sec_outdir = tk.Entry(c, width=30)
        self.sec_outdir.insert(0, "secondary_metabolism")
        self.sec_outdir.grid(row=6, column=1, sticky='w', padx=5, pady=4)
        tk.Button(c, text="Run Analysis", command=self.run_secondary,
                  font=('Segoe UI', 12, 'bold'), padx=30, pady=8).grid(row=7, column=0, columnspan=3, pady=15)
        self.sec_log = scrolledtext.ScrolledText(c, height=8, bg='#1e293b', fg='#e2e8f0', font=('Consolas', 9))
        self.sec_log.grid(row=8, column=0, columnspan=3, padx=10, pady=10, sticky='we')

    def build_pharma_tab(self):
        tab = tk.Frame(self.notebook)
        self.notebook.add(tab, text="Pharma Analysis")
        c = tk.Frame(tab, padx=20, pady=20)
        c.pack(fill=tk.BOTH, expand=True)
        tk.Label(c, text="Pharmaceutical Analysis", font=('Segoe UI', 14, 'bold')).grid(
            row=0, column=0, columnspan=3, pady=(0, 20), sticky='w')
        self._row(c, 1, "eggNOG file:", 'ph_eggnog')
        self._row(c, 2, "InterProScan file:", 'ph_interpro')
        self._row(c, 3, "GFF file:", 'ph_gff')
        self._row(c, 4, "Genome FASTA:", 'ph_genome')
        self._row(c, 5, "HMM profiles dir:", 'ph_hmm', browse='dir')
        tk.Label(c, text="Threshold:").grid(row=6, column=0, sticky='e', padx=5, pady=4)
        self.ph_threshold = tk.Entry(c, width=10)
        self.ph_threshold.insert(0, "4")
        self.ph_threshold.grid(row=6, column=1, sticky='w', padx=5, pady=4)
        tk.Label(c, text="Network nodes:").grid(row=7, column=0, sticky='e', padx=5, pady=4)
        self.ph_nodes = tk.Entry(c, width=10)
        self.ph_nodes.insert(0, "15")
        self.ph_nodes.grid(row=7, column=1, sticky='w', padx=5, pady=4)
        tk.Label(c, text="Network method:").grid(row=8, column=0, sticky='e', padx=5, pady=4)
        self.ph_method = ttk.Combobox(c, values=['score', 'enrichment', 'balanced'], width=15)
        self.ph_method.set('enrichment')
        self.ph_method.grid(row=8, column=1, sticky='w', padx=5, pady=4)
        tk.Label(c, text="Output folder:").grid(row=9, column=0, sticky='e', padx=5, pady=4)
        self.ph_outdir = tk.Entry(c, width=30)
        self.ph_outdir.insert(0, "pharma_analysis")
        self.ph_outdir.grid(row=9, column=1, sticky='w', padx=5, pady=4)
        tk.Button(c, text="Run Analysis", command=self.run_pharma,
                  font=('Segoe UI', 12, 'bold'), padx=30, pady=8).grid(row=10, column=0, columnspan=3, pady=15)
        self.ph_log = scrolledtext.ScrolledText(c, height=8, bg='#1e293b', fg='#e2e8f0', font=('Consolas', 9))
        self.ph_log.grid(row=11, column=0, columnspan=3, padx=10, pady=10, sticky='we')

    def build_resistance_tab(self):
        tab = tk.Frame(self.notebook)
        self.notebook.add(tab, text="Disease Resistance")
        c = tk.Frame(tab, padx=20, pady=20)
        c.pack(fill=tk.BOTH, expand=True)
        tk.Label(c, text="Disease Resistance Analysis", font=('Segoe UI', 14, 'bold')).grid(
            row=0, column=0, columnspan=3, pady=(0, 20), sticky='w')
        self._row(c, 1, "eggNOG file:", 'res_eggnog')
        self._row(c, 2, "InterProScan file:", 'res_interpro')
        self._row(c, 3, "GFF file:", 'res_gff')
        self._row(c, 4, "Genome FASTA:", 'res_genome')
        self._row(c, 5, "HMM profiles dir:", 'res_hmm', browse='dir')
        tk.Label(c, text="Threshold:").grid(row=6, column=0, sticky='e', padx=5, pady=4)
        self.res_threshold = tk.Entry(c, width=10)
        self.res_threshold.insert(0, "5")
        self.res_threshold.grid(row=6, column=1, sticky='w', padx=5, pady=4)
        tk.Label(c, text="Network nodes:").grid(row=7, column=0, sticky='e', padx=5, pady=4)
        self.res_nodes = tk.Entry(c, width=10)
        self.res_nodes.insert(0, "15")
        self.res_nodes.grid(row=7, column=1, sticky='w', padx=5, pady=4)
        tk.Label(c, text="Network method:").grid(row=8, column=0, sticky='e', padx=5, pady=4)
        self.res_method = ttk.Combobox(c, values=['score', 'enrichment', 'balanced'], width=15)
        self.res_method.set('enrichment')
        self.res_method.grid(row=8, column=1, sticky='w', padx=5, pady=4)
        tk.Label(c, text="Output folder:").grid(row=9, column=0, sticky='e', padx=5, pady=4)
        self.res_outdir = tk.Entry(c, width=30)
        self.res_outdir.insert(0, "resistance_analysis")
        self.res_outdir.grid(row=9, column=1, sticky='w', padx=5, pady=4)
        tk.Button(c, text="Run Analysis", command=self.run_resistance,
                  font=('Segoe UI', 12, 'bold'), padx=30, pady=8).grid(row=10, column=0, columnspan=3, pady=15)
        self.res_log = scrolledtext.ScrolledText(c, height=8, bg='#1e293b', fg='#e2e8f0', font=('Consolas', 9))
        self.res_log.grid(row=11, column=0, columnspan=3, padx=10, pady=10, sticky='we')

    def build_eggnog_tab(self):
        tab = tk.Frame(self.notebook)
        self.notebook.add(tab, text="eggnog2kegg")
        c = tk.Frame(tab, padx=20, pady=20)
        c.pack(fill=tk.BOTH, expand=True)
        tk.Label(c, text="eggNOG to KEGG Converter", font=('Segoe UI', 14, 'bold')).grid(
            row=0, column=0, columnspan=3, pady=(0, 20), sticky='w')
        self._row(c, 1, "eggNOG file:", 'ek_input')
        tk.Label(c, text="Output file:").grid(row=2, column=0, sticky='e', padx=5, pady=4)
        self.ek_output = tk.Entry(c, width=30)
        self.ek_output.insert(0, "kegg_mapper.txt")
        self.ek_output.grid(row=2, column=1, sticky='w', padx=5, pady=4)
        tk.Label(c, text="Output folder:").grid(row=3, column=0, sticky='e', padx=5, pady=4)
        self.ek_outdir = tk.Entry(c, width=30)
        self.ek_outdir.insert(0, "eggnog2kegg_output")
        self.ek_outdir.grid(row=3, column=1, sticky='w', padx=5, pady=4)
        tk.Button(c, text="Convert", command=self.run_ek,
                  font=('Segoe UI', 12, 'bold'), padx=30, pady=8).grid(row=4, column=0, columnspan=3, pady=15)
        self.ek_log = scrolledtext.ScrolledText(c, height=8, bg='#1e293b', fg='#e2e8f0', font=('Consolas', 9))
        self.ek_log.grid(row=5, column=0, columnspan=3, padx=10, pady=10, sticky='we')

    def run_secondary(self):
        eggnog = self.sec_eggnog.get().strip()
        interpro = self.sec_interpro.get().strip()
        gff = self.sec_gff.get().strip()
        genome = self.sec_genome.get().strip()
        hmm = self.sec_hmm.get().strip()
        outdir = self.sec_outdir.get().strip()
        if not (eggnog or interpro or gff):
            messagebox.showerror("Error", "Provide at least eggNOG, InterPro or GFF file")
            return
        if (gff and not genome) or (genome and not gff):
            messagebox.showerror("Error", "Both GFF and Genome are required together")
            return
        self.run_thread(lambda: run_secondary_analysis(eggnog, interpro, gff, genome, hmm, outdir),
                        self.sec_log, f"Secondary analysis done: {outdir}")

    def run_pharma(self):
        eggnog = self.ph_eggnog.get().strip()
        interpro = self.ph_interpro.get().strip()
        gff = self.ph_gff.get().strip()
        genome = self.ph_genome.get().strip()
        hmm = self.ph_hmm.get().strip()
        outdir = self.ph_outdir.get().strip()
        try:
            threshold = int(self.ph_threshold.get().strip() or 4)
            nodes = int(self.ph_nodes.get().strip() or 15)
        except ValueError:
            messagebox.showerror("Error", "Threshold and nodes must be integers")
            return
        method = self.ph_method.get().strip() or 'enrichment'
        if not (eggnog or interpro or gff):
            messagebox.showerror("Error", "Provide at least eggNOG, InterPro or GFF file")
            return
        if (gff and not genome) or (genome and not gff):
            messagebox.showerror("Error", "Both GFF and Genome are required together")
            return
        self.run_thread(lambda: run_pharma_analysis(eggnog, interpro, gff, genome, hmm, outdir,
                                                     threshold, nodes, method),
                        self.ph_log, f"Pharma analysis done: {outdir}")

    def run_resistance(self):
        eggnog = self.res_eggnog.get().strip()
        interpro = self.res_interpro.get().strip()
        gff = self.res_gff.get().strip()
        genome = self.res_genome.get().strip()
        hmm = self.res_hmm.get().strip()
        outdir = self.res_outdir.get().strip()
        try:
            threshold = int(self.res_threshold.get().strip() or 5)
            nodes = int(self.res_nodes.get().strip() or 15)
        except ValueError:
            messagebox.showerror("Error", "Threshold and nodes must be integers")
            return
        method = self.res_method.get().strip() or 'enrichment'
        if not (eggnog or interpro or gff):
            messagebox.showerror("Error", "Provide at least eggNOG, InterPro or GFF file")
            return
        if (gff and not genome) or (genome and not gff):
            messagebox.showerror("Error", "Both GFF and Genome are required together")
            return
        self.run_thread(lambda: run_resistance_analysis(eggnog, interpro, gff, genome, hmm, outdir,
                                                         threshold, nodes, method),
                        self.res_log, f"Resistance analysis done: {outdir}")

    def run_ek(self):
        infile = self.ek_input.get().strip()
        outfile = self.ek_output.get().strip() or 'kegg_mapper.txt'
        outdir = self.ek_outdir.get().strip() or 'eggnog2kegg_output'
        if not infile or not os.path.exists(infile):
            messagebox.showerror("Error", "Invalid input file")
            return
        self.run_thread(lambda: run_eggnog2kegg(infile, outfile, outdir),
                        self.ek_log, f"Conversion done: {os.path.join(outdir, outfile)}")


def main():
    configure_style()
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == '__main__':
    main()
