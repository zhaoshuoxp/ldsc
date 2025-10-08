#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import argparse
import gzip
from pybedtools import BedTool

def read_gene_coords(gene_coord_file, gene_set_file=None):
    print("Reading gene coordinates...")
    df_gene = pd.read_csv(gene_coord_file, delim_whitespace=True)
    if gene_set_file:
        with open(gene_set_file) as f:
            gene_set = set(line.strip() for line in f)
        df_gene = df_gene[df_gene['GENE'].isin(gene_set)]
    bed_dict = {}
    for _, row in df_gene.iterrows():
        chrom = str(row['CHR']).replace('chr','')
        start = int(row['START'])
        end = int(row['END'])
        bed_dict[row['GENE']] = BedTool([[chrom, start, end, row['GENE']]])
    return bed_dict


def make_annot_file(bimfile, bed_dict=None, annot_file=None, windowsize=0, merge=True, bed_file=None):
    print("Reading BIM file...")
    df_bim = pd.read_csv(bimfile, sep="\t", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    df_bim['CHR'] = df_bim['CHR'].astype(str).str.replace('chr','')
    df_bim['base'] = 1

    print("Building SNP BedTool...")
    snp_bed = BedTool([
        [str(int(x[0])), int(x[1]), int(x[1])+1, str(x[3])]
        for x in df_bim[['CHR','BP','BP','SNP']].values
    ])

    if bed_file:
        print("Using BED file for intersection...")
        bedtool_all = BedTool(bed_file)
    else:
        all_bed_list = []
        for gene, bed in bed_dict.items():
            for f in bed:
                start = max(0, f.start - windowsize)
                end = f.end + windowsize
                all_bed_list.append([f.chrom, start, end, gene])
        bedtool_all = BedTool(all_bed_list)

    genes = list(set(f.fields[3] for f in bedtool_all))

    if merge:
        print("Merge mode: combining overlapping genes...")
        all_bed = bedtool_all.sort().merge(c=[4], o=['distinct'])
        intersect = snp_bed.intersect(all_bed, wa=True)
        snp_hits = set([x.name for x in intersect])
        annot_df = pd.DataFrame({
            'CHR': df_bim['CHR'],
            'BP': df_bim['BP'],
            'SNP': df_bim['SNP'],
            'CM': df_bim['CM'],
            'base': df_bim['base'],
            'ANNOT': df_bim['SNP'].apply(lambda x: 1 if x in snp_hits else 0)
        })
    else:
        print("Nomerged mode: calculating proportion per gene based on overlap length...")
        annot_dict = {gene: np.zeros(len(df_bim)) for gene in genes}
        snp_idx_map = {snp:i for i,snp in enumerate(df_bim['SNP'])}

        intersect = snp_bed.intersect(bedtool_all, wa=True, wb=True)

        snp_to_gene_lengths = {}
        for f in intersect:
            snp_name = f.name
            gene_name = f.fields[7]
            start_overlap = max(int(f.start), int(f.fields[4]))
            end_overlap = min(int(f.end), int(f.fields[5]))
            length = max(0, end_overlap - start_overlap)
            snp_to_gene_lengths.setdefault(snp_name, []).append((gene_name, length))

        for snp_name, idx in snp_idx_map.items():
            overlaps = snp_to_gene_lengths.get(snp_name, [])
            if not overlaps:
                continue
            if len(overlaps) == 1:
                gene_name, _ = overlaps[0]
                annot_dict[gene_name][idx] = 1.0
            else:
                total_len = sum(l for _, l in overlaps)
                if total_len == 0:
                    continue
                for gene_name, l in overlaps:
                    annot_dict[gene_name][idx] = l / total_len

        base_df = df_bim[['CHR','BP','SNP','CM','base']].copy()
        gene_df = pd.DataFrame(annot_dict)
        annot_df = pd.concat([base_df, gene_df], axis=1)

    print("Writing annot file...")
    if annot_file.endswith('.gz'):
        with gzip.open(annot_file, 'wt') as f:
            annot_df.to_csv(f, sep='\t', index=False)
    else:
        annot_df.to_csv(annot_file, sep='\t', index=False)


    print("Fixing column order...")
    expected_cols = ['CHR', 'BP', 'SNP', 'CM', 'base']
    df = pd.read_csv(annot_file, sep='\t', compression='gzip' if annot_file.endswith('.gz') else None)
    annot_cols = [c for c in df.columns if c not in expected_cols]
    df = df[expected_cols + annot_cols]

    if annot_file.endswith('.gz'):
        df.to_csv(annot_file, sep='\t', index=False, compression='gzip')
    else:
        df.to_csv(annot_file, sep='\t', index=False)

    print(f"✅ Output file updated with correct column order: {annot_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--bimfile', required=True)
    parser.add_argument('--gene-coord-file', default=None)
    parser.add_argument('--gene-set-file', default=None)
    parser.add_argument('--windowsize', type=int, default=0)
    parser.add_argument('--annot-file', required=True)
    parser.add_argument('--nomerge', action='store_true')
    parser.add_argument('--bed-file', default=None, help="Direct BED file for intersection, keeps merge/nomerge logic")
    args = parser.parse_args()

    bed_dict = None
    if args.bed_file is None and args.gene_coord_file:
        bed_dict = read_gene_coords(args.gene_coord_file, args.gene_set_file)

    make_annot_file(
        args.bimfile,
        bed_dict,
        args.annot_file,
        windowsize=args.windowsize,
        merge=not args.nomerge,
        bed_file=args.bed_file
    )
