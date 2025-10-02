#!/usr/bin/env python3
# coding: utf-8

"""
@author      CS Lim
@create date 2025-08-18 17:23:00
@modify date 2025-08-18 18:07:41
@desc        RIBOSS module for generating flatten, ORF-centric GTF files
"""


import pandas as pd
import pyranges as pr
import csv
from .orfs import orf_finder


def flatten_orfs(df_or_path, gtf_path):
    """
    Flatten overlapping ORF annotations into non-redundant BED6 intervals.

    Args:
        df_or_path (Union[pd.DataFrame, str]): Either a pandas DataFrame or a path to a tab-delimited or pickle file.

    Returns:
        pd.DataFrame: BED6-format DataFrame with columns:
            Chromosome, Start, End, Name, Score, Strand
    """
        
    if isinstance(df_or_path, str):
        df = pd.read_pickle(df_or_path) if df_or_path.endswith(".pkl") else pd.read_csv(df_or_path, sep="\t")
    else:
        df = df_or_path.copy()

    gr = df.iloc[:, [0, 7, 8, 1, 9, 2]].copy()
    gr.columns = ["Chromosome", "Start", "End", "transcript_id", "Score", "Strand"]

    gtf = pd.read_csv(gtf_path, sep="\t", header=None, comment="#")
    gtf = gtf[gtf[2] == "transcript"].copy()
    gtf["transcript_id"] = gtf[8].str.replace('"', '').str.replace(";", "").str.split().str[3]
    gtf["Name"] = gtf[8].str.replace('"', '').str.replace(";", "").str.split().str[1]
    gt_map = gtf[["Name", "transcript_id"]].drop_duplicates()

    gr = gr.merge(gt_map, on="transcript_id", how="left").drop("transcript_id", axis=1)
    gr.drop_duplicates(inplace=True)
    gr.sort_values(by=["Chromosome", "Start", "Name"], inplace=True)
    gr = pr.PyRanges(gr)

    merged = gr.merge(strand=True)
    joined = merged.join(gr, strandedness="same")

    agg = joined.df.groupby(["Chromosome", "Start", "End", "Strand"]).agg({
        "Name": lambda x: sorted(set(x))[0],
        "Score": lambda x: ",".join(sorted(set(x)))
    }).reset_index()

    agg = agg[~agg["Score"].str.contains("dORF,uORF")]

    morfs = pr.PyRanges(agg[agg.Score.str.contains("mORF")])
    xorfs = pr.PyRanges(agg[~agg.Score.str.contains("mORF")])
    remove = xorfs.overlap(morfs, strandedness=False)
    xorfs = pd.concat([xorfs.df, remove.df]).drop_duplicates(keep=False)
    
    morfs = morfs.df
    morfs["Score"] = "mORF"
    
    orf_df = pd.concat([morfs, xorfs]).sort_values(by=["Chromosome", "Start", "Name"])
    orf_df = orf_df[["Chromosome", "Start", "End", "Name", "Score", "Strand"]].copy()

    gt_map.columns = ["name", "transcript_id"]
    gt_map.sort_values(by=["name", "transcript_id"], inplace=True)

    return orf_df, gt_map


def bed6_to_bed12(df):
    """
    Convert BED6-format DataFrame to BED12 format.

    Args:
        df (pd.DataFrame): BED6-format DataFrame.

    Returns:
        pd.DataFrame: BED12-format DataFrame.
    """
    
    bed12_rows = []
    for name, group in df.groupby("Name"): # df.groupby(["Name", "Strand"], observed=True):
        group = group.sort_values("Start")
        chrom = group.iloc[0]["Chromosome"]
        strand = group.iloc[0]["Strand"]
        tx_start = group["Start"].min()
        tx_end = group["End"].max()
        block_sizes = (group["End"] - group["Start"]).tolist()
        block_starts = (group["Start"] - tx_start).tolist()
        block_count = len(block_sizes)

        bed12_rows.append([
            chrom, tx_start, tx_end, name, 0, strand,
            tx_start, tx_end, "255,0,0", block_count,
            ",".join(map(str, block_sizes)) + ",",
            ",".join(map(str, block_starts)) + ","
        ])

    return pd.DataFrame(bed12_rows, columns=[
        "chrom", "chromStart", "chromEnd", "name", "score", "strand",
        "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"
    ])


def bed12_to_gtf(bed12_df, gene_to_transcript_map=None):
    """
    Convert BED12-format DataFrame to DOTSeq/DEXSeq-compatible GTF format.

    Args:
        bed12_df (pd.DataFrame): BED12-format DataFrame.
        gene_to_transcript_map (pd.DataFrame): Gene to transcript ID mapping DataFrame.
        
    Returns:
        pd.DataFrame: GTF-format DataFrame.
    """
        
    gff_rows = []

    if isinstance(gene_to_transcript_map, pd.DataFrame):
        gt_map = gene_to_transcript_map.drop_duplicates("name")
        bed12_df = pd.merge(bed12_df, gt_map)
    else:
        bed12_df["transcript_id"] = bed12_df["name"]
        bed12_df["name"] = bed12_df["name"].str.split(".").str[0]

    for _, row in bed12_df.iterrows():
        chrom = row["chrom"]
        strand = row["strand"]
        transcript_id = row["transcript_id"]
        gene_id = row["name"]
        start = int(row["chromStart"])
        end = int(row["chromEnd"])
        block_sizes = list(map(int, row["blockSizes"].strip(",").split(",")))
        block_starts = list(map(int, row["blockStarts"].strip(",").split(",")))

        gff_rows.append([chrom, "orf_finder", "gene", start + 1, end, ".", strand, ".", f'gene_id "{gene_id}";'])

        for i, (bs, rs) in enumerate(zip(block_sizes, block_starts)):
            exon_start = start + rs
            exon_end = exon_start + bs
            attribute = f'gene_id "{gene_id}"; transcripts "{transcript_id}"; exon_number "{i+1:03d}";'
            gff_rows.append([chrom, "orf_finder", "exon", exon_start + 1, exon_end, ".", strand, ".", attribute])

    return pd.DataFrame(gff_rows, columns=[
        "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"
    ])

