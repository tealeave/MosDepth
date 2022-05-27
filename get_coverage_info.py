#!/usr/bin/python

import argparse
from pathlib import Path
import pandas as pd
import logging
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess
import shutil

"""
    This script is used to generate coverage analysis from bam files and a target bed file

"""
__version__ = "1.0"
log_level = logging.DEBUG
logging.basicConfig(
    filename=None,
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=log_level,
)


def get_args():
    """Parse Arguments"""

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--bam_list",
        action="store",
        dest="bam_list",
        required=True,
        help="[Input] list of bam fils for cov analysis",
    )
    parser.add_argument(
        "--bed",
        action="store",
        dest="bed",
        required=True,
        help="[Input] target bed file where 4th column has the $gene_exon info",
    )

    parser.add_argument(
        "--heatmap",
        action="store_true",
        default=False,
        help="[Optional] to genrate correlation heatmap",
    )

    parser.add_argument(
        "--outdir",
        action="store",
        dest="outdir",
        type=Path,
        default=Path.cwd(),
        help="[Optional] Path to output directory",
    )

    return parser.parse_args()


def read_lst(txt):
    lst = []
    with open(txt) as f:
        for line in f:
            lst.append(line.strip())
    return lst


def run_mosdepth(bam_list, bed, inter_folder):

    """This function runs mosdepth jobs per sample and produce preliminary results in inter_folder"""

    bam_lst = read_lst(bam_list)
    work_dir = Path.cwd()

    # # Multi-processing
    cmds_list = [
        [
            "/home/dlin/miniconda3/envs/mose/bin/mosdepth",
            "--by",
            f"{bed}",
            "-T",
            "10,20,50,100,500",
            "-x",
            f"{Path(bam).stem}",
            f"{bam}",
        ]
        for bam in bam_lst
    ]
    logging.debug(cmds_list)
    procs_list = [
        subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )
        for cmd in cmds_list
    ]

    for proc in procs_list:
        proc.wait()

    # Write to log files
    for i in range(len(bam_lst)):
        mosdepth_log = work_dir.joinpath(f"{Path(bam_lst[i]).stem}_mosdepth.log")
        mosdepth_log.touch()
        with open(f"{mosdepth_log}", "w") as log_file:
            log_file.write(procs_list[i].stdout.read())

    file2move = []
    file_types = ("*.gz", "*.txt", "*.csi", "*_mosdepth.log")
    for file in file_types:
        file2move.extend(work_dir.glob(file))
    logging.debug(file2move)

    for f in file2move:
        source = f
        destinaiton = inter_folder.joinpath(f"{Path(f).name}")
        shutil.move(source, destinaiton)


def read_region(region_lst, summary_lst):

    """this function reads and compiles per target coverage output from MosDepth"""

    df_lst_for_avg = []
    df_lst_for_all_sample = []
    for r in region_lst:
        sample_name = r.stem.split(".")[0]
        # logging.debug(sample_name)
        df = pd.read_csv(
            r,
            compression="gzip",
            sep="\t",
            header=None,
            names=["chr", "start", "stop", "target", "average"],
        )
        tdf = df[["target", "average"]]
        df_lst_for_avg.append(tdf)
        sdf = df[["target", "average"]]
        sdf.columns = ["target", f"{sample_name}"]
        sdf.set_index("target", inplace=True)
        df_lst_for_all_sample.append(sdf)

    avg_mdf = pd.concat(df_lst_for_avg, axis=0)
    gdf = avg_mdf.groupby("target").mean()
    gdf.reset_index(inplace=True)

    sample_mdf = pd.concat(df_lst_for_all_sample, axis=1)

    sdf_lst = []
    for s in summary_lst:
        sample_name = s.stem.split(".")[0]
        df = pd.read_csv(s, sep="\t", header=0)
        df["sample_name"] = sample_name
        df = df[df["chrom"] == "total_region"]
        df = df[["sample_name", "length", "mean", "min", "max"]]
        sdf_lst.append(df)
    smdf = pd.concat(sdf_lst, axis=0)

    return sample_mdf, gdf, smdf


def read_threshold(threshold_lst):

    """This function reads and compiles .threshold.bed.gz file from the mosdepth output"""

    df_lst = []
    sdf_lst = []
    for t in threshold_lst:
        sample_name = t.stem.split(".")[0]
        df = pd.read_csv(t, compression="gzip", sep="\t", header=0)
        df["sample_name"] = sample_name
        df["region_size"] = df["end"] - df["start"]

        sdf = df.groupby("sample_name").sum()
        sdf["%Over10X"] = sdf["10X"] / sdf["region_size"]
        sdf["%Over20X"] = sdf["20X"] / sdf["region_size"]
        sdf["%Over50X"] = sdf["50X"] / sdf["region_size"]
        sdf["%Over100X"] = sdf["100X"] / sdf["region_size"]
        sdf["%Over500X"] = sdf["500X"] / sdf["region_size"]
        sdf = sdf[["%Over10X", "%Over20X", "%Over50X", "%Over100X", "%Over500X"]]
        # print(sdf)
        sdf_lst.append(sdf)

        df["%Over10X"] = df["10X"] / df["region_size"]
        df["%Over20X"] = df["20X"] / df["region_size"]
        df["%Over50X"] = df["50X"] / df["region_size"]
        df["%Over100X"] = df["100X"] / df["region_size"]
        df["%Over500X"] = df["500X"] / df["region_size"]

        df = df[
            [
                "region",
                "%Over10X",
                "%Over20X",
                "%Over50X",
                "%Over100X",
                "%Over500X",
            ]
        ]
        # print(df.head())

        df_lst.append(df)

        # To get to one sample per row
        # des_df = df.describe().reset_index()
        # des_df = des_df[des_df['index'] == 'mean']
        # des_df['sample'] = sample
        # des_df = des_df[['sample', '10Xp', '20Xp', '50Xp', '100Xp', '150Xp']]
        # df_lst.append(des_df)

    mdf = pd.concat(df_lst, axis=0)
    gdf = mdf.groupby("region").mean()
    gdf.reset_index(inplace=True)

    msdf = pd.concat(sdf_lst, axis=0)
    msdf.reset_index(inplace=True)

    return gdf, msdf


def pwise_corr(outdir, df):

    """this function plots pwise correlation. figure size 100 for Nextseq, 200 for Nova384,
    300 for Nova768, 400 will cause memory error"""

    sns.set(rc={"figure.figsize": (100, 100)})

    sns_plot = sns.heatmap(
        df.corr(),
        #  annot = True,
        cmap="coolwarm",
        square=True,
    )

    plt.savefig(outdir.joinpath("correlation_heatmap.png"))


def run():
    args = get_args()
    # logging.debug(args.outdir)
    intermediate_folder = args.outdir.joinpath("intermediate_files")
    intermediate_folder.mkdir(parents=True, exist_ok=True)
    logging.info("running mosdepth!")
    run_mosdepth(args.bam_list, args.bed, intermediate_folder)

    logging.info("compiling per target and threshold result!")
    mdf, gdf, smdf = read_region(
        intermediate_folder.glob("*.regions.bed.gz"),
        intermediate_folder.glob("*.mosdepth.summary.txt"),
    )
    writer = pd.ExcelWriter(
        args.outdir.joinpath("PerTargetMeanCov.xlsx"), engine="xlsxwriter"
    )
    mdf.to_excel(writer, sheet_name="All_Samples", index=True)
    gdf.to_excel(writer, sheet_name="AvgAcrossSamples", index=False)
    smdf.to_excel(writer, sheet_name="SampleLevel", index=False)
    writer.save()

    tgdf, tmsdf = read_threshold(intermediate_folder.glob("*.thresholds.bed.gz"))
    writer = pd.ExcelWriter(
        args.outdir.joinpath("PerTargetThresholdCov.xlsx"), engine="xlsxwriter"
    )
    tgdf.to_excel(writer, sheet_name="AvgAcrossSamples", index=False)
    tmsdf.to_excel(writer, sheet_name="SampleLevel", index=False)
    writer.save()

    if args.heatmap:
        mdf = mdf.reindex(sorted(mdf.columns), axis=1)
        pwise_corr(args.outdir, mdf)

    

    logging.debug("All done!")


if __name__ == "__main__":
    run()
