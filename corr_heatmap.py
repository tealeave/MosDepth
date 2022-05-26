import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import sys

"""This script helps to plot pwise correlation heatmaps from per target coverage xlsx"""


def read_df(cov_xlsx):

    """this read per target coverage file to df"""

    df = pd.read_excel(cov_xlsx, sheet_name="All_Samples")
    # if dealing with pwise correlation in 2 df then sort to make sure columns are in order
    df = df.reindex(sorted(df.columns), axis=1)
    df.set_index("target", inplace=True)

    return df


def pwise_corr(df):

    """this function plots pwise correlation. figure size 100 for Nextseq, 200 for Nova384,
     300 for Nova768, 400 will cause memory error"""

    sns.set(rc={"figure.figsize": (100, 100)})

    sns_plot = sns.heatmap(
        df.corr(),
        #  annot = True,
        cmap="coolwarm",
        square=True,
    )

    plt.savefig("correlation_heatmap.png")


def pwise_corr_2runs(nova_df, next_df):

    """ "this function take two pt_cov df and compare to sample cov vector between them"""

    series_lst = []
    for noi, noc in nova_df.iteritems():
        new_series_name = noc.name
        pearson_lst = []
        # storing row names
        index_lst = []
        #     iterate oveer next_df columns for pearson
        for nei, nec in next_df.iteritems():
            # only take a series if it is also in Novaseq
            if nec.name in nova_df.columns:
                pearson = nec.corr(noc)
                pearson_lst.append(pearson)
                index_lst.append(nec.name)
        new_series = pd.Series(pearson_lst, index=index_lst, name=new_series_name)
        series_lst.append(new_series)
    # this Novaseq as series and Nextseq as row indices
    NovaVsNextdf = pd.concat(series_lst, axis=1)
    # Sort the df by mean value of a column
    NovaVsNextdfSorted = NovaVsNextdf.reindex(
        NovaVsNextdf.mean().sort_values().index, axis=1
    )
    print(NovaVsNextdfSorted.head())
    sns.set(rc={"figure.figsize": (300, 300)})
    # Novaseq as x axis and nextseq as y axis
    sns_plot = sns.heatmap(NovaVsNextdfSorted, annot=True, cmap="coolwarm", square=True)

    plt.savefig("NovaXVsNextY.png")


def run():
    pt_cov_xlsx = sys.argv[1]
    df = read_df(pt_cov_xlsx)
    pwise_corr(df)

    # if comparing same set of sample in two diff runs
    # pt_cov_xlsx1 = sys.argv[1]
    # pt_cov_xlsx2 = sys.argv[2]
    # df1 = read_df(pt_cov_xlsx1)
    # df2 = read_df(pt_cov_xlsx2)
    # pwise_corr_2runs(df1, df2)


if __name__ == "__main__":
    run()
