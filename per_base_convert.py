import pandas as pd
import sys, os

"""This script converts per base coverage output from Mosdepth to the format 
    that is used in GATK"""

bedtools_output = sys.argv[1]

def read_bedtools_output(pb):

    """This function reads bedtools ouput and do the inial parsing and cleaning"""

    df = pd.read_csv(pb, sep="\t", names=["chr", "start", "stop", "cov"])
    df["length"] = df["stop"] - df["start"]
    cleaned_df1 = df[
        ((df["chr"] != "artifact_seq") & (df["length"] <= 300))
        # & (df['cov']<20)
        # & (df['length'] > 1)
    ]

    return cleaned_df1


def adding_per_base(df):

    """This function iterate through MosDepth per base cov range result and generate a per base resolutoin coverage file"""

    cord_lst = []
    for i, row in df.iterrows():
        chr = row["chr"]
        start = row["start"]
        stop = row["stop"]
        cov = row["cov"]
        length = row["length"]
        if length == 1:
            cord_lst.append([chr, start + 1, cov])
        if length > 1:
            for i in range(length):
                cord_lst.append([chr, start + 1 + i, cov])
    updated_df = pd.DataFrame(cord_lst, columns=["chr", "position", "coverage"])

    return updated_df


def run():
    
    sample_name = os.path.basename(bedtools_output).split(".")[0]
    output_df = read_bedtools_output(bedtools_output)
    df = adding_per_base(output_df)
    df.to_csv(f"{sample_name}.coverage", sep="\t", index=False, header=False)


if __name__ == "__main__":
    run()
