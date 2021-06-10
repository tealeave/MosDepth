import pandas as pd
import sys, os, glob

bedtools_output = sys.argv[1]
# bedtools_output = '/mnt/dlin/share/REPO/common_scripts_dl/MosDepth/210528_NB501613_0394_AHK5WFBGXH/RD_21_127030_S1A.bedtool.out'
def read_bedtools_output(pb):

    df = pd.read_csv(pb, sep='\t', names=['chr', 'start', 'stop', 'cov'])
    df['length'] = df['stop'] - df['start']
    cleaned_df1 = df[((df['chr'] != 'artifact_seq') & ( df['length']<=300)) 
        # & (df['cov']<20)
        # & (df['length'] > 1)
        ]
    
    print(cleaned_df1)
    return cleaned_df1

def adding_per_base(df):
    cord_lst = []
    for i, row in df.iterrows():
        chr = row['chr']
        start = row['start']
        stop = row['stop']
        cov = row['cov']
        length = row['length']
        if length == 1:
            cord_lst.append([chr, start+1, cov])
        if length > 1:
            for i in range(length):
                cord_lst.append([chr, start+1+i, cov])
    updated_df = pd.DataFrame(cord_lst,columns=['chr', 'position', 'coverage'])
    return updated_df

def run():
    sample_name = os.path.basename(bedtools_output).split('.')[0]
    output_df = read_bedtools_output(bedtools_output)

    df = adding_per_base(output_df)
    # print(df)
    df.to_csv(f'{sample_name}.coverage', sep='\t', index=False, header=False)

if __name__=="__main__":
    run()