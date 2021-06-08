import pandas as pd
import sys, os, glob

# This will glob all the region and threshold output from mosdepth
region_glob = glob.glob('*.regions.bed.gz')
threshold_glob = glob.glob('*.thresholds.bed.gz')
summary_glob = glob.glob('*.mosdepth.summary.txt')
cov_threshold = sys.argv[1]

def read_region(region_lst, summary_lst):
    df_lst_for_avg = []
    df_lst_for_all_sample = []
    for r in region_lst:
        sample_name = r.split('/')[-1].split('.')[0]
        df = pd.read_csv( r, compression='gzip', sep='\t' , header=None, names=['chr', 'start', 'stop', 'target', 'average'])
        tdf = df[['target', 'average']]
        df_lst_for_avg.append(tdf)
        sdf = df[['target', 'average']]
        sdf.columns = ['target', f'{sample_name}']
        sdf.set_index('target', inplace=True)
        df_lst_for_all_sample.append(sdf)

    avg_mdf = pd.concat(df_lst_for_avg, axis=0)
    # print(mdf)
    gdf = avg_mdf.groupby('target').mean()
    gdf.reset_index(inplace=True)
    # print(gdf)

    sample_mdf = pd.concat(df_lst_for_all_sample, axis=1)

    sdf_lst = []
    for s in summary_lst:
        sample_name = s.split('/')[-1].split('.')[0]
        df = pd.read_csv( s, sep='\t' , header=0)
        df['sample_name'] = sample_name
        df = df[df['chrom'] == "total_region"]
        df = df[['sample_name', 'length', 'mean', 'min', 'max']]
        # print(df)
        sdf_lst.append(df)
    smdf = pd.concat(sdf_lst, axis=0)
    # print(smdf)



    return sample_mdf, gdf, smdf

def read_threshold(threshold_lst):
    """This reads and compiles .threshold.bed.gz file from the mosdepth output"""

    df_lst = []
    sdf_lst = []
    for t in threshold_lst:
        sample_name = t.split('/')[-1].split('.')[0]
        df = pd.read_csv( t, compression='gzip', sep='\t' , header=0)
        df['sample_name'] = sample_name
        df['region_size'] = df['end'] - df['start']

        sdf = df.groupby('sample_name').sum()
        sdf['%Over10X'] = sdf['10X'] / sdf['region_size']
        sdf['%Over20X'] = sdf['20X'] / sdf['region_size']
        sdf['%Over50X'] = sdf['50X'] / sdf['region_size']
        sdf['%Over100X'] = sdf['100X'] / sdf['region_size']
        sdf['%Over500X'] = sdf['500X'] / sdf['region_size']
        sdf['%Over1000X'] = sdf['1000X'] / sdf['region_size']
        sdf = sdf[['%Over10X', '%Over20X', '%Over50X', '%Over100X', '%Over500X','%Over1000X']]
        # print(sdf)
        sdf_lst.append(sdf)

        df['%Over10X'] = df['10X'] / df['region_size']
        df['%Over20X'] = df['20X'] / df['region_size']
        df['%Over50X'] = df['50X'] / df['region_size']
        df['%Over100X'] = df['100X'] / df['region_size']
        df['%Over500X'] = df['500X'] / df['region_size']
        df['%Over1000X'] = df['1000X'] / df['region_size']
        
        df = df[['region', '%Over10X', '%Over20X', '%Over50X', '%Over100X', '%Over500X','%Over1000X']]
        # print(df.head())

        df_lst.append(df)

        # To get to one sample per row
        # des_df = df.describe().reset_index()
        # des_df = des_df[des_df['index'] == 'mean']
        # des_df['sample'] = sample
        # des_df = des_df[['sample', '10Xp', '20Xp', '50Xp', '100Xp', '150Xp']]
        # df_lst.append(des_df)

    mdf = pd.concat(df_lst, axis=0)
    gdf = mdf.groupby('region').mean()
    gdf.reset_index(inplace=True)
    # print(gdf)

    msdf = pd.concat(sdf_lst, axis=0)
    msdf.reset_index(inplace=True)
    # print(msdf)

    return gdf, msdf

def read_lst(txt):
    lst = []
    with open(txt) as f:
        for line in f:
            lst.append(line.strip())
    return lst 
    
def compile_exonStat(run_lst):
    df_lst = []
    for run in run_lst:
        exonStat_glob = glob.glob(f'{run}/Aligned_Panel_9000_SGE/Project_9000_HGMDNext/Sample_RD_*/RD_*.coverage.exonStat')
        exonStat_glob = [g for g in exonStat_glob if 'NEG' not in g]
        # print(len(exonStat_glob))
        for e in exonStat_glob:
            df = pd.read_csv( e, skiprows=5)
            print(df.head())
            df_lst.append(df)
    mdf = pd.concat(df_lst, axis=0)
    gdf = mdf.groupby(['Gene','Isoform','Coding exon(CDS)']).mean()
    print(gdf)
    gdf.reset_index(inplace=True)
    gdf.to_excel('exonStat.xlsx', index=False)

def run():
    print('starting to compile')
    writer = pd.ExcelWriter(f'PerTargetMeanCov.xlsx', engine='xlsxwriter')
    mdf, gdf, smdf = read_region(region_glob, summary_glob)
    cols = [ c for c in mdf.columns if 'neg' not in c.lower() ]
    no_neg_mdf = mdf[cols]
    # df < 20 will turn any cells => 20 to NaN, and dropna.all will remove the row if it is all NaN( all > 20), previous cov threshold was set to 20
    print(f'filtering targets below {cov_threshold}X')
    low_cov_df = no_neg_mdf[no_neg_mdf < int(cov_threshold)].dropna(how='all')
    low_cov_df.reset_index(inplace=True)
    low_cov_df_target_names = low_cov_df[['target']]
    below_cov_threshold_df = pd.merge(no_neg_mdf, low_cov_df_target_names, on='target', how='inner') 
    mdf.to_excel(writer, sheet_name = 'All_Samples' ,index = True)
    below_cov_threshold_df.to_excel(writer, sheet_name = f'TargetsBelow_{cov_threshold}X' ,index = False)

    gdf.to_excel(writer, sheet_name = 'AvgAcrossSamples' ,index = False)
    smdf.to_excel(writer, sheet_name = 'SampleLevel' ,index = False)
    writer.save()
    print('done with per target coverage')


    writer = pd.ExcelWriter(f'PerTargetThresholdCov.xlsx', engine='xlsxwriter')
    tgdf, tmsdf = read_threshold(threshold_glob)
    tgdf.to_excel(writer, sheet_name = 'AvgAcrossSamples' ,index = False)
    tmsdf.to_excel(writer, sheet_name = 'SampleLevel' ,index = False)
    writer.save()
    print('done with threshold coverage')
    
    
if __name__=="__main__":
    run()