# Filename: pubfetch.py
# Author: Thomas H. Smith 2017

import pandas as pd
import os
from numpy import mean

def import_reported_sites(psites_file, save_pickle=True):
    """
    Convert tsv files downloaded from PhosphoSite.org database to
    Pandas DataFrames, which will be used to retrieve info from
    later in the analysis.  If save_pickle option is enabled, the
    DataFrames will be saved in their pre-processed form for later
    use.

    Parameters
    ----------
    ps_file_path : str
        Path to unzipped and unaltered tsv file downloaded from
        PhosphoSitePlus.
    save_pickle : bool
        Option to save the generated DataFrame as a .pickle for
        easy use in later analyses.
    pickle_name : str
        Name of pickle file to be saved, if save_pickle is chosen

    Returns
    -------
    df_db : pandas.DataFrame
        DataFrame containing all information from the original tsv
        file in a format more convenient for data access in later
        steps of the analysis workflow.
    """

    f = open(psites_file, 'r')
    L = f.readlines()
    f.close()
    row_lens = [len(line.split('\t')) for line in L]
    if (mean(row_lens) % 10) >= 0.5:
        ncols = int(mean(row_lens))+1
    else: ncols = int(mean(row_lens))
    temp_fname = '%s_temp' % psites_file
    f = open(temp_fname, 'w')
    for line in L:
        if len(line.split('\t')) == ncols:
            f.write(line)
    f.close()
    df = pd.read_csv(temp_fname, delimiter='\t')
    os.remove(temp_fname)
    # Drop empty column that is included sometimes
    for col in df.columns:
        if df[col].isnull().all():
            df.drop(col, axis=1, inplace=True)
    df = df[df['ORGANISM'] == 'human'].copy()
    df.reset_index(inplace=True, drop=True)
    if save_pickle:
        pickle_fname = '%s.pickle' % psites_file.split('.')[0]
        df.to_pickle(pickle_fname)
    return df

def import_fasta_seqs(fasta_file, save_pickle=True):
    header = True
    firstEntry = True
    f = open(fasta_file)
    seqs_dict = {}
    results = []
    for line in f.readlines():
        if line[0] == '>':
            if firstEntry==False:
                seqs_dict[entry] = seq
            firstEntry, header = False, False
            words = line.split('|')
            species = words[2].rstrip()
            acc_id = words[3].rstrip()
            gene = words[0].split('GN:')[1].rstrip()
            entry = '%s_%s_%s' % (acc_id, species, gene)
            seq = ''
        else:
            if header == False:
                seq = seq+line.rstrip()
    df = pd.DataFrame(data=seqs_dict.values(), index=seqs_dict.keys())
    df.reset_index(inplace=True)
    df.columns = ['entry', 'Sequence']
    df['species'] = df['entry'].apply(lambda x: x.split('_')[1])
    df['Uniprot_ID'] = df['entry'].apply(lambda x: x.split('_')[0])
    df['Gene'] = df['entry'].apply(lambda x: x.split('_')[2])
    df = df[df['species'] == 'human'].copy()
    df.reset_index(inplace=True, drop=True)
    df.drop(['species', 'entry'], inplace=True, axis=1)
    if save_pickle:
        pickle_fname = '%s.pickle' % fasta_file.split('.')[0]
        df.to_pickle(pickle_fname)
    return df

