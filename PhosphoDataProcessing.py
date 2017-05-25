#!/usr/bin/python
# Filename: PhosphoDataProcessing
# Author: Thomas H. Smith 2017

"""
These are a collection of functions I found useful for
processing phosphoproteomics datasets that contain
phosphopeptide sequences with associated intensity values for
each condition, RS score, num phospho sites and unitprot ACC
identifier in a csv or excel format
v0.4
"""
import re
import xml.dom.minidom
import urllib2
from urllib2 import HTTPError, URLError
from operator import itemgetter
from multiprocessing.dummy import Pool as ThreadPool
from tqdm import tqdm_notebook
import numpy as np
import imp
import pandas as pd
from pyshorteners import Shortener
from requests.exceptions import ReadTimeout
pmt = imp.load_source('PubmedTools','/scripts/PubmedTools.py')

# Path to directory containing libraries
LIB_PATH = '/scripts/lib'

# Import required libraries
#DF_SYNS = pd.read_pickle(LIB_PATH+'/prot_genes_syns_dict.pickle')
DF_DOMAINS = pd.read_pickle(LIB_PATH+'/Uniprot_domains_df_20170318.pickle')
DF_SEQS = pd.read_pickle(LIB_PATH+'/phosphosite_ref_seqs.pickle')
DF_SEQS.drop('SPECIES', axis=1, inplace=True)
DF_SEQS.columns = ['Gene', 'Gene_Isoform', 'Uniprot_ID', 'Sequence']
DF_REG = pd.read_excel(LIB_PATH+'/phosphosite_regulatory_sites.xlsx')
DF_REG = DF_REG[DF_REG['ORGANISM'] == 'human'].copy()
DF_REG.drop('Unnamed: 20', axis=1, inplace=True)
STR_COLS = ['ON_PROCESS', 'MOD_RSD', 'DOMAIN', 'ON_PROT_INTERACT', 'PROT_TYPE', 'ON_FUNCTION',
            'ON_OTHER_INTERACT', 'NOTES']
DF_REG[STR_COLS] = DF_REG[STR_COLS].astype(str)
DF_REPORTED = pd.read_csv(LIB_PATH+'/phosphosite_reported_site_dataset', delimiter='\t')
DF_REPORTED = DF_REPORTED[DF_REPORTED['ORGANISM'] == 'human'].copy()
DF_PUBS = pd.read_pickle(LIB_PATH+'/Gene_complete_pub_counts_wLinks_20170318.pickle')

def fix_gene_names(desc, curr_name):
    """
    Fix gene names that Excel has converted to dates
    ie Excel automatically changes SEPT2 to Sep-2

    Parameters
    ----------
    desc : str
        'Description' data, which contains gene names folloed by "GN=".
    curr_name : str
        Current Gene name from dataset, in case Description does not
        contain "GN="

    Returns
    -------
    gene_name : str
        The corrected gene name.

    """
    words = desc.split('GN=')
    if len(words) > 1:
        gene_name = words[1].split()[0]
        return gene_name
    else:
        return curr_name

def get_gene_name(desc):
    """
    Fix gene names that Excel has converted to dates
    ie Excel automatically changes SEPT2 to Sep-2

    Parameters
    ----------
    desc : str
        'Description' data, which contains gene names folloed by "GN=".

    Returns
    -------
    gene_name : str
        The corrected gene name.

    """
    words = desc.split('GN=')
    if len(words) > 1:
        gene_name = words[1].split()[0]
        return gene_name

def build_simple_seq_string(concatamer):
    """

    Simple function to simplify concatamer sequence to string

    Parameters
    ----------
    concatamer_col : str
        Contacamer string, which is in the form [K].MWGRSTLLYSRX.[L]

    Returns
    -------
    seq : str
        Simplified sequence.

    """
    return concatamer.split('.')[1]

def build_syns_df(_df):
    """

    Retrieve synonyms associated with a given protein. Synonyms
    db was retrieved from Uniprot.org

    Parameters
    ----------
    protein : str
        Protein uniprot ACC ID

    Returns
    -------
    syns : str
        Synonyms seperated by commas

    """
    #syns = DF_SYNS.ix[protein].Synonyms
    df_in = _df.copy()
    uniprot_ids = list(df_in.Protein.unique())
    pool = ThreadPool(32)
    results = pool.map(_fetchUniprotSynonyms, uniprot_ids)
    pool.close()
    pool.join()
    df_syns = pd.DataFrame(results)
    df_syns['Synonyms'] = df_syns['Synonyms'].apply(lambda x: _conv_syns(x))
    df_syns.set_index('UniprotID', inplace=True)
    return df_syns

def _conv_syns(syns):
    syns = ','.join([str(x) for x in syns])
    return syns


def _fetchUniprotSynonyms(uniprot_ID):
    url = 'http://www.uniprot.org/uniprot/%s.xml' % uniprot_ID
    try:
        fp = urllib2.urlopen(url)
        doc = xml.dom.minidom.parse(fp)
        fp.close()
    except HTTPError:
        fp.close()
        print 'Caught HTTPError for %s' % uniprot_ID
        return [0]
    syns = []
    protein_elements = doc.getElementsByTagName('protein')
    if len(protein_elements) > 0:
        for prot_elem in protein_elements:
            full_names = prot_elem.getElementsByTagName('fullName')
            short_names = prot_elem.getElementsByTagName('shortName')
            if len(full_names) > 0:
                for node in full_names:
                    syns.append(node.childNodes[0].data)
            if len(short_names) > 0:
                for node in short_names:
                    syns.append(node.childNodes[0].data)
    gene_elements = doc.getElementsByTagName('gene')
    if len(gene_elements) > 0:
        for gene_elem in gene_elements:
            names = gene_elem.getElementsByTagName('name')
            if len(names) > 0:
                for node in names:
                    syns.append(node.childNodes[0].data)
    return {'UniprotID':uniprot_ID, 'Synonyms':syns}


def _populate_seqs_dict(protein):
    global SEQS_DICT
    SEQS_DICT[protein] = _get_full_protein_sequence(protein)

def _build_seqs_db(missing):
    global SEQS_DICT
    SEQS_DICT = {}
    pool = ThreadPool(32)
    results = pool.map(_populate_seqs_dict, missing)
    #close the pool and wait for the work to finish
    pool.close()
    pool.join()

def _get_class1_sites(df_in, row_ix, rs_col_name, threshold,
                      include_uncertain=True):
    # Helper function
    # Given df and RS score column name, extract p-sites
    # above given threshold and return as a list of sites
    # in format [S1, T4, y5]
    sites = []
    mysplit = df_in.iloc[row_ix][rs_col_name].split(';')
    for item in mysplit:
        words = item.split('(Phospho):')
        residue = words[0].strip()
        prob = float(words[1])
        if prob >= threshold:
            sites.append(residue)
        else:
            if include_uncertain:
                unsure_res = residue[0].lower() + ''.join(residue[1:])
                sites.append(unsure_res)
    return sites

def _webfetch_uniprot_seq(protein):
    # Helper function
    # Retrieve reference sequence for a protein from uniprot
    # given the uniprot accession ID, return 0 if not found

    # Construct URL string pointing to uniprot fasta file
    #global SEQS_DICT
    url = 'http://www.uniprot.org/uniprot/%s.fasta' % protein
    try: # Catch errors arising from invalid URL
        data = urllib2.urlopen(url)
        seq = ''
        for line in data.readlines()[1:]:
            seq = seq + line.rstrip()
        #SEQS_DICT[protein] = seq
        return seq
    except HTTPError, URLError:
        print '%s: caught HTTPError' % protein
        return ''

def _get_full_protein_sequence(uniprot_acc):
    global SEQS_DICT
    if DF_SEQS.Uniprot_ID.str.contains(uniprot_acc).any():
        uniprot_seq = DF_SEQS[DF_SEQS.Uniprot_ID == uniprot_acc].Sequence.values[0]
        return uniprot_seq
    else:
        if SEQS_DICT.has_key(uniprot_acc):
            return SEQS_DICT[uniprot_acc]
        else: return _webfetch_uniprot_seq(uniprot_acc)

# Helper function - search full-length protein seq for subseq and index phosph-sites in this context
def _identify_site_in_seq(protein, subseq, site):
    uniprot_seq = _get_full_protein_sequence(protein)
    if len(uniprot_seq) < 5:
        return 0
    seq_pos_index = uniprot_seq.find(str.upper(subseq))
    site_pos = int(site[1:]) + seq_pos_index
    full_annot = '%s%d' % (site[0], site_pos)
    return full_annot

def _populate_site_annotation_cols(_df, site):
    df_in = _df.copy()
    j = len(df_in) -1
    protein = df_in.loc[j]['Protein']
    subseq = df_in.loc[j]['Sequence']
    df_in.loc[j, 'LocalSite'] = site
    annot_site = _identify_site_in_seq(protein, subseq, site)
    if annot_site == 0:
        df_in.loc[j, 'AnnotatedSite'] = 'na'
        df_in.loc[j, 'Residue'] = 'na'
        df_in.loc[j, 'Position'] = 0
        return df_in
    else:
        residue = annot_site[0]
        position = annot_site[1:]
        df_in.loc[j, 'AnnotatedSite'] = annot_site
        df_in.loc[j, 'Residue'] = residue
        df_in.loc[j, 'Position'] = position
        return df_in


def identify_phosphosites(_df, phos_rs_col, threshold):
    """

    Annotate phosphosites over a specified threshold phosRS value
    in the context of corresponding full-length protein sequences.
    Rows with multiple sites are split into multiple rows, one for
    each phosphosite. Sites under threshold are dropped.

    Parameters
    ----------
    _df : pandas.DataFrame
        DataFrame containing data
    phos_rs_col : str
        Name of column containing phosRS scores
    threshold : int
        Minimum phosRS score to include

    Returns
    -------
    df_new : pandas.DataFrame
        New DataFrame containing additional rows expanded from peptides
        containing multiple phosphosites, and additional columns containing
        phosphosite annotation data.

    """

    df_in = _df.copy()
    df_new = pd.DataFrame()
    pbar1 = tqdm_notebook(range(len(df_in)), total=len(df_in))
    for i in pbar1:
        protein = df_in.iloc[i]['Protein']
        subseq = df_in.iloc[i]['Sequence']
        sites = _get_class1_sites(df_in, i, phos_rs_col, threshold, include_uncertain=False)
        if len(sites) > 0:
            for site in sites:
                df_new = df_new.append(df_in.iloc[i], ignore_index=True)
                j = len(df_new) - 1
                df_new.loc[j, 'LocalSite'] = site
                annot_site = _identify_site_in_seq(protein, subseq, site)
                residue = annot_site[0]
                position = annot_site[1:]
                df_new.loc[j, 'AnnotatedSite'] = annot_site
                df_new.loc[j, 'Residue'] = residue
                df_new.loc[j, 'Position'] = position
    return df_new

def identify_phosphosites_one_run(_df, phos_rs_col1, val_cols1, threshold, missing_values=np.nan):
    """

    Annotate phosphosites over a specified threshold phosRS value
    in the context of corresponding full-length protein sequences.
    Use data acquired from one run. Rows with multiple sites are split
    into multiple rows, one for each phosphosite. Sites under threshold
    are dropped.

    Parameters
    ----------
    _df : pandas.DataFrame
        DataFrame containing data
    phos_rs_col1 : str
        Name of column containing phosRS scores from run 1
    val_cols1 : list of str
        Names of columns containing data values from run 1
    threshold : int
        Minimum phosRS score to include
    missing_values : int, str, or NaN
        Value to impute for missing data

    Returns
    -------
    df_new : pandas.DataFrame
        New DataFrame containing additional rows expanded from peptides
        containing multiple phosphosites, and additional columns containing
        phosphosite annotation data.

    """

    global SEQS_DICT
    num_rows_dropped = 0
    df_in = _df.copy()
    df_new = pd.DataFrame()

    # Build dict of missing sequences using thread pool
    # this significantly decreases processing time
    SEQS_DICT = {}
    known = list(DF_SEQS.Uniprot_ID)
    targets = df_in.Protein.unique()
    missing = set(targets) - set(known)
    missing = list(missing)
    _build_seqs_db(missing)

    # Annotation
    pbar1 = tqdm_notebook(range(len(df_in)), total=len(df_in))
    for i in pbar1:
        protein = df_in.iloc[i]['Protein']
        subseq = df_in.iloc[i]['Sequence']
        sites1 = _get_class1_sites(df_in, i, phos_rs_col1, threshold, include_uncertain=False)
        if len(sites1) > 0: # If any sites met threshold for run1
            for site in sites1:
                df_new = df_new.append(df_in.iloc[i], ignore_index=True)
                df_new = _populate_site_annotation_cols(df_new, site)
        else: # if sites1 does not contain any sites
            num_rows_dropped += 1
    print 'Dropped '+str(num_rows_dropped)+' rows that failed to meet phosRS'\
    'threshold for at least one run'
    return df_new


def identify_phosphosites_two_runs(_df, phos_rs_col1, val_cols1, phos_rs_col2,
                                   val_cols2, threshold, missing_values=0):
    """

    Annotate phosphosites over a specified threshold phosRS value
    in the context of corresponding full-length protein sequences.
    Use data acquired from two runs, where each run has a seperate
    phosRS score. Rows with multiple sites are split into multiple
    rows, one for each phosphosite. Sites under threshold are dropped.
    For sites under threshold for only one run, data columns are only
    conserved for the run meeting the threshold, and imputed with
    zeroes or NaN for the other run.

    Parameters
    ----------
    _df : pandas.DataFrame
        DataFrame containing data
    phos_rs_col1 : str
        Name of column containing phosRS scores from run 1
    val_cols1 : list of str
        Names of columns containing data values from run 1
    phos_rs_col2 : str
        Name of column containing phosRS scores from run 2
    val_cols2 : list of str
        Names of columns containing data values from run 2
    threshold : int
        Minimum phosRS score to include
    missing_values : int, str, or NaN
        Value to impute for missing data

    Returns
    -------
    df_new : pandas.DataFrame
        New DataFrame containing additional rows expanded from peptides
        containing multiple phosphosites, and additional columns containing
        phosphosite annotation data.

    """

    global SEQS_DICT
    num_rows_dropped = 0
    df_in = _df.copy()
    df_new = pd.DataFrame()

    # Build dict of missing sequences using thread pool
    # this significantly decreases processing time
    SEQS_DICT = {}
    known = list(DF_SEQS.Uniprot_ID)
    targets = df_in.Protein.unique()
    missing = set(targets) - set(known)
    missing = list(missing)
    _build_seqs_db(missing)

    # Annotation
    pbar1 = tqdm_notebook(range(len(df_in)), total=len(df_in))
    for i in pbar1:
        protein = df_in.iloc[i]['Protein']
        subseq = df_in.iloc[i]['Sequence']
        if df_in.iloc[i][phos_rs_col1] != 'na':
            sites1 = _get_class1_sites(df_in, i, phos_rs_col1, threshold, include_uncertain=False)
        else: sites1 = []
        if df_in.iloc[i][phos_rs_col2] != 'na':
            sites2 = _get_class1_sites(df_in, i, phos_rs_col2, threshold, include_uncertain=False)
        else: sites2 = []

        if len(sites1) > 0: # If any sites met threshold for run1
            if len(sites2) > 0: # If both runs have at least one site that met threshold
                for site in sites1:
                    if site in sites2: # For a site in both runs, just copy over the data
                        df_new = df_new.append(df_in.iloc[i], ignore_index=True)
                        sites2.remove(site) # Remove site from sites2 to avoid counting twice
                        df_new = _populate_site_annotation_cols(df_new, site)
                    else: # for a site only in run1, copy over data and set val_cols2 = NaN
                        df_new = df_new.append(df_in.iloc[i], ignore_index=True)
                        j = len(df_new) - 1
                        df_new.loc[j, val_cols2] = missing_values
                        df_new = _populate_site_annotation_cols(df_new, site)
                for site in sites2: # finish up by dealing with whatever is left in sites2 list
                    df_new = df_new.append(df_in.iloc[i], ignore_index=True)
                    j = len(df_new) - 1
                    df_new.loc[j, val_cols1] = missing_values
                    df_new = _populate_site_annotation_cols(df_new, site)
            else: # if sites1 has at least one site,
                  # but sites2 is empty, copy over data, set vals2=NaN
                for site in sites1:
                    df_new = df_new.append(df_in.iloc[i], ignore_index=True)
                    j = len(df_new) - 1
                    df_new.loc[j, val_cols2] = missing_values
                    df_new = _populate_site_annotation_cols(df_new, site)
        else: # if sites1 does not contain any sites
            if len(sites2) > 0: # if sites2 contains at least one site
                for site in sites2:
                    df_new = df_new.append(df_in.iloc[i], ignore_index=True)
                    j = len(df_new) - 1
                    df_new.loc[j, val_cols1] = missing_values
                    df_new = _populate_site_annotation_cols(df_new, site)
            else: num_rows_dropped += 1
    print 'Dropped '+str(num_rows_dropped)+' rows that failed to meet phosRS'\
    'threshold for at least one run'
    return df_new


def phosphosite_funct_annot(_df):
    """

    Annotate functional phosphosites using
    PhosphoSitePlus database as reference

    Parameters
    ----------
    _df : pandas.DataFrame
        DataFrame containing data

    Returns
    -------
    df_new : pandas.DataFrame
        New DataFrame with additional columns
        containing functional site data.

    """

    df_new = _df.copy()
    pbar1 = tqdm_notebook(range(len(df_new)), total=len(df_new))
    for i in pbar1:
        protein = df_new.iloc[i]['Protein']
        gene = df_new.iloc[i]['Gene']
        site = df_new.iloc[i]['AnnotatedSite'] + '-p'
        df_match = DF_REG[(DF_REG['ACC_ID'] == protein) & (DF_REG['MOD_RSD'] == site)]
        if len(df_match > 0):
            df_new.loc[i, 'Functional_site'] = '+'
            df_new.loc[i, 'DOMAIN'] = str(df_match['DOMAIN'].values[0])
            df_new.loc[i, 'ON_FUNCTION'] = str(df_match['ON_FUNCTION'].values[0])
            df_new.loc[i, 'ON_PROCESS'] = str(df_match['ON_PROCESS'].values[0])
            df_new.loc[i, 'ON_PROT_INTERACT'] = str(df_match['ON_PROT_INTERACT'].values[0])
            df_new.loc[i, 'ON_OTHER_INTERACT'] = str(df_match['ON_OTHER_INTERACT'].values[0])
            df_new.loc[i, 'NOTES'] = str(df_match['NOTES'].values[0])
        else:
            df_new.loc[i, 'Functional_site'] = '-'
            df_new.loc[i, 'DOMAIN'] = 0
            df_new.loc[i, 'ON_FUNCTION'] = 0
            df_new.loc[i, 'ON_PROCESS'] = 0
            df_new.loc[i, 'ON_PROT_INTERACT'] = 0
            df_new.loc[i, 'ON_OTHER_INTERACT'] = 0
            df_new.loc[i, 'NOTES'] = 0
        if len(DF_REPORTED[(DF_REPORTED['ACC_ID'] == protein) &
                           (DF_REPORTED['MOD_RSD'] == site)]) > 0:
            df_new.loc[i, 'Known_site'] = '+'
        else:
            if len(DF_REPORTED[(DF_REPORTED['GENE'] == gene) &
                               (DF_REPORTED['MOD_RSD'] == site)]) > 0:
                df_new.loc[i, 'Known_site'] = '+'
            else:
                df_new.loc[i, 'Known_site'] = '-'
    return df_new


def add_uniprot_psites_string_col(_df, protein_col, annot_seq_col, site_sep='_'):
    """
    Assemble a string conatining the uniprot accession ID
    with phosphosites appended ie A014012_S123_T130

    Parameters
    ----------
    _df : pandas.DataFrame
        DataFrame from which to pull data.
    protein_col : str
        Name of df_in column that contains uniprot ACC IDs.
    annot_seq_col : str
        Name of df_in column that contains phosphorylation
        sites in form S312_t314_y315.
    site_set : str, optional
        Character(s) that seperate site annotations for
        peptides containing multiple p-sites ie S312_t314_y315
        Default character is "_" which is generated by the
        annotate_phosphosites function.

    Returns
    -------
    df : pandas.DataFrame
        A copy of the original input DataFrame but with the
        addition of an extra column named 'Uniprot_pSites"

    Examples
    --------
    >>> from AnalysisTools import add_uniprot_psites_column
    >>> df = add_uniprot_psites_column(df, 'Protein', 'Full_annot')

    """

    df_in = _df.copy()
    pbar1 = tqdm_notebook(list(df_in.index), total=len(df_in))
    for i in pbar1:
        protein_str = df_in.loc[i, protein_col]
        annot_str = df_in.loc[i, annot_seq_col]
        uniprot_annot = '%s_%s' % (protein_str, annot_str)
        df_in.loc[i, 'Uniprot_pSites'] = uniprot_annot
    return df_in


def count_sty_sites(_df):
    """

    Simple function to count number of
    S/T/Y phosphosites in dataset.

    Parameters
    ----------
    _df : pandas.DataFrame
        DataFrame containing data

    Returns
    -------
    counts_dict : dict
        Dict containing value counts for S/T/Y
        phosphosites.

    """

    counts_dict = {'S':0, 'T':0, 'Y':0}
    indices = range(len(df))
    for index in indices:
        sites = list(df.iloc[index].Full_annot)
        for site in sites:
            mod = site[0]
            counts_dict[mod] = counts_dict[mod] + 1
    return counts_dict


def get_flanking_sites(protein, gene, site, webfetch_unknown_seqs=True):
    """

    Given a site, return the n flanking residues before
    and after the site within the full protein sequence.

    Parameters
    ----------
    protein : str
        Uniprot ACC ID of protein of interest
    gene : str
        Gene symbol ID
    site : str
        Phosphosite of interest
    webfetch_unknown_seqs : boolean
        Whether to attempt to retrieve unknown protein
        sequences from the web

    Returns
    -------
    motif : str
        Subsequence of full protein sequence corresponding
        to region n residues before and after specified
        phosphosite.

    """

    site_mod = site.lower()
    site_loc = int(site[1:])-1
    if DF_SEQS.Uniprot_ID.isin([protein]).any():
        seq = DF_SEQS[DF_SEQS.Uniprot_ID == protein].iloc[0]['Sequence']
    else:
        if DF_SEQS.Gene.isin([gene]).any():
            seq = DF_SEQS[DF_SEQS.Gene == gene].iloc[0]['Sequence']
        else:
            if webfetch_unknown_seqs:
                seq = _webfetch_uniprot_seq(protein)
                if seq == '':
                    print 'no refseq for %s and webfetch failed' % protein
                    return ''
            else:
                return ''
    if site_mod.upper() == seq[site_loc]:
        if (site_loc - num_residues >= 0) & (site_loc+num_residues < len(seq)):
            start = (site_loc - num_residues)
            stop = (site_loc + num_residues + 1)
            motif = seq[start:site_loc] + site_mod + seq[site_loc+1:stop]
            return motif
    else:
        if webfetch_unknown_seqs:
            seq = _webfetch_uniprot_seq(protein)
            if site_mod.upper() == seq[site_loc]:
                if (site_loc - num_residues >= 0) & (site_loc+num_residues < len(seq)):
                    start = (site_loc - num_residues)
                    stop = (site_loc + num_residues + 1)
                    motif = seq[start:site_loc] + site_mod + seq[site_loc+1:stop]
                    return motif
            else:
                return ''
        else:
            print 'no refseq for %s and webfetch disabled' % protein
            return ''


def _fetch_uniprot_func_domains(uniprot_id, sites, wiggle):
    """

    Helper function to check for sites that fall within
    uniprot functional domains for proteins not already
    in functional domain local db. Pulls domains from
    Uniprot.org.

    Parameters
    ----------
    protein : str
        Uniprot ACC ID of protein of interest
    sites : list of str
        List of phosphosites with each item in format
        'X###' ie 'R420'
    wiggle : int
        For identifying sites that occur near the start
        or stop of a functional domain

    Returns
    -------
    matches : list of str
        List where each item is a string corresponding to
        an input site that falls within a functional domain.
        Formatted as (<name of domain>, <start position>,
        <site>, <stop position>)

    """

    matches = []
    ignore_list = ['chain', 'splice variant', 'sequence conflict']
    url = 'http://www.uniprot.org/uniprot/%s.xml' % uniprot_id
    try:
        _fp = urllib2.urlopen(url)
        doc = xml.dom.minidom.parse(_fp)
        _fp.close()
        attributes_dict = {}
        features = doc.getElementsByTagName('feature')
    except HTTPError:
        return 0
    #print 'Checking %d features...' % (len(features))
    for feature in features:
        feat_name = feature.getAttributeNode('type').value
        if str(feat_name) in ignore_list:
            continue
        try:
            feat_name = feature.getAttributeNode('description').value
        except AttributeError:
            feat_desc = ''
        loc = feature.getElementsByTagName('location')[0]
        # Skip anything with no range (ie modification sites)
        if len(loc.getElementsByTagName('begin')) > 0:
            begin = loc.getElementsByTagName('begin')[0]
            end = loc.getElementsByTagName('end')[0]
            start_pos = int(begin.attributes['position'].value)
            stop_pos = int(end.attributes['position'].value)
        else:
            continue
        for site in sites:
            site_pos = int(site[1:])
            if (site_pos >= start_pos-wiggle) & (site_pos <= stop_pos+wiggle):
                match = '%s:[%d-(%s)-%d]' % (feat_name, start_pos, site, stop_pos)
                matches.append(match)
    return matches


def annotate_sites_in_domains(protein, sites, wiggle):
    """

    Annotate phosphosites that occur within functional
    domains as defined by Uniprot.org.

    Parameters
    ----------
    protein : str
        Uniprot ACC ID of protein of interest
    sites : list of str
        List of phosphosites with each item in format
        'X###' ie 'R420'
    wiggle : int
        For identifying sites that occur near the start
        or stop of a functional domain

    Returns
    -------
    matches : list of str
        List where each item is a string corresponding to
        an input site that falls within a functional domain.
        Formatted as (<name of domain>, <start position>,
        <site>, <stop position>)

    """

    matches = []
    try:
        func_domains = DF_DOMAINS.ix[protein][0]
    except KeyError:
        print 'couldnt find %s..trying webfetch' % protein
        matches = _fetch_uniprot_functional_domains(protein, sites, wiggle)
        if matches != 0:
            return matches
        else:
            print 'webfetch FAILED for %s' % protein
            return [0]
    if len(func_domains) == 0:
        return [0]
    for domain in func_domains:
        name = domain.split(':')[0:-1]
        locs = domain.split(':')[-1].replace('[', '').replace(']', '').split('-')
        start_pos = int(locs[0])
        stop_pos = int(locs[1])
        for site in sites:
            site_pos = int(site[1:])
            if site_pos in range(start_pos, stop_pos):
                match = '%s:[%d-(%s)-%d]' % (name[0], start_pos, site, stop_pos)
                matches.append(match)
    return matches


def pubmed_pub_count(term1_list, term2_list, get_pmids_list=False):
    base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term='
    base_link_url = 'https://www.ncbi.nlm.nih.gov/pubmed/?term='
    term1_string = '(%s)' % term1_list[0]
    if len(term1_list) > 1:
        for term in term1_list[1:]:
            term1_string = '%s OR (%s)' % (term1_string, term)
    term2_string = '(%s)' % term2_list[0]
    if len(term2_list) > 1:
        for term in term2_list[1:]:
            term2_string = '%s OR (%s)' % (term2_string, term)

    search_string = '(%s) AND (%s)' % (term1_string, term2_string)
    search_string = search_string.replace(' ', '+')
    url = '%s%s&retmax=100000' % (base_url, search_string)
    link_url = base_link_url + search_string
    try:
        fp = urllib2.urlopen(url)
        doc = xml.dom.minidom.parse(fp)
    except URLError, HTTPError:
        print 'Caught exception for %s' % term1_string
        return ['na', 0]

    res = doc.getElementsByTagName('eSearchResult')
    c = res[0].getElementsByTagName('Count')
    count = int(c[0].childNodes[0].data)
    print search_string
    print count
    print link_url
    if get_pmids_list == False:
        return [count, link_url]
    else:
        idList = doc.getElementsByTagName('IdList')
        id_nodes = idList[0].getElementsByTagName('Id')
        ids_list = []
        for pmid in id_nodes:
            ids_list.append(str(pmid.childNodes[0].data))
        return ids_list



    """
    Add number of publications associated with gene names
    and synonyms.

    Parameters
    ----------
    _df : pandas.DataFrame
        DataFrame containing dataset.
    as_links : boolean, optional
        Whether to format data as a string that will be interpreted
        as a hyperlink when exported to Excel

    Returns
    -------
    df_in : pandas.DataFrame
        DataFrame containing original dataset with additional
        columns containing publication numbers

    """

def _get_pub_count(term1_list, term2_list, title_abs_only=True,
                  return_IDs_list=False, return_link=False):
    if title_abs_only:
        suffix = '[TIAB]'
    else:
        suffix = ''

    term1_string = '(%s%s)' % (term1_list[0], suffix)
    if len(term1_list) > 1:
        for term in term1_list[1:]:
            term1_string = '%s OR (%s%s)' % (term1_string, term, suffix)
    if term2_list[0] == '':
        search_string = '%s' % term1_string
    else:
        term2_string = '(%s%s)' % (term2_list[0], suffix)
        if len(term2_list) > 1:
            for term in term2_list[1:]:
                term2_string = '%s OR (%s%s)' % (term2_string, term, suffix)
        search_string = '(%s) AND (%s)' % (term1_string, term2_string)
    search_string = search_string.replace(' ', '+')
    return _get_esearch_data(search_string, return_IDs_list, return_link)



def _build_args_tuples(df_in, search_terms):
    # Build list of gene-specific search terms for each row then pass to get_pub_count function
    # save results either in place in original df or in a copy to be returned

    #build new df to serve as db for fetched counts data
    unique_genes = list(df_in.Gene.unique())
    print 'Building list of search strings for %d unique genes...' % len(unique_genes)
    args_tuples = []
    for gene in unique_genes:
        row = df_in[df_in.Gene == gene].reset_index().loc[0]
        desc = str(row.Description)
        try:
            term1_list = [str(x) for x in row.Synonyms.split(',')]
            if len(term1_list) > 0:
                term1_list.extend([gene, desc])
                args_tuples.append( (gene, term1_list, search_terms))
        except AttributeError:
            print 'Error getting synonyms for gene: %s' % gene
            args_tuples.append( (gene, [''], search_terms))
    return args_tuples

def _get_esearch_data(search_string, return_IDs_list=False, return_link=True):
    URL_PREFIX = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term='
    LINK_URL_PREFIX = 'https://www.ncbi.nlm.nih.gov/pubmed/?term='
    if return_IDs_list == False:
        URL_SUFFIX = '&rettype=Count&retmode=json'
    else:
        URL_SUFFIX = '&retmax=100000'
    search_string = search_string.replace(' ', '+')
    url = '%s%s%s' % (URL_PREFIX, search_string, URL_SUFFIX)

    if return_IDs_list == False:
        s = urllib2.urlopen(url).read()
        i = s.find('count')
        count = int(s[i:].split(':')[1].split('"')[1])
        if return_link:
            link_url = '=HYPERLINK("%s%s",%d)' % (LINK_URL_PREFIX, search_string, count)
#            if len(link_url) > 255:
#                url = '%s%s' % (LINK_URL_PREFIX, search_string)
#                api_key = 'AIzaSyBRdk9gGnyvMAe_QNiwQfRDW_1jILh_3D8'
#                shortener = Shortener('Google', api_key=api_key)
            return (count, link_url)
        else:
            return count


def shorten_large_links(link):
    url = link.split('"')[1]
    count = link.split('"')[2].split(',')[1].split(')')[0]
    api_key = 'AIzaSyBRdk9gGnyvMAe_QNiwQfRDW_1jILh_3D8'
    shortener = Shortener('Google', api_key=api_key)
    try:
        short_url = str(shortener.short(url))
        new_link = '=HYPERLINK("%s", %s)' % (short_url, count)
        return new_link
    except ReadTimeout:
        print 'Caught readtimeout error for %s' % link
        return link


# called by each thread
def _thread_helper_func(args_tuple):
    gene, term1_list, search_terms = args_tuple
    count, link = _get_pub_count(term1_list, search_terms, return_link=True)
    return {'Gene': gene, 'Count_link':link}


def build_pubcount_df(df_in, search_terms, col_name):
    args_tuples = _build_args_tuples(df_in, search_terms)
    print 'Getting data from PubMed...'
    pool = ThreadPool(32)
    results = pool.map(_thread_helper_func, args_tuples)
    df_db = pd.DataFrame(results)
    df_db.set_index('Gene', inplace=True)
    df_db[col_name] = df_db['Count_link']
    df_db.drop(['Count_link'], inplace=True, axis=1)
    #df_db[col_name].apply(lambda x: _shorten_large_links(x))
    return df_db


