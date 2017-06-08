# Filename: synonymsfetch.py
# Author: Thomas H. Smith 2017

from multiprocessing.dummy import Pool as ThreadPool
import pandas as pd
import xml.dom.minidom
import urllib2
from urllib2 import HTTPError, URLError
from requests.exceptions import ReadTimeout
import pkg_resources, os

def build_syns_df(_df, save_pickle=True, load_pickle=True, add_col_inplace=False):
    """
    Retrieve synonyms associated with a given protein. Synonyms
    are retrieved from Uniprot.org. Option to save the db or
    begin with a db generated from a previous analysis.

    Parameters
    ----------
    _df : pandas.DataFrame
        Input df containing protein IDs

    save_pickle : bool
        Option to save results to pickle file. If True, results will
        be saved in package /data directory as saved_synonyms.pickle

    load_pickle : bool
        If True, will attempt to load a previously generated synonyms
        DataFrame from data/saved_synonyms.pickle. If successful, then
        synonyms will only be retrieved for proteins not contained
        in existing DataFrame, and an appended DataFrame will be returned.

    add_col_inplace : bool
        If True, rather than returning DataFrame containing synonyms, the
        input DataFrame will be returned with a 'Synonyms' column appended.

    Returns
    -------
    df_syns : pandas.DataFrame
        Output df containing synonyms associated with input
        df unique protein IDs

    """
    have_existing=False
    df_in = _df.copy()
    uniprot_ids = list(df_in.Protein.unique())

    if load_pickle==True:
        PICKLE_FILE = pkg_resources.resource_filename(__name__, 'data/saved_synonyms.pickle')
        print 'Attempting to load saved synonyms file...'
        if os.path.isfile(PICKLE_FILE):
            df_syns_loaded = pd.read_pickle(PICKLE_FILE)
            have_existing=True
            print 'Loaded saved synonyms for %d unique proteins' % len(df_syns_loaded)
            known_ids = list(df_syns_loaded.index)
            uniprot_ids = set(uniprot_ids) - set(known_ids)
            uniprot_ids = list(uniprot_ids)
        else:
            print 'Saved synonyms file not found.'

    if len(uniprot_ids) > 0:
        print 'Getting synonyms from Uniprot for %d unique proteins.' % len(uniprot_ids)
        print 'This may take several minutes...'
        pool = ThreadPool(32)
        results = pool.map(_fetchUniprotSynonyms, uniprot_ids)
        pool.close()
        pool.join()
        df_syns = pd.DataFrame(results)
        df_syns.set_index('UniprotID', inplace=True)
        if have_existing==True:
            df_syns = df_syns_loaded.append(df_syns)
        print 'Finished building synonyms table.'
        if save_pickle==True:
            PICKLE_FILE = pkg_resources.resource_filename(__name__, 'data/saved_synonyms.pickle')
            df_syns.to_pickle(PICKLE_FILE)
            print 'Saved synonyms for %d unique proteins as %s' % (len(df_syns), PICKLE_FILE)
        if add_col_inplace:
            df_in['Synonyms'] = df_in['Protein'].apply(lambda x: df_syns.loc[x, 'Synonyms'])
            print "Added 'Synonyms' column to input DataFrame"
            return df_in
        else:
            return df_syns
    else:
        print 'All proteins were found in saved synonyms file'
        if add_col_inplace:
            df_in['Synonyms'] = df_in['Protein'].apply(lambda x: df_syns_loaded.loc[x, 'Synonyms'])
            print "Added 'Synonyms' column to input DataFrame"
            return df_in
        else:
            return df_syns_loaded

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
    syns = set(syns)
    syns = list(syns)
    return {'UniprotID':uniprot_ID, 'Synonyms':syns}

