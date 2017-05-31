# Filename: pubfetch.py
# Author: Thomas H. Smith 2017

import xml.dom.minidom
import urllib2
from urllib2 import HTTPError, URLError
#from operator import itemgetter
from multiprocessing.dummy import Pool as ThreadPool
from tqdm import tqdm_notebook
import pandas as pd
from requests.exceptions import ReadTimeout



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
            return (count, link_url)
        else:
            return count


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

