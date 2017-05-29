# Filename: preprocessing.py
# Author: Thomas H. Smith 2017

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

