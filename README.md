# PhosphoproteomicsAnalysis
Collection of functions for processing phosphoproteomics data, including:

 - Annotate phospho-sites using PhosRS cutoff value
 - Retrieve Gene/Protein synonyms from Uniprot
 - Annotate previously reported phospho-sites (as reported by PhosphoSitePlus database) 
 - Annotate functional phospho-sites (as reported by PhosphoSitePlus database)
 - Statistical analysis by pi-score method for differential phosphopeptide abundance between two groups
 - Retrieve publication counts for genes with desired keywords
               
## Getting Started

#### Python Requirements
- Python 2.7 with the following packages:
- `Pandas` 
- `Tqdm`
- `numpy`
- `scipy.stats` 

Aside from Tqdm, all of these are included in the Anaconda distribution of Python 2.7 (https://www.continuum.io/downloads).  Missing packages should be installed automatically as part of the pip install of the phosphoprotools package.

#### Installation
`pip install phosphoprotools`

#### Other Requirements
 For annotation of known/functional phosphosites, you will need to download three files from the PhosphoSitePlus database (http://www.phosphosite.org/staticDownloads.action).  After registering for a free account, download the following files to your Desktop.
 
 - Phosphorylation_site_dataset.gz (all reported functional sites)
 - Phosphosite_seq.fasta.gz (reference sequences)
 - Regulatory_sites.gz (all reported functional sites)

When phosphoprotools is imported for the first time, the above three compessed files will be automaticallly extracted to the package `/data` directory and saved as a .pickle.  After first import the files downloaded to your Desktop can be deleted.
	
	PhosphoProTools
	├── LICENSE
	├── MANIFEST.in
	├── README.md
	├── setup.py
	└── src
	    └── phosphoprotools
	        ├── __init__.py
	        ├── data
	        │   ├── Phosphorylation_site_dataset.pickle
	        │   ├── Phosphosite_seq.fasta.pickle
	        │   └── Regulatory_sites.pickle
	        ├── piscoreanalysis.py
	        ├── preprocessing.py
	        ├── pubfetch.py
	        ├── siteannotation.py
	        └── synonymsfetch.py


## TODO
- Finish up comments
- Add example ipynb
- Add sample data for testing
- Add easy command to update db pickels
- Add ref seq pickle that is permanently saved, updates self
- Add synonynms df pickle that is permanently saved, updates self
- Fix long link problem
- Add pubcount column with number only
- fix publication fetch error when only searching 3 terms
	- remove [TIAB] from each term!!!!!!

## Authors

* **Thomas Smith** - [ThomasHSmith](https://github.com/ThomasHSmith)


## License
This work is licensed under the MIT License - see LICENSE.md for details.
