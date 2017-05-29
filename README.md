# PhosphoproteomicsAnalysis

This is a set of simple scripts that I wrote for the purpose of parsing/aggregating, filtering, and visualizing open-access RNA-Seq gene expression data obtained from the National Cancer Institute Genomic Data Commons (NCI-GDC, <https://gdc-portal.nci.nih.gov>).  I have specifically used it for the TCGA project datasets.

## Getting Started


#### Python Requirements
- Python 2.7 with the following packages:
- `Pandas` 
- `Tqdm`
- `numpy`
- `scipy.stats` 

Aside from Tqdm, all of these are included in the Anaconda distribution of Python 2.7 (https://www.continuum.io/downloads)

#### Installation


#### Other Requirements
 For functional site annotation, you will need to download three files from the PhosphoSitePlus database (http://www.phosphosite.org/staticDownloads.action).  After registering for a free account, download the following files:
 
 - Phosphorylation_site_dataset.gz (all reported functional sites)
 - Phosphosite_seq.fasta.gz (reference sequences)
 - Regulatory_sites.gz (all reported functional sites)

Unzip each file and  save to `PhosphoProTools/src/phosphoprotools/data`, make sure file names match those shown below:
	
	PhosphoProTools
	├── LICENSE
	├── MANIFEST.in
	├── README.md
	├── setup.py
	└── src
	    └── phosphoprotools
	        ├── __init__.py
	        ├── data
	        │   ├── Phosphorylation_site_dataset
	        │   ├── Phosphosite_seq.fasta
	        │   └── Regulatory_sites
	        ├── piscoreanalysis.py
	        ├── preprocessing.py
	        ├── pubfetch.py
	        ├── siteannotation.py
	        └── synonymsfetch.py


## Authors

* **Thomas Smith** - [ThomasHSmith](https://github.com/ThomasHSmith)


## License
This work is licensed under the MIT License - see LICENSE.md for details.
