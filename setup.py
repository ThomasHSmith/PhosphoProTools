from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='phosphoprotools',
      version='0.14',
      description='Tools for phosphoproteomics data analysis',
      long_description=('Collection of functions for processing phosphoproteomics data,'\
                       'including annotation of phospho-sites using PhosRS cutoff value,'\
                       'annotation of functional sites, retrieval of synonyms, publication'\
                       'counts, statistical analysis (by pi-score method), and more.'),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Text Processing :: Linguistic',
      ],
      url='https://github.com/ThomasHSmith/phosphoprotools',
      author='Thomas H. Smith',
      author_email='ThomasHoraceSmith@gmail.com',
      license='MIT',
      packages= find_packages('src'),
      package_dir={'': 'src', 'refdbprocessing':'data'},
      package_data={'': ['*.pickle']},
      install_requires=['tqdm',],
      zip_safe=False)
