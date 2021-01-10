# Biosynthetic Pathway Generation Repository 

### Introduction
This repository contains code for effective synthetic circuit models written in the [Python](https://www.python.org/downloads/) programming language. 
The algorithm and examples are described in the publication:


### Installation and Requirements
There are several required [Python](https://www.python.org/downloads/) packages that must be installed to use the biochemical pathway generation code. 

package | version | download | documentation 
---: | ---: | --- | ---
  Rdkit | 2019.09.3.0 | [Installation Guide](https://www.rdkit.org/docs/Install.html) | [Documentation](https://buildmedia.readthedocs.org/media/pdf/rdkit/latest/rdkit.pdf).  
  Numpy | 1.17.2  | [Installation Guide](https://numpy.org/install/) | [Documentation](https://numpy.org/doc/stable/) 
  Pandas | 0.25.1 | [Installation Guide](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html) | [Documentation](https://pandas.pydata.org/docs/) 
  Keras | 2.3.1 | [Installation Guide](https://keras.io/)|[Documentation](https://keras.io/guides/)
  urllib3 | 1.25.8 | [Installation Guide](https://pypi.org/project/urllib3/) |[Documentation](https://urllib3.readthedocs.io/en/latest/user-guide.html) 
 scikit-learn | 0.22.2.post1 | [Installation Guide](https://scikit-learn.org/stable/install.html) | [Documentation](https://scikit-learn.org/stable/_downloads/scikit-learn-docs.pdf)
scipy | 1.4.1  | [Installation Guide](https://www.scipy.org/install.html) | [Documentation](https://docs.scipy.org/doc/scipy/reference/)

After the above packages are installed the notebook has to opened in the rdkit environment. Terminal commands for finding the rdkit environment are:

    > conda info --envs conda info -e
    > conda activate my-rdkit-env
    > jupyter notebook 

Lastly, some files required to run the notebook are too big to be maintained in this repository. We have to store these files externally.
Please download the [folder](https://drive.google.com/drive/folders/14nG2eAxLNvol8CD6sGzZzKVpzlixenFM?usp=sharing). The link directs to two sub-folders "Required to run code" and "reaction rules database (2.36 GB)". The prior is required to run the notebook with example pathways. If the transformation is not found, the latter should be downloaded, and rules added.

### How do I execute the biochemical route planning job?
The function required to run a pathway generation job is the:
    
    result = final_main(initialCompound, finalCompound, threshold, hydrogenAdd, types)
    
Argument | Type | Description
--: | -- | -- 
`initial`, `final` | Smile String | The product from which you would want to iterate back to target (this is the reaction's precursor.) For example, when creating a Tyrosine pathway to Morphine, the initial compound and final target compound are Morphine and Tyrosine, respectively|
`threshold` | int | This parameter controls the maximum [Tanimoto threshold](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2914517/) cutoff accepted for any two successive compounds in the pathway. The range of the Tanimoto similarity is [0,1].  Default value: 0.1|
`hydrogen` | Boolean | Depending on whether the reaction rules have hydrogen implicitly or explicitly added, this parameter controls whether we add hydrogen or not.
`types` | string |This parameter controls which deep learning is being chosen for candidate ranking. This has only three accepted values: "reaction", "molecular", "atomic_level", "atomic_spectator_model"

### Funding
The work described was supported by the [Center on the Physics of Cancer Metabolism at Cornell University](https://psoc.engineering.cornell.edu) through Award Number 1U54CA210184-01 from the [National Cancer Institute](https://www.cancer.gov). The content is solely the responsibility of the authors and does not necessarily
represent the official views of the [National Cancer Institute](https://www.cancer.gov) or the [National Institutes of Health](https://www.nih.gov).
