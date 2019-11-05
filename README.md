# AnaCore-Utils

## Description
Scripts, rules/components and resources for workflow managers and web components and JS libraries developped for NGS in IUCT Oncopole.

### AnaCore-Utils folder
The application folder has the following structure:

    <APP_DIR>/
    ├── bin/                  # Scripts
    ├── workflows_managers/
    │   ├── jflow/
    │   │   └── components/   # Components used in jFlow
    │   └── snakemake/        # Rules used in snakemake      
    ├── README.md
    └── web/
        ├── components/       # Re-usable Vue.js components for bioinformatic's results
        └── models/           # JS objects for manipulate biological/analysis entities (variants, MSI, ...)

## Installation
AnaCore-Utils can be installed in two ways.

### From conda
Install binaries in environments. Resources for workflows managers and web will not be downloaded.

    conda install AnaCore-utils

### From sources
All the resources will be downloaded and binaries can be added to the environment.

    git clone https://github.com/bialimed/AnaCore-utils.git
    cd anacore-utils
    git checkout $RELEASE

    # Install bin in $PATH
    python setup.py install --user

## License
GNU GPL v3

## Copyright
2017 Laboratoire d'Anatomo-Cytopathologie de l'Institut Universitaire du Cancer Toulouse - Oncopole

## Contact
escudie.frederic@iuct-oncopole.fr
