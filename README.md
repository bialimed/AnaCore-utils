# AnaCore

## Description
Anapath Core contains binaries, libraries and resources developped for
Anatomo-Cytopathologie department of IUCT Oncopole.

## Installation
### Download

* [user way] Downloads the latest released versions from `https://bitbucket.org/fescudie/anacore/downloads/?tab=tags`.
* [developper way] Clones the repository from the latest unreleased version: `git clone https://bitbucket.org/fescudie/anacore.git`.

### AnaCore folder
The application folder has the following structure:

    <APP_DIR>/
    ├── bin/               # Scripts
    ├── jflow/
    │   └── components/    # Wrappers used in workflow manager (jflow)
    ├── lib/               # Librairies used in binaries
    ├── README.md
    └── web/
        ├── components/    # Re-usable Vue.js components for bioinformatic's results
        └── models/        # JS objects for manipulate biological/analysis entities (variants, MSI, ...)

## License
GNU GPL v3

## Copyright
2017 Laboratoire d'Anatomo-Cytopathologie de l'Institut Universitaire du Cancer
Toulouse - Oncopole

## Contact
escudie.frederic@iuct-oncopole.fr
