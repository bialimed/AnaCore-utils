# Release 3.2.0 [DEV]

### Improvements:
  * Manage new ID format from COSMIC in `bin/VCFToJSON.py`: COSM* -> COS*.
  * Add scripts to prepare data mining of fusions from breakends annotations:
  `bin/ensemblInterProToGFF.py` to generate the proteins domains annotations file
  and `bin/inspectBND.py` to produce a detailed annotation of each breakend.
  * Add management for included or near full overlapped amplicons in
   `bin/addAmpliRG.py`.
  * Add the shallow gene file as output from `bin/shallowsAreas.py`. It summarizes
  by gene the shallow areas and potentially missed variants.

### Bug fixes:
  * Fix bug in `bin/fixVEPAnnot.py` when alleles are in lower case.
  * Fix bug in writeJSON of `bin/shallowsAnalysis.py` when optional parameters
  input-annotations and inputs-variants are not provided.

# Release 3.1.0 [2020-05-13]

### Improvements:
  * Add count by target in `--output-summary` file from `bin/addAmpliRG.py`.
  * Add `--single-mode` for process single-end alignments in `bin/addAmpliRG.py`.
  * Replace `bin/fusionCatcherToVCF.py` by `bin/fusionsToVCF.py` to add management
  of Arriba and STAR-Fusion.

### Bug fixes:
  * Fix bug in `bin/splitBAMByRG.py` with change in pysam.AlignmentHeader v1.5:
  `reader.header.copy()` to `reader.header.to_dict()`.

# Release 3.0.0 [2020-01-03]
Split AnaCore project in:
* AnaCore (python libraries),
* AnaCore-utils (scripts),
* AnaCore-web (JS libraries and web components).

Replace section workflows managers utilities by the project AnaCore-snakemake
(AnaCore-jflow is dropped).

# Release 1.0.0 [2017‑10‑07]
First release
