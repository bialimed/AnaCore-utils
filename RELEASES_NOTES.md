# Release 3.7.0 [DEV]

### Improvements:
  * Use user selected annotation tag in header of empty variants file in
  `bin/VEPWrapper.py`.


# Release 3.6.0 [2022-10-01]

### Improvements:
  * Add scripts to manage microsatellites classification:
    * `bin/microsatCreateModel.py`
    * `bin/microsatLenDistrib.py`
    * `bin/microsatMergeResults.py`
    * `bin/microsatMsingsClassif.py`
    * `bin/microsatMSIsensorproProClassif.py`
    * `bin/microsatSklearnClassify.py`
    * `bin/microsatStatusToAnnot.py`
  * Add storage for flag of main alternative annotation in `bin/VCFToJSON.py`.
  This flag is optional and is retrieved from PICK field in annotation (example:
  VCF coming from VEP with flag_pick option).

### Bug fixes:
  * Fix bug with secondary mapping for aligners producing extra lines in the BAM
  (like STAR). These secondary mapping were taken into account and will no longer
  be taken into account.
  * Fix bug in annotating precise breakends with CIPOS > 0 on opposite strands
  and no overlap with an exon edge in `bin/annotBND.py`.


# Release 3.5.1 [2022-03-19]

### Improvements:
  * Update scipy from `1.2.1` to `1.7.3`.

### Bug fixes:
  * Fix bug with missing type of quality parameters in `bin/inspectBND.py` and
  `bin/shallowsAnalysis.py`.


# Release 3.5.0 [2022-03-10]

### Improvements:
  * Update pysam from `0.15.3` to `0.18.0`.
  * Increase execution speed for `bin/inspectBND.py`.
  * Change JSON structure from `bin/fusionsToJSON.py` to improve usage in
  document databases: replacement of the dict feature_by_id by a list with _id
  in each feature.
  * Now `bin/evalVariantControl.py` takes into account false positives variants
  by default. The previous behavior can be recovered with the `--only-expected`
  option: exclude FP from the analysis.
  * Manage annotations out of transcript in `bin/fixHGVSMutalyzer.py`.
  * Add script to add UMI sequence in SAM RX tag from reads ID:
  `bin/setUMITagFromID.py`.

### Bug fixes:
  * Fix bug with refskip for spliced alignment in `bin/mergeCoOccurVar.py`.
  * Fix bug with insertion and substitution of same nucleotid at same position
  in `bin/evalVariantControl.py`. Example: with *1:587874=A/T* and *1:587874=./T*
  only the last in VCF was taking into account.
  * Fix bug on depths in `bin/depthsPanel.py`, `bin/inspectBND.py`,
  `bin/mergeVCF.py`, `bin/mergeVCFAmpli.py` and `bin/shallowsAnalysis.py`. In
  overlapping read pairs, overlapping bases were only taken into account for one
  of the reads.


# Release 3.4.0 [2021-10-13]

### Improvements:
  * Add script to report depths on panel: `bin/depthsPanel.py`.
  * Add script to create known fusions partners database: `bin/buildKnownBNDDb.py`.
  * Add script to create data for saturation curve `bin/saturationCurve.py`.

### Bug fixes:
  * Fix bug in `bin/filterBND.py` when a gene SYMBOL is None (in some lncRNA).


# Release 3.3.0 [2020-04-28]

### Improvements:
  * Replace non-allele specific COSMIC annotations produced by VEP to
  allelle-specific annotations in `bin/fixVEPAnnot.py`.


# Release 3.2.0 [2020-11-11]

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
  * Fix critical bug in `bin/shallowsAnalysis.py`: missing shallow area when
  there is no read on the area.
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
