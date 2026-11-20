# Proteosim Proteomic Analysis Package
The Proteosin Package can be used to simulate proteomics experiments based on protein sequences. It includes functions for simulating protein digestion, liquid chromatography and ionisation of the peptides, fragmentation and mass spectrometry analysis.
The code is under active development.

## installation steps
`pip install -r requirements.txt` followed by `pip install .`

## In- and Outputs
Expected inputs: FASTA proteins
Intermediate processing steps: digestion, chromatography, MS simulation), and typical outputs (mz tables, fragment spectra

## Functions
### file_handling
digest_protein_sequence: cleave a protein sequence using an enzyme cleavage pattern
digest_protein_collection:
compute_sequence_coverage:

### liquid_chromatography
predict_lc_retention_times,
plot_retention_time,
select_retention_time_window,

### mass_spectra_simulation
    calculate_mol_mass,
    calculate_mol_mass_collection,
    calculate_mz_collection,
    plot_spectrum,
    fragment_peptide

### end-to-end example
`tutorials/ms_experiment_final.ipynb` notebook