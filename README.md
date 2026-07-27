# 2025_NDBS_ASDsubtypes_dmrichr_ml_George

Methylation profiles of sex-stratified clinical subtypes of autism spectrum disorder using newborn dried blood spot WGBS data.

---

## Project overview

This project investigates sex-stratified clinical and epigenetic subtypes of autism spectrum disorder (ASD) using neurodevelopmental outcome measures and whole-genome bisulfite sequencing (WGBS) methylation data.

The workflow has two main components:

1. **Clinical subtype discovery** using sex-stratified latent profile analysis (LPA)
2. **Differential methylation analysis** using DMRichR based on the subtype-derived comparisons

Machine learning analyses will be added in a later phase.

---

## Project summary

| Feature                    | Detail                                                       |
| -------------------------- | ------------------------------------------------------------ |
| Data type                  | Whole Genome Bisulfite Sequencing                            |
| Cohorts                    | CHARGE and CHDS                                              |
| Genome build               | hg38                                                         |
| Biological sample          | Newborn dried blood spots                                    |
| Study design               | Case-control                                                 |
| Total sample size          | 456 child participants                                       |
| Primary clinical method    | Latent profile analysis                                      |
| Primary methylation method | DMRichR                                                      |
| Main outputs               | Clinical subtypes, DMRs, DMR annotations, enrichment results |

---

## Cohort summary

| Cohort    |     ASD | Typically Developing / Control |   Total |
| --------- | ------: | -----------------------------: | ------: |
| CHARGE    |     161 |                             99 |     260 |
| CHDS      |      84 |                            112 |     196 |
| **Total** | **245** |                        **211** | **456** |

---

## Outcomes and measures

**Primary outcome**

* ASD diagnosis

**Neurodevelopmental measures**

* Mullen Scales of Early Learning
* Vineland Adaptive Behavior Scales
* Social Communication Questionnaire
* Autism Diagnostic Observation Schedule
* Autism Diagnostic Interview-Revised

**Methylation outputs**

* Differentially methylated regions
* Annotated DMRs
* Enrichment results
* Subtype-specific methylation signatures

---

## Workflow

### 1. Latent profile analysis

LPA scripts are located in:

```text
scripts/LPA/
```

Main outputs are written to:

```text
results/LPA/
```

This step identifies sex-stratified ASD and neurotypical clinical profiles and prepares DMRichR comparison files.

### 2. DMRichR analysis

DMRichR scripts are located in:

```text
scripts/DMRichR/
```

Main outputs are written to:

```text
results/DMRichR/
```

This step links Bismark cytosine reports, runs DMRichR comparisons, and generates DMR tables, annotations, enrichment outputs, and reports.

### 3. Future machine learning

Machine learning analyses will be developed later. Current DMRichR-generated machine learning folders are not considered the final project-specific ML workflow.

---

## Repository structure

```text
2025_NDBS_ASDsubtypes_dmrichr_ml_George/
├── analysis/         # notebooks and analysis notes
├── data/             # metadata and processed local data; sensitive files ignored
├── docs/             # detailed workflow and methods documentation
├── reproducibility/  # environment files and Pixi setup
├── results/          # generated LPA and DMRichR outputs; large files ignored
├── scripts/          # LPA and DMRichR scripts
├── tests/            # smoke tests and example testing scripts
├── CHANGELOG.md
├── CITATION.cff
├── LICENSE
└── README.md
```

---

## Quick start

### Clone repository

```bash
git clone https://github.com/dreusebio/2025_NDBS_ASDsubtypes_dmrichr_ml.git
cd 2025_NDBS_ASDsubtypes_dmrichr_ml
```

---

## Reproducible environment

Environment files are located in:

```text
reproducibility/
├── env/
└── pixi/
```

To install the Pixi environment:

```bash
cd reproducibility/pixi
pixi install --concurrent-downloads 1 --concurrent-solves 1
```

To install DMRichR:

```bash
pixi run install-dmrichr
```

To check available Pixi tasks:

```bash
pixi task list
```

---

## Running the analysis

Run LPA scripts from the project root:

```bash
Rscript scripts/LPA/01_data_cleaning.R
Rscript scripts/LPA/02_recode_variables.R
Rscript scripts/LPA/03_subset_neuromeasures.R
Rscript scripts/LPA/10_prepare_DMRichR_inputs.R
```

Run DMRichR setup:

```bash
bash scripts/DMRichR/01_link_bismark_reports.sh
```

Additional DMRichR settings are stored in:

```text
scripts/DMRichR/config_dmrichr.R
```

For detailed workflow instructions, see:

```text
docs/analysis_workflow.md
```

---

## Data and file tracking

Raw WGBS data, cytosine reports, `.rds`, `.RData`, DMRichR intermediate objects, and participant-level metadata are not intended to be committed to GitHub.

These files should be excluded using `.gitignore`.

The repository is intended to track:

* scripts
* configuration files
* documentation
* environment files
* small non-sensitive summary outputs

---

## Team members

* Janine LaSalle
* George Eusebio Kuodza
* Giovani Pimentel-Solorio

---

## Citation

If using this workflow, please cite or acknowledge:

* DMRichR
* Bismark
* DSS / bsseq and relevant Bioconductor tools
* CHARGE and CHDS cohort resources where appropriate

---

## License

MIT License. See `LICENSE`.

---

## Contact

**George Eusebio Kuodza**
Postdoctoral Scholar
University of California, Davis
[gekuodza@health.ucdavis.edu](mailto:gekuodza@health.ucdavis.edu)
