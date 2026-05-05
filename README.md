# Methylation profiles of Sex-stratified Clinical Subtypes of ASD

## Project overview

This project investigates sex-stratified clinical and epigenetic subtypes of autism spectrum disorder (ASD) using neurodevelopmental outcome measures and whole-genome bisulfite sequencing (WGBS) methylation data. In our first aim, we apply sex-stratified clustering of clinical measures to identify distinct ASD subtypes within each sex. Building on these clinically defined subtypes, our second aim characterizes their epigenetic architecture through differential methylation analysis, then evaluates whether these methylation signatures can predict ASD diagnosis in an independent cohort. 

---
## Project summary

**Data type(s):** WGBS

**Study/cohort:** CHARGE(Childhood Autism Risk from Genetics and Environments) and CHDS(Child Health Developmental Studies)

**Genome build:** hg38

**Primary outputs:** Clinical subtypes, DMRs

---

## Study design

* **Design:** Case-Control
* **Sample size:** 260 child participants from CHARGE and 196 from CHDS
* **Biological sample:** Newborn Dried Blood

**Outcomes:**

* ASD Diagnosis
* Neurodevelopmental Outcomes(Mullen,Vineland,SCQ,ADOS, and ADIR)

**Statistical framework:** Clinical subtypes will be identified using Latent profile Analysis. DMRs will be identified using DMRichr. Machine Learning will be performed with random forrest.

---

## Analysis workflow

---

## Quick start

### Clone
```bash
git clone https://github.com/dreusebio/2025_NDBS_ASDsubtypes_dmrichr_ml.git
cd 2025_NDBS_ASDsubtypes_dmrichr_ml
```
---
## Repository Structure

- `config/` — configuration + sample sheets
- `metadata/` — codebooks, data dictionaries, provenance notes (**committed**)
- `data/` — raw data (**not committed**)
- `scripts/` — entrypoints + SLURM submit scripts
- `analysis/` — downstream statistics/figures (optionally `renv`)
- `docs/` — methods + workflow documentation
- `results/` — generated outputs (**not committed**)
- `reproducibility/` — how to generate environment (** committed**)
- `test/` — examine test data  (** committed**)

---

## Team Members
Janine LaSalle 

George Eusebio Kuodza

Giovani Pimentel-Solorio 

---

## License

MIT License (see `LICENSE`)
