# Methylation profiles of Sex-stratified Clinical Subtypes of ASD

This project investigates sex-stratified clinical and epigenetic subtypes of autism spectrum disorder (ASD) using neurodevelopmental outcome measures and whole-genome bisulfite sequencing (WGBS) methylation data. In our first aim, we apply sex-stratified clustering of clinical measures to identify distinct ASD subtypes within each sex. Building on these clinically defined subtypes, our second aim characterizes their epigenetic architecture through differential methylation analysis, then evaluates whether these methylation signatures can predict ASD diagnosis in an independent cohort. 

---

**Data type(s):** WGBS
**Study/cohort:** CHARGE(Childhood Autism Risk from Genetics and Environments) and CHDS(Child Health Developmental Studies)
**Genome build:** hg38
**Primary outputs:** Clinical subtypes, DMRs

---

## Quick start

### Clone
```bash
git clone https://github.com/dreusebio/2025_NDBS_ASDsubtypes_dmrichr_ml_George.git
cd 2025_NDBS_ASDsubtypes_dmrichr_ml_George
```

---

```

---

## Repository structure

- `config/` — configuration + sample sheets
- `metadata/` — codebooks, data dictionaries, provenance notes (**committed**)
- `resources/` — reference versions + checksums (**committed**)
- `workflow/` — orchestration workflows (**optional**)
- `scripts/` — entrypoints + SLURM submit scripts
- `analysis/` — downstream statistics/figures (optionally `renv`)
- `docs/` — methods + workflow documentation
- `results/` — generated outputs (**not committed**)
- `logs/` — generated logs (**not committed**)

---

## Reproducibility standard

See `docs/lab_reproducibility_standard.md`.

---

## License

MIT License (see `LICENSE`)
