# 1000G EUR Rare PTV Burden Reference

Build a reproducible 1000 Genomes unrelated European reference distribution for rare protein-truncating variant burden on GRCh38.

The goal is an auditable research reference package: the same frozen resources,
variant definitions, and burden formula can be reused by downstream analyses
without changing the reference-panel computation.

The workflow is intentionally idempotent: each stage writes a completion sentinel under `work/.done/`, so rerunning the same command skips completed work unless `--force` is supplied.

## Motivation

Large biobank exome studies suggest that rare damaging coding variants add
phenotypic information beyond common-variant polygenic scores. Chen et al.
showed that cognitive-function polygenic scores and carrier status for rare PTVs
or damaging missense variants in pLI > 0.9 genes have additive associations with
educational attainment and verbal-numerical reasoning.

<p align="center">
  <img src="https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41588-023-01398-8/MediaObjects/41588_2023_1398_Fig5_HTML.png" alt="Chen et al. figure showing additive contributions of cognitive-function polygenic score and rare damaging coding variant carrier status to educational attainment and verbal-numerical reasoning" width="850">
</p>

<p align="center"><em>Chen et al. 2023, Fig. 5: common-variant PRS and rare PTV or damaging missense carrier status contribute additively to EDU and VNR.</em></p>

Tian et al. similarly showed that depression prevalence increases with
polygenic risk score and is higher among carriers of rare PTVs or damaging
missense variants, supporting a joint common-plus-rare architecture.

<p align="center">
  <img src="https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-024-45774-2/MediaObjects/41467_2024_45774_Fig3_HTML.png" alt="Tian et al. figure showing depression prevalence across PRS percentiles stratified by rare PTV or damaging missense carrier status" width="850">
</p>

<p align="center"><em>Tian et al. 2024, Fig. 3: depression risk varies with PRS percentile and rare damaging coding-variant carrier status.</em></p>

This repository builds a reproducible 1000G EUR reference distribution for rare
PTV burden so downstream analyses can compare rare-variant burden against a
stable public reference panel.

## Quick Start

Run a chr22 smoke test against an existing local gnomAD v4.1 joint VCF directory:

```bash
bash main.sh \
  --chromosomes 22 \
  --gnomad-dir /path/to/gnomad_v4.1_joint \
  --kgp-dir /path/to/1000g_cache
```

The gnomAD directory must contain files named:

```text
gnomad.joint.v4.1.sites.chr22.vcf.bgz
gnomad.joint.v4.1.sites.chr22.vcf.bgz.tbi
```

For all autosomes, prefer the sequential runner on machines where raw 1000G
VCFs will not all fit comfortably on the same filesystem:

```bash
bash scripts/run_chromosomes_sequential.sh \
  --chromosomes 1-22 \
  --gnomad-dir /path/to/gnomad_v4.1_joint \
  --kgp-dir /path/to/1000g_cache \
  --download-missing \
  --download-cadd \
  --enable-gtex
```

The plain `main.sh` all-autosome command is useful when all raw chromosome
VCFs are already present:

```bash
bash main.sh \
  --chromosomes 1-22 \
  --gnomad-dir /path/to/gnomad_v4.1_joint \
  --kgp-dir /path/to/1000g_cache
```

Use explicit download flags when public resources are absent:

```bash
bash main.sh \
  --chromosomes 22 \
  --gnomad-dir /path/to/gnomad_v4.1_joint \
  --kgp-dir /path/to/1000g_cache \
  --download-missing \
  --enable-gtex \
  --stream-1000g
```

Use `--force` to rerun completed stages.

## Setup

From a fresh checkout, the only things expected to exist before setup are:

- `bash`
- `git`
- `curl`
- `micromamba`, `mamba`, `conda`, or another way to install Bioconda tools
- enough local or mounted storage for the ignored resource directories

The pipeline itself uses these runtime tools:

- `bcftools`
- `tabix`
- `bgzip`
- `samtools`
- VEP 115 with a GRCh38 cache
- Python packages from `requirements.txt`

Those runtime tools do not need to be preinstalled system-wide. Install them
locally into ignored paths, then rerun the pipeline. The pipeline automatically
uses `tools/envs/ptv/bin` when that environment exists.

For example:

```bash
micromamba create -y -p tools/envs/ptv \
  -c conda-forge -c bioconda \
  python=3.11 pip git curl ensembl-vep=115.2 bcftools samtools htslib tabix \
  perl-list-moreutils perl-set-intervaltree perl-dbi perl-dbd-mysql \
  perl-html-tableextract

export PATH="${PWD}/tools/envs/ptv/bin:${PATH}"

python -m pip install -r requirements.txt

perl tools/envs/ptv/bin/vep_install \
  -a cf \
  -s homo_sapiens \
  -y GRCh38 \
  -c resources/vep/cache \
  --CACHE_VERSION 115 \
  --CONVERT \
  --NO_UPDATE \
  --USE_HTTPS_PROTO

bash scripts/setup_vep_loftee.sh
bash scripts/download_loftee_aux.sh
```

VEP, its GRCh38 cache, LOFTEE, and LOFTEE auxiliary resources are large
third-party resources. They are intentionally stored under ignored `tools/` and
`resources/` paths rather than committed to this repository.

## Main Definitions

The default rare PTV definitions use VEP/LOFTEE annotations, QC-filtered 1000G unrelated EUR allele count, gnomAD v4.1 external allele frequency, CADD v1.7 PHRED, pLI, and optional GTEx brain/CNS expression percentiles.

Key tiers include:

- `broad_ptv_gnomad1e3_ac1`: stop gained, frameshift, splice donor, or splice acceptor; gnomAD AF < 1e-3; 1000G unrelated EUR AC <= 1.
- `loftee_hc_ptv_gnomad1e3_ac1`: broad PTV plus LOFTEE high confidence.
- `chen_broad_ptv_gnomad1e5_ac1_pli90`: broad PTV, gnomAD AF < 1e-5, 1000G unrelated EUR AC <= 1, and pLI >= 0.9.
- `hc_ptv_gnomad1e5_ac1_pli90`: LOFTEE high-confidence PTV, gnomAD AF < 1e-5, 1000G unrelated EUR AC <= 1, and pLI >= 0.9.
- `gardner_private_proxy_ptv_gnomad_absent_ac1`: Gardner primary-style private proxy: LOFTEE high confidence, CADD v1.7 PHRED > 25, not last exon/intron, absent from gnomAD v4.1, and 1000G unrelated EUR AC <= 1.
- `gardner_core_ptv_gnomad1e3_ac1`, `gardner_core_ptv_gnomad1e4_ac1`, `gardner_core_ptv_gnomad1e5_ac1`, `gardner_core_ptv_gnomad2e6_ac1`: Gardner-style strict PTV masks with gnomAD AF < 1e-3, < 1e-4, < 1e-5, and < 2e-6. Gardner et al.'s primary burden was private/singleton; the gnomAD AF thresholds are rare-threshold sensitivity tiers.
- `rare_pli30_ptv_gnomad1e3_ac1`: broad rare PTV with pLI > 0.3.
- `modern_pli90_ptv_gnomad1e3_ac1`: Gardner-style constrained-gene sensitivity tier with pLI >= 0.9.

The `s_het` score is computed per sample and tier as:

```text
1 - product(1 - shet_gene)
```

Multiple qualifying variants in the same gene count once for `s_het` burden.

1000G carrier calls are filtered before allele counts and carrier rows are
computed. By default, a site must have `FILTER=PASS`, each counted genotype
must have `GQ >= 20` and `DP >= 10`, and heterozygous `AB` must be between
0.2 and 0.8 when the `AB` field is present. Missing `AB` does not fail a call
unless `KGP_REQUIRE_HET_AB=1` is set.

## PTV Definitions

The pipeline separates the base molecular definition of a protein-truncating
variant from extra filters intended to enrich for clinical or complex-trait
impact. The base broad PTV consequence set is:

```text
stop_gained, frameshift_variant, splice_acceptor_variant, splice_donor_variant
```

The implemented tiers map to the literature masks as follows:

| Definition | Paper motivation | Implemented tier | Criteria |
| --- | --- | --- | --- |
| Broad rare PTV | Baseline consequence-only mask used before stricter impact filters | `broad_ptv_gnomad1e3_ac1`, `broad_ptv_gnomad1e4_ac1`, `broad_ptv_gnomad1e5_ac1` | Broad PTV consequences, 1000G unrelated EUR AC <= 1, and gnomAD AF below the tier threshold. |
| Ganna-style constrained disruptive variant | Ganna et al. studied ultra-rare disruptive variants in highly constrained genes, with brain-expression sensitivity analyses. | `chen_broad_ptv_gnomad1e5_ac1_pli90` is the closest current constrained broad-PTV mask. | Broad PTV consequences, gnomAD AF < 1e-5, 1000G unrelated EUR AC <= 1, pLI >= 0.9. GTEx CNS expression percentiles are reported rather than used as a hard filter. |
| Chen-style constrained rare PTV | Chen et al. defined rare PTVs from VEP stop-gain/splice-disruptive/frameshift consequences, used UKB MAF < 1e-5, and stratified genes by pLI >= 0.9. | `chen_broad_ptv_gnomad1e5_ac1_pli90` | Broad PTV consequences, gnomAD AF < 1e-5 as the public-reference proxy for UKB MAF < 1e-5, 1000G unrelated EUR AC <= 1, pLI >= 0.9, no CADD requirement. |
| LOFTEE high-confidence constrained rare PTV | Tian et al. and van den Berg et al. use rare PTV masks with LOFTEE high-confidence filtering and pLI-stratified burden analyses. | `hc_ptv_gnomad1e5_ac1_pli90` | Broad PTV consequences, LOFTEE HC, gnomAD AF < 1e-5, 1000G unrelated EUR AC <= 1, pLI >= 0.9, no CADD requirement. |
| Gardner primary-style private PTV proxy | Gardner et al.'s primary burden used private/singleton deleterious variants with stricter PTV filters and `s_het` burden aggregation. | `gardner_private_proxy_ptv_gnomad_absent_ac1` | Broad PTV consequences, LOFTEE HC, CADD v1.7 PHRED > 25, not last exon/intron, absent from gnomAD v4.1, and 1000G unrelated EUR AC <= 1. This is a public-reference proxy, not an exact UK Biobank-private mask. |
| Gardner-style rare-threshold sensitivity PTV | Gardner et al. also reported rare-threshold sensitivity analyses at MAF < 1e-3, < 1e-4, and < 1e-5; this pipeline adds a stricter < 2e-6 sensitivity tier. | `gardner_core_ptv_gnomad1e3_ac1`, `gardner_core_ptv_gnomad1e4_ac1`, `gardner_core_ptv_gnomad1e5_ac1`, `gardner_core_ptv_gnomad2e6_ac1` | Broad PTV consequences, LOFTEE HC, CADD v1.7 PHRED > 25, not last exon/intron, 1000G unrelated EUR AC <= 1, and gnomAD AF below the tier threshold. These tiers are expected to produce larger burden distributions than the private proxy. |
| pLI 0.3 sensitivity tier | A more permissive constrained-gene screen requested for broader phenotype relevance. | `rare_pli30_ptv_gnomad1e3_ac1` | Broad PTV consequences, gnomAD AF < 1e-3, 1000G unrelated EUR AC <= 1, pLI > 0.3. |

Two implementation details differ from exact paper reproductions. First, this
reference panel uses public 1000G and gnomAD frequencies, so `gnomAD AF < 1e-5`
is a proxy for study-internal UK Biobank MAF thresholds. Second, the parser
retains PTV transcript consequences emitted by VEP/LOFTEE and records canonical
and MANE annotations, but the tier table does not currently enforce a
canonical-only or MANE-only transcript mask.

## Data Sources

The main 1000G input is the raw 20201028 30x high-coverage genotype VCF
collection:

```text
20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr*.recalibrated_variants.vcf.gz
```

The phased high-coverage panel is intentionally not used for singleton tiers,
because it removes singleton variants.

1000G is supplied by path at runtime:

```bash
--kgp-dir /path/to/1000g_cache
```

The pipeline stores raw 1000G files and metadata under:

```text
/path/to/1000g_cache/raw/
```

This directory should be on a filesystem with enough space for the requested raw
chromosome VCFs. `--download-missing` downloads absent files idempotently into
that cache, so interrupted downloads can be resumed and existing chromosomes are
reused. For full reference builds on finite local storage, use
`scripts/run_chromosomes_sequential.sh --download-missing`; it processes one
chromosome at a time and can optionally remove each raw chromosome VCF after a
successful chromosome pass with `--delete-raw-after-success`.

`--stream-1000g` reads the indexed public VCFs remotely and extracts only
coding/splice intervals, avoiding whole-chromosome downloads. This requires a
`bcftools`/htslib build with HTTPS support and can be slower than local
sequential mode on large all-autosome runs.

gnomAD is supplied by path at runtime:

```bash
--gnomad-dir /path/to/gnomad_v4.1_joint
```

If the files are not already present, download them into a directory on a
filesystem with enough space, then pass that same directory as `--gnomad-dir`:

```bash
bash scripts/download_gnomad.sh \
  --chromosomes 1-22 \
  --out-dir /path/to/gnomad_v4.1_joint
```

Expected files:

```text
gnomad.joint.v4.1.sites.chr1.vcf.bgz
gnomad.joint.v4.1.sites.chr1.vcf.bgz.tbi
...
gnomad.joint.v4.1.sites.chr22.vcf.bgz
gnomad.joint.v4.1.sites.chr22.vcf.bgz.tbi
```

The Weghorn `s_het` resource is taken from the public
`HurlesGroupSanger/UKBBFertility` raw data bundle and normalized from
`rawdata/genelists/shet.weghorn.txt`, using the `s_het_drift` column.

CADD defaults to GRCh38 v1.7. For production runs, use `--download-cadd` to
download local CADD v1.7 SNV and indel lookup tables. If local tables are
absent, the pipeline can query the CADD GRCh38-v1.7 API for post-filter SNVs
only. That fallback is meant for small smoke tests; the CADD API is
experimental and not intended for thousands of lookups.

LOFTEE is run with explicit GRCh38 auxiliary resources:

```text
resources/loftee/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
resources/loftee/loftee_data/GRCh38/human_ancestor.fa.gz
resources/loftee/loftee_data/GRCh38/human_ancestor.fa.gz.fai
resources/loftee/loftee_data/GRCh38/human_ancestor.fa.gz.gzi
```

GTEx CNS/brain expression is optional and uses GTEx v8 median tissue TPM. CNS
expression percentile is computed across genes using all GTEx tissue columns
whose names contain `Brain`.

## Pipeline Stages

`main.sh` runs these stages:

1. Write tool versions.
2. Optionally download public resources.
3. Validate gnomAD, the user-supplied 1000G cache, reference, VEP, LOFTEE auxiliary files, pLI, and `s_het` inputs.
4. Build unrelated EUR sample lists from 1000G metadata.
5. Build a protein-coding CDS/splice BED from GENCODE.
6. Extract and normalize candidate chr-level variants from raw 1000G VCFs.
7. Run VEP and LOFTEE.
8. Parse PTV transcript consequences.
9. Compute QC-filtered unrelated EUR allele counts and carriers.
10. Query gnomAD only for PTVs with 1000G unrelated EUR AC <= 1.
11. Query CADD v1.7 only for variants passing earlier rarity filters.
12. Apply tier definitions.
13. Compute per-sample burden and empirical percentiles.
14. Build QC and reference distribution summaries.
15. Build a frozen reference manifest.

`scripts/run_chromosomes_sequential.sh` wraps the same stages one chromosome
at a time. It is intended for full autosome builds where the raw 1000G cache is
kept on a large external or mounted filesystem specified by `--kgp-dir`. With
`--stream-1000g`, raw 1000G chromosome VCFs are not downloaded at all; in
local-raw mode, each raw VCF can be deleted after the chromosome has been
reduced to the small intermediate files needed for final merging.

Each stage writes a sentinel under `work/.done/`. Logs are in `logs/`.

## Outputs

The generated reference package is written to:

```text
results/reference_package/
```

The repository also includes a precomputed public example package at:

```text
example-outputs/reference_package/
```

This example package was generated from the full autosomal 1000G unrelated EUR
reference run and can be used immediately for inspection, tutorials, or
downstream comparison without rerunning the full workflow.

Important files:

- `burden_by_sample.tsv.gz`
- `qualifying_variants_by_tier.tsv.gz`
- `gene_scores_shet_constraint_expression.tsv`
- `reference_distribution.tier_summary.tsv`
- `reference_distribution.top_carriers.tsv`
- `reference_distribution.summary.json`
- `reference_manifest.json`
- `tool_versions.tsv`
- `tiers.tsv`
- `checksums.tsv`

QC outputs are written under `results/qc/`.

## Interpreting Results

The main burden table has one row per 1000G unrelated EUR sample per tier. Key
columns:

- `n_qualifying_variants`: unique qualifying variants carried by the sample.
- `n_qualifying_genes`: unique genes hit by qualifying variants.
- `shet_burden`: `1 - product(1 - shet_gene)` across hit genes.
- `empirical_percentile_within_1000G_EUR`: percentile of `shet_burden` within
  the reference panel for the same tier.
- `mean_cns_expression_percentile` and `max_cns_expression_percentile`: GTEx
  CNS/brain expression summaries for hit genes when GTEx is enabled.

Sparse singleton tiers are expected to have many zero-burden samples.

## Resource Notes

- Use the raw 20201028 1000G high-coverage genotype VCFs for singleton/AC<=1 work.
- Do not use the phased high-coverage panel for singleton analyses; that panel removes singleton variants.
- 1000G genotype QC defaults are controlled by `KGP_REQUIRE_SITE_PASS`,
  `KGP_MIN_GQ`, `KGP_MIN_DP`, `KGP_HET_AB_MIN`, `KGP_HET_AB_MAX`, and
  `KGP_REQUIRE_HET_AB` in `config/default.env`.
- CADD defaults to v1.7 GRCh38. Lookup is performed only after earlier filters have reduced the variant set. Local CADD tables are recommended for full autosome runs.
- LOFTEE high-confidence calls require the GRCh38 GERP bigWig and ancestral
  FASTA; `scripts/download_loftee_aux.sh` downloads those files.
- The Weghorn/Gardner `s_het` table is required for production `s_het` burden. The setup script fails clearly if no usable table is provided.

## Notes And Limitations

- This is a research reference distribution, not a validated clinical score.
- The reference panel size is modest; use empirical percentiles rather than
  normal approximations.
- gnomAD v4.1 differs from the historical resources used in older publications;
  outputs are Gardner-style, not exact reproductions.
- CADD API fallback is intended for small post-filter SNV smoke tests. Full
  local CADD mirroring is recommended for production all-autosome runs,
  especially if frameshift indels need complete CADD coverage.

## References

This workflow is a modernized implementation of rare PTV burden ideas from:

- Ganna A, Genovese G, Howrigan DP, et al. Ultra-rare disruptive and damaging
  mutations influence educational attainment in the general population. Nature
  Neuroscience. 2016;19:1563-1565. doi:10.1038/nn.4404.
- Gardner EJ, Neville MDC, Samocha KE, et al. Reduced reproductive success is
  associated with selective constraint on human genes. Nature.
  2022;603:858-863. doi:10.1038/s41586-022-04549-9.
- Chen C-Y, et al. The impact of rare protein coding genetic variation on adult
  cognitive function. Nature Genetics. 2023;55:927-938.
  doi:10.1038/s41588-023-01398-8.
- Tian R, Ge T, Kweon H, et al. Whole-exome sequencing in UK Biobank reveals
  rare genetic architecture for depression. Nature Communications.
  2024;15:1755. doi:10.1038/s41467-024-45774-2.
- van den Berg DM, Huang W, Malawsky DS, et al. Imputation of fluid
  intelligence scores reduces ascertainment bias and increases power for
  analyses of common and rare variants. medRxiv. 2025.
  doi:10.1101/2025.06.18.25329418.
