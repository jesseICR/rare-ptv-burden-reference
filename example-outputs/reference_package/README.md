# 1000G EUR Reference Package

This directory contains a precomputed public example reference package from the
full autosomal 1000G unrelated EUR run.

Key files:

- `burden_by_sample.tsv.gz`: one row per 1000G unrelated EUR sample per PTV tier.
- `qualifying_variants_by_tier.tsv.gz`: qualifying sample-variant-tier rows.
- `reference_distribution.tier_summary.tsv`: distribution summary by tier.
- `reference_distribution.top_carriers.tsv`: highest-burden reference samples by tier.
- `reference_distribution.summary.json`: machine-readable distribution summary.
- `gene_scores_shet_constraint_expression.tsv`: gene-level `s_het`, constraint, and GTEx CNS expression annotations.
- `tiers.tsv`: frozen tier definitions used for this package.
- `checksums.tsv`: SHA-256 checksums for the package files.
- `reference_manifest.json`: resource/version metadata for the package.

The package is intended as a reproducible public reference output and as a small
example of the files produced by the pipeline.
