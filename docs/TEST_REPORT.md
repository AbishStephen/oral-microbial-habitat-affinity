# Verification of the Methods-organised release

Verified on 17 October 2026 with Python 3.12 and the pinned scientific environment.

| Check | Result |
|---|---|
| Scientific integrity and organisation tests | 13 passed |
| Retained workflow sources | 133 checksums matched; Python and shell syntax passed |
| Methods layout | All workflow sources catalogued once; eleven section guides; unique runtime destinations |
| Runtime preparation | Organised sources retain exactly the original code after anonymous path substitution, including paths containing spaces |
| Cohort reconstruction | 47 HMTs, 13 genera and 1,121 genomes |
| Full 4,096-permutation feature screen | Identical 14-candidate set and candidate P values; 119 raw nominal, 23 jointly nominal; no all-profile FWER hits |
| Four observed nested prediction models | Every HMT prediction matched the archived output exactly |
| Seed-0 complete prediction null refits | Adjusted and unadjusted results matched exactly for both outcomes |
| Stored full null distributions | All 200/100 replicates per outcome checked; empirical P values recalculated |
| Supplementary coordinate extraction | 30 loci and 206 gene records reconstructed |
| Reference Python figures | Five main figures rebuilt successfully |
| Native R scripts and delivered figures | All six figures rebuilt after author text edits; ten input checks passed; retained PDF text matched all six edited originals |

The full `python run.py check` command passed. After adding the two organisation-specific tests, `python run.py test` passed all 13 tests. The separate `python run.py cohort` command also passed. Machine-readable evidence is in `test_evidence/reorganised_primary_reproduction.json`, `reorganised_permutation_smoke.json` and `reorganisation.json`.

The figure-text revision was rendered with R 4.2.0 on 17 October 2026. All six PDFs were inspected visually; their whitespace-normalised text-token counts matched the author-edited PDFs. All 17 exported source tables remained byte-identical to the preceding release. The R export check verified 206 gene records and 47 held-out predictions; the 13 Python integrity tests passed again. Evidence and exact package versions are in `test_evidence/author_figure_text_reconciliation.json`. Earlier verification remains documented in `PRIOR_VERIFICATION.md`. Figure generation requires a UTF-8 locale for consistent Unicode text in exported tables.

Upstream raw annotation, full phylogenetic inference, external structure services and all supplementary analyses were not rerun. Complete prediction null distributions were checked, not fully recomputed. Hosted GitHub CI has not yet run. The full analytical route is tested on Linux; native Windows execution is not claimed.
