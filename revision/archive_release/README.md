# Archive release

Files for the permanent archive of the revised analysis, built by `build_release.sh` from a commit of this
repository. Archive files and their checksums are not committed (`.gitignore`); each archive records the commit it
was built from in `REVISED_COMMIT.txt`.

## Versions

| Version | Commit | Content |
|---|---|---|
| Analysis version | `55e41b3c030b95c965a6085c28e11e830006ada9` (2026-11-01) | Code of the final analysis run (run 7, `revision/verification/final_run.md`) |
| Archive commit | Recorded in `REVISED_COMMIT.txt` of each archive | The analysis version plus later licence and documentation commits. `REVISED_COMMIT.txt` lists the files changed since the analysis version; no file under `pipeline/` or `tests/`, `run_all.sh`, or any file under `revision/results/`, `revision/figures/`, `revision/tables/` or `revision/source_data/` changes |
| Original submission | tag `v1-submitted` = `e0d6bb527369d8e70b63356d149a579db82567fe` | Scripts posted at initial submission (`legacy/`) |

## Files

| File | Contents |
|---|---|
| `GPR153-GBM-Analysis-<commit>-source.tar.gz` | `git archive` of the archive commit: pipeline, tests, `run_all.sh`, READMEs, licences, legacy scripts, environment specifications, manifests, results, figures, tables, source data, verification records |
| `GPR153-GBM-Analysis-<commit>-submission-package.zip` | The submission package (`revision/submission_package/`): results, tables, figures, source data, manifests, environment, verification, literature, plan, licences; `PACKAGE_CONTENTS.md` and `SHA256SUMS.txt` inside |
| `SHA256SUMS.txt` | SHA-256 of both files |
| `REVISED_COMMIT.txt` | Analysis version, archive commit, files changed between them, original-submission tag |

## Zenodo record

- DOI 10.5281/zenodo.23872434: **reserved for a draft record, not yet published**. The DOI resolves only after
  the record is published.
- Licences chosen by the authors: MIT for author-created software; CC BY 4.0 for author-created figures, result
  tables and non-software documentation; third-party materials keep their original terms (`LICENSING.md`).
- Upload the two archive files of the current `SHA256SUMS.txt`, replacing any earlier files in the draft, and
  check the checksums Zenodo displays (MD5) against the files uploaded.
- Related identifiers: the GitHub repository at the archive commit (`https://github.com/reachdky/GPR153-GBM-Analysis/tree/<commit>`)
  and the analysis plan, https://osf.io/qpd2f/.
- After publication, record the publication date in `revision/data_availability_draft.md`; the DOI is already
  stated there.
