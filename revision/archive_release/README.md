# Archive release (prepared, not deposited)

Files for a permanent archive of the revised analysis. They are built from the revised commit
(`REVISED_COMMIT.txt`) and are not yet deposited: a deposit needs the authors' Zenodo (or OSF) account, and no
DOI exists until it is made.

| File | Contents |
|---|---|
| `GPR153-GBM-Analysis-<commit>-source.tar.gz` | `git archive` of the revised commit: pipeline, tests, `run_all.sh`, READMEs, legacy scripts, environment specifications, manifests, results, figures, tables, source data, verification records |
| `GPR153-GBM-Analysis-<commit>-submission-package.zip` | The submission package (`revision/submission_package/`) |
| `SHA256SUMS.txt` | SHA-256 of both files |
| `REVISED_COMMIT.txt` | Full commit hash and date; original-submission tag |

## Deposit (authors' account)

Either route gives a DOI; use one.

1. **Zenodo upload.** zenodo.org, New upload: add both archive files; resource type Software (or Dataset);
   title, authors, affiliations and licence as the authors decide; related identifiers: the GitHub repository
   URL and https://osf.io/qpd2f/; version: the revised commit hash. Publish, then record the DOI in the Data
   Availability statement (`{REVISED_COMMIT}` and the DOI placeholder in `revision/data_availability_draft.md`).
2. **Zenodo GitHub integration.** zenodo.org, Account, GitHub: enable `reachdky/GPR153-GBM-Analysis`; then create
   a GitHub release from the revised commit (for example tag `v2-revision`). Zenodo archives the release and
   assigns a DOI. Creating a release and tag needs the authors' approval.

The creator and licence metadata are not filled in here; they are for the authors to state.
