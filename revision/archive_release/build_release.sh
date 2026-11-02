#!/usr/bin/env bash
# Build the archive files from the current commit (HEAD).
#
# Writes revision/submission_package/ and revision/submission_package.zip, then in revision/archive_release/:
#   GPR153-GBM-Analysis-<HEAD>-source.tar.gz        git archive of HEAD
#   GPR153-GBM-Analysis-<HEAD>-submission-package.zip
#   REVISED_COMMIT.txt, SHA256SUMS.txt
# Stops if the working tree is not clean, or if any analysis file differs from the analysis version
# (analysis code, tests, results, figures, tables, source data). Earlier archive files are moved to
# revision/archive/release_<commit>/, never deleted.
#
# Run from the repository root:
#     bash revision/archive_release/build_release.sh
set -euo pipefail

ANALYSIS=55e41b3c030b95c965a6085c28e11e830006ada9
H=$(git rev-parse HEAD)
A=revision/archive_release
P=revision/submission_package

[ -z "$(git status --porcelain --untracked-files=no)" ] || { echo "working tree not clean" >&2; exit 1; }
changed=$(git diff --name-only "$ANALYSIS" "$H")
if echo "$changed" | grep -E '^(pipeline/|tests/|run_all\.sh$|revision/(results|figures|tables|source_data)/)'; then
    echo "analysis files changed since $ANALYSIS (listed above); not a documentation-only release" >&2
    exit 1
fi

# Keep earlier archive files
for f in "$A"/GPR153-GBM-Analysis-*-source.tar.gz; do
    [ -e "$f" ] || continue
    old=$(basename "$f" | sed -E 's/^GPR153-GBM-Analysis-([0-9a-f]{40})-source\.tar\.gz$/\1/')
    [ "$old" = "$H" ] && continue
    dest=revision/archive/release_$old
    mkdir -p "$dest"
    mv "$A"/GPR153-GBM-Analysis-"$old"-* "$dest"/
    for g in SHA256SUMS.txt REVISED_COMMIT.txt; do [ -e "$A/$g" ] && cp "$A/$g" "$dest/"; done
    [ -e revision/submission_package.zip ] && mv revision/submission_package.zip "$dest/submission_package.zip"
done

T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
git archive "$H" revision LICENSE LICENSE-CC-BY-4.0.txt LICENSING.md | tar -x -C "$T"
R=$T/revision
rm -rf "$P" revision/submission_package.zip
mkdir -p "$P"/{environment,manifests/download_logs,logs,plan,claims_ledger,figures,source_data}

cp "$T"/{LICENSE,LICENSE-CC-BY-4.0.txt,LICENSING.md} "$P"/
cp "$R"/environment/{pip_freeze_2026-10-28.txt,pip_freeze_venv_rsf_2026-10-28.txt,R_sessionInfo_2026-10-28.txt,R_packages_2026-10-28.csv,os_info.txt,python_version.txt,requirements_analysis.txt,requirements_rsf.txt} "$P"/environment/
cp "$R"/environment/r00_*.log "$P"/manifests/download_logs/
cp -r "$R"/manifests/{tcga_samples.csv,external_cohort_samples.csv,gdc_queries,DOWNLOAD_LOG.csv,download_dates.md,raw_inputs_manifest.csv,raw_inputs_expected_md5.csv,provenance.md} "$P"/manifests/
cp -r "$R"/results "$R"/tables "$R"/verification "$R"/literature "$P"/
cp "$R"/figures/Fig[1-6].tif "$R"/figures/S*_Fig.tif "$P"/figures/
cp "$R"/source_data/Fig[1-6].csv "$R"/source_data/S*_Fig.csv "$P"/source_data/
cp "$R"/{figure_script_map.csv,methods_parameters.json,methods_concordance.csv,DEVIATIONS.md,editor_checklist.md,numbers_for_manuscript.md,claims_not_reproduced.md} "$P"/
sed "s/{ARCHIVE_COMMIT}/$H/" "$R"/data_availability_draft.md > "$P"/data_availability_draft.md
cp "$R"/logs/{run_analysis_2026-11-01_run7.log,run_analysis_start_run7.txt,r10_fresh_download_2026-11-01.log} "$P"/logs/
cp "$R"/{PREREGISTRATION.md,PREREGISTRATION.sha256,PREREGISTRATION_REGISTRATION.md} "$P"/plan/
cp "$R"/tables/S9_Table.csv "$R"/results/claims_trace_v1.csv "$R"/results/claims_trace_v2.csv "$P"/claims_ledger/
(cd "$P"/tables && zip -q -X -9 S3_Table.zip S3_Table.csv)

{
    printf "Analysis version: %s (%s; code of the final run, run 7)\n" "$ANALYSIS" "$(git log -1 --format=%cI "$ANALYSIS")"
    printf "Archive commit: %s (%s)\n" "$H" "$(git log -1 --format=%cI "$H")"
    printf "Files changed between them (licence files, documentation, archive build script; no analysis code or output):\n"
    if [ -n "$changed" ]; then echo "$changed" | sed 's/^/  /'; else echo "  none"; fi
    printf "Repository: https://github.com/reachdky/GPR153-GBM-Analysis\n"
    printf "Original submission tag: v1-submitted -> %s\n" "$(git rev-parse 'v1-submitted^{commit}')"
    printf "Zenodo: DOI 10.5281/zenodo.23872434 (reserved; record not published when this archive was built)\n"
    printf "Package built: %s from git archive of the archive commit\n" "$(date -u +%FT%TZ)"
} > "$P"/REVISED_COMMIT.txt

n() { find "$P/$1" -type f | wc -l; }
cat > "$P"/PACKAGE_CONTENTS.md <<EOF
# Submission package: PONE-D-26-34494 revised analysis

Built from the archive commit \`$H\` with \`git archive\` (\`REVISED_COMMIT.txt\`). The analysis version is commit
\`$ANALYSIS\` (final run 7); the archive commit adds licence files, documentation and the archive build script only. Every file is a copy of
a file in the archive commit (under \`revision/\`, or the licence files at the repository root), except
\`data_availability_draft.md\` (archive commit filled in), \`tables/S3_Table.zip\` (zip of \`S3_Table.csv\`, for
upload under the 10 MB recommendation for supporting files) and the files written at packaging
(\`REVISED_COMMIT.txt\`, this file, \`SHA256SUMS.txt\`). Paths in the documents that start with \`revision/\` refer
to the repository; here the files sit at the path without that prefix (logs in \`logs/\`, plan in \`plan/\`).

Licences: MIT for author-created software (here only \`verification/scripts/\`); CC BY 4.0 for the authors' figures, result
tables and documentation; third-party materials keep their original terms (\`LICENSING.md\`).

| Folder or file | Contents | Files |
|---|---|---|
| \`LICENSE\`, \`LICENSE-CC-BY-4.0.txt\`, \`LICENSING.md\` | Licence texts and scope, including third-party terms | 3 |
| \`environment/\` | pip freezes of both Python environments, requirement files, R sessionInfo and package list, OS and Python version | $(n environment) |
| \`manifests/\` | TCGA and external-cohort manifests (S1, S2 Table), GDC query JSON, download log with URLs, dates, md5; raw-input checksum records; \`provenance.md\`; download logs | $(n manifests) |
| \`results/\` | All final result CSVs, untruncated (incl. \`screen/screen_gene_level.csv\`, tests, independent check) | $(n results) |
| \`tables/\` | Tables 1-6, S1-S9 Table (CSV); \`S3_Table.zip\` | $(n tables) |
| \`figures/\` | Fig1-Fig6.tif, S1_Fig-S8_Fig.tif: TIFF, LZW, 300 dpi, RGB, width 789-2250 px, height at most 2625 px, Arial 8-12 pt | $(n figures) |
| \`source_data/\` | One CSV per figure (panel column) | $(n source_data) |
| \`claims_ledger/\` | S9 Table and the claims traces | $(n claims_ledger) |
| \`verification/\` | Final-run identification, run comparisons, regeneration check, md5 checks of all inputs, audit against the decision letter, scripts used | $(n verification) |
| \`literature/\` | PubMed and Europe PMC searches (exact queries, dates, counts, records), citation audit, summary | $(n literature) |
| \`logs/\` | Final run log and start/end times; fresh-download log | $(n logs) |
| \`plan/\` | Registered analysis plan v1.5, SHA-256, registration record (https://osf.io/qpd2f/) | $(n plan) |
| \`figure_script_map.csv\`, \`methods_parameters.json\`, \`methods_concordance.csv\`, \`DEVIATIONS.md\` | As in the repository | 4 |
| \`editor_checklist.md\`, \`numbers_for_manuscript.md\`, \`claims_not_reproduced.md\`, \`data_availability_draft.md\` | Manuscript-support documents | 4 |

## Omitted or unavailable

| Item | Reason |
|---|---|
| \`FigS1_cohort_flow.png/.csv\` | Written by \`pipeline/r01_tcga_cohort.py\`; not a numbered figure (the cohort flow is Fig 1). In the repository |
| Figures 7 and 8 of the submitted manuscript | No counterpart: the revised set is Fig 1-6 and S1-S8 Fig |
| Raw downloads (\`data/\`) | Open access; rebuilt by \`r00\` and \`r10\` with md5 checks (\`manifests/\`) |
| Intermediate matrices (\`revision/work/\`, about 1.2 GB) | Derived by the loader scripts; not distributed |
| Earlier runs (\`revision/archive/\`, runs 1-6) | Superseded; comparisons in \`verification/\` |
| Earlier environment snapshots and R build logs | Not the final-run environment |
| Arial font files | Licensed fonts, used locally to render the figures; not distributed |
| \`revision/ai_tool_log.md\` | Excluded as instructed |
| Code | In the source archive and the GitHub repository, not in this package |
| Manuscript and submission letters | Prepared separately for the journal submission; not part of this package |

## Figure format

Checked by \`pipeline/plotstyle.py\` at creation: dimensions, RGB, 300 dpi, LZW, 8-12 pt, and every text element
rendered from an Arial font file (recorded in each TIFF). NAAS has not been run.
EOF

(cd "$P" && find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS.txt)
(cd revision && zip -qr -X submission_package.zip submission_package)

git archive --format=tar.gz --prefix="GPR153-GBM-Analysis-$H/" -o "$A/GPR153-GBM-Analysis-$H-source.tar.gz" "$H"
cp revision/submission_package.zip "$A/GPR153-GBM-Analysis-$H-submission-package.zip"
cp "$P"/REVISED_COMMIT.txt "$A"/
(cd "$A" && sha256sum "GPR153-GBM-Analysis-$H-source.tar.gz" "GPR153-GBM-Analysis-$H-submission-package.zip" > SHA256SUMS.txt && cat SHA256SUMS.txt)
