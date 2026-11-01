# Literature search protocol

Status: protocol only. No search has been run. Fill in each field when the search is run;
do not edit a completed row afterwards (add a new row with the new date instead).

## 1. PubMed searches

Database: PubMed (https://pubmed.ncbi.nlm.nih.gov/). Fields: all fields, no filters, no date limit.
Record the exact string as entered, the search date (UTC) and the result count.

| # | Search string | Search date (UTC) | Result count | Records screened | Records relevant | Notes |
|---|---|---|---|---|---|---|
| 1 | `GPR153 AND glioblastoma` | | | | | |
| 2 | `GPR153 AND glioma` | | | | | |
| 3 | `GPR153 AND cancer` | | | | | |
| 4 | `GPR153 AND MGMT` | | | | | |
| 5 | `GPR153 AND (GLI OR GLI1)` | | | | | |
| 6 | `GPR153 AND WNT` | | | | | |
| 7 | `GPR153 AND temozolomide` | | | | | |
| 8 | `GPR153 AND (leukemia OR AML)` | | | | | |
| 9 | `orphan GPCR AND glioblastoma` | | | | | |

Export the PMID list for each search to `revision/literature/pubmed_<#>_<date>.txt`.

## 2. Modulator searches (ligands, agonists, antagonists)

| Step | Resource | Action | Date (UTC) | Release / version | Result |
|---|---|---|---|---|---|
| 1 | IUPHAR/BPS Guide to Pharmacology (https://www.guidetopharmacology.org) | Search target "GPR153"; record target ID, class, and whether it is listed as orphan | | | |
| 2 | IUPHAR/BPS Guide to Pharmacology | Open the target's ligand tab; record every listed ligand with type (agonist, antagonist, other) and reference | | | |
| 3 | ChEMBL (https://www.ebi.ac.uk/chembl) | Search target "GPR153" (human); record the ChEMBL target ID | | | |
| 4 | ChEMBL | Record the number of bioactivity records and compounds for the target, with assay types | | | |
| 5 | Both | Record "no modulator listed" explicitly if none is found | | | |

## 3. Citation audit

One row per sentence in the manuscript (including the response letter) that cites a reference.
"Supports" is judged from the full text of the cited source, not the abstract alone.

| # | Manuscript location (section, paragraph) | Sentence (verbatim) | Reference cited | Full text checked (yes/no) | Supports the sentence (yes/no/partly) | Action (keep / reword / replace / remove) |
|---|---|---|---|---|---|---|
| 1 | | | | | | |
