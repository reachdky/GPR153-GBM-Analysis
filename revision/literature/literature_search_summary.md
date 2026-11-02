# Literature search and citation audit: evidence record

Evidence collected for the rewrite; not manuscript text. Searches were run on 2026-11-01 (UTC). Scripts:
`pipeline/x05_pubmed_search.py` (PubMed), `pipeline/x06_citation_audit.py` (citation audit). Raw records:
`pubmed_searches.csv` (query, run time, count, PMIDs, PubMed query translation), `pubmed_records.csv` (PMID, PubMed
link, title, journal, year, DOI and DOI link, term flags), `pubmed_raw/` (esearch responses: PMIDs and query
translations), `europepmc_fulltext_searches.csv`, `citation_audit.csv`.

## 1. PubMed searches (NCBI E-utilities esearch, all fields, no date limit; run 2026-11-01 18:04:31-18:04:38 UTC)

| # | Query (exact) | Purpose | Records |
|---|---|---|---|
| 1 | GPR153 AND glioblastoma | decision letter Section 10 | 0 |
| 2 | GPR153 AND glioma | decision letter Section 10 | 0 |
| 3 | GPR153 AND cancer | decision letter Section 10 | 3 |
| 4 | GPR153 AND MGMT | decision letter Section 10 | 0 |
| 5 | GPR153 AND (GLI OR GLI1) | decision letter Section 10 | 0 |
| 6 | GPR153 AND WNT | decision letter Section 10 | 0 |
| 7 | GPR153 AND temozolomide | decision letter Section 10 | 0 |
| 8 | GPR153 AND (leukemia OR AML) | decision letter Section 10 | 3 |
| 9 | orphan GPCR AND glioblastoma | decision letter Section 10 | 5 |
| 10 | GPR153 | every GPR153 record (screening set) | 18 |
| 11 | GPR153 AND (agonist OR antagonist OR ligand OR inhibitor OR modulator OR inverse agonist) | GPR153 modulators | 8 |
| 12 | GPR153 AND (small molecule OR compound OR drug OR screening) | GPR153 modulators | 3 |
| 13 | GPR153 AND (cAMP OR signaling OR signalling) | GPR153 signaling | 7 |

23 unique records (`pubmed_records.csv`). PubMed searches titles, abstracts and indexing terms only; a GPR153
mention in the full text alone is not found (see Section 2).

## 2. Europe PMC full-text searches (run 2026-11-01 18:06:58 UTC)

| Query (exact) | Records |
|---|---|
| "GPR153" AND (glioblastoma OR glioma) | 24 |
| "GPR153" AND MGMT | 8 |
| "GPR153" AND (agonist OR antagonist OR ligand) | 109 |

These match the terms anywhere in the full text (including supplementary gene lists and pharmacology guides), so
most hits mention GPR153 only in passing. The glioma/glioblastoma set includes Fève et al. 2014 (PMID 24662753,
the study the editor names), GPCR expression in medulloblastoma subgroups (PMID 24252460), and an SHH/GLI target
study in medulloblastoma (PMID 42298590). None of the 8 MGMT full-text hits is a glioma study (titles in
`europepmc_fulltext_searches.csv`). The full texts were not read; relevance beyond the title needs author review.

## 3. Relevant findings from the GPR153 records (from titles and abstracts as retrieved)

| PMID | Year | Finding relevant to the revision |
|---|---|---|
| 24662753 | 2014 | Fève et al., PLoS One: GPCR repertoire in glioblastoma stem-like cells (named by the editor). Not returned by PubMed "GPR153" searches (GPR153 not in the indexed abstract); found by DOI and by full-text search |
| 35021931 | 2022 | Expression of GPR21, GPR39, GPR135 and GPR153 in several tumour tissues, including astrocytoma; the abstract reports results for GPR21, GPR39 and GPR135 only. The full text is needed to see the GPR153 result in astrocytoma |
| 40623996 | 2025 | Shao et al., Nat Commun: GPR153 negatively regulates cellular cAMP; SMC-specific deficiency reduces neointima formation; YAP/TAZ and NF-kB (named by the editor) |
| 42095237 | 2026 | Yang et al., Haematologica: fibrinogen upregulates miR-486-5p, which targets GPR153, in AML (named by the editor); with a published comment (42719965) and reply (42817842) |
| 21981325 | 2011 | Sreedharan et al., FEBS J: GPR153 shares origin with GPR162; high expression in thalamus, cerebellum, arcuate nucleus (cited in the submission as ref 9) |
| 15777626 | 2005 | Identification of nine new human rhodopsin-family GPCRs, including GPR153 |
| 42269381 | 2026 | IUPHAR review of class A orphan GPCRs and the criteria for deorphanization (status of GPR153: see full text) |
| 40886928 | 2025 | S100A9 regulates GPR153 and Kcnk16 in DRG neurons (neuropathic pain) |
| 32643557 | 2020 | GPR153 expression altered in kidney of diabetic rats |
| 28960688, 32249236 | 2018, 2020 | GPR153 in bovine anterior pituitary gonadotrophs |

**Modulators.** No retrieved record reports an agonist, antagonist, inverse agonist or other ligand for GPR153;
the records found by searches 11-12 describe GPR153 as an orphan receptor. The editor's statement that GPR153 is
a class A orphan without pharmacology (IUPHAR/BPS Guide to Pharmacology) is consistent with these records; the
Guide to Pharmacology page itself was not checked.

**GBM, MGMT, GLI, WNT, temozolomide.** No PubMed record links GPR153 with glioblastoma, glioma, MGMT, GLI, WNT or
temozolomide in its title or abstract.

## 4. Recommended references (decision letter) verified

All 8 exist with the DOI given in the letter; first author, year, volume and first page match Crossref
(`citation_audit.csv`, rows R1-R8). PMIDs: R1 24662753, R2 20460306, R3 23446420, R4 40623996, R5 42095237
(Haematologica 2026, volume and pages not yet assigned in Crossref), R6 34185076, R7 33293629, R8 40827022.

## 5. Audit of the submitted manuscript's references (16)

| Ref | Result | Detail |
|---|---|---|
| 1-4, 6, 9, 12, 16 | Match | Author, title, journal, year, volume, first page and DOI agree with Crossref |
| 5 (Bady 2012) | Title truncated | Registered title continues "..., yielding a prediction model for comparisons across datasets, tumor grades, and CIMP-status" |
| 7 (Touat) | Year wrong | Nature 2020;580:517-523 (cited as 2019) |
| 8 (Wu, onco-GPCRome) | DOI and details wrong | Cited DOI 10.1016/j.jbc.2021.100601 resolves to a different article (Omata et al., ADAR1). The paper is Wu V et al., J Biol Chem 2019;294(29):11062-11086, doi:10.1074/jbc.REV119.005601, PMID 31171722 |
| 10 (Human Protein Atlas) | URL wrong gene | ENSG00000170925 resolves to TEX13B (request 2026-11-01). GPR153 is ENSG00000158292 |
| 11 (Atwood 2013) | Author list wrong | Fifth author cited as "Bhatt RN"; PubMed lists Atwood SX, Li M, Lee A, Tang JY, Oro AE (PMID 23446420) |
| 13 (Altaner 2008) | Match | No DOI cited; PubMed PMID 18665745, Neoplasma 2008;55(5):369-74 (title search by exact phrase failed; found by author, journal and year) |
| 14 (Minniti 2009) | Match | No DOI cited; PMID 20044633, Anticancer Res 2009;29(12):5171-84 |
| 15 (Wickström) | Journal, DOI and title wrong | Cited as Front Oncol 2015;5:108, doi:10.3389/fonc.2015.00108, which resolves to a different article (Sun et al., histone acetyltransferases). The paper is Wickström M et al., "Wnt/β-catenin pathway regulates MGMT gene expression in cancer and inhibition of Wnt signalling prevents chemoresistance", Nat Commun 2015;6:8904, doi:10.1038/ncomms9904, PMID 26603103 |

Only the first author was checked automatically; the full author list was checked by hand for ref 11 only.

## 6. Changes to the distributed records (2026-11-02)

- The `abstract` column was removed from `pubmed_records.csv` because PubMed abstracts may be protected by the
  publishers' copyright. The file keeps identifiers, bibliographic metadata, PubMed and DOI links, and the term
  flags. The term flags were computed from titles and abstracts when the search was run; rerunning
  `pipeline/x05_pubmed_search.py` retrieves the records again, including abstracts, from PubMed. Earlier repository
  commits still contain the column; history was not rewritten.
- Class A orphan status of GPR153 was checked on 2026-11-02 in the full text of the Concise Guide to PHARMACOLOGY
  2023/24: G protein-coupled receptors (Europe PMC PMC13324819): GPR153 is listed among the class A orphans with no
  ligand listed. The IUPHAR/BPS Guide to Pharmacology website requires login and was not consulted.
- Lee et al. 2010 (PNAS 107:9736, PMID 20460306), suggested in the decision letter as the primary source for Gpr153
  as a Gli1 target: the abstract does not name Gpr153, and the full text and supplementary material could not be
  retrieved (PMC lists the article as not open access; publisher and repository requests returned access
  challenges, which were not bypassed). The claim could not be verified and is not made in the revised manuscript.
- Atwood et al. 2013 (Nature 494:484), cited in the submission for GPR153 as a GLI1 target: its main text (PMC3761364)
  does not mention GPR153; it is not cited in the revision.
