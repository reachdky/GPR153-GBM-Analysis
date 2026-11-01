#!/usr/bin/env python3
"""
r00_download_external.py
Download the open-access external resources used by the revised pipeline and
record URL, retrieval time (UTC), size and md5 for each file in
data/downloads/DOWNLOAD_LOG.csv.

Files already present are not downloaded again; their log entry is kept.

Run from the repository root:
    python pipeline/r00_download_external.py [step ...]
Steps: expression_meta maf clinical cbioportal static cna methylation cptac depmap
cellxgene (default: all, in that order).
"""

import csv
import fcntl
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

OUT = Path("data/downloads")
LOG = OUT / "DOWNLOAD_LOG.csv"
GDC = "https://api.gdc.cancer.gov"
CBIO = "https://www.cbioportal.org/api"
CGGA = "https://www.cgga.org.cn/download?file=download/{date}/{name}&type={type}&time={date}"
TCGA_EXPR = Path("data/raw/tcga/expression/TCGA-GBM/Transcriptome_Profiling/Gene_Expression_Quantification")
TOIL = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/"
GDSC = "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/"
GEO_SER = "https://ftp.ncbi.nlm.nih.gov/geo/series/{p}/{g}/"
GEO_PLATFORM = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={p}&targ=self&form=text&view=data"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/{}/fullTextXML"
MSIGDB = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.2.Hs/"
ENRICHR = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="
PMC_BLOB = "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/dd12/3444709/"
GDC_DATA = GDC + "/data/"

STATIC = [
    # (subdir, filename, url, description)
    ("toil", "TcgaTargetGtex_gene_expected_count.gz", TOIL + "TcgaTargetGtex_gene_expected_count.gz",
     "UCSC Toil RSEM expected counts, log2(x+1), GENCODE v23"),
    ("toil", "TcgaTargetGTEX_phenotype.txt.gz", TOIL + "TcgaTargetGTEX_phenotype.txt.gz",
     "UCSC Toil sample phenotype (study, sample type, GTEx tissue)"),
    ("toil", "gencode.v23.annotation.gene.probemap", TOIL + "probeMap/gencode.v23.annotation.gene.probemap",
     "UCSC Toil gene ID map, GENCODE v23"),
    ("cgga", "CGGA.mRNAseq_693.RSEM-genes.20200506.txt.zip",
     CGGA.format(date="20200506", name="CGGA.mRNAseq_693.RSEM-genes.20200506.txt.zip", type="mRNAseq_693"),
     "CGGA mRNAseq_693 RSEM FPKM"),
    ("cgga", "CGGA.mRNAseq_693_clinical.20200506.txt.zip",
     CGGA.format(date="20200506", name="CGGA.mRNAseq_693_clinical.20200506.txt.zip", type="mRNAseq_693_clinical"),
     "CGGA mRNAseq_693 clinical"),
    ("cgga", "CGGA.mRNAseq_325.RSEM-genes.20200506.txt.zip",
     CGGA.format(date="20200506", name="CGGA.mRNAseq_325.RSEM-genes.20200506.txt.zip", type="mRNAseq_325"),
     "CGGA mRNAseq_325 RSEM FPKM"),
    ("cgga", "CGGA.mRNAseq_325_clinical.20200506.txt.zip",
     CGGA.format(date="20200506", name="CGGA.mRNAseq_325_clinical.20200506.txt.zip", type="mRNAseq_325_clinical"),
     "CGGA mRNAseq_325 clinical"),
    ("estimate", "tidyestimate_1.1.1.tar.gz", "https://cran.r-project.org/src/contrib/tidyestimate_1.1.1.tar.gz",
     "ESTIMATE stromal/immune gene sets (tidyestimate 1.1.1, CRAN)"),
    # GDSC release 8.5 (drug response, compound list, cell line annotation)
    *[("gdsc", n, GDSC + n, d) for n, d in [
        ("screened_compounds_rel_8.5.csv", "GDSC release 8.5 screened compounds"),
        ("GDSC1_fitted_dose_response_27Oct23.xlsx", "GDSC1 fitted dose response (LN_IC50, AUC), release 8.5"),
        ("GDSC2_fitted_dose_response_27Oct23.xlsx", "GDSC2 fitted dose response (LN_IC50, AUC), release 8.5"),
        ("Cell_Lines_Details.xlsx", "GDSC cell line details, release 8.5")]],
    # GEO: single-cell sample annotation and survival-replication candidates
    ("geo_sc", "GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx",
     GEO_SER.format(p="GSE131nnn", g="GSE131928") + "suppl/GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx",
     "GSE131928 (Neftel 2019) cell-to-tumor annotation"),
    ("geo_replication", "GSE16011_series_matrix.txt.gz",
     GEO_SER.format(p="GSE16nnn", g="GSE16011") + "matrix/GSE16011_series_matrix.txt.gz", "GSE16011 series matrix"),
    ("geo_replication", "GSE68848_series_matrix.txt.gz",
     GEO_SER.format(p="GSE68nnn", g="GSE68848") + "matrix/GSE68848_series_matrix.txt.gz",
     "GSE68848 (Rembrandt subset) series matrix"),
    *[("geo_replication", n, GEO_SER.format(p="GSE108nnn", g="GSE108474") + "suppl/" + n, d) for n, d in [
        ("GSE108474_REMBRANDT_clinical.data.txt.gz", "Rembrandt clinical data"),
        ("GSE108474_REMBRANDT_biospecimen_mapping_GEO.txt.gz", "Rembrandt subject to GSM mapping"),
        ("GSE108474_REMBRANDT_data.dictionary.xlsx", "Rembrandt clinical data dictionary"),
        ("GSE108474_REMBRANDT_GeneExpression.txt.gz", "Rembrandt GPL570 expression")]],
    # GEO platform tables (annotation only)
    *[("platforms", f"{p}_table.txt", GEO_PLATFORM.format(p=p), d) for p, d in [
        ("GPL8490", "Illumina HumanMethylation27 annotation"), ("GPL13534", "Illumina HumanMethylation450 annotation"),
        ("GPL8542", "GSE16011 platform (Brainarray HGU133Plus2_Hs_ENTREZG v10) annotation"),
        ("GPL570", "Affymetrix HG-U133 Plus 2.0 annotation")]],
    # Gene sets: MSigDB 2026.2.Hs (Hallmark, Reactome) and Enrichr libraries
    *[("genesets", n, MSIGDB + n, d) for n, d in [
        ("h.all.v2026.2.Hs.symbols.gmt", "MSigDB 2026.2.Hs Hallmark gene sets"),
        ("c2.cp.reactome.v2026.2.Hs.symbols.gmt", "MSigDB 2026.2.Hs Reactome gene sets")]],
    *[("genesets", f"{lib}.gmt", ENRICHR + lib, f"Enrichr library {lib}")
      for lib in ["Reactome_2022", "KEGG_2021_Human", "GO_Biological_Process_2023"]],
    # MGMT-STP27 source publication (Bady et al. 2012) and erratum, open access
    ("literature", "PMC3444709.xml", EPMC.format("PMC3444709"), "Bady 2012 Acta Neuropathol 124:547, full text"),
    ("literature", "PMC4079675.xml", EPMC.format("PMC4079675"), "Bady 2012 erratum (Acta Neuropathol 2013)"),
    ("literature", "PMC3444709_M1.gif", PMC_BLOB + "ddbdbfe5067e/M1.gif", "Bady 2012 M-value equation"),
    ("literature", "PMC3444709_M2.gif", PMC_BLOB + "762e08f1f6ff/M2.gif", "Bady 2012 MGMT-STP27 model equation"),
    # TCGA publication-page files (gdc.cancer.gov/about-data/publications): extent of resection and MGMT call method
    ("tcga_pub", "TCGA-CDR-SupplementalTableS1.xlsx", GDC_DATA + "1b5f413e-a8d1-4d10-92eb-7c4ae739ed81",
     "TCGA Pan-Cancer Clinical Data Resource (Liu 2018, Cell 173:400), Table S1"),
    ("tcga_pub", "clinical_PANCAN_patient_with_followup.tsv", GDC_DATA + "0fc78496-818b-4896-bd83-52db1f533c5c",
     "PanCanAtlas merged BCR clinical table, patient with follow-up"),
    ("tcga_pub", "lgggbm_2016_AWG_clinical_supplementary_info", GDC_DATA + "a0886693-4c6e-40de-9cb1-8bfcb4bb30a7",
     "Ceccarelli 2016 publication page: AWG Clinical Supplementary Info"),
    ("tcga_pub", "lgggbm_2016_GBM_clinical_data", GDC_DATA + "e988cf7c-135c-4770-8b0b-d44160a03ac5",
     "Ceccarelli 2016 publication page: GBM Clinical data"),
    ("tcga_pub", "TCGA_GBM_dnameth_call_categories.txt", GDC_DATA + "a4c616e1-711c-4230-9da1-73135245a42d",
     "Brennan 2013 publication page: DNA methylation call categories"),
]

# DepMap 24Q4 Public (figshare 10.25452/figshare.plus.27993248.v1); md5 from figshare
DEPMAP_ARTICLE = 27993248
DEPMAP_FILES = ["Model.csv", "OmicsProfiles.csv", "README.txt", "OmicsExpressionProteinCodingGenesTPMLogp1.csv"]
CELLXGENE = "https://api.cellxgene.cziscience.com/curation/v1"
GBMAP_COLLECTION = "999f2a15-3d7e-440b-96ae-2c806799c08c"
GBMAP_CORE = "c888b684-6c51-431f-972a-6c963044cef0"

CBIO_STUDIES = ["gbm_tcga_pub2013", "gbm_tcga_pan_can_atlas_2018", "lgggbm_tcga_pub"]


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_log():
    if not LOG.exists():
        return {}
    with open(LOG) as f:
        return {r["path"]: r for r in csv.DictReader(f)}


def save_log(log):
    fields = ["path", "url", "description", "retrieved_utc", "bytes", "md5"]
    with open(LOG, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for k in sorted(log):
            w.writerow(log[k])


def record(log, path, url, desc):
    """Add one entry; merge with the log on disk under a lock so parallel runs do not lose entries."""
    entry = {"path": str(path), "url": url, "description": desc,
             "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "bytes": path.stat().st_size, "md5": md5(path)}
    with open(LOG.with_suffix(".lock"), "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        log.update(load_log())
        log[str(path)] = entry
        save_log(log)


def fetch(url, path, method="GET", payload=None, retries=3):
    tmp = path.with_suffix(path.suffix + ".part")
    for attempt in range(1, retries + 1):
        try:
            with requests.request(method, url, json=payload, stream=True, timeout=300) as r:
                r.raise_for_status()
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
            tmp.rename(path)
            return
        except requests.RequestException as e:
            print(f"  attempt {attempt} failed: {e}")
            time.sleep(10 * attempt)
    sys.exit(f"Download failed: {url}")


def get_static(log):
    for sub, name, url, desc in STATIC:
        path = OUT / sub / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and str(path) in log:
            print(f"present: {path}")
            continue
        print(f"downloading: {path}")
        fetch(url, path)
        record(log, path, url, desc)


def get_cbioportal(log):
    d = OUT / "cbioportal"
    d.mkdir(parents=True, exist_ok=True)
    for study in CBIO_STUDIES:
        for level in ["PATIENT", "SAMPLE"]:
            path = d / f"{study}_{level.lower()}_clinical.json"
            if path.exists() and str(path) in log:
                print(f"present: {path}")
                continue
            url = f"{CBIO}/studies/{study}/clinical-data?clinicalDataType={level}&projection=DETAILED&pageSize=100000"
            print(f"downloading: {path}")
            fetch(url, path)
            record(log, path, url, f"cBioPortal {study} {level.lower()} clinical data")
        path = d / f"{study}_attributes.json"
        if not (path.exists() and str(path) in log):
            url = f"{CBIO}/studies/{study}/clinical-attributes?projection=DETAILED"
            fetch(url, path)
            record(log, path, url, f"cBioPortal {study} clinical attribute definitions")


def get_gdc_maf(log):
    """Open-access masked somatic mutation MAFs for TCGA-GBM (one per tumor aliquot)."""
    d = OUT / "gdc_maf"
    d.mkdir(parents=True, exist_ok=True)
    filt = {"op": "and", "content": [
        {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-GBM"}},
        {"op": "=", "content": {"field": "data_type", "value": "Masked Somatic Mutation"}},
        {"op": "=", "content": {"field": "access", "value": "open"}}]}
    fields = ("file_id,file_name,md5sum,file_size,analysis.workflow_type,data_format,"
              "cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type")
    q = {"filters": filt, "fields": fields, "format": "JSON", "size": "2000"}
    meta_path = d / "maf_files_meta.json"
    r = requests.post(f"{GDC}/files", json=q, timeout=300)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    status = requests.get(f"{GDC}/status", timeout=60).json()
    with open(meta_path, "w") as f:
        json.dump({"query": q, "gdc_status": status,
                   "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                   "hits": hits}, f, indent=1)
    record(log, meta_path, f"{GDC}/files", "GDC query and metadata for TCGA-GBM open masked MAFs")
    for h in hits:
        path = d / h["file_name"]
        if path.exists() and md5(path) == h["md5sum"]:
            continue
        fetch(f"{GDC}/data/{h['file_id']}", path)
        if md5(path) != h["md5sum"]:
            sys.exit(f"md5 mismatch: {path}")
    print(f"MAFs: {len(hits)} files, GDC {status.get('data_release')}")


def get_gdc_expression_meta(log):
    """GDC metadata (aliquot, sample type, center, workflow, md5) for the STAR count files on disk."""
    d = OUT / "gdc_meta"
    d.mkdir(parents=True, exist_ok=True)
    ids = sorted(p.name for p in TCGA_EXPR.iterdir() if p.is_dir())
    fields = ("file_id,file_name,md5sum,file_size,created_datetime,analysis.workflow_type,experimental_strategy,"
              "platform,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type,"
              "cases.samples.portions.analytes.aliquots.submitter_id,"
              "cases.samples.portions.analytes.aliquots.center.short_name")
    q = {"filters": {"op": "in", "content": {"field": "file_id", "value": ids}},
         "fields": fields, "format": "JSON", "size": str(len(ids) + 10)}
    r = requests.post(f"{GDC}/files", json=q, timeout=300)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    path = d / "expression_files_meta.json"
    with open(path, "w") as f:
        json.dump({"query": q, "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                   "gdc_status": requests.get(f"{GDC}/status", timeout=60).json(), "hits": hits}, f, indent=1)
    record(log, path, f"{GDC}/files", "GDC metadata for the TCGA-GBM STAR count files on disk")
    rows = []
    for h in hits:
        s = h["cases"][0]["samples"][0]
        a = s["portions"][0]["analytes"][0]["aliquots"][0]
        local = next((TCGA_EXPR / h["file_id"]).glob("*.tsv"))
        rows.append({"file_id": h["file_id"], "file_name": h["file_name"], "aliquot": a["submitter_id"],
                     "patient": h["cases"][0]["submitter_id"], "sample_type": s["sample_type"],
                     "center": a.get("center", {}).get("short_name", ""), "platform": h.get("platform", ""),
                     "workflow": h["analysis"]["workflow_type"], "gdc_created": h["created_datetime"][:10],
                     "md5_gdc": h["md5sum"], "md5_local": md5(local)})
    csv_path = d / "expression_files_meta.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(sorted(rows, key=lambda x: x["aliquot"]))
    bad = [x["aliquot"] for x in rows if x["md5_gdc"] != x["md5_local"]]
    print(f"Expression files: {len(rows)} in GDC metadata; md5 mismatches: {len(bad)}")


def get_gdc_clinical(log):
    """Case-level GDC clinical and per-file metadata for the TCGA-GBM STAR count files on disk."""
    d = OUT / "gdc_meta"
    d.mkdir(parents=True, exist_ok=True)
    path = d / "tcga_gbm_cases.json"
    filt = {"op": "=", "content": {"field": "project.project_id", "value": "TCGA-GBM"}}
    q = {"filters": filt, "expand": "demographic,diagnoses,diagnoses.treatments,follow_ups",
         "format": "JSON", "size": "2000"}
    r = requests.post(f"{GDC}/cases", json=q, timeout=300)
    r.raise_for_status()
    with open(path, "w") as f:
        json.dump({"query": q, "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                   "gdc_status": requests.get(f"{GDC}/status", timeout=60).json(),
                   "hits": r.json()["data"]["hits"]}, f, indent=1)
    record(log, path, f"{GDC}/cases", "GDC TCGA-GBM case clinical records")


def gdc_files(log, filt, fields, subdir, desc, download=True):
    """Query GDC /files, save query + metadata JSON, and download the files (md5-verified)."""
    d = OUT / subdir
    d.mkdir(parents=True, exist_ok=True)
    q = {"filters": filt, "fields": fields, "format": "JSON", "size": "10000"}
    r = requests.post(f"{GDC}/files", json=q, timeout=300)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    meta = d / "files_meta.json"
    with open(meta, "w") as f:
        json.dump({"query": q, "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                   "gdc_status": requests.get(f"{GDC}/status", timeout=60).json(), "hits": hits}, f, indent=1)
    record(log, meta, f"{GDC}/files", f"GDC query and metadata: {desc}")
    if download:
        for h in hits:
            path = d / h["file_id"] / h["file_name"]
            path.parent.mkdir(exist_ok=True)
            if path.exists() and md5(path) == h["md5sum"]:
                continue
            fetch(f"{GDC}/data/{h['file_id']}", path)
            if md5(path) != h["md5sum"]:
                sys.exit(f"md5 mismatch: {path}")
    print(f"{desc}: {len(hits)} files")
    return hits


def get_gdc_methylation(log):
    """TCGA-GBM DNA methylation beta values (SeSAMe), 27K and 450K, open access."""
    filt = {"op": "and", "content": [
        {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-GBM"}},
        {"op": "=", "content": {"field": "data_type", "value": "Methylation Beta Value"}},
        {"op": "=", "content": {"field": "access", "value": "open"}}]}
    fields = ("file_id,file_name,md5sum,file_size,platform,analysis.workflow_type,cases.submitter_id,"
              "cases.samples.submitter_id,cases.samples.sample_type,"
              "cases.samples.portions.analytes.aliquots.submitter_id")
    gdc_files(log, filt, fields, "gdc_methylation", "TCGA-GBM methylation beta values (27K, 450K)")


def get_cptac(log):
    """CPTAC-3 brain tumors: STAR counts, open masked MAFs, case clinical records."""
    base = [{"op": "=", "content": {"field": "cases.project.project_id", "value": "CPTAC-3"}},
            {"op": "=", "content": {"field": "cases.primary_site", "value": "Brain"}},
            {"op": "=", "content": {"field": "access", "value": "open"}}]
    fields = ("file_id,file_name,md5sum,file_size,data_type,analysis.workflow_type,cases.submitter_id,"
              "cases.samples.submitter_id,cases.samples.sample_type")
    gdc_files(log, {"op": "and", "content": base + [
        {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}}]},
        fields, "cptac/expression", "CPTAC-3 brain STAR counts")
    gdc_files(log, {"op": "and", "content": base + [
        {"op": "=", "content": {"field": "data_type", "value": "Masked Somatic Mutation"}}]},
        fields, "cptac/maf", "CPTAC-3 brain open masked MAFs")
    path = OUT / "cptac" / "cptac_brain_cases.json"
    q = {"filters": {"op": "and", "content": [
            {"op": "=", "content": {"field": "project.project_id", "value": "CPTAC-3"}},
            {"op": "=", "content": {"field": "primary_site", "value": "Brain"}}]},
         "expand": "demographic,diagnoses,diagnoses.treatments,follow_ups", "format": "JSON", "size": "2000"}
    r = requests.post(f"{GDC}/cases", json=q, timeout=300)
    r.raise_for_status()
    with open(path, "w") as f:
        json.dump({"query": q, "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                   "gdc_status": requests.get(f"{GDC}/status", timeout=60).json(),
                   "hits": r.json()["data"]["hits"]}, f, indent=1)
    record(log, path, f"{GDC}/cases", "GDC CPTAC-3 brain case clinical records")


def get_gdc_bcr_xml(log):
    """TCGA-GBM clinical supplements (original BCR XML), open access."""
    filt = {"op": "and", "content": [
        {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-GBM"}},
        {"op": "=", "content": {"field": "data_type", "value": "Clinical Supplement"}},
        {"op": "=", "content": {"field": "data_format", "value": "BCR XML"}}]}
    gdc_files(log, filt, "file_id,file_name,md5sum,file_size,cases.submitter_id", "gdc_bcr_xml",
              "TCGA-GBM clinical supplements (BCR XML)")


def get_cbioportal_cna(log):
    """EGFR (Entrez 1956) GISTIC copy-number calls, TCGA-GBM."""
    d = OUT / "cbioportal"
    for study in ["gbm_tcga_pan_can_atlas_2018", "gbm_tcga_pub2013"]:
        path = d / f"{study}_EGFR_gistic.json"
        url = f"{CBIO}/molecular-profiles/{study}_gistic/molecular-data/fetch?projection=SUMMARY"
        fetch(url, path, method="POST", payload={"entrezGeneIds": [1956], "sampleListId": f"{study}_cna"})
        record(log, path, url, f"cBioPortal {study} EGFR GISTIC calls")


def get_depmap(log):
    d = OUT / "depmap"
    d.mkdir(parents=True, exist_ok=True)
    art = requests.get(f"https://api.figshare.com/v2/articles/{DEPMAP_ARTICLE}", timeout=60).json()
    with open(d / "figshare_article.json", "w") as f:
        json.dump(art, f, indent=1)
    record(log, d / "figshare_article.json", f"https://api.figshare.com/v2/articles/{DEPMAP_ARTICLE}",
           f"DepMap figshare record: {art['title']}, {art.get('doi')}")
    for fi in [x for x in art["files"] if x["name"] in DEPMAP_FILES]:
        path = d / fi["name"]
        if not (path.exists() and str(path) in log):
            fetch(fi["download_url"], path)
            record(log, path, fi["download_url"], f"{art['title']}: {fi['name']}")
        if md5(path) != fi["computed_md5"]:
            sys.exit(f"md5 mismatch: {path}")


def get_cellxgene(log):
    """Core GBmap (Ruiz-Moreno et al.) H5AD from CELLxGENE Discover."""
    d = OUT / "cellxgene"
    d.mkdir(parents=True, exist_ok=True)
    meta = requests.get(f"{CELLXGENE}/collections/{GBMAP_COLLECTION}", timeout=120).json()
    with open(d / "gbmap_collection.json", "w") as f:
        json.dump(meta, f, indent=1)
    record(log, d / "gbmap_collection.json", f"{CELLXGENE}/collections/{GBMAP_COLLECTION}", "GBmap collection metadata")
    ds = next(x for x in meta["datasets"] if x["dataset_id"] == GBMAP_CORE)
    asset = next(a for a in ds["assets"] if a["filetype"] == "H5AD")
    path = d / f"core_GBmap_{ds['dataset_version_id']}.h5ad"
    if not (path.exists() and str(path) in log):
        fetch(asset["url"], path)
        record(log, path, asset["url"], f"Core GBmap H5AD, dataset version {ds['dataset_version_id']}")


STEPS = {"expression_meta": get_gdc_expression_meta, "maf": get_gdc_maf, "clinical": get_gdc_clinical,
         "cbioportal": get_cbioportal, "static": get_static, "cna": get_cbioportal_cna,
         "methylation": get_gdc_methylation, "cptac": get_cptac, "depmap": get_depmap, "cellxgene": get_cellxgene,
         "bcr_xml": get_gdc_bcr_xml}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    log = load_log()
    for name in (sys.argv[1:] or STEPS):
        STEPS[name](log)
    print(f"Done. Log: {LOG}")


if __name__ == "__main__":
    main()
