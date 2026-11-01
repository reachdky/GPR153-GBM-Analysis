#!/usr/bin/env python3
"""
s00_make_sim_data.py
Simulated inputs for testing the analysis scripts end to end without reading any real expression or
methylation value (GBM_MODE=sim). Writes the same files, in the same formats, as the loaders
a10_tcga_matrix.py, a19_methylation_inputs.py, a29_external_inputs.py and a59_de_inputs.py, under
revision/sim/work/. Expression values are random draws from a latent-factor model (seed 42).

Clinical tables are copied from the real manifests (sample IDs, cohort flags, survival, covariates) so that the
simulated runs have the real sample sizes; no expression or methylation value is copied. Gene identifiers come
from the GENCODE v23 probe map (identifiers and symbols only).

Outputs (revision/sim/)
  manifests/tcga_samples.csv, manifests/external_cohort_samples.csv  (copies of the real manifests)
  work/tcga_expr.npz, work/methylation_inputs.csv, work/ext_*.npz, work/ext_*_clinical.csv,
  work/de_toil.npz, work/de_toil_samples.csv
  results/cohort/ (copy of the real cohort results, clinical only)

Run from the repository root:
    python pipeline/sim/s00_make_sim_data.py
"""

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import os  # noqa: E402

os.environ["GBM_MODE"] = "sim"
import common  # noqa: E402

rng = np.random.default_rng(common.SEED)
REAL_MAN = Path("revision/manifests")


def genes_v23():
    pm = pd.read_csv(common.DL / "toil/gencode.v23.annotation.gene.probemap", sep="\t", usecols=[0, 1])
    pm.columns = ["id", "gene"]
    pm["gid"] = pm["id"].str.replace(r"\.\d+(_PAR_Y)?$", "", regex=True)
    pm = pm[~pm["id"].str.endswith("_PAR_Y")].drop_duplicates("gid")
    return pm


def latent_matrix(n_genes, factors, special, batch=None, on_frac=0.62):
    """log2(TPM + 1) = mu + loadings . factors (+ batch) + noise; `special` maps row -> loadings on factor 0."""
    k, n = factors.shape
    on = rng.random(n_genes) < on_frac
    mu = np.where(on, np.clip(rng.normal(3.5, 2.0, n_genes), 0.6, 14), rng.uniform(0, 0.25, n_genes))
    load = rng.normal(0, 0.35, (n_genes, k)) * (rng.random((n_genes, k)) < 0.3)
    for r, l0 in special.items():
        mu[r], on[r] = 4.0, True
        load[r, 0] = l0
    x = mu[:, None] + load @ factors + rng.normal(0, 0.5, (n_genes, n))
    if batch is not None:
        x += (rng.normal(0, 0.3, n_genes) * (rng.random(n_genes) < 0.1))[:, None] * batch[None, :]
    x[~on] = np.clip(x[~on], 0, None) * (rng.random(((~on).sum(), n)) < 0.3)
    return np.clip(x, 0, None)


def tcga(pm):
    m = pd.read_csv(REAL_MAN / "tcga_samples.csv")
    pts = m["patient"].unique()
    F = dict(zip(pts, rng.normal(size=(len(pts), 6))))
    factors = np.column_stack([F[p] for p in m["patient"]])
    batch = (m["batch"] == "GDC2024").to_numpy(float)
    gid = pm["gid"].to_numpy()
    row = {g: i for i, g in enumerate(gid)}
    special = {row[common.MGMT]: 0.7, row[common.GPR153]: 0.45}
    for s in common.NEURONAL + common.FOREST_GENES:
        hit = np.where(pm["gene"].to_numpy() == s)[0]
        if len(hit):
            special.setdefault(hit[0], rng.normal(0, 0.2))
    L = latent_matrix(len(gid), factors, special, batch)
    tpm = np.round(2 ** L - 1, 4)
    counts = rng.poisson(tpm * 25).astype(np.int64)
    fpkm = np.round(tpm * rng.uniform(0.25, 0.35, tpm.shape[1])[None, :], 4)
    common.save_matrix("tcga_expr.npz", tpm=tpm, fpkm=fpkm, counts=counts, gene_id=gid.astype(str),
                       symbol=pm["gene"].fillna("").to_numpy(str), gene_type=np.repeat("simulated", len(gid)),
                       aliquot=m["aliquot"].to_numpy(str))
    print(f"tcga_expr.npz: {tpm.shape}")
    return m


def methylation(m):
    """Betas drawn conditional on the Ceccarelli call where one exists, so that kappa is high."""
    p = m[(m["sample_type"] == "Primary Tumor") & m["methylation_platforms"].notna()].drop_duplicates("patient")
    meth = np.where(p["mgmt_promoter"].eq("methylated"), True,
                    np.where(p["mgmt_promoter"].eq("unmethylated"), False, rng.random(len(p)) < 0.4))
    b = np.where(meth[:, None], rng.beta(6, 3, (len(p), 2)), rng.beta(1, 12, (len(p), 2)))
    b[rng.random(b.shape) < 0.03] = np.nan
    d = pd.DataFrame({"patient": p["patient"].values,
                      "platform_used": np.where(p["methylation_platforms"].str.contains("HM450"), "HM450", "HM27"),
                      "aliquot": p["aliquot"].str[:15].values + "-sim", "file_id": "sim",
                      "platforms_available": p["methylation_platforms"].values,
                      "beta_cg12434587": b[:, 0], "beta_cg12981137": b[:, 1]})
    d.to_csv(common.work_path("methylation_inputs.csv"), index=False)


def ext_matrix(features, n, special, linear):
    f = rng.normal(size=(1, n))
    X = rng.normal(7, 1.5, (len(features), 1)) + 0.4 * rng.normal(size=(len(features), n))
    for name, l0 in special.items():
        X[features.index(name)] = 7 + l0 * f[0] + 0.8 * rng.normal(size=n)
    return 2 ** X if linear else X


def external(pm):
    man = pd.read_csv(REAL_MAN / "external_cohort_samples.csv")
    probes = [f"sim_{i}_at" for i in range(22000)] + common.GPR153_PROBES + [common.MGMT_PROBE]
    sp = {"221902_at": 0.3, "64942_at": 0.1, common.MGMT_PROBE: 0.5}
    for plat, name in (("GPL96", "GSE13041"), ("GPL570", "GSE13041_GPL570")):
        g = man[(man.cohort == "GSE13041_GPL96") & (man.platform == plat)]
        common.save_matrix(f"ext_{name}.npz", X=ext_matrix(probes, len(g), sp, True), feature=np.array(probes),
                           sample=g.sample_id.to_numpy(str))
        g.to_csv(common.work_path(f"ext_{name}_clinical.csv"), index=False)
    g = man[man.cohort == "GSE4290_GBM"]
    common.save_matrix("ext_GSE4290.npz", X=ext_matrix(probes, len(g), sp, True), feature=np.array(probes),
                       sample=g.sample_id.to_numpy(str))
    g.to_csv(common.work_path("ext_GSE4290_clinical.csv"), index=False)
    syms = list(dict.fromkeys(list(pm["gene"].dropna().drop_duplicates().iloc[:23000]) + ["GPR153", "MGMT"]))
    for n in ("693", "325"):
        g = man[(man.cohort == f"CGGA_{n}") & man.included].copy()
        g["sex"] = g["sex"].str.strip()
        X = np.clip(ext_matrix(syms, len(g), {"GPR153": 0.5, "MGMT": 0.6}, False) - 5, 0, None)
        common.save_matrix(f"ext_CGGA_{n}.npz", X=2 ** X - 1, feature=np.array(syms), sample=g.sample_id.to_numpy(str))
        g.to_csv(common.work_path(f"ext_CGGA_{n}_clinical.csv"), index=False)
    n = 220
    t = rng.exponential(450, n)
    t[rng.choice(n, 10, replace=False)] = np.nan
    lo = rng.choice(np.arange(20, 85, 5), n)
    clin = pd.DataFrame({"sample_id": [f"sim{i:05d}_U133P2" for i in range(n)], "subject": [f"S{i}" for i in range(n)],
                         "os_days": t, "os_event": (rng.random(n) < 0.89).astype(float), "age_range": [f"{a}-{a + 4}" for a in lo],
                         "age": lo + 2.0, "sex": rng.choice(["Male", "Female"], n)})
    common.save_matrix("ext_Rembrandt.npz", X=ext_matrix(probes, n, sp, False), feature=np.array(probes),
                       sample=clin.sample_id.to_numpy(str))
    clin.to_csv(common.work_path("ext_Rembrandt_clinical.csv"), index=False)
    n = 191
    t = rng.exponential(500, n)
    t[:3] = [-8, 0, 0]
    clin = pd.DataFrame({"sample_id": [f"C3L-sim{i:03d}" for i in range(n)], "file_id": "sim", "os_days": t,
                         "os_event": (rng.random(n) < 0.74).astype(float), "vital_status": "sim",
                         "age": rng.normal(60, 12, n), "sex": rng.choice(["Male", "Female"], n),
                         "idh_status": np.where(rng.random(n) < 0.92, "wildtype", "mutant")})
    gid = list(pm["gid"].iloc[:30000]) + [common.GPR153, common.MGMT]
    gid = list(dict.fromkeys(gid))
    X = np.clip(ext_matrix(gid, n, {common.GPR153: 0.5, common.MGMT: 0.6}, False) - 5, 0, None)
    common.save_matrix("ext_CPTAC.npz", X=2 ** X - 1, feature=np.array(gid), sample=clin.sample_id.to_numpy(str))
    clin.to_csv(common.work_path("ext_CPTAC_clinical.csv"), index=False)


def de(pm):
    """Toil-format counts: 144 tumors (primary-cohort BI barcodes) and 206 normals from 151 donors."""
    m = pd.read_csv(REAL_MAN / "tcga_samples.csv")
    tum = m[m.in_screen_primary & m.selected_aliquot & (m.batch == "BI")]["aliquot"].str[:15].iloc[:144].tolist()
    donors = [f"GTEX-SIM{i:03d}" for i in range(151)]
    nor = [(f"{d}-0011-R10a-SM-SIM", "Brain - Cortex", d) for d in donors[:105]]
    nor += [(f"{d}-0011-R10b-SM-SIM", "Brain - Frontal Cortex (Ba9)", d) for d in donors[54:155]]
    samples = pd.concat([
        pd.DataFrame({"sample": tum, "group": "tumor", "tissue": "TCGA-GBM primary", "donor": [s[:12] for s in tum],
                      "normal_one_per_donor": False}),
        pd.DataFrame({"sample": [x[0] for x in nor], "group": "normal", "tissue": [x[1] for x in nor],
                      "donor": [x[2] for x in nor]})])
    first = samples[samples.group == "normal"].assign(pref=lambda d: (d.tissue != "Brain - Cortex").astype(int)).sort_values(
        ["donor", "pref", "sample"]).groupby("donor").head(1)["sample"]
    samples["normal_one_per_donor"] = samples["sample"].isin(first)
    gid = pm["gid"].to_numpy()
    base = np.exp(rng.normal(3, 2.5, len(gid)))
    lfc = rng.normal(0, 0.4, len(gid)) * (rng.random(len(gid)) < 0.2)
    is_t = (samples.group == "tumor").to_numpy()
    mu = base[:, None] * np.where(is_t[None, :], 2 ** lfc[:, None], 1) * rng.uniform(0.6, 1.6, len(samples))[None, :]
    counts = rng.negative_binomial(5, 5 / (5 + mu)).astype(np.int64)
    common.save_matrix("de_toil.npz", counts=counts, gene_id=gid.astype(str), sample=samples["sample"].to_numpy(str))
    samples.to_csv(common.work_path("de_toil_samples.csv"), index=False)


def main():
    for d in (common.MAN, common.WORK, common.RES):
        d.mkdir(parents=True, exist_ok=True)
    for f in ("tcga_samples.csv", "external_cohort_samples.csv"):
        shutil.copy(REAL_MAN / f, common.MAN / f)
    shutil.copytree("revision/results/cohort", common.RES / "cohort", dirs_exist_ok=True)
    pm = genes_v23()
    m = tcga(pm)
    methylation(m)
    external(pm)
    de(pm)


if __name__ == "__main__":
    main()
