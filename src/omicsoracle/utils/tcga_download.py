import requests
import json
import os
import tarfile
import shutil
from pathlib import Path
from tqdm import tqdm

GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
GDC_DATA_ENDPOINT  = "https://api.gdc.cancer.gov/data"

CANCER_TYPES = {
    "BRCA": "Breast Invasive Carcinoma",
    "LUAD": "Lung Adenocarcinoma",
    "GBM":  "Glioblastoma Multiforme",
}

def build_query(cancer_type: str, n_files: int = 20) -> dict:
    return {
        "filters": json.dumps({
            "op": "and",
            "content": [
                {"op": "=", "content": {"field": "cases.project.project_id",
                                         "value": f"TCGA-{cancer_type}"}},
                {"op": "=", "content": {"field": "data_type",
                                         "value": "Gene Expression Quantification"}},
                {"op": "=", "content": {"field": "analysis.workflow_type",
                                         "value": "STAR - Counts"}},
                {"op": "=", "content": {"field": "data_format",
                                         "value": "TSV"}},
            ]
        }),
        "fields": "file_id,file_name,cases.submitter_id,cases.case_id,"
                  "cases.demographic.gender,cases.diagnoses.tumor_stage,"
                  "cases.samples.sample_type",
        "format": "JSON",
        "size": n_files,
    }

def fetch_file_manifest(cancer_type: str, n_files: int = 20) -> list[dict]:
    resp = requests.get(GDC_FILES_ENDPOINT, params=build_query(cancer_type, n_files))
    resp.raise_for_status()
    hits = resp.json()["data"]["hits"]
    print(f"  {cancer_type}: found {len(hits)} files")
    return hits

def download_files(file_ids: list[str], out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    archive = out_dir / "download.tar.gz"
    payload = {"ids": file_ids}
    with requests.post(GDC_DATA_ENDPOINT, json=payload, stream=True) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with open(archive, "wb") as f, tqdm(
            total=total, unit="B", unit_scale=True, desc="Downloading"
        ) as bar:
            for chunk in r.iter_content(chunk_size=1024 * 256):
                f.write(chunk)
                bar.update(len(chunk))
    return archive

def extract_archive(archive: Path, out_dir: Path) -> None:
    with tarfile.open(archive, "r:gz") as tar:
        tar.extractall(out_dir)
    archive.unlink()
    print(f"  Extracted to {out_dir}")

def download_tcga(cancer_types: list[str] = None,
                  n_files_per_cancer: int = 20,
                  raw_dir: str = "data/raw/tcga") -> None:
    if cancer_types is None:
        cancer_types = list(CANCER_TYPES.keys())
    raw_dir = Path(raw_dir)

    for ct in cancer_types:
        print(f"\n── {ct} ({CANCER_TYPES[ct]}) ──")
        out = raw_dir / ct
        if out.exists() and any(out.iterdir()):
            print(f"  Already downloaded, skipping.")
            continue
        hits     = fetch_file_manifest(ct, n_files_per_cancer)
        file_ids = [h["file_id"] for h in hits]
        manifest = {h["file_id"]: h for h in hits}
        archive  = download_files(file_ids, out)
        extract_archive(archive, out)
        manifest_path = out / "manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2)
        print(f"  Manifest saved → {manifest_path}")

if __name__ == "__main__":
    download_tcga(n_files_per_cancer=20)