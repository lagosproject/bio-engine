import time
import logging
import json
from services.job_manager import JobManager
from tasks.worker import process_job_background
from data.models import JobStatus

# Enhanced logging to see what's happening internally
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_test_job():
    jm = JobManager()
    
    # Configuration matches what the API would send
    job_name = "Test NG_008866.1 Real Read"
    reference = {"type": "ncbi", "value": "NG_008866.1"}
    patients = [{
        "id": "patient1",
        "name": "Patient 1",
        "reads": [{
            "file": "/home/vant/Documentos/Trabajo/Secuencias/Sanger/18-138498_101_102___R-_B3_2_20260121_163604.ab1"
        }]
    }]
    hgvs_config = {
        "transcript": "NG_008866.1",
        "auto_hgvs": True,
        "auto_vep": True
    }
    # Standard tracy config
    config = {
        "trimLeft": 30,
        "trimRight": 30
    }

    # 1. Create Job
    print("\n--- Creating Job ---")
    job = jm.create_job(
        name=job_name,
        reference=reference,
        patients=patients,
        hgvs_config=hgvs_config,
        config=config
    )
    job_id = job.id
    print(f"Job created with ID: {job_id}")

    # 2. Run Job
    # We call it directly here instead of background_tasks.add_task
    # to wait for it synchronously in this test script
    print("\n--- Starting Analysis Pipeline ---")
    try:
        process_job_background(job_id)
    except Exception as e:
        print(f"Fatal error in background process: {e}")

    # 3. Final Check
    job = jm.get_job(job_id)
    print("\n--- FINAL JOB STATUS ---")
    print(f"Status: {job.status}")
    print(f"Progress: {job.progress}%")
    print(f"Status Message: {job.status_message}")

    # 4. Show Analysis results
    print("\n--- ANALYSIS RESULTS ---")
    print(f"HGVS Alternatives Count: {len(job.hgvs_alternatives)}")
    print(f"VEP Annotations Count: {len(job.vep_annotations)}")
    
    if job.vep_annotations:
        print("\nVEP Annotations Summary:")
        for hgvs, ann in job.vep_annotations.items():
            print(f"  - {hgvs}: {ann.gene_symbol} ({ann.consequence})")
    else:
        print("\n[!] WARNING: vep_annotations is EMPTY!")

    # Check one result read to see the 'hgvs' column
    if job.results:
        res = job.results[0]
        variants = res.get("alignment", {}).get("variants", {})
        header = variants.get("columns", [])
        rows = variants.get("rows", [])
        if "hgvs" in header:
            h_idx = header.index("hgvs")
            print("\nVariants found (primary HGVS):")
            for row in rows:
                print(f"  - {row[h_idx]}")

if __name__ == "__main__":
    run_test_job()
