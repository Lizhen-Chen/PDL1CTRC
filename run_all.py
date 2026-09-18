"""Reproduce the analyses from bundled public processed inputs."""
from pathlib import Path
import subprocess,sys,os
W=Path(__file__).resolve().parent;os.chdir(W)
steps=['analyze_singlecell.py','analyze_clinical_followup.py','analyze_spatial.py','analyze_additional_sensitivity.py','prepare_bulk_expression.py','analyze_bulk.py','analyze_mr.py']
for name in steps:subprocess.run([sys.executable,name],check=True)
subprocess.run(['Rscript','--vanilla','analyze_nhanes.R','.'],check=True)
subprocess.run(['Rscript','--vanilla','check_survival.R'],check=True)
subprocess.run([sys.executable,'prepare_tables.py'],check=True)
subprocess.run([sys.executable,'draw_figures.py'],check=True)
