# UPDATED: Education level vs Tau (tTau & pTau) — creates a grouped bar chart and stops there

import csv
import warnings
import matplotlib.pyplot as plt
from collections import defaultdict
import statistics
import numpy as np

warnings.filterwarnings("ignore")

# ---------------------------
# 1) Patient class
# ---------------------------
class Patient: 
    all_patients = []
    education_lvl = {}

    def __init__(self, DonorID, ABeta40: float , ABeta42: float, tTau: float, pTau: float):
        self.DonorID = DonorID
        self.ABeta40 = ABeta40
        self.ABeta42 = ABeta42
        self.tTau = tTau
        self.pTau = pTau
        self.sex = None
        self.death_age = None
        self.ed_lvl = None
        self.cog_stat = None
        self.age_symp_on = None
        self.age_diag = None 
        self.head_inj = None
        self.thal_score = None
        Patient.all_patients.append(self)

    def __repr__(self):
        return (f"{self.DonorID} | sex: {self.sex} | ABeta40 {self.ABeta40} | "
                f"tTau {self.tTau} | pTau {self.pTau} | Death Age {self.death_age} | "
                f"Thal Score {self.thal_score}")

    def get_id(self):
        return self.DonorID

    # --- Metadata join ---
    @classmethod
    def combine_data(cls, filename: str):
        """Augment Patient objects with metadata from the metadata CSV."""
        with open(filename, encoding="utf8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            meta_by_id = {row["Donor ID"]: row for row in rows}

        for p in Patient.all_patients:
            row = meta_by_id.get(p.DonorID)
            if not row:
                warnings.warn(f"ID {p.DonorID} not found in metadata.")
                continue

            if row.get("Sex", ""):
                p.sex = row["Sex"]

            if row.get("Age at Death", ""):
                try:
                    p.death_age = int(row["Age at Death"])
                except ValueError:
                    pass

            if row.get("Highest level of education", ""):
                p.ed_lvl = row["Highest level of education"]

            if row.get("Cognitive Status", ""):
                p.cog_stat = row["Cognitive Status"]

            if row.get("Age of onset cognitive symptoms", ""):
                try:
                    p.age_symp_on = int(row["Age of onset cognitive symptoms"])
                except ValueError:
                    pass

            if row.get("Age of Dementia diagnosis", ""):
                try:
                    p.age_diag = int(row["Age of Dementia diagnosis"])
                except ValueError:
                    pass

            if row.get("Known head injury", ""):
                p.head_inj = row["Known head injury"]

            if row.get("Thal", ""):
                try:
                    p.thal_score = int(row["Thal"].split()[1])
                except Exception:
                    pass

    @classmethod
    def instantiate_from_csv(cls, luminex_file: str, meta_file: str):
        """Create Patient objects from the Luminex CSV, then enrich them with metadata."""
        with open(luminex_file, encoding="utf8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                Patient(
                    DonorID = row['Donor ID'],
                    ABeta40 = float(row['ABeta40 pg/ug']),
                    ABeta42 = float(row['ABeta42 pg/ug']),
                    tTau = float(row['tTAU pg/ug']),
                    pTau = float(row['pTAU pg/ug'])
                )
        # sort for consistency then enrich
        Patient.all_patients.sort(key=Patient.get_id)
        Patient.combine_data(meta_file)

# ---------------------------
# 2) Load your CSVs
# ---------------------------
luminex_file = "UpdatedLuminex.csv"   # exact filename you provided
meta_file    = "UpdatedMetaData.csv"  # exact filename you provided

Patient.all_patients.clear()
Patient.education_lvl.clear()
Patient.instantiate_from_csv(luminex_file, meta_file)

# ---------------------------
# 3) Compute mean tTau & pTau by education level
#    Set alzheimer_only = True to focus on dementia-labeled patients.
# ---------------------------
alzheimer_only = True

if alzheimer_only:
    patients_use = [p for p in Patient.all_patients
                    if (p.cog_stat is None) or (p.cog_stat.strip().lower() == "dementia")]
else:
    patients_use = Patient.all_patients

tTau_by_ed = defaultdict(list)
pTau_by_ed = defaultdict(list)

for p in patients_use:
    ed = p.ed_lvl if (p.ed_lvl is not None and p.ed_lvl != "") else "Unknown"
    if p.tTau is not None:
        tTau_by_ed[ed].append(p.tTau)
    if p.pTau is not None:
        pTau_by_ed[ed].append(p.pTau)

ed_levels = sorted(set(list(tTau_by_ed.keys()) + list(pTau_by_ed.keys())))

mean_tTau = [statistics.mean(tTau_by_ed[ed]) if tTau_by_ed[ed] else 0.0 for ed in ed_levels]
mean_pTau = [statistics.mean(pTau_by_ed[ed]) if pTau_by_ed[ed] else 0.0 for ed in ed_levels]

# ---------------------------
# 4) Bar graph (STOP HERE)
# ---------------------------
x = np.arange(len(ed_levels))
width = 0.35

plt.figure()
plt.bar(x - width/2, mean_tTau, width=width, label="Total Tau (tTau)")
plt.bar(x + width/2, mean_pTau, width=width, label="Phosphorylated Tau (pTau)")

plt.xticks(x, ed_levels, rotation=45, ha="right")
plt.ylabel("Average Tau Level (pg/μg)")
plt.title("Average Tau Levels by Education Level (Alzheimer’s Patients)")
plt.legend()
plt.tight_layout()
plt.show()
