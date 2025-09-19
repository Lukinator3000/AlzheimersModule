#1) BUILD OUR CLASS BY COMBINING TWO .csv FILES OF DATA
 
import csv
import warnings
import re
import matplotlib.pyplot as plt

class Patient: 

    all_patients = []

    death_age = []

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
        self.oadc = None  # Overall AD neuropathological Change (Not AD / Low / Intermediate / High)
        Patient.all_patients.append(self)

    def __repr__(self):
        return f"{self.DonorID} | sex: {self.sex} | ABeta40 {self.ABeta40} | tTau {self.tTau} | pTau {self.pTau} | Death Age {self.death_age} | Thal Score {self.thal_score}"

    def get_id(self):
        return self.DonorID

    def get_ABeta42(self):
        return self.ABeta42
    
    def get_thal(self):
        return self.thal_score
    
    def get_death_age(self):
        return self.death_age


    @classmethod
    def combine_data(cls, filename: str):
        # Read metadata
        with open(filename, encoding="utf8") as f:
            reader = csv.DictReader(f)
            meta_rows = list(reader)
            fieldnames = reader.fieldnames or []

        # Build a fast lookup by Donor ID
        meta_by_id = {}
        for r in meta_rows:
            did = (r.get("Donor ID") or "").strip()
            if did:
                meta_by_id[did] = r

        # Identify the OADC/ADNC column, robust to naming variants
        oadc_col = None
        candidates = []
        for fn in fieldnames:
            low = fn.lower()
            if (
                "oadc" in low
                or "adnc" in low
                or (
                    "overall" in low
                    and ("ad" in low or "alzheimer" in low)
                    and ("neuro" in low or "neuropath" in low)
                    and ("change" in low or "nc" in low)
                )
            ):
                candidates.append(fn)
        if candidates:
            # Prefer the shortest, most specific header
            oadc_col = sorted(candidates, key=len)[0]

        def _norm_oadc_val(v):
            if v is None:
                return None
            s = str(v).strip()
            if not s or s.lower() in {"na", "nan", "none"}:
                return None
            sl = s.lower()

            # Try numeric codes first (0–3 sometimes stored as "0", "1.0", etc.)
            try:
                n = int(float(s))
                mapping = {0: "Not AD", 1: "Low", 2: "Intermediate", 3: "High"}
                if n in mapping:
                    return mapping[n]
            except Exception:
                pass

            # Keyword mapping (also catches forms like "ADNC-High", "OADC: low", etc.)
            if "not" in sl and "ad" in sl:
                return "Not AD"
            if "low" in sl:
                return "Low"
            if "inter" in sl:
                return "Intermediate"
            if "high" in sl:
                return "High"

            # Fallback: title-case whatever is there
            return s.title()

        # Update instantiated patients
        for p in Patient.all_patients:
            r = meta_by_id.get(p.DonorID)
            if not r:
                warnings.warn(f"No metadata found for Donor ID {p.DonorID}")
                continue

            if r.get("Sex"):
                p.sex = r["Sex"]

            if r.get("Age at Death"):
                try:
                    p.death_age = int(r["Age at Death"])
                except ValueError:
                    p.death_age = None

            if r.get("Highest level of education"):
                p.ed_lvl = r["Highest level of education"]

            if r.get("Cognitive Status"):
                p.cog_stat = r["Cognitive Status"]

            if r.get("Age of onset cognitive symptoms"):
                try:
                    p.age_symp_on = int(r["Age of onset cognitive symptoms"])
                except ValueError:
                    p.age_symp_on = None

            if r.get("Age of Dementia diagnosis"):
                try:
                    p.age_diag = int(r["Age of Dementia diagnosis"])
                except ValueError:
                    p.age_diag = None

            if r.get("Known head injury"):
                p.head_inj = r["Known head injury"]

            if r.get("Thal"):
                m = re.search(r"(\d+)", r["Thal"])
                p.thal_score = int(m.group(1)) if m else None

            # Robust OADC extraction/normalization
            raw_oadc = r.get(oadc_col) if oadc_col else None
            p.oadc = _norm_oadc_val(raw_oadc)

    @classmethod
    def instantiate_from_csv(cls, filename: str, other_file: str):
        # open csv and create list of all rows
        with open(filename, encoding="utf8") as f:
            reader = csv.DictReader(f)
            rows_of_patients = list(reader)

            # for line in csv create object
            for row in rows_of_patients:
                Patient(
                    DonorID = (row['Donor ID'] or '').strip(),
                    ABeta40 = float(row['ABeta40 pg/ug']),
                    ABeta42 = float(row['ABeta42 pg/ug']),
                    tTau    = float(row['tTAU pg/ug']),
                    pTau    = float(row['pTAU pg/ug'])
                )

            # sort by Donor ID so printing is stable; metadata join does not rely on row order
            Patient.all_patients.sort(key = Patient.get_id)
            Patient.combine_data(other_file)
