# conda install -c conda-forge termcolor
# optional: pip/conda install scipy

import csv
import os
import warnings
import statistics
from typing import List, Dict, Any, Optional

import matplotlib.pyplot as plt
import numpy as np

# Optional color output
try:
    from termcolor import colored
except Exception:
    def colored(x, *_args, **_kwargs):
        return x

# Optional stats
try:
    from scipy import stats as scipy_stats  # type: ignore
except Exception:
    scipy_stats = None


class Patient:
    """
    Represents one subject with Luminex and Meta data merged by Donor ID.
    """

    all_patients: List["Patient"] = []
    death_age: List[int] = []          # kept from your original, not used directly
    education_lvl: Dict[str, List["Patient"]] = {}

    def __init__(self, DonorID: str, ABeta40: float, ABeta42: float, tTau: float, pTau: float):
        self.DonorID: str = DonorID
        self.ABeta40: float = ABeta40
        self.ABeta42: float = ABeta42
        self.tTau: float = tTau
        self.pTau: float = pTau

        # Will be filled from MetaData
        self.sex: Optional[str] = None
        self.death_age: Optional[int] = None
        self.ed_lvl: Optional[str] = None
        self.cog_stat: Optional[str] = None
        self.age_symp_on: Optional[int] = None
        self.age_diag: Optional[int] = None
        self.head_inj: Optional[str] = None
        self.thal_score: Optional[int] = None

        Patient.all_patients.append(self)

    def __repr__(self):
        return (f"{self.DonorID} | sex: {self.sex} | "
                f"ABeta40 {self.ABeta40} | ABeta42 {self.ABeta42} | "
                f"tTau {self.tTau} | pTau {self.pTau} | "
                f"Death Age {self.death_age} | Thal Score {self.thal_score}")

    # ----- getters used for sorting -----
    def get_id(self) -> str:
        return self.DonorID

    def get_ABeta42(self) -> float:
        return self.ABeta42

    def get_thal(self) -> int:
        # Fallback: treat None as very large so Nones sort to the end
        return self.thal_score if self.thal_score is not None else 10**9

    def get_death_age(self) -> int:
        return self.death_age if self.death_age is not None else -1

    # ---------- utilities ----------
    @staticmethod
    def _script_dir() -> str:
        """Folder where this script (patient.py) lives."""
        return os.path.dirname(__file__)

    @staticmethod
    def _full_path(filename: str) -> str:
        """Build absolute path to a file that sits next to patient.py."""
        return os.path.join(Patient._script_dir(), filename)

    @staticmethod
    def _require_file(fp: str) -> None:
        if not os.path.isfile(fp):
            raise FileNotFoundError(
                f"Could not find file:\n  {fp}\n"
                f"Make sure it sits in the same folder as patient.py."
            )

    @staticmethod
    def _safe_float(s: str) -> Optional[float]:
        try:
            s = s.strip()
            if s == "":
                return None
            return float(s)
        except Exception:
            return None

    @staticmethod
    def _safe_int(s: str) -> Optional[int]:
        try:
            s = s.strip()
            if s == "":
                return None
            return int(s)
        except Exception:
            return None

    # ---------- core loaders/mergers ----------
    @classmethod
    def instantiate_from_csv(cls, luminex_csv: str, meta_csv: str):
        """
        Load the Luminex file first (creates Patient objects), then join MetaData by Donor ID.
        File paths are resolved relative to this script's directory.
        """
        # Resolve absolute paths next to patient.py
        file1 = cls._full_path(luminex_csv)
        file2 = cls._full_path(meta_csv)
        cls._require_file(file1)
        cls._require_file(file2)

        # Clear any previous data if this is re-run
        cls.all_patients.clear()
        cls.education_lvl.clear()

        # ---- Read Luminex ----
        with open(file1, encoding="utf8", newline="") as f:
            reader = csv.DictReader(f)
            required_cols = {'Donor ID', 'ABeta40 pg/ug', 'ABeta42 pg/ug', 'tTAU pg/ug', 'pTAU pg/ug'}
            if not required_cols.issubset(set(reader.fieldnames or [])):
                raise KeyError(
                    f"Luminex CSV missing required columns. Found columns:\n{reader.fieldnames}\n"
                    f"Required: {sorted(required_cols)}"
                )

            for row in reader:
                did = (row.get('Donor ID') or "").strip()
                if not did:
                    warnings.warn("Skipping row with empty Donor ID in Luminex.")
                    continue

                a40 = cls._safe_float(row.get('ABeta40 pg/ug', ''))
                a42 = cls._safe_float(row.get('ABeta42 pg/ug', ''))
                ttau = cls._safe_float(row.get('tTAU pg/ug', ''))
                ptau = cls._safe_float(row.get('pTAU pg/ug', ''))

                # If any key value is missing, skip with warning
                if None in (a40, a42, ttau, ptau):
                    warnings.warn(f"Skipping Donor {did}: missing numeric Luminex value(s).")
                    continue

                Patient(
                    DonorID=did,
                    ABeta40=float(a40),
                    ABeta42=float(a42),
                    tTau=float(ttau),
                    pTau=float(ptau)
                )

        # Sort by ID before merging (not strictly required)
        cls.all_patients.sort(key=Patient.get_id)

        # Merge MetaData
        cls.combine_data(file2)

    @classmethod
    def combine_data(cls, meta_csv_fullpath: str):
        """
        Merge MetaData fields into existing Patient objects by Donor ID.
        Works even if the row orders are different between CSVs.
        """
        cls._require_file(meta_csv_fullpath)

        # Build quick lookup from DonorID -> Patient
        index: Dict[str, Patient] = {p.DonorID: p for p in cls.all_patients}

        with open(meta_csv_fullpath, encoding="utf8", newline="") as f:
            reader = csv.DictReader(f)
            # do not hard-fail on column names; attempt to read leniently
            for row in reader:
                did = (row.get("Donor ID") or "").strip()
                if not did or did not in index:
                    if did:
                        warnings.warn(f"Meta row for Donor ID '{did}' has no match in Luminex; skipping.")
                    continue

                p = index[did]

                sex = (row.get("Sex") or "").strip() or None
                p.sex = sex

                p.death_age = cls._safe_int(row.get("Age at Death", ""))
                p.ed_lvl = (row.get("Highest level of education") or "").strip() or None
                p.cog_stat = (row.get("Cognitive Status") or "").strip() or None
                p.age_symp_on = cls._safe_int(row.get("Age of onset cognitive symptoms", ""))
                p.age_diag = cls._safe_int(row.get("Age of Dementia diagnosis", ""))
                p.head_inj = (row.get("Known head injury") or "").strip() or None

                # Thal can appear like "Thal 2" or just "2"; handle both
                thal_raw = (row.get("Thal") or "").strip()
                if thal_raw:
                    # try "Thal 2"
                    parts = thal_raw.split()
                    last = parts[-1] if parts else thal_raw
                    p.thal_score = cls._safe_int(last)

    # ---------- grouping/sorting ----------
    @classmethod
    def sort_ed(cls):
        """
        Build dict {education level: [patients with that level]}.
        """
        cls.education_lvl.clear()
        for patient in cls.all_patients:
            key = patient.ed_lvl or "Unknown"
            if key not in cls.education_lvl:
                cls.education_lvl[key] = []
            cls.education_lvl[key].append(patient)

    @classmethod
    def subsort_thal(cls):
        """
        Sort each education group's list by Thal score (None values go to end).
        """
        for key, values in cls.education_lvl.items():
            values.sort(key=Patient.get_thal)
            cls.education_lvl[key] = values

    # ---------- filtering ----------
    @classmethod
    def filter(cls,
               patients: List["Patient"],
               ABeta40: Any = "any",
               ABeta42: Any = "any",
               tTau: Any = "any",
               pTau: Any = "any",
               sex: Any = "any",
               death_age: Any = "any",
               ed_lvl: Any = "any",
               cog_stat: Any = "any",
               age_symp_on: Any = "any",
               age_diag: Any = "any",
               head_inj: Any = "any",
               thal_score: Any = "any") -> List["Patient"]:
        """
        Generic equality filter across provided attributes.
        Usage stays compatible with your original calls.
        """
        attr_list = (
            ABeta40, ABeta42, tTau, pTau, sex, death_age, ed_lvl,
            cog_stat, age_symp_on, age_diag, head_inj, thal_score
        )
        attr_name = (
            "ABeta40", "ABeta42", "tTau", "pTau", "sex", "death_age",
            "ed_lvl", "cog_stat", "age_symp_on", "age_diag", "head_inj",
            "thal_score"
        )

        result = patients[:]
        for value, name in zip(attr_list, attr_name):
            if value == "any":
                continue
            # Keep patients where attribute equals the desired value
            keep = []
            for p in result:
                if getattr(p, name) == value:
                    keep.append(p)
            result = keep
        return result


# -------------------- convenience functions for plotting --------------------

def _mean_and_stdev(values: List[float]) -> (float, float):
    if len(values) == 0:
        return float("nan"), float("nan")
    if len(values) == 1:
        # stdev undefined for n=1; use 0 for errorbar aesthetics
        return float(values[0]), 0.0
    return statistics.mean(values), statistics.stdev(values)


def _ttest_ind_safe(a: List[float], b: List[float]) -> Optional[tuple]:
    """
    Two-sample t-test if scipy is available and both groups have n>=2; else None.
    """
    if scipy_stats is None:
        return None
    if len(a) < 2 or len(b) < 2:
        return None
    try:
        return scipy_stats.ttest_ind(a, b, equal_var=True)  # you can change to Welch if needed
    except Exception:
        return None


# -------------------- main script (only runs if executed directly) --------------------

if __name__ == "__main__":
    # Filenames expected to sit next to patient.py
    LUMINEX_FILE = "UpdatedLuminex.csv"
    META_FILE = "UpdatedMetaData.csv"

    # 1) Load and merge
    Patient.instantiate_from_csv(LUMINEX_FILE, META_FILE)

    # 2) Print merged patients
    print(colored("\n=== Patients (merged) ===", "green"))
    for patient in Patient.all_patients:
        print(patient)

    # 3) Sort by ABeta42 (ascending) and print
    Patient.all_patients.sort(key=Patient.get_ABeta42, reverse=False)
    print(colored("\n=== Patients sorted by ABeta42 (asc) ===", "green"))
    for patient in Patient.all_patients:
        print(patient)

    # 4) Build education groups and sub-sort by Thal
    Patient.sort_ed()
    Patient.subsort_thal()

    print(colored("\n=== Groups by Education (sub-sorted by Thal) ===", "green"))
    for key in Patient.education_lvl:
        print(colored(key, "red"))
        for patient in Patient.education_lvl.get(key, []):
            print(patient)
        print()

    # 5) Use filter to get counts by sex + cog status
    fem_healthy = Patient.filter(Patient.all_patients, sex="Female", cog_stat="No dementia")
    male_healthy = Patient.filter(Patient.all_patients, sex="Male", cog_stat="No dementia")
    fem_diseased = Patient.filter(Patient.all_patients, sex="Female", cog_stat="Dementia")
    male_diseased = Patient.filter(Patient.all_patients, sex="Male", cog_stat="Dementia")

    print(colored("=== Counts ===", "green"))
    print(f"Female Healthy Patients = {len(fem_healthy)} | Male Healthy Patients = {len(male_healthy)}")
    print(f"Female Diseased Patients = {len(fem_diseased)} | Male Diseased Patients = {len(male_diseased)}")

    # 6) Bar graph: ABeta42 mean ± stdev by sex
    fem_vals = [p.ABeta42 for p in Patient.filter(Patient.all_patients, sex="Female")]
    male_vals = [p.ABeta42 for p in Patient.filter(Patient.all_patients, sex="Male")]

    x_fem_bar, fem_sd = _mean_and_stdev(fem_vals)
    x_male_bar, male_sd = _mean_and_stdev(male_vals)

    print(colored("\n=== ABeta42 summary ===", "green"))
    print(f"x_fem_bar = {x_fem_bar}, ABeta_fem_stdev = {fem_sd}")
    print(f"x_male_bar = {x_male_bar}, ABeta_male_stdev = {male_sd}")

    # Optional t-test
    ttest = _ttest_ind_safe(fem_vals, male_vals)
    if ttest is not None:
        t_stat, p_val = ttest
        print(f"t-statistic: {t_stat}")
        print(f"p-value: {p_val}")
    else:
        print("(t-test skipped: need scipy and at least 2 values in each group)")

    # Plot
    sex_cols = ['Female', 'Male']
    means = [x_fem_bar, x_male_bar]
    yerr = [fem_sd, male_sd]

    plt.bar(sex_cols, means, yerr=yerr, capsize=10)
    plt.title("ABeta42 Levels")
    plt.xlabel("Sex")
    plt.ylabel("ABeta42 (pg/ug)")
    plt.tight_layout()
    plt.show()
