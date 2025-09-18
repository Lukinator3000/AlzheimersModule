#1) BUILD OUR CLASS BY COMBINING TWO .csv FILES OF DATA

import csv
import os
import warnings
import statistics
import matplotlib.pyplot as plt


class Patient:

    all_patients = []
    death_age = []
    education_lvl = {}

    def __init__(self, DonorID, ABeta40: float, ABeta42: float, tTau: float, pTau: float):
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

    def get_ABeta42(self):
        return self.ABeta42

    def get_thal(self):
        return self.thal_score if self.thal_score is not None else 10**9

    def get_death_age(self):
        return self.death_age if self.death_age is not None else -1

    @staticmethod
    def _here_path(filename: str) -> str:
        base_dir = os.path.dirname(__file__)
        return os.path.join(base_dir, filename)

    @classmethod
    def combine_data(cls, filename: str):
        meta_path = cls._here_path(filename)
        with open(meta_path, encoding="utf8", newline="") as f:
            reader = csv.DictReader(f)
            rows_of_meta = list(reader)

        index = {p.DonorID: p for p in Patient.all_patients}

        for row in rows_of_meta:
            did = (row.get("Donor ID") or "").strip()
            if not did:
                continue
            if did not in index:
                warnings.warn(f"Meta row Donor ID '{did}' not in Luminex; skipping.")
                continue

            p = index[did]

            if row.get("Sex", "") != "":
                p.sex = row["Sex"]

            if row.get("Age at Death", "") != "":
                try:
                    p.death_age = int(row["Age at Death"])
                except:
                    pass

            if row.get("Highest level of education", "") != "":
                p.ed_lvl = row["Highest level of education"]

            if row.get("Cognitive Status", "") != "":
                p.cog_stat = row["Cognitive Status"]

            if row.get("Age of onset cognitive symptoms", "") != "":
                try:
                    p.age_symp_on = int(row["Age of onset cognitive symptoms"])
                except:
                    pass

            if row.get("Age of Dementia diagnosis", "") != "":
                try:
                    p.age_diag = int(row["Age of Dementia diagnosis"])
                except:
                    pass

            if row.get("Known head injury", "") != "":
                p.head_inj = row["Known head injury"]

            if row.get("Thal", "") != "":
                raw = row["Thal"].strip()
                try:
                    last = raw.split()[-1]
                    p.thal_score = int(last)
                except:
                    pass

    @classmethod
    def instantiate_from_csv(cls, filename: str, other_file: str):
        luminex_path = cls._here_path(filename)
        with open(luminex_path, encoding="utf8", newline="") as f:
            reader = csv.DictReader(f)
            rows_of_patients = list(reader)
            for row in rows_of_patients:
                try:
                    Patient(
                        DonorID=row['Donor ID'],
                        ABeta40=float(row['ABeta40 pg/ug']),
                        ABeta42=float(row['ABeta42 pg/ug']),
                        tTau=float(row['tTAU pg/ug']),
                        pTau=float(row['pTAU pg/ug'])
                    )
                except:
                    warnings.warn(f"Skipping row with missing numeric values for Donor ID {row.get('Donor ID')}")

        Patient.all_patients.sort(key=Patient.get_id)
        Patient.combine_data(other_file)

#4) BUILD DICTIONARIES TO DO SORTING AND SUB-SORTING
    @classmethod
    def sort_ed(cls):
        Patient.education_lvl.clear()
        for patient in Patient.all_patients:
            key = patient.ed_lvl
            if key is None:
                continue
            if key not in Patient.education_lvl:
                Patient.education_lvl[key] = []
            Patient.education_lvl[key].append(patient)

    @classmethod
    def subsort_thal(cls):
        for key in Patient.education_lvl:
            values = Patient.education_lvl.get(key, [])
            values.sort(key=Patient.get_thal)
            Patient.education_lvl[key] = values

#6) MAKE A FILTER TO PULL OUT PATIENTS WITH SPECIFIC ATTRIBUTES
    @classmethod
    def filter(cls, plist, ABeta40: float = "any", ABeta42: float = "any", tTau: float = "any", pTau: float = "any",
               sex: str = "any", death_age: int = "any", ed_lvl: str = "any", cog_stat: str = "any",
               age_symp_on: int = "any", age_diag: int = "any", head_inj: str = "any", thal_score: int = "any"):
        all_patients = list(plist)
        remove_list = []
        attr_list = (
            ABeta40,
            ABeta42,
            tTau,
            pTau,
            sex,
            death_age,
            ed_lvl,
            cog_stat,
            age_symp_on,
            age_diag,
            head_inj,
            thal_score
        )
        attr_name = (
            "ABeta40",
            "ABeta42",
            "tTau",
            "pTau",
            "sex",
            "death_age",
            "ed_lvl",
            "cog_stat",
            "age_symp_on",
            "age_diag",
            "head_inj",
            "thal_score"
        )
        for i in range(len(attr_list)):
            if attr_list[i] != "any":
                for patient in all_patients:
                    if getattr(patient, attr_name[i]) != attr_list[i]:
                        remove_list.append(patient)
                all_patients = [patient for patient in all_patients if patient not in remove_list]
                remove_list.clear()

        return all_patients

#2) PRINT OUR LIST OF PATIENTS / RUN SCRIPT

if __name__ == "__main__":
    Patient.instantiate_from_csv('UpdatedLuminex.csv', 'UpdatedMetaData.csv')

    for patient in Patient.all_patients:
        print(patient)

#3) SORT OUR LIST OF PATIENTS
    Patient.all_patients.sort(key=Patient.get_ABeta42, reverse=False)

    for patient in Patient.all_patients:
        print(patient)

#5) SORT THE SUB-LISTS WITHIN THE LIST AND PRINT IT
    Patient.sort_ed()
    Patient.subsort_thal()

    for key in Patient.education_lvl:
        print(key)
        for patient in Patient.education_lvl.get(key, []):
            print(patient)
        print()

#7) USE OUR FILTER TO PULL OUT AND COUNT SUB-SETS OF PATIENTS THAT ARE FILTERED
    fem_healty_patients = range(len(Patient.filter(Patient.all_patients, sex="Female", cog_stat="No dementia")))
    male_healthy_patients = range(len(Patient.filter(Patient.all_patients, sex="Male", cog_stat="No dementia")))
    fem_diseased_patients = range(len(Patient.filter(Patient.all_patients, sex="Female", cog_stat="Dementia")))
    male_diseased_patients = range(len(Patient.filter(Patient.all_patients, sex="Male", cog_stat="Dementia")))

    print(f'Female Healthy Patients = {len(fem_healty_patients)} | Male Healthy Patients = {len(male_healthy_patients)}')
    print(f'Female Diseased Patients = {len(fem_diseased_patients)} | Male Diseased Patients = {len(male_diseased_patients)}')

#8) PLOT BAR GRAPHS

    # Build groups: {education: [values]}
    tTau_groups = {}
    pTau_groups = {}

    for patient in Patient.all_patients:
        lvl = patient.ed_lvl
        if not lvl:
            continue
        if patient.tTau is not None:
            tTau_groups.setdefault(lvl, []).append(float(patient.tTau))
        if patient.pTau is not None:
            pTau_groups.setdefault(lvl, []).append(float(patient.pTau))

    # Sort labels alphabetically
    tTau_labels = sorted(tTau_groups.keys())
    pTau_labels = sorted(pTau_groups.keys())

    # Compute means and SDs
    tTau_means = [statistics.mean(tTau_groups[k]) for k in tTau_labels]
    tTau_sds = [statistics.stdev(tTau_groups[k]) if len(tTau_groups[k]) > 1 else 0 for k in tTau_labels]

    pTau_means = [statistics.mean(pTau_groups[k]) for k in pTau_labels]
    pTau_sds = [statistics.stdev(pTau_groups[k]) if len(pTau_groups[k]) > 1 else 0 for k in pTau_labels]

    # Plot tTau by education
    plt.bar(tTau_labels, tTau_means, yerr=tTau_sds, capsize=10)
    plt.title("tTau by Education Level")
    plt.xlabel("Education Level")
    plt.ylabel("tTau (pg/mL)")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.show()

    # Plot pTau by education
    plt.bar(pTau_labels, pTau_means, yerr=pTau_sds, capsize=10)
    plt.title("pTau by Education Level")
    plt.xlabel("Education Level")
    plt.ylabel("pTau (pg/mL)")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.show()