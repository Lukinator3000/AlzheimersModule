#2) PRINT OUR LIST OF PATIENTS

from patient import Patient

Patient.instantiate_from_csv('UpdatedLuminex.csv', 'UpdatedMetaData.csv')

for patient in Patient.all_patients:
    print(patient)

#3) SORT OUR LIST OF PATIENTS

Patient.all_patients.sort(key=Patient.get_ABeta42, reverse=False)

for patient in Patient.all_patients:
    print(patient)

#5) SORT THE SUB-LISTS WITHIN THE LIST AND PRINT IT

from termcolor import colored

Patient.sort_ed()
Patient.subsort_thal()
#change to fit Tau score function above!!!!!

for key in Patient.education_lvl:
    print(colored(key,"red"))
    for patient in Patient.education_lvl.get(key):
        print(patient)
    print()

#7) USE OUR FILTER TO PULL OUT AND COUNT SUB-SETS OF PATIENTS THAT ARE FILTERED

fem_healty_patients = range(len(Patient.filter(Patient.all_patients, sex = "Female", cog_stat = "No dementia")))
male_healthy_patients = range(len(Patient.filter(Patient.all_patients, sex = "Male", cog_stat = "No dementia")))
fem_diseased_patients = range(len(Patient.filter(Patient.all_patients, sex = "Female", cog_stat = "Dementia")))
male_diseased_patients = range(len(Patient.filter(Patient.all_patients, sex = "Male", cog_stat = "Dementia")))

print(f'Female Healthy Patients = {len(fem_healty_patients)} | Male Healthy Patients = {len(male_healthy_patients)}')
print(f'Female Diseased Patients = {len(fem_diseased_patients)} | Male Diseased Patients = {len(male_diseased_patients)}')

#8) PLOT A BAR GRAPH OF THE SORTED DATA

from patient import Patient
from termcolor import colored
from matplotlib import pyplot as plt
from scipy import stats
import numpy as np
import statistics 


# Patient.instantiate_from_csv('UpdatedLuminex.csv', 'UpdatedMetaData.csv')

ABeta42_fem_vals = []
ABeta42_male_vals = []

for patient in Patient.filter(Patient.all_patients, sex = "Female"):
     ABeta42_fem_vals.append(patient.ABeta42)
for patient in Patient.filter(Patient.all_patients, sex = "Male"):
     ABeta42_male_vals.append(patient.ABeta42)

x_fem_bar = (statistics.mean(ABeta42_fem_vals))
x_male_bar = (statistics.mean(ABeta42_male_vals))

ABeta_fem_stdev = (statistics.stdev(ABeta42_fem_vals))
ABeta_male_stdev = (statistics.stdev(ABeta42_male_vals))

print(f'x_fem_bar = {x_fem_bar}, ABeta_fem_stdev {ABeta_fem_stdev}')
print(f'x_male_bar = {x_male_bar}, ABeta_male_stdev {ABeta_male_stdev}')

x_fem_vals = range(len(Patient.filter(Patient.all_patients, sex = "Female")))
x_male_vals = range(len(Patient.filter(Patient.all_patients, sex = "Male")))

sex_cols = ['Female', 'Male']
mean_sex_ABeta42 = [x_fem_bar, x_male_bar]
stdev_sex_ABeta42 = [ABeta_fem_stdev, ABeta_male_stdev]
colors = ["pink", "blue"]
yerr = [np.zeros(len(mean_sex_ABeta42)), stdev_sex_ABeta42]

# Equal variance assumed
t_stat, p_val = stats.ttest_ind(ABeta42_fem_vals, ABeta42_male_vals)
print("t-statistic:", t_stat)
print("p-value:", p_val)

plt.bar(sex_cols, mean_sex_ABeta42, yerr=yerr, capsize=10, color=["red", "blue"])
plt.title("ABeta42 Levels")
plt.xlabel("Sex")
plt.ylabel("Abeta42")
plt.show()