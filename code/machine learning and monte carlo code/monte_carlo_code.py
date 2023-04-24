import matlab.engine
import csv
import random
import numpy as np
from pymatgen.io.vasp.inputs import Poscar
import pandas as pd
from sklearn.svm import SVR
from sklearn import linear_model
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import cross_validate
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
import os

"""Temperature"""
T = 473

"""Load data"""
df = pd.read_csv("MgAl2O4-dataset.csv")
X = df.iloc[:, 1:24]
y = df.iloc[:, 24]

""" Machine learning part """
"""  Support vector machine  """
Model = SVR(C=6.2, epsilon=0.001, kernel='poly')
Model.fit(X, y)


""" Matlab """
eng = matlab.engine.start_matlab()  # load matlab code.

"""Monte Carlo Simulation"""


def swap_atom(poscar_infile):
    """
    swap the positins of two ions in the structure.
    :param poscar_infile: the POSCAR file of structure before swapping positions.
    :return: the bond numbers counted by matlab code (eng.count_bond_in_single_spinel). [num1, num2, ..., num23]
    """
    with open(poscar_infile, "r") as f:
        content = csv.reader(f)
        data = [i for i in content]
        f.close()
    f.close()
    pos1 = random.randint(1, 24)
    pos2 = random.randint(1, 24)
    coord = data[8:]
    coord[pos1 - 1], coord[pos2 - 1] = coord[pos2 - 1], coord[pos1 - 1]
    with open("POSCAR-new", "w") as f:
        for i in range(8):
            f.write(data[i][0])
            f.write("\n")
        for i in range(len(coord)):
            f.write(coord[i][0])
            f.write('\n')
        f.close()
    f.close()
    bond_num = eng.count_bond_in_single_spinel("POSCAR-new") # for single spinel, use matlab code of count_bond_in_single_spinel("POSCAR-new")
    # bond_num = eng.count_bond_in_double_spinel("POSCAR-new") # for double spinel, use matlab code of count_bond_in_double_spinel("POSCAR-new")
    return bond_num


def cal_energy(data, model):
    """
    predict the energy of structure using machine learning model.
    :param data: structural features, [num1, num2, ..., num23]
    :param model: machine learning model
    :return: a number.
    """
    return model.predict(data)[0]


def calX(file):
    """
    calculate the degree of inversion of a given structure.
    :param file: POSCAR file of s structure.
    :return: for example, in MgFe2O4, return fraction of tetrahedral sites occupied by Fe ion.
    """
    with open(file, "r") as f:
        content = csv.reader(f)
        data = [i for i in content]
        f.close()
    f.close()
    site_info_Mg = data[8: 16]
    site_info_Fe = data[16: 32]
    sites_Mg = []
    for i in range(len(site_info_Mg)):
        sites_Mg.append(site_info_Mg[i][0].split())
    sites_float_Mg = []
    for i in range(len(sites_Mg)):
        tmp = []
        for j in range(len(sites_Mg[i])):
            tmp.append(float(sites_Mg[i][j]))
        sites_float_Mg.append(tmp)
    count_Mg = 0
    for i in range(len(sites_float_Mg)):
        if sites_float_Mg[i][0] in [0, 0.25, 0.5, 0.75, 1]:
            count_Mg += 1

    sites_Fe = []
    for i in range(len(site_info_Fe)):
        sites_Fe.append(site_info_Fe[i][0].split())
    sites_float_Fe = []
    for i in range(len(sites_Fe)):
        tmp = []
        for j in range(len(sites_Fe[i])):
            tmp.append(float(sites_Fe[i][j]))
        sites_float_Fe.append(tmp)
    count_Fe = 0
    for i in range(len(sites_float_Fe)):
        if sites_float_Fe[i][0] in [0, 0.25, 0.5, 0.75, 1]:
            count_Fe += 1
    return count_Fe / 8


iteration = 100000 # in Monte Carlo simulation, we run 100,000 steps.
initial_file = "POSCAR-initial"
print(eng.count_bond_in_single_spinel(initial_file))
energy_start = cal_energy(eng.count_bond_in_single_spinel(initial_file), model)
print(energy_start)
energy_set = []  # save the energies of 100,000 steps in MC simulation.
x_set = []   # save the degree of inversion of 100,000 steps in MC simulation.
bond_set = []   # save the structure features of 100,000 steps in MC simulation.
for i in range(iteration):
    bond_set.append(eng.count_bond_in_single_spinel(initial_file)[0])
    energy_set.append(energy_start)
    x_set.append(calX(initial_file))
    bond_num = swap_atom(initial_file)
    energy_new = cal_energy(bond_num, model)
    dE = energy_new - energy_start
    kbT = 0.0000860 * (T)  # 0.00008617333262145
    p = np.exp(-dE / kbT)
    if min(1, p) >= 1:
        f = open("POSCAR-new", "rb")
        f2 = open("POSCAR-initial", "wb")
        for line in f.readlines():
            f2.write(line)
        f.close()
        f2.close()
        energy_start = energy_new
    if min(1, p) < 1:
        random_p = random.random()
        if random_p < p:
            f = open("POSCAR-new", "rb")
            f2 = open("POSCAR-initial", "wb")
            for line in f.readlines():
                f2.write(line)
            f.close()
            f2.close()
            energy_start = energy_new


df_x = pd.DataFrame(x_set)
df_x.to_csv("x-{}K.csv".format(T))  # save the degree of inversion into a csv file.
df_energy = pd.DataFrame(energy_set)
df_energy.to_csv("energy-{}K.csv".format(T))   # save the energies into a csv file.