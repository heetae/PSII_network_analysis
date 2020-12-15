#Computational analysis of FRET rates between Chlorophylls

#import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

#load coordinates of Chlorophylls
    #data set should be like below
    #Molecule   Label   Element X   Y   Z
    #example)
    #CLA	A 405	NB	244.604	280.744	222.325
    #CLA	A 405	ND	247.173	283.695	223.102
    #CLA	A 405	NA	247.45	281.078	221.843
    #CLA	A 405	NC	244.455	283.404	223.519
    #CLA	A 405	MG	245.864	282.315	222.423

DataFiles = ['database_all_a.csv', 'database_major_a.csv', 'database_minor_a.csv', 'database_natural.csv', 'database_minor_b.csv', 'database_major_b.csv', 'database_all_b.csv']

for database in DataFiles:
    temp = pd.read_csv(database)
    temp.replace("\t", " ")
    labeltemp = temp["Label"]
    label2 = labeltemp.drop_duplicates(keep='first')
    label3 = label2.reset_index()
    label = []
    label = label3["Label"]
    data = temp.set_index(["Label","Element"])
    result=np.array(["Source","Target","Weight"])
    n = len(label)
    
    #set n_index (Reference "Gradinaru et al., Biophysics Journal 1998")
    n_index = 1.55
    
    #Calculate FRET rates
    for i in range(0,n):
        Donor = label[i]
        for j in range(0,n):
            Acceptor = label[j]
    
            if Donor != Acceptor:
    
                Rx = data.loc[(Donor, "MG"), "X"] - data.loc[(Acceptor, "MG"), "X"]
                Ry = data.loc[(Donor, "MG"), "Y"] - data.loc[(Acceptor, "MG"), "Y"]
                Rz = data.loc[(Donor, "MG"), "Z"] - data.loc[(Acceptor, "MG"), "Z"]
                R = (Rx ** 2 + Ry ** 2 + Rz ** 2) ** (1 / 2)
                Rux = Rx / R
                Ruy = Ry / R
                Ruz = Rz / R
                Dx = data.loc[(Donor, "ND"), "X"] - data.loc[(Donor, "NB"), "X"]
                Dy = data.loc[(Donor, "ND"), "Y"] - data.loc[(Donor, "NB"), "Y"]
                Dz = data.loc[(Donor, "ND"), "Z"] - data.loc[(Donor, "NB"), "Z"]
                D = (Dx ** 2 + Dy ** 2 + Dz ** 2) ** (1 / 2)
                Dux = Dx / D
                Duy = Dy / D
                Duz = Dz / D

                Ax = data.loc[(Acceptor, "ND"), "X"] - data.loc[(Acceptor, "NB"), "X"]
                Ay = data.loc[(Acceptor, "ND"), "Y"] - data.loc[(Acceptor, "NB"), "Y"]
                Az = data.loc[(Acceptor, "ND"), "Z"] - data.loc[(Acceptor, "NB"), "Z"]
                A = (Ax ** 2 + Ay ** 2 + Az ** 2) ** (1 / 2)
                Aux = Ax / A
                Auy = Ay / A
                Auz = Az / A

                DA = Dux * Aux + Duy * Auy + Duz * Auz
                DR = Dux * Rux + Duy * Ruy + Duz * Ruz
                AR = Aux * Rux + Auy * Ruy + Auz * Ruz

                Kfs = (DA - 3 * DR * AR) ** 2

                # constants (nm^6 ps^-1) (Reference "Gradinaru et al., Biophysics Journal 1998")
                DChl = data.loc[(Donor, "MG"), "Mol"]
                AChl = data.loc[(Acceptor, "MG"), "Mol"]

                if DChl == "CLA":

                    if AChl == "CLA":

                        C_DA = 32.26
                        # Chl a -> Chl a
                    else:

                        C_DA = 1.11
                        # Chl a -> Chl b
                else:
                    if AChl == "CLA":

                        C_DA = 9.61
                        # Chl b -> Chl a
                    else:

                        C_DA = 14.45
                        # Chl b -> Chl b

                FRETrate = (C_DA * Kfs) / ((n_index ** 4) * ((R * 0.1) ** 6))
                Dprotein = Donor.split(" ")[0]
                Dmol = Donor.split(" ")[1]
                Aprotein = Acceptor.split(" ")[0]
                Amol = Acceptor.split(" ")[1]
                result = np.vstack((result, [Dprotein + "_" + Dmol, Aprotein + "_" + Amol, FRETrate]))
        print(i/n)
    #print and export results
    np.savetxt("FRET_" + database[9:], result, delimiter=",", fmt='%s')
