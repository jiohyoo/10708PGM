import sys
import pickle

fileName = 'C:/Users/HyunAh/Desktop/human_liver_cohort/expression/expression.txt'

feat_idx = 0
protein_ID = []
protein_name = []
Exp = []

for line in open(fileName,'r'):
    line = line.strip()
    if line:
        arr = line.split('\t')
        
        if feat_idx == 0:
            N_pat = len(arr) - 4
            patient_ID = arr[4:]
            feat_idx = feat_idx + 1
        else:
            protein_ID.append(arr[0])
            protein_name.append(arr[1])
            Exp.append(arr[4:])






# Write out in txt file
# fileOut = open('Exp_Out.txt', 'w')
# pickle.dump(Exp, fileOut)
# fileOut.close()

fileOut = open('Exp_Out.txt', 'w')
for x in range(0, len(Exp)):
    for y in range(0, len(Exp[x])):
        fileOut.write(str(Exp[x][y]) + " ")
    fileOut.write("\n")
fileOut.close()


fileOut = open('PatientID_Out.txt', 'w')
for x in range(0, len(patient_ID)):
    fileOut.write(str(patient_ID[x]) + "\n")
fileOut.close()


fileOut = open('ProteinID_Out.txt', 'w')
for x in range(0, len(protein_ID)):
    fileOut.write(str(protein_ID[x]) + "\n")
fileOut.close()

fileOut = open('ProteinName_Out.txt', 'w')
for x in range(0, len(protein_name)):
    fileOut.write(str(protein_name[x]) + "\n")
fileOut.close()


fileOut = open('Data_info.txt', 'w')
fileOut.write("# of patients: " + str(N_pat) + "\n")
fileOut.write("# of features: " + str(len(Exp)) + "\n")
fileOut.close()