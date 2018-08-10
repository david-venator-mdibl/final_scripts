#!/usr/bin/env python3

import numpy as np
import pandas as pd
import scipy
from scipy import stats

#inputs

table = "steinmetz_paGene.txt" 
design = "Steinmetz_design.txt"
detail = "_AvgL3pUTR"  #if name in design file differs from column titles in table, write the addendem here
control = "wt" #name the control, t-tests are run against it
column_with_difference = "genetics"

#end inputs

sample_names = []
t_test_array = []
genetics = []
Associative_Array = {}
AArray2 = {}
#AArray3 = {}

db = pd.read_csv(table, delimiter = "\t") #creates a pandas table from the text file for table
df = pd.read_csv(design, delimiter = "\t") #creates a pandas table from the text tile for design

if "paGene" in table :
  dg = db.query('UTR3count > 1')
  extra_array_name = "UTR3count"
elif "paProb" in table :
  dg = db.query('all > 0.1')
  extra_array_name = "position"

gene_array = dg["#gene"].tolist()      	
extra_array = dg[extra_array_name].tolist()    #print any added information, such as position on gene
genetic_array = df[column_with_difference].tolist()  #the differentiator, such as genetics, as it appears as a column title in the design file
sample_array = df["sample"].tolist()

for n in range(len(genetic_array)) :
  if genetic_array[n] not in genetics and genetic_array[n] != control :
    genetics.append(genetic_array[n])

genetics.append(control)

#this loop finds which samples belong to each genetic type, and puts the list of names into a .data file containing a list of samples for the type. They are stored using pickle
for s in range(0, len(genetics)):
  sample_names = []
  for i in range(0, len(genetic_array)):
    if genetic_array[i] == genetics[s]:
      sample_names.append(sample_array[i])
  Associative_Array[genetics[s]] = sample_names

for i in range(0, len(genetics)) :
   mutation = Associative_Array[genetics[i]]
   bigarray = []
   array = []
   print(mutation)
   for s in range(0, len(mutation)) :
     array = dg[mutation[s]+detail].tolist()
     AArray2[mutation[s]] = array
     with open('means_'+table, 'a+') as f:
       print(mutation[s], "	", np.mean(array), file=f)
     bigarray  = np.hstack((bigarray,array)) 
   with open('means_'+table, 'a+') as f:
     print("	", file=f)
   with open('means_'+table, 'a+') as f:
     print(genetics[i], "mean - ", np.mean(bigarray), file=f)
   with open('means_'+table, 'a+') as f:
     print("	", file=f)
   #AArray3[genetics[i]] = bigarray
   
with open('unfiltered_t-test_'+detail+'_'+table, 'a+') as f:
  print("Gene", extra_array_name, "name", "p_value", "mutant_avg", "mutant_stdev", "wt_avg", "wt_stdev", "difference", file = f) 

for n in range(0, len(gene_array)) :
  print(n) 
  t_test_array = []
  control_array = []
  for i in range(0, len(genetics)) :
    type = Associative_Array[genetics[i]]   
    if genetics[i] == control:
      control_genetic_array = []
      for s in range(0, len(type)):
        Control = AArray2[type[s]]
        control_genetic_array.append(Control[n])
      control_array.append(control_genetic_array)
    else:
      gene_genetic_array = []
      for s in range(0, len(type)):
        array = AArray2[type[s]]
        gene_genetic_array.append(array[n])
      t_test_array.append(gene_genetic_array)
  for a in range(0, len(t_test_array)) :
     genetics.remove(control)
     mutant_avg = np.mean(t_test_array[a])
     mutant_stdev = np.std(t_test_array[a])
     wt_avg = np.mean(control_array[0])
     wt_stdev = np.std(control_array[0])
     answer = scipy.stats.ttest_ind(control_array[0], t_test_array[a], equal_var = False)
     with open('unfiltered_t-test_'+detail+'_'+table, 'a+') as f:
        print(gene_array[n], extra_array[n], genetics[a], answer[1], mutant_avg, mutant_stdev, wt_avg, wt_stdev, mutant_avg-wt_avg, file=f)  
     genetics.append(control)

dv = pd.read_csv('unfiltered_t-test_'+detail+'_'+table, delimiter = " ") #creates a pandas table from the text file for table

signif = dv.query('p_value < 0.05')

dj = pd.DataFrame(signif)

dj.to_csv('filtered_t-test'+detail+'_'+table, header = True, index = True, sep = "\t")
