#!/usr/bin/env python
# coding: utf-8
#-------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
from xlsxwriter import Workbook
import os, glob, sys, datetime, scipy
import argparse
#-------------------------------------------------------------------------------
# Author: Aidan Shands
# USAGE
# python Calculate_EC50.py -i input.csv -v 0,5,25,100
# Outputs:
# EC50.Results.date.xlsx
# EC50.Summary.date.xlsx
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input file')
parser.add_argument('-v', '--values', help='Concentration values separated by commas (4 max and include 0 first)')
args = parser.parse_args()
input_file = args.input
value_list = args.values.split(',')
dt = datetime.datetime.now()
date = dt.strftime("%m%d%y")
#-------------------------------------------------------------------------------
df = pd.read_csv(input_file)
df["Average"] = ((df[['M1', 'M2']].mean(axis=1)-5.3)/2)
cols = ["Isolate", "Replicate"]
df[cols] = df[cols].astype(str)
df["Isolate_w_Rep"] = df["Isolate"] + '-' + df["Replicate"]
cols2drop = ["Day", 'M1', 'M2', 'Isolate', 'Replicate']
Isolates = df['Isolate_w_Rep'].unique().tolist()
df.drop(cols2drop, axis=1, inplace = True)
concentrations = value_list
#-------------------------------------------------------------------------------
writer = pd.ExcelWriter('EC50.Results.'+date+'.xlsx', engine = 'xlsxwriter')
writer2 = pd.ExcelWriter('EC50.Summary.'+date+'.xlsx', engine = 'xlsxwriter')
Summary = []
#-------------------------------------------------------------------------------
for i in Isolates:
    df3 = df[(df['Isolate_w_Rep'] == i)]
    Control = df3['Average'].values[0]
    PI5= 100-df3['Average'].values[1]/Control*100
    PI25= 100-df3['Average'].values[2]/Control*100
    PI100= 100-df3['Average'].values[3]/Control*100
    PI_List = [0, PI5, PI25,PI100]
    Dict = dict(zip(concentrations, PI_List))
    PI_DF = pd.DataFrame.from_dict(Dict, orient='index')
    PI_DF.reset_index(inplace = True)
    PI_DF.columns =["Concentration", "Percent_Inhibition"]
    PI_DF[["Concentration", "Percent_Inhibition"]] = PI_DF[["Concentration", "Percent_Inhibition"]].apply(pd.to_numeric)
    PI_DF['Decimal']=PI_DF['Percent_Inhibition']/100
    PI_DF = PI_DF.iloc[1:]
    PI_DF.drop(PI_DF[(PI_DF['Percent_Inhibition'] < 0) | (PI_DF['Percent_Inhibition'] > 100) ].index, inplace = True)
    PI_DF['LnLOGIT']=-(np.log((1/PI_DF['Decimal'])-1))
    concentrations2 = PI_DF['Concentration'].unique().tolist()
    Con_LnLOGIT = []
    for conc in concentrations2:
        Value = np.log(conc)
        Con_LnLOGIT.append(Value)
    PI_DF["LnLOGIT_Conc"] = Con_LnLOGIT
    x = Con_LnLOGIT
    y = PI_DF.LnLOGIT
    curve = np.polyfit(x, y, 1)
    p = np.poly1d(curve)
    plt.figure()
    plt.scatter(x, y)
    plt.plot(x, p(x), "r--")
    plt.title(i+ " % Inhibition and Concentration")
    EQ = "y=%.2fx-%.2f"%(curve[0],np.absolute(curve[1]))
    EC50 = np.exp(((0-curve[1]))/curve[0])
    EC = "EC50=" + str(EC50)
    res = stats.linregress(x, y)
    R2 = f"R-squared: {res.rvalue**2:.3f}"
    Sum_Data = [str(i), EC50]
    Summary.append(Sum_Data)
    plt.savefig(i+'_EC50_plot.png')
    PI_DF.to_excel(writer, sheet_name=i)
    worksheet = writer.sheets[i]
    df3.to_excel(writer, sheet_name=i, startcol=0 , startrow = 7)
    worksheet.write(0, 7, EQ)
    worksheet.write(1, 7, EC)
    worksheet.write(2, 7, R2)
    worksheet.insert_image('J2',i+'_EC50_plot.png')
Summary_DF = pd.DataFrame(Summary, columns = ['Isolate', 'EC50'], index=None)
Summary_DF.to_excel(writer2, sheet_name="Summarized_Data")
writer.save()
writer2.save()
#-------------------------------------------------------------------------------
current_directory = os.getcwd()
files = glob.glob(current_directory+'/*.png')
for f in files:
    try:
        os.remove(f)
    except OSError as e:
        print("Error: %s : %s" % (f, e.strerror))
