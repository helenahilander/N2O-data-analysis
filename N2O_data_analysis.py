#!/usr/bin/env python
# coding: utf-8

# # File for N2O data analysis
# 
# By: Helena Hilander
# 

#first in this file comes all the functions.
#Those you do not necessarily need to edit if you do not want to. But feel free to make improvements if you want! :) 


################################## NO NEED TO EDIT BELOW #####################################################################################################
#import needed packages and functions
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import datetime
import scipy.stats
import seaborn as sns
import matplotlib.dates as mdates
import warnings
from scipy.optimize import minimize

#function for reading single N2O data file
#returns N2O data time series as a data frame
def read_N2O_data(filename, plant):
    #initialize constants
    pN = 101325 
    TN = 298 
    R = 8.314 
    M_N2O = 44.013 #g/mol
    #read data
    file = pd.read_excel(filename,header = 0, index_col = None, parse_dates=[["Date", "Time"]])
    N2O_data  = file.loc[:,["Date_Time", "Line", "Compartment","Nitrous oxide"]] #fetch relevant data (can add:  "Cell temperature", "Pressure") 
    N2O_data = N2O_data.dropna(axis=0, how='any') #remove rows with no values
    N2O_data.loc[:,"Plant"] = N2O_data.shape[0]*[plant]  # add plant info
    N2O_data = N2O_data.astype({"Compartment": 'int', "Line":'int'}) #"Cell temperature":'float', 
    N2O_data.columns = ["Time", "Line", "Compartment", "N2O (ppmv)",  "Plant" ] #edit column names ("Cell temperature (K)" , "Pressure (Pa)",)
    N2O_data.index = N2O_data.loc[:,'Time'] #add timestamp as index
    N2O_data.loc[:,"Weekday"] = N2O_data.index.dayofweek #add weekday number (Monday=0, Sunday=6.)
    N2O_data = N2O_data.loc[:,["Time", "Weekday", "Plant", "Line", "Compartment", "N2O (ppmv)"]] #rearrange columns (,"Cell temperature (K)" , "Pressure (Pa)")
    #Caclulate N2O concentration (g/Nm3)
    N2O_data.loc[:,"N2O (g/Nm3)"] = M_N2O*pN/(10**6*TN*R) * N2O_data.loc[:,"N2O (ppmv)"] 
    return N2O_data


#function for reading multiple N2O data files to single dataframe:
#returns N2O data time series as a data frame
def N2O_data_to_df(path , plant):
    data_list = []
    for root, dirs, files in os.walk(path):
        for filename in files:
            if filename.endswith((".xlsx")):
                filepath = os.path.join(root, filename)
                df = read_N2O_data(filepath, plant)
                data_list.append(df)
                print("Reading file: ")
                print(filepath)
    frame = pd.concat(data_list, axis=0, ignore_index=True) #combine several files
    frame.set_index("Time", inplace=True) #Set datetime as index
    return frame


#function for reading nitrite data
#returns nitrite data time series as a data frame
def read_nitrite(filepath):
    Nitrite_data = pd.read_excel(filepath,  header = 0) #read file
    Nitrite_data.dropna(axis=0, how="any", thresh = 3, inplace = True) #remove rows with no values
    #Merge time and date cols
    Nitrite_data.loc[:,'Date'] = Nitrite_data.loc[:,'Date'].dt.date 
    for i, row in Nitrite_data.iterrows():
        Nitrite_data.loc[i,'Date_time'] = datetime.datetime.combine(Nitrite_data.loc[i,'Date'] ,  Nitrite_data.loc[i,'Time'])
    Nitrite_data.loc[:,"Compartment"] = Nitrite_data.loc[:,"Compartment"].astype('int32') #Change datatype
    Nitrite_data.set_index("Date_time", inplace=True)
    return Nitrite_data


#objective function for solving the standard cubic meters (Nm3/h) based on volumetric flow (m3/h)
def objective(x, args):
    P2 = args[0]
    P1 = args[1]
    T1 = args[2]
    V2_measured = args[3]
    
    V1 = x[0]
    T2 = x[1]
    
    V2 = V1*(P1/P2)*((T2+ 273.15)/(T1 + 273.15))
    #print("V2: ", V2)
    return abs(V2-V2_measured)

#function for converting m3/h to normal units (Nm3/h)
def find_Nm3(aeration_data):
    new_aeration_data = {}
    for key in aeration_data.keys():
        aeration_df = aeration_data[key]
        new_aeration_df = aeration_df.copy()
        standard_volumes = []
        if int(key[-1]) in [1,2,7,8]:
            dP = 54 #pressure change in the old aeration side
        else:
            dP = 64.5 #pressure change in the new aeration side
        for i, row in aeration_df.iterrows():
            
            #search for data value
            #print(row)
            V2_measured = row[0]
            #print(V2_measured)
            
            #set initial values
            P1 = 101.32 #kPa
            P2 = P1 + dP #kPa
            T2 = 64 #first guess
            T1 = 0 
            V1 = V2_measured  #first guess
            x= [V1,T2] #variables to solve
                    
            #constraints
            cons=({'type': 'ineq', 'fun': lambda x: x[0]},{'type': 'ineq','fun': lambda x: x[1]})
            
            #solve opitimization
            result = minimize(objective, x, args=([P2,P1,T1, V2_measured]), constraints = cons)
            #print(result.x)
            standard_volumes.append(float(result.x[0])) #Save standard volume
        new_aeration_df.iloc[:,0] = standard_volumes
        old_col_name = new_aeration_df.columns[0]
        new_aeration_df.rename(columns = { old_col_name: old_col_name[0:-7] + " (Nm3/h)"}, inplace = True)
        new_aeration_data[key] = new_aeration_df
        #new_aeration_df.to_excel("Aeration_Nm3_" + key.replace(" ", "") + ".xlsx") 
    return new_aeration_data


#function for reading aeration data
#returns aeration data as hourly averages for every weekday. 
#The returned data is in the form of dictionary: each key and value pair belongs to one compartment from a single line
def read_aeration_data(aeration_path, in_n_units = True, unit = "hour", plot= False, saving_path = None,  fig_name = None):
    hourly_aeration_data = {}
    #read data
    aeration_data = pd.read_excel(aeration_path ,sheet_name = None,  header = 0, index_col = None )#, parse_dates=[["Date", "Time"]]
    for line in aeration_data.keys():
        aeration_data_line = aeration_data[line]
        n_cols = aeration_data_line.shape[1]
        for col in range(0, n_cols, 2):
            aeration_data_column = aeration_data_line.iloc[:, col : col +2 ] #data from one compartment + its time column
            aeration_data_column.dropna(axis = 0, how= "any", inplace=True)
            aeration_data_column.iloc[:,0] = pd.to_datetime(aeration_data_column.iloc[:,0])
            aeration_data_column.set_index(aeration_data_column.iloc[:,0],  inplace= True)
            aeration_data_column.loc[:,"Weekday"] = aeration_data_column.index.dayofweek #add weekday
            if pd.api.types.is_string_dtype(aeration_data_column.iloc[:,1]):
                aeration_data_column.iloc[:,1] = [np.char.replace(value, ",", ".") for value in aeration_data_column.iloc[:,1]]
            aeration_data_column.iloc[:,1] = aeration_data_column.iloc[:,1].astype(float)
            #aeration_data_column.to_excel(saving_path + line.replace(" ", "")  + str(col) +  "_turku_aeration_spring.xlsx") 
            if plot:
                plot_aeration(aeration_data_column, line +" col "+ str(col), saving_path = saving_path,  fig_name = fig_name)
            #print(aeration_data_column)mean()
            hourly_aeration_data_column = aeration_data_column.groupby([aeration_data_column.index.weekday, aeration_data_column.index.hour ]).mean()
            hourly_aeration_data_column.index.rename(["Weekday", "Hours"], inplace=True)
            #print(hourly_aeration_data_column)
            column_number = hourly_aeration_data_column.columns[0][6]
            hourly_aeration_data["line " + str(line[-1]) + " compartment " + str(column_number)] = hourly_aeration_data_column
    #print(hourly_aeration_data.keys())
    
    #if aeration data is not in Nm3/h unit, convert
    if in_n_units == False:
        hourly_aeration_data= find_Nm3(hourly_aeration_data)
    return hourly_aeration_data


#function for reading inflow data with time series for each day in separate columns
#returns inflow timeseries as dataframe
def read_inflow_days_separate_cols( file, inflow_sheet , total_N_sheet, read_treated = False):
    inflow_df = read_inflow_sheet(file, inflow_sheet)
    inflow_df.columns = ["Inflow (m3/h)"]
    totalN_df = read_inflow_sheet(file, total_N_sheet)
    if read_treated:
        totalN_df.columns = ["Total N treated (gN/m3)"]
    else:
        totalN_df.columns = ["Total N (gN/m3)"]
    df_out = pd.concat( [inflow_df ,totalN_df ], axis = 1, ignore_index=False)
    return df_out


#function for reading one sheet of inflow excel file
#returns inflow timeseries as dataframe
def read_inflow_sheet(file, sheetname):
    inflow = pd.read_excel(file, sheet_name = sheetname)
    inflow.drop(labels = "Hour", axis = 1, inplace = True)
    inflow = pd.DataFrame(inflow.stack())
    inflow.reset_index(inplace = True)
    inflow.loc[:,"level_0"] = pd.to_timedelta(inflow.loc[:,"level_0"], unit="hour")
    time_str = inflow.loc[:,"level_0"].astype("string")
    timestamps = [w[7:] for w in time_str]
    inflow.loc[:,"level_1"] = inflow.loc[:,"level_1"].astype("string") + " " +  timestamps
    inflow.loc[:,"level_1"] = pd.to_datetime(inflow.loc[:,"level_1"])
    inflow.drop("level_0", axis = 1, inplace=True)
    inflow.set_index("level_1", inplace= True)
    inflow.index.rename("Time", inplace = True)
    return inflow


#function that calculates hourly mean N2O data for every weekday
#returns N2O data as hourly averages for every weekday as dataframe
#missing values are filled with average of other week/weekend days
def hourly_mean_data_fill_with_average_days(N2O_data, line_fig, compartment, plot = False, saving_path = None,  fig_name = None):
    data_list= []
    all_weekdays = [0,1,2,3,4,5,6]
    #print(N2O_data.columns)
    plant = N2O_data.loc[:,"Plant"][0]
    line = N2O_data.loc[:,"Line"][0]
    N2O_data_comp = N2O_data[N2O_data["Compartment"] == compartment]
    N2O_data_comp.drop(columns = ["N2O (ppmv)"], inplace = True)
    weekdays = np.unique(N2O_data_comp.loc[:,"Weekday"])
    missing_weekdays = list(set(all_weekdays) - set(weekdays))
    print("missing weekdays: ", missing_weekdays)
    avg_template_week = avg_template(N2O_data_comp, [0,1,2,3,4], plant)
    avg_template_weekend = avg_template(N2O_data_comp, [5,6], plant)
    #interpolate < 3 missing point gaps linearly but larger with average values from week/weekend
    for weekday in weekdays:
        N2O_data_sub = N2O_data_comp[N2O_data_comp["Weekday"] == weekday]
        if plot: 
            plot_weekday_N20(N2O_data_sub , line_fig, [compartment], saving_path, "Weekday_" + str(weekday) + "_orig_" + fig_name)
        N2O_data_inter = average_weekday_data(N2O_data_sub, weekday, line, compartment, plant)
        N2O_data_inter.loc[:,"Hours"] = list(range(0,24))
        #print(N2O_data_inter)
        ##################################################################################################
        #if more than 2 missing values, fill with weekday/weekend average values
        if weekday in [5,6]:
            N2O_data_inter.loc[:,"N2O (g/Nm3)"] = N2O_data_inter.loc[:,"N2O (g/Nm3)"].fillna(avg_template_weekend.loc[:,"N2O (g/Nm3)"])
        else:
            N2O_data_inter.loc[:,"N2O (g/Nm3)"] = N2O_data_inter.loc[:,"N2O (g/Nm3)"].fillna(avg_template_week.loc[:,"N2O (g/Nm3)"])
        #print(N2O_data_inter)
        if plot:
            plot_weekday_N20(N2O_data_inter , line_fig, [compartment], saving_path, "Weekday_" + str(weekday) + "_interpol_" + fig_name)
        data_list.append(N2O_data_inter)
    #add data from all weekdays to a single dataframe
    New_data = pd.concat(data_list, axis=0, ignore_index=True)
    #print(New_data)
    #print(New_data.columns)
    
    #add average week/weekend template for missing days
    if len(missing_weekdays) > 0:
        for day in missing_weekdays:
            if day in [5,6]: #Weekend
                #addd day to weekday template
                missing_day_data =avg_template_weekend.copy()
            else: #Weekday
                missing_day_data =avg_template_week.copy()
            missing_day_data.loc[:,"Weekday"] = np.ones(missing_day_data.shape[0])*day
            #print(missing_day_data)
            missing_day_data.loc[:,"Hours"] = list(range(0,24))
            #N2O_data_comp = N2O_data_comp.reindex(columns=
            plot_weekday_N20(missing_day_data , line_fig, [compartment], saving_path, "MISSING_Weekday_" + str(weekday) + "_interpol_" + fig_name)
            New_data = pd.concat([New_data, missing_day_data], axis=0, ignore_index=True)
            #print(New_data)
    if plot:
        plot_weekday_N20(New_data , line_fig,  [compartment], saving_path, fig_name)
    New_data.set_index(["Weekday", "Hours"], inplace=True)
    New_data.sort_values(by=['Weekday', "Hours"],inplace = True)
    #print(New_data)
    return New_data   


#calculate average hourly averages for given days
#will be used to fill missing data
def avg_template(N2O_data_comp_orig, days, plant):
    all_weekdays = [0,1,2,3,4,5,6]
    other_days = list(set(all_weekdays) - set(days))
    N2O_data_comp = N2O_data_comp_orig.copy()
    #calculate average day for given days
    avg_template_days = create_avg_template(N2O_data_comp, days, plant)
    if avg_template_days.empty:
        #if there was no data at all from the wanted weekdays, make template based on all the other days
        avg_template_other_days = create_avg_template(N2O_data_comp, other_days, plant)
        return avg_template_other_days
    #if given days are missing some values, fill with average from other weekdays
    avg_template_other_days = create_avg_template(N2O_data_comp, other_days, plant)
    avg_template_days = avg_template_days.reindex(list(range(0,24)))
    avg_template_days.fillna(avg_template_other_days , inplace=True)
    return avg_template_days
   
 
#calculate average hourly averages for given days
def create_avg_template(N2O_data_comp, days , plant):
    avg_template = N2O_data_comp[N2O_data_comp.loc[:,"Weekday"].isin(days)]
    if avg_template.empty:
        return avg_template
    avg_template =  avg_template.groupby( avg_template.index.hour ).mean()
    avg_template = avg_template.loc[:,["Line", "Compartment",  "N2O (g/Nm3)"]]
    avg_template.loc[:,"Plant"] = avg_template.shape[0]*[plant]
    return avg_template


#function that calculates hourly averages for one weekday based on N2O data
#missing values are left with NaN
def average_weekday_data(N2O_data_sub, weekday, line, compartment, plant):
    N2O_data_sub = N2O_data_sub.groupby( N2O_data_sub.index.hour).mean()
    N2O_data_sub = N2O_data_sub.reindex(list(range(0,24)))
    #print(N2O_data_sub)
    N2O_data_sub.loc[:,"Weekday"] = np.ones(N2O_data_sub.shape[0])*weekday #add correct weekday
    N2O_data_sub.loc[:,"Line"] = np.ones(N2O_data_sub.shape[0])*line #add correct weekday
    N2O_data_sub.loc[:,"Compartment"] = np.ones(N2O_data_sub.shape[0])*compartment #add correct weekday
    N2O_data_sub.loc[:,"Plant"] = N2O_data_sub.shape[0]*[plant]
    return N2O_data_sub


#function for calculating N2O_N load [gN/h]
def calculate_N2ON_load(N2O_data, aeration_data, plot= False , saving_path = None,  fig_name = None,  data_saving_path= None):
    # N2O_N_load [gN/h] = N2O_N concentration off-gas * air flow rate in aeration tank [gN/Nm3]*[Nm3/h]
    # Nm3 refers to standard cubic meters per 
    comp_data_dict = {}
    for key in aeration_data.keys():
        aeration_data_comp = aeration_data[key]
        N2O_data_comp = N2O_data[key]
        print(aeration_data_comp)
        print(N2O_data_comp)
        data_comp = pd.concat( [aeration_data_comp , N2O_data_comp], axis = 1, ignore_index=False)
        data_comp.loc[:, "N2O-N load (gN/h)"] = data_comp.iloc[:,0]*28.02/44.013 *data_comp.loc[:,"N2O (g/Nm3)"]
        data_comp.to_excel( data_saving_path + "aeration_and_N2O_" + key.replace(" ", "") + ".xlsx") 
        #print(data_comp.columns)
        if plot:
            plot_N2O_N_load(data_comp, key= key, saving_path = saving_path , fig_name = fig_name)
        comp_data_dict[key] = data_comp
    return comp_data_dict


def calculate_day_ef(N2O_N_load, inflow_data):
    inflow_data_mes = inflow_data.dropna(axis =0)
    return None
        

#calculates emission factor based on weekday hourly average N2O-load and inflow timeseries (inflow not averaged, original timeseries)
#if  calculate_to_treated_N is set = True, the emission factor is calculated as the ratio of N2O-N load and treated total nitrogen
def calculate_ef_whole_period(total_N2O_N_load_data, inflow_data,  plant, season, data_saving_path,calculate_to_treated_N = False):
    N2O_N_load_mean = total_N2O_N_load_data.loc[:,"N2O-N load (gN/h)"].mean()
    N2O_N_load_median = total_N2O_N_load_data.loc[:,"N2O-N load (gN/h)"].median()
    if calculate_to_treated_N:
        n_flow_mean, n_flow_median = calculate_avg_and_median_treated_N(inflow_data, plant, season, data_saving_path)
    else:
        n_flow_mean, n_flow_median = calculate_avg_and_median_N_inflow(inflow_data, plant, season, data_saving_path)
    ef_mean = N2O_N_load_mean/n_flow_mean*100
    ef_median = N2O_N_load_median/n_flow_median*100
    
    return ef_mean, ef_median


#calculates average and median nitrogen inflow
def calculate_avg_and_median_N_inflow(inflow_data, plant ,season, data_saving_path):
    #assume missing Total N concentrations to be average of measured days
    inflow_data_avg = inflow_data.fillna(inflow_data.mean())
    inflow_data_avg.loc[:,"Total N inflow (gN/h)"] = inflow_data_avg.loc[:,"Total N (gN/m3)"]* inflow_data_avg.loc[:,"Inflow (m3/h)"]
    inflow_data_avg.to_excel(data_saving_path + plant + "_" + season + "_filled_inflow.xlsx") 
    return inflow_data_avg.loc[:,"Total N inflow (gN/h)"].mean(), inflow_data_avg.loc[:,"Total N inflow (gN/h)"].median()


#calculates average and median treated nitrogen
def calculate_avg_and_median_treated_N(data, plant ,season, data_saving_path):
    #assume missing Total N concentrations to be average of measured days
    data_avg = data.fillna(data.mean())
    data_avg.loc[:,"Total N treated (gN/h)"] = data_avg.loc[:,"Total N treated (gN/m3)"]*data_avg.loc[:,"Inflow (m3/h)"]
    data_avg.to_excel(data_saving_path + plant + "_" + season + "_filled_treated.xlsx") 
    return data_avg.loc[:,"Total N treated (gN/h)"].mean(), data_avg.loc[:,"Total N treated (gN/h)"].median()
       

#Function for correlating hourly mean N2O to
#(obs. these hourly mean N2O values are not weekday averages but separate for each day in measurement data)
def correlate_hourly_mean_N2O_to_nitrite(N2O_data_orig, Nitrite_data_orig,line, linefig,  compartment, saving_path):
    N2O_data = N2O_data_orig.copy()
    Nitrite_data = Nitrite_data_orig.copy()
    #Select data from the right line and compartment
    Nitrite_data = Nitrite_data[ Nitrite_data["Line"] == line ]
    Nitrite_data = Nitrite_data[Nitrite_data["Compartment"] == compartment]
    N2O_data = N2O_data[N2O_data["Line"] == line]
    N2O_data = N2O_data[N2O_data["Compartment"] == compartment]
    #take hourly mean
    N2O_data = N2O_data.groupby([N2O_data.index.date, N2O_data.index.hour ]).mean()
    #index manipulation
    N2O_data.reset_index(level=[0,1], inplace=True)
    N2O_data.loc[:,"Time"] = pd.to_timedelta(N2O_data.loc[:,"Time"], unit="hour") #pd.to_datetime(N2O_data.loc[:,"Time"], format='%H')#p
    time_str = N2O_data.loc[:,"Time"].astype("string")
    timestamps = [w[7:] for w in time_str]
    N2O_data.loc[:,"level_0"] = N2O_data.loc[:,"level_0"].astype("string") + " " +  timestamps
    N2O_data.loc[:,"level_0"] = pd.to_datetime(N2O_data.loc[:,"level_0"])
    N2O_data.drop("Time", axis = 1, inplace=True)
    N2O_data.set_index("level_0", inplace= True)
    N2O_data.index.rename("Time", inplace = True)
    #initialize N2O columns to nitrite data
    Nitrite_data.loc[:,"N2O (g/Nm3)"] = np.zeros(Nitrite_data.shape[0])
    Nitrite_data.loc[:,"N2O timestamp"] = np.zeros(Nitrite_data.shape[0])
    #find closest N2O values
    for i, row in Nitrite_data.iterrows():
        timepoint = i #.to_pydatetime()
        timepoint_round = timepoint.replace(minute=0, second=0) #round time point to hour
        indexer = N2O_data.index.get_indexer([timepoint_round], method='nearest', tolerance = datetime.timedelta(seconds=3600))[0]
        if indexer != -1:
            N2O_value = N2O_data.iloc[indexer, N2O_data.columns.get_loc("N2O (g/Nm3)")]
            N2O_timestamp = N2O_data.index[indexer]
            Nitrite_data.loc[timepoint,"N2O (g/Nm3)"] = N2O_value
            Nitrite_data.loc[timepoint,"N2O timestamp"] = N2O_timestamp
        else:
            Nitrite_data = Nitrite_data.drop( timepoint, axis= 0 ) #remove row if no corrersponding N2O point found
    #calculate pearson correlation
    y = Nitrite_data.loc[:,"N2O (g/Nm3)"]
    x= Nitrite_data.loc[:,"N02- sample (mg/l)"]
    corr_coeff = scipy.stats.pearsonr(x, y)
    print("Pearson correlation at " + Nitrite_data.loc[:,"Plant"][0] + " WWTP line " + str(linefig) + " compartment " + str(compartment) + " is: \n" )
    print(corr_coeff)
    print("\n With N= ")
    print(np.shape(y))
    r_2 = corr_coeff[0]
    p_val = corr_coeff[1]
    N_val = np.shape(y)[0]
    plot_N2O_and_nitrite_correlation(Nitrite_data, compartment, linefig, saving_path, r_2, p_val, N_val)
    return Nitrite_data #return nitrite data with N2O values added  


#Function for correlating hourly mean N2O to a general variable
#(obs. these hourly mean N2O values are not weekday averages but separate for each day in measurement data)
def correlate_hourly_mean_general(N2O_data_orig, gen_data_orig, gen_data_name, hrs_time_lag_N2O_from_gen, line, linefig,  compartment, saving_path, fig_name):
    N2O_data = N2O_data_orig.copy()
    gen_data = gen_data_orig.copy()
    plant = N2O_data.loc[:,"Plant"][0]
    #Select data from the right line and compartment
    #Nitrite_data = Nitrite_data[ Nitrite_data["Line"] == line ]
    #Nitrite_data = Nitrite_data[Nitrite_data["Compartment"] == compartment]
    N2O_data = N2O_data[N2O_data["Line"] == line]
    N2O_data = N2O_data[N2O_data["Compartment"] == compartment]
    #take hourly mean
    N2O_data = N2O_data.groupby([N2O_data.index.date, N2O_data.index.hour ]).mean()
    #index manipulation
    N2O_data.reset_index(level=[0,1], inplace=True)
    N2O_data.loc[:,"Time"] = pd.to_timedelta(N2O_data.loc[:,"Time"], unit="hour") #pd.to_datetime(N2O_data.loc[:,"Time"], format='%H')#p
    time_str = N2O_data.loc[:,"Time"].astype("string")
    timestamps = [w[7:] for w in time_str]
    N2O_data.loc[:,"level_0"] = N2O_data.loc[:,"level_0"].astype("string") + " " +  timestamps
    N2O_data.loc[:,"level_0"] = pd.to_datetime(N2O_data.loc[:,"level_0"])
    N2O_data.drop("Time", axis = 1, inplace=True)
    N2O_data.set_index("level_0", inplace= True)
    N2O_data.index.rename("Time", inplace = True)
    #initialize N2O columns to nitrite data
    gen_data.loc[:,"N2O (g/Nm3)"] = np.zeros(gen_data.shape[0])
    gen_data.loc[:,"N2O timestamp"] = np.zeros(gen_data.shape[0])
    #find closest N2O values
    for i, row in gen_data.iterrows():
        timepoint = i 
        timepoint_round = timepoint.replace(minute=0, second=0) #round time point to hour
        timepoint_lag = timepoint_round + datetime.timedelta(hours= hrs_time_lag_N2O_from_gen) #add timelag to corrrelation
        indexer = N2O_data.index.get_indexer([timepoint_lag], method='nearest', tolerance = datetime.timedelta(seconds=60))[0]
        if indexer != -1:
            N2O_value = N2O_data.iloc[indexer, N2O_data.columns.get_loc("N2O (g/Nm3)")]
            N2O_timestamp = N2O_data.index[indexer]
            gen_data.loc[timepoint,"N2O (g/Nm3)"] = N2O_value
            gen_data.loc[timepoint,"N2O timestamp"] = N2O_timestamp
        else:
            gen_data = gen_data.drop( timepoint, axis= 0 ) #remove row if no corrersponding N2O point found
    #calculate pearson correlation
    y = gen_data.loc[:,"N2O (g/Nm3)"]
    x= gen_data.loc[:,gen_data_name]
    corr_coeff = scipy.stats.pearsonr(x, y)
    r_2 = corr_coeff[0]
    p_val = corr_coeff[1]
    N_val = np.shape(y)[0]
    print("Pearson correlation at " + plant + " WWTP line " + str(linefig) + " compartment " + str(compartment) + " is: \n" )
    print(corr_coeff)
    print("\n With N= ")
    print(np.shape(y))
    print("With time lag (hrs):",  hrs_time_lag_N2O_from_gen)
    plot_N2O_and_gen_correlation(gen_data, gen_data_name, hrs_time_lag_N2O_from_gen, plant,  compartment, linefig, saving_path, fig_name, r_2, p_val, N_val)
    return gen_data #return gen data with N2O values added  


#function that estimates missing compartments' N2O by linear regression from 2 compartment data.
#If you have data from more than 2 compartments, you should create a new function
def estimate_missing_compartments(first_comp_orig, sec_comp_orig, compartments_to_estimate, line,  plot= False,  saving_path =None, fig_name =None):
    first_comp = first_comp_orig.copy()
    sec_comp = sec_comp_orig.copy()
    first_comp.reset_index(inplace =True)
    sec_comp.reset_index(inplace =True)
    if "N2O (ppmv)" in first_comp.columns:
        first_comp.drop(columns = ["N2O (ppmv)"], inplace = True)
    if "N2O (ppmv)" in sec_comp.columns:  
        sec_comp.drop(columns = ["N2O (ppmv)"], inplace = True)
    weekdays = np.unique(first_comp.loc[:,"Weekday"])
    hours  = np.unique(first_comp.loc[:,"Hours"])

    first_comp_numb = int(first_comp.loc[0,"Compartment"])
    sec_comp_numb = int(sec_comp.loc[0,"Compartment"])
    first_comp_to_dic = first_comp_orig.copy()
    sec_comp_to_dic = sec_comp_orig.copy()
    if "N2O (ppmv)" in first_comp_to_dic.columns:
        first_comp_to_dic.drop(columns = ["N2O (ppmv)"], inplace = True)
    if "N2O (ppmv)" in sec_comp_to_dic.columns:
        sec_comp_to_dic.drop(columns = ["N2O (ppmv)"], inplace = True)
    
    comp_dict = {}
    comp_list = []
    if first_comp_numb < compartments_to_estimate[0]:
        comp_dict["Compartment " + str(first_comp_numb)] = first_comp_to_dic
    for comp in compartments_to_estimate:
        df_comp = first_comp_orig.copy()
        if "N2O (ppmv)" in df_comp.columns:  
            df_comp.drop(columns = ["N2O (ppmv)"], inplace = True)
        df_comp.loc[:,"Compartment"] = comp*np.ones(df_comp.shape[0])
        for day in weekdays:
            day_data1 = first_comp[first_comp.loc[:,"Weekday"] == day]
            day_data2 = sec_comp[sec_comp.loc[:,"Weekday"] == day]
            for hour in hours:
                h_data1 = day_data1[day_data1.loc[:,"Hours"] == hour]
                h_data2 = day_data2[day_data2.loc[:,"Hours"] == hour]
                x1 = h_data1.loc[:,"Compartment"]
                #print(x1.index)
                x1 = x1[x1.index[0]]
                y1 = h_data1.loc[:,"N2O (g/Nm3)"]
                y1 = y1[y1.index[0]]
                x2 = h_data2.loc[:,"Compartment"]
                x2 = x2[x2.index[0]]
                y2 = h_data2.loc[:,"N2O (g/Nm3)"]
                y2 = y2[y2.index[0]]
                k = (y2-y1)/(x2-x1)
                b = y1 - k*x1
                df_comp.loc[(day, hour),"N2O (g/Nm3)"] = max(k*comp + b,0)
        comp_dict["Compartment " + str(comp)] = df_comp
        comp_list.append(df_comp.reset_index())
        if first_comp_numb == comp +1:
            comp_dict["Compartment " + str(first_comp_numb)] = first_comp_to_dic
        if sec_comp_numb == comp+1:
            comp_dict["Compartment " + str(sec_comp_numb)] = sec_comp_to_dic
    #plot
    if plot:
        comp_list.append(first_comp)
        comp_list.append(sec_comp)
        df_plot = pd.concat(comp_list, axis=0, ignore_index=True)
        #print(df_plot)
        for day in weekdays:
            df_plot_day = df_plot[df_plot.loc[:,"Weekday"] == day]
            plot_weekday_N20_comp_to_same(df_plot_day , line, saving_path, fig_name)
    return comp_dict

#BELOW ALL FUNCTIONS ARE FOR DRAWING FIGURES

# function to convert to subscript (for plotting)
def get_sub(x):
    normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
    sub_s = "ₐ₈CDₑբGₕᵢⱼₖₗₘₙₒₚQᵣₛₜᵤᵥwₓᵧZₐ♭꜀ᑯₑբ₉ₕᵢⱼₖₗₘₙₒₚ૧ᵣₛₜᵤᵥwₓᵧ₂₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎"
    res = x.maketrans(''.join(normal), ''.join(sub_s))
    return x.translate(res)


# function to convert to superscript (for plotting)
def get_super(x):
    normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
    super_s = "ᴬᴮᶜᴰᴱᶠᴳᴴᴵᴶᴷᴸᴹᴺᴼᴾQᴿˢᵀᵁⱽᵂˣʸᶻᵃᵇᶜᵈᵉᶠᵍʰᶦʲᵏˡᵐⁿᵒᵖ۹ʳˢᵗᵘᵛʷˣʸᶻ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾"
    res = x.maketrans(''.join(normal), ''.join(super_s))
    return x.translate(res)


def plot_N2O(N2O_data, line, compartments, saving_path, fig_name ):
    #create titles and axis names
    N2O_str= 'N{}O'.format(get_sub('2')) #Create subscripted N2O for visualization
    m3_str = 'm{}'.format(get_super('3')) #Create subscripted m3 for visualization
    title = "Off-gas "+ N2O_str +" concentration at " + N2O_data.loc[:,"Plant"][0] + " WWTP line " + str(line)
    y_axis_label = "Off-gas " + N2O_str + " concentration (g/N"+ m3_str + ")" 
    x_axis_label = "Time" 
    #plot compartment data 
    fig = plt.figure(figsize=(10, 9))
    #if len(compartments) > 1:
    for compartment in compartments:
        #take subset only at given compartment
        N2O_data_sub = N2O_data[N2O_data["Compartment"] == compartment]
        #plot subset
        plt.scatter( x= N2O_data_sub.index, y = N2O_data_sub.loc[:,"N2O (g/Nm3)"],s=4, label = "Compartment " + str(compartment))
    plt.title(title, fontsize=14)
    plt.xlabel(x_axis_label, fontsize=14)
    plt.ylabel(y_axis_label, fontsize=14)
    plt.xticks(fontsize=14, rotation = 45)
    plt.yticks(fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(saving_path + fig_name) #save plot
    plt.show()
    #plot daily data
    dates = N2O_data.index.date
    dates = np.unique(dates)
    if len(dates) > 1:
        for date in dates:
            #select daily data
            N2O_data_daily = N2O_data[N2O_data.index.date == date]
            #plot
            fig = plt.figure(figsize=(10, 7))
            for compartment in compartments:
                N2O_data_sub = N2O_data_daily[N2O_data_daily["Compartment"] == compartment]
                if N2O_data_sub.shape[0] != 0:
                    plt.scatter( x= N2O_data_sub.index, y = N2O_data_sub.loc[:,"N2O (g/Nm3)"],s=4, label = "Compartment " + str(compartment))
            plt.title(title, fontsize=14)
            plt.xlabel(x_axis_label, fontsize=14)
            plt.ylabel(y_axis_label, fontsize=14)
            plt.xticks(fontsize=14, rotation = 45)
            plt.yticks(fontsize=14)
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(saving_path + str(date) + fig_name) #save plot
            plt.show()
        

def plot_weekday_N20(N2O_data_orig , line ,compartments, saving_path, fig_name):
    N2O_data = N2O_data_orig.copy() #make a copy of data so that original is not edited
    N2O_data.astype({"Weekday": int, "Compartment": int, "Line": int }) #check that datatype is correct
    #create titles and axis names
    N2O_str= 'N{}O'.format(get_sub('2')) #Create subscripted N2O for visualization
    m3_str = 'm{}'.format(get_super('3')) #Create subscripted m3 for visualization
    y_axis_label = "Off-gas " + N2O_str + " concentration (g/N"+ m3_str + ")" 
    x_axis_label = "Time" 
    title = "Hourly mean off-gas "+ N2O_str +" concentration at " + N2O_data.loc[:,"Plant"][N2O_data.index[0]] + " WWTP line " + str(line)
    weekday_labels = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
    N2O_data.fillna(0, inplace=True)
    for compartment in compartments:
        fig = plt.figure(figsize=(10, 6))
        N2O_data_comp = N2O_data[N2O_data["Compartment"] == compartment]
        weekdays = np.unique(N2O_data_comp.loc[:,"Weekday"])
        for weekday in weekdays:
            weekday = int(weekday)
            N2O_data_sub = N2O_data_comp[N2O_data_comp["Weekday"] == weekday]
            if isinstance(N2O_data_sub.index, pd.core.indexes.datetimes.DatetimeIndex):
                N2O_data_sub = N2O_data_sub.groupby( N2O_data_sub.index.hour).mean() #take mean if more than one same timestamps
                N2O_data_sub = N2O_data_sub.reindex(list(range(0,24)),fill_value=0)
                #plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            else:
                N2O_data_sub.loc[:,"Hours"] = list(range(0,24))
                N2O_data_sub.set_index("Hours", inplace =True)
            #plot subset
            plt.plot( N2O_data_sub.index, N2O_data_sub.loc[:,"N2O (g/Nm3)"],"o--", linewidth=4, label = weekday_labels[weekday])
        plt.title(title + " compartment " + str(compartment), fontsize=14)
        plt.xlabel(x_axis_label, fontsize=14)
        plt.ylabel(y_axis_label, fontsize=14)
        plt.xticks(fontsize=14, rotation = 45)
        plt.yticks(fontsize=14)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(saving_path + str(compartment)+"_"+ fig_name) #save plot
        plt.show()
        
def plot_weekday_N20_comp_to_same(N2O_data_orig , line,  saving_path, fig_name):
    N2O_data = N2O_data_orig.copy() #make a copy of data so that original is not edited
    N2O_data.astype({"Weekday": int, "Compartment": int, "Line": int }) #check that datatype is correct
    #create titles and axis names
    N2O_str= 'N{}O'.format(get_sub('2')) #Create subscripted N2O for visualization
    m3_str = 'm{}'.format(get_super('3')) #Create subscripted m3 for visualization
    y_axis_label = "Off-gas " + N2O_str + " concentration (g/N"+ m3_str + ")" 
    x_axis_label = "Time" 
    title = "Hourly mean off-gas "+ N2O_str +" concentration at " + N2O_data.loc[:,"Plant"][N2O_data.index[0]] + " WWTP line " + str(line)
    weekday_labels = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
    N2O_data.fillna(0, inplace=True)
    fig = plt.figure(figsize=(10, 6))
    compartments = np.unique(N2O_data.loc[:,"Compartment"])
    for compartment in compartments:
        N2O_data_comp = N2O_data[N2O_data["Compartment"] == compartment]
        weekdays = np.unique(N2O_data_comp.loc[:,"Weekday"])
        for weekday in weekdays:
            weekday = int(weekday)
            N2O_data_sub = N2O_data_comp[N2O_data_comp["Weekday"] == weekday]
            if isinstance(N2O_data_sub.index, pd.core.indexes.datetimes.DatetimeIndex):
                N2O_data_sub = N2O_data_sub.groupby( N2O_data_sub.index.hour).mean() #take mean if more than one same timestamps
                N2O_data_sub = N2O_data_sub.reindex(list(range(0,24)),fill_value=0)
                #plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            else:
                N2O_data_sub.loc[:,"Hours"] = list(range(0,24))
                N2O_data_sub.set_index("Hours", inplace =True)
            #plot subset
            plt.plot( N2O_data_sub.index, N2O_data_sub.loc[:,"N2O (g/Nm3)"],"o--", linewidth=4, label = "Compartment " + str(int(compartment)))
    plt.title(title + "\nWeekday: " + weekday_labels[weekday], fontsize=14)
    plt.xlabel(x_axis_label, fontsize=14)
    plt.ylabel(y_axis_label, fontsize=14)
    plt.xticks(fontsize=14, rotation = 45)
    plt.yticks(fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(saving_path + weekday_labels[weekday] +"_"+ fig_name) #save plot
    plt.show()


def plot_N2O_N_load(df_orig,  key, saving_path, fig_name):
    df = df_orig.copy()
    weekday_labels = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
    #create titles and axis names
    N2O_str= 'N{}O'.format(get_sub('2')) 
    m3_str = 'm{}'.format(get_super('3')) #Create subscripted m3 for visualization
    plant = df.loc[:,"Plant"][0][0]
    #print(plant)
    title = N2O_str + "-N load (gN/h) at " + plant + " WWTP " + key
    y_axis_label =  N2O_str + "-N load (gN/h)" #Edit: change figure y-axis name if needed
    x_axis_label = "Time (hour)" #Edit: change figure x-axis name if needed
    if "Weekday" in df.columns:
        df.drop(columns = "Weekday",  inplace=True)
    df.reset_index(inplace= True)
    weekdays = np.unique(df.loc[:,"Weekday"])
    fig = plt.figure(figsize=(10, 6))
    for weekday in weekdays:
        weekday = int(weekday)
        df_sub = df[df["Weekday"] == weekday]
        #print(df_sub.loc[:,"EF (%)"])
        #print(df_sub.loc[:,"Hours"])
        plt.plot( df_sub.loc[:,"Hours"], df_sub.loc[:,"N2O-N load (gN/h)"],"o--", linewidth=4, label = weekday_labels[weekday])
    #print(title)
    plt.title(title, fontsize=14)
    plt.xlabel(x_axis_label, fontsize=14)
    plt.ylabel(y_axis_label, fontsize=14)
    plt.xticks(fontsize=14, rotation = 45)
    plt.yticks(fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(saving_path + key.replace(" ", "") +fig_name) #save plot
    plt.show()
    
    
def plot_total_load(df_orig, saving_path, fig_name):
    df = df_orig.copy()
    weekday_labels = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
    #create titles and axis names
    plant = df.loc[:,"Plant"][0][0]
    N2O_str= 'N{}O'.format(get_sub('2')) 
    m3_str = 'm{}'.format(get_super('3')) #Create subscripted m3 for visualization
    title = N2O_str + "-N load at " + plant + " WWTP"
    y_axis_label =  N2O_str + "-N load (gN/h)" #Edit: change figure y-axis name if needed
    x_axis_label = "Time (hour)" #Edit: change figure x-axis name if needed
    df.reset_index(inplace= True)
    weekdays = np.unique(df.loc[:,"Weekday"])
    fig = plt.figure(figsize=(10, 6))
    for weekday in weekdays:
        weekday = int(weekday)
        df_sub = df[df["Weekday"] == weekday]
        #print(df_sub.loc[:,"EF (%)"])
        #print(df_sub.loc[:,"Hours"])
        plt.plot( df_sub.loc[:,"Hours"], df_sub.loc[:,"N2O-N load (gN/h)"],"o--", linewidth=4, label = weekday_labels[weekday])
    #print(title)
    plt.title(title, fontsize=14)
    plt.xlabel(x_axis_label, fontsize=14)
    plt.ylabel(y_axis_label, fontsize=14)
    plt.xticks(fontsize=14, rotation = 45)
    plt.yticks(fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(saving_path + fig_name) #save plot
    plt.show()
    

def plot_aeration(data, label, saving_path = None,  fig_name = None):
    title = "Aeration data " + label
    y_axis_label =  "Aeration" #Edit: change figure y-axis name if needed
    x_axis_label = "Time" #Edit: change figure x-axis name if needed
    plt.scatter( x= data.iloc[:,0], y= data.iloc[:,1] , s= 4)
    plt.title(title, fontsize=14)
    plt.xlabel(x_axis_label, fontsize=14)
    plt.ylabel(y_axis_label, fontsize=14)
    plt.xticks(fontsize=14, rotation = 45)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(saving_path + label.replace(" ", "") + fig_name) #save plot
    plt.show()
    

def plot_nitrite(Nitrite_data, saving_path, fig_name):
    #create titles and axis names
    NO2_str= 'NO{}-'.format(get_sub('2')) #Create subscripted NO2- for visualization
    m3_str = 'm{}'.format(get_super('3')) #Create subscripted m3 for visualization
    title = "Dissolved " + NO2_str + " concentrations at " + Nitrite_data.loc[:,"Plant"][0] + " WWTP line " + str(Nitrite_data.loc[:,"Line"][0])
    y_axis_label =  NO2_str + " concentration (mg/l)" #Edit: change figure y-axis name if needed
    x_axis_label = "Time" #Edit: change figure x-axis name if needed

    fig = plt.figure(figsize=(10, 7))
    Nitrite_data.groupby('Compartment')['N02- sample (mg/l)'].plot(legend=True, kind = "line", style = "o--", lw = 3)
    plt.title(title, fontsize=14)
    plt.xlabel(x_axis_label, fontsize=14)
    plt.ylabel(y_axis_label, fontsize=14)
    plt.xticks(fontsize=14, rotation = 45)
    plt.yticks(fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(saving_path + fig_name) #save plot
    plt.show()


# function for plotting N2O and nitrite
def plot_N2O_and_nitrite(N2O_data_orig, Nitrite_data_orig, compartment, line, saving_path ):
    N2O_data = N2O_data_orig.copy()
    Nitrite_data = Nitrite_data_orig.copy()
    #take subset only at given compartment
    N2O_data = N2O_data[N2O_data["Compartment"] == compartment]
    Nitrite_data = Nitrite_data[Nitrite_data["Compartment"] == compartment]
    #create titles and axis names
    N2O_str= 'N{}O'.format(get_sub('2')) #Create subscripted N2O for visualization
    NO2_str= 'NO{}-'.format(get_sub('2')) #Create subscripted NO2- for visualization
    m3_str = 'm{}'.format(get_super('3')) #Create subscripted m3 for visualization
    title = "Off-gas "+ N2O_str +" and dissolved " + NO2_str + " concentrations at " + N2O_data.loc[:,"Plant"][0] + " WWTP line " + str(line) + " compartment " + str(compartment)
    N2O_y_axis_label = "Off-gas " + N2O_str + " concentration (g/N"+ m3_str + ")" 
    Nitrite_y_axis_label = NO2_str + " concentration (g/"+ m3_str + ")" 
    x_axis_label = "Time" 
    #plot
    fig, ax1 = plt.subplots( figsize=(10,8))
    ax1.plot( N2O_data.index, N2O_data.loc[:,"N2O (g/Nm3)"], "b.", label = "Off-gas " + N2O_str)
    ax1.set_xlabel( x_axis_label)
    plt.xticks(fontsize=14, rotation = 45)
    #mn, mx = ax1.set_ylim(mean-amp, mean+amp)
    ax1.set_ylabel(N2O_y_axis_label, fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc='upper left' )
    plt.tight_layout()
    ax2 = ax1.twinx()
    #ax2.set_ylim(mn*km3yearToSv, mx*km3yearToSv)
    ax2.set_ylabel(Nitrite_y_axis_label , fontsize=14)
    ax2.set_xlabel( x_axis_label)
    plt.yticks(fontsize=14) 
    ax2.plot( Nitrite_data.index, Nitrite_data.loc[:,'N02- sample (mg/l)'] , "ro--", ms = 4, label = NO2_str)
    plt.title(title, fontsize=14)
    plt.legend(loc='upper right' )
    plt.tight_layout()
    plt.savefig(saving_path) #save plot 
    plt.show()
    

def plot_N2O_and_flow(N2O_data_orig, flow_data_orig, compartment, saving_path ):
    N2O_data = N2O_data_orig.copy()
    flow_data = flow_data_orig.copy()
    #take subset only at given compartment
    N2O_data = N2O_data[N2O_data["Compartment"] == compartment]
    #Nitrite_data = Nitrite_data[Nitrite_data["Compartment"] == compartment]
    #create titles and axis names
    N2O_str= 'N{}O'.format(get_sub('2')) #Create subscripted N2O for visualization
    #NO2_str= 'NO{}-'.format(get_sub('2')) #Create subscripted NO2- for visualization
    m3_str = 'm{}'.format(get_super('3')) #Create subscripted m3 for visualization
    title = "Off-gas "+ N2O_str +" concentration and inflow at " + N2O_data.loc[:,"Plant"][0] + " WWTP"
    N2O_y_axis_label = "Off-gas " + N2O_str + " concentration (g/N"+ m3_str + ")" 
    flow_y_axis_label = "Inflow (" +m3_str+  "/h)"
    x_axis_label = "Time" 
    #plot
    fig, ax1 = plt.subplots( figsize=(12,8))
    ax1.plot( N2O_data.index, N2O_data.loc[:,"N2O (g/Nm3)"], "b.", label = "Off-gas " + N2O_str )
    ax1.set_xlabel( x_axis_label)
    plt.xticks(fontsize=14, rotation = 45)
    #mn, mx = ax1.set_ylim(mean-amp, mean+amp)
    ax1.set_ylabel(N2O_y_axis_label, fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc='upper left' )
    ax2 = ax1.twinx()
    #ax2.set_ylim(mn*km3yearToSv, mx*km3yearToSv)
    ax2.set_ylabel(flow_y_axis_label , fontsize=14)
    plt.yticks(fontsize=14) 
    ax2.plot( flow_data.index, flow_data.loc[:,"Inflow (m3/h)"] , "ro--", ms = 4, label ="Inflow" )
    plt.title(title, fontsize=14)
    plt.legend(loc='upper right' )
    plt.tight_layout()
    plt.savefig(saving_path) #save plot 
    plt.show()


#function for plotting correlation line
def plot_N2O_and_nitrite_correlation(Nitrite_data, compartment, line, saving_path,  r_2, p_val, N_val):
    #create titles and axis names
    N2O_str= 'N{}O'.format(get_sub('2')) #Create subscripted N2O for visualization
    m3_str = 'm{}'.format(get_super('3')) #Create subscripted m3 for visualization
    NO2_str= 'NO{}-'.format(get_sub('2')) #Create subscripted NO2- for visualization
    title = "Off-gas "+ N2O_str +" and dissolved " + NO2_str + " correlation at \n" + Nitrite_data.loc[:,"Plant"][0] + " WWTP line " + str(line) + " compartment " + str(compartment)
    y_axis_label = "Off-gas " + N2O_str + " concentration (g/N"+ m3_str + ")" 
    x_axis_label = NO2_str + " concentration (g/"+ m3_str + ")" 
    #plot  
    #fig = plt.figure(figsize=(20, 10))
    g= sns.lmplot(x= "N02- sample (mg/l)", y= "N2O (g/Nm3)", data = Nitrite_data , fit_reg =True)
    title_add = "\n Pearson correlation = {:.2f}, p= {:.3f}, N = {:}".format(r_2, p_val, N_val)
    plt.title(title + title_add, fontsize=14)
    plt.xlabel(x_axis_label, fontsize=14)
    plt.ylabel(y_axis_label, fontsize=14)
    plt.xticks(fontsize=14, rotation = 45)
    plt.yticks(fontsize=14)
    #plt.legend()
    #plt.grid(True)
    plt.savefig(saving_path,bbox_inches="tight", pad_inches = 0.5) #save plot     
    plt.show()

#function for plotting correlation line for all data from several different plants
def plot_N2O_and_nitrite_correlation_all_data(Nitrite_data, saving_path,  r_2, p_val, N_val):
    #create titles and axis names
    N2O_str= 'N{}O'.format(get_sub('2')) #Create subscripted N2O for visualization
    m3_str = 'm{}'.format(get_super('3')) #Create subscripted m3 for visualization
    NO2_str= 'NO{}-'.format(get_sub('2')) #Create subscripted NO2- for visualization
    title = "Off-gas "+ N2O_str +" and dissolved " + NO2_str + " correlation"
    y_axis_label = "Off-gas " + N2O_str + " concentration (g/N"+ m3_str + ")" 
    x_axis_label = NO2_str + " concentration (g/"+ m3_str + ")" 
    #plot  
    #fig = plt.figure(figsize=(20, 10))
    g= sns.lmplot(x= "N02- sample (mg/l)", y= "N2O (g/Nm3)", data = Nitrite_data , fit_reg =True)
    title_add = "\n Pearson correlation = {:.2f}, p= {:.3f}, N = {:}".format(r_2, p_val, N_val)
    plt.title(title + title_add, fontsize=14)
    plt.xlabel(x_axis_label, fontsize=14)
    plt.ylabel(y_axis_label, fontsize=14)
    plt.xticks(fontsize=14, rotation = 45)
    plt.yticks(fontsize=14)
    #plt.legend()
    #plt.grid(True)
    plt.savefig(saving_path,bbox_inches="tight", pad_inches = 0.5) #save plot     
    plt.show()


#function for plotting correlation line
def plot_N2O_and_gen_correlation(gen_data, gen_data_name, time_lag, plant,  compartment, line,saving_path, fig_name, r_2, p_val, N_val):
    #create titles and axis names
    N2O_str= 'N{}O'.format(get_sub('2')) #Create subscripted N2O for visualization
    m3_str = 'm{}'.format(get_super('3')) #Create subscripted m3 for visualization
    #NO2_str= 'NO{}-'.format(get_sub('2')) #Create subscripted NO2- for visualization
    title = "Off-gas "+ N2O_str +" and " + gen_data_name + " correlation at \n" + plant + " WWTP" 
    if time_lag !=0:
        title_add1 =  " with time lag of " + str(time_lag) + " hours"
    else:
        title_add1 = ""
    y_axis_label = "Off-gas " + N2O_str + " concentration (g/N"+ m3_str + ")" 
    x_axis_label = gen_data_name
    #plot  
    #fig = plt.figure(figsize=(20, 10))
    g= sns.lmplot(x= gen_data_name, y= "N2O (g/Nm3)", data = gen_data , fit_reg =True)
    #ax = g.axes[0,0]
    title_add2 = "\n Pearson correlation = {:.2f}, p = {:.3f}, N = {:}".format(r_2, p_val, N_val)
    plt.title(title + title_add1 + title_add2, fontsize=14)
    plt.xlabel(x_axis_label, fontsize=14)
    plt.ylabel(y_axis_label, fontsize=14)
    plt.xticks(fontsize=14, rotation = 45)
    plt.yticks(fontsize=14)
    #ax.annotate(, xy= (x_coord,y_coord) )
    #plt.legend()
    #plt.grid(True)
    plt.savefig(saving_path + "lag_hrs_" + str(time_lag) + fig_name,bbox_inches="tight", pad_inches = 0.5) #save plot     
    plt.show()
################################## NO NEED TO EDIT ABOVE #####################################################################################################


#This is the main function that calls all the functions above.
def main():
    warnings.filterwarnings("ignore")
    ################## EDIT BELOW ####################################################################################################################
    
    #each measurement period is its own element in all the lists. For example 0 element in all lists includes information on first measurement period.
    
    #ADD INFO ABOUT EACH MEASUREMENT CAMPAIGN
    plants = ["Kakolanmäki", "Nenäinniemi", "Kakolanmäki", "Nenäinniemi"] #plant names
    seasons = ["spring", "spring", "summer", "summer"] #season or month when measurement was conducted
    measured_line_fig_names = [2, "D", 2, "D",] #line name for figures. Give name of line in which measurements were conducted
    measured_line_numbers = [2,4 ,2,4] #line number since numbers are easier in data analysis than strings.
    compartments = [[6], [1,8], [4,6],[4,8]] #give all compartments from which there is data
    all_line_names = [[1,2,3,4], ["A","B","C","D"],[1,2,3,4], ["A","B","C","D"]] #give names for all aeration lines in the plants
    aerated_compartments = [[3,4,5,6], [1,2,3,4,5,6,7,8],[3,4,5,6], [3,4,5,6,7,8]] #give numbers for all aerated compartments
    inflow_sep_cols = [True, False, False, False] #Insert True if inflow data has separate columns for each day, False if all inflow data is in one column
    aeration_in_N_units = [True, False, True, False] #Insert True if aearation data is in units Nm3/h, false if unit is m3/h.
    
    #ADD DATA PATHS FOR EACH MEASUREMENT CAMPAIGN
    N2O_paths =  [  ] #give FOLDER in which N2O data subfolders exist
    fig_saving_paths = [  ] #give FOLDER to which save figures
    output_data_paths = [   ] #give FOLDER where to save analyzed data
    nitrite_paths = [ ] #give nitrite data FILE
    aeration_data_paths = [] #give aeration data FILE
    inflow_paths = [ ] #give inflow data FILE
    
    #SELECT CAMPAIGNS TO ANALYZE
    analyze_paths = [0,1,2,3] #give measurement periods to analyze (0=first)
    
    #uncomment line below if you want to print full dataframes
    #pd.set_option("display.max_rows", None, "display.max_columns", None)
    
    ################ EDIT ABOVE #################################################################################################################################
    
    #THIS LOOP DOES ALL THE ANALYSIS. NO NEED TO EDIT IF YOU DO NOT WANT TO. 
    
    #initialize list to store all nitrite data 
    nitrite_list = []
 
    
    for i in analyze_paths:
        #define variables
        plant = plants[i]
        season = seasons[i]
        line_fig_name = measured_line_fig_names[i]
        line_number = measured_line_numbers[i]
        compartment_list = compartments[i]
        all_line_names_list = all_line_names[i]
        aerated_compartments_list = aerated_compartments[i]
        
        #1. read, plot and correlate N2O and nitrite data
        N2O_data = N2O_data_to_df(path = N2O_paths[i] , plant = plant) #read N2O data (values recorded ~ every minute)
        plot_N2O(N2O_data, line = line_fig_name, compartments = compartment_list, saving_path=fig_saving_paths[i] , fig_name = "N2O_timeseries_" + plant + "_" +season + ".pdf")
        if nitrite_paths[i] != None:
            nitrite_data = read_nitrite(nitrite_paths[i]) #read nitrite data
        for compartment in compartment_list:
            if nitrite_paths[i] != None:
                plot_N2O_and_nitrite(N2O_data, nitrite_data ,compartment= compartment, line= line_fig_name, saving_path= fig_saving_paths[i]+ "N2O_Nitrite_compartment_" + str(compartment) + "_" + plant+ "_" +season+ ".pdf" )
                #correlate hourly mean N2O and ntirite samples
                nitrite_N2O_df = correlate_hourly_mean_N2O_to_nitrite(N2O_data, nitrite_data , line= line_number, linefig=line_fig_name ,compartment= compartment,  saving_path= fig_saving_paths[i]+ "N2O_Nitrite_correlation_compartment_" + str(compartment) + "_" + plant + "_" +season + ".pdf" )
                nitrite_list.append(nitrite_N2O_df)
            print("Compartment {:} N2O data shape: {:}".format(compartment, N2O_data[N2O_data.loc[:,"Compartment"] == compartment].shape))
     
        #2. Calculate hourly average N2O data
        #calculate hourly average N2O for every day of week
        hourly_N2O_list = []
        for compartment in compartment_list:
            hourly_N2O = hourly_mean_data_fill_with_average_days(N2O_data, line_fig = line_fig_name , compartment= compartment, plot = True, saving_path = fig_saving_paths[i],  fig_name = "Weekday_hourly_avg_N2O_compartment_" + str(compartment) + ".pdf" )
            hourly_N2O.describe().to_csv(output_data_paths[i] + plant + "_compartment_" + str(compartment) +"_" + seasons[i] + "_SUMMARY_hourly_N2O.txt", header = True, index = True, sep = " ", mode= "w")
            hourly_N2O.to_excel(output_data_paths[i] + plant + "_compartment_" + str(compartment) +"_" + seasons[i] + "_hourly_N2O.xlsx") 
            #print(hourly_N2O)
            hourly_N2O_list.append(hourly_N2O)
        
        #3. Read inflow data
        #read hourly inflow data from the whole measuring period (2 weeks)
        if inflow_sep_cols[i]:
            inflow_data = read_inflow_days_separate_cols( inflow_paths[i], inflow_sheet = "Inflow" , total_N_sheet = "Total N")
            treated_data = read_inflow_days_separate_cols( inflow_paths[i], inflow_sheet = "Inflow" , total_N_sheet = "Total N treated", read_treated = True) #read treatment results separately
        else:
            inflow_data = pd.read_excel(inflow_paths[i],header = 0, index_col = 0)
            treated_data = inflow_data  #treatment results are stored in the same dataframe
        #print(inflow_data)
        
        #4. Plot and correlate  N2O and inflow
        #correlation done for the whole measuring period for hourly averages
        for compartment in compartment_list:
            plot_N2O_and_flow(N2O_data , inflow_data , compartment = compartment ,saving_path= fig_saving_paths[i] +"N2O_inflow_compartment_" + str(compartment) + "_" + plant+ "_" +season+ ".pdf"  )
            #correlate inflow flow and N2O
            for lag in [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]: #study all these timelags
                correlation_data = correlate_hourly_mean_general(N2O_data, inflow_data, "Inflow (m3/h)",hrs_time_lag_N2O_from_gen= lag, line = line_number, linefig = line_fig_name, compartment = compartment, saving_path = fig_saving_paths[i] , fig_name = "N2O_inflow_correlation_compartment_" + str(compartment) + "_" + plant+ "_" +season+ ".jpg" )
                #print(correlation_data)
        
        #5. read aeration data
        #read aeration from whole measurement period (data ~ every X mins)
        aeration_data = read_aeration_data(aeration_data_paths[i], in_n_units =  aeration_in_N_units[i] , unit = "hour", plot= True , saving_path = fig_saving_paths[i],  fig_name = "aeration_" + plant +"_" + season + ".pdf")
        #print(aeration_data)
    
        #6. calculate N2O-N load
        hourly_N2O_all_lines = {} #initialize for data storage
        if len(compartment_list) >1:
            #if there is data from 2 comaprtments estimate missing compartments with linear regression
            compartments_estimate = list(set(aerated_compartments_list) - set(compartment_list))
            hourly_N2O_all_comp = estimate_missing_compartments(hourly_N2O_list[0], hourly_N2O_list[1], line = line_fig_name, compartments_to_estimate = compartments_estimate, plot= True,  saving_path = fig_saving_paths[i] , fig_name = "missing_comp_linear_regression_Weekday_hourly_avg_N2O.pdf")
            for line in all_line_names_list: 
                #assume all lines are equal. Compartments estimated linearly above
                for key in hourly_N2O_all_comp.keys():
                    label = "line " + str(line) + " compartment " + key[-1]
                    hourly_N2O_all_lines[label] = hourly_N2O_all_comp[key]
        else:
            #if data only from 1 compartment, assume all compartments' N2O concentration equal
            for key in aeration_data.keys(): #loop all compartments and lines in plant 
                hourly_N2O_all_lines[key] = hourly_N2O_list[0] 
        #calculate N2O-N load for all lines and compartments based on weekday hourly average data
        N2O_N_load = calculate_N2ON_load(hourly_N2O_all_lines, aeration_data, plot= True , saving_path = fig_saving_paths[i],  fig_name = "N2O_N_load_" + "_weekday_hourly_avg_" + plant +"_" + season + ".pdf", data_saving_path = output_data_paths[i])
        #calculate total N2O-load of the plant
        total_N2O_N_load = list(N2O_N_load.items())[0][1].copy() #initialize variable total N2O-N load
        total_N2O_N_load = total_N2O_N_load.loc[:,["Plant", "N2O-N load (gN/h)"]] #initialize variable total N2O-N load
        total_N2O_N_load.loc[:,"N2O-N load (gN/h)"] = np.zeros(total_N2O_N_load.shape[0]) #initialize variable total N2O-N load
        for key in N2O_N_load: #loop all compartments and lines in plant 
            N2O_N_load_single_compartment = N2O_N_load[key]
            total_N2O_N_load.loc[:,"N2O-N load (gN/h)"] = total_N2O_N_load.loc[:,"N2O-N load (gN/h)"] + N2O_N_load_single_compartment.loc[:,"N2O-N load (gN/h)"] #sum N2O-N loads from all compartments from all lines
        #print(total_N_load)
        plot_total_load(total_N2O_N_load, saving_path = fig_saving_paths[i],  fig_name = "TOTAL_N2O_N_load_" + "_weekday_hourly_avg_" + plant +"_" + season + ".pdf")
        total_N2O_N_load.to_excel(output_data_paths[i] + plant + "_" + seasons[i] + "_weekday_hourly_N2O_N_load.xlsx") 
        total_N2O_N_load.describe().to_csv(output_data_paths[i] + plant + "_" + seasons[i] + "_SUMMARY_weekday_hourly_N2O_N_load.txt", header = True, index = True, sep = " ", mode= "w")
        
        #7. Calculate emission factor (EF) as a fraction of incoming nitrogen emitted as N2O
        ef_mean, ef_median = calculate_ef_whole_period(total_N2O_N_load, inflow_data, plant, season, output_data_paths[i])
        print("Average EF at " + plant + " during " + season + " campaign: ")
        print(ef_mean)
        print("\n Median EF at " + plant + " during " + season +" campaign: ")
        print(ef_median)       
        ef_df = pd.DataFrame({"Plant": plant, "Season": season, "EF (%) mean": ef_mean, "EF (%) median": ef_median} , index=[0])
        ef_df.to_csv(output_data_paths[i] + plant + "_" + seasons[i] + "_EF.txt", header = True, index = None, sep = " ", mode= "w") 
    
        #8. Calculate emission factor as the ratio of emitted N2O to treated total nitrogen
        ef_mean, ef_median = calculate_ef_whole_period(total_N2O_N_load, treated_data, plant, season, output_data_paths[i], calculate_to_treated_N = True)
        print("Average EF to treated N at " + plant + " during " + season + " campaign: ")
        print(ef_mean)
        print("\n Median EF to treated N at " + plant + " during " + season +" campaign: ")
        print(ef_median)       
        ef_df = pd.DataFrame({"Plant": plant, "Season": season, "EF treated (%) mean": ef_mean, "EF treated (%) median": ef_median} , index=[0])
        ef_df.to_csv(output_data_paths[i] + plant + "_" + seasons[i] + "_EF_treated.txt", header = True, index = None, sep = " ", mode= "w")
       
    
     

main()

