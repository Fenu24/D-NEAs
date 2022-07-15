import numpy as np
import pandas as pd
import csv
import json # JavaScript Object Notation
import requests
import os.path
from os import path
import csv

# Reading JPL SBDB (database of NEOs with defined A2)
jpl = pd.read_csv('jpl.csv', sep=',', keep_default_na=False, na_values=[''], quoting = csv.QUOTE_NONNUMERIC)

jpl_spkid=np.array(jpl['spkid']).astype(int)

jpl_number=[]
for i in range(len(jpl_spkid)): # objects with spkid>3000000 are unnumbered
    if jpl_spkid[i]<3000000:
        jpl_number.append(str(jpl_spkid[i]-2000000))
    else:
        jpl_number.append("NaN")
jpl_number=np.array(jpl_number)


# np.arrays of data
jpl_fn=np.array(jpl['full_name']).astype(str)
jpl_name=np.array(jpl['name']).astype(str)
jpl_H=np.array(jpl['H']).astype(float)
jpl_D=np.array(jpl['diameter'])
jpl_albedo=np.array(jpl['albedo'])
jpl_P=np.array(jpl['rot_per'])
jpl_BV=np.array(jpl['BV'])
jpl_UB=np.array(jpl['UB'])
jpl_IR=np.array(jpl['IR'])
jpl_H_sigma=np.array(jpl['H_sigma'])
jpl_D_sigma=np.array(jpl['diameter_sigma'])
jpl_e=np.array(jpl['e'])
jpl_a=np.array(jpl['a'])
jpl_i=np.array(jpl['i'])
jpl_A2=np.array(jpl['A2'])
jpl_n=np.array(jpl['n'])
jpl_a=np.array(jpl['a'])
jpl_e=np.array(jpl['e'])


# taking identification from the full name (e.g. taking '1988 XB' from '  7753 (1988 XB)')
for i in range(len(jpl_fn)):
    fn=jpl_fn[i]
    start=fn.find('(')+1
    stop=fn.find(')')
    if jpl_name[i]=='nan':
        jpl_name[i]=fn[start:stop]
        
#deleting spaces from IDs (in order to be compared with Marco's database)
for i in range(len(jpl_fn)):
    jpl_name[i]=jpl_name[i].replace(' ','')

# *****************************************************************************
# retrieving sigma_A2 from jpl SBDB (takes few minutes)
    
# This procedure lasts approx 1 second per object. It is executed only if file
# jpl_A2_sigma.txt does not exist already
if path.exists("jpl_A2_sigma.txt"): 
    jpl_A2_sigma=np.loadtxt("jpl_A2_sigma.txt")
    
else: # API
     
    jpl_A2_sigma=[]
    for i in range(len(jpl)):
        print(i)
        spikd=int(jpl['spkid'][jpl.index[i]]) 
        page='https://ssd-api.jpl.nasa.gov/sbdb.api?spk='+str(spikd) #page for specific object
        data = requests.get(page).json() # getting data from the page (JavaScript Object Notation)
        
        for j in range(len(data['orbit']['model_pars'])): # searching model_pars section
            if data['orbit']['model_pars'][j]['name']=='A2': # searching A2 inside model_pars
                ind=j # index which corresponds to A2 (here is also A1 and other model parameters)
        
        # getting only A2_sigma (A2 is already taken directly from JPL SBDB)
        jpl_A2_sigma.append(data['orbit']['model_pars'][ind]['sigma']) 
    
    # converting to np.array
    jpl_A2_sigma=np.array(jpl_A2_sigma).astype(float)
    
    np.savetxt("jpl_A2_sigma.txt", jpl_A2_sigma)
# *****************************************************************************

# Calculating dadt and dadt_sigma
jpl_dadt=2*jpl_A2/jpl_a**2/(1-jpl_e**2)/np.deg2rad(jpl_n)*365250000
jpl_dadt_sigma=2*jpl_A2_sigma/jpl_a**2/(1-jpl_e**2)/np.deg2rad(jpl_n)*365250000


# *****************************************************************************
# reading Marcos's database
marco=pd.read_csv('YarkoDB.csv', sep=',', keep_default_na=False, na_values=[''], quoting = csv.QUOTE_NONNUMERIC)
marco_ref=np.array(marco['Ref']).astype(str)

# deleting entries from Marco's base which have data from jpl
remove=[]
for i in range(len(marco)):
    if marco_ref[i]=='D':
        remove.append(marco.index[i])
    
marco=marco.drop(remove)


                 
marco_number=np.array(marco['ID']).astype(str)
marco_name=np.array(marco['Name']).astype(str)
marco_snr=np.array(marco['SNR']).astype(float)
        
# checking if there is more then one entry for same object in Marco's base
remove=[]
for i in range(len(marco_number)):
    for j in range(i,len(marco_number)):
        if marco_name[i]==marco_name[j] and i!=j:
            if marco_snr[i]>=marco_snr[j]:
                remove.append(marco.index[j])
            else:
                remove.append(marco.index[i])

# deleting less relevant entries
marco=marco.drop(remove)

# taking data again after deleting entries

marco_number=np.array(marco['ID']).astype(str)
marco_number[np.where(marco_number=='0')]=''
marco_name=np.array(marco['Name']).astype(str)
marco_dadt=np.array(marco['dadt']).astype(float)
marco_dadt_sigma=np.array(marco['dadt_std']).astype(float)

marco_a=np.array(marco['a']).astype(float)
marco_e=np.array(marco['e']).astype(float)
marco_i=np.array(marco['i']).astype(float)
marco_H=np.array(marco['H']).astype(float)
marco_H_sigma=np.array(marco['H_std']).astype(float)
marco_D=np.array(marco['D']).astype(float)
marco_D_sigma=np.array(marco['D_std']).astype(float)
marco_P=np.array(marco['P']).astype(float)
marco_P_sigma=np.array(marco['P_std']).astype(float)
marco_ref=np.array(marco['Ref']).astype(str)

# Comparing JPL and cleared Marco's database
new=[] # all objects in "new" are not in Marco's base
exist_jpl=[]
exist_marco=[]
for i in range(len(jpl_number)):
    not_exist=-1 # not in Marco's database
    for j in range(len(marco_number)):
        if (jpl_number[i]==marco_number[j] and jpl_number[i]!='NaN') or jpl_name[i]==marco_name[j]: # checking if the object is in Marco's database
            not_exist=i
            exist=j

    if not_exist==-1: # object is NOT in Marco's database
        new.append(i) 
    else: # object IS in Marco's database
        exist_jpl.append(i) # ordered number on JPL database
        exist_marco.append(exist) # ordered number of the same object in Marco's database

# converting to np.array         
exist_marco=np.array(exist_marco)
exist_jpl=np.array(exist_jpl)
   
# Merging JPL and Marco's database (Marco first, JPL after)     
number=np.concatenate((marco_number, jpl_number[new]))
name=np.concatenate((marco_name, jpl_name[new]))
a=np.concatenate((marco_a, jpl_a[new]))
e=np.concatenate((marco_e, jpl_e[new]))
inc=np.concatenate((marco_i, jpl_i[new]))
H=np.concatenate((marco_H, jpl_H[new]))
H_sigma=np.concatenate((marco_H_sigma, jpl_H_sigma[new]))
D=np.concatenate((marco_D, jpl_D[new]))
D_sigma=np.concatenate((marco_D_sigma, jpl_D_sigma[new]))
P=np.concatenate((marco_P, jpl_P[new]))

# JPL does not provide P_sigma
P_sigma=np.empty(len(number))
P_sigma[:] = np.NaN
P_sigma[:len(marco_number)]=marco_P_sigma


# Initializing arrays of NaNs
dadt_jpl=np.array([np.nan for _ in range(len(number))]) # from JPL (some objects have data frm both JPL and Marco's database)
dadt_papers=np.array([np.nan for _ in range(len(number))]) # from Marco's database

dadt_sigma_jpl=np.array([np.nan for _ in range(len(number))]) # from JPL (some objects have data frm both JPL and Marco's database)
dadt_sigma_papers=np.array([np.nan for _ in range(len(number))]) # from Marco's database


for i in range(len(marco_number)):
    if i in exist_marco: # if object from JPL is also in Marco's database
        ind=np.where(i==exist_marco)[0][0] # index in JPL
        dadt_jpl[i]=jpl_dadt[exist_jpl[ind]] # adding that value from the object (this object have both values)
        dadt_sigma_jpl[i]=jpl_dadt_sigma[exist_jpl[ind]] # adding that value from the object (this object have both values)
        
        
dadt_jpl[len(marco_number):]=jpl_dadt[new] # after Marco's objects are JPL objects
dadt_sigma_jpl[len(marco_number):]=jpl_dadt_sigma[new] # after Marco's objects are JPL objects

dadt_papers[:len(marco_number)]=marco_dadt # Marco's objects
dadt_sigma_papers[:len(marco_number)]=marco_dadt_sigma # Marco's objects

# selectingonly those with SNR>3
selection_papers=np.abs(dadt_papers)/dadt_sigma_papers>3
selection_jpl=np.abs(dadt_jpl)/dadt_sigma_jpl>3

# ******************************************************************************
# Checking periods on LCDB

with open("lc_summary_pub.txt", "r") as f:
    data = f.readlines() # reading whole lines as strings

data=data[5:] # skipping header

# Initialiying array of rotation periods
lcdb_P=np.array([np.nan for _ in range(len(data))])

# making arrays for comparison with JPL and Marco
lcdb_number=[]
lcdb_name=[]
lcdb_min=[]
lcdb_max=[]


for i in range(len(data)):
    period1=(data[i][145:158]).replace(' ','') # rotation period
    lcdb_min.append((data[i][177:181]).replace(' ','')) # minimum magnitude
    lcdb_max.append((data[i][182:186]).replace(' ','')) # maximum magnitude
    number1=(data[i][0:7]).replace(' ','')
    name1=(data[i][10:40]).replace(' ','')
    lcdb_number.append(number1)
    lcdb_name.append(name1)
    if len(period1)>0: # "0" is flag that there is no period in LCDB
        lcdb_P[i]=period1
        

# rouding to 3 decimals in order to compare values
lcdb_P=np.round(lcdb_P,3)
P=np.round(P,3)

# initializing final array of periods from LCDB. Here are written only those which have periods in LCDB which are different from JPL or Marco
P_lcdb=np.array([np.nan for _ in range(len(number))]) # final arraz of LCDB periods which is written to our database
min_mag=np.array([np.nan for _ in range(len(number))]) # minimum reported magnitude
max_mag=np.array([np.nan for _ in range(len(number))]) # maximum reported magnitude

spin_axes_flag=np.array([' ' for _ in range(len(number))]) # "Y" means that object has determined axes pole

br=0
for i in range(len(name)):
    if name[i] in lcdb_name: # checking name
        spin_axes_flag[i]=(data[lcdb_name.index(name[i])][200:203]).replace(' ', '')
        try:
            min_mag[i]=lcdb_min[lcdb_name.index(name[i])]
            max_mag[i]=lcdb_max[lcdb_name.index(name[i])]
        except:
            pass
        
        if lcdb_P[lcdb_name.index(name[i])]!=P[i]: # onlz if the period different from othat we alreadz have from JPL or Marco
            P_lcdb[i]=lcdb_P[lcdb_name.index(name[i])]
            
    
    elif number[i] in lcdb_number: # checking number
        spin_axes_flag[i]=(data[lcdb_number.index(number[i])][200:203]).replace(' ', '')
        if lcdb_P[lcdb_number.index(number[i])]!=P[i]:
            P_lcdb[i]=lcdb_P[lcdb_number.index(number[i])]
            
# ******************************************************************************
# spin axes in LCDB       
spin_axes_lcdb=[] 
with open("lc_spinaxis_pub.txt", "r") as f1:
    data1 = f1.readlines() # reading whole line as a string

# this database has specific structure so it is neccessary to identify the line at
# which the relevant data are written (we take only the most recent estimation of pole axes direction)
spin_number=[]
spin_name=[]
spin_start=[] # line at which start data from specific object (oldest data)
spin_stop=[] # line at which ends data from specific object (mos recent data)


for i in range(len(data1)):
    if data1[i]=='\n': # serching for empty line
        try:
            spin_number.append(int(data1[i+1][:7]))
            spin_name.append((data1[i+1][10:40]).replace(' ', ''))
            spin_start.append(i)
            spin_stop.append(i)
        except:
            pass

spin_start=np.array(spin_start[:-1])+2
spin_stop=np.array(spin_stop[1:])

start=[]
stop=[]
last=[]

spin_number1=[] #  number in our database
spin_name1=[] # name in our database

# 4 solutions for axes direction
l1=np.array([np.nan for _ in range(len(number))]) 
b1=np.array([np.nan for _ in range(len(number))])
l2=np.array([np.nan for _ in range(len(number))])
b2=np.array([np.nan for _ in range(len(number))])
l3=np.array([np.nan for _ in range(len(number))])
b3=np.array([np.nan for _ in range(len(number))])
l4=np.array([np.nan for _ in range(len(number))])
b4=np.array([np.nan for _ in range(len(number))])


for i in range(len(name)):
    if spin_axes_flag[i]=='Y':
        if name[i] in spin_name:

            last=spin_stop[spin_name.index(name[i])]-1
            spin_number1.append(number[i])
            spin_name1.append(name[i])
 
            try:
                x=np.float(data1[last][64:70])
                try:
                    if x!=0:
                        l1[i]=x
                except:
                    pass
            except:
                pass
            
            try:
                x=np.float(data1[last][70:76])
                try:
                    if x!=0:
                        b1[i]=x
                except:
                    pass
            except:
                pass
            
            try:
                x=np.float(data1[last][76:82])
                try:
                    if x!=0:
                        l2[i]=x
                except:
                    pass
            except:
                pass
            
            try:
                x=np.float(data1[last][82:88])
                try:
                    if x!=0:
                        b2[i]=x
                except:
                    pass
            except:
                pass
            
            try:
                x=np.float(data1[last][88:94])
                try:
                    if x!=0:
                        l3[i]=x
                except:
                    pass
            except:
                pass
            
            try:
                x=np.float(data1[last][94:100])
                try:
                    if x!=0:
                        b3[i]=x
                except:
                    pass
            except:
                pass
            
            try:
                x=np.float(data1[last][100:106])
                try:
                    if x!=0:
                        l4[i]=x
                except:
                    pass
            except:
                pass
            
            try:
                x=np.float(data1[last][106:112])
                try:
                    if x!=0:
                        b4[i]=x
                except:
                    pass
            except:
                pass
            
# ******************************************************************************
# colors in SDSS
sdss_name=[] # for comparison with our database
sdss_number=[]  # for comparison with our database  

u=np.array([np.nan for _ in range(len(number))])
u_sigma=np.array([np.nan for _ in range(len(number))])
g=np.array([np.nan for _ in range(len(number))])
g_sigma=np.array([np.nan for _ in range(len(number))])
r=np.array([np.nan for _ in range(len(number))])
r_sigma=np.array([np.nan for _ in range(len(number))])
ii=np.array([np.nan for _ in range(len(number))])
ii_sigma=np.array([np.nan for _ in range(len(number))])
z=np.array([np.nan for _ in range(len(number))])
z_sigma=np.array([np.nan for _ in range(len(number))])
       
with open("ADR3.dat", "r") as f1:
    data2 = f1.readlines()      
            
for i in range(len(data2)):
    sdss_number.append(int(data2[i][244:250]))
    ime=(data2[i][250:271]).replace('_','')
    ime=ime.replace(' ', '')
    sdss_name.append(ime)

br=0
for i in range(len(number)):
    ind=-99 # our object is not found in SDSS
    try:
        if name[i] in sdss_name:
            br+=1
            ind=sdss_name.index(name[i])
        elif int(number[i]) in sdss_number:
            br+=1
            ind=sdss_number.index(int(number([i])))
   
    except:
        pass
    
    if ind!=-99:
            
        u[i]=data2[ind][163:169]
        u_sigma[i]=data2[ind][169:174]
        g[i]=data2[ind][174:180]
        g_sigma[i]=data2[ind][180:185]
        r[i]=data2[ind][185:191]
        r_sigma[i]=data2[ind][191:196]
        ii[i]=data2[ind][196:202]
        ii_sigma[i]=data2[ind][202:207]
        z[i]=data2[ind][207:213]
        z_sigma[i]=data2[ind][213:218]

# Ordering according to SNR

snr_jpl=abs(dadt_jpl/dadt_sigma_jpl)
snr_papers=abs(dadt_papers/dadt_sigma_papers)

snr_max=np.zeros(len(snr_papers)) # array of maximum values of SNR (some object have more than one value (papers and JPL) so we take the better one)

# for object which have data both form the papers and JPL we take better value of SNR
for i in range(len(number)):
    if snr_jpl[i]>snr_papers[i] or np.isnan(snr_papers[i]):
        snr_max[i]=snr_jpl[i]
    else:
        snr_max[i]=snr_papers[i]

order = snr_max.argsort()[::-1] # sorting according to SNR

# making CSV database
# open the file in the write mode
f = open('database.csv', 'w')

# create the csv writer
writer = csv.writer(f)
# write a row to the csv file
first_row=['Number', 'Name', 'a', 'e', 'i', 
           'H', 'H_std', 'D', 'D_std', 'P', 'P_std', 'P (LCDB)', 
           'l1', 'b1', 'l2', 'b2','l3', 'b3','l4', 'b4', 'Mag(min)', 'Mag(max)',
           'u', 'u_sigma', 'g', 'g_sigma','r', 'r_sigma','i', 'i_sigma','z', 'z_sigma',
           'dadt (papers)', 'dadt_std (papers)', 'snr_papers',
           'dadt (JPL)', 'dadt_sigma (JPL)', 'snr_jpl']

writer.writerow(first_row)

for i in range(len(number)):
    
    writer.writerow([number[order][i], name[order][i], np.round(a[order][i],9), np.round(e[order][i],9), np.round(inc[order][i],9), 
                     H[order][i], H_sigma[order][i], D[order][i], D_sigma[order][i], P[order][i], P_sigma[order][i], P_lcdb[order][i], 
                     l1[order][i], b1[order][i], l2[order][i], b2[order][i],l3[order][i], b3[order][i],l4[order][i], b4[order][i], min_mag[order][i], max_mag[order][i],
                     u[order][i], u_sigma[order][i], g[order][i], g_sigma[order][i], r[order][i], r_sigma[order][i], ii[order][i], ii_sigma[order][i], z[order][i], z_sigma[order][i],
                     dadt_papers[order][i], dadt_sigma_papers[order][i], np.round(snr_papers[order][i],2),
                     np.round(dadt_jpl[order][i],6), np.round(dadt_sigma_jpl[order][i],6), np.round(snr_jpl[order][i],2)])

# close the file
f.close()
#

            

    


