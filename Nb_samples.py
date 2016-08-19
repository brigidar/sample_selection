#!/usr/bin/env python
#########################################################################################
#											#
# Name	      :	nb_samples.py								#
# Version     : 0.4								#
# Project     : targeted Metabolomics						#
# Description : Script to select stool samples based on weight, sample id, time or combination of them		#
# Author      : Brigida Rusconi								#
# Date        : August 17th, 2016							#
#											#
#########################################################################################


import argparse, os, sys, csv,pandas,pdb, numpy
from pandas import *
from numpy import *
import math
from math import floor
import re
from re import search

parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="table of nb of samples above cutoff")
parser.add_argument('-i', '--table', help="stool master list as txt file")
parser.add_argument('-w', '--weight', help="min sample weight",const=2400,nargs='?')
parser.add_argument('-b', '--DOB', help="list with DOB. Required: header: Infant date of birth, file format: txt")
parser.add_argument('-s', '--sample', help="list with samples. Required: header: Sample, file format: txt ")
parser.add_argument('-d', '--days', help="list with days.Required: header: Days, file format: txt")
parser.add_argument('-t', '--DOL', help="make table with DOL of stool sample list (True/False)",default=False)
parser.add_argument('-m', '--matched', help="list with sample and days to look up")
parser.add_argument('-a','--all',help="output with all matching time points and samples")

args = parser.parse_args()
output_file = args.output
input_file = args.table
DOB = args.DOB
dol_table=args.DOL
#pdb.set_trace()
#-------------------read tables and parse them---------------------------------



data=read_csv(input_file,sep='\t', dtype=object)
dob=read_csv(DOB, dtype=object)
data.rename(columns={'Time Collected':'time_collected'},inplace=True)
dob.rename(columns={'Infant date of birth':'DOB'},inplace=True)
dob.DOB=to_datetime(dob.DOB)
data.dropna(subset=['time_collected'], inplace=True)
data[data.time_collected != 'D/T Unknown']
data.time_collected=to_datetime(data.time_collected, errors='coerce')
data.Weight=to_numeric(data.Weight,errors='coerce')
data.dropna(subset=['Weight'], inplace=True)

#sum weight of duplicate values
#df1=data.groupby(['Study','time_collected'])['Weight'].sum().reset_index()
#http://pandas.pydata.org/pandas-docs/stable/groupby.html

data['total_weight']=data.groupby(['Study','time_collected'])['Weight'].transform('sum')

#pdb.set_trace()
#reorder dataframe without duplicates of study and time collected

#replace with sum values of weight



#----------------optional------------- weight of sample--------------------------------
if args.weight:
    w=args.weight
    df5=DataFrame()
    
#get samples with enough weight
# & has higher precedence than == BE CAREFUL http://pandas.pydata.org/pandas-docs/stable/indexing.html#boolean-indexing
# with unique each item is only iterated once
#have to remove any characters so that it will not fail on the float conversion

    uniq=list(data.Patient.unique())
    for i,item in enumerate(uniq):
        if bool(search(r'\d',item)) is not True:
            uniq.pop(i)
    for item in uniq:
            if float(item)>=float(100.01):
                df5=df5.append(data[(data.Patient==item) & (data.total_weight >=float(w))])
            else:
                df5=df5.append(data[(data.Patient==item) & (data.total_weight >=float(250))])

else:
    df5=data[data.Weight>0]


#pdb.set_trace()
#----------------optional-------------selected samples--------------------------------

if args.sample:
    samp=args.sample
    samp1=read_csv(samp,sep='\t',dtype=object)
    samp1.sort('Sample',inplace=True)
    samp1.set_index('Sample',inplace=True)
    df5.set_index('Patient',inplace=True)
    df5=df5[df5.index.isin(samp1.index)]
    df5.reset_index(inplace=True)
    #group samples by patient
    nb_sampl=df5.groupby('Patient').size()
else:
    nb_sampl=DataFrame(df5.groupby('Patient').size(),columns=['count'])

#--------------------count of samples-----------------------------
with open(output_file ,'w') as output:
    nb_sampl.to_csv(output, sep='\t')

#----------------------------------------------------------------------------------------
#-----------------DOL information if wanted--------------------------------------------------------
if dol_table=='True':
    output2=args.all
    
    #make list with DOL value
    df5.set_index(['Patient'],inplace=True)
    dob.set_index(['Study ID'],inplace=True)

    df6=concat([df5,dob['DOB']],axis=1,join_axes=[df5.index])
    df6['DOL']=df6['time_collected']-df6['DOB']

    #convert to float days
    new=list()
    for item in df6['DOL']:
        new.append(float("{0:.2f}".format(item.total_seconds()/86400)))

    df6.insert(df6.columns.size, "DOL_dec",new)


    #--------------------optional---------select time point-------------------------------------------
    if args.days and (not args.sample) and (not args.matched):
        day=args.days
        dy1=read_csv(day,sep='\t',dtype=object)
        dy1=list(to_numeric(dy1.Days))
        rd=list()
        for item in df6['DOL_dec']:
            rd.append(floor(item))
        df6.insert(df6.columns.size,'rounded',rd)
        df6.reset_index(inplace=True)
        pos=list()
        for i in df6.index:
            if float(df6['rounded'][i]) in dy1:
                pos.append(i)
        df8=df6[df6.index.isin(pos)]
        with open(output2 ,'w') as output:
            df8.to_csv(output, sep='\t')

    #--------------------optional---------select time point for select samples-------------------------------------------
    elif args.days and args.sample and (not args.matched):
        day=args.days
        dy1=read_csv(day,sep='\t',dtype=object)
        dy1=list(to_numeric(dy1.Days))
        samp=args.sample
        samp1=read_csv(samp,sep='\t',dtype=object)
        samp1.sort('Sample',inplace=True)
        samp1.set_index('Sample',inplace=True)
        rd=list()
        for item in df6['DOL_dec']:
            rd.append(floor(item))
        df6.insert(df6.columns.size,'rounded',rd)
        df6.reset_index(inplace=True)
        pos=list()
        for i in df6.index:
            if float(df6['rounded'][i]) in dy1:
                pos.append(i)
        df8=df6[df6.index.isin(pos)]
        df8.set_index('Patient',inplace=True)
        df8=df8[df8.index.isin(samp1.index)]
        df8.reset_index(inplace=True)
        with open(output2 ,'w') as output:
            df8.to_csv(output, sep='\t')
#--------------------optional---------select matched controls/cases-------------------------------------------
    elif args.matched and (not args.days) and (not args.sample):
        match=args.matched
        matched1=read_csv(match,sep='\t',dtype=object)
        matched1.Days=to_numeric(matched1.Days)
        
        cases=matched1[matched1.Phenotype=='Cases']
        #find closest day to onset
        cases.reset_index(inplace=True)
        rd=list()
        for item in df6['DOL_dec']:
            rd.append(floor(item))
        df6.insert(df6.columns.size,'rounded',rd)
        df6.reset_index(inplace=True)
        closest=list()

        for i,item in enumerate(cases.Sample.unique()):
            l=list(df6[(df6.Patient==item) & (df6.rounded<cases.Days[i])]['rounded'].values)
            if len(l)>0:
                closest.append(min(l, key=lambda x:cases.Days[i]-x))
            else:
                closest.append(cases.Days[i])
        cases.insert(cases.columns.size,'closest',closest)
        cases.set_index(['Sample','closest'],inplace=True)
        df6.set_index(['Patient','rounded'],inplace=True)
        sel_cases=df6[df6.index.isin(cases.index)]
        sel_cases.reset_index(inplace=True)
        #replace closest sample time point for controls to match
        controls=matched1[matched1.Phenotype=='Controls']
        controls.reset_index(inplace=True)

        cases.reset_index(inplace=True)
        for i,item in enumerate(cases.Sample):
            for x,item1 in enumerate(controls.Matching_Case):
                if item==item1:
                    controls.Days[x]=cases.closest[i]

        #make a range from the days that are selected to have more control hits, maybe modify so that interval can be changed and not fixed currently only 3 days before onset
        re=list()
        id=list()
        for i,item in enumerate(controls.Days):
            re.append(range((int(controls.Days[i])-3),(int(controls.Days[i]+1))))
            id.append(repeat(controls.Sample[i],len(re[i])).tolist())
        flat_re=[n for item in re for n in item]
        flat_id=[n for item in id for n in item]
        ri=zip(flat_id,flat_re)
        tes=concatenate([z for z in ri])
        tes=tes.reshape(-1,2)
        matched2=DataFrame(tes,columns=['Sample','Days'])
        matched2.Days=matched2.Days.astype(float)

        df6.reset_index(inplace=True)
        df7=DataFrame()
        for item in matched2.Sample.unique():
            l=df6[df6.Patient==item]
            g=list(matched2[matched2.Sample==item]['Days'])
            df7=df7.append(l[l['rounded'].isin(g)])
        df8=concat([sel_cases,df7],axis=0)
        # append onset day if you use index.value you get the index of the match
        onset=list()
        for item in df8.Patient:
            onset.append(matched1.Days[matched1[matched1.Sample==item].index.values].values[0])
#pdb.set_trace()
        df8.reset_index(inplace=True)
        df8.insert(df8.columns.size,'onset',onset)
                    #pdb.set_trace()
        df8.sort(columns=['NEC','Patient','time_collected'],inplace=True)
        df8.drop('index',inplace=True,axis=1)
        missing=list()
        for item in matched1.Sample:
            if item in df8.Patient.values:
                pass
            else:
                missing.append(item)
#pdb.set_trace()
        with open(output2 ,'w') as output:
            df8.to_csv(output, sep='\t',index=False)

        with open("missing.txt" ,'w') as output:
            for item in missing:
                output.write('%s\n' % item)
    else:
        samp=args.sample
        samp1=read_csv(samp,sep='\t',dtype=object)
        samp1.sort('Sample',inplace=True)
        samp1.set_index('Sample',inplace=True)
        df6=df5[df5.index.isin(samp1.index)]
        df7=concat([df5,dob['DOB']],axis=1,join_axes=[df5.index])
        df7.reset_index(inplace=True)
        df7['DOL']=df7['time_collected']-df7['DOB']
        new=list()
        for item in df7['DOL']:
            new.append(float("{0:.2f}".format(item.total_seconds()/86400)))
    
        df7.insert(df7.columns.size, "DOL_dec",new)
        df8=df7.sort(columns=['NEC','Patient','time_collected'])
#pdb.set_trace()
        with open(output2 ,'w') as output:
            df8.to_csv(output, sep='\t',index=False)