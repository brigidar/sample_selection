#!/usr/bin/env python
#########################################################################################
#											#
# Name	      :	nb_samples.py								#
# Version     : 0.2								#
# Project     : targeted Metabolomics						#
# Description : Script to select stool samples based on weight, sample id, time or combination of them		#
# Author      : Brigida Rusconi								#
# Date        : August 15th, 2016							#
#											#
#########################################################################################


import argparse, os, sys, csv,pandas,pdb, numpy
from pandas import *
from numpy import *
import math
from math import floor

parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="table of nb of samples above cutoff")
parser.add_argument('-i', '--table', help="stool master list as txt file")
parser.add_argument('-w', '--weight', help="min sample weight",const=2400,nargs='?')
parser.add_argument('-b', '--DOB', help="list with DOB. Required: header: Infant date of birth, file format: txt")
parser.add_argument('-s', '--sample', help="list with samples. Required: header: Sample, file format: txt ")
parser.add_argument('-d', '--days', help="list with days.Required: header: Days, file format: txt")
parser.add_argument('-t', '--DOL', help="make table with DOL of stool sample list (True/False)",default=False)
parser.add_argument('-m', '--matched', help="list with sample and days to look up")

args = parser.parse_args()
output_file = args.output
input_file = args.table
DOB = args.DOB
dol_table=args.DOL
#pdb.set_trace()
#-------------------read tables and parse them---------------------------------



test=read_csv(input_file,sep='\t', dtype=object)
dob=read_csv(DOB, dtype=object)

test.rename(columns={'Time Collected':'time_collected'},inplace=True)
dob.rename(columns={'Infant date of birth':'DOB'},inplace=True)
dob.DOB=to_datetime(dob.DOB)
test.dropna(subset=['time_collected'], inplace=True)
test[test.time_collected != 'D/T Unknown']
test.time_collected=to_datetime(test.time_collected, errors='coerce')
test.Weight=to_numeric(test.Weight,errors='coerce')
test.dropna(subset=['Weight'], inplace=True)

#sum weight of duplicate values
df1=test.groupby(['Study','time_collected'])['Weight'].sum().reset_index()
df1.set_index(['Study','time_collected'], inplace=True)

#reorder dataframe without duplicates of study and time collected
test.drop_duplicates(subset=['Study','time_collected'], inplace=True)
test.set_index(['Study','time_collected'],inplace=True)

test.drop(['Weight'],axis=1,inplace=True)
#replace with sum values of weight

df3=concat([test,df1['Weight']],axis=1,join_axes=[df1.index])

df3.reset_index(inplace=True)

#----------------optional------------- weight of sample--------------------------------
if args.weight:
    w=args.weight
    df5=DataFrame()
    
#get samples with enough weight
    for i in df3.index:
        if df3.Patient[i]>='100.01':
            df4=df3[df3.Weight >=w]
        else:
            df5=df3[df3.Weight >=250]

    if df5.index.size>0:
        df6=concat([df4,df5],axis=0)
    else:
        df6=df4
else:
    df6=df3


#----------------optional-------------selected samples--------------------------------

if args.sample:
    samp=args.sample
    samp1=read_csv(samp,sep='\t',dtype=object)
    samp1.sort('Sample',inplace=True)
    samp1.set_index('Sample',inplace=True)
    df6.set_index('Patient',inplace=True)
    df6[df6.index.isin(samp1.index)]
    df6.reset_index(inplace=True)
    #group samples by patient
    nb_sampl=df6.groupby('Patient').size()
else:
    nb_sampl=DataFrame(df6.groupby('Patient').size(),columns=['count'])

#--------------------count of samples-----------------------------
with open(output_file ,'w') as output:
    nb_sampl.to_csv(output, sep='\t')

#----------------------------------------------------------------------------------------
#-----------------DOL information if wanted--------------------------------------------------------
if dol_table=='True':
    
    #make list with DOL value
    df6.set_index(['Patient'],inplace=True)
    dob.set_index(['Study ID'],inplace=True)

    df7=concat([df6,dob['DOB']],axis=1,join_axes=[df6.index])
    df7['DOL']=df7['time_collected']-df7['DOB']

    #convert to float days
    new=list()
    for item in df7['DOL']:
        new.append(float("{0:.2f}".format(item.total_seconds()/86400)))

    df7.insert(df7.columns.size, "DOL_dec",new)
#    df7.rename(columns={'Unnamed: 11':'notes'},inplace=True)
#    df7.drop(['notes'],axis=1,inplace=True)
    df7.reset_index()

    #--------------------optional---------select time point-------------------------------------------
    if args.days and (not args.sample) and (not args.matched):
        day=args.days
        dy1=read_csv(day,sep='\t',dtype=object)
        dy1=list(to_numeric(dy1.Days))
        rd=list()
        for item in df7['DOL_dec']:
            rd.append(floor(item))
        df7.insert(df7.columns.size,'rounded',rd)
        df7.reset_index(inplace=True)
        pos=list()
        for i in df7.index:
            if float(df7['rounded'][i]) in dy1:
                pos.append(i)
        df8=df7[df7.index.isin(pos)]
        with open("stool_timepoints.txt" ,'w') as output:
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
        for item in df7['DOL_dec']:
            rd.append(floor(item))
        df7.insert(df7.columns.size,'rounded',rd)
        df7.reset_index(inplace=True)
        pos=list()
        for i in df7.index:
            if float(df7['rounded'][i]) in dy1:
                pos.append(i)
        df8=df7[df7.index.isin(pos)]
        df8.set_index('Patient',inplace=True)
        df8=df8[df8.index.isin(samp1.index)]
        df8.reset_index(inplace=True)
        with open("stool_timepoints.txt" ,'w') as output:
            df8.to_csv(output, sep='\t')

    elif args.matched and (not args.days) and (not args.sample):
        match=args.matched
        matched1=read_csv(match,sep='\t',dtype=object)
        matched1.Days=to_numeric(matched1.Days)
        #make a range from the days that are selected to have more control hits, maybe modify so that interval can be changed and not fixed
        re=list()
        id=list()
        for i,item in enumerate(matched1.Days):
            re.append(range((int(item)-3),(int(item)+3)))
            id.append(repeat(matched1.Sample[i],len(re)).tolist())

        flat_re=[n for item in re for n in item]
        flat_id=[n for item in id for n in item]
        ri=zip(flat_id,flat_re)
        matched2=DataFrame(ri, columns=['Sample','Days'])
        matched2.Days=matched2.Days.astype(float)
        matched2.set_index(['Sample','Days'], inplace=True)
        rd=list()
        for item in df7['DOL_dec']:
            rd.append(floor(item))
        df7.insert(df7.columns.size,'rounded',rd)
        df7.reset_index(inplace=True)
        df7.set_index(['Patient','rounded'],inplace=True)
        pdb.set_trace()
        df8=df7[df7.index.isin(matched2.index)]
#pdb.set_trace()


        with open("stool_timepoints.txt" ,'w') as output:
            df8.to_csv(output, sep='\t')

