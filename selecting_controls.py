#!/usr/bin/env python

#########################################################################################
#											#
# Name	      :	selecting_controls.py								#
# Version     : 0.2								#
# Project     : targeted Metabolomics						#
# Description : Script to select matched controls based on gestational age, brithweight, DOB		#
# Author      : Brigida Rusconi								#
# Date        : August 17th, 2016							#
#											#
#########################################################################################


import argparse, os, sys, csv,pandas,pdb, numpy
from pandas import *
from numpy import *

parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="matched_controls")
parser.add_argument('-i', '--table', help="cases and controls")
parser.add_argument('-b', '--birthweight', help="difference in birthweight",default=100)
parser.add_argument('-g', '--gestational', help="difference in gestational age in weeks", default=1)
parser.add_argument('-e', '--exclusion', help="exclusion set")

args = parser.parse_args()
output_file = args.output
input_file = args.table
w=args.birthweight
g=args.gestational
e=args.exclusion

#reading in table and formatting
df =read_csv(input_file,sep='\t', dtype=object)
exc=read_csv(e,sep='\t', dtype=object)
df.Gestational=to_numeric(df.Gestational)
df.rename(columns={'Necrotizing enterocolitis':'NEC'},inplace=True)
df.rename(columns={'Birthweight':'BW'},inplace=True)
df.rename(columns={'Gestational age:  days':'GA_days'},inplace=True)
df.rename(columns={'Gestational age:  weeks':'GA_weeks'},inplace=True)
df.rename(columns={'Study ID':'Sample'},inplace=True)
df.rename(columns={'Infant date of birth':'DOB'},inplace=True)
df.BW=to_numeric(df.BW)
df.DOB=to_datetime(df.DOB)

#df1=df[df.Gestational < 27] #depends if we want to only look at less than 27 weeks cases
#removes samples that have sepsis and other exclusion criteria
uniq=np.unique(exc[['SIP','CHD','no_samples']])
df2=df[~df['Sample'].isin(uniq)]
cases1=df2[df2.NEC=="Yes"]
cases=cases1[(cases1.Sample).astype(float) >= float(33.01)]
controls1=df2[df2.NEC!="Yes"]
cases.sort(columns=['Sample'],inplace=True)

#sepsis is only excluded from controls
controls=controls1[~controls1['Sample'].isin(exc['Sepsis'])]
controls.reset_index(inplace=True)
cases.reset_index(inplace=True)

#find all controls that match gestational age +/- 1 week
li=list()
id=list()
for i in cases.index:
    li.append(list(arange((cases.Gestational[i]-g),(cases.Gestational[i]+g+0.01),0.1)))
    id.append(repeat(cases.Sample[i],len(li[i])).tolist())
flat_re=[n for item in li for n in item]
flat_id=[n for item in id for n in item]
ri=zip(flat_id,flat_re)
tes1=concatenate([z for z in ri])
tes1=tes1.reshape(-1,2)
gesta=DataFrame(tes1,columns=['Sample','Gestational'])
gesta.Gestational=gesta.Gestational.astype(str)


dat=list()
dat1=list()
for x,item in enumerate(controls.Gestational):
    if str(item) in gesta.Gestational.values:
        dat.append(gesta[gesta.Gestational==str(item)]['Sample'].values.tolist())
        dat1.append(repeat(controls.Sample[x],len(dat[-1])).tolist())

flat_re1=[n for item in dat for n in item]
flat_id1=[n for item in dat1 for n in item]
dat2=list(set(zip(flat_re1,flat_id1)))


#find all controls that match brithweight +/- 100g

wei=list()
id=list()
for i in cases.index:
    w1=list()
    for x in controls.index:
        if controls.BW[x] in range((cases.BW[i]-w),(cases.BW[i]+(w+1)),1):
            w1.append(controls.Sample[x])
    wei.append(w1)
    id.append(repeat(cases.Sample[i],len(w1)).tolist())

flat_w=[n for item in wei for n in item]
flat_id2=[n for item in id for n in item]
dat3=list(set(zip(flat_id2,flat_w)))

#find overlapping controls
m=list()
for item in dat2:
    if item in dat3:
       m.append(item)

tes=concatenate([z for z in m])
tes=tes.reshape(-1,2)
matched=DataFrame(tes,columns=['Cases','Controls'])
matched.sort(columns=['Cases'],inplace=True)
m2=list()
for item in matched.Cases.unique():
    m2.append('/'.join([str(n) for n in matched[matched.Cases==item]['Controls']]))

cases.insert(cases.columns.size,"matched",m2)

m=list()
for item in m2:
    m.append(item.split('/'))

#find distance to DOB of cases
data1=list()
for i,v in enumerate(cases.index):
    dat=list()
    for item in m[i]:
        tim=controls.DOB[controls.Sample==str(item)]-cases.DOB[v]
        test=tim.tolist()
        t=int(str(str(str(test[0]).split("'")).split(" days")[0])[2:])
        dat.append(t)
    data1.append(dat)

#select closest by DOB
closest=list()
for i,item in enumerate(data1):
    if len(item)>1:
        gi=item.index(min(item, key=abs))
        closest.append(str(m[i][gi]))
    elif len(item)==1:
        closest.append('unique')
    else:
        closest.append('no_match')


#drop the closest from the date and redo minimum
def remov(list_match,to_remove):
    sec=list()
    for i,item in enumerate(list_match):
        if to_remove[i]=='no_match' or to_remove[i]=='unique':
            sec.append(to_remove[i])
        else:
            for s,x in enumerate(item):
                if str(x)==to_remove[i]:
                    sec.append(s)
    return sec
sec=remov(m,closest)

for i,item in enumerate(closest):
    if item=='unique':
        closest[i]=str(m[i][0])
#pdb.set_trace()
for i,item in enumerate(data1):
    if len(item)>1:
        item.pop(sec[i])
        m[i].pop(sec[i])
    else:
        m[i]='unique'

sec_closest=list()
for i,item in enumerate(data1):
    if len(item)>1:
        gi=item.index(min(item, key=abs))
        sec_closest.append(str(m[i][gi]))
    elif len(item)==1:
        sec_closest.append('unique')
    else:
        sec_closest.append('no_match')

sec2=remov(m,sec_closest)

for i,item in enumerate(sec_closest):
    if item=='unique':
        sec_closest[i]=str(m[i][0])
for i,item in enumerate(data1):
    if len(item)>1:
        item.pop(sec2[i])
        m[i].pop(sec2[i])
    else:
        m[i]='unique'


tir_closest=list()
for i,item in enumerate(data1):
    if len(item)>1:
        gi=item.index(min(item, key=abs))
        tir_closest.append(str(m[i][gi]))
    elif len(item)==1:
        tir_closest.append('unique')
    else:
        tir_closest.append('no_match')
sec3=remov(m,tir_closest)

for i,item in enumerate(tir_closest):
    if item=='unique':
        tir_closest[i]=str(m[i][0])
for i,item in enumerate(data1):
    #del item[sec[i]]
    if len(item)>1:
        item.pop(sec3[i])
        m[i].pop(sec3[i])
    else:
        m[i]='unique'

four_closest=list()
for i,item in enumerate(data1):
    if len(item)>0:
        gi=item.index(min(item, key=abs))
        four_closest.append(str(m[i][gi]))
    elif len(item)==1:
        four_closest.append('unique')
    else:
        four_closest.append('no_match')

for i,item in enumerate(four_closest):
    if item=='unique':
        four_closest[i]=str(m[i][0])
#identify closest control that are duplicated

def dup(thelist):
    d=list()
    seen = set()
    for x in thelist:
        if x in seen:
            d.append('yes')
        else:
            d.append('no')
            seen.add(x)
    return d

#check if second best is duplicated with first or within secondary hits

def dup2(thelist,thelist2):
    d=list()
    seen = set(thelist2)
    for x in thelist:
        if x in seen:
            d.append('yes')
        else:
            d.append('no')
            seen.add(x)
    return d




#duplicated within the column
dupli=dup(closest)
duplit2=dup(sec_closest)


#values of first present in second
dup_c=dup2(sec_closest,closest)
#all together
merg=list()
for i,item in enumerate(closest):
    mer=list()
    mer.extend([closest[i],sec_closest[i],tir_closest[i],four_closest[i]])
    merg.append(mer)

d2=list()
seen = set()
for item in merg:
    d=list()
    for x in item:
        if x in seen:
            d.append('yes')
        else:
            d.append('no')
            seen.add(x)
    d2.append(list(d))
dupl=list()
for item in d2:
    b=set(item)
    if len(b)==1:
        dupl.extend(list(b))
    else:
        dupl.append('yes')

#pdb.set_trace()
#combine the two
def combin(list1,list2):
    comb=list()
    for i,item in enumerate(list1):
        if item=='yes':
            comb.append(item)
        elif dupli[i]=='yes':
            comb.append(list2[i])
        else:
                comb.append(item)
    return comb
#duplicated between the first two
comb_2=combin(dupli,dup_c)



#append columns control

cases.insert(cases.columns.size, "closest",closest)
cases.insert(cases.columns.size, "2nd_closest",sec_closest)
cases.insert(cases.columns.size, "3rd_closest",tir_closest)
cases.insert(cases.columns.size, "4th_closest",four_closest)
cases.insert(cases.columns.size, "duplicated_first_two",comb_2)
cases.insert(cases.columns.size, "duplicated_all",dupl)
cases=cases.drop("fract",1)
cases=cases.drop("GA_weeks",1)
cases=cases.drop("GA_days",1)
#pdb.set_trace()

#write table
with open(output_file ,'w') as output:
    cases.to_csv(output, sep='\t', index=False)

#all possible matches
with open("all_matches.txt",'w') as output:
    unique_matches=set()
    for item in m2:
        if '/' in item:
            for n in item.split('/'):
                unique_matches.add(n)
        else:
            unique_matches.add(item)
    uniq_m=DataFrame({'matches':list(unique_matches)})
    uniq_m=uniq_m[uniq_m['matches']!='no_match']
    uniq_m.to_csv(output, sep='\t',index=False)


              #        g=list(gesta[gesta.Sample==item]['Gestational'])
              #        cont_sel=df7.append(l[l['rounded'].isin(g)])
              
              
              #l=list()
              #for i in cases.index:
              #    l1=list()
              #
              #    for x in controls.index:
              #        if controls.Gestational[x] in linspace((cases.Gestational[i]-g),(cases.Gestational[i]+g),(g*10+1)):
              #            l1.append(controls.Sample[x])
              #    l.append(l1)

