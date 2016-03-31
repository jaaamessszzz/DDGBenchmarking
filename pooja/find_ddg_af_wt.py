import os
import sys
import re
import numpy

#jobs_1
mut=os.getcwd()[38:43]
pdb=os.getcwd()[33:37]

#jobs_2
#pdb=sys.argv[1]
#mut=sys.argv[1]
#pub=sys.argv[3]

mut_complexes=[]
wt_complexes=[]
mut_monomers=[]
wt_monomers=[]
mut_partners=[]
wt_partners=[]

for filename in os.listdir(os.getcwd()):
    if filename.endswith("last_mutate_0001_0001_natro_0001.pdb"):
        mut_complexes.append(filename)
    if filename.endswith("last_repack_0001_0001_natro_0001.pdb"):
        wt_complexes.append(filename)
    if filename.endswith("last_mutate_0001_0001_monomer_natro_0001.pdb"):
        mut_monomers.append(filename)
    if filename.endswith("last_repack_0001_0001_monomer_natro_0001.pdb"):
        wt_monomers.append(filename)
    if filename.endswith("last_mutate_0001_0001_partner_natro_0001.pdb"):
        mut_partners.append(filename)
    if filename.endswith("last_repack_0001_0001_partner_natro_0001.pdb"):
        wt_partners.append(filename)

mut_complexes.sort()
wt_complexes.sort()
mut_monomers.sort()
wt_monomers.sort()
mut_partners.sort()
wt_partners.sort()

mut_energy=[]
wt_energy=[]
mut_mon_energy=[]
wt_mon_energy=[]
mut_part_energy=[]
wt_part_energy=[]

for mut_complex in mut_complexes:
    open_mut=open(mut_complex, "r")
    for line in open_mut:
        if "pose" in line:
            mut_energy.append(float(line.split()[17]))
    open_mut.close()

for wt_complex in wt_complexes:
    open_wt=open(wt_complex, "r")
    for line in open_wt:
        if "pose" in line:
            wt_energy.append(float(line.split()[17]))
    open_wt.close()

for mut_monomer in mut_monomers:
    open_mut_mon=open(mut_monomer, "r")
    for line in open_mut_mon:
        if "pose" in line:
            mut_mon_energy.append(float(line.split()[17]))
    open_mut_mon.close()

for wt_monomer in wt_monomers:
    open_wt_mon=open(wt_monomer, "r")
    for line in open_wt_mon:
        if "pose" in line:
            wt_mon_energy.append(float(line.split()[17]))
    open_wt_mon.close()

#print len(wt_mon_energy), len(mut_mon_energy)

for mut_partner in mut_partners:
    open_mut_part=open(mut_partner, "r")
    for line in open_mut_part:
        if "pose" in line:
            mut_part_energy.append(float(line.split()[17]))
    open_mut_part.close()

for wt_partner in wt_partners:
    open_wt_part=open(wt_partner, "r")
    for line in open_wt_part:
        if "pose" in line:
            wt_part_energy.append(float(line.split()[17]))
    open_wt_part.close()

print len(mut_energy), len(mut_part_energy), len(mut_mon_energy), len(wt_energy), len(wt_part_energy), len(wt_mon_energy)


dg=[]
for i in range(0,len(mut_energy)):
    dg.append(mut_energy[i]-mut_part_energy[i]-mut_mon_energy[i]-wt_energy[i]+wt_part_energy[i]+wt_mon_energy[i])

print numpy.mean(dg), numpy.std(dg)
ddg_boltzmann_averaged=numpy.mean(dg) #kT=infinity


'''
boltzmann_factor=[]

for i in range(len(wtmut_energy_complex)):
    boltzmann_factor.append(numpy.exp(-(wtmut_energy_complex[i]-wt_energy_complex[i])/0.5))

boltzmann_sum= sum(boltzmann_factor)

ddg_boltzmann_weighted=[]

for i in range(len(wtmut_energy_complex)):
    ddg_boltzmann_weighted.append(float((wtmut_energy_complex[i]-wt_energy_complex[i])*boltzmann_factor[i]))

ddg_boltzmann_averaged=sum(ddg_boltzmann_weighted)/boltzmann_sum
'''

print ddg_boltzmann_averaged


ddg_wt=open("/netapp/home/psuresh2/DDG/jobs_1/ddg_val_wt_ktlow", "a")
ddg_wt.write(str(pdb)+"\t"+str(mut)+"\t"+str(ddg_boltzmann_averaged)+"\n")
ddg_wt.close()
