import numpy as np

file_name = "ams_fermion.data"

# read level degrees of system
norbs = np.loadtxt(file_name, max_rows=1, dtype = int)
# read spin degrees of system
nspin = np.loadtxt(file_name, skiprows=1, max_rows=1, dtype = int)

nrho = (nspin*2)**norbs
nrho2 = nrho*nrho

ncount = 0

ams = np.zeros((norbs,nspin,nrho,nrho), dtype = float)

for i in range(norbs):
  for j in range(nspin):
    tmp1 = 2 + ncount * (nrho2 + 2)
    tmp2 = tmp1 + 1
    tmp3 = tmp2 + 1
    # read orbs of ams
    orbs = np.loadtxt(file_name, skiprows=tmp1, max_rows=1, dtype = int)-1
    # read spin of ams
    spin = np.loadtxt(file_name, skiprows=tmp2, max_rows=1, dtype = int)-1
    # read ams reshaped in the nrho * nrho form
    tmpmat = np.loadtxt(file_name, skiprows= tmp3, max_rows=nrho2, dtype = float).reshape(nrho,nrho)
    ncount = ncount + 1
    ams[orbs][spin] = tmpmat

print(ams)
