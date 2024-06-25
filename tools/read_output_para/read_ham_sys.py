import numpy as np

# read the dimension of Hs
nrho = np.loadtxt("ham_sys.data", max_rows=1, dtype = int)
# read Hs reshaped in the nrho * nrho form
hs = np.loadtxt("ham_sys.data", skiprows=1, dtype = float).reshape(nrho,nrho)

