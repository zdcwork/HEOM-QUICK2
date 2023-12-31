import numpy as np

# read the dimension of relative parameters
nsgn = np.loadtxt("res_corr.data", max_rows=1, dtype = int)
nalf = np.loadtxt("res_corr.data", skiprows=1, max_rows=1, dtype = int)
ncor = np.loadtxt("res_corr.data", skiprows=2, max_rows=1, dtype = int)
nspin = np.loadtxt("res_corr.data", skiprows=3, max_rows=1, dtype = int)
nvar2 = np.loadtxt("res_corr.data", skiprows=4, max_rows=1, dtype = int)


# read cb(nvar2,nspin,ncor,nalf,nsgn)
cb_real = np.loadtxt("res_corr.data", skiprows=5, usecols=6, dtype = np.float64)
cb_img = np.loadtxt("res_corr.data", skiprows=5, usecols=7, dtype = np.float64)
cb = np.zeros((cb_real.shape), dtype=complex)
for i in range(cb_real.size):
  cb[i] = complex(cb_real[i], cb_img[i])


# read cd(nvar2,nspin,ncor,nalf,nsgn)
cd_real = np.loadtxt("res_corr.data", skiprows=5, usecols=8, dtype = np.float64)
cd_img = np.loadtxt("res_corr.data", skiprows=5, usecols=9, dtype = np.float64)
cd = np.zeros((cd_real.shape), dtype=complex)
for i in range(cd_real.size):
  cd[i] = complex(cd_real[i], cd_img[i])


# read cgama(nvar2,nspin,ncor,nalf,nsgn)
cgama_real = np.loadtxt("res_corr.data", skiprows=5, usecols=10, dtype = np.float64)
cgama_img = np.loadtxt("res_corr.data", skiprows=5, usecols=11, dtype = np.float64)
cgama = np.zeros((cd_real.shape), dtype=complex)
for i in range(cd_real.size):
  cgama[i] = complex(cgama_real[i], cgama_img[i])


# read ifff(nvar2,nspin,ncor,nalf,nsgn)
ifff = np.loadtxt("res_corr.data", skiprows=5, usecols=12, dtype = int)
