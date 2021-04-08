import sys
import cf
import numpy as np

day = int(sys.argv[1])

ol = cf.read("/neodc/esacci/sst/data/CDR_v2/Climatology/L4/v2.1/D%03d-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc"%(day))
nl = cf.read("clim/D%03d-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc"%(day))
diffs = nl[0].array - ol[0].array
print(np.min(diffs),np.max(diffs),np.mean(diffs))
mdiffs = np.bitwise_xor(nl[0].array.mask,ol[0].array.mask)
np.any(mdiffs)

print(np.histogram(diffs))
