import numpy as np
import yt
from yt.extensions.astro_analysis.halo_analysis import HaloCatalog
import yt.units as u
import os
import sys
import glob
yt.enable_parallelism()

sim_dir = sys.argv[-1]
all_snaps = list(glob.glob(sim_dir+"/RD00*"))
for snap in all_snaps:
    snap_num = snap[-2:]

    fname = snap+"/RedshiftOutput00"+snap_num
    print("Loading ", fname)
    data_ds = yt.load(fname)

    print("Creating HaloCatalog()")
    hc=HaloCatalog(data_ds = data_ds, finder_method="fof", output_dir=sim_dir+"/fof_halos")
    print("Running")
    hc.create()
