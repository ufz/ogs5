mHM2OGS is a command line tool to process the recharge data of mHM for the groundwater modelling by using OGS5

It takes the recharge data of mHM, finds face elements on the top surface of 3D domain, and performs face integration to transform the mHM recharge data to the nodal flux for the finite element analysis in OGS.

In addition  to the data files from mHM, two more input files with the same base file name are needed to run the command. One input file is for mesh, and the other file is for the recharge data files from mHM. The later has an entension  name of 'pcp', and it contains data of time unit, a factor of ratio, and the names of recharge data files. The syntax of pcp file looks like:

---
TimeUnit: Day
Ratio: 0.8
mHM_sim_Percolation_1991_1.asc
mHM_sim_Percolation_1991_2.asc
mHM_sim_Percolation_1991_3.asc
---

In pcp file, Ratio is a factor that is mutilplied to the recharge data. One can use it as unit convertion factor. Note: in OGS5, the unit of velocity is  m/[time unit].

The tool is run in command line as
 mHM2OGS [file name] [option]

options are:
---
    --version:  display version number.
    --help:     display help info.
    --output-directory: set output directory.
---


