# JTB_CSPcode_CIM_dePillis
Fortran code for manuscript of J Theor Biol, Patsatzis 2021, under review


The present Fortran code is used for running the simulations of the above manuscript, so that
the user acquires all the CSP data for the desired simulation. The model introduced in de Pillis et al. J Theor Biol 238(4) 2006 

The code includes the following external packages, which run independantly without prior installation:
1) ODEPACK by LLNL for numerical integration of ODEs
2) EISPACK by Netlib for solving the eigenproblem

The code runs through *sh script.sh*.

Before running:
1) Edit *Makefile* to match the OS installed libraries.
2) Edit *defs.h* to select the operations you want to perform (numerical/analytic Jacobian, use/no use of the CSP tools)
3) Provide the parameter set, the initial conditions and the ending time for your simulation in *paramet.i*.

After running the code outputs the following:
1) Solution of cell populations at *Asol.dat*, rhs at *Arhs.dat*, process rates at *Arates.dat* 
2) Timescales at *ATmscl.dat*, CSP amplitudes at *AFi_PosNeg.dat* and their absolute values at *AFi_Pos.dat*
3) The CSP data: (i) number of exhausted CSP modes at *ANofExhMod.dat*, (ii) API for all modes at *AAPI.dat*, (iii) TPI for all modes *ATPI.dat*, (iv) Po for all modes *APointers.dat*, (v) II for all cell populations *AII.dat*
4) The folder ADiag, which contains all the above indices API, TPI, Po and II from (ii)-(v) sorted in descedning order and devided per mode. In additions, each *Mod1.dat* - *Mod4.dat* files contain all the CSP data sorted in descedning order per mode, but printed at selected timepoints, as defined in SUB_main.f.
5) The absolute error of the constraint solution with respect to the solution at *ASIMAbsErr.dat*

At the end of the run you can see the ending tumor cell population of the simulation at the terminal.

Please use with care and contact me for further assistance at patsatzisdim@gmail.com
