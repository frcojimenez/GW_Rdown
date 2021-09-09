### Use of KerrQNMOvertones.f90 and the .nml parameter files:
(See the comments within the .f90 code itself for more information regarding the aims and assumptions of the code, relevant references, and the meaning of each of the parameters provided in the parameter files.)

1. Compile the .f90 code to generate the executable.
2. Create or edit a parameter file following the format of the example namelist (.nml) files provided. All 16 parameters must be provided and given some value, in the same order as in the example namelists.
3. Run the executable with the chosen parameter file passed as a command-line argument: <br> \<path_to_executable\> \<path_to_parameter_file\>

If using a compiler compatible with Fortran 90 or 95 but not 2003 or higher, the get_command_argument() command used to read the path to the parameter file as a command-line argument may prevent compilation. In this case, create the parameter file first, comment out lines 183 to 187 of the .f90 code, then uncomment line 193 and provide the path to the parameter file in that line, before compiling. The executable is then to be run directly with no argument. In this scenario however, each time a different parameter file is to be used (different path/name), the .f90 code (line 193) needs to be edited again accordingly and recompiled.

<br>

Two example parameter files are provided, KerrQNMOvertones_params_example_n3.nml and KerrQNMOvertones_params_example_n8.nml. They set parameters for the computation of two s=-2, (l,m) = (2,2) Kerr quasi-normal mode (QNM) frequency and angular separation constant values: respectively an arbitrary typical case, the n=3 overtone with a dimensionless spin of 0.69, and a "difficult" case, the n=8 overtone with a low dimensionless spin of 5\*10<sup>-3</sup>. This latter case requires many more iterations to converge due to the low efficiency of Leaver's method near the (l,m,n) = (2,2,8) Scharzschild QNM (which is algebraically special). <br>
The results and computing times may be compared to those obtained with the KerrQNMOvertones.ipynb and KerrQNMOvertones.nb notebooks provided in other folders; the example parameters provided in those notebooks already correspond to those of KerrQNMOvertones_params_example_n3.nml. The Fortran code should be found to be much faster.

The code may be switched between working with double-precision (as in the version provided) or quadruple-precision real and complex numbers by editing accordingly line 78 of the .f90 code and recompiling. Quadruple precision might be needed in some cases where convergence is difficult to reach; however it would make the computation substantially slower, with execution times comparable to the Python and Mathematica codes.
