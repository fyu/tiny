    **************************************************************
                                 sparseLM 
                              version 1.3
                          By Manolis Lourakis

                     Institute of Computer Science
            Foundation for Research and Technology - Hellas
                       Heraklion, Crete, Greece
    **************************************************************


=================================== GENERAL ===================================
This is sparseLM, a copylefted C/C++ software package for arbitrarily sparse
non-linear least squares. IT implements a variant of the Levenberg-Marquardt
optimization algorithm that can efficiently deal with non-linear least squares
problems giving rise to high-dimensional Jacobians with arbitrary sparseness.
Jacobians can be encoded either in compressed row storage (CRS) or compressed
column storage (CCS) format. sparseLM can be downloaded from
http://www.ics.forth.gr/~lourakis/sparseLM

For solving the augmented normal equations arising in the course of the LM
algorithm, sparseLM includes interfaces to a wide variety of direct sparse
solvers, the list of which currently consists of HSL's MA57, MA77, MA47 & MA27,
PARDISO, SuperLU, TAUCS, UMFPACK, CHOLMOD, LDL, CSparse, SPOOLES and MUMPS.
While no direct sparse solver is universally superior in terms of performance,
the open-source sparse Cholesky factorization package CHOLMOD is reasonably
efficient in practice and is hence chosen as sparseLM's default solver. PARDISO
(which is integrated into Intel's MKL), is another good commercial alternative.

A high-level overview regarding sparseLM can be found in M.I.A. Lourakis,
"Sparse Non-linear Least Squares Optimization for Geometric Vision", Proceedings
of European Conference on Computer Vision, vol. 2, pages 43-56, 2010.
More detailed descriptions will follow in future publications.

In case that you use sparseLM in your published work, you are requested to
include a reference to the above article:

@inproceedings{lourakis10,
    author    = {Manolis I.A. Lourakis},
    title     = {Sparse Non-linear Least Squares Optimization for Geometric Vision},
    booktitle = {European Conference on Computer Vision},
    volume    = {2},
    year      = {2010},
    pages     = {43-56},
    doi       = {http://dx.doi.org/10.1007/978-3-642-15552-9_4},
}


==================================== FILES ====================================
splm.c: core sparse Levenberg-Marquardt implementation
splm_ccsm.c: CCS sparse matrix manipulation routines
splm_crsm.c: CRS sparse matrix manipulation routines
splm_hessian.c: approximate sparse Hessian J^t*J computation
splm_stm.c: Sparse triplet matrix manipulation routines
splm_time.c: time measurement utility functions
splm_misc.c: miscellaneous functions
splm.h: Function prototypes & related data structures
splm_priv.h: Private prototypes and definitions for sparseLM
splm_cholmod.c: interface to CHOLMOD
splm_ldlp.c: interface to LDL
splm_csparse.c: interface to CSparse
splm_ma27.c: interface to MA27
splm_ma47.c: interface to MA47
splm_ma57.c: interface to MA57
splm_mumps.c: interface to MUMPS
splm_pardiso.c:interface to PARDISO
splm_dss.c: interface to MKL's DSS interface to PARDISO
splm_spooles.c: interface to SPOOLES
splm_superlu.c: interface to SuperLU
splm_taucs.c: interface to TAUCS
splm_umfpack.c: interface to UMFPACK
ma77wrap.f90 splm_ma77.c: interface to MA77

splmdemo.c: demo program with examples of sparseLM use
matlab/*: Matlab MEX-file interface and demo

============================= COMPILING & LINKING =============================
sparseLM has been developed under Linux and Cygwin with gcc, therefore most of
the instructions that follow apply to these environments. Compiling on any other
Unix-like system should be straightforward. The software also compiles under
MS Windows, however its deployment in this case becomes more complicated due
to the fact that most of the required direct solvers do not support compilation
in such a setting.
Before compiling sparseLM itself, a number of other packages on which it is
dependent upon should be compiled. More details are provided below:

 - First, the desired direct sparse solvers should be specified in solvers.mk.
   If not already installed, the specified solvers should be compiled following
   the appropriate steps for each, as described in their documentation. As
   mentioned above, the suggested solver is CHOLMOD which is part of the
   SuiteSparse bundle that can be obtained from
   http://www.cise.ufl.edu/research/sparse/SuiteSparse

   Unfortunately, the official distribution of SuiteSparse supports compilation
   only on Linux/Cygwin/Unix systems, not including any support for MS Windows.
   Nevertheless, compilation instructions for SuiteSparse with MS Visual Studio
   based on a workaround using Cygwin are at
   http://matrixprogramming.com/2008/05/umfpack-vc
   This site has also instructions for compiling MUMPS and TAUCS on Windows.
   Another option is to download CHOLMOD precompiled for WinXP from
   http://sites.google.com/site/kevinkaixu/publications/fast/code/cholmod.zip

 - Any additional packages required by the chosen solver(s) should be installed
   as well. For instance, both CHOLMOD and HSL's MA57 solvers require METIS from
   http://www.cs.umn.edu/~metis
   Please refer to each solver's documentation for more details on required
   additional software.

 - LAPACK and BLAS (or equivalent vendor libraries such as MKL) should also be
   installed since they are needed by all solvers. If not already available, they
   can be obtained from http://www.netlib.org/clapack and http://www.netlib.org/cblas,
   respectively.

 - sparseLM includes makefiles for GNU make (i.e., Makefile/Makefile_demo) that
   can be used on Linux/Cygwin/Unix and MS nmake (i.e., Makefile.vc/Makefile_demo.vc)
   that can be used with MS Visual C under Windows. The compilation instructions
   that follow should be applied on the makefile that is appropriate for your
   target system. Alternatively, the included CMake configuration file can be
   used to generate suitable build files; more details can be found below.


 - After the installation of the sparse solvers, the include paths for their
   header files should be specified in the configuration section at the top of
   Makefile. Note that only the paths corresponding to the chosen solver(s) are 
   relevant here; make variables corresponding to unused solvers can point to any
   path.

 - In order to link your code against sparseLM, the link paths to the libraries
   corresponding to the chosen solver(s) and any additional packages should be
   specified. An example of doing this is provided in the configuration section
   of the demo program in Makefile_demo.
   
 - Finally, typing "make" ("nmake /f Makefile.vc" on Windows) will build both
   sparseLM and the demo program.
 
 - sparseLM can also be built using the CMake cross-platform build system.
   The included CMakeLists.txt file can be used to generate makefiles for 
   Unix systems or project files for Windows systems that build both sparseLM
   and the demo program. CMakeLists.txt defines several configuration variables
   that control which sparse direct solvers are available and the paths where
   they have been installed. These variables can be modified from CMake's user
   interface. For example, to compile support for the CHOLMOD linear solver,
   variable HAVE_CHOLMOD should be set to TRUE and CHOLMOD_INCDIR should be made
   to point to the directory where CHOLMOD's header files are installed. Defining
   these two variables suffices to compile the sparseLM library, however linking
   against it requires two additional variables to be defined: CHOLMOD_LIBDIR,
   which should contain the path to CHOLMOD's libraries and CHOLMOD_LIBS which
   should contain their actual names. Support for any other linear solver can
   be compiled in a similar manner.
   More information on how to use CMake can be found at http://www.cmake.org

============================= MATLAB INTERFACE =============================
A matlab interface to sparseLM is included in the 'matlab' subdirectory of
the distribution. Please refer to it for more information and examples of use.


Send your comments/bug reports to lourakis (at) ics (dot) forth (dot) gr
