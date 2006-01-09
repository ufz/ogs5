README.TXT

                                PHREEQC

       A program for speciation, batch-reaction, one-dimensional 
            transport, and inverse geochemical calculations

NOTE: This file describes the batch version of PHREEQC. However, a
graphical user interface, PhreeqcI, is available that runs under
Windows operating systems. This interactive program can be downloaded
from http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/index.html. 
PhreeqcI has internal documentation for all PHREEQC data input and is a
complete interface to PHREEQC. A brief description of the program has
been published as a Fact Sheet:

Charlton, S.R., and Parkhurst, D.L., 2002, PhreeqcI--A graphical user
   interface to the geochemical model PHREEQC: U.S. Geological Survey 
   Fact Sheet FS-031-02, 2 p.


                  PHREEQC - Version 2.15.0  February  5, 2008

Instructions for installation, execution, and testing are provided
below.  After installation, see the phreeqc.txt file in the doc directory
of the PHREEQC installation for summary information on PHREEQC.

For assistance, enhancement requests, or to report bugs, contact 
dlpark@usgs.gov.


                          TABLE OF CONTENTS

                          A. DISTRIBUTION FILES
                          B. DOCUMENTATION
                          C. EXTRACTING FILES
                          D  COMPILING
                          E. INSTALLING
                          F. RUNNING THE PROGRAM
                          G. TESTING
                          H. CONTACTS


A. DISTRIBUTION FILES

The following distribution packages (containing the software, test data sets,
and information files) are currently available for UNIX systems:

  phreeqc-2.15.0-2697.Linux.tar.gz  - Compiled with gcc 4.1.2 20070302 (prerelease) (4.1.2-1mdv2007.1)
                                           kernel 2.6.17-16mdv
  phreeqc-2.15.0-2697.source.tar.gz - Source code but no executable


B. DOCUMENTATION

Parkhurst, D.L., and Appelo, C.A.J., 1999, User's guide to PHREEQC
   (Version 2)--A computer program for speciation, batch-reaction,
   one-dimensional transport, and inverse geochemical calculations:
   U.S. Geological Survey Water-Resources Investigations Report 99-4259,
   312 p.

   This user's guide is available in electronic format. Portable Document
   Format (PDF) files are included in the doc subdirectory of the
   PHREEQC program distribution. An online HTML version of the User's
   Guide can also be found at
   http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/

Thorstenson, D.C., and Parkhurst, D.L., 2002, Calculation of
   individual isotope equilibrium constants for implementation in
   geochemical models: U.S. Geological Survey Water-Resources
   Investigations Report 02-4172, 129 p.

   This report documents the theory and implementation of isotopes 
   in PHREEQC. Portable Document Format (PDF) is included in the 
   doc subdirectory of the PHREEQC program distribution. 

C. EXTRACTING FILES

Compressed tar files are used to distribute the source code and
versions of the software compiled for selected UNIX operating systems.
All of the PHREEQC files needed to install and test the program are
contained in the file phreeqc-2.15.0.2697-OS.tar.gz (where OS is a
string indicating the intended operating system.)  If there is not a
tar file for your operating system or you want to compile the
software, the source version of the tar file contains all of the
PHREEQC files needed to compile and install the program on a
UNIX-based computer.  For all of these distributions, the directory
phreeqc-2.15.0 is created (or overwritten) when the files are extracted
from the tar file; if this directory already exists, you may want to
delete or rename it before extracting the files.

Follow the steps below to extract the files from a distribution tar
file. 

  Steps in extracting files                 explanation
  ----------------------------------------  -----------------------------------
  gunzip phreeqc-2.15.0-2697.OS.tar.gz       Uncompress the distribution file.

  tar -xvpof phreeqc-2.15.0-2697.OS.tar      Extract files from the tar file.

The following directory structure is created (the contents of each
directory are shown to the right):

   phreeqc-2.15.0      files NOTICE.TXT, RELEASE.TXT, and this README.TXT
     `-----bin         compiled executable, template for execution script
     `-----database    database files required during execution
     `-----doc         documentation files 
     `-----examples    examples from user's guide--used in verification tests
     `-----src         Makefile and source code
     `-----test        scripts to run verification tests

Notes:  (a) A compiled executable is not included in the source
            distribution.
        (b) It is recommended that no user files be kept in the PHREEQC
            directory structure.  If you plan to put files in the PHREEQC
            directory structure, do so only by creating subdirectories.


D. COMPILING

If you have retrieved a pre-compiled distribution of the software, skip to
the Installing section below.

If a compiled version of the software is not available for your computer
or you want to build the executable yourself, follow the instructions in
this section.  The source distribution is provided for those users who
want the source code.  Little or no support is provided for users
generating their own versions of the software.  In general, to compile a
new version of the software, you will need:

       (a) a C compiler, and
       (b) a minimal level of knowledge of Make, the compiler, and
           the UNIX operating system.

As provided in the source distribution, the software is set up to be
compiled under Linux with gcc.

To generate a new executable, do the following:

1.  Change directory to the source directory:
        cd phreeqc-2.15.0/src

2.  Modify the beginning of the file named Makefile to correctly specify 
    system-dependent variables:
     
        CC             C compiler name
        CCFLAGS        C compiler flags
        LOADFLAGS      linker flags
        INVERSE_CL1MP  remove definition if gmp package 
                       is not available on your system

3.  Run make in the src directory to compile the source:

        make

    make will:

        a.  compile the source code, and
        b.  place the program executable in the PHREEQC bin directory.


E. INSTALLING

To make PHREEQC easy to use, a link to the script that runs PHREEQC
should be placed in a directory that is included in each user's search
path.  Run make in the PHREEQC src subdirectory to create the link:

    make install BINDIR=directory_for_links

A script for running the executable is placed in the PHREEQC installation
directory and a link to that script is placed in the directory defined
by BINDIR.  If each user's search path consists of

    /usr/bin:/usr/opt/bin:/usr/local/bin

using the command

    make install BINDIR=/usr/local/bin

will make PHREEQC accessible from any directory without requiring the full
pathname of the executable.  Note that to create and delete links to the
PHREEQC script, the installer must have sufficient rights to the BINDIR
directory. By default, BINDIR is $HOME/bin.


F. RUNNING THE PROGRAM

If PHREEQC has been installed in a directory included in the users'
PATH (as described above), the program can be executed with any of the
commands below. 

  command to execute PHREEQC     explanation
  -----------------------------  -----------------------------------------
  phreeqc                        The program will query for each of the
                                 needed files.

  phreeqc input                  The input file is named input, the output
                                 file will be named input.out and the
                                 default database file will be used.
  
  phreeqc input output           The input file is named input, the output
                                 file is named output, and the default
                                 database file will be used.
  
  phreeqc input output database  All file names are specified explicitly.
  
  phreeqc input output database screen_output     

                                 All file names are specified explicitly,
                                 and screen output is directed to the
                                 file screen_output.

If the "make install" command described above has not been executed,
the executable must be copied to the current directory. The executable
is normally phreeqc-2.15.0/bin/phreeqc. In addition the default database
file will be phreeqc.dat in the current directory. Phreeqc.dat is
installed in phreeqc-2.15.0/database/phreeqc.dat.

The shell variable PHREEQC_DATABASE can be used to specify the default
database. In the C shell, this variable can be set with the command:

  setenv PHREEQC_DATABASE /home/jdoe/local/project/mydata.dat

In the Bourne or Korn shell, this variable can be set with the command:

  export PHREEQC_DATABASE=/home/jdoe/local/project/mydata.dat

The shell variable can be set permanently by including the appropriate
command in a file that is read when the shell is initiated, frequently
$HOME/.login or $HOME/.profile.  If this shell variable is not set,
the default database is set in the script in the installation directory
to phreeqc-2.15.0/database/phreeqc.dat.  It is possible to specify a
different default database by editing the script.


G. TESTING

Test data sets are provided to verify that the program is correctly installed
and running on the system.  The tests execute the examples presented in
the user's guide, which demonstrate most of the program's capabilities.
The PHREEQC test directory contains the script to run the tests.  The PHREEQC
examples directory contains the input data and the expected results for
each test.

To test the installation, change to the phreeqc-2.15.0/test directory and type
the command:

     ./test.sh [start [stop]]

where:  start = the number of the first test to perform, default = 1
        stop  = the number of the last test to perform, default = 18

For example:

     command                             what happens
     ----------------------------------  --------------------------------
     ./test.sh                           runs all of the tests
     ./test.sh 1 1                       runs the first test
     ./test.sh 2 3                       runs tests 2 and 3
     ./test.sh 4                         runs test 4 through the last test


After the tests are completed, the results are compared to the expected
results (found in the PHREEQC examples directory).  See the file check.log; 
if all goes well, there should be no differences.  

To clean up after the tests, type the command:

     ./clean.sh

Notes: The output data files for phreeqc-2.15.0.2697.source.tar.gz
       were created on a Linux system. You may notice slight numeric differences 
       in the results. These are generally due to different round-off algorithms
       and the different architecture of the central processing unit chip.
       Slight differences in output formats may occur on other computers,
       particularly for the value 0.0.

       Problem 9 generates a warning message that indicates negative 
       concentrations were generated in a kinetic run. The numerical
       method automatically reduces the step size until negative 
       concentrations are eliminated and an accurate kinetic integration
       is obtained.

       Problems 11, 12, 13, and 15 generate warning messages. The
       messages simply indicate a shorthand method for defining cell
       lengths and cell dispersivities was used.

The tests are described in the table below, where 'test' is the test
number and the 'usage' column indicates how a file is used, with i for
input and o for output.

test  description of test and files                  file name & usage
----  --------------------------------------------   -----------------
  1   Add uranium and speciate seawater

      PHREEQC input file                             ex1             i
      Printed results of calculation                 ex1.out         o
                                                                      
  2   Temperature dependence of solubility                            
      of gypsum and anhydrite                                         
                                                                      
      PHREEQC input file                             ex2             i
      Printed results of calculation                 ex2.out         o
      Results written to the selected-output file    ex2.sel         o

  3   A. Calcite equilibrium at log Pco2 = -2.0 and 25C
      B. Definition of seawater
      C. Mix 70% ground water, 30% seawater
      D. Equilibrate mixture with calcite and dolomite
      E. Equilibrate mixture with calcite only

      PHREEQC input file                             ex3             i
      Printed results of calculation                 ex3.out         o
                                                                      
  4   A. Rain water evaporation                                       
      B. Factor of 20 more solution                                   
                                                                      
      PHREEQC input file                             ex4             i
      Printed results of calculation                 ex4.out         o

  5   Add oxygen, equilibrate with pyrite, calcite, and goethite.

      PHREEQC input file                             ex5             i
      Printed results of calculation                 ex5.out         o
      Results written to the selected-output file    ex5.sel         o
                                                                      
  6   6A. React to phase boundaries                                    
      6A1.--Find amount of k-feldspar dissolved to                         
            reach gibbsite saturation                                  
      6A2.--Find amount of k-feldspar dissolved to                         
            reach kaolinite saturation                                 
      6A3.--Find amount of k-feldspar dissolved to                         
            reach K-mica saturation                                 
      6A4.--Find amount of k-feldspar dissolved to                         
            reach k-feldspar saturation                                    
      6A5.--Find point with kaolinite present,                         
            but no gibbsite                                            
      6A6.--Find point with K-mica present,                         
            but no kaolinite                                           
      6B. Path between phase boundaries                                
      6C. Kinetic calculation
                                                                      
      PHREEQC input file                             ex6             i
      Printed results of calculation                 ex6.out         o
      Results written to the selected-output file    ex6A-B.sel      o
      Results written to the selected-output file    ex6C.sel        o

  7   Organic decomposition with fixed-pressure and
      fixed-volume gas phases  

      PHREEQC input file                             ex7             i
      Printed results of calculation                 ex7.out         o
      Results written to the selected-output file    ex7.sel         o
                                                                      
  8   Sorption of zinc on hydrous iron oxides                         
                                                                      
      PHREEQC input file                             ex8             i
      Printed results of calculation                 ex8.out         o
      Results written to the selected-output file    ex8.sel         o
                                                                      
  9   Kinetically controlled oxidation of ferrous iron. 
      Decoupled valence states of iron.                                      
                                                                      
      PHREEQC input file                             ex9             i
      Printed results of calculation                 ex9.out         o
      Results written to the selected-output file    ex9.sel         o
                                                                      
 10   Solid solution of strontianite and aragonite.
                                                                      
      PHREEQC input file                             ex10            i
      Printed results of calculation                 ex10.out        o
      Results written to the selected-output file    ex10.sel        o
                                                                      
 11   Transport and ion exchange.
                                                                      
      PHREEQC input file                             ex11            i
      Printed results of calculation                 ex11.out        o
      Results written to the selected-output file    ex11adv.sel     o
      Results written to the selected-output file    ex11trn.sel     o

 12   Advective and diffusive transport of heat and solutes. Constant
      boundary condition at one end, closed at other. The problem is 
      designed so that temperature should equal Na-conc (in mmol/kgw) 
      after diffusion.

      PHREEQC input file                             ex12            i
      Printed results of calculation                 ex12.out        o
      Results written to the selected-output file    ex12.sel        o

 13   A. 1 mmol/l NaCl/NO3 enters column with stagnant zones. Implicit 
         definition of first-order exchange model.
      B. 1 mmol/l NaCl/NO3 enters column with stagnant zones. Explicit 
         definition of first-order exchange factors.
      C. 1 mmol/l NaCl/NO3 enters column with stagnant zones. 5 layer 
         stagnant zone with finite differences.

      PHREEQC input file                             ex13a           i
      Printed results of calculation                 ex13a.out       o
      Results written to the selected-output file    ex13a.sel       o
      PHREEQC input file                             ex13b           i
      Printed results of calculation                 ex13b.out       o
      Results written to the selected-output file    ex13b.sel       o
      PHREEQC input file                             ex13c           i
      Printed results of calculation                 ex13c.out       o
      Results written to the selected-output file    ex13c.sel       o

 14   Transport with equilibrium_phases, exchange, and surface reactions

      PHREEQC input file                             ex14            i
      Printed results of calculation                 ex14.out        o
      Results written to the selected-output file    ex14.sel        o

 15   1D Transport: Kinetic Biodegradation, Cell Growth, and Sorption 

      PHREEQC input file                             ex15            i
      database file                                  ex15.dat        i
      Printed results of calculation                 ex15.out        o
      Results written to the selected-output file    ex15.sel        o

 16   Inverse modeling of Sierra springs

      PHREEQC input file                             ex16            i
      Printed results of calculation                 ex16.out        o

 17   Inverse modeling of Black Sea water evaporation

      PHREEQC input file                             ex17            i
      Printed results of calculation                 ex17.out        o

 18   Inverse modeling of Madison aquifer

      PHREEQC input file                             ex18            i
      Printed results of calculation                 ex18.out        o


H. CONTACTS

Inquiries about this software distribution should be directed to:

e-mail:  h2osoft@usgs.gov or dlpark@usgs.gov

