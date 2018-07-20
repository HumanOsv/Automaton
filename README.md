# Automaton

We introduce a structural search algorithm implemented in the new **AUTOMATON** program. The program consists of two main procedures: in the first one, a discrete population is generated combining rules of a simplified Celular Automaton model with geometrical optimizations (to the nearest stationary point) using an ab initio method. In the second one, this population is evolved through genetic operations followed by geometrical optimizations (to the nearest stationary point) using an *ab initio* method. Additionally, **AUTOMATON** includes a structure-recognition routine, which is used in different stages of the search process to identify and eliminate duplicates.

# Getting Started

**1)	Prerequisites**

AUTOMATON is written in Perl. The program has only been tested on Mac OS X, Linux and Windows, so it can’t be guaranteed to run on other operating systems.

This topic lists library and software that must be installed prior to installing AUTOMATON.

-Install CPAN modules (http://www.cpan.org/modules/INSTALL.html)

    user$ sudo cpan Parallel::ForkManager
      
    user$ sudo cpan Math::Matrix

-Install External softwares

  •	Mopac (http://openmopac.net/Download_MOPAC_Executable_Step2.html)

  •	Gaussian (http://gaussian.com/)

  •	Lammps (https://lammps.sandia.gov/download.html#ubuntu)
  
  NOTE: mpiexec for Lammps (https://www.mpich.org/static/docs/v3.1/www1/mpiexec.html)

**2)	Running AUTOMATON**

The program does not have a graphical user interface, it has a command line interface that is very simple to use with some instruction. AUTOMATON program interfaces with a computational program in the background, thus the program to be used has to be available. The program allows energy calculations to be performed using a wide variety of external quantum chemistry programs including Gaussian, Mopac and Lammps (ReaxFF).

To download the AUTOMATON you need Git installed on your computer. If Git is installed use the following command to download the AUTOMATON: 

    user$ git clone https://github.com/HumanOsv/Automaton.git

    user$ cd ./Automaton

The following necessary files should appear in the working directory:

    • Config.in              : The AUTOMATON input file

    • AUTOMATON.pl           : The executable file for structure prediction

    • ReaxFF file (optional) : Reactive MD-force field file of Lammps

Now use the following commands to execute this program:

    user$ setsid perl AUTOMATON.pl Config.in >out.log

After a successful run of the program, several output files named as: 01Final_coords.xyz, 02Duplicate_coords.xyz, 03Report_Output.txt and 04Files_not_converge.txt (only Gaussian) will be generated in your working directory.

	01Final_coords.xyz       : Final coordinates XYZ file format of each species ordered less energy at higher energy.
	02Duplicate_coords.xyz   : Candidates that are a duplicate of a candidate in the population, in XYZ file format.
	03Report_Output.txt      : Print summary information after each set of this many iterations.
	04Files_not_converge.txt : Summary information of files that have problems SCF convergence (only Gaussian program).


**3)	Input File**

The main input file named as input.dat, contains all necessary parameters for the structure
prediction.

Number of structures (3N or 5N, N = Atoms number)

    numb_conf = 60

Genetic operations Most genetic algorithms implement several genetic operators; mating
and mutation operator.

    mutations = YES
    crossing_over = YES

The size of the box (in Angstroms) length, width, and height. AUTOMATON build an automatic box.

    box_size = 

Different atomic (or chemical) species in the system.( example: H 02 Pb 03 Ca 04 )

    chemical_formula = H 06 C 06
    
*NOTE: Respect the spaces of separation.*

Software mopac and gaussian (mopac/gaussian/lammps)

    software = gaussian

*Configuring the program for chemistry packages*

The number of processors to use in the run (the value may be used to create the input file) # and memory to be used in GB.

    core_mem = 8,8

The charge and multiplicity of the candidate.

    charge_multi = 0,1

keywords for gaussian, mopac, or lammps

*Gaussian*

    header = PBE1PBE/SDDAll scf=(maxcycle=512) opt=(cartesian,maxcycle=512)

*Mopac*

    header = AUX LARGE PM6

*Lammps (ReaxFF file)*

    header = reaxxFF.Co

General Note: Respect the spaces of separation between the symbol "=".

    Correct : software = gaussian
    Wrong   : software=gaussian
	
	


