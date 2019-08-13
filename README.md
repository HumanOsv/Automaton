# Automaton

![alt text](https://github.com/HumanOsv/Logos/blob/master/DesignEvo.jpg)

AUTOMATON is a hybrid program that combines a probabilistic cellular automata and a genetic algorithm for the global minimum search of clusters and molecules. The main procedure consists of two major steps: first, a discrete population is generated combining rules of a simplified Cellular Automata model with geometrical optimizations (to the nearest stationary point) using an ab initio or DFT method; second, this population is primarily evolved through genetic operations, and then followed by geometrical optimizations, towards the nearest stationary point. Additionally, AUTOMATON includes a structure-recognition routine, which is used at different stages of the search process to identify and eliminate duplicate structures.

• Osvaldo Yañez, Rodrigo Báez-Grez, Diego Inostroza, Walter A. Rabanal-León, Ricardo Pino-Rios, Jorge Garza, and W. Tiznado
Journal of Chemical Theory and Computation 2019 15 (2), 1463-1475. **DOI: 10.1021/acs.jctc.8b00772**
https://pubs.acs.org/doi/abs/10.1021/acs.jctc.8b00772

# Getting Started

**1)	Step Zero**

Before starting the installation, it is important to know that AUTOMATON is not a full or autonomous software; instead, it needs an assisting program to calculate energies and perform optimizations like Gaussian, Mopac or Lammps. AUTOMATON uses these pre-installed tools for the minima local search on the potential energy surface (PES).

**1. Installing external softwares.**

  •	Mopac (http://openmopac.net/Download_MOPAC_Executable_Step2.html)

  •	Gaussian (http://gaussian.com/)

  •	Lammps (https://lammps.sandia.gov/download.html#ubuntu)
  
  NOTE: mpiexec for Lammps (https://www.mpich.org/static/docs/v3.1/www1/mpiexec.html)

**2. Installing Perl environment.**

Once the programs for both energy calculations and geometry optimizations are working correctly, the Perl environment needs to be installed as well. Perl is a highly capable, feature-rich programming language that runs on many platforms from portable to mainframes.
It can be installed from:
- https://www.perl.org/get.html

There are some additional libraries and softwares that must also be installed to allow AUTOMATON to work:

-Install CPAN modules (http://www.cpan.org/modules/INSTALL.html or https://egoleo.wordpress.com/2008/05/19/how-to-install-perl-modules-through-cpan-on-ubuntu-hardy-server/)

    user$ cpan Parallel::ForkManager
      
    user$ cpan Math::Matrix

**2)	Downloading and Installing AUTOMATON**

AUTOMATON can be directly downloaded as a zip file from the page:

-https://github.com/HumanOsv/Automaton

Alternatively, it can be downloaded using the Git tools using the following command:

    user$ git clone https://github.com/HumanOsv/Automaton.git

    user$ cd ./Automaton

**Note: before downloading using Git tools, make sure to be in your final installation path.**

We recommend to install using Git tools to update future AUTOMATON software easily. To update the program, use the following command:

	user$ git pull master
	
Alternatively, AUTOMATON could be installed as follows: choose a final installation path, and then extract the ZIP file (containing the software). Provide all the basic permissions for use and, optionally, set AUTOMATON.pl file as a system call.

**3)	Running AUTOMATON**

To run AUTOMATON the following files are necessary in the working directory:

    • Input.dat              : The AUTOMATON input file, see below for more information.

    • AUTOMATON.pl           : The executable file for structure prediction. **

    • ReaxFF file (optional) : Reactive MD-force field file of Lammps.

**Note: AUTOMATON.pl can be called from another path if correctly set**

Now, use the following commands to execute this program:

    user$  perl AUTOMATON.pl Config.in > out.log

Alternatively, the user can set AUTOMATON to run in the background using one of the following methods:

	user$ nohup perl AUTOMATON.pl Config.in > out.log
	user$ setsid perl AUTOMATON.pl Config.in > out.log

**4)	Input File**

The main input file, known as input.dat, contains all the necessary parameters for a correct calculation. Each variable is explained below.

Number of structures (3N or 5N, N = Atoms number)

    numb_conf = 60


Use of genetic operators; mating and mutation operator.

    mutations = YES
    crossing_over = YES


Rules of Cellular Automata using Chemistry intelligence (Others Elements > Hydrogen > Alkali > Alkaline Earth > Halogens)

    clever_automata = NO


Maximun Cycles and Energy (kcal/mol) AUTOMATON

    maximun_cycles = 30
    maximun_energy = 100


Chemical formula for the system ( example: H 02 Pb 03 Ca 04 ).

    chemical_formula = H 06 C 06

**NOTE: Respect the spaces.**


The size of the box (Angstroms) with format "length, width, height". AUTOMATON builds an automatic cuadricular box using the sum of all covalents radii of the system.

    box_size = 
    

Software that will be used (mopac/gaussian/lammps).

    software = gaussian

*Configuring the program for chemistry packages*


Procesor and memory (GB) that will be used for each calculation.

    core_mem = 8,8


The charge and multiplicity of the candidate.

    charge_multi = 0,1


Keywords for gaussian, mopac, or lammps

*1. Gaussian*

    header = PBE1PBE/SDDAll scf=(maxcycle=512) opt=(cartesian,maxcycle=512)

*2. Mopac*

    header = AUX LARGE PM6

*3. Lammps (ReaxFF file)*

    header = reaxxFF.Co

**General Note:** Respect the spaces of separation between the symbol "=".

    Correct : software = gaussian
    Wrong   : software=gaussian

**5) AUTOMATON outputs**

After a successful run of the program, several output files will be generated in your working directory.

	01Final_coords.xyz         : Final coordinates XYZ file format of each species ordered less energy at higher energy.
	02Duplicate_coords.xyz     : Candidates that are a duplicate of a candidate in the population, in XYZ file format.
	Output.log                 : Print summary information after each set of this many iterations.
	PostCoords-NumberCycle.xyz : Local optimization for each structure.
	PreCoords-NumberCycle.xyz  : Non-optimized structures.
	
	
	
