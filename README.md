# Automaton

We introduce a structural search algorithm implemented in the new **AUTOMATON** program. The program consists of two main procedures: in the first one, a discrete population is generated combining rules of a simplified Celular Automaton model with geometrical optimizations (to the nearest stationary point) using an *ab initio* method. In the second one, this population is evolved through genetic operations followed by geometrical optimizations (to the nearest stationary point) using an *ab initio* method. Additionally, **AUTOMATON** includes a structure-recognition routine, which is used in different stages of the search process to identify and eliminate duplicates.

# Getting Started

**1)	Step Zero**

Before starting to install, copying or calculating stuff is important to know that AUTOMATON isn´t a full or autosuficienct software, it needs an installed calculation software like Gaussian, Mopac or Lammps. AUTOMATON use these pre-installed tools for the minima local search.

**1. Install external softwares.**

  •	Mopac (http://openmopac.net/Download_MOPAC_Executable_Step2.html)

  •	Gaussian (http://gaussian.com/)

  •	Lammps (https://lammps.sandia.gov/download.html#ubuntu)
  
  NOTE: mpiexec for Lammps (https://www.mpich.org/static/docs/v3.1/www1/mpiexec.html)

**2. Install Perl environment.**

Once the calculation tools are correctly installed and working, the Perl environment needs to be installed as well. Perl is highly capable, feature-rich programming language that runs on many plataforms from portable to mainframes.
Can be installed from:
- https://www.perl.org/get.html

There is some libraries and software that must be installed as well fro AUTOMATON to works:

-Install CPAN modules (http://www.cpan.org/modules/INSTALL.html or https://egoleo.wordpress.com/2008/05/19/how-to-install-perl-modules-through-cpan-on-ubuntu-hardy-server/)

    user$ cpan Parallel::ForkManager
      
    user$ cpan Math::Matrix

**2)	Downloading and Installing AUTOMATON**

The program does not have a graphical user interface, it has a command line interface that is very simple to use with some instruction. AUTOMATON program interfaces with a computational program in the background, thus the program to be used has to be available. The program allows energy calculations to be performed using a wide variety of external quantum chemistry programs including Gaussian, Mopac and Lammps (ReaxFF).

The first step is always "Download the software", AUTOMATON can be directly downloaded as a zip file from the page:

-https://github.com/HumanOsv/Automaton

Alternatively can be downloaded using the Git tools using the following command:

    user$ git clone https://github.com/HumanOsv/Automaton.git

    user$ cd ./Automaton

**Note: before downloading using Git make sure to be in your final installation path.**

We recommend to install using Git for easy futures updates to the AUTOMATON software. To update the progrma simply use

	user$ git pull master
	
The next logical step is installing AUTOMATON, for thiss just choose a final installation path and extract/git the software there. Give all the minimal permission to use and, optionally set AUTOMATON.pl file as a system call. That's it, no further instruction needed.

**3)	Running AUTOMATON**

The following necessary files needs to be in the working directory:

    • Input.dat              : The AUTOMATON input file, see below for more information.

    • AUTOMATON.pl           : The executable file for structure prediction. **

    • ReaxFF file (optional) : Reactive MD-force field file of Lammps.

**Note: AUTOMATON.pl can be call from another path if correctly set**

Now use the following commands to execute this program:

    user$  AUTOMATON.pl Config.in > out.log

alternatively, the user can set AUTOMATON to run in the background using one of the following methods:

	user$ nohup perl NICSall.pl Config.in > out.log
	user$ setsid perl NICSall.pl Config.in > out.log

**4)	Input File**

The main input file known as input.dat contains all the necessary parameters for a correct calculation. Each variable is explained below.

Number of structures (3N or 5N, N = Atoms number)

    numb_conf = 60

Use of genetic operators; mating and mutation operator.

    mutations = YES
    crossing_over = YES

Chemical formula for the system ( example: H 02 Pb 03 Ca 04 ).

    chemical_formula = H 06 C 06

**NOTE: Respect the spaces of separation.**

The size of the box (Angstroms) with format "length, width, height". AUTOMATON builds an automatic cuadricular box using the sum of all covalents radii of the system.

    box_size = 
    
Software that will be used (mopac/gaussian/lammps).

    software = gaussian

*Configuring the program for chemistry packages*

Procesor and memory (GB) that will be used for each calculation.

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

**5) Outputs of AUTOMATON**

After a successful run of the program, several output files will be generated in your working directory.

	01Final_coords.xyz       : Final coordinates XYZ file format of each species ordered less energy at higher energy.
	02Duplicate_coords.xyz   : Candidates that are a duplicate of a candidate in the population, in XYZ file format.
	03Report_Output.txt      : Print summary information after each set of this many iterations.
	04Files_not_converge.txt : Summary information of files that have problems SCF convergence (only Gaussian program).
	out.log			 : Log file with all information.
