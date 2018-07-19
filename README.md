# Automaton

We introduce a structural search algorithm implemented in the new AUTOMATON program. The program consists of two main procedures: in the first one, a discrete population is generated combining rules of a simplified Celular Automaton model with geometrical optimizations (to the nearest stationary point) using an ab initio method. In the second one, this population is evolved through genetic operations followed by geometrical optimizations (to the nearest stationary point) using an ab initio method. Additionally, AUTOMATON includes a structure-recognition routine, which is used in different stages of the search process to identify and eliminate duplicates.

# Getting Started

1)	Prerequisites

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

2)	Running AUTOMATON

To download the AUTOMATON you need Git installed on your computer. If Git is installed use the following command to download the AUTOMATON: 

    user$ git clone https://github.com/HumanOsv/Automaton.git

    user$ cd ./Automaton

The following necessary files should appear in the working directory:

    • Config.in              : The AUTOMATON input file

    • AUTOMATON.pl           : The executable file for structure prediction

    • ReaxFF file (optional) : Reactive MD-force field file of Lammps

Now use the following commands to execute this program:

    user$ setsid perl AUTOMATON.pl Config.in >out.log

