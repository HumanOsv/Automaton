#!/usr/bin/perl


#
# Write code: 
#            Diego Inostrosa & Osvaldo Yañez Osses
#            contact: osvyanezosses@gmail.com

use strict;
use warnings;
use Benchmark;
use List::Util qw( min max );
use Math::Trig;

#
# install: sudo cpan Parallel::ForkManager
use Parallel::ForkManager;
#
# install: sudo cpan Math::Matrix
use Math::Matrix;









# Numero de colas para correr en 
# paralelo (Ej: Si mandas un calculo de 4 procesadores y tu 
#               computador es de 16, por lo tanto son 4 colas 
#                puesto que 4*4 = 16)
my $numero_colas     = 3;
my $gaussian_version = "g09";
#
my $exec_bin_g09     = "Gaussian16.b01";
#
# Percentage of initial population 
my $Percentage_1D = 10;
my $Percentage_2D = 30;
my $Percentage_3D = 60;

# # # 
# Check cluster send jobs
# Local = 0 y queue = 1 only Gaussian 
my $local_cluster   = 0;
#
my ($file_name,$queue);
if ($local_cluster == 1) { 
	($file_name,$queue) = @ARGV;
} else {
	($file_name) = @ARGV;
	$queue = 0;
}
#
# # #
# Config Automaton
my $MaximunEnergy            = 0;
my $Maximun_Cycles_Automaton = 0;
# Fix Value 0.3
my $cell_size_w = 0.3;
# For connect (Delta mayor al Grid)
my $delta       = 0.4;
# Mutations type of Grid 1D, 2D and 3D
my $option_1D = 0;
my $option_2D = 1;
my $option_3D = 1;
# porcentaje de atomos mutados
my $PCENTATOMMUTATED   = 0.3;
my $kick_punch         = 0.4;
# New species (1/5)*N
my $New_Pobl_Mut      = (1/5);
my $New_Pobl_Cell     = (1/5);
# Restart file
my $file_restart      = "Fit_Cell.xyz";
# Number of maximum cycles in automaton 
my $max_numb_convergence = 9;
# Check similarity
my $threshold_duplicate  = 0.006;
# Number of process
my $nprocess             = 100;

# # #
# For gaussian
my $max_numb_conver  = 3;
my $max_numb_NImag   = 3;

# # # 
# convergence criteria Lammps
my $criteria        = "1e-06";
# specify the maximum numer of steps 
my $steps           = "100";
# For linux
my $path_bin_lammps = "mpiexec -np 2 lammps-daily";
# For Windows
#my $path_bin_lammps = "lmp_mpi.exe";

# # # 
# Path for software MOPAC
my $path_bin_mopac       = "/opt/mopac/MOPAC2016.exe";




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
my $num_atoms_xyz;
# Global array
my @tmp_new_neighbors_atoms = ();
#
my $count_struc_global   = 0;
my $count_struc_mut      = 0;
my $count_struc_mut_Ex   = 0;
my $count_struc_crossing = 0;
#
my $Global_Port_Normal  = 0;
my $Global_Port_Error   = 0;
my $Global_Port_Corrupt = 0;
my $Global_No_Coords    = 0;
#
my @files_Port_Error   = ();
my @files_No_Coords    = ();
#
my $Geom_Exchange_Option;
#
# Vector normal del plano.
my @v = ();
# Vector ortogonal 
my @v_orth;
# Indica separacion entre fragmento y plano.
my $spacingInPlane = 0.6; 
# Hash contador de atomos x tipo
my %atom_counting = (default=>0);
# tmp
my $separation = "";

##############
# Hashs 

#  Cordero, B.; Gómez, V.; Platero-Prats, A. E.; Revés, M.; Echeverría, J.; Cremades, E.; Barragán, F.; Alvarez, S. 
#  Covalent Radii Revisited. J. Chem. Soc. Dalt. Trans. 2008, No. 21, 2832–2838. DOI: 10.1039/B801115J 
my $other_element = 0.8;

my %Atomic_radii = ( 'H'  => '0.31', 'He' => '0.28', 'Li' => '1.28', 'Be' => '0.96',
                     'B'  => '0.84', 'C'  => '0.76', 'N'  => '0.71', 'O'  => '0.66',
                     'F'  => '0.57', 'Ne' => '0.58', 'Na' => '1.66', 'Mg' => '1.41',
                     'Al' => '1.21', 'Si' => '1.11', 'P'  => '1.07', 'S'  => '1.05',
                     'Cl' => '1.02', 'Ar' => '1.06', 'K'  => '2.03', 'Ca' => '1.77',
                     'Sc' => '1.70', 'Ti' => '1.60', 'V'  => '1.53', 'Cr' => '1.39',
                     'Mn' => '1.39', 'Fe' => '1.32', 'Co' => '1.26', 'Ni' => '1.24',
                     'Cu' => '1.32', 'Zn' => '1.22', 'Ga' => '1.22', 'Ge' => '1.20',
                     'As' => '1.19', 'Se' => '1.20', 'Br' => '1.20', 'Kr' => '1.16',
                     'Rb' => '2.20', 'Sr' => '1.95', 'Y'  => '1.90', 'Zr' => '1.75',
                     'Nb' => '1.64', 'Mo' => '1.54', 'Tc' => '1.47', 'Ru' => '1.46',
                     'Rh' => '1.42', 'Pd' => '1.39', 'Ag' => '1.45', 'Cd' => '1.44',
                     'In' => '1.42', 'Sn' => '1.39', 'Sb' => '1.39', 'Te' => '1.38',
                     'I'  => '1.39', 'Xe' => '1.40', 'Cs' => '2.44', 'Ba' => '2.16',
                     'La' => '2.07', 'Ce' => '2.04', 'Pr' => '2.03', 'Nd' => '2.01',
                     'Pm' => '1.99', 'Sm' => '1.98', 'Eu' => '1.98', 'Gd' => '1.96',
                     'Tb' => '1.94', 'Dy' => '1.92', 'Ho' => '1.92', 'Er' => '1.89',
                     'Tm' => '1.90', 'Yb' => '1.87', 'Lu' => '1.87', 'Hf' => '1.75',
                     'Ta' => '1.70', 'W'  => '1.62', 'Re' => '1.51', 'Os' => '1.44',
                     'Ir' => '1.41', 'Pt' => '1.36', 'Au' => '1.36', 'Hg' => '1.32',
                     'Tl' => '1.45', 'Pb' => '1.46', 'Bi' => '1.48', 'Po' => '1.40',
                     'At' => '1.50', 'Rn' => '1.50', 'Fr' => '2.60', 'Ra' => '2.21',
                     'Ac' => '2.15', 'Th' => '2.06', 'Pa' => '2.00', 'U'  => '1.96',
                     'Np' => '1.90', 'Pu' => '1.87', 'Am' => '1.80', 'Cm' => '1.69' );

my %Atomic_number = ( '89'  => 'Ac', '13'  => 'Al', '95'  => 'Am', '51'  => 'Sb',
		      '18'  => 'Ar', '33'  => 'As', '85'  => 'At', '16'  => 'S',  
		      '56'  => 'Ba', '4'   => 'Be', '97'  => 'Bk', '83'  => 'Bi',	
		      '107' => 'Bh', '5'   => 'B',  '35'  => 'Br', '48'  => 'Cd',	
		      '20'  => 'Ca', '98'  => 'Cf', '6'   => 'C',  '58'  => 'Ce',	
		      '55'  => 'Cs', '17'  => 'Cl', '27'  => 'Co', '29'  => 'Cu',	
		      '24'  => 'Cr', '96'  => 'Cm', '110' => 'Ds', '66'  => 'Dy',
		      '105' => 'Db', '99'  => 'Es', '68'  => 'Er', '21'  => 'Sc',	
		      '50'  => 'Sn', '38'  => 'Sr', '63'  => 'Eu', '100' => 'Fm',	
		      '9'   => 'F',  '15'  => 'P',  '87'  => 'Fr', '64'  => 'Gd',	
		      '31'  => 'Ga', '32'  => 'Ge', '72'  => 'Hf', '108' => 'Hs',	
		      '2'   => 'He', '1'   => 'H',  '26'  => 'Fe', '67'  => 'Ho',	
		      '49'  => 'In', '53'  => 'I',  '77'  => 'Ir', '70'  => 'Yb',
		      '39'  => 'Y',  '36'  => 'Kr', '57'  => 'La', '103' => 'Lr',	
		      '3'   => 'Li', '71'  => 'Lu', '12'  => 'Mg', '25'  => 'Mn',	
		      '109' => 'Mt', '101' => 'Md', '80'  => 'Hg', '42'  => 'Mo',	
		      '60'  => 'Nd', '10'  => 'Ne', '93'  => 'Np', '41'  => 'Nb',	
		      '28'  => 'Ni', '7'   => 'N',  '102' => 'No', '79'  => 'Au',	
		      '76'  => 'Os', '8'   => 'O',  '46'  => 'Pd', '47'  => 'Ag',	
		      '78'  => 'Pt', '82'  => 'Pb', '94'  => 'Pu', '84'  => 'Po',	
		      '19'  => 'K',  '59'  => 'Pr', '61'  => 'Pm', '91'  => 'Pa',	
		      '88'  => 'Ra', '86'  => 'Rn', '75'  => 'Re', '45'  => 'Rh',	
		      '37'  => 'Rb', '44'  => 'Ru', '104' => 'Rf', '62'  => 'Sm',
		      '106' => 'Sg', '34'  => 'Se', '14'  => 'Si', '11'  => 'Na',
		      '81'  => 'Tl', '73'  => 'Ta', '43'  => 'Tc', '52'  => 'Te',	
		      '65'  => 'Tb', '22'  => 'Ti', '90'  => 'Th', '69'  => 'Tm',	
		      '112' => 'Uub','116' => 'Uuh','111' => 'Uuu','118' => 'Uuo',	
		      '115' => 'Uup','114' => 'Uuq','117' => 'Uus','113' => 'Uut',
		      '92'  => 'U',  '23'  => 'V',  '74'  => 'W',  '54'  => 'Xe',
		      '30'  => 'Zn', '40'  => 'Zr' );

my %Atomic_mass   = ( 'H'   => '1.0079'  ,'He' => '4.003'   ,'Li'  => '6.941'   ,'Be'  => '9.0122',
                      'B'   => '10.811'  ,'C'  => '12.018'  ,'N'   => '14.0067' ,'O'   => '15.9994', 
                      'F'   => '18.998'  ,'Ne' => '20.179'  ,'Na'  => '22.9897' ,'Mg'  => '24.305',
                      'Al'  => '26.981'  ,'Si' => '28.085'  ,'P'   => '30.9738' ,'Cl'  => '35.453',
                      'K'   => '39.098'  ,'Ar' => '39.948'  ,'Ca'  => '40.078'  ,'Sc'  => '44.9559',
                      'Ti'  => '47.867'  ,'V'  => '50.942'  ,'Cr'  => '51.9961' ,'Mn'  => '54.938',
                      'Fe'  => '55.845'  ,'Ni' => '58.693'  ,'Co'  => '58.9332' ,'Cu'  => '63.546',
                      'Zn'  => '65.390'  ,'Ga' => '69.723'  ,'Ge'  => '72.64'   ,'As'  => '74.9216', 
                      'Se'  => '78.960'  ,'Br' => '79.904'  ,'Kr'  => '83.8'    ,'Rb'  => '85.4678', 
                      'Sr'  => '87.620'  ,'Y'  => '88.906'  , 'Zr' => '91.224'  ,'Nb'  => '92.9064',
                      'Mo'  => '95.940'  ,'Tc' => '98.000'  ,'Ru'  => '101.07'  ,'Rh'  => '102.9055',
                      'Pd'  => '106.420' ,'Ag' => '107.868' , 'Cd' => '112.411' ,'In'  => '114.818',
                      'Sn'  => '118.710' ,'Sb' => '121.760' ,'I'   => '126.9045','Te'  => '127.6',
                      'Xe'  => '131.290' ,'Cs' => '132.906' ,'Ba'  => '137.327' ,'La'  => '138.9055',
                      'Ce'  => '140.116' ,'Pr' => '140.908' ,'Nd'  => '144.24'  ,'Pm'  => '145',
                      'Sm'  => '150.360' ,'Eu' => '151.964' ,'Gd'  => '157.25'  ,'Tb'  => '158.9253' ,
                      'Dy'  => '162.500' ,'Ho' => '164.930' , 'Er' => '167.259' ,'Tm'  => '168.9342',
                      'Yb'  => '173.040' ,'Lu' => '174.967' ,'Hf'  => '178.49'  ,'Ta'  => '180.9479',
                      'W'   => '183.840' ,'Re' => '186.207' ,'Os'  => '190.23'  ,'Ir'  => '192.217',
                      'Pt'  => '195.078' ,'Au' => '196.967' ,'Hg'  => '200.59'  ,'Tl'  => '204.3833',
                      'Pb'  => '207.200' ,'Bi' => '208.980' ,'Po'  => '209'     ,'At'  => '210',
                      'Rn'  => '222.000' ,'Fr' => '223.000' ,'Ra'  => '226'     ,'Ac'  => '227',
                      'Pa'  => '231.035' ,'Th' => '232.038' ,'Np'  => '237'     ,'U'   => '238.0289',
                      'Am'  => '243.000' ,'Pu' => '244'     ,'Cm'  => '247'     ,'Bk'  => '247', 
                      'Cf'  => '251.000' ,'Es' => '252'     ,'Fm'  => '257'     ,'Md'  => '258',
                      'No'  => '259.000' ,'Rf' => '261'     ,'Lr'  => '262'     ,'Db'  => '262',
                      'Bh'  => '264.000' ,'Sg' => '266'     ,'Mt'  => '268'     ,'Hs'  => '277' );


#
sub ReadXYZFile {
	my ($file, $tmpOption)=@_;
	my @at = ();
	my @cx = ();
	my @cy = ();
	my @cz = ();
	my @lines = @{$file};
	# Las demas lineas tienen las coordedanas.
	foreach my $line(@lines){
		my ($a, $x, $y, $z)=split(" ", $line);
		# This option works only for 1 recombination, more than 1 will duplicate the 
		# amount of maximmum atoms in memory
		if ($tmpOption == 1) {
			$atom_counting{$a}++;
			$atom_counting{"all"}++;
		}
		push @at, $a;
		push @cx, $x;
		push @cy, $y;
		push @cz, $z;
	}
	my @coords=([@at],
				[@cx],
				[@cy],
				[@cz]);
	return @coords;
}
# The geometries are cut by a random cutting plane
sub CreateRandomPlane {
	my $a1=(rand(3))-(rand(3));
	my $b1=(rand(3))-(rand(3));
	my $c1=(rand(3))-(rand(3));

	my $a2=(rand(3))-(rand(3));
	my $b2=(rand(3))-(rand(3));
	my $c2=(rand(3))-(rand(3));	
	# vector en plano
	@v = ($b1*$c2-$c1*$b2, 
		  $c1*$a2-$a1*$c2,
		  $a1*$b2-$b1*$a2,0);
	# Separacion del plano
	my $module=sqrt(($v[0]**2)+($v[1]**2)+($v[2]**2)); 
	my $newX=$v[0]*($spacingInPlane/$module);
	my $newY=$v[1]*($spacingInPlane/$module);	
	my $newZ=$v[2]*($spacingInPlane/$module);
	@v_orth=($newX,$newY,$newZ);	
}					  
#
sub PlanePosition {
	my ($type, @coords) = @_ ;
	my @cx       = @{$coords[1]}; 
	my @cy       = @{$coords[2]}; 
	my @cz       = @{$coords[3]}; 
	my @a        = @{$coords[0]};	
	my @dots     = ();
	my %dotHash;
	for (my $var = 0; $var <= $#cx; $var++) {
		my @v_coordinates = ($cx[$var],$cy[$var],$cz[$var],1);
		my $dot           = dotprod(\@v_coordinates,\@v);
		#security measure		
		$dot              = $dot.$var;
		# separacion
		if( $type == 1 ) {
			$v_coordinates[0]=$v_coordinates[0]+$v_orth[0];
			$v_coordinates[1]=$v_coordinates[1]+$v_orth[1];
			$v_coordinates[2]=$v_coordinates[2]+$v_orth[2];
			#@v_coordinates=@v_coordinates+$spacingInPlane;
		} else {
			$v_coordinates[0]=$v_coordinates[0]-$v_orth[0];
			$v_coordinates[1]=$v_coordinates[1]-$v_orth[1];
			$v_coordinates[2]=$v_coordinates[2]-$v_orth[2];
		}		
		# hash(Dot product)=Atom, x, y, z	
		$dotHash{$dot}=[$a[$var],@v_coordinates];
		#array of dots products		
		push @dots, $dot;			
	}
	# hash type DotProduct->Coordinates, Array of DotsProduct
	return (\%dotHash, \@dots);
}  
#
sub CombineMolecule {
	# This function combine 2 molecules taking the dotproduct and coordinates.
	my ($hash1, $hash2, $sort1 , $sort2)=@_;
	# $hash is the hash variable (by reference) with dotProdut->coordinates; 
	# $sort is the dotproduct sorted in both sizes, molecule 1 will be from 
	# greater to least and molecule 2 will be from least to greater value	
	# atom_counting is a global variable, has the information of how many atoms 
	# from an specific type there are
	my $natom = $atom_counting{"all"};
	#
	delete $atom_counting{"all"};
	# hash tmp for later use
	my %verifyCorrectNumberOfAtoms;		
	my @finalCoordinates = ();
	# desreference arrays
	my @dotSortPlusMinus = @{$sort1};
	my @dotSortMinusPlus = @{$sort2};
	# desrereference hashs
	my %coordsByHash1    = %{$hash1};
	my %coordsByHash2    = %{$hash2};
	# get al the jeys from the global hash
	my @keys_atoms       = keys %atom_counting;
	# constructor for the tmp hash, 
	# all values are set to 0
	foreach my $keyInMemory(@keys_atoms){
		$verifyCorrectNumberOfAtoms{$keyInMemory} = 0;
	}
	#
	my ($i, $minor, $big) = (0, 0, 0);	
	for (my $i = 0; $i <= $natom*2; $i++) {
		my @tmpCoord;
		# greater dot value from molecule 1
		if(($i%2) == 0){
			@tmpCoord = $coordsByHash1{$dotSortPlusMinus[$big]};
			$big++;
			# lesser dot value from molecule 2
		} else {
			@tmpCoord = $coordsByHash2{$dotSortMinusPlus[$minor]};
			$minor++;
		}
		my $typeatom  = $tmpCoord[0][0];
		# if the number of a type of atoms is lesser than the maximum 
		# value for the specific atom, 
		# then the coordinates are save to the final array
		if($verifyCorrectNumberOfAtoms{$typeatom} < $atom_counting{$typeatom}){
			$verifyCorrectNumberOfAtoms{$typeatom}++;
			push @finalCoordinates, @tmpCoord;
			#
			if(($#finalCoordinates+1) == $natom){
				last;
			}
		}
	}
	# Simple print to show coordinates
	my $string_recom;
	for (my $var = 0; $var <= $#finalCoordinates; $var++) {
		#
		my $name_atom = $finalCoordinates[$var][0];
		#
		my $c_x = sprintf '%.5f', $finalCoordinates[$var][1];
		my $c_y = sprintf '%.5f', $finalCoordinates[$var][2];	
		my $c_z = sprintf '%.5f', $finalCoordinates[$var][3];
		#
		$string_recom.= "$name_atom  $c_x  $c_y  $c_z\n";
	}
	delete @atom_counting{@keys_atoms};
	return $string_recom;
}
#
sub SortDotMolecule {
	# Simple function to sort an array
	my ($dot, $option) = @_;
	my @dots   = @{$dot};
	my @sorted =();
	if($option == 1){
		#mayor a menor
		@sorted = sort {$b <=> $a} @dots
	}else{
		#menor a mayor
		@sorted = sort {$a <=> $b} @dots
	}
	return @sorted;
}
#
sub dotprod {
	# Dot Product, takes 2 vector and returns dotProduct
    my($vec_a, $vec_b) = @_;
    die "they must have the same size\n" unless @$vec_a == @$vec_b;
    my $sum = 0;
    $sum += $vec_a->[$_] * $vec_b->[$_] for 0..$#$vec_a;
    return $sum;
}
# 
sub CenterMolecule {
	my ($input)=@_;
	my @problematic_molecule=@{$input};
	my @coord_x=@{$problematic_molecule[1]};
	my @total_array = ();
	#
	my @array_center_mass = measure_center(\@{$problematic_molecule[1]},\@{$problematic_molecule[2]},\@{$problematic_molecule[3]});
	my @array_vecinvert   = vecinvert(\@array_center_mass);
	my @array_catersian   = vecadd (\@{$problematic_molecule[1]},\@{$problematic_molecule[2]},\@{$problematic_molecule[3]},\@array_vecinvert);
	@total_array=(${problematic_molecule[0]},
					$array_catersian[0],
					$array_catersian[1],
					$array_catersian[2]);
	return (@total_array);	
}
###################################
# phi, theta, psi
sub gen_ptp {
	my $pi     = 3.14159265;
	my $phi    = sprintf '%.6f', rand()*2*$pi;
	my $theta  = sprintf '%.6f', rand()*2*$pi;
	my $psi    = sprintf '%.6f', rand()*2*$pi;
	my @angles = ($phi, $theta, $psi);
	return @angles;
}
#
sub cluster_rotation {	
	my (@coords_atom) = @_;
	my @atom          = @{$coords_atom[0]}; 
	my @crds_x        = @{$coords_atom[1]}; 
	my @crds_y        = @{$coords_atom[2]}; 
	my @crds_z        = @{$coords_atom[3]}; 
	# get a set of angles
	my @base_angles = gen_ptp();
	my $phi         = $base_angles[0];
	my $theta       = $base_angles[1];
	my $psi         = $base_angles[2];
	# do the trig
	my $cos_phi     = sprintf '%.6f', cos($phi);
	my $cos_theta   = sprintf '%.6f', cos($theta);
	my $cos_psi     = sprintf '%.6f', cos($psi);
	my $sin_phi     = sprintf '%.6f', sin($phi);
	my $sin_theta   = sprintf '%.6f', sin($theta);
	my $sin_psi     = sprintf '%.6f', sin($psi);
	# make the rotation matrix
	my $D = new Math::Matrix ([$cos_phi,$sin_phi,0],[-$sin_phi,$cos_phi,0],[0,0,1]);
	my $C = new Math::Matrix ([1,0,0],[0,$cos_theta,$sin_theta],[0,-$sin_theta,$cos_theta]);
	my $B = new Math::Matrix ([$cos_psi,$sin_psi,0],[-$sin_psi,$cos_psi,0],[0,0,1]);
	my $A = $B->multiply($C)->multiply($D);
	#
	my @new_element_rot  = ();
	my @new_coords_x_rot = ();
	my @new_coords_y_rot = ();
	my @new_coords_z_rot = ();
	#
	my @FragLines     = ();
	for (my $i = 0; $i < scalar (@atom); $i++) {
		push (@FragLines,"$atom[$i]  $crds_x[$i]  $crds_y[$i]  $crds_z[$i]");
	}
	#
	my @total_array      = ();
	while (my $Fline = shift (@FragLines)) {
		my @Cartesians              = split '\s+', $Fline;
		my ($Atom_label, @orig_xyz) = @Cartesians;
		my $matrix_xyz  = new Math::Matrix ([$orig_xyz[0],$orig_xyz[1],$orig_xyz[2]]);
		my $trans_xyz   = ($matrix_xyz->transpose);
		my $rotated_xyz = $A->multiply($trans_xyz);
		my @new_xyz     = split '\n+',$rotated_xyz;
		push (@new_element_rot ,$Atom_label);		
		push (@new_coords_x_rot,$new_xyz[0]);
		push (@new_coords_y_rot,$new_xyz[1]);		
		push (@new_coords_z_rot,$new_xyz[2]);
	}
	@total_array = ( [@new_element_rot], 
                     [@new_coords_x_rot],
                     [@new_coords_y_rot],					 
                     [@new_coords_z_rot]);
	return (@total_array);
}
#
sub Crossing_Over {
	my ($numberOfRec,$coords_opt) = @_;
	my @files      = @{$coords_opt};
	my $string_mol;
	my @array_xyz_coords = ();
	#
	if ( $numberOfRec > 1 ) {
		for (my $h = 0; $h < $numberOfRec; $h++) {
			for (my $y = 0; $y < $numberOfRec; $y++) {		
				if ( $h < $y ){
					my $pos1  = $h;
					my $pos2  = $y;
					#
					my $file1 = $files[$pos1];
					my $file2 = $files[$pos2];
					#
					my @array_F1 = split (/\n/,$file1);
					my @array_F2 = split (/\n/,$file2);		
					#
					my $boolean      = 0;
					my $count_repeat = 0;
					while ( $boolean < 1 ) {
						my @frag1     = ReadXYZFile (\@array_F1,1);
						my @frag2     = ReadXYZFile (\@array_F2,2);
						#
						my @centered  = CenterMolecule (\@frag1);
						my @centered2 = CenterMolecule (\@frag2);
						# Rotacion
						my @mol_rot_1 = cluster_rotation (@centered);
						my @mol_rot_2 = cluster_rotation (@centered2);
						#
						# Create random plane with a vector
						CreateRandomPlane();
						# Get coordinates and DotProduct for each 
						# molecule for the plane
						my ($hashWithCoordinates_1, $PP_coords1) = PlanePosition (1,@mol_rot_1);
						my ($hashWithCoordinates_2, $PP_coords2) = PlanePosition (2,@mol_rot_2);
						#
						# Sort molecules, 1 is decreasing order; 2 is increasing order.		
						my @sort1 = SortDotMolecule($PP_coords1, 1);
						my @sort2 = SortDotMolecule($PP_coords2, 2);
						#
						# Use data to combine 2 molecules. returns array with coordinates
						$string_mol    = CombineMolecule ($hashWithCoordinates_1, $hashWithCoordinates_2, \@sort1, \@sort2);
						my @arr_crds   = split (/\n/,$string_mol);		
						my $option_boo = steric_impediment_atoms (\@arr_crds,1);
						#
						if ( $option_boo == 1 ) {
							#print "$string_mol";
							$count_repeat++;
							if ( $count_repeat == 20 ) {
								#print "Crossing Over $count_repeat\n";
								#print "$string_mol\n";
								$boolean = 1;
							}
						}
						#
						if ( $option_boo == 0 ) { $boolean = 1;}
					}
					push (@array_xyz_coords,$string_mol);
				}
			}
		}
		#
		return @array_xyz_coords;
	} else {
		return @files;
	}
}
#
sub restart_species {
	my ($InputFile) = @_;
	#
	my @info_array = read_file ($InputFile);
	my $count      = 0;
	#
	my @first_arr  = ();
	foreach my $a_1 (@info_array){
		#
		$a_1 =~ s/^\s+|\s+$//g;
		$info_array[0] =~ s/^\s+|\s+$//g;
		if ( $a_1 eq $info_array[0] ){
			my $numb = $count + 2;
			push (@first_arr,$numb);
		}
		$count++;
	}
	#
	my @total_coords = ();
	my @total_id     = ();
	my $uniq_count   = 0;
	foreach my $a_1 (@first_arr){
		my $string;	
		my $sum = ($a_1 + $info_array[0] - 1);  
	    foreach my $i ($a_1..$sum) {
			$string.= "$info_array[$i]\n";
		}
		my $id               = sprintf '%.6d', $uniq_count;
		my $filebase         = "Adapt$id";
		push (@total_id,$filebase);
		push (@total_coords,$string);
		$uniq_count++;
	}
	#
	my $num_species = scalar (@total_id);
	return (\@total_id,\@total_coords,$num_species);	
}
#
sub final_opt_coords {
	my ($energy_opt,$coords_opt,$file_opt,$num_element,$option_unit,$name_file) = @_;
	# sort, same thing in reversed order
	my @array_energy_zero = @{$energy_opt}; 
	my @array_coords      = @{$coords_opt};
	my @array_files       = @{$file_opt};	
	#
	my @value_energy_sort = ();
	my @value_coords_sort = ();
	my @value_files_sort  = ();	
	my @idx = sort { $array_energy_zero[$a] <=> $array_energy_zero[$b] } 0 .. $#array_energy_zero;
	@value_energy_sort = @array_energy_zero[@idx];
	@value_coords_sort = @array_coords[@idx];
	@value_files_sort  = @array_files[@idx];	
	#
	my $filename = $name_file;
	open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
	for (my $i=0; $i < scalar (@value_energy_sort); $i++){
		my $resta = abs($value_energy_sort[0]) - abs($value_energy_sort[$i]);	
        # 1 Hartree = 27,2114 ev
        # 1 Hartree = 627,509 Kcal/mol
		my $eV      = 0;
		my $Kcalmol = 0;
		my $Hartree = 0;
		if ($option_unit == 0) {
			$eV      = sprintf("%.4f",(27.2114 * $resta ));
			$Kcalmol = sprintf("%.4f",(627.509 * $resta ));
			$Hartree = sprintf("%.4f",$value_energy_sort[$i]);
		} elsif ( $option_unit == 1 ) {
			$eV      = sprintf("%.4f",$value_energy_sort[$i]);
			$Kcalmol = sprintf("%.4f",( 23.060542   * $resta ));
			$Hartree = sprintf("%.4f",( 0.036749305 * $resta ));		
		} else {
			$eV      = sprintf("%.4f", 0.04336412 * $resta );
			$Kcalmol = sprintf("%.4f", $resta);
			$Hartree = sprintf("%.4f",( 0.001593601 * $resta ));
		}
        print $fh "$num_element\n";
        print $fh "$Kcalmol Kcal/mol $eV eV $Hartree H $value_files_sort[$i]\n";
		#print "$Hartree\t$value_files_sort[$i]\n";
        print $fh "$value_coords_sort[$i]";
	}
	close $fh;
}
#
sub coords_energy_file {
	my ($info_array) = @_;
	my @Secuencias   = @{$info_array};
	# coodenadas
	my @coords;
	# numero de atomos
	my $atom_numb = 0;
	# energia
	my $energy      = 0;
	#
	my @columns_1N  = ();
	my @columns_2N  = ();
	my @columns_3N  = ();
	my $count_lines = 0;
	#
	foreach my $a_1 (@Secuencias){
		# SCF Done:  E(RPBE1PBE) =  -56.7829127857     A.U. after   40 cycles
		if ( ($a_1=~/SCF/gi ) && ($a_1=~/Done/gi ) && ($a_1=~/after/gi ) ){
			my @array_tabs = ();
			@array_tabs = split (/ /,$a_1);
			push (@columns_1N  ,$array_tabs[7]);
		}
		# Standard orientation:
		if ( ($a_1=~/Standard/gi ) && ($a_1=~/orientation/gi ) && ($a_1=~/:/gi ) ){
			push (@columns_2N  ,$count_lines);
		}
		# Rotational constants (GHZ):
		if ( ($a_1=~/Rotational/gi ) && ($a_1=~/constants/gi ) && ($a_1=~/GHZ/gi ) ){
			push (@columns_3N  ,$count_lines);
		}
		$count_lines++;
	}
	#
	if ( scalar (@columns_1N) > 0 ){
		for (my $i=0; $i < scalar (@columns_1N); $i++){
			my $start = $columns_2N[$i] + 5;
			my $end   = $columns_3N[$i] - 2;
			$atom_numb = $end - $start + 1;
			$energy = $columns_1N[$i];
			@coords = ();
			foreach my $j (@Secuencias[$start..$end]){
				push (@coords,$j);		
			}
		}
		my $total_coords;
		foreach my $i (@coords){
			my @tmp = ();
			@tmp =  split (/\s+/,$i);
			$total_coords.= "$tmp[2]  $tmp[4]  $tmp[5]  $tmp[6]\n" 
			#push (@total_coords,"$tmp[2]\t$tmp[4]\t$tmp[5]\t$tmp[6]");
		}
		#
		return ($energy,$atom_numb,$total_coords);
	}
}
#
sub frequencies_coords {
	my ($info_array,$num_atoms) = @_;
	my @Secuencias   = @{$info_array};
	# coodenadas
	my @coords;
	my @columns_1N  = ();
	my $count_lines = 0;
	#
	foreach my $a_1 (@Secuencias){
		# Atom  AN      X      Y      Z	
		if ( ($a_1=~/Atom/gi ) && ($a_1=~/AN/gi ) && ($a_1=~/X/gi ) && ($a_1=~/Y/gi ) && ($a_1=~/Z/gi ) ){
			push (@columns_1N  ,$count_lines);
		}
		$count_lines++;
	}
	my $start = $columns_1N[0];
	my $end   = $columns_1N[1];
	foreach my $j (@Secuencias[$start..$end]){
		push (@coords,$j);		
	}
	# coords frecuencies
	my $string;
	for (my $i = 1; $i <= $num_atoms ; $i++) {
		my @array_tabs  = split (/\s+/,$coords[$i]);
		$string.= "$array_tabs[3]  $array_tabs[4]  $array_tabs[5]\n"
	}
	return $string;
}
#
sub frequencies_sum {
	my ($info_coords,$info_freq_coords,$atom_numb) = @_;
	my @array_coords = split (/\n/,$info_coords);
	my @array_freq   = split (/\n/,$info_freq_coords);
	#
	my $string;
	for (my $i = 0; $i < $atom_numb ; $i++) {
		my ($element,$cx,$cy,$cz) = split (/\s+/,$array_coords[$i]);
		my ($fx,$fy,$fz)          = split (/\s+/,$array_freq[$i]);
		my $sum_x = ($fx + $cx );
		my $sum_y = ($fy + $cy );
		my $sum_z = ($fz + $cz );
		$string.="$element  $sum_x  $sum_y  $sum_z\n";
	}
	return $string;
}
#
sub verification_of_termination {
	my ($info_array) = @_;
	my @Secuencias   = @{$info_array};
	#
	my $boolean = 3;
	#
	foreach my $a_1 (@Secuencias){
		# Normal termination
		if ( ($a_1=~/Normal/gi ) && ($a_1=~/termination/gi ) ){
			$boolean = 1;
		}
		# The combination of multiplicity or charge is impossible
		if (($a_1=~/combination/gi ) &&($a_1=~/multiplicity/gi ) && ($a_1=~/impossible/gi ) ){
			print ("ERROR  The Combination of Multiplicity or Charge is Impossible\n");
			exit(1);
		}
		# Error termination
		if ( ($a_1=~/Error/gi ) && ($a_1=~/termination/gi ) ){
			$boolean = 0;
		}
		# traceback not available
		if ( ($a_1=~/traceback/gi ) ){
			$boolean = 0;
		}
		# open-new-file
		if ( ($a_1=~/open-new-file/gi ) ){
			$boolean = 2;
		}		
		# Error: segmentation violation
	#	if ( ($a_1=~/Error:/gi ) && ($a_1=~/segmentation/gi ) && ($a_1=~/violation/gi ) ){
	#		$boolean = 2;
	#	}
	}
	return $boolean;
}
#
sub frequencies_NImag {
	my ($info_array) = @_;
	my @Secuencias   = @{$info_array};
	# energia
	my $energy      = 0;
	my $energy_zpe  = 0;
	my $numb_NImag  = 0;
	#
	my @columns_1N  = ();
	my @columns_2N  = ();
	my @columns_3N  = ();
	my @columns_4N  = ();	
	my $count_lines = 0;
	#
	foreach my $a_1 (@Secuencias){	
		# #
		# Item Value Threshold Converged
		if ( ($a_1=~/Item/gi ) && ($a_1=~/Value/gi ) && ($a_1=~/Threshold/gi ) && ($a_1=~/Converged/gi ) ){
			push (@columns_1N  ,$count_lines);
		}
		# Normal termination
		if ( ($a_1=~/Normal/gi ) && ($a_1=~/termination/gi ) ){
			push (@columns_2N  ,$count_lines);
		}
		# #
        # SCF Done:  E(RPBE1PBE) =  -56.7829127857     A.U. after   40 cycles
        if ( ($a_1=~/SCF/gi ) && ($a_1=~/Done/gi ) && ($a_1=~/after/gi ) ){
            my @array_tabs = ();
            @array_tabs = split (/ /,$a_1);
            push (@columns_3N  ,$array_tabs[7]);
        }
		# Zero-point correction= 
        if ( ($a_1=~/Zero-point/gi ) && ($a_1=~/correction/gi ) ){
            my @array_tabs = ();
            @array_tabs = split '\s+',$a_1;
            push (@columns_4N  ,$array_tabs[3]);
        }		
		$count_lines++;
	}
	my $start = $columns_1N[$#columns_1N];
	my $end   = $columns_2N[$#columns_2N];
	my $string;
	foreach my $j (@Secuencias[$start..$end]){
		chomp ($j);
		$string.=$j;		
	}
	# remove all whitespace
	$string =~ s/\s+//g;
	my @NImag_array = ();
	foreach my $pw (split(/\\/,$string)) {
		if ( ($pw=~/NImag/gi ) && ($pw=~/=/gi ) ){
			@NImag_array = split(/=/,$pw);
		}
	}
	#
	$energy     = $columns_3N[$#columns_3N];
	$energy_zpe = $columns_4N[$#columns_4N];
	$numb_NImag = $NImag_array[$#NImag_array]; 
	return ($energy,$energy_zpe,$numb_NImag);
}
#
sub Normal_check_Freq {
	my ($info_array) = @_;
	my @Secuencias   = @{$info_array};
	#
	my $option = 0;
	#
	my $count_lines  = 0;
	foreach my $a_1 (@Secuencias){	
		# Zero-point correction= 
        if ( ($a_1=~/Zero-point/gi ) && ($a_1=~/correction/gi ) ){
			$count_lines++;
         }		
	}
	if ( $count_lines == 1 ) {
		$option = 1;
	}
	#
	return $option;
}
#
sub distribute {
    my ($n, $array) = @_;

    my @parts;
    my $i = 0;
    foreach my $elem (@$array) {
        push @{ $parts[$i++ % $n] }, $elem;
    };
    return \@parts;
}
#
sub parallel_cpu_local {
	my ($info_file,$file_one,$count,$option,$software_option)  = @_;
	my @tmp_arr = @{$info_file};
	#
	my @Outputs_lammps = ();
	#
	my $slrm = "Lanz-tmp_$count.sh";
	open (SLURMFILE, ">$slrm");
	#
	print SLURMFILE "#!/bin/sh \n";
	print SLURMFILE "\n";
	#
	if ( $option == 0 ) {
		foreach my $i (@tmp_arr) {
			(my $without_ext = $i) =~ s/\.[^.]+$//;
            if ($software_option == 0) { print SLURMFILE "$gaussian_version $without_ext.com\n"; }
            if ($software_option == 1) { print SLURMFILE "$path_bin_mopac $without_ext.mop >tmp_mopac_1.txt 2>tmp_mopac_2.txt\n"; }
			if ($software_option == 2) { print SLURMFILE "$path_bin_lammps -in $without_ext.in >$without_ext.out\n"; push (@Outputs_lammps,"$without_ext.out"); }
		}
	} else {
		if ($software_option == 0) { print SLURMFILE "$gaussian_version $file_one.com\n"; }
        if ($software_option == 1) { print SLURMFILE "$path_bin_mopac $file_one.mop >tmp_mopac_1.txt 2>tmp_mopac_2.txt\n"; }
        if ($software_option == 2) { print SLURMFILE "$path_bin_lammps -in $file_one.in >$file_one.out\n"; }
	}    
	close (SLURMFILE);
	#
	system ("sh $slrm &");
	sleep (3);
}
#
sub slurm_cluster {
	#
	my ($info_file,$file_one,$count,$option)  = @_;
	my @tmp_arr                               = @{$info_file};
	#
	my $slrm = "Lanz-tmp_$count.slrm";
	open (SLURMFILE, ">$slrm");
	#
	print SLURMFILE "#!/bin/bash \n";
	print SLURMFILE "#SBATCH --job-name=AutomataCelular\n";
	print SLURMFILE "#SBATCH --output=AutomataCelular.out\n";
	print SLURMFILE "#SBATCH --error=error_AutomataCelular\n";
	my $days     = 7;
	my $word     = $days."d";
	my $new_days = $word;
	print SLURMFILE "#SBATCH --partition=q$new_days-20p\n";
	print SLURMFILE "#SBATCH --nodes=1\n";
	print SLURMFILE "#SBATCH --ntasks=20\n";
	print SLURMFILE "#SBATCH -t 0$days-00:00:00\n";
	print SLURMFILE "##SBATCH --exclude=nc[15,17,18]\n";
	print SLURMFILE "\n";
	print SLURMFILE "module load gaussian/09-sse4\n";
	print SLURMFILE "\n";
	print SLURMFILE "cd \$SLURM_SUBMIT_DIR\n";
	print SLURMFILE "\n";
	print SLURMFILE "scontrol show hostname \$SLURM_NODELIST > hostlist.dat\n";
	print SLURMFILE "\n";
	print SLURMFILE "export GAUSS_SCRDIR=/tmp\n";
	print SLURMFILE "\n";
	print SLURMFILE "# #\n";
	print SLURMFILE "#\n";
	if ( $option == 0 ) {
		foreach my $i (@tmp_arr) {
			(my $without_ext = $i) =~ s/\.[^.]+$//;
			print SLURMFILE "$gaussian_version $without_ext.com\n";
		}
	} else {
		print SLURMFILE "$gaussian_version $file_one.com\n";	
	}	
	print SLURMFILE "# #\n\n";
	print SLURMFILE "rm /tmp/* \n";
	#
	close (SLURMFILE);
	#
	system ("sbatch $slrm");
}
###################################
# submit Mopac && Lammps
sub submit_Mopac_Lammps {
	my ($arrayInputs,$path_mopac,$software_option) = @_;
	#
	my @files  = @{$arrayInputs};
    #
	my $number_paral     = ($numero_colas - 1);
	#
	my %Info_files_dir   = ();
	my @array_keys       = ();
	#
	my @total_files_qm   = ();
	# # # # 
	# Local node
	#
	my @arrfile = ();
	my $arr     = distribute($number_paral,\@files);
	my @loarr   = @{$arr};
	my $div     = int ( scalar (@files) / $number_paral);		
	for (my $i = 0; $i < scalar (@loarr) ; $i++) {
		my @tmp_arr = ();
		for (my $j = 0; $j < $div ; $j++) {
			push (@tmp_arr,$loarr[$i][$j]);
			push (@arrfile,$loarr[$i][$j]);
		}
		sleep(3);		
		parallel_cpu_local (\@tmp_arr,"NOFILE",$i,0,$software_option);
	}
	my @delete_arr    = ();
	my @element_files = ();
	foreach my $i (@files) {
		foreach my $j (@arrfile) {
			if ( ($i=~/$j/gi ) ){
				#print "$j\n";
				push (@element_files,$j);
				@delete_arr = ();
			} else {
				push (@delete_arr,$i);
			}
		}
	}
	my @filtered   = uniq(@delete_arr);
	my $data_value = scalar (@loarr);
    	#
	sleep(3);    
	parallel_cpu_local (\@filtered,"NOFILE",$data_value,0,$software_option);
	foreach my $i (@filtered) {
		push (@element_files,$i);
	}
	if ( ( scalar (@files) ) != ( scalar (@element_files) ) ) { print "Problem Number Files Local\n"; exit;}
	# For Lammps
	my @Outputs_lammps = ();
	for ( my $i=0 ; $i < scalar(@files); $i++) {
		(my $without_ext = $files[$i]) =~ s/\.[^.]+$//;
		push (@Outputs_lammps,"$without_ext.out");
	}
	return @Outputs_lammps	
}
###################################
# submit gaussian
sub submit_queue_gaussian { 
	my ($arrayInputs,$number_ele,$header,$ncpus,$mem,$Charge,$Multiplicity,$queue,$exec_bin_g09,$software_option) = @_;
	#
	my @files  = @{$arrayInputs}; 
	# numero de atomos 
	my $atom_numb        = $number_ele;
	my $number_paral     = ($numero_colas - 1);
	# Corrupted File
	my @name_files_corr       = ();
	my @report_structure_corr = ();
	# Convergence Failure
	my @name_files       = ();
	my @report_structure = ();
	#
	my @report_NImag     = ();
	#
	my @NImag_files      = ();
	#
	my %Info_files_dir   = ();
	my @array_keys       = ();
	#
	my @total_files_qm   = ();
	#
	my $Port_Corrupt     = 0;
	#
	for (my $i = 0; $i < scalar(@files) ; $i++) {
		#
		(my $without_ext = $files[$i]) =~ s/\.[^.]+$//;
		#
		my $id               = sprintf("%06d",$i);
		$Info_files_dir{$id} = $files[$i];
		push (@array_keys,$id);
		push (@total_files_qm,$files[$i]);
		#
		if ($local_cluster == 1) {
			my $env_a = `$exec_bin_g09 $without_ext $without_ext.com $ncpus $queue`;
		}
	}
	# # # # 
	# Cluster Slurm
	#
	if ($local_cluster == 0) {
		my @arrfile = ();
		my $arr     = distribute($number_paral,\@files);
		my @loarr   = @{$arr};
		my $div     = int ( scalar (@files) / $number_paral);		
		for (my $i = 0; $i < scalar (@loarr) ; $i++) {
			my @tmp_arr = ();
			for (my $j = 0; $j < $div ; $j++) {
				push (@tmp_arr,$loarr[$i][$j]);
				push (@arrfile,$loarr[$i][$j]);
			}
			#slurm_cluster (\@tmp_arr,"NOFILE",$i,0);
			parallel_cpu_local (\@tmp_arr,"NOFILE",$i,0,$software_option);
		}
		my @delete_arr    = ();
		my @element_files = ();
		foreach my $i (@files) {
			foreach my $j (@arrfile) {
				if ( ($i=~/$j/gi ) ){
					#print "$j\n";
					push (@element_files,$j);
					@delete_arr = ();
				} else {
					push (@delete_arr,$i);
				}
			}
		}
		my @filtered   = uniq(@delete_arr);
		my $data_value = scalar (@loarr);
		#slurm_cluster (\@filtered,"NOFILE",$data_value,0);
		parallel_cpu_local (\@filtered,"NOFILE",$data_value,0,$software_option);
		foreach my $i (@filtered) {
			push (@element_files,$i);
		}
		if ( ( scalar (@files) ) != ( scalar (@element_files) ) ) { print "Problem Number Files Slurm, PBS & Local\n"; exit;}
	}
	#
	#
	#
	# # # # # # # # # # #
	my $count = 0;
	while ( $count < 1 ) {
		while (my ($key, $value) = each %Info_files_dir) {
			# Corrupted File
			my %Info_count_files_Corr = ();
			# Convergence Failure
			my %Info_count_files = ();
			#	
			my %Info_NImag_files = ();
			#
			my $input_file = $value;
			############# 		
			# Main
			(my $without_extension = $input_file) =~ s/\.[^.]+$//;
            
			if( ( -e "$without_extension.out" ) || ( -e "$without_extension.log" ) ) {
				my @Secuencias   = ();
				if( -e "$without_extension.out" ) {
					sleep(3);
					@Secuencias      = read_file ("$without_extension.out");
				}
				if( -e "$without_extension.log" ) {
					sleep(3);
					@Secuencias      = read_file ("$without_extension.log");
				}
				#
				my $option       = verification_of_termination(\@Secuencias);
							
				# # # # # # # # # # # # #
				# Convergence Failure
				if ( $option == 0 ) {
					#print "Convergence Failure -> File: $input_file\n";	
					my $energy = 0;
					my $total_coords;
					my $deci   = 0;
					#
					push (@name_files,$without_extension); 
					foreach my $element( @name_files ) {
						++$Info_count_files{$element};
					}
					if ( $Info_count_files{$without_extension} > $max_numb_conver ) {
						push (@report_structure,$Info_files_dir{$key});
						delete $Info_files_dir{$key};
						$deci = 1;
					}
					#
					if ( $deci == 0 ) {
						# funcion
						if( -e "$without_extension.com" ){ unlink ("$without_extension.com");}
						($energy,$atom_numb,$total_coords) = coords_energy_file (\@Secuencias);
						#
						# If the coordinates are empty
						if (!defined($total_coords)) {
							#
							#print "No Coordinates\n";
							#
							delete $Info_files_dir{$key};					
						} else {
							my $G09name = G03Input ($without_extension,$header,$ncpus,$mem,$Charge,$Multiplicity,$total_coords,$energy);
							if( -e "$without_extension.chk" ){ unlink ("$without_extension.chk");}
							if( -e "$without_extension.out" ){ unlink ("$without_extension.out");}
							if( -e "$without_extension.log" ){ unlink ("$without_extension.log");}
							#
							my @empty_arr = ();
							#
							if ($local_cluster == 0) {
								sleep(2);
								#slurm_cluster (\@empty_arr,$without_extension,100,1);
								parallel_cpu_local (\@empty_arr,$without_extension,100,1,$software_option);
							} else {
								sleep(2);
								my $env_b = `$exec_bin_g09 $without_extension $without_extension.com $ncpus $queue`;
							}
						}
					}				
				}
				# # # # # # # # # # # # #
				# Corrupted File
				if ( $option == 2 ) {
					#print "Corrupted File -> File: $input_file\n";	
					my $deci   = 0;
					#
					push (@name_files_corr,$without_extension); 
					foreach my $element( @name_files_corr ) {
						++$Info_count_files_Corr{$element};
					}
					if ( $Info_count_files_Corr{$without_extension} > $max_numb_conver ) {
						push (@report_structure_corr,$Info_files_dir{$key});
						delete $Info_files_dir{$key};
						$deci = 1;
					}
					#
					if ( $deci == 0 ) {
						# funcion
						#if( -e "$without_extension.com" ){ unlink ("$without_extension.com");}
						if( -e "$without_extension.chk" ){ unlink ("$without_extension.chk");}
						if( -e "$without_extension.out" ){ unlink ("$without_extension.out");}
						if( -e "$without_extension.log" ){ unlink ("$without_extension.log");}
						#
						my @empty_arr = ();
						#
						if ($local_cluster == 0) {
							sleep(2);
							#slurm_cluster (\@empty_arr,$without_extension,100,1);
							parallel_cpu_local (\@empty_arr,$without_extension,100,1,$software_option);
						} else {
							sleep(2);
							my $env_b = `$exec_bin_g09 $without_extension $without_extension.com $ncpus $queue`;
						}
					}
					$Port_Corrupt++;
				}	
				# # # # # # # # # # # # #
				# Normal termination
				my $energy = 0;
				my $total_coords;
				#print "$Info_files_dir{$key}\n";
				my $string = "Freq";
				if (index($header, $string) != -1) {			
					#
					# Verificar si existen frecuencias
					my $verif_normal = Normal_check_Freq (\@Secuencias);
					if ( $verif_normal == 1 ) { # XXXX
						my $energy      = 0;
						my $energy_zpe  = 0;
						my $numb_NImag  = 0;
						#
						($energy,$energy_zpe,$numb_NImag) = frequencies_NImag (\@Secuencias);
						if ( $numb_NImag == 0 ) {
							#print "Normal Termination -> File : $input_file\n";
							delete $Info_files_dir{$key};
						} else {
							#print "Imaginary Frequencies -> File : $input_file\n";
							my $dec_nimag = 0;
							push (@NImag_files,$without_extension); 
							foreach my $element( @NImag_files ) {
								++$Info_NImag_files{$element};
							}
							if ( $Info_NImag_files{$without_extension} > $max_numb_NImag ) {
								push (@report_NImag,"$Info_files_dir{$key} = $numb_NImag");
								delete $Info_files_dir{$key};
								$dec_nimag = 1;
							}
							#
							if ( $dec_nimag == 0 ) {
								# funcion
								($energy,$atom_numb,$total_coords) = coords_energy_file (\@Secuencias);
								my $all_freq_coords = frequencies_coords (\@Secuencias,$atom_numb);
								# Aqui le sumamos las frecuencias del susodicho							
								my $final_freq      = frequencies_sum ($total_coords,$all_freq_coords,$atom_numb);
								if( -e "$without_extension.com" ){ unlink ("$without_extension.com");}
								my $G09name = G03Input ($without_extension,$header,$ncpus,$mem,$Charge,$Multiplicity,$final_freq,$energy);
								if( -e "$without_extension.chk" ){ unlink ("$without_extension.chk");}
								if( -e "$without_extension.out" ){ unlink ("$without_extension.out");}
								if( -e "$without_extension.log" ){ unlink ("$without_extension.log");}
								my @empty_arr = ();
								#	
								if ($local_cluster == 0) {
									sleep(2);
									#slurm_cluster (\@empty_arr,$without_extension,100,1);
									parallel_cpu_local (\@empty_arr,$without_extension,100,1,$software_option);
								} else {
									sleep(2);
									my $env_c = `$exec_bin_g09 $without_extension $without_extension.com $ncpus $queue`;
								}	
								#
							}
						}
					} # XXXX
				} else {
					if ( $option == 1 ) { #-----
						#print "Normal Termination -> File : $input_file\n";
						delete $Info_files_dir{$key};
					}	#-----	
				}
			}
		}
		if (!%Info_files_dir) {
			#print "Normal Termination CELL\n";
			$count = 1;
		}
	}
	#
	my @coords_opt       = ();
	my @energy_opt       = ();	
	my @file_opt         = ();
	my @num_element      = ();
	#
	my $Port_Normal    = 0;
	my $Port_Error     = 0;
	my $Port_No_Coords = 0;
	#
	@files_Port_Error = ();
	# No coords
	@files_No_Coords  = ();	
	#
	foreach my $file (@total_files_qm) {
		(my $without_extension = $file) =~ s/\.[^.]+$//;
        #
		if( ( -e "$without_extension.out" ) || ( -e "$without_extension.log" ) ) {
			my @Secuencias   = ();
			if( -e "$without_extension.out" ) {
				sleep(3);
				@Secuencias      = read_file ("$without_extension.out");
			}
			if( -e "$without_extension.log" ) {
				sleep(3);
				@Secuencias      = read_file ("$without_extension.log");
			}
			my $option       = verification_of_termination(\@Secuencias);
			# # # # # # # # # # # # #
			# Normaltermination
			if ( $option == 1 ) {
				my $string = "Freq";
				if (index($header, $string) != -1) {
					my $energy      = 0;
					my $energy_zpe  = 0;
					my $numb_NImag  = 0;
					#
					($energy,$energy_zpe,$numb_NImag)     = frequencies_NImag (\@Secuencias);
					my ($energy_o,$atom_numb,$coords_opt) = coords_energy_file (\@Secuencias);
					my $sum_zpe = ($energy + $energy_zpe);
					push (@energy_opt,$sum_zpe);
					push (@coords_opt,$coords_opt);
					push (@file_opt,$without_extension);
					push (@num_element,$atom_numb);			
				} else {
					my ($energy_scf,$atom_numb,$coords_opt) = coords_energy_file (\@Secuencias);
					push (@energy_opt ,$energy_scf);
					push (@coords_opt ,$coords_opt);
					push (@file_opt   ,$without_extension);
					push (@num_element,$atom_numb);			
				}
			} else {
				my ($energy_scf,$atom_numb,$coords_opt) = coords_energy_file (\@Secuencias);
				#
				if (!defined($coords_opt)) {	
					push (@files_No_Coords,$without_extension);	
					$Port_No_Coords++;				
				} else {
					push (@energy_opt ,$energy_scf);
					push (@coords_opt ,$coords_opt);
					push (@file_opt   ,$without_extension);
					push (@num_element,$atom_numb);
					#
					push (@files_Port_Error,$without_extension);
				}
			}
			# Proportion of Normal and Error Termination
			if ( $option == 1 ) {
				$Port_Normal++;
			} else {
				$Port_Error++;
			}
		}
	}
	#
	$Global_Port_Normal  = $Port_Normal;
	$Global_Port_Error   = $Port_Error;
	$Global_Port_Corrupt = $Port_Corrupt;
	#
	$Global_No_Coords    = $Port_No_Coords;
	#
	my @value_energy_sort = ();
	my @value_coords_sort = ();
	my @value_files_sort  = ();
	my @idx               = sort { $energy_opt[$a] <=> $energy_opt[$b] } 0 .. $#energy_opt;
	@value_energy_sort    = @energy_opt[@idx];
	@value_coords_sort    = @coords_opt[@idx];
	@value_files_sort     = @file_opt[@idx];
	#
	return (\@value_energy_sort,\@value_coords_sort,\@report_structure,\@report_NImag,\@value_files_sort);
}
###################################
# delete repeat data
sub uniq {
	my %seen;
	grep !$seen{$_}++, @_;
}
###################################
# index duplicate data
sub index_elements {
	my ($duplicate_name,$files_name) = @_;
	# reference arrays	
	my @array_1     = @{$duplicate_name}; 
	my @array_2     = @{$files_name};
	my @array_index = ();
	#
	my @filtered = uniq(@array_1);
	foreach my $u (@filtered){
		my @del_indexes = reverse( grep { $array_2[$_] eq "$u" } 0..$#array_2);
		foreach my $k (@del_indexes) {
			push (@array_index,$k);
		}
	}
	return @array_index;
}
###################################
# promedio
sub promedio {
	my ($num,$data) = @_;
	# write file
	my $sum = 0;
	for ( my $i = 0 ; $i < $num ; $i = $i + 1 ){
		$sum+= @$data[$i];
	}
	my $div = $sum / $num;
	return $div; 
}
###################################
# verify similar structure
sub info_duplicate_structures {
	my ($number_cycle,$numb_atoms,$coords_xyz,$file_energy,$name_file,$ncpus,$threshold_duplicate) = @_;
	my @array_coords = @{$coords_xyz};
	my @array_energy = @{$file_energy};
	my @array_name   = @{$name_file};
	#
	my $add = $number_cycle;
	#
	my %Info_Coords = ();
	#
	my %Info_Axis   = ();
	my %Info_Energy = ();
	my %Info_Name   = ();
	#
	my @array_keys  = ();
	for (my $i=0; $i < $add ; $i++) { 
		my $id            = sprintf("%.5d",$i);
		#
		$Info_Energy{$id} = $array_energy[$i];
		$Info_Name{$id}   = $array_name[$i];
		$Info_Axis{$id}   = $array_coords[$i];
		#
		my @abc           = split (/\n/,$array_coords[$i]);
		$Info_Coords{$id} = \@abc;
		push(@array_keys,$id);
	}
	my $pm        = new Parallel::ForkManager($ncpus);
	my $iteration = 0;
	#
	my $file_tmp = "Dupli.tmp";
	open (FILE, ">$file_tmp") or die "Unable to open XYZ file: $file_tmp";
	for ( my $x = 0 ; $x < scalar (@array_keys); $x = $x + 1 ) {
		$pm->start($iteration) and next;
		# All children process havee their own random.			
		srand();
		for ( my $y = 0 ; $y < scalar (@array_keys); $y = $y + 1 ) {
			if ( $x < $y ){
				#
				my @matrix_1 = @{$Info_Coords{$array_keys[$x]}};
				my @matrix_2 = @{$Info_Coords{$array_keys[$y]}};
				# # # # # # # # # # # # # # # # #
				#
				my @array_name_atoms_1 = ();
				my @array_coord_x_1    = ();
				my @array_coord_y_1    = ();
				my @array_coord_z_1    = ();
				#
				my @array_name_atoms_2 = ();
				my @array_coord_x_2    = ();
				my @array_coord_y_2    = ();
				my @array_coord_z_2    = ();	
				#
				for ( my $i = 0 ; $i < $numb_atoms ; $i = $i + 1 ){
					my @array_tabs_1  = split (/\s+/,$matrix_1[$i]);
					#
					my $radii_val;
					my $other_element = 0;
					if ( exists $Atomic_number{$array_tabs_1[0]} ) {
						# exists
						$radii_val = $Atomic_number{$array_tabs_1[0]};
						$array_name_atoms_1[++$#array_name_atoms_1] = $radii_val;
					} else {
						# not exists
						$radii_val = $array_tabs_1[0] ;
						$array_name_atoms_1[++$#array_name_atoms_1] = $radii_val;
					}
					$array_coord_x_1[++$#array_coord_x_1]   = $array_tabs_1[1];
					$array_coord_y_1[++$#array_coord_y_1]   = $array_tabs_1[2];
					$array_coord_z_1[++$#array_coord_z_1]   = $array_tabs_1[3];
				}
				#
				for ( my $i = 0 ; $i < $numb_atoms ; $i = $i + 1 ){
					my @array_tabs_2 = split (/\s+/,$matrix_2[$i]);
					#
					my $radii_val;
					my $other_element = 0;
					if ( exists $Atomic_number{$array_tabs_2[0]} ) {
						# exists
						$radii_val = $Atomic_number{$array_tabs_2[0]};
						$array_name_atoms_2[++$#array_name_atoms_2] = $radii_val;
					} else {
						# not exists
						$radii_val = $array_tabs_2[0] ;
						$array_name_atoms_2[++$#array_name_atoms_2] = $radii_val;
					}
					$array_coord_x_2[++$#array_coord_x_2]   = $array_tabs_2[1];
					$array_coord_y_2[++$#array_coord_y_2]   = $array_tabs_2[2];
					$array_coord_z_2[++$#array_coord_z_2]   = $array_tabs_2[3];
				}
				my $Springborg = Grigoryan_Springborg ($numb_atoms,\@array_coord_x_1 ,\@array_coord_y_1 ,\@array_coord_z_1 
																,\@array_coord_x_2 ,\@array_coord_y_2 ,\@array_coord_z_2 );												  
				#
				if ( $Springborg < $threshold_duplicate ) {
					my $number      = sprintf '%.6f', $Springborg;
					print FILE "$array_keys[$y]\n";
					print FILE "Value = $number\n";
				}
				#
			}
			$iteration++;
		}
		$pm->finish;	
	}
	close (FILE);
	# Paralel
	$pm->wait_all_children;
	# # #
	my @data_tmp = read_file ($file_tmp);	
	my @duplicates_name = ();
	my @Value_simi      = ();
	foreach my $info (@data_tmp) {
		if ( ($info =~ m/Value/) ) {
			my @array_tabs = ();
			@array_tabs    = split ('\s+',$info);
			push (@Value_simi,$array_tabs[2]);
		} else {
			push (@duplicates_name,$info);
		}
	}
	my @array_similar_coords = ();
	my @array_similar_energy = ();
	my @array_similar_files  = ();
	# Delete similar structures
	my @index_files = index_elements (\@duplicates_name,\@array_keys);
	#
	# Delete similar structures	
	for my $k (@index_files) {
		push (@array_similar_coords,$Info_Axis{$array_keys[$k]});
		push (@array_similar_energy,$Info_Energy{$array_keys[$k]});
		push (@array_similar_files ,$Info_Name{$array_keys[$k]});
		delete $Info_Coords{$array_keys[$k]};
	}
	#
	unlink ($file_tmp);
	#
	my @keys_arr    = keys %Info_Coords;
	my @axis_all    = ();
	my @energy_all 	= ();
	my @name_all 	= ();
	for my $co (@keys_arr) {
		push (@axis_all  ,$Info_Axis{$co});
		push (@energy_all,$Info_Energy{$co});
		push (@name_all  ,$Info_Name{$co});		
	}
	#
	return (\@axis_all,\@energy_all,\@name_all,\@array_similar_coords,\@array_similar_energy,\@array_similar_files);
}
###################################
# Grigoryan Springborg similitud
sub Grigoryan_Springborg {
	my ($numb_atoms,$array_coord_x_1,$array_coord_y_1, $array_coord_z_1,
	                $array_coord_x_2,$array_coord_y_2, $array_coord_z_2) = @_;
	#
	my @distance_alpha = ();
	my @distance_beta  = ();
	#
	my $sum_1 = 0;
	for ( my $i = 0 ; $i < $numb_atoms ; $i = $i + 1 ){
		for ( my $j = 0 ; $j < $numb_atoms ; $j = $j + 1 ){
			if ( $i < $j ){
				my $distance = Euclidean_distance (@$array_coord_x_1[$i],@$array_coord_y_1[$i],@$array_coord_z_1[$i],
												@$array_coord_x_1[$j],@$array_coord_y_1[$j],@$array_coord_z_1[$j]);
				push (@distance_alpha,$distance);
			}
		}
	}
	#
	my $sum_2 = 0;
	for ( my $i = 0 ; $i < $numb_atoms ; $i = $i + 1 ){
		for ( my $j = 0 ; $j < $numb_atoms ; $j = $j + 1 ){
			if ( $i < $j ){
				my $distance = Euclidean_distance (@$array_coord_x_2[$i],@$array_coord_y_2[$i],@$array_coord_z_2[$i],
												@$array_coord_x_2[$j],@$array_coord_y_2[$j],@$array_coord_z_2[$j]);
				push (@distance_beta,$distance);
			}
		}
	}
	#
	my $InterDist_1 = (2/($numb_atoms*($numb_atoms-1)));
	my $InterDist_2 = (($numb_atoms*($numb_atoms-1))/2);  
	#
	my @mol_alpha = ();
	my @mol_beta  = ();
	my @idx_1 = sort { $distance_alpha[$a] <=> $distance_alpha[$b] } 0 .. $#distance_alpha;
	my @idx_2 = sort { $distance_beta[$a]  <=> $distance_beta[$b]  } 0 .. $#distance_beta;
	@mol_alpha = @distance_alpha[@idx_1];
	@mol_beta  = @distance_beta[@idx_2];
	#
	my $num_1 = scalar (@mol_alpha);
	my $num_2 = scalar (@mol_beta);
	my $dim_alpha =  promedio ($num_1,\@mol_alpha);
	my $dim_beta  =  promedio ($num_2,\@mol_beta);
	#
	my $sumX;
	my $sumY;
	# Sin normalizar
	for ( my $i = 0 ; $i < $InterDist_2 ; $i = $i + 1 ){
		my $mult = ( $mol_alpha[$i] - $mol_beta[$i] )**2;
		$sumX+=$mult;
	}
	my $Springborg_1 = sqrt( $InterDist_1 * $sumX );
	# Normalizado
	for ( my $i = 0 ; $i < $InterDist_2 ; $i = $i + 1 ){
		my $mult = ( ($mol_alpha[$i]/$dim_alpha) - ($mol_beta[$i]/$dim_beta) )**2;
		$sumY+=$mult;
	}
	my $Springborg_2 = sqrt( $InterDist_1 * $sumY );
	#
	return $Springborg_2; 
}
###################################
# read files
sub read_file {
	# filename
	my ($filename) = @_;
	(my $input_file = $filename) =~ s/\s//g;
	my @array          = ();
	# open file
	open(FILE, "<", $input_file ) || die "Can't open $input_file: $!";
	while (my $row = <FILE>) {
		chomp($row);
		push (@array,$row);
	}
	close (FILE);
	# return array	
	return @array;
}
###################################
# drawing a box around a molecule 
sub drawing_box_molecule_VMD {
	my ($coordsmin, $coordsmax) = @_;
	#
	my $minx = @$coordsmin[0];
	my $maxx = @$coordsmax[0];
	my $miny = @$coordsmin[1];
	my $maxy = @$coordsmax[1];
	my $minz = @$coordsmin[2];
	my $maxz = @$coordsmax[2];
	# raw the lines
	
	my $filename = '04BOX_kick.vmd';
	open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";	
	print $fh "draw delete all\n";
	print $fh "draw materials off\n";
	print $fh "draw color green\n";
	#
	print $fh "draw line \"$minx $miny $minz\" \"$maxx $miny $minz\" \n";
	print $fh "draw line \"$minx $miny $minz\" \"$minx $maxy $minz\" \n";
	print $fh "draw line \"$minx $miny $minz\" \"$minx $miny $maxz\" \n";
	#
	print $fh "draw line \"$maxx $miny $minz\" \"$maxx $maxy $minz\" \n";
	print $fh "draw line \"$maxx $miny $minz\" \"$maxx $miny $maxz\" \n";
	#
	print $fh "draw line \"$minx $maxy $minz\" \"$maxx $maxy $minz\" \n";
	print $fh "draw line \"$minx $maxy $minz\" \"$minx $maxy $maxz\" \n";
	#
	print $fh "draw line \"$minx $miny $maxz\" \"$maxx $miny $maxz\" \n";
	print $fh "draw line \"$minx $miny $maxz\" \"$minx $maxy $maxz\" \n";
	#
	print $fh "draw line \"$maxx $maxy $maxz\" \"$maxx $maxy $minz\" \n";
	print $fh "draw line \"$maxx $maxy $maxz\" \"$minx $maxy $maxz\" \n";
	print $fh "draw line \"$maxx $maxy $maxz\" \"$maxx $miny $maxz\" \n";
	close $fh;
}					  
###################################
# input file gaussian 
sub G03Input {
	#
	my ($filebase,$Header,$ncpus,$mem,$Charge,$Multiplicity,$coordsMM,$iteration) = @_;
	my $G03Input     = "$filebase.com";
	#
	if (!defined($coordsMM)) {
		print "$filebase = Without Coordinates\n";
	}
	#
	open (COMFILE, ">$G03Input");
	#print COMFILE "%chk=$filebase.chk\n";
	if ( $ncpus > 0 ) {
		print COMFILE "%NProc=$ncpus\n";
	}	
	(my $word_nospaces = $mem) =~ s/\s//g;
	print COMFILE "%mem=$word_nospaces"."GB\n";
	print COMFILE "# $Header \n";
	print COMFILE "\nCellular Automata job $iteration\n";
	print COMFILE "\n";
	print COMFILE "$Charge $Multiplicity\n";
	print COMFILE "$coordsMM";
	print COMFILE "\n";
	close COMFILE;
	#
	return $G03Input;
}
###################################
# submit mopac
sub submit_queue_mopac {
	my ($arrayInputs,$path_mopac) = @_;
	my @Inputs_mopac = @{$arrayInputs};
	for ( my $i=0 ; $i < scalar(@Inputs_mopac); $i++) {
		my $env = `$path_mopac $Inputs_mopac[$i] >tmp_mopac_1.txt 2>tmp_mopac_2.txt`;	
	}
	unlink ("tmp_mopac_1.txt");
	unlink ("tmp_mopac_2.txt");
}
###################################
# input file Mopac
sub MopacInput {
	#
	my $filebase     = $_[0];
	my $coordsMM     = $_[1];
	my $MopacInput   = "$filebase.mop";
	my $iteration    = $_[2];
	my $Headerfile   = $_[3];
	my $Charge       = $_[4];
	my $Multiplicity = $_[5];
	#
	my $mem          = $_[7];
	#
	my $tmp   = 1;
	my @words = split (/\n/,$coordsMM);
	#
	open (COMFILE, ">$MopacInput");
	#
	my $word;
	# Spin multiplicity:
	if ( $Multiplicity == 0 ) { $word = "NONET"   };			
	# singlet	- 0 unpaired electrons
	if ( $Multiplicity == 1 ) { $word = "SINGLET" };
	# doublet	- 1 unpaired electrons
	if ( $Multiplicity == 2 ) { $word = "DOUBLET" };
	# triplet	- 2 unpaired electrons
	if ( $Multiplicity == 3 ) { $word = "TRIPLET" };
	# quartet	- 3 unpaired electrons
	if ( $Multiplicity == 4 ) { $word = "QUARTET" };
	# quintet	- 4 unpaired electrons			
	if ( $Multiplicity == 5 ) { $word = "QUINTET" };
	# sextet	- 5 unpaired electrons
	if ( $Multiplicity == 6 ) { $word = "SEXTET"  };
	# septet	- 6 unpaired electrons
	if ( $Multiplicity == 7 ) { $word = "SEPTET"  };
	# octet	- 7 unpaired electrons
	if ( $Multiplicity == 8 ) { $word = "OCTET"   };
	#
	my $ncpus        = ($_[6] * 2);
	if ( $ncpus == 0 ) {
		print COMFILE "$Headerfile $word CHARGE=$Charge";
	} else {
		# The maximum number of threads is normally equal to the number of cores, 
		# even if each core supports two threads.
		# In the special case of THREADS=1, parallelization is switched off.
		print COMFILE "$Headerfile $word CHARGE=$Charge THREADS=$ncpus";
	}	
	print COMFILE "\n";
	print COMFILE "Cell job $iteration\n";
	print COMFILE "\n";	
	foreach my $i (@words){
		my @axis    = split (" ",$i);
		#
		my $label  = $axis[0];  
		my $axis_x = $axis[1];
		my $axis_y = $axis[2];
		my $axis_z = $axis[3];
		#
		print COMFILE "$label\t$axis_x\t$tmp\t$axis_y\t$tmp\t$axis_z\t$tmp\n";
	}
	print COMFILE "\n";
	print COMFILE "\n";
	close (COMFILE);
	#
	return $MopacInput;
}
###################################
# Energy ouputsmopac
sub energy_mopac {
	my ($num_atoms_xyz,$arrayInputs) = @_;	
	#
	my @array              = ();
	my @array_coords_mopac = ();	
	foreach my $files ( @{$arrayInputs} ) {
		(my $without_extension = $files) =~ s/\.[^.]+$//;	
		if ( -e "$without_extension.arc" ) {
			push( @array, "$without_extension.arc");
		} else {
			print "WARNING: Mopac file $files error termination\n";		
		}
	}
	#
	my $tam_esc = scalar (@array);
	if ($tam_esc == 0) { print "ERROR: Problem MOPAC $tam_esc files .arc, check .out\n"; exit(0);}
	#
	my @HeaderLines   = ();
	my @ZeroPoint     = ();
	my @energyy       = ();
	my @energy_hartre = ();
	my $energy        = '';
	my $number_atoms  = $num_atoms_xyz;
	foreach my $i (@array) {
		(my $without_extension = $i) =~ s/\.[^.]+$//;
		open(HEADER,"$without_extension.arc") or die "Unable to open $i";
		@HeaderLines  = <HEADER>;
		close HEADER;
		while (my $HLine = shift (@HeaderLines)) {
			chomp ($HLine);
			my $hatfield_1 = "TOTAL ENERGY";
			if ( $HLine =~/$hatfield_1/ ){
				$energy = $HLine;
				my @words_1 = split (" ",$energy);
				push (@energyy,$words_1[3]);
			}
		}
	}
	my @idx = sort { $energyy[$a] <=> $energyy[$b] } 0 .. $#energyy;
	my @energyy_1     = @energyy[@idx];
	my @array_1       = @array[@idx];
	for ( my $i = 0; $i < scalar(@array_1); $i++) {
		open (HEADER,"$array_1[$i]") or die "Unable to open $i";
		my @HeaderLines = <HEADER>;
		close HEADER;
		my $count_lines = 0;
		my $first_line  = 0;
		my @array_lines = ();
		while (my $HLine = shift (@HeaderLines)) {
			my $hatfield = "FINAL GEOMETRY OBTAINED";
			if ( $HLine =~/$hatfield/ ){
				$first_line = $count_lines;			
			}
			$count_lines++;
			push (@array_lines,$HLine);
		}
		my $tmp_rest = $first_line + 3;	
		my $concat;
		for ( my $i = $tmp_rest; $i < $count_lines; $i++) {
			my $wordlength = length ($array_lines[$i]);
			if ( $wordlength > 3) {
				chomp ($array_lines[$i]);
				my @words = split (" ",$array_lines[$i]);
				$concat.= "$words[0]\t$words[1]\t$words[3]\t$words[5]\n";
			}
		}
		push (@array_coords_mopac,$concat);
	}
	return (\@array_coords_mopac,\@array_1,\@energyy_1);
}
###################################
# Keywords Errors
sub errors_config {
	my ($data) = @_;
	my $bolean = 1;
	if ( ( @$data[0]  =~/numb_conf/gi )        ){ } else { print "ERROR Correct Keywords: numb_conf\n";         $bolean = 0;};
	if ( ( @$data[1]  =~/box_size/gi )         ){ } else { print "ERROR Correct Keywords: box_size\n";          $bolean = 0;};
	if ( ( @$data[2]  =~/chemical_formula/gi ) ){ } else { print "ERROR Correct Keywords: chemical_formula\n";  $bolean = 0;};
	if ( ( @$data[3]  =~/mutations/gi )        ){ } else { print "ERROR Correct Keywords: mutations\n";         $bolean = 0;};
	if ( ( @$data[4]  =~/core_mem/gi )         ){ } else { print "ERROR Correct Keywords: core_mem\n";          $bolean = 0;};
	if ( ( @$data[5]  =~/charge_multi/gi)      ){ } else { print "ERROR Correct Keywords: charge_multi\n";      $bolean = 0;};
	if ( ( @$data[6]  =~/header/gi)            ){ } else { print "ERROR Correct Keywords: header\n";            $bolean = 0;};
	if ( ( @$data[8]  =~/software/gi)          ){ } else { print "ERROR Correct Keywords: software\n";          $bolean = 0;};
	if ( ( @$data[7]  =~/crossing_over/gi)     ){ } else { print "ERROR Correct Keywords: crossing_over\n";     $bolean = 0;};		
	return $bolean;
}
###################################
# join string
sub string_tmp {
	my ($array_input) = @_;
	#
	my $concat_string;
	for ( my $i=2 ; $i < scalar (@{$array_input}); $i++) {
		$concat_string.="@$array_input[$i] ";
	}
	return $concat_string;
}
###################################
# Verification
sub verification_no_opt {
	my ($a1, $a2, $dist)=@_;
	# hash values	
	my $v1  = $Atomic_radii{$a1} || $other_element; 
	my $v2  = $Atomic_radii{$a2} || $other_element;
	my $sum = $v1 + $v2;  
	my $resultado;
	# steric effects if radio1+radio2 < distance
	if($dist <= $sum){
		# Steric problem	
		$resultado = 1; 
	}else{
		$resultado = 0;
	}
	return $resultado;
}
###################################
# Verification
sub verification_opt {
	my ($a1, $a2, $dist)=@_;
	# hash values
	my $sum = 0.7;  
	my $resultado;
	# steric effects if radio1+radio2 < distance
	if($dist <= $sum){
		# Steric problem	
		$resultado = 1; 
	}else{
		$resultado = 0;
	}
	return $resultado;
}
###################################
# Steric Impediment for atoms
sub steric_impediment_atoms {
	# array are send by reference
	my ($hash_tmp,$option) = @_;
	# get size
	my $final_trial = 0;
	my $resultado   = 0;
	my @total       = @{$hash_tmp};
	#
	for (my $i=0; $i < scalar(@total);$i++){
		for (my $j=0; $j < scalar(@total); $j++){
			if ( $i < $j ){
				my ($atom_1,$axis_x_1,$axis_y_1,$axis_z_1) = split '\s+', $total[$i];
				my ($atom_2,$axis_x_2,$axis_y_2,$axis_z_2) = split '\s+', $total[$j];
                #
                #                
				my $distance    = Euclidean_distance($axis_x_1,$axis_y_1,$axis_z_1,$axis_x_2,$axis_y_2,$axis_z_2);
				if ( $option == 0) {
					$final_trial = $final_trial + verification_no_opt($atom_1,$atom_2,$distance);
				} else {
					$final_trial = $final_trial + verification_opt($atom_1,$atom_2,$distance);
				}
                #
                #
				if( $final_trial ==	 1 ) {
					$resultado = 1;
					last;
				}
                #
			}
		}
	}
	# verify for steric impediment, 1 yes, 0 no;
	return $resultado;		
}
###################################
# compute the center of mass
sub measure_center {
	my ($coord_x,$coord_y,$coord_z) = @_;
	my $num_data = scalar (@{$coord_x});
	my @array  = ();
	my $weight = 1;
	# variable sum
	my $sum_weight = 0;
	my $sum_x = 0;
	my $sum_y = 0;
	my $sum_z = 0;
	for ( my $j = 0 ; $j < $num_data ; $j = $j + 1 ){
		$sum_weight+= $weight;
		$sum_x+= $weight * @$coord_x[$j];
		$sum_y+= $weight * @$coord_y[$j];
		$sum_z+= $weight * @$coord_z[$j];		
	}
	my $com_x = $sum_x / $sum_weight;
	my $com_y = $sum_y / $sum_weight;
	my $com_z = $sum_z / $sum_weight;
	# array
	@array = ($com_x,$com_y,$com_z);
	# return array	
	return @array;
}
###################################
# Returns the additive inverse of v(-v)
sub vecinvert {
	my ($center_mass) = @_;
	my @array         = ();
	foreach my $i (@$center_mass) {
		my $invert        = $i * -1;
		$array[++$#array] = $invert; 
	}	
	# return array	
	return @array;
}
###################################
# Returns the vector sum of all the terms.
sub vecadd {
	my ($coord_x,$coord_y,$coord_z,$vecinvert_cm ) = @_;
	my $num_data = scalar (@{$coord_x});
	my @array   = ();
	my $sum_coord_x;
	my $sum_coord_y;
	my $sum_coord_z;
	# array 
	my @array_x = ();
	my @array_y = ();
	my @array_z = ();
	for ( my $i = 0 ; $i < $num_data ; $i = $i + 1 ){	
		$sum_coord_x = @$coord_x[$i]+@$vecinvert_cm[0] ; 
		$sum_coord_y = @$coord_y[$i]+@$vecinvert_cm[1] ;
		$sum_coord_z = @$coord_z[$i]+@$vecinvert_cm[2] ;
		# save array
		$array_x[++$#array_x] = $sum_coord_x;
		$array_y[++$#array_y] = $sum_coord_y;
		$array_z[++$#array_z] = $sum_coord_z;
	}
	@array = ( [@array_x], 
              [@array_y], 
              [@array_z] ); 
	# return array	
	return @array;
}
###################################
# value covalent radii
sub value_covalent_radii {
	my ($element) = @_;
	my $radii_val = 0;
	if ( exists $Atomic_radii{$element} ) {
		$radii_val = $Atomic_radii{$element};
	} else {
		$radii_val = $other_element ;
	}	
	return $radii_val
}
###################################
# Automatic box length
sub automatic_box_length {
	my ($input_array) = @_;
	my $sum = 0;
	#
	foreach my $i (@{$input_array}) {
		my $radii_val;
		if ( exists $Atomic_radii{$i} ) {
			# exists
			$radii_val = $Atomic_radii{$i};
		} else {
			# not exists
			$radii_val = $other_element ;
		}
		$sum+=$radii_val;
	}
	return $sum;
}
###################################
# Euclidean distance between points
sub Euclidean_distance {
	# array coords basin 1 and basin 2
	my ($p1,$p2,$p3, $axis_x, $axis_y, $axis_z) = @_;
	# variables
	my $x1 = $axis_x;
	my $y1 = $axis_y;
	my $z1 = $axis_z;
	# measure distance between two point
	my $dist = sqrt(
					($x1-$p1)**2 +
					($y1-$p2)**2 +
					($z1-$p3)**2
					); 
	return $dist;
}
###################################
# function to determine if a string is numeric
sub looks_like_number {
	my ($array_input) = @_;	
	#
	my $number;
	my $element;
	my @array_elements = ();
	for ( my $i=2 ; $i < scalar (@{$array_input}); $i++) {
		$number  = 0;
		if ( @$array_input[$i] =~ /^[0-9,.E]+$/ ) {
			$number = @$array_input[$i];
		} else {
			$element = @$array_input[$i];
		}
		for ( my $j = 0 ; $j < $number ; $j++) {
			push (@array_elements,$element);
		}
	}
	return @array_elements;
}
###################################
# fisher yates shuffle
sub fisher_yates_shuffle {
	my $array = shift;
	my $i = @$array;
	while ( --$i ) {
		my $j = int rand( $i+1 );
		@$array[$i,$j] = @$array[$j,$i];
	}
	return @$array; 
}
###################################
# Construct Discrete Search Space Cube
sub Construct_Discrete_Search_Space_Cube {
	my ($length_side_box_x,$length_side_box_y,$length_side_box_z,$cell_size_w, $sphereType) = @_;
	#
	my $max_coord_x = $length_side_box_x;
	my $max_coord_y = $length_side_box_y;
	my $max_coord_z = $length_side_box_z;	
	#
	my $cell_center          = ( $cell_size_w / 2 ); 
	#
	# For spheres
	my $middleX = $length_side_box_x/2;	
	my $middleY = $length_side_box_y/2;
	my $middleZ = $length_side_box_z/2;
	#
	my $total_number_of_cell_x = ( $max_coord_x / $cell_size_w );
	my $total_number_of_cell_y = ( $max_coord_y / $cell_size_w );
	my $total_number_of_cell_z = ( $max_coord_z / $cell_size_w );
	#
	my @discretized_search_space = ();
	#
	my $pm = Parallel::ForkManager->new($nprocess);
	$pm->run_on_finish(sub {
						my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $string_coords) = @_;
							#
							my @string_array = ();
							#
							if (!defined($$string_coords)) {
							} else {
								@string_array = split ( /\n/, $$string_coords );
								foreach (@string_array) {
									push (@discretized_search_space,$_);
								}							
							}
						});	
	#
	for (my $x=0; $x < $total_number_of_cell_x ; $x++) {
		$pm->start and next;
		# All children process havee their own random.
		srand();
		my $string_coords;	
		for (my $y=0; $y < $total_number_of_cell_y ; $y++) { 
			for (my $z=0; $z < $total_number_of_cell_z ; $z++) { 
				my $coord_x = ( ( $x * $cell_size_w ) + $cell_center );
				my $coord_y = ( ( $y * $cell_size_w ) + $cell_center );
				my $coord_z = ( ( $z * $cell_size_w ) + $cell_center );
				# Si existe radio de esfera entra 
				if($sphereType != -1){
					my $distance_sphere = Euclidean_distance($coord_x, $coord_y, $coord_z, $middleX, $middleY, $middleZ);
					if ( ($distance_sphere < $sphereType) || ($distance_sphere > ($sphereType + 2*$delta)) ) {
						# Delete all cells inside the radio of the sphere
					}else{
						#print "X\t$coord_x\t$coord_y\t$coord_z\n";
						$string_coords.= "$coord_x\t$coord_y\t$coord_z\n";	
					}					
				} else {
					$string_coords.= "$coord_x\t$coord_y\t$coord_z\n";
				}
			}
		}
		$pm->finish(0,\$string_coords);
	}
	$pm->wait_all_children;	
	#
	return \@discretized_search_space;
}
###################################
# Construct Solutions
sub Construct_Solutions {
	my ($number_atoms,$Atoms_set,$discretized_search_space) = @_;
	#
	my @array_all_atoms      = @{$Atoms_set};
	# Delete construct solutions for atoms
	@tmp_new_neighbors_atoms = ();
	#
	for (my $i=0; $i < $number_atoms ; $i++) {
		my $element_neighbors = $i;
		# Este >= y <= va si o si
		if ( $element_neighbors >= 1) {
			my @array_neighbors = ();
			for (my $j=0; $j <= $element_neighbors ; $j++) {
				push (@array_neighbors,$array_all_atoms[$j]);
			}
			Rules_Cell_Automata (\@array_neighbors,\@{$discretized_search_space});
			#
		}
	}
	# return final cell
	return @tmp_new_neighbors_atoms;
}
###################################
# Rules Cell Automata
sub Rules_Cell_Automata {
	my ($array_neighbors,$discretized_search_space) = @_;
	#
	my @tmp_neighbors      = @{$array_neighbors};
	my $length_neighbors   = scalar (@{$array_neighbors});
	#
	my @search_space_atoms = @{$discretized_search_space};
	my $length_array_dis   = scalar (@{$discretized_search_space});
	###############################
	# my first cell in the space
	#
	if ( scalar (@tmp_new_neighbors_atoms) == 0 ) {
		my $int_coord_rand           = int rand($length_array_dis);
		my ($axis_x,$axis_y,$axis_z) = split '\s+', $search_space_atoms[$int_coord_rand];
		#
		my $element_1    = $tmp_neighbors[0];
		my $element_2    = $tmp_neighbors[1];
		my $radii_atom_1 = value_covalent_radii ($element_1);
		my $radii_atom_2 = value_covalent_radii ($element_2);
		my $total_radii  = $radii_atom_1 + $radii_atom_2;
		my $delta_radii  = $delta + $total_radii;
		#
		my $line_1 = "$element_1\t$axis_x\t$axis_y\t$axis_z";
		#
		my @array_tmp = ();
		foreach ( @search_space_atoms ) {
			my ($tmp_x,$tmp_y,$tmp_z) = split '\s+', $_;
			my $distance = Euclidean_distance ($axis_x,$axis_y,$axis_z,$tmp_x,$tmp_y,$tmp_z);
			if ( $distance > $total_radii && $distance <= $delta_radii ) {
				push (@array_tmp,"$element_2\t$tmp_x\t$tmp_y\t$tmp_z");
			}
		}
		my $tmp_length  = scalar (@array_tmp);
		my $tmp_int 	= int rand($tmp_length);
		my ($tmp_element,$tmp_x,$tmp_y,$tmp_z) = split '\s+', $array_tmp[$tmp_int];
		#
		my $line_2 = "$tmp_element\t$tmp_x\t$tmp_y\t$tmp_z";
		#
		@tmp_new_neighbors_atoms = ($line_1,$line_2);		
	} else {
		# choose rand new neighbors
		my $length_new_neighbors = scalar (@tmp_new_neighbors_atoms);
		#
		my $line;
		my $boolean = 0;
		# do...while loop execution			
		do {
			my $int_new_neighbors_rand            = int rand($length_new_neighbors);
			my ($element,$axis_x,$axis_y,$axis_z) = split '\s+', $tmp_new_neighbors_atoms[$int_new_neighbors_rand];
			# The last element of an array
			my $last_element = $tmp_neighbors[$#tmp_neighbors];
            #
            #
            my $radii_atom_1 = value_covalent_radii ($element);
            my $radii_atom_2 = value_covalent_radii ($last_element);
            my $total_radii  = $radii_atom_1 + $radii_atom_2;
            my $delta_radii  = $delta + $total_radii;
            #
            my @array_tmp = ();        
            foreach ( @search_space_atoms ) {
                my ($tmp_x,$tmp_y,$tmp_z) = split '\s+', $_;
                my $distance = Euclidean_distance ($axis_x,$axis_y,$axis_z,$tmp_x,$tmp_y,$tmp_z);
                if ( $distance > $total_radii && $distance <= $delta_radii ) {            
                    push (@array_tmp,"$last_element\t$tmp_x\t$tmp_y\t$tmp_z");
                }
            }
            #
            my $tmp_length  = scalar (@array_tmp);
            #
            my $tmp_int 	= int rand($tmp_length);
            my ($tmp_element,$tmp_x,$tmp_y,$tmp_z) = split '\s+', $array_tmp[$tmp_int];
            #
            my @arr = ();
            foreach (@tmp_new_neighbors_atoms) {
                push (@arr,$_);
            }
            $line = "$tmp_element\t$tmp_x\t$tmp_y\t$tmp_z";
            push (@arr,$line);
            #
            # Rules
            $boolean = steric_impediment_atoms (\@arr,0);
        #    
		} while ( $boolean > 0 );
		push (@tmp_new_neighbors_atoms,$line);
	}
}
###################################
# Center system
sub Center_All_Cell {
	my ($array_CS,$option) = @_;
	#
	my @Atoms   = ();
	my @coord_x = ();
	my @coord_y = ();	
	my @coord_z = ();
	#
	my @total_array = ();
	foreach my $line (@{$array_CS}) {
		if ( $option == 1 ) {
			my ($element,$axis_x,$axis_y,$axis_z) = split '\s+', $line;
			$Atoms[++$#Atoms]       = $element;
			$coord_x[++$#coord_x]   = $axis_x;
			$coord_y[++$#coord_y]   = $axis_y;
			$coord_z[++$#coord_z]   = $axis_z;
			#
		} else {
			my ($axis_x,$axis_y,$axis_z) = split '\s+', $line;
			$coord_x[++$#coord_x]   = $axis_x;
			$coord_y[++$#coord_y]   = $axis_y;
			$coord_z[++$#coord_z]   = $axis_z;
		} 	
	}
	#
	# for coords xyz molecules, moveby {x y z} (translate selected atoms to origin)
	my @array_center_mass = measure_center(\@coord_x,\@coord_y,\@coord_z);
	my @array_vecinvert   = vecinvert(\@array_center_mass);
	my @array_catersian   = vecadd (\@coord_x,\@coord_y,\@coord_z,\@array_vecinvert);
	if ( $option == 1 ) {		
		for ( my $i = 0 ; $i < scalar (@Atoms) ; $i = $i + 1 ){
			push (@total_array,"$Atoms[$i]\t$array_catersian[0][$i]\t$array_catersian[1][$i]\t$array_catersian[2][$i]");
		}
	} else {
		for ( my $i = 0 ; $i < scalar (@coord_x) ; $i = $i + 1 ){
			push (@total_array,"$array_catersian[0][$i]\t$array_catersian[1][$i]\t$array_catersian[2][$i]");
		}
	}	
	#
	return (@total_array);
}
###################################
# Rules Automata Hydrogen and Halogens
sub New_Rules_Automaton {
	my ($Atoms) = @_;
    my @Array_Elements = @{$Atoms};
    #
    my @arr_Hydrogen          = ();
    my @arr_Halogens          = ();
    my @arr_Alcalinos         = ();
    my @arr_Alcalinos_Terreos = ();
    my @arr_Others            = ();
    #
    foreach my $i (@Array_Elements) {
        #
        if ( $i eq 'H') {
            push (@arr_Hydrogen, $i);
        } elsif ( ($i eq 'F') or ($i eq 'Cl') or ($i eq 'Br') or ($i eq 'I') or ($i eq 'At') ) {
            push (@arr_Halogens, $i);
        } elsif ( ($i eq 'Li') or ($i eq 'Na') or ($i eq 'K') or ($i eq 'Rb') or ($i eq 'Cs') or ($i eq 'Fr') ) {
            push (@arr_Alcalinos, $i);
        } elsif ( ($i eq 'Be') or ($i eq 'Mg') or ($i eq 'Ca') or ($i eq 'Sr') or ($i eq 'Ba') or ($i eq 'Ra') ) {
            push (@arr_Alcalinos_Terreos, $i);
        } else {
            push (@arr_Others, $i);
        } 
    }
    my @shuffled_atoms = fisher_yates_shuffle(\@arr_Others);
    my @Total_Atoms    = (@shuffled_atoms,@arr_Hydrogen,@arr_Alcalinos_Terreos,@arr_Alcalinos,@arr_Halogens);
    #
    return @Total_Atoms;
}
###################################
# Automata 1D
sub z_axis_path {
	my ($largeAxisZ) = @_;
	#
	my $window = 1.8;
	#
    my $cell_center = 0; 
    my $total_number_of_cell_z = ( (4*$largeAxisZ) / $window );
    $total_number_of_cell_z=~ s/\.\d+$//;
    #
    my @array = ();
    #
    my $div_points_z = ($total_number_of_cell_z / 2);
    $div_points_z=~ s/\.\d+$//;
    my $neg_points_z = $div_points_z * -1;
    my $pos_points_z = $div_points_z;
    #
    for (my $z=$neg_points_z; $z <= ($pos_points_z + 0.0001); $z++) { 
        my $coord_z = ( ( $z * $window ) + $cell_center );
        #
        my $c_x = sprintf '%.6f', 0;
        my $c_y = sprintf '%.6f', 0;
        my $c_z = sprintf '%.6f', $coord_z;
		#
        my $xyz_axis = "$c_x\t$c_y\t$c_z";
        #
        push (@array,$xyz_axis);
    }
    return @array;
}
###################################
# Construct Solution in 1 Dimension
sub Construct_Solution_1D {
	my ($number_atoms,$Atoms_set,$discretized_search_space) = @_;
	#
	my @array_all_atoms  = @{$Atoms_set};
	my @array_coords_xyz = @{$discretized_search_space};
	#
	my @total_specie_1D = (); 
	#
	for (my $i=0; $i < $number_atoms ; $i++) {		
		my $line = "$array_all_atoms[$i]\t$array_coords_xyz[$i]";
		#
		push (@total_specie_1D, $line);
	}
	return @total_specie_1D;
}
###################################
# Generar poblacion Automata celular
sub Opt_Cell_Automata {
	my ($Num_of_geometries_input,$Atoms,$discretized_search_space,$Num_of_atoms,$id_name,$switchAU,$dimension) = @_;
	#
	my @cartesian_coords = ();
	my @filename         = ();
	#
	for (my $iteration = 0 ; $iteration < $Num_of_geometries_input ; $iteration++) {
		my $id               = sprintf '%.6d', $count_struc_global;
		my $filebase         = "$id_name\_$id";
		push (@filename,$filebase);
		$count_struc_global++;
	}
	#
	my $pm = Parallel::ForkManager->new($nprocess);
	#my $pm = Parallel::ForkManager->new(10);
	$pm->run_on_finish(sub {
						my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $string_atoms_coords) = @_;
							push (@cartesian_coords,$$string_atoms_coords);
						});	
	#
	for (my $iteration = 0 ; $iteration < $Num_of_geometries_input ; $iteration++) {
		$pm->start and next;
		# All children process havee their own random.			
		srand();
        #
        my @shuffled_atoms = ();
		# New Rules Automaton
        if ( $switchAU == 1 ) {
			@shuffled_atoms = New_Rules_Automaton (\@{$Atoms});
		} else {
			# Old Rules Automaton
			@shuffled_atoms = fisher_yates_shuffle(\@{$Atoms});
		}
		#
        #
        my @array_CS  = ();
        if ( $dimension == 1 ) {
			@array_CS  = Construct_Solution_1D ($Num_of_atoms,\@shuffled_atoms,\@{$discretized_search_space});
		} else {
			@array_CS  = Construct_Solutions ($Num_of_atoms,\@shuffled_atoms,\@{$discretized_search_space});
		}
        #
		my @final_construct_molecule = Center_All_Cell (\@array_CS,1);
		#	
		my $string_atoms_coords;
		for ( my $i = 0 ; $i < scalar (@final_construct_molecule) ; $i = $i + 1 ) {
			my ($element,$axis_x,$axis_y,$axis_z) = split '\s+', $final_construct_molecule[$i];
			my $axis_x_oxs  = sprintf '%.4f', $axis_x;
			my $axis_y_oxs  = sprintf '%.4f', $axis_y;
			my $axis_z_oxs  = sprintf '%.4f', $axis_z; 
			$string_atoms_coords.="$element  $axis_x_oxs  $axis_y_oxs  $axis_z_oxs\n";
		}
        #print "$string_atoms_coords\n\n";
		#push (@cartesian_coords,$string_atoms_coords);
		$pm->finish(0,\$string_atoms_coords);
	}
	# Paralel
	$pm->wait_all_children;
	#
	return (\@filename,\@cartesian_coords);
}
###################################
# Generar mutaciones en la poblacion
sub Mutations_Automata {
	my ($Atoms,$Num_of_atoms) = @_;
	#
	my @tmp           = split (/\n/,$Atoms);
	my @recombination = fisher_yates_shuffle (\@tmp);
	#
	my @name          = ();
	my @coords_x      = ();
	my @coords_y      = ();
	my @coords_z	  = ();
	#
	my $delta_kick = $kick_punch;
	#
	my $atoms_to_change = int($Num_of_atoms * $PCENTATOMMUTATED);
	for (my $i = 0; $i < $Num_of_atoms; $i++) {
		my($element, $cx, $cy, $cz) = split (/\s+/,$recombination[$i]);
		push (@name, $element);
		if($i < $atoms_to_change){
			my $radii_atom  = value_covalent_radii ($element);
			my $total_radii = ($radii_atom + $radii_atom);
			my $delta_radii = $delta_kick + $total_radii ;
			#
			my $mult = int(rand(2));
			if ( $mult == 0 ) { $mult=-1; }
			my $num     = ($radii_atom + rand($delta_radii - $radii_atom))*$mult;
			#
			# print "dX $num\n";
			my $sum_x   = $num + $cx ;
			$mult       = int(rand(2));
			if ( $mult == 0 ){ $mult = -1; }
			$num        = ($radii_atom + rand($delta_radii - $radii_atom))*$mult;
		    # print "dY $num\n";
			my $sum_y = $num + $cy ;
			$mult       =int(rand(2));
			if ( $mult == 0 ){ $mult = -1;}
			$num        = ($radii_atom + rand($delta_radii - $radii_atom))*$mult;
			# print "dZ $num\n";
			my $sum_z = $num + $cz;
			#
			#print "$element\t$sum_x\t$sum_y\t$sum_z\n";
			push @coords_x, $sum_x;
			push @coords_y, $sum_y;
			push @coords_z, $sum_z;
		}else{
			#print "$element\t$cx\t$cy\t$cz\n";
			push @coords_x, $cx;
			push @coords_y, $cy;
			push @coords_z, $cz;
		}
	}
	my $string;
	my @total_mol = ();
	for (my $i = 0; $i < $Num_of_atoms; $i++) {
		my $mi_x = sprintf '%.6f',$coords_x[$i];
		my $mi_y = sprintf '%.6f',$coords_y[$i];
		my $mi_z = sprintf '%.6f',$coords_z[$i];
		push (@total_mol,"$name[$i]  $mi_x  $mi_y  $mi_z");
		$string.="$name[$i]  $mi_x  $mi_y  $mi_z\n";
	}
	return (\@total_mol,$string);
}
###################################
# Generation of random coordinates (Type KICK)
sub gen_xyz {
	my ($Box_x, $Box_y, $Box_z) = @_;
	# generate a random number in perl in the range box size
	my $lower_limit_x = ($Box_x * -1);
	my $upper_limit_x = $Box_x;
	my $lower_limit_y = ($Box_y * -1);
	my $upper_limit_y = $Box_y;
	my $lower_limit_z = ($Box_z * -1);
	my $upper_limit_z = $Box_z;
	#
	my $x = (rand($upper_limit_x-$lower_limit_x) + $lower_limit_x);
	my $y = (rand($upper_limit_y-$lower_limit_y) + $lower_limit_y);
	my $z = (rand($upper_limit_z-$lower_limit_z) + $lower_limit_z);
	#
	my $x_coord = sprintf '%.6f', $x;
	my $y_coord = sprintf '%.6f', $y;
	my $z_coord = sprintf '%.6f', $z;
	my @coords  = ($x_coord, $y_coord, $z_coord);
	return @coords;
}
###################################
# Generation of random pob (Type KICK)
sub pob_type_kick {
	my ($side_plus_x,$side_plus_y,$side_plus_z,$Atoms,$num_cycle,$id_name) = @_;
	#
	my @tmp             = ();
	my @file_base_array = ();
	#
	for (my $i = 0 ; $i < $num_cycle; $i++) {
		my @array_CS = ();
		my $boolean  = 0;
		while ( $boolean < 1 ) {
			@array_CS   = ();
			foreach my $species (@{$Atoms}) {
				my @xyz = gen_xyz ($side_plus_x,$side_plus_y,$side_plus_z);
				push (@array_CS,"$species  $xyz[0]  $xyz[1]  $xyz[2]");
			}
			my $option_boo = steric_impediment_atoms (\@array_CS,0);
			if ( $option_boo == 0 ) { $boolean = 1;}
		}
		#
		my @final_construct_molecule = Center_All_Cell     (\@array_CS,1);
		my $string;
		foreach (@final_construct_molecule) {
			$string.="$_\n";
		}
		my $id               = sprintf '%.6d', $count_struc_global;
		my $filebase         = "$id_name$id";
		push (@tmp,$string);
		push (@file_base_array,$filebase);
		$count_struc_global++;
	}
	#
	return (\@file_base_array,\@tmp);
}
# 
sub LammpsInput {
	#
	my $filebase         = $_[0];
	my $coordsMM         = $_[1];
	my $LammpsInput      = "$filebase.in";
	my $LammpsCoords     = "$filebase.dat";
	my $LammpsOutputAxis = "$filebase.xyz";
	my $iteration        = $_[2];
	my $Headerfile       = $_[3];
	my $Box_Length       = $_[4];
	#
	my @words = split (/\n/,$coordsMM);
	my %num_atoms_lammps = ();
	#
	open (COORDSFILE, ">$LammpsCoords");
	print COORDSFILE "# $LammpsCoords file format coords\n";
	print COORDSFILE "\n";
	my $count    = 0;
	my @elements = ();
	foreach my $i (@words){
		my @axis    = split (" ",$i);
		my $label  = $axis[0];  
		push (@elements,$label);
		$count++;
	}
	print COORDSFILE "$count atoms\n";
	my @unique_elements   = uniq @elements;
	my $num_uniq_elements = scalar (@unique_elements);
	print COORDSFILE "$num_uniq_elements atom types\n";
	print COORDSFILE "\n";
	my ($mi_x,$mi_y,$mi_z,$ma_x,$ma_y,$ma_z) = split (" ",$Box_Length);
	print COORDSFILE " $mi_x   $ma_x     xlo xhi\n";
	print COORDSFILE " $mi_y   $ma_y     ylo yhi\n";
	print COORDSFILE " $mi_z   $ma_z     zlo zhi\n";
	print COORDSFILE "\n";
	print COORDSFILE " Masses\n";
	print COORDSFILE "\n";
	my $numb_at = 1;
	for (my $i = 0 ; $i < $num_uniq_elements; $i++) {
		my $mass_val  = 0;
		my $element   = $unique_elements[$i];
		if ( exists $Atomic_mass{$element} ) {
			$mass_val = $Atomic_mass{$element};
		} else {
			$mass_val = $other_element ;
		}
		$num_atoms_lammps{$element} = $numb_at;
		print COORDSFILE " $numb_at $mass_val\n";
		$numb_at++;		
	}
	print COORDSFILE "\n";
	print COORDSFILE " Atoms\n";
	print COORDSFILE "\n";	
	my $count_atoms = 1;
	my $tmp   = "0.0";
	foreach my $i (@words){
		my @axis    = split (" ",$i);
		#
		my $label  = $axis[0];  
		my $axis_x = $axis[1];
		my $axis_y = $axis[2];
		my $axis_z = $axis[3];
		#
		print COORDSFILE "   $count_atoms  $num_atoms_lammps{$label}  $tmp  $axis_x  $axis_y  $axis_z\n";
		$count_atoms++;
	}
	print COORDSFILE "\n";
	close (COORDSFILE);
	#
	open (LAMMPSFILE, ">$LammpsInput");
	print LAMMPSFILE "# REAX FF parameters\n";
	print LAMMPSFILE "#\n";
	print LAMMPSFILE "dimension       3\n";
	print LAMMPSFILE "boundary        p p p\n";
	print LAMMPSFILE "units		    real\n";
	print LAMMPSFILE "\n";
	print LAMMPSFILE "atom_style      charge\n";
	print LAMMPSFILE "atom_modify     map array sort 0 0.0\n";
	print LAMMPSFILE "\n";
	print LAMMPSFILE "read_data	    $LammpsCoords\n";
	print LAMMPSFILE "\n";
	print LAMMPSFILE "pair_style	    reax/c NULL\n";
	#
	print LAMMPSFILE "pair_coeff	    * * $Headerfile";
	for (my $i = 0 ; $i < $num_uniq_elements; $i++) {
		my $mass_val  = 0;
		my $element   = $unique_elements[$i];
		print LAMMPSFILE "$element ";
	}		
	print LAMMPSFILE "\n";	
	#
	print LAMMPSFILE "neighbor	    2 bin\n";
	print LAMMPSFILE "neigh_modify	every 10 check yes\n";
	print LAMMPSFILE "fix             2 all qeq/reax 1 0.0 10.0 1e-6 reax/c\n";
	print LAMMPSFILE "\n";
	print LAMMPSFILE "# should equilibrate much longer in practice\n";
	print LAMMPSFILE "#\n";
	print LAMMPSFILE "fix             1 all nvt temp 300 300 10\n";
	print LAMMPSFILE "#fix		        1 all npt temp 273.0 273.0 10.0 iso 1.0 1. 2000.0\n";
	print LAMMPSFILE "timestep        0.2\n";
	print LAMMPSFILE "thermo_style    custom step temp epair pe ke etotal \n";
	print LAMMPSFILE "thermo          1\n";
	print LAMMPSFILE "#\n";
	print LAMMPSFILE "dump             1 all custom 100 $LammpsOutputAxis element x y z\n";
	#
	print LAMMPSFILE "dump_modify      1 element ";
	for (my $i = 0 ; $i < $num_uniq_elements; $i++) {
		my $mass_val  = 0;
		my $element   = $unique_elements[$i];
		print LAMMPSFILE "$element ";
	}		
	print LAMMPSFILE "\n";		
	#
	print LAMMPSFILE "\n";
	print LAMMPSFILE "# $steps step of minimize\n";	
	print LAMMPSFILE "minimize           $criteria $criteria $steps $steps\n";
	print LAMMPSFILE "\n";
	#
	return $filebase;
}
# 
sub submit_queue_lammps {
	my ($arrayInputs,$path_lammps) = @_;
	my @Inputs_lammps  = @{$arrayInputs};
	my @Outputs_lammps = ();
	for ( my $i=0 ; $i < scalar(@Inputs_lammps); $i++) {
		my $env = `$path_lammps -in $Inputs_lammps[$i].in >$Inputs_lammps[$i].out`;
		push (@Outputs_lammps,"$Inputs_lammps[$i].out");
	}
	return @Outputs_lammps
}
#
sub coords_energy_lammps {
	my ($num_atoms_xyz,$arrayInputs) = @_;	
	#
	my @array              = ();
	my @array_coords_mopac = ();
	my @energy_lmp         = ();
	my @coords_lmp         = ();
	foreach my $files ( @{$arrayInputs} ) {
		(my $without_extension = $files) =~ s/\.[^.]+$//;
		open(HEADER,"$without_extension.out") or die "Unable to open $files";
		my @HeaderLines  = <HEADER>;
		close HEADER;
		#
		my $count_lines = 0; 
		for ( my $i=0; $i < scalar(@HeaderLines); $i++) {		
			if ( ($HeaderLines[$i]=~/Energy/gi ) && ($HeaderLines[$i]=~/initial/gi ) 
			  && ($HeaderLines[$i]=~/next-to-last/gi ) && ($HeaderLines[$i]=~/final/gi ) ){
				$count_lines = $i;
			}		
		}
		my $data_E = $count_lines + 1;
		my ($E_initial,$E_next_to_last,$E_final) = split (" ",$HeaderLines[$data_E]);
		push (@energy_lmp,$E_final);
	}
	#
	foreach my $files ( @{$arrayInputs} ) {
		(my $without_extension = $files) =~ s/\.[^.]+$//;
		open(HEADER,"$without_extension.xyz") or die "Unable to open $files";
		my @HeaderLines  = <HEADER>;
		close HEADER;
		#
		my $string_coords;
		my $count_lines = 0; 
		for ( my $i=0; $i < scalar(@HeaderLines); $i++) {		
			if ( ($HeaderLines[$i]=~/ITEM/gi ) && ($HeaderLines[$i]=~/ATOMS/gi ) 
			  && ($HeaderLines[$i]=~/element/gi ) && ($HeaderLines[$i]=~/x/gi ) ){
				$count_lines = $i;
			}		
		}
		my $data_xyz = $count_lines + $num_atoms_xyz;		
		for ( my $i = ($count_lines + 1); $i <= $data_xyz; $i++) {
			chomp ($HeaderLines[$i]);
			#print "$HeaderLines[$i]\n";
			$string_coords.="$HeaderLines[$i]\n";
		}
		push (@coords_lmp,$string_coords);
	}
	return (\@coords_lmp,\@energy_lmp);
}
###################################
#
sub EnergyFilter {
	my ($energy_opt,$coords_opt,$file_opt, $option_unit, $kcalmolfilt) = @_;
	# sort, same thing in reversed order
	my @array_energy_zero = @{$energy_opt}; 
	my @array_coords      = @{$coords_opt};
	my @array_files       = @{$file_opt};	
	#
	my @value_energy_sort = ();
	my @value_coords_sort = ();
	my @value_files_sort  = ();	
	my @idx = sort { $array_energy_zero[$a] <=> $array_energy_zero[$b] } 0 .. $#array_energy_zero;
	@value_energy_sort = @array_energy_zero[@idx];
	@value_coords_sort = @array_coords[@idx];
	@value_files_sort  = @array_files[@idx];	
    #
    my @arr_energy = ();  
    my @arr_coords = ();
    my @arr_files  = ();
    #
	for (my $i=0; $i < scalar (@value_energy_sort); $i++){
		my $resta = abs($value_energy_sort[0]) - abs($value_energy_sort[$i]);	
        # 1 Hartree = 27,2114 ev
        # 1 Hartree = 627,509 Kcal/mol
		my $Kcalmol = 0;
		if ($option_unit == 0) {
			$Kcalmol = sprintf("%.4f",( 627.509   * $resta ));
            #
            if ( $Kcalmol <= $kcalmolfilt ) {    
                push (@arr_energy, $value_energy_sort[$i]);
                push (@arr_coords, $value_coords_sort[$i]);
                push (@arr_files , $value_files_sort[$i]);
            }
		} elsif ( $option_unit == 1 ) {
			$Kcalmol = sprintf("%.4f",( 23.060542 * $resta ));
            #
            if ( $Kcalmol <= $kcalmolfilt ) {    
                push (@arr_energy, $value_energy_sort[$i]);
                push (@arr_coords, $value_coords_sort[$i]);
                push (@arr_files , $value_files_sort[$i]);
            }            
		} else {
			$Kcalmol = sprintf("%.4f", $resta);
            #
            if ( $Kcalmol <= $kcalmolfilt ) {    
                push (@arr_energy, $value_energy_sort[$i]);
                push (@arr_coords, $value_coords_sort[$i]);
                push (@arr_files , $value_files_sort[$i]);
            }            
		}
	}
    return (\@arr_energy,\@arr_coords,\@arr_files);
}
###################################
#
sub CycleCoords {
	my ($InputCoords,$InputFiles,$numb,$NameCycleFile) = @_;
    #
    my @arrCoords = @{$InputCoords};
    my @arrFiles  = @{$InputFiles};
    #
    my @Files_sort  = ();
	my @Coords_sort = ();
	my @idx = sort { $arrFiles[$a] cmp $arrFiles[$b] } 0 .. $#arrFiles;
	@Files_sort  = @arrFiles[@idx];
	@Coords_sort = @arrCoords[@idx];
    #
    open (FILECYCLE, ">$NameCycleFile.xyz");
    for ( my $i = 0; $i < scalar (@Files_sort); $i++) {
        print FILECYCLE "$numb\n";
        (my $without_extension = $Files_sort[$i]) =~ s/\.[^.]+$//;
        print FILECYCLE "$without_extension\n";
        print FILECYCLE "$Coords_sort[$i]";
    }
    close (FILECYCLE);
}
###################################
# print logo
sub print_logo {
	print "\n";
	print "    ___  _   _ _____ ________  ___  ___ _____ _____ _   _             \n";
	print "   / _ \\| | | |_   _|  _  |  \\/  | / _ \\_   _|  _  | \\ | |        \n";
	print "  / /_\\ \\ | | | | | | | | | .  . |/ /_\\ \\| | | | | |  \\| |       \n";
	print "  |  _  | | | | | | | | | | |\\/| ||  _  || | | | | | . ` |           \n";
	print "  | | | | |_| | | | \\ \\_/ / |  | || | | || | \\ \\_/ / |\\  |       \n";
	print "  \\_| |_/\\___/  \\_/  \\___/\\_|  |_/\\_| |_/\\_/  \\___/\\_| \\_/  \n";
	print "\n";                                                        
	print "       A probabilistic cellular automaton method for a           \n";
	print "        novel program for the search of global minimum            \n";
	print "         structures of atomic clusters and molecules              \n";
	print "                    in the gas phase.                           \n\n\n";
	print "                       TiznadoLab\n";
	print "\n";
	my $datestring = localtime();
	print   "              $datestring\n\n";
}







# # # # # # # # #
# Main
if (not defined $file_name) {
	print_logo ();
	die "\nAUTOMATON must be run with:\n\nUsage:\n\t perl AUTOMATON.pl [configure-file]\n\n\n";
	exit;  
}
#
if (not defined $queue) {
	print_logo ();
	die " Queue\n\n\n";
	exit;
}


#
my @delete_rot;
my @arrayOutputs = ();
my @arrayOptMM   = ();
my @EnergyOptMM  = ();
# read and parse files
my @data          = read_file($file_name);
my @arrays_errors = ();
# data parse
my $Num_of_geometries_input  = 0;
#
my $Box;
my @Box_dimensions = ();
my ($Box_x,$Box_y,$Box_z);
my $option_box;
#
my @Atoms        = ();
my $Num_of_atoms = 0;
#
my $Submit_guff;
my @Submit_parameters = ();
my $ncpus;
my $mem;
#
my $charge_multi;
my @charge_multi_parameters = ();
my $Charge; 
my $Multiplicity;
#
my $header;
my $software;
#
my $sphereType=-1;
#
my $mutations;
my $crossing_over;
#
my $DeepAU;
#
my %Info_files_dir_global = ();
my @array_keys            = ();
#
my $file_tmp = "Output.log";
open (FILEREP, ">$file_tmp") or die "Unable to open XYZ file: $file_tmp";
#
print FILEREP "\n";
print FILEREP "    ___  _   _ _____ ________  ___  ___ _____ _____ _   _             \n";
print FILEREP "   / _ \\| | | |_   _|  _  |  \\/  | / _ \\_   _|  _  | \\ | |        \n";
print FILEREP "  / /_\\ \\ | | | | | | | | | .  . |/ /_\\ \\| | | | | |  \\| |       \n";
print FILEREP "  |  _  | | | | | | | | | | |\\/| ||  _  || | | | | | . ` |           \n";
print FILEREP "  | | | | |_| | | | \\ \\_/ / |  | || | | || | \\ \\_/ / |\\  |       \n";
print FILEREP "  \\_| |_/\\___/  \\_/  \\___/\\_|  |_/\\_| |_/\\_/  \\___/\\_| \\_/  \n";
print FILEREP "\n";                                                        
print FILEREP "        A probabilistic cellular automaton method for a           \n";
print FILEREP "        novel program for the search of global minimum            \n";
print FILEREP "         structures of atomic clusters and molecules              \n";
print FILEREP "                    in the gas phase.                           \n\n\n";
print FILEREP "                       TiznadoLab\n";
print FILEREP "\n";
my $datestring = localtime();
print FILEREP "              $datestring\n\n";
#
foreach my $a_1 (@data){
	if ( ($a_1=~/#/gi ) ){
	#	print "$a_1\n";
	} else {
		if ( ($a_1=~/numb_conf/gi ) ){
			my @tmp = ();
			@tmp    = split (/\s+/,$a_1);		
			# Identify empty string
			if (!defined($tmp[2])) {			
				print "ERROR input population numbers empty\n";
				exit;
			} else {
				$Num_of_geometries_input = $tmp[2];
			}	
			#
			$arrays_errors[0] = "numb_conf";
		}		
		if ( ($a_1=~/box_size/gi ) ){
			my @tmp = ();
			@tmp    = split (/\s+/,$a_1);
			# Identify empty string			
			if (!defined($tmp[2])) {
				print FILEREP "MESSAGE Automatic size box\n";
				$option_box = 0;
			} else {
				my $var_tmp = string_tmp (\@tmp);
				$Box = $var_tmp;
				@Box_dimensions = split(/,/, $Box);
				$Box_x = $Box_dimensions[0];
				$Box_y = $Box_dimensions[1];
				$Box_z = $Box_dimensions[2];
				$option_box = 1;
			}
			#
			$arrays_errors[1] = "box_size";			
		}
		if ( ($a_1=~/chemical_formula/gi ) ){
			my @tmp = ();
			@tmp    = split (/\s+/,$a_1);
			# Identify empty string
			if (!defined($tmp[2])) {
				print "ERROR chemical formula empty\n";
				exit;				
			} else {
				@Atoms = looks_like_number (\@tmp);
				$Num_of_atoms = scalar (@Atoms);
			}
			my ($val1,$val2,@chem_form) = @tmp;
			print FILEREP "MESSAGE Chemical Formula @chem_form\n";
			#
			$arrays_errors[2] = "chemical_formula";			
		}
		if ( ($a_1=~/mutations/gi ) ){
			my @tmp = ();
			@tmp    = split (/\s+/,$a_1);
			# Identify empty string
			if (!defined($tmp[2])) {
				print "ERROR option mutation empty\n";
				exit;
			} else {			
				$mutations = $tmp[2];
			}
			#
			$arrays_errors[3] = "mutations";			
		}
		if ( ($a_1=~/core_mem/gi ) ){
			my @tmp = ();
			@tmp    = split (/\s+/,$a_1);
			# Identify empty string
			if (!defined($tmp[2])) {
				# default cpus 1 and memory 1GB
				$ncpus             = 1;
				$mem               = 1;
			} else {
				my $var_tmp        = string_tmp (\@tmp);
				$Submit_guff       = $var_tmp;
				@Submit_parameters = split(/,/, $Submit_guff);
				$ncpus             = $Submit_parameters[0];
				$mem               = $Submit_parameters[1];
			}
			#
			$arrays_errors[4] = "core_mem";			
		}
		if ( ($a_1=~/charge_multi/gi ) ){
			my @tmp = ();
			@tmp    = split (/\s+/,$a_1);
			# Identify empty string
			if (!defined($tmp[2])) {
				# default multiplicity 1 and charge 0
				$Charge       = 0; 
				$Multiplicity = 1;
			} else {
				my $var_tmp   = string_tmp (\@tmp);
				$charge_multi = $var_tmp;
				@charge_multi_parameters = split(/,/, $charge_multi);
				$Charge       = $charge_multi_parameters[0]; 
				$Multiplicity = $charge_multi_parameters[1];
			}
			#
			$arrays_errors[5] = "charge_multi";			
		}
		if ( ($a_1=~/header/gi ) ){
			my @tmp = ();
			@tmp    = split (/\s+/,$a_1);
			# Identify empty string
			if (!defined($tmp[2])) {
				print "ERROR theory level empty\n";
				exit;
			} else {
				my $var_tmp = string_tmp (\@tmp);			
				$header     = $var_tmp;
			}
			#
			$arrays_errors[6] = "header";			
		}
		if ( ($a_1=~/crossing_over/gi ) ){
			my @tmp = ();
			@tmp    = split (/\s+/,$a_1);
			# Identify empty string
			if (!defined($tmp[2])) {
				print "ERROR option crossing over empty\n";
				exit;
			} else {			
				$crossing_over = $tmp[2];
			}
			#
			$arrays_errors[7] = "crossing_over";			
		}		
		if ( ($a_1=~/software/gi ) ){
			my @tmp = ();
			@tmp    = split (/\s+/,$a_1);
			# Identify empty string
			if (!defined($tmp[2])) {
				print "ERROR software empty\n";
				exit;
			} else {			
				$software = $tmp[2];
			}
			#
			$arrays_errors[8] = "software";			
		}
		if ( ($a_1=~/radius_sphere/gi ) ){
			my @tmp = ();
			@tmp    = split (/\s+/,$a_1);
			# Identify empty string
			if (!defined($tmp[2])) {
				print "ERROR Radius of the Sphere empty\n";
				exit;
			} else {			
				$sphereType = $tmp[2];
			}
			#
			$arrays_errors[9] = "radius_sphere";			
		}
		if ( ($a_1=~/maximun_cycles/gi ) ){
			my @tmp = ();
			@tmp    = split (/\s+/,$a_1);
			# Identify empty string
			if (!defined($tmp[2])) {
				print "ERROR Maximun Cycles AUTOMATON empty\n";
				exit;
			} else {			
				$Maximun_Cycles_Automaton = $tmp[2];
			}
			#
			$arrays_errors[10] = "maximun_cycles";			
		}
		if ( ($a_1=~/maximun_energy/gi ) ){
			my @tmp = ();
			@tmp    = split (/\s+/,$a_1);
			# Identify empty string
			if (!defined($tmp[2])) {
				print "ERROR Maximun Energy (kcal/mol) AUTOMATON empty\n";
				exit;
			} else {			
				$MaximunEnergy = $tmp[2];
			}
			#
			$arrays_errors[11] = "maximun_energy";			
		}
		if ( ($a_1=~/clever_automata/gi ) ){
			my @tmp = ();
			@tmp    = split (/\s+/,$a_1);
			# Identify empty string
			if (!defined($tmp[2])) {
				print "ERROR Clever Automata empty\n";
				exit;
			} else {			
				$DeepAU = $tmp[2];
			}
			#
			$arrays_errors[12] = "clever_automata";			
		}		
	}
}
# El numero de atomos debe ser mayor o igual a 3
my $atom_numb_condition = scalar (@Atoms);
if ( $atom_numb_condition < 3 )  { 
	print "\nERROR The number of atoms must be equal to or greater than three\n\n";
	exit(1);
}
#
my @UniqAtoms     = uniq (@Atoms);
my $num_type_atom = scalar (@UniqAtoms);
if ( $num_type_atom == 1) {
	$Geom_Exchange_Option = 0;
} else {
	$Geom_Exchange_Option = 1;
}
# Required
if (!defined($arrays_errors[0]))  { $arrays_errors[0]  = "NO"; }
if (!defined($arrays_errors[1]))  { $arrays_errors[1]  = "NO"; }
if (!defined($arrays_errors[2]))  { $arrays_errors[2]  = "NO"; }
if (!defined($arrays_errors[3]))  { $arrays_errors[3]  = "NO"; }
if (!defined($arrays_errors[4]))  { $arrays_errors[4]  = "NO"; }
if (!defined($arrays_errors[5]))  { $arrays_errors[5]  = "NO"; }
if (!defined($arrays_errors[6]))  { $arrays_errors[6]  = "NO"; }
if (!defined($arrays_errors[7]))  { $arrays_errors[7]  = "NO"; }
if (!defined($arrays_errors[8]))  { $arrays_errors[8]  = "NO"; }
# Optional
if (!defined($arrays_errors[9]))  { $arrays_errors[9]  = "NO"; }
if (!defined($arrays_errors[10])) { $arrays_errors[10] = "NO"; }
if (!defined($arrays_errors[11])) { $arrays_errors[11] = "NO"; }
if (!defined($arrays_errors[12])) { $arrays_errors[12] = "NO"; }
#
my $bolean = errors_config (\@arrays_errors);
if ( $bolean == 0) { exit; }
# Inputs for Gaussian and Mopac
my $option_software;
if (($software=~/gaussian/gi )) {
	$option_software = 0;
	print FILEREP "MESSAGE Choose software Gaussian\n";
} elsif (($software=~/mopac/gi )) {
	$option_software = 1;
	print FILEREP "MESSAGE Choose software Mopac \n";
} elsif (($software=~/lammps/gi )) {
	$option_software = 2;
	print FILEREP "MESSAGE Choose software Lammps\n";	
} else {
	print "ERROR Choose software Gaussian, Mopac or Lammps\n";
	exit (1);
}



###################################
# MAIN
#
# Funcion para el tiempo de ejecucion del programa
my $tiempo_inicial  = new Benchmark;
# automatic box
my @min_coords = ();
my @max_coords = ();
my ($side_plus_x,$side_plus_y,$side_plus_z);
my $side_box = 0;
if ($option_box == 0) {
	# measure side cube
	if ( $Num_of_atoms > 0 ) {	
		my $sides   = automatic_box_length (\@Atoms);;
		$side_box   = sprintf '%.3f',($sides * 2);
	}
	$Box_x = $side_box;
	$Box_y = $side_box;
	$Box_z = $side_box;	
	#
	$side_plus_x  = ($side_box / 2);
	$side_plus_y  = ($side_box / 2);
	$side_plus_z  = ($side_box / 2);		
	#
	my $side_minus = (-1 * $side_plus_x);
	@min_coords    = ($side_minus,$side_minus,$side_minus);
	@max_coords    = ($side_plus_x,$side_plus_y,$side_plus_z);
} else {
	$side_plus_x  = ($Box_x / 2);
	my $side_minus_x = (-1 * $side_plus_x);
	$side_plus_y  = ($Box_y / 2);
	my $side_minus_y = (-1 * $side_plus_y);
	$side_plus_z  = ($Box_z / 2);
	my $side_minus_z = (-1 * $side_plus_z);
	#
	@min_coords    = ($side_minus_x,$side_minus_y,$side_minus_z);
	@max_coords    = ($side_plus_x,$side_plus_y,$side_plus_z);
}
#
my $mi_x = sprintf '%.3f',$min_coords[0];
my $mi_y = sprintf '%.3f',$min_coords[1];
my $mi_z = sprintf '%.3f',$min_coords[2];
#
my $ma_x = sprintf '%.3f',$max_coords[0];
my $ma_y = sprintf '%.3f',$max_coords[1];
my $ma_z = sprintf '%.3f',$max_coords[2];
#
print FILEREP "MESSAGE Box size Min = $mi_x $mi_y $mi_z\n"; 
print FILEREP "MESSAGE Box size Max = $ma_x $ma_y $ma_z\n";
#drawing_box_molecule_VMD (\@min_coords,\@max_coords);
# # # # 
# Genetic operators
my $operator_crossing_over = 0;
if (($crossing_over=~/YES/gi )) {
	$operator_crossing_over = 1;
}
#
my $operator_mutations     = 0;
if (($mutations=~/YES/gi )) {
	$operator_mutations   = 1;
}
#
my @discretized_search_space_1D = ();
my @discretized_search_space_2D = ();
my @discretized_search_space_3D = ();
#
# 1 = YES (type automata) &&  0 = NO (type kick)
my $type_Cellular_automaton = 1;
#
# Automata Inteligente
my $operator_DeepAU = 0;
if (($DeepAU=~/YES/gi )) {
	print FILEREP "MESSAGE Rules of Cellular Automata using Chemistry Intelligence\n";	
	$operator_DeepAU = 1;
}
#
if ( $type_Cellular_automaton == 1 ) {
	# Atomata Cubico
	if( $sphereType == -1) {
		@discretized_search_space_1D = z_axis_path($Box_z);			
		@discretized_search_space_2D = @{Construct_Discrete_Search_Space_Cube ($Box_x,$Box_y,"1.0","0.2",$sphereType)};
		@discretized_search_space_3D = @{Construct_Discrete_Search_Space_Cube ($Box_x,$Box_y,$Box_z,$cell_size_w,$sphereType)};
	} else {
	# Automata Esferico
		$option_1D = 0;
		$option_2D = 0;
		$option_3D = 0;
		@discretized_search_space_3D = @{Construct_Discrete_Search_Space_Cube ($Box_x,$Box_y,$Box_z,$cell_size_w,$sphereType)};
	}
	print FILEREP "MESSAGE Construct discrete search space is done\n";
	print FILEREP "        Grid Size = $cell_size_w\n";
}
# Obtener el numero total de atomos para un formato XYZ
my $NumAtoms    = $Num_of_atoms;
$num_atoms_xyz  = $NumAtoms;
#
my @num_energy_minimun = ();
my $convergence_count  = 0;
my @files_not_converge = ();
my @files_NImag        = ();
my @files_energy       = ();
my @files_coords       = ();
my @files_name         = ();
#
my @array_similar_coords = ();
my @array_similar_energy = ();
my @array_similar_files  = ();
#
my $count_cycle = 1;
while ( $convergence_count < 1 ) {
	my @arrayInputs      = (); 
	my %Info_conver_cell = ();
	#
	my @file_qm          = ();
	my @cartesian_coords = ();
	#
	my @total_D          = ();
	#
	my ($filebase,$string_atoms_coords);
	my ($filebase_2d,$string_atoms_coords_2d);
	my ($filebase_3d,$string_atoms_coords_3d);
	#
	if ($count_cycle == 1) {
		if ( $Num_of_geometries_input > 0 ) {
			print FILEREP "MESSAGE Initial Species (1D, 2D and 3D)\n";			
			if ( $type_Cellular_automaton == 1 ) {
				# Generar poblacion Automata celular
				my $Pob_1D = int ( $Num_of_geometries_input * ( $Percentage_1D / 100) );
				my $Pob_2D = int ( $Num_of_geometries_input * ( $Percentage_2D / 100) );
				my $Pob_3D = int ( $Num_of_geometries_input * ( $Percentage_3D / 100) );
				#
				my $rest_3D = $Num_of_geometries_input - ($Pob_1D + $Pob_2D);
				#
				if( $sphereType == -1) {
					my $Porc_1D = $Pob_1D;				
					($filebase,$string_atoms_coords)       = Opt_Cell_Automata ($Porc_1D,\@Atoms,\@discretized_search_space_1D,$Num_of_atoms,"Cell1D_$count_cycle", 0, 1);
					print FILEREP "        Species 1D = $Porc_1D\n";
					#
					my $Porc_2D = $Pob_2D;
					($filebase_2d,$string_atoms_coords_2d) = Opt_Cell_Automata ($Porc_2D,\@Atoms,\@discretized_search_space_2D,$Num_of_atoms,"Cell2D_$count_cycle", $operator_DeepAU, 2);
					print FILEREP "        Species 2D = $Porc_2D\n";
					#
					my $Porc_3D = $rest_3D;
					($filebase_3d,$string_atoms_coords_3d) = Opt_Cell_Automata ($Porc_3D,\@Atoms, \@discretized_search_space_3D,$Num_of_atoms,"Cell3D_$count_cycle", $operator_DeepAU, 3);
					print FILEREP "        Species 3D = $Porc_3D\n";
				} else {
					($filebase_3d,$string_atoms_coords_3d) = Opt_Cell_Automata ($Num_of_geometries_input,\@Atoms, \@discretized_search_space_3D,$Num_of_atoms,"Cages3D_$count_cycle", 0, 3);
					print FILEREP "        Species 3D Cages = $Num_of_geometries_input\n";					
				}
                #
                #
			} else {			
				# Generar poblacion tipo kick
				#($filebase,$string_atoms_coords) = pob_type_kick ($ma_x,$ma_y,$ma_z,\@Atoms,$Num_of_geometries_input,"Kick3D");
			}
			#
			if( $sphereType == -1) {
				#
				for (my $i = 0 ; $i < scalar (@{$filebase}); $i++) {
					push (@file_qm         ,${$filebase}[$i]);
					push (@cartesian_coords,${$string_atoms_coords}[$i]);
				}
				for (my $i = 0 ; $i < scalar (@{$filebase_2d}); $i++) {
					push (@file_qm         ,${$filebase_2d}[$i]);
					push (@cartesian_coords,${$string_atoms_coords_2d}[$i]);
				}
				for (my $i = 0 ; $i < scalar (@{$filebase_3d}); $i++) {
					push (@file_qm         ,${$filebase_3d}[$i]);
					push (@cartesian_coords,${$string_atoms_coords_3d}[$i]);
				}
				#
			} else {
				for (my $i = 0 ; $i < scalar (@{$filebase_3d}); $i++) {
					push (@file_qm         ,${$filebase_3d}[$i]);
					push (@cartesian_coords,${$string_atoms_coords_3d}[$i]);
				}
			}
		} else {
			my $species_fit   = $file_restart;
			my $new_numb_pob;
			($filebase,$string_atoms_coords,$new_numb_pob) = restart_species ($species_fit);
			print FILEREP "MESSAGE Genetically-fit species\n";
			$Num_of_geometries_input = $new_numb_pob;
			@file_qm          = @{$filebase};
			@cartesian_coords = @{$string_atoms_coords};
		}
		#
		print FILEREP "MESSAGE Input Files:\n";
		for (my $iteration = 0 ; $iteration < $Num_of_geometries_input ; $iteration++) {	
			my $Input_QM;
			if ( $option_software == 0) {
				$Input_QM = G03Input ($file_qm[$iteration],$header,$ncpus,$mem,$Charge,$Multiplicity,$cartesian_coords[$iteration],$file_qm[$iteration]);			
			}  elsif ( $option_software == 1) {
				$Input_QM = MopacInput ($file_qm[$iteration],$cartesian_coords[$iteration],$file_qm[$iteration],$header,$Charge,$Multiplicity,$ncpus,$mem);
			} else {		
				my $tmp_box = "$mi_x $mi_y $mi_z $ma_x $ma_y $ma_z";
				$Input_QM = LammpsInput ($file_qm[$iteration],$cartesian_coords[$iteration],$file_qm[$iteration],$header,$tmp_box);
			}
			print FILEREP "        $Input_QM\n";
			push (@arrayInputs,$Input_QM);
		}
        #
        CycleCoords (\@cartesian_coords,\@file_qm, $Num_of_atoms,"PreCoords-$count_cycle" );     
        #
	} else {
		# # # # # # # # # # # # #
		# Ver duplicados 
		my $Num_of_geometries = scalar (@files_energy);
		print FILEREP "MESSAGE Find Similar 3-dimensional Structures\n";
		print FILEREP "        Threshold Duplicates = $threshold_duplicate\n";
		#
		my ($not_struc_dup,$energy_dup,$files_dup,$coords_sim,$energy_sim,$files_sim) =  
		info_duplicate_structures ($Num_of_geometries,$Num_of_atoms,\@files_coords,\@files_energy,\@files_name,$nprocess, $threshold_duplicate);		
		# Guardar Duplicados
		my @arr_1 = @{$coords_sim};
		my @arr_2 = @{$energy_sim};
		my @arr_3 = @{$files_sim};
		for (my $dup = 0 ; $dup < scalar (@arr_1); $dup++) {
			push (@array_similar_coords,$arr_1[$dup]);
			push (@array_similar_energy,$arr_2[$dup]);
			push (@array_similar_files ,$arr_3[$dup]);
		}
		# # # # # # # # # # # # # 
		# Evaluar funcion fitness
		my @fitness_values = ();
		my $min_energy     = min @{$energy_dup};
		my $max_energy     = max @{$energy_dup};
		foreach my $energy ( @{$energy_dup} ) {
			my $prob              = ($energy - $min_energy)/( $max_energy - $min_energy );
			my $alpha             = 3;
			my $fitness_functions = exp (-($alpha*$prob));
			#my $fitness_functions = 0.5*( 1 - tanh ((2*$prob) - 1));
			#my $fitness_functions = 1 - (0.7*$prob);
			push (@fitness_values,$fitness_functions);
		}
		my @fitness_ener_sort = ();
		my @val_ener_sort     = ();		
		my @coord_xyz_sort    = ();
		my @name_sort         = ();		
		my @idx = sort { $fitness_values[$a] <=> $fitness_values[$b] } 0 .. $#fitness_values;
		@fitness_ener_sort = @fitness_values[@idx];
		@val_ener_sort     = @{$energy_dup}[@idx];
		@coord_xyz_sort    = @{$not_struc_dup}[@idx];
		@name_sort         = @{$files_dup}[@idx];
		#
		my $string_axis;
		my $array_axis;	
		my @new_elements_atoms = ();
		my @new_id             = ();
		#
		my @the_best_coords_fitness = ();
		#
		for (my $i = 0 ; $i < scalar (@fitness_ener_sort); $i++) {
			if ( $fitness_ener_sort[$i] >= 0.5 ) {
				push (@the_best_coords_fitness,$coord_xyz_sort[$i]);
			}
		}
		# # # # # # # # # # # # #
		# Crossing Over
		# Recombinacion del material genetico de los individuos	
		if ($operator_crossing_over == 1) {
			print FILEREP "MESSAGE Crossing Over\n";
			print FILEREP "        Spacing Between Planes = $spacingInPlane\n";
			
			my @Crossing_coords     = ();
			my @Crossing_id         = ();
			#
			my $numberOfRec    = scalar (@the_best_coords_fitness);
			my @tmp_array_pair = Crossing_Over ($numberOfRec,\@the_best_coords_fitness);
			# Se escogen de forma aleatoria los hijos, para pasar a la proxima 	
			# generacion	
			my @shuffle        = fisher_yates_shuffle (\@tmp_array_pair);
			for (my $i = 0 ; $i < scalar (@shuffle) ; $i++) {
				my 	$string_axis_pair = ($shuffle[$i]);
				push (@Crossing_coords,$string_axis_pair);
				my $id               = sprintf '%.6d', $count_struc_crossing;
				my $filebase         = "Child_$count_cycle\_$id";
				push (@Crossing_id,$filebase);
				$count_struc_crossing++;
			}
			#
			if ( scalar (@tmp_array_pair) < $Num_of_geometries_input) {
				for (my $i = 0 ; $i < scalar (@tmp_array_pair) ; $i++) {
					push (@new_elements_atoms,$Crossing_coords[$i]);
					push (@new_id            ,$Crossing_id[$i]);
				}
			} else {
				for (my $i = 0 ; $i < $Num_of_geometries_input; $i++) {
					push (@new_elements_atoms,$Crossing_coords[$i]);
					push (@new_id            ,$Crossing_id[$i]);
				}
			}
		}
		# # # # # # # # # # # # #
		# Mutaciones
		# Creacion de la nueva poblacion mediante mutaciones o movimientos de los atomos
		my $fix_numb_new_cell = int ( $New_Pobl_Mut * $Num_of_geometries_input);
		if ($operator_mutations    == 1) {
			print FILEREP "MESSAGE Mutations\n";
			print FILEREP "        Kick Factor = $kick_punch\n";
			my @Mutants_coords     = ();
			my @Mutants_id         = ();
			#
			my $range = scalar (@the_best_coords_fitness);			
			for (my $i = 0 ; $i < $fix_numb_new_cell; $i++) {
				my $random_number = int(rand($range));
				my $boolean = 0;
				while ( $boolean < 1 ) {
					($array_axis,$string_axis) = Mutations_Automata ($the_best_coords_fitness[$random_number],$Num_of_atoms);
					my $option_boo             = steric_impediment_atoms (\@{$array_axis},1);
					if ( $option_boo == 0 ) { $boolean = 1;}
				}
				push (@Mutants_coords,$string_axis);
				my $id               = sprintf '%.6d', $count_struc_mut;
				my $filebase         = "Mut_$count_cycle\_$id";
				push (@Mutants_id,$filebase);
				$count_struc_mut++;
			}
			#
			for (my $i = 0 ; $i < scalar (@Mutants_id) ; $i++) {
				push (@new_elements_atoms,$Mutants_coords[$i]);
				push (@new_id            ,$Mutants_id[$i]);
			}
			# # # # # # # # # # # # # # # #
			# Geometric exchange operator
			my @MutExtCoord = ();
			my @MutExtId    = ();
			#
			my @tmp_FY      = split (/\n/,$the_best_coords_fitness[0]);
			if ( $Geom_Exchange_Option == 1) {
				for (my $i = 0; $i < $fix_numb_new_cell; $i++) {
					my @recombination_FY = fisher_yates_shuffle (\@tmp_FY);
					my @cooTMP           = split (/\n/,$the_best_coords_fitness[0]);
					my $formatxyz;
					for (my $i = 0; $i < $Num_of_atoms; $i++) {
						my($element_best, $cx_best, $cy_best, $cz_best) = split (/\s+/,$cooTMP[$i]);
						my($element, $cx, $cy, $cz)                     = split (/\s+/,$recombination_FY[$i]);
						$formatxyz.= "$element  $cx_best  $cy_best  $cz_best\n";
					}
					my $id               = sprintf '%.6d', $count_struc_mut_Ex;
					my $filebase         = "MutEx_$count_cycle\_$id";
					#
					push (@MutExtCoord,$formatxyz);
					push (@MutExtId   ,$filebase);
					$count_struc_mut_Ex++;
				}
				#
				for (my $i = 0 ; $i < scalar (@MutExtId) ; $i++) {
					push (@new_elements_atoms,$MutExtCoord[$i]);
					push (@new_id            ,$MutExtId[$i]);
				}
			}
		}
		# # # # # # # # # # # # # #
		# La nueva poblacion generada por el automata celular
		print FILEREP "MESSAGE New species (1D, 2D and 3D)\n";
		my $numb_1D_2D      = int ( $New_Pobl_Cell * $Num_of_geometries_input);
		my @new_name_base   = ();
		my @new_coords_base = ();
		if ( $type_Cellular_automaton == 1 ) {
			my $new_numb = int ($numb_1D_2D/2);
			if ( $option_1D == 1) {	
				# Generar poblacion Automata celular 1D
				($filebase,$string_atoms_coords) = Opt_Cell_Automata ($new_numb,\@Atoms, \@discretized_search_space_1D,$Num_of_atoms,"Cell1D_$count_cycle", 0, 1);
			
				for (my $i = 0 ; $i < scalar (@{$filebase}) ; $i++) {
					push (@new_name_base  ,$$filebase[$i]);
					push (@new_coords_base,$$string_atoms_coords[$i]);
				}
				#
				print FILEREP "        New species 1D Build\n";
			}
			if ( $option_2D == 1) {	
				$filebase = ();
				$string_atoms_coords = ();	
				# Generar poblacion Automata celular 2D
				($filebase,$string_atoms_coords) = Opt_Cell_Automata ($new_numb,\@Atoms,\@discretized_search_space_2D,$Num_of_atoms,"Cell2D_$count_cycle",$operator_DeepAU, 2);

				for (my $i = 0 ; $i < scalar (@{$filebase}) ; $i++) {
					push (@new_name_base  ,$$filebase[$i]);
					push (@new_coords_base,$$string_atoms_coords[$i]);
				}
				#
				print FILEREP "        New species 2D Build\n";
			}
			if ( $option_3D == 1) {	
				$filebase = ();
				$string_atoms_coords = ();			
				# Generar poblacion Automata celular 3D
				($filebase,$string_atoms_coords) = Opt_Cell_Automata ($new_numb,\@Atoms,\@discretized_search_space_3D,$Num_of_atoms,"Cell3D_$count_cycle",$operator_DeepAU, 3);

				for (my $i = 0 ; $i < scalar (@{$filebase}) ; $i++) {
					push (@new_name_base  ,$$filebase[$i]);
					push (@new_coords_base,$$string_atoms_coords[$i]);
				}
				#
				print FILEREP "        New species 3D Build\n";
			}	
		} else {
			# Generar poblacion tipo kick
			# Generar poblacion 1D
			# Generar poblacion 2D			
			# Generar poblacion 3D
		}
		#
		my @name_new_file  = @new_name_base; 
		my @name_new_atoms = @new_coords_base;
		for (my $i = 0 ; $i < scalar (@name_new_atoms); $i++) {
			push (@new_elements_atoms,$name_new_atoms[$i]);
			push (@new_id,$name_new_file[$i]);
		}
		# # # # # # # # # # # # # # # # # 
		# Enviar a hacer calculos cuanticos
		my @file_qm          = @new_id;
		my @cartesian_coords = @new_elements_atoms;
        #
        CycleCoords (\@cartesian_coords,\@file_qm,$Num_of_atoms,"PreCoords-$count_cycle" );
        #
		print FILEREP "MESSAGE Input Files:\n";
		for (my $iteration = 0 ; $iteration < scalar(@file_qm) ; $iteration++) {	
			my $Input_QM;
			if ( $option_software == 0) {
				$Input_QM = G03Input ($file_qm[$iteration],$header,$ncpus,$mem,$Charge,$Multiplicity,$cartesian_coords[$iteration],$file_qm[$iteration]);
			} elsif ( $option_software == 1) {
				$Input_QM = MopacInput ($file_qm[$iteration],$cartesian_coords[$iteration],$file_qm[$iteration],$header,$Charge,$Multiplicity,$ncpus,$mem);
			} else {
				my $tmp_box = "$mi_x $mi_y $mi_z $ma_x $ma_y $ma_z";
				$Input_QM = LammpsInput ($file_qm[$iteration],$cartesian_coords[$iteration],$file_qm[$iteration],$header,$tmp_box);	
			}
			#
			print FILEREP "        $Input_QM\n";
			push (@arrayInputs,$Input_QM);
		}
	}
	#
	my @opt_name      = ();
	my @value_energy  = ();
	my @value_coords  = ();
	# Estructura con errores
	my @rep_structure = ();
	# Estructura con frecuencias imaginarias
	my @rep_NImag	  = ();
	if ( $option_software == 0) {
		print FILEREP "MESSAGE Send Gaussian\n";
		my ($value_energy_sort,$value_coords_sort,
			$report_structure,$report_NImag,$file_opt) = submit_queue_gaussian (\@arrayInputs,$num_atoms_xyz,$header,$ncpus,$mem,$Charge,$Multiplicity,$queue,$exec_bin_g09,$option_software);
	
		#
		print FILEREP "MESSAGE End Gaussian\n";
		#
		my $sumPort = $Global_Port_Normal + $Global_Port_Error + $Global_Port_Corrupt + $Global_No_Coords;
		print FILEREP "        Normal Termination Number        = $Global_Port_Normal/$sumPort\n";
		print FILEREP "        Convergence Failure Number       = $Global_Port_Error/$sumPort\n";		
		print FILEREP "        Corrupted Files Number           = $Global_Port_Corrupt/$sumPort\n";
		print FILEREP "        Files Without Coordinates Number = $Global_No_Coords/$sumPort\n";
		#
		print FILEREP "        * Convergence Failure Files: \n";
		foreach my $outx (@files_Port_Error) {
			print FILEREP "        $outx.log\n";
		}
		#
		print FILEREP "        * Files Without Coordinates: \n";
		foreach my $outx (@files_No_Coords) {
			print FILEREP "        $outx.log\n";
		}		
		#
		@opt_name      = @{$file_opt};
		@value_energy  = @{$value_energy_sort};
		@value_coords  = @{$value_coords_sort};
		# Estructura con errores
		@rep_structure = @{$report_structure};
		# Estructura con frecuencias imaginarias
		@rep_NImag	   = @{$report_NImag};
		#
		for ( my $i=0; $i < scalar(@arrayInputs); $i++) {
			(my $without_extension = $arrayInputs[$i]) =~ s/\.[^.]+$//;
			unlink ("$without_extension.com");
			unlink ("$without_extension.log");
			unlink ("$without_extension.out");
		}
	} elsif ( $option_software == 1) {
		print FILEREP "MESSAGE Send MOPAC\n";
		submit_queue_mopac (\@arrayInputs,$path_bin_mopac);        
        #submit_Mopac_Lammps (\@arrayInputs,$path_bin_mopac, $option_software);
		print FILEREP "MESSAGE End MOPAC\n";
		my ($value_coords_sort,$file_opt,$value_energy_sort) = energy_mopac ($num_atoms_xyz,\@arrayInputs);
		@opt_name      = @{$file_opt};
		@value_energy  = @{$value_energy_sort};
		@value_coords  = @{$value_coords_sort};
		for ( my $i=0; $i < scalar(@arrayInputs); $i++) {
			(my $without_extension = $arrayInputs[$i]) =~ s/\.[^.]+$//;
			unlink ("$without_extension.mop");
			unlink ("$without_extension.aux");
			unlink ("$without_extension.out");	
			unlink ("$without_extension.arc");			
		}
		unlink ("tmp_mopac_1.txt");
		unlink ("tmp_mopac_2.txt");
	} else {
		print FILEREP "MESSAGE Send LAMMPS\n";
		#my @Outputs_lammps = submit_queue_lammps (\@arrayInputs,$path_bin_lammps);
		my @Outputs_lammps = submit_Mopac_Lammps (\@arrayInputs,$path_bin_lammps,$option_software);
		print FILEREP "MESSAGE End LAMMPS\n";
		my ($value_coords_sort,$value_energy_sort) = coords_energy_lammps ($num_atoms_xyz,\@Outputs_lammps);
		@opt_name      = @Outputs_lammps;
		@value_energy  = @{$value_energy_sort};
		@value_coords  = @{$value_coords_sort};		
		for ( my $i=0; $i < scalar(@Outputs_lammps); $i++) {
			(my $without_extension = $Outputs_lammps[$i]) =~ s/\.[^.]+$//;
			unlink ("$without_extension.dat");
			unlink ("$without_extension.out");	
			unlink ("$without_extension.xyz");
			unlink ("$without_extension.in");
		}	
	}
    #
    #
    CycleCoords (\@value_coords,\@opt_name,$Num_of_atoms,"PostCoords-$count_cycle" );
    #
	#
	foreach (@opt_name)      { push (@files_name,        $_);}	
	foreach (@value_energy)  { push (@files_energy,      $_);}
	foreach (@value_coords)  { push (@files_coords,      $_);}
	if ( $option_software == 0) {
		# Estructura con errores
		foreach (@rep_structure) { push (@files_not_converge,$_);}
		# Estructura con frecuencias imaginarias
		foreach (@rep_NImag)     { push (@files_NImag,   $_);}	
	}
	#
	#
	#
	#
	# En esta parte se ve el numero de veces que se repite 
	# la energia para que termine el ciclo while
	#my @val_ener_sort  = ();
	#my @coord_xyz_sort = ();
	#my @name_sort      = ();
	#my @idx = sort { $files_energy[$a] <=> $files_energy[$b] } 0 .. $#files_energy;
	#@val_ener_sort  = @files_energy[@idx];
	#@coord_xyz_sort = @files_coords[@idx];
	#@name_sort      = @files_name[@idx];
    #
    #
    my ($files_energy_tmp,$files_coords_tmp,$files_name_tmp ) = EnergyFilter (\@files_energy,\@files_coords,\@files_name, $option_software, $MaximunEnergy);
    #
    my @val_ener_tmp  = @{$files_energy_tmp};
	my @coord_xyz_tmp = @{$files_coords_tmp};
	my @name_tmp      = @{$files_name_tmp};
    #
    my @val_ener_sort  = ();
	my @coord_xyz_sort = ();
	my @name_sort      = ();
    my @idx = sort { $val_ener_tmp[$a] <=> $val_ener_tmp[$b] } 0 .. $#val_ener_tmp;
	@val_ener_sort  = @val_ener_tmp[@idx];
	@coord_xyz_sort = @coord_xyz_tmp[@idx];
	@name_sort      = @name_tmp[@idx];
    #
    #
	my $minimun_energy = sprintf '%.5f', $val_ener_sort[0];
    #
	push (@num_energy_minimun,$minimun_energy); 
	foreach my $element( @num_energy_minimun ) {
		++$Info_conver_cell{$element};
	}
	print FILEREP "\n";
	print FILEREP "Energy = $minimun_energy H\n";
	print FILEREP "* Lead Species Coords $name_sort[0]\n";	
	print FILEREP "$coord_xyz_sort[0]";
	print FILEREP "\n";
	# File for restart
	# Imprimir duplicados y no duplicados archivos
	#my $Num_of_geometries_output = scalar (@files_energy);
    #
    #
    my $Num_of_geometries_output = scalar (@val_ener_sort);
    #
    #
    #my ($not_struc_dup,$energy_dup,$files_dup,$coords_sim,$energy_sim,$files_sim) = info_duplicate_structures ($Num_of_geometries_output,$num_atoms_xyz,\@files_coords,\@files_energy,\@files_name,$nprocess, $threshold_duplicate);
    #
    #
	my ($not_struc_dup,$energy_dup,$files_dup,$coords_sim,$energy_sim,$files_sim) = info_duplicate_structures ($Num_of_geometries_output,$num_atoms_xyz,\@coord_xyz_sort,\@val_ener_sort,\@name_sort,$nprocess, $threshold_duplicate);
	#
    #
	# Duplicados
	final_opt_coords (\@{$energy_sim},\@{$coords_sim},\@{$files_sim},$num_atoms_xyz,$option_software,"02Duplicate_coords.xyz");
    #
    #
	# No duplicados
	final_opt_coords (\@{$energy_dup},\@{$not_struc_dup},\@{$files_dup},$num_atoms_xyz,$option_software,"01Final_coords.xyz");	
	print FILEREP "MESSAGE Cycle $count_cycle Completed\n\n\n";
    #
    #
#	final_opt_coords (\@files_energy,\@files_coords,\@files_name,$num_atoms_xyz,$option_software,"03Restart_coords.xyz");
	if ( $Info_conver_cell{$minimun_energy} > $max_numb_convergence ) {
		# Finish Cell algorith
		print FILEREP "\n\nNormal Termination AUTOMATON\n\n";
		print FILEREP "Final-Energy = $minimun_energy H\n";
		print FILEREP "* Final Lead Species Coords $name_sort[0]\n";
		print FILEREP "$coord_xyz_sort[0]";
		print FILEREP "\n";		
		$convergence_count = 2;
	}
	$count_cycle++;
    #
    # Maximun cycles
    if ($count_cycle > $Maximun_Cycles_Automaton ) {
        my $tiempo_final  = new Benchmark;
        my $tiempo_total  = timediff($tiempo_final, $tiempo_inicial);
        print FILEREP "\n\tExecution Time: ",timestr($tiempo_total),"\n";
        print FILEREP "\n";
        my $datesFinal = localtime();
        print FILEREP "              $datesFinal\n\n";
        #
        close (FILEREP);
        exit (1);
    }
}
#
#
my $tiempo_final  = new Benchmark;
my $tiempo_total  = timediff($tiempo_final, $tiempo_inicial);
print FILEREP "\n\tExecution Time: ",timestr($tiempo_total),"\n";
print FILEREP "\n";
my $datesFinal = localtime();
print FILEREP "              $datesFinal\n\n";
#
close (FILEREP);
#
#
unlink ("log.cite");
unlink ("log.lammps");
exit(0);
