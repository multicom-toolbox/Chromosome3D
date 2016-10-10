#!/usr/bin/perl -w
# version 1.2, Badri Adhikari, 2/18/2016

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

# Location of CNSsuite
my $cns_suite      = "/home/bap54/chr3D_projects/project01/programs/cns_solve_1.3";
my $cns_executable = "$cns_suite/intel-x86_64bit-linux/bin/cns_solve";
my $dg_sa_log      = "/dev/null"; # replace with dgsa.log if space is NOT a concern

# Chromosome Model building parameters
my $KSCALING   = 11;
my $ALPHA      = 0.5;
my $SEPARATION = 5;
my $MODELCOUNT = 20;

# User inputs 
my $help;
my $dir_out;
my $file_if;

GetOptions(
	"h"		=> \$help,
	"o=s"	=> \$dir_out,
	"k=i"	=> \$KSCALING,
	"a=s"	=> \$ALPHA,
	"m=i"	=> \$MODELCOUNT,
	"i=s"	=> \$file_if)
or confess "ERROR! Error in command line arguments!";

print_usage() if defined $help;
print_usage("Input IF matrix not found!") if not defined $file_if;
print_usage("Output directory not defined!") if not defined $dir_out;
print_usage("Input IF file $file_if does not exist!") if not -f $file_if;
$KSCALING = 11 if not defined $KSCALING;
$ALPHA = 0.5 if not defined $ALPHA;
$MODELCOUNT = 20 if not defined $MODELCOUNT;

mkdir $dir_out or confess "ERROR! Could not create output directory $dir_out!" if not -d $dir_out;
print_usage("Output directory $dir_out does not exist!") if not -d $dir_out;

print "Start Time : ".(localtime)." [$0]\n";
print "Input      : $file_if\n";
print "Output Dir : $dir_out\n";
print "Scaling(K) : $KSCALING\n";
print "Alpha      : $ALPHA\n";
print "Effective Conversion Equation is : D = $KSCALING/[(IF^$ALPHA)*(avg of IF^$ALPHA)]\n";

my $ID = basename($file_if, ".txt");
system_cmd("rm -f $dir_out/*");
system_cmd("cp $file_if $dir_out/") if not -f $ID;
chdir $dir_out or confess $!;
$file_if = $ID.".txt";
$dir_out = ".";
my $L = calc_len_IF($file_if);
print "L          : $L\n";

# CNS Suite parameters
my $min_sep     = $SEPARATION;
my $con_wt      = 10;
my $model_count = $MODELCOUNT;
my $rep2        = 0.85;
my $atomselect  = 2;
my $rep1        = 1.0;
my $mini        = 15000;
my $model_info_log =  "model_info.log";
my $mode        = "trial";
my $DISTRELAX   = 0.5;

# http://proteopedia.org/wiki/index.php/Standard_Residues and http://prody.csb.pitt.edu/manual/reference/atomic/flags.html
our %AA3TO1 = qw(ALA A ASN N CYS C GLN Q HIS H LEU L MET M PRO P THR T TYR Y ARG R ASP D GLU E GLY G ILE I LYS K PHE F SER S TRP W VAL V);
our %AA1TO3 = reverse %AA3TO1;

# Verify CNS Suite installation
confess "ERROR! Cannot find CNS suite folder at $cns_suite! Check CNS installation!" if not -d $cns_suite;
confess "ERROR! cns_solve_env.sh not found inside $cns_suite/! Check CNS installation!" if not -f "$cns_suite/cns_solve_env.sh";
confess "ERROR! CNS executable $cns_executable not found! Check CNS installation!" if not -f $cns_executable;
confess "ERROR! protein.param not found inside $cns_suite/libraries/toppar/! Check CNS installation!" if not -f "$cns_suite/libraries/toppar/protein.param";

# Generate restraints from IF matrix
IF2dist_new("$ID.txt", "$ID.dist", $KSCALING);
dist2rr("$ID.dist", "$ID.rr");
carr2tbl("$ID.rr", "contact.tbl");
print "Restraints : ".count_lines("contact.tbl")." lines in tbl file\n";

# Make Fasta file for CNS Suite input
my $REFSEQUENCE = "RSEDWQCPRTPYAASRDFDVKYVVPSFSAGGLVQAMVTYEGDRNESAVFVAIRNRLHVLGPDLKSVQSLATGPAGDPGCQTCAACGPGPHGPPGDTDTKVLVLDPALPALVSCGSSLQGRCFLHDLEPQGTAVHLAAPACLFSAHHNRPDDCPDCVASPLGTRVTVVEQGQASYFYVASSLDAAVAGSFSPRSVSIRRLKADASGFAPGFVALSVLPKHLVSYSIEYVHSFHTGAFVYFLTVQPASVTDDPSALHTRLARLSATEPELGDYRELVLDCRFAPKLVPRGSPEGGQPYPVLQVAHSAPVGAQLATELSIAEGQEVLFGVFVTGKDGGPGVGPNSVVCAFPIDLLDTLIDEGVERCCESPVHPGLRRGLDFFQSPSFCPNPPGLEALSPNTSCRHFPLLVSSSFSRVDLFNGLLGPVQVTALYVTRLDNVTVAHMGTMDGRILQVELVRSLNYLLYVSNFSLGDSGQPVQRDVSRLGDHLLFASGDQVFQVPIQGPGCRHFLTCGRCLRAWHFMGCGWCGNMCGQQKECPGSWQQDHCPPKLTEFHPHSGPLRGSTRLTLCGSNFYLHPSGLVPEGTHQVTVGQSPCRPLPKDSSKLRPVPRKDFVEEFECELEPLGTQAVGPTNVSLTVTNMPPGKHFRVDGTSVLRGFSFMETG";
my $sequence = substr $REFSEQUENCE, 0, $L;
my $file_seq = "$ID.fasta";
system_cmd("rm -f $file_seq");
print2file($file_seq, ">$ID");
print2file($file_seq, "$sequence");

print "\n";
print "(A) Build extended mtf and pdb..\n";
build_extended();
print "(B) Build models using CNS dgsa..\n";
build_models();
print "(C) Assess models..\n";
assess_dgsa();

print "\nFinished [$0]: ".(localtime)."\n";

sub IF2dist_new{
	my $matrix = shift;
	my $file_dist = shift;
	my $a = shift;
	system_cmd("rm -f $file_dist");
	my %input = ();
	open MATRIX, $matrix or confess $!;
	my $i = 0;
	while(<MATRIX>){
		chomp $_;
		$_ =~ s/^\s+//;
		confess "??" if length($_) < 3;
		my @R = split /\s+/, $_;
		for(my $j = 0; $j <= $#R; $j++){
			confess "ERROR! dist not defined!" if not defined $R[$j];
			$input{$i+1}{$j+1} = $R[$j];
		}
		$i++;
	}
	close MATRIX;
	my %output = ();
	my $mean_sqrt = 0;
	for(my $i = 1; $i <= $L; $i++){
		for(my $j = 1; $j <= $L; $j++){
			$output{$i}{$j} = $input{$i}{$j} ** $ALPHA;
			$mean_sqrt += $output{$i}{$j};
		}
	}
	# divide by mean
	$mean_sqrt = $mean_sqrt/($L * $L);
	for(my $i = 1; $i <= $L; $i++){
		for(my $j = 1; $j <= $L; $j++){
			$output{$i}{$j} = $output{$i}{$j}/$mean_sqrt;
		}
	}
	# perform conversion
	for(my $i = 1; $i <= $L; $i++){
		for(my $j = 1; $j <= $L; $j++){
			if ($output{$i}{$j} == 0){
				$output{$i}{$j} = -1;
				next;
			}
			$output{$i}{$j} = $a / $output{$i}{$j};
		}
	}
	# write output
	for(my $i = 1; $i <= $L; $i++){
		for(my $j = 1; $j <= $L; $j++){
			print2line($file_dist, sprintf "%.1f ", $output{$i}{$j});
		}
		print2file($file_dist, "");
	}
}

sub calc_len_IF{
	my $matrix = shift;
	my $length;
	open MATRIX, $matrix or confess $!;
	while(<MATRIX>){
		chomp $_;
		$_ =~ s/^\s+//;
		die "??" if length($_) < 2;
		my @R = split /\s+/, $_;
		$length = $#R + 1;
		last;
	}
	close MATRIX;
	confess ":(" if not defined $length;
	return $length;
}

sub dist2rr{
	my $matrix = shift;
	my $rr = shift;
	my $SEP  = shift;
	confess "NEED matrix and rr" if not ($matrix && $rr); 
	my %distances = ();
	open MATRIX, $matrix or confess $!;
	my $i = 0;
	while(<MATRIX>){
		chomp $_;
		$_ =~ s/^\s+//;
		die "??" if length($_) < 2;
		my @R = split /\s+/, $_;
		for(my $j = $i+1; $j <= $#R; $j++){
			next if abs($j - $i) < $min_sep;
			next if $R[$j] <= 0;
			$distances{($i+1)." ".($j+1)} = $R[$j];
		}
		$i++;
	}
	close MATRIX;
	system_cmd("rm -f $rr");
	foreach (sort keys %distances){
		print2file($rr, (sprintf "$_ %.2f %.2f 1.0", $distances{$_}, $distances{$_}));
	}
}

sub add_connect_rows{
	my $in_pdb = shift;
	my $L = length(seq_chain($in_pdb));
	for(my $i = 1; $i < $L; $i++){
		print2file($in_pdb, sprintf "CONECT%5s%5s", $i, $i+1);
	}
	print2file($in_pdb, "END");
}

sub seq_chain{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my $seq = "";
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		next unless (parse_pdb_row($_,"aname") eq "CA");
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		my $res = $AA3TO1{parse_pdb_row($_,"rname")};
		$seq .= $res;
	}
	close CHAIN;
	confess "$chain has less than 1 residue!" if (length($seq) < 1);
	return $seq;
}

sub build_extended{
	confess "ERROR! expected fasta file not found in the extended directory!" if not -f $file_seq;
	write_cns_seq($file_seq, "input.seq");
	write_cns_generate_seq_file();
	write_cns_generate_extended_file();
	open  JOB, ">job.sh" or confess "ERROR! Could not open job.sh $!";
	print JOB "#!/bin/bash                                     \n";
	print JOB "# CNS-CONFIGURATION                             \n";
	print JOB "source $cns_suite/cns_solve_env.sh        \n";
	print JOB "export KMP_AFFINITY=none                        \n";
	print JOB "$cns_executable < gseq.inp > gseq.log \n";
	print JOB "$cns_executable < extn.inp > extn.log \n";
	close JOB;
	system_cmd("chmod +x job.sh");
	system_cmd("./job.sh", "job.log");
	confess "FAILED! extended.mtf not found!" if not -f "extended.pdb";
	system_cmd("rm -f gseq.log");
	system_cmd("rm -f extn.log");
}

sub build_models{
	# prepare CNS task files
	system_cmd("rm -f iam.*");
	write_cns_dgsa_file();
	open  JOB, ">job.sh" or confess "ERROR! Could not open job.sh $!";
	print JOB "#!/bin/bash                                       \n";
	print JOB "echo \"starting cns..\"                           \n";
	print JOB "touch iam.running                                 \n";
	print JOB "# CNS-CONFIGURATION                               \n";
	print JOB "source $cns_suite/cns_solve_env.sh                \n";
	print JOB "export KMP_AFFINITY=none                          \n";
	print JOB "$cns_suite/intel-x86_64bit-linux/bin/cns_solve < dgsa.inp > $dg_sa_log \n";
	print JOB "if [ -f \"${ID}_${model_count}.pdb\" ]; then      \n";
	print JOB "   rm iam.running                                 \n";
	print JOB "   echo \"trial structures written.\"             \n";
	print JOB "   rm *embed*                                     \n";
	print JOB "   exit                                           \n";
	print JOB "fi                                                \n";
	print JOB "if [ -f \"${ID}a_${model_count}.pdb\" ]; then \n";
	print JOB "   rm iam.running                                 \n";
	print JOB "   echo \"accepted structures written.\"          \n";
	print JOB "   rm *embed*                                     \n";
	print JOB "   exit                                           \n";
	print JOB "fi                                                \n";
	if ($dg_sa_log ne "/dev/null"){
		print JOB "tail -n 30 $dg_sa_log                         \n";
	}
	print JOB "echo \"ERROR! Final structures not found!\"       \n";
	print JOB "echo \"CNS FAILED!\"                              \n";
	print JOB "mv iam.running iam.failed                         \n";
	close JOB;
	print "Starting job [$dir_out/job.sh > job.log]\n";
	system_cmd("chmod +x job.sh");
	system_cmd("./job.sh > job.log");
	confess "ERROR! Something went wrong while running CNS! Check job.log and dgsa.log!" if -f "iam.failed";
}

sub system_cmd{
	my $command = shift;
	my $log = shift;
	confess "EXECUTE [$command]?\n" if (length($command) < 5  and $command =~ m/^rm/);
	if(defined $log){
		system("$command &> $log");
	}
	else{
		system($command);
	}
	if($? != 0){
		my $exit_code  = $? >> 8;
		confess "ERROR!! Could not execute [$command]! \nError message: [$!]";
	}
}

sub count_lines{
	my $file = shift;
	my $lines = 0;
	return 0 if not -f $file;
	open FILE, $file or confess "ERROR! Could not open $file! $!";
	while (<FILE>){
		chomp $_;
		$_ =~ tr/\r//d; # chomp does not remove \r
		next if not defined $_;
		next if length($_) < 1;
		$lines ++;
	}
	close FILE;
	return $lines;
}

sub print2file{
	my $file = shift;
	my $message = shift;
	my $newline = shift;
	$newline = "\n" if not defined $newline;
	if (-f $file){
		open  FILE, ">>$file" or confess $!;
		print FILE $message.$newline;
		close FILE;
	}
	else{
		open  FILE, ">$file" or confess $!;
		print FILE $message.$newline;
		close FILE;
	}
}

sub carr2tbl{
	my $file_rr  = shift;
	my $file_tbl = shift;
	confess ":(" if not -f $file_rr;
	open RR, $file_rr or confess $!;
	while(<RR>){
		next unless ($_ =~ /[0-9]/);
		$_ =~ s/^\s+//;
		next unless ($_ =~ /^[0-9]/);
		my @C = split(/\s+/, $_);
		confess "ERROR! Expecting a pair in row [".$_."]!\n" if (not defined $C[0] || not defined $C[1]);
		confess "ERROR! Confidence column not defined in row [".$_."] in file $file_rr!\n" if not defined $C[4];
		my $distance = sprintf("%.2f", ($C[3]+$C[2])/2);
		my $negdev   = sprintf("%.2f", ($C[3]-$C[2])/2);
		my $posdev   = $negdev;
		if ($C[2] eq "0"){
			$distance = sprintf("%.2f", 3.6);
			$negdev   = sprintf("%.2f", 0.1);
			$posdev   = sprintf("%.2f", ($C[3] - 3.6));
		}
		print2file($file_tbl, (sprintf "assign45 (resid %3d and name %2s) (resid %3d and name %2s) %.2f %.2f %.2f", $C[0], "ca", $C[1], "ca", $distance, $negdev, $posdev));
	}
}

sub print2line{
	my $file = shift;
	my $message = shift;
	if (-f $file){
		open  FILE, ">>$file" or confess $!;
		print FILE $message;
		close FILE;
	}
	else{
		open  FILE, ">$file" or confess $!;
		print FILE $message;
		close FILE;
	}
}

sub tbl2rows_hash{
	my $file_tbl = shift;
	confess "ERROR! tbl file $file_tbl does not exist!" if not -f $file_tbl;
	my %noe = ();
	open NOETBL, $file_tbl or confess $!;
	while (<NOETBL>){
		$_ =~ s/^\s+//;
		next if $_ !~ m/^assign/;
		$_ =~ s/^\s+//;
		$_ =~ s/\(/ /g;
		$_ =~ s/\)/ /g;
		$noe{$_} = 1;
	}
	close NOETBL;
	confess "$file_tbl seems empty!" if not %noe;
	return %noe;
}

sub coverage_tbl{
	my $file_fasta = shift;
	my $file_tbl = shift;
	my $flag_dihe = shift;
	confess "ERROR! fasta file $file_fasta does not exist!" if not -f $file_fasta;
	confess "ERROR! tbl file $file_tbl does not exist!" if not -f $file_tbl;
	$flag_dihe = 0 if not defined $flag_dihe;
	my $seq = seq_fasta($file_fasta);
	my $L = length $seq;
	my $cov = $seq;
	$cov =~ s/[A-Z]/-/g;
	my %noe = tbl2rows_hash($file_tbl);
	foreach (keys %noe){
		#assign (resid 123 and name CA) (resid  58 and name CA) 3.60 0.10 3.40
		my @C = split /\s+/, $_;
		my $r1 = $C[2]; my $r2 = $C[7];
		# the case when "ca or cb" restraints are provided
		#assign ((resid 123 and name ca) or (resid 123 and name cb)) ((resid 58 and name ca) or (resid 58 and name cb)) 3.60 0.10 3.40
		$r2 = $C[13] if $C[6] eq "or";
		$r2 = $C[17] if $flag_dihe;
		my $c1 = substr $cov, ($r1 - 1), 1;
		my $c2 = substr $cov, ($r2 - 1), 1;
		if ($c1 eq "-" ){
			$c1 = 1;
		}
		elsif ($c1 eq "*" ){
			$c1 = "*";
		}
		else{
			$c1++;
			$c1 = "*" if ($c1 > 9);
		}
		if ($c2 eq "-" ){
			$c2 = 1;
		}
		elsif ($c2 eq "*" ){
			$c2 = "*";
		}
		else{
			$c2++;
			$c2 = "*" if ($c2 > 9);
		}
		substr $cov, ($r1 - 1), 1, $c1;
		substr $cov, ($r2 - 1), 1, $c2;
	}
	my $cov2 = $cov;
	$cov2 =~ s/-//g;
	return sprintf "$cov [%12s : %3s restraints touching %s residues]", basename($file_tbl), (scalar keys %noe), length($cov2);
}

sub count_satisfied_tbl_rows{
	my $file_pdb = shift;
	my $file_tbl = shift;
	my $type    = shift;
	confess "ERROR! file pdb $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! file tbl $file_tbl does not exist!" if not -f $file_tbl;
	confess "ERROR! Invalid type!" if not ($type eq "noe" or $type eq "dihedral");
	my $count = 0;
	my $total = 0;
	my %log_rows = ();
	my %tbl_hash = ssnoe_tbl_min_pdb_dist($file_tbl, $file_pdb);
	foreach (sort keys %tbl_hash){
		my $viol_flag = 1;
		my $distance = $tbl_hash{$_}{"distance"};
		my @D = split /\s+/, $distance;
		my $pdb_distance = $tbl_hash{$_}{"pdb_distance"};
		my $deviation = $pdb_distance - ( $D[0] + $D[2] );
		if ($pdb_distance < ( $D[0] + $D[2] + $DISTRELAX) ){
			$count ++;
			$viol_flag = 0;
			$deviation = 0.0;
		}
		if ($pdb_distance < ( $D[0] - $D[1] - $DISTRELAX) ){
			$count--;
			$viol_flag = 1;
			$deviation = -($D[0] - $D[1] - $pdb_distance);
		}
		
		$log_rows{(sprintf "%3s\t%.2f\t%.2f # $_", $viol_flag, $deviation, $pdb_distance)} = $viol_flag;
		$total++;
	}
	my $file_viol = "".basename($file_tbl,".tbl")."_violation.txt";
	print2file($file_viol, "#NOE violation check; $file_pdb against $file_tbl");
	print2file($file_viol, "#violation-flag, deviation, actual-measurement, Input-NOE-restraint");
	foreach (sort {$log_rows{$b} <=> $log_rows{$a}} keys %log_rows){
		print2file($file_viol, $_);
	}
	return $count."/".$total;
}

sub ssnoe_tbl_min_pdb_dist{
	my $file_tbl = shift;
	my $file_pdb = shift;
	confess ":(" if not -f $file_tbl;
	confess ":(" if not -f $file_pdb;
	my %noe_hash = ();
	open NOETBL, $file_tbl or confess $!;
	while (<NOETBL>){
		chomp $_;
		$_ =~ s/^\s+//;
		$_ =~ s/\)/ /g;
		$_ =~ s/\(/ /g;
		confess "$_" if $_ !~ m/^assign/;
		my @C = split /\s+/, $_;
		if ($C[6] eq "or" and $C[17] eq "or"){
			$noe_hash{$_}{"left"}     = $C[2]." ".$C[5]." ".$C[8]." ".$C[11];
			$noe_hash{$_}{"right"}    = $C[13]." ".$C[16]." ".$C[19]." ".$C[22];
			$noe_hash{$_}{"distance"} = $C[23]." ".$C[24]." ".$C[25];
		}
		elsif ($C[6] eq "or" and $C[17] ne "or"){
			$noe_hash{$_}{"left"}     = $C[2]." ".$C[5]." ".$C[8]." ".$C[11];
			$noe_hash{$_}{"right"}    = $C[13]." ".$C[16];
			$noe_hash{$_}{"distance"} = $C[17]." ".$C[18]." ".$C[19];
		}
		elsif ($C[6] ne "or" and $C[11] eq "or"){
			$noe_hash{$_}{"left"}     = $C[2]." ".$C[5];
			$noe_hash{$_}{"right"}    = $C[7]." ".$C[10]." ".$C[13]." ".$C[16];
			$noe_hash{$_}{"distance"} = $C[17]." ".$C[18]." ".$C[19];
		}
		else{
			$noe_hash{$_}{"left"}     = $C[2]." ".$C[5];
			$noe_hash{$_}{"right"}    = $C[7]." ".$C[10];
			$noe_hash{$_}{"distance"} = $C[11]." ".$C[12]." ".$C[13];
		}
	}
	close NOETBL;
	confess "$file_tbl seems empty!" if not %noe_hash;
	my %xyzPDB = xyz_pdb($file_pdb, "all");
	foreach (sort keys %noe_hash){
		my $left = $noe_hash{$_}{"left"};
		my $right = $noe_hash{$_}{"right"};
		my $distance = $noe_hash{$_}{"distance"};
		my @L = split /\s+/, $left;
		my @R = split /\s+/, $right;
		my @D = split /\s+/, $distance;
		my %left_list = ();
		my %right_list = ();
		for(my $i = 0; $i <= $#L; $i = $i+2){
			$left_list{$L[$i]." ".$L[$i+1]} = 1;
		}
		for(my $i = 0; $i <= $#R; $i = $i+2){
			$right_list{$R[$i]." ".$R[$i+1]} = 1;
		}
		my $distance_pdb = 1000.0;
		foreach my $le (keys %left_list){
			foreach my $ri (keys %right_list){
				my @L = split /\s+/, $le;
				my @R = split /\s+/, $ri;
				confess "$file_pdb does not have ".$L[0]." ".uc($L[1])."\n" if not defined $xyzPDB{$L[0]." ".uc($L[1])};
				confess "$file_pdb does not have ".$R[0]." ".uc($R[1])."\n" if not defined $xyzPDB{$R[0]." ".uc($R[1])};
				my $d = calc_dist($xyzPDB{$L[0]." ".uc($L[1])}, $xyzPDB{$R[0]." ".uc($R[1])});
				$distance_pdb = $d if $d < $distance_pdb;
			}
		}
		$noe_hash{$_}{"pdb_distance"} = $distance_pdb;
	}
	return %noe_hash;
}

sub noe_tbl_violation_coverage{
	my $file_pdb = shift;
	my $file_tbl = shift;
	confess "ERROR! file pdb $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! file tbl $file_tbl does not exist!" if not -f $file_tbl;
	my $cov = seq_chain($file_pdb);
	$cov =~ s/[A-Z]/-/g;
	my %tbl_hash = ssnoe_tbl_min_pdb_dist($file_tbl, $file_pdb);
	my %xyz = xyz_pdb($file_pdb, "all");
	foreach (sort keys %tbl_hash){
		my $left = $tbl_hash{$_}{"left"};
		my $right = $tbl_hash{$_}{"right"};
		my $distance = $tbl_hash{$_}{"distance"};
		my @L = split /\s+/, $left;
		my @R = split /\s+/, $right;
		my @D = split /\s+/, $distance;
		my $pdb_distance = $tbl_hash{$_}{"pdb_distance"};
		substr $cov, $L[0] - 1, 1, "x" if $pdb_distance > $D[0] + $D[2] + 0.2;
		substr $cov, $R[0] - 1, 1, "x" if $pdb_distance > $D[0] + $D[2] + 0.2;
		substr $cov, $L[0] - 1, 1, "x" if $pdb_distance < $D[0] - $D[1] - 0.2;
		substr $cov, $R[0] - 1, 1, "x" if $pdb_distance < $D[0] - $D[1] - 0.2;
	}
	return $cov;
}

sub sum_noe_dev{
	my $file_pdb = shift;
	my $file_tbl = shift;
	confess "ERROR! file pdb $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! file tbl $file_tbl does not exist!" if not -f $file_tbl;
	my $sum_dev = 0.0;
	my %tbl_hash = ssnoe_tbl_min_pdb_dist($file_tbl, $file_pdb);
	foreach (sort keys %tbl_hash){
		my $viol_flag = 1;
		my @D = split /\s+/, $tbl_hash{$_}{"distance"};
		my $pdb_distance = $tbl_hash{$_}{"pdb_distance"};
		if ($pdb_distance > ( $D[0] + $D[2] + 0.2) ){
			$sum_dev += $pdb_distance - ( $D[0] + $D[2] );
		}
		if ($pdb_distance < ( $D[0] - $D[1] - 0.2) ){
			$sum_dev += ( $D[0] - $D[1] ) - $pdb_distance;
		}
	}
	return sprintf "%.2f", $sum_dev;
}

sub get_cns_energy{
	my $cns_pdb = shift;
	my $energy_term = shift; # "overall", "vdw", "bon", "noe"
	confess "ERROR! file $cns_pdb does not exist!" if not -f $cns_pdb;
	confess "ERROR! energy term must be overall vdw bon or noe" if not ($energy_term eq "overall" or $energy_term eq "bon" or $energy_term eq "vdw" or $energy_term eq "noe");
	my $value = "X";
	open CHAIN, $cns_pdb or confess $!;
	while(<CHAIN>){
		chomp $_;
		next if $_ !~ m/^REMARK\ $energy_term/;
		$_ =~ s/\s+//g;
		my @C = split /=/, $_;
		$value = $C[1];
	}
	close CHAIN;
	return int($value);
}

sub load_pdb{
	my $dir_chains = shift;
	confess ":( directory $dir_chains does not exist!" if not -d $dir_chains;
	my @pdb_list = <$dir_chains/*.pdb>;
	if(not (@pdb_list)){
		@pdb_list = <$dir_chains/*.ent>;
	}
	confess "ERROR! Directory $dir_chains has no pdb files!\n" unless(@pdb_list);
	return @pdb_list;
}

sub pdb2rnum_rname{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my %rnum_rname = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$rnum_rname{parse_pdb_row($_,"rnum")} = parse_pdb_row($_,"rname");
	}
	close CHAIN;
	confess ":(" if not scalar keys %rnum_rname;
	return %rnum_rname;
}

sub xyz_pdb{
	my $chain = shift;
	my $atom_selection = shift; # ca or cb or all
	$atom_selection = uc($atom_selection);
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	confess "ERROR! Selection must be ca or cb or all" if not ($atom_selection eq "CA" or $atom_selection eq "ALL" or $atom_selection eq "CB");
	my %xyz_pdb = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$xyz_pdb{"".parse_pdb_row($_,"rnum")." ".parse_pdb_row($_,"aname")} = "".parse_pdb_row($_,"x")." ".parse_pdb_row($_,"y")." ".parse_pdb_row($_,"z");
	}
	close CHAIN;
	confess "ERROR!: xyz_pdb is empty\n" if (not scalar keys %xyz_pdb);
	return %xyz_pdb if $atom_selection eq "ALL";

	my %rnum_rname = pdb2rnum_rname($chain);
	my %selected_xyz = ();
	foreach (sort keys %xyz_pdb){
		my @C = split /\s+/, $_;
		my $this_atom = $atom_selection;
		$this_atom = "CA" if ($atom_selection eq "CB" and $rnum_rname{$C[0]} eq "GLY");
		next if $C[1] ne $this_atom;
		$selected_xyz{$C[0]} = $xyz_pdb{$_};
	}
	confess ":(" if not scalar keys %selected_xyz;
	return %selected_xyz;
}

sub parse_pdb_row{
	my $row = shift;
	my $param = shift;
	my $result;
	$result = substr($row,6,5) if ($param eq "anum");
	$result = substr($row,12,4) if ($param eq "aname");
	$result = substr($row,16,1) if ($param eq "altloc");
	$result = substr($row,17,3) if ($param eq "rname");
	$result = substr($row,22,5) if ($param eq "rnum");
	$result = substr($row,26,1) if ($param eq "insertion");
	$result = substr($row,21,1) if ($param eq "chain");
	$result = substr($row,30,8) if ($param eq "x");
	$result = substr($row,38,8) if ($param eq "y");
	$result = substr($row,46,8) if ($param eq "z");
	confess "Invalid row[$row] or parameter[$param]" if (not defined $result);
	$result =~ s/\s+//g;
	return $result;
}

sub clash_count{
	my $file_pdb = shift;
	my $threshold = shift;
	confess "ERROR! pdb file $file_pdb not defined!" if not -f $file_pdb;
	confess "ERROR! threshold not defined!" if not defined $threshold;
	my $count = 0;
	my %ca_xyz = xyz_pdb($file_pdb, "ca");
	foreach my $r1 (sort keys %ca_xyz){
		my @R1 = split /\s+/, $ca_xyz{$r1};
		my $x1 = $R1[0]; my $y1 = $R1[1]; my $z1 = $R1[2];
		foreach my $r2 (sort keys %ca_xyz){
			next if $r1 >= $r2;
			my @R2 = split /\s+/, $ca_xyz{$r2};
			my $x2 = $R2[0]; my $y2 = $R2[1]; my $z2 = $R2[2];
			my $d = sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
			if ($d <= $threshold){
				$count++;
			}
		}
	}
	return $count;
}

sub calc_dist{
	my $x1y1z1 = shift;
	my $x2y2z2 = shift;
	confess "ERROR! x1y1z1 not defined!" if !$x1y1z1;
	confess "ERROR! x2y2z2 not defined!" if !$x2y2z2;
	my @row1 = split(/\s+/, $x1y1z1);
	my $x1 = $row1[0]; my $y1 = $row1[1]; my $z1 = $row1[2];
	confess "ERROR! One of the xyz in $x1y1z1 is not a number!" if (not looks_like_number($x1) or not looks_like_number($y1) or not looks_like_number($z1));
	my @row2 = split(/\s+/, $x2y2z2);
	my $x2 = $row2[0]; my $y2 = $row2[1]; my $z2 = $row2[2];
	confess "ERROR! One of the xyz in $x2y2z2 is not a number!" if (not looks_like_number($x2) or not looks_like_number($y2) or not looks_like_number($z2));
	my $d = sprintf "%.3f", sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
	return $d;
}

sub seq_fasta{
	my $file_fasta = shift;
	confess "ERROR! Fasta file $file_fasta does not exist!" if not -f $file_fasta;
	my $seq = "";
	open FASTA, $file_fasta or confess $!;
	while (<FASTA>){
		next if (substr($_,0,1) eq ">"); 
		chomp $_;
		$_ =~ tr/\r//d; # chomp does not remove \r
		$seq .= $_;
	}
	close FASTA;
	return $seq;
}

sub write_cns_seq{
	my $file_fasta = shift;
	my $file_cns_seq = shift;
	confess "ERROR! Fasta file $file_fasta does not exist!" if not -f $file_fasta;
	confess "ERROR! Output not defined!" if not defined $file_cns_seq;
	my @seq = split //, seq_fasta($file_fasta);
	my $three_letter_seq;
	foreach (@seq) {
		$three_letter_seq .= $AA1TO3{$_}." ";
	}
	system_cmd("rm -f $file_cns_seq");
	while($three_letter_seq){
		if(length ($three_letter_seq) <= 64 ){
			print2file($file_cns_seq, $three_letter_seq);
			$three_letter_seq = ();
		}
		else{
			print2file($file_cns_seq, substr($three_letter_seq, 0, 64));
			$three_letter_seq = substr($three_letter_seq, 64);
		}
	}
}

sub assess_dgsa{
	my @seq_files = <./*.fasta>;
	confess ":( fasta not found!" if not defined $file_seq;
	my $ID = basename($file_seq, ".fasta");
	my @pdbList = load_pdb($dir_out);
	# also verify if CNS actually accepted all the distance restraints provided
	my $search_string = "N1";
	if(not -f $dg_sa_log){
		# PROBABLY THE SERVER DOES NOT HAVE SPACE!
		system_cmd("touch dgsa.log.noverification");
	}
	else{
		my $result = `grep NOEPRI $dg_sa_log | grep $search_string | head -n 1`;
		my @C = split /\s+/, $result;
		my $count = $C[($#C)-1];
		if ($count != count_lines("contact.tbl")){
			system_cmd("touch assess.failed");
			confess ":( CNS did not accept all restraints of contact.tlb! Something wrong somewhere! Only $count accepted";
		}
	}

	# remove "trial" structure of corresponding "accepted" structure because they are same
	for(my $i = 0; $i <= 1000; $i++){
		next if not -f "${ID}a_$i.pdb";
		print "\ndeleting ${ID}_$i.pdb because ${ID}a_$i.pdb exists!";
		system_cmd("rm ./${ID}_$i.pdb");
	}
	my %energyNoe = ();
	foreach my $pdb (@pdbList) {
		next if $pdb =~ m/sub_embed/;
		next if $pdb =~ m/extended/;
		$energyNoe{$pdb} = get_cns_energy($pdb, "noe");
	}
	my $bestPdb = (sort {$energyNoe{$a} <=> $energyNoe{$b}} keys %energyNoe)[0];
	print "\n";
	print "NOE_SATISFIED(±${DISTRELAX}A)  SUM_OF_DEVIATIONS>= 0.2  PDB\n";
	foreach my $pdb (sort {$energyNoe{$b} <=> $energyNoe{$a}} keys %energyNoe){
		my ($n1, $s1);
		$n1 = count_satisfied_tbl_rows($pdb, "contact.tbl", "noe");
		$s1 = sum_noe_dev($pdb, "contact.tbl");
		printf "%-9s             %-9s                %-25s\n", $n1, $s1, basename($pdb, ".pdb"); 
	}
	print "\n";
	print "removing non-CA ATOM rows and backing up REMARK rows..\n";
	foreach my $pdb (sort {$energyNoe{$b} <=> $energyNoe{$a}} keys %energyNoe){
		system_cmd("rm -f ca_filtered.pdb");
		filter_nonCA($pdb, "ca_filtered.pdb", $model_info_log);
		reindex_chain("ca_filtered.pdb", 1, $pdb); # to reorder the atom numbers as well
		system_cmd("rm -f ca_filtered.pdb");
		system_cmd("sed -i \"s/END//g\" $pdb");
		add_connect_rows($pdb);
	}
	print "\n";
	my $i = 1;
	foreach( sort {$energyNoe{$a} <=> $energyNoe{$b}} keys %energyNoe){
		print "model$i.pdb <= $_\n";
		system_cmd("mv $_ ${ID}_model$i.pdb");
		$i++;
		last if $i > 5;
	}
}

sub reindex_chain{
	my $file_pdb = shift;
	my $index = shift;
	my $out_pdb = shift;
	confess "ERROR! file $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! index $index is invalied!" if not defined $index;
	my $res_counter = $index - 1;
	my $atom_counter = 0;
	my $prev_res_num = "XX";
	open OUTPDB, ">$out_pdb" or confess $!;
	open PDBFILE, $file_pdb or confess $!;
	foreach (<PDBFILE>) {
		last if $_ =~ /^TER/;
		last if $_ =~ /^ENDMDL/;
		next if $_ !~ m/^ATOM/;
		next if not ((parse_pdb_row($_,"altloc") eq "") or (parse_pdb_row($_,"altloc") eq "A"));
		next if not defined $AA3TO1{parse_pdb_row($_,"rname")};
		my $this_rnum = parse_pdb_row($_,"rnum");
		if ($prev_res_num ne $this_rnum) {
			$prev_res_num = $this_rnum;
			$res_counter++;
		}
		$atom_counter++;
		my $rnum_string = sprintf("%4s", $res_counter);
		my $anum_string = sprintf("%5s", $atom_counter);
		my $row = substr($_,0,6).$anum_string.substr($_,11,5)." ".substr($_,17,3)." "." ".$rnum_string." ".substr($_,27);
		print OUTPDB $row;
	}
	close PDBFILE;
	print OUTPDB "END\n";
	close OUTPDB;
}

sub filter_nonCA{
	my $in_pdb = shift;
	my $out_pdb = shift;
	my $log_file = shift;
	confess "ERROR! file does not exist!" if not -f $in_pdb;
	print2line($log_file, $in_pdb);
	open FILEPDB, $in_pdb or confess "ERROR! Could not open $in_pdb! $!";
	while (<FILEPDB>){
		my $r = $_;
		print2line($log_file, $r) if ($r =~ /^REMARK/);
		next if $r !~ /^ATOM/;
		next if ($r =~ /^ATOM/ && $r !~ /CA/);
		print2line($out_pdb, $r);
	}
	close FILEPDB;
	print2line($log_file, "");
}

sub write_cns_dgsa_file{
	print2file("dgsa.inp", "{+ file: dgsa.inp +}");
	print2file("dgsa.inp", "{+ directory: nmr_calc +}");
	print2file("dgsa.inp", "{+ description: distance geometry, full or substructure, with ");
	print2file("dgsa.inp", "                simulated annealing regularization starting from ");
	print2file("dgsa.inp", "                extended strand or pre-folded structures. +}");
	print2file("dgsa.inp", "{+ authors: Gregory Warren, Michael Nilges, John Kuszewski, ");
	print2file("dgsa.inp", "	    Marius Clore and Axel Brunger +}");
	print2file("dgsa.inp", "{+ copyright: Yale University +}");
	print2file("dgsa.inp", "{+ reference: Clore GM, Gronenborn AM, Tjandra N, Direct structure refinement ");
	print2file("dgsa.inp", "              against residual dipolar couplings in the presence of rhombicity");
	print2file("dgsa.inp", "              of unknown magnitude., J. Magn. Reson., 131, In press, (1998) +}");
	print2file("dgsa.inp", "{+ reference: Clore GM, Gronenborn AM, Bax A, A robust method for determining ");
	print2file("dgsa.inp", "              the magnitude of the fully asymmetric alignment tensor of");
	print2file("dgsa.inp", "              oriented macromolecules in the absence of structural");
	print2file("dgsa.inp", "              information., J. Magn. Reson., In press (1998) +}");
	print2file("dgsa.inp", "{+ reference: Garrett DS, Kuszewski J, Hancock TJ, Lodi PJ, Vuister GW,");
	print2file("dgsa.inp", "              Gronenborn AM, Clore GM, The impact of direct refinement against ");
	print2file("dgsa.inp", "              three-bond HN-C alpha H coupling constants on protein structure");
	print2file("dgsa.inp", "              determination by NMR., J. Magn. Reson. Ser. B, 104(1), ");
	print2file("dgsa.inp", "              99-103, (1994) May +}");
	print2file("dgsa.inp", "{+ reference: Kuszewski J, Nilges M, Brunger AT,   Sampling and efficiency ");
	print2file("dgsa.inp", "              of metric matrix distance geometry:  A novel partial metrization ");
	print2file("dgsa.inp", "              algorithm.  J. Biomol. NMR 2, 33-56, (1992). +} ");
	print2file("dgsa.inp", "{+ reference: Kuszewski J, Qin J, Gronenborn AM, Clore GM, The impact of direct");
	print2file("dgsa.inp", "              refinement against 13C alpha and 13C beta chemical shifts on ");
	print2file("dgsa.inp", "              protein structure determination by NMR., J. Magn. Reson. Ser. B,");
	print2file("dgsa.inp", "              106(1), 92-6, (1995) Jan +}");
	print2file("dgsa.inp", "{+ reference: Kuszewski J, Gronenborn AM, Clore GM, The impact of direct");
	print2file("dgsa.inp", "              refinement against proton chemical shifts on protein structure ");
	print2file("dgsa.inp", "              determination by NMR., J. Magn. Reson. Ser. B, 107(3), 293-7, ");
	print2file("dgsa.inp", "              (1995) Jun +}");
	print2file("dgsa.inp", "{+ reference: Kuszewski J, Gronenborn AM, Clore GM, A potential involving ");
	print2file("dgsa.inp", "              multiple proton chemical-shift restraints for ");
	print2file("dgsa.inp", "              nonstereospecifically assigned methyl and methylene protons.");
	print2file("dgsa.inp", "              J. Magn. Reson. Ser. B, 112(1), 79-81, (1996) Jul. +}");
	print2file("dgsa.inp", "{+ reference: Nilges M, Clore GM, Gronenborn AM, Determination of ");
	print2file("dgsa.inp", "              three-dimensional structures of proteins from interproton ");
	print2file("dgsa.inp", "              distance data by hybrid distance geometry-dynamical simulated ");
	print2file("dgsa.inp", "              annealing calculations. FEBS Lett. 229, 317-324 (1988). +}");
	print2file("dgsa.inp", "{+ reference: Nilges M, Clore GM, Gronenborn AM,  Determination of ");
	print2file("dgsa.inp", "              three-dimensional structures of proteins from interproton ");
	print2file("dgsa.inp", "              distance data by dynamical simulated annealing from a random ");
	print2file("dgsa.inp", "              array of atoms. FEBS LEtt. 239, 129-136 (1988). +}");
	print2file("dgsa.inp", "{+ reference: Nilges M, Kuszewski J, Brunger AT, In: Computational Aspects ");
	print2file("dgsa.inp", "              of the Study of Biological Macromolecules by NMR. ");
	print2file("dgsa.inp", "              (J.C. Hoch, ed.),  New York: Plenum Press, (1991). +}");
	print2file("dgsa.inp", "{+ reference: Tjandra N, Garrett DS, Gronenborn AM, Bax A, Clore GM, Defining");
	print2file("dgsa.inp", "              long range order in NMR structure determination from the ");
	print2file("dgsa.inp", "              dependence of heteronuclear relaxation times on rotational ");
	print2file("dgsa.inp", "              diffusion anisotropy. Nature Struct. Biol., 4(6), 443-9,");
	print2file("dgsa.inp", "              (1997) June +}");
	print2file("dgsa.inp", "{+ reference: Tjandra N, Omichinski JG, Gronenborn AM, Clore GM, Bax A, Use of");
	print2file("dgsa.inp", "              dipolar 1H-15N and 1H-13C couplings in the structure");
	print2file("dgsa.inp", "              determination of magnetically oriented macromolecules in");
	print2file("dgsa.inp", "              solution. Nature Struct. Biol., 4(9), 732-8, (1997) Sept +} ");
	print2file("dgsa.inp", "              ");
	print2file("dgsa.inp", "{- Guidelines for using this file:");
	print2file("dgsa.inp", "   - all strings must be quoted by double-quotes");
	print2file("dgsa.inp", "   - logical variables (true/false) are not quoted");
	print2file("dgsa.inp", "   - do not remove any evaluate statements from the file -}");
	print2file("dgsa.inp", "{- begin block parameter definition -} define(");
	print2file("dgsa.inp", "{======================= molecular structure =========================}");
	print2file("dgsa.inp", "{* parameter file(s) *}");
	# CNS_TOPPAR:protein-allhdg5-4.param
	print2file("dgsa.inp", "{===>} par.1=\"CNS_TOPPAR:protein.param\";");
	print2file("dgsa.inp", "{===>} par.2=\"\";");
	print2file("dgsa.inp", "{===>} par.3=\"\";");
	print2file("dgsa.inp", "{===>} par.4=\"\";");
	print2file("dgsa.inp", "{===>} par.5=\"\";");
	print2file("dgsa.inp", "{* structure file(s) *}");
	print2file("dgsa.inp", "{===>} struct.1=\"extended.mtf\";");
	print2file("dgsa.inp", "{===>} struct.2=\"\";");
	print2file("dgsa.inp", "{===>} struct.3=\"\";");
	print2file("dgsa.inp", "{===>} struct.4=\"\";");
	print2file("dgsa.inp", "{===>} struct.5=\"\";");
	print2file("dgsa.inp", "{* input coordinate file(s) *}");
	print2file("dgsa.inp", "{===>} pdb.in.file.1=\"extended.pdb\";");
	print2file("dgsa.inp", "{===>} pdb.in.file.2=\"\";");
	print2file("dgsa.inp", "{===>} pdb.in.file.3=\"\";");
	print2file("dgsa.inp", "{========================== atom selection ===========================}");
	print2file("dgsa.inp", "{* input \"backbone\" selection criteria for average structure generation *}");
	print2file("dgsa.inp", "{* for protein      (name n or name ca or name c)");
	print2file("dgsa.inp", "   for nucleic acid (name O5' or name C5' or name C4' or name C3' ");
	print2file("dgsa.inp", "                     or name O3' or name P) *}");
	print2file("dgsa.inp", "{===>} pdb.atom.select=(name n or name ca or name c);");
	print2file("dgsa.inp", "{======================= refinement parameters ========================}");
	print2file("dgsa.inp", "{* distance geometry *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} flg.dg.flag=true;");
	print2file("dgsa.inp", "{* distance geometry/simualted annealing regularization (DGSA) *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} flg.dgsa.flag=true;");
	print2file("dgsa.inp", "{* if only regularizing coordinate files (no DG) then ");
	print2file("dgsa.inp", "   enter the number of coordinate files to be regularized (DGSA) *}");
	print2file("dgsa.inp", "{===>} pdb.dg.count=$model_count;");
	print2file("dgsa.inp", "{* seed for random number generator *}");
	print2file("dgsa.inp", "{* change to get different initial velocities *}");
	print2file("dgsa.inp", "{===>} md.seed=82364;");
	print2file("dgsa.inp", "{* select whether the number of structures will be either trial or 	");
	print2file("dgsa.inp", "   accepted structures and whether to print only the trial, accepted, 	");
	print2file("dgsa.inp", "   both sets of structures. The printing format is as follows:");
	print2file("dgsa.inp", "   trial = pdb.out.name + _#.pdb , accepted = pdb.out.name + a_#.pdb *} ");
	print2file("dgsa.inp", "{* are the number of structures to be trials or accepted? *}");
	print2file("dgsa.inp", "{+ choice: \"trial\" \"accept\" +}");
	print2file("dgsa.inp", "{===>} flg.trial.struc=\"$mode\";");
	print2file("dgsa.inp", "{* number of trial or accepted structures *}");
	print2file("dgsa.inp", "{===>} pdb.end.count=$model_count;");
	print2file("dgsa.inp", "{* print accepted structures *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} flg.print.accept=true;");
	print2file("dgsa.inp", "{* print trial structures *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} flg.print.trial=true;");
	print2file("dgsa.inp", "{* calculate an average structure for either the trial or 	");
	print2file("dgsa.inp", "   accepted structure.  If calculate accepted average is false then ");
	print2file("dgsa.inp", "   an average for the trial structures will be calculated. *}");
	print2file("dgsa.inp", "{* calculate an average structure? *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} flg.calc.ave.struct=false;");
	print2file("dgsa.inp", "{* calculate an average structure for the accepted structures? *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} flg.calc.ave.accpt=false;");
	print2file("dgsa.inp", "{* minimize average coordinates? *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} flg.min.ave.coor=false;");
	print2file("dgsa.inp", "{============ parameters for the distance geometry stage ==============}");
	print2file("dgsa.inp", "{* shortest path algorithm *}");
	print2file("dgsa.inp", "{+ choice: \"auto\" \"full\" \"sparse\" +}");
	print2file("dgsa.inp", "{===>} md.dg.algo=\"auto\";");
	print2file("dgsa.inp", "{* distance geometry on substructure or complete structure? *}");
	print2file("dgsa.inp", "{* proteins: \"sub\" or \"complete\"; dna/rna: \"complete\" *}");
	print2file("dgsa.inp", "{+ choice: \"sub\" \"complete\" +}");
	print2file("dgsa.inp", "{===>} md.dg.type=\"sub\";");
	print2file("dgsa.inp", "{* input atom selection for substructure  *}");
	# Atom selection scheme for distance geometry
	# (1) as-is / ca + ha + n + hn + c + cb + cg (2) as-is + o (3) as-is + o + h 
	# (2) as-is + o (3) as-is + o + h (4) c + cα + n + o (backbone atoms)
	# (5) bkbone + cβ (6) bkbone + cβ + h (7)  bkbone + cβ + h + cg                              
	# Notes: 
	# (a) CNS recommended is 6, but with option 6, error "%MATROT error encountered: Error in internal consistency check" was observed for some pdbs
	# (b) with option 7, chirality issues were faced on target 1QJP
	# (c) with option 2 (best), maximum reconstruction accuracy and sheet_detection results for 150 fragfold proteins.
	if ($atomselect == 1){
		# as-is
		print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name ha or name n or name hn");
		print2file("dgsa.inp", "		        or name c or name cb* or name cg*);");
	}
	if ($atomselect == 2){
		# add oxygen
		print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name ha or name n or name hn");
		print2file("dgsa.inp", "		        or name c or name cb* or name cg* or name o);");
	}
	if ($atomselect == 3){
		# add oxygen and hydrogen
		print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name ha or name n or name hn or name h");
		print2file("dgsa.inp", "		        or name c or name cb* or name cg* or name o);");
	}
	if ($atomselect == 4){
		# backbone atoms only
		print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name c or name n or name o);");
	}
	if ($atomselect == 5){
		# backbone atoms with cb
		print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name c or name n or name o");
		print2file("dgsa.inp", "		        or name cb);");
	}
	if ($atomselect == 6){
		# backbone atoms with cb and hydrogen
		# atom selection according to instructions in the NIH-XPLORE manual
		print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name c or name n or name o");
		print2file("dgsa.inp", "		        or name cb or name h);");
	}
	if ($atomselect == 7){
		# backbone atoms with cb, hydrogen and cg
		print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name c or name n or name o");
		print2file("dgsa.inp", "		        or name cb or name h or name cg*);");
	}
	print2file("dgsa.inp", "{* when using \"complete\" input the rigid group atom selection *}");
	print2file("dgsa.inp", "{===>} md.dg.group.slct=(known);");
	print2file("dgsa.inp", "{* group interatomic error value in angstroms *}");
	print2file("dgsa.inp", "{===>} md.dg.group.err=0.5;");
	print2file("dgsa.inp", "{* use metrization for complete distance geometry? *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} md.dg.metr.flag=false;");
	print2file("dgsa.inp", "{* ordered or random metrization *}");
	print2file("dgsa.inp", "{+ choice: \"ordered\" \"random\" +}");
	print2file("dgsa.inp", "{===>} md.dg.ord=\"random\";");
	print2file("dgsa.inp", "{* input metrization atom selection *}");
	print2file("dgsa.inp", "{===>} md.dg.metr.atom=(all);");
	print2file("dgsa.inp", "{* input number of atoms from selection used during retightening *}");
	print2file("dgsa.inp", "{===>} md.dg.metr.num=4;");
	print2file("dgsa.inp", "{* reference for building the reference data base *}");
	print2file("dgsa.inp", "{+ choice: \"parameter\" \"coordinate\" +}");
	print2file("dgsa.inp", "{===>} md.dg.ref=\"parameter\";");
	print2file("dgsa.inp", "{* scale factor for distance geometry restraint term *}");
	print2file("dgsa.inp", "{===>} md.dg.scale=100.;");
	print2file("dgsa.inp", "{* exponent for distance geometry restraint term *}");
	print2file("dgsa.inp", "{===>} md.dg.expo=2;");
	print2file("dgsa.inp", "{* bond length (in angstroms) error value *}");
	print2file("dgsa.inp", "{===>} md.dg.bacc=0.01;");
	print2file("dgsa.inp", "{* angle (in degrees) error value *}");
	print2file("dgsa.inp", "{===>} md.dg.tacc=2.;");
	print2file("dgsa.inp", "{* improper (in degrees) error value *}");
	print2file("dgsa.inp", "{===>} md.dg.iacc=2.;");
	print2file("dgsa.inp", "{* dihedral (in degrees) error value *}");
	print2file("dgsa.inp", "{===>} md.dg.pacc=2.;");
	print2file("dgsa.inp", "{* number of steps of minimization *}");
	print2file("dgsa.inp", "{===>} md.dg.step=200;");
	print2file("dgsa.inp", "{=== parameters for the distance geometry/simulated annealing stage ===}");
	print2file("dgsa.inp", "{* starting temperature *}");
	print2file("dgsa.inp", "{===>} md.hot.temp=2000;");
	print2file("dgsa.inp", "{* number of steps for high temperature dyanmics *}");
	print2file("dgsa.inp", "{===>} md.hot.step=1000;");
	print2file("dgsa.inp", "{* number of steps for slow-cool annealing *}");
	print2file("dgsa.inp", "{===>} md.cool.step=1000;");
	print2file("dgsa.inp", "{* hot molecular dynamics timestep *}");
	print2file("dgsa.inp", "{===>} md.hot.ss=0.003;");
	print2file("dgsa.inp", "{* slow-cool molecular dynamics timestep *}");
	print2file("dgsa.inp", "{===>} md.cool.ss=0.005;");
	print2file("dgsa.inp", "{* initial scale factor for van der Waals (repel) energy term *}");
	print2file("dgsa.inp", "{===>} md.cool.vdw.init=0.003;");
	print2file("dgsa.inp", "{* final scale factor for van der Waals (repel) energy term *}");
	print2file("dgsa.inp", "{===>} md.cool.vdw.finl=4.0;");
	print2file("dgsa.inp", "{* initial van der Waals repel radius *}");
	print2file("dgsa.inp", "{===>} md.cool.init.rad=$rep1;"); # 0.9 originally
	print2file("dgsa.inp", "{* final van der Waals repel radius *}");
	print2file("dgsa.inp", "{===>} md.cool.fina.rad=$rep2;"); # 0.8 originally
	print2file("dgsa.inp", "{* scale factor for NOE energy term *}");
	print2file("dgsa.inp", "{===>} md.cool.noe=$con_wt;"); # 50 originally
	print2file("dgsa.inp", "{* high temperature scale factor for dihedral angle energy term *}");
	print2file("dgsa.inp", "{===>} md.hot.cdih=5;");
	print2file("dgsa.inp", "{* slow-cooling scale factor for dihedral angle energy term *}");
	print2file("dgsa.inp", "{===>} md.cool.cdih=200;"); # 200 originally
	print2file("dgsa.inp", "{* slow-cool annealing temperature step *}");
	print2file("dgsa.inp", "{===>} md.cool.tmpstp=25.;");
	print2file("dgsa.inp", "{=============== parameters for final minimization stage ==============}");
	print2file("dgsa.inp", "{* scale factor for NOE energy term *}");
	print2file("dgsa.inp", "{===>} md.pow.noe=$con_wt;"); # 50 originally
	print2file("dgsa.inp", "{* scale factor for dihedral angle energy term *}");
	print2file("dgsa.inp", "{===>} md.pow.cdih=400;"); # 400 originally
	print2file("dgsa.inp", "{* number of minimization steps *}");
	print2file("dgsa.inp", "{===>} md.pow.step=$mini;"); # 200 originally
	print2file("dgsa.inp", "{* number of cycles of minimization *}");
	print2file("dgsa.inp", "{===>} md.pow.cycl=10;");
	print2file("dgsa.inp", "      ");
	print2file("dgsa.inp", "{============================= noe data ===============================}");
	print2file("dgsa.inp", "{- Important - if you do not have a particular data set then");
	print2file("dgsa.inp", "   set the file name to null (\"\") -}");
	print2file("dgsa.inp", "{* NOE distance restraints files. *}");
	print2file("dgsa.inp", "{* restraint set 1 file *}");
	print2file("dgsa.inp", "{===>} nmr.noe.file.1=\"contact.tbl\";");
	print2file("dgsa.inp", "{* restraint set 2 file *}");
	print2file("dgsa.inp", "{===>} nmr.noe.file.2=\"\";");
	print2file("dgsa.inp", "{* restraint set 3 file *}");
	print2file("dgsa.inp", "{===>} nmr.noe.file.3=\"\";");
	print2file("dgsa.inp", "{* restraint set 4 file *}");
	print2file("dgsa.inp", "{===>} nmr.noe.file.4=\"\";");
	print2file("dgsa.inp", "{* restraint set 5 file *}");
	print2file("dgsa.inp", "{===>} nmr.noe.file.5=\"\";");
	print2file("dgsa.inp", "{* NOE averaging modes *}");
	print2file("dgsa.inp", "{* restraint set 1 *}");
	print2file("dgsa.inp", "{+ choice: \"sum\" \"cent\" \"R-6\" \"R-3\" \"symm\" +}");
	print2file("dgsa.inp", "{===>} nmr.noe.ave.mode.1=\"cent\";");
	print2file("dgsa.inp", "{* restraint set 2 *}");
	print2file("dgsa.inp", "{+ choice: \"sum\" \"cent\" \"R-6\" \"R-3\" \"symm\" +}");
	print2file("dgsa.inp", "{===>} nmr.noe.ave.mode.2=\"sum\";");
	print2file("dgsa.inp", "{* restraint set 3 *}");
	print2file("dgsa.inp", "{+ choice: \"sum\" \"cent\" \"R-6\" \"R-3\" \"symm\" +}");
	print2file("dgsa.inp", "{===>} nmr.noe.ave.mode.3=\"R-6\";");
	print2file("dgsa.inp", "{* restraint set 4 *}");
	print2file("dgsa.inp", "{+ choice: \"sum\" \"cent\" \"R-6\" \"R-3\" \"symm\" +}");
	print2file("dgsa.inp", "{===>} nmr.noe.ave.mode.4=\"\";");
	print2file("dgsa.inp", "{* restraint set 5 *}");
	print2file("dgsa.inp", "{+ choice: \"sum\" \"cent\" \"R-6\" \"R-3\" \"symm\" +}");
	print2file("dgsa.inp", "{===>} nmr.noe.ave.mode.5=\"\";");
	print2file("dgsa.inp", "{======================== hydrogen bond data ==========================}");
	print2file("dgsa.inp", "{* hydrogen-bond distance restraints file. *}");
	print2file("dgsa.inp", "{===>} nmr.noe.hbnd.file=\"\";");
	print2file("dgsa.inp", "{* enter hydrogen-bond distance averaging mode *}");
	print2file("dgsa.inp", "{+ choice: \"sum\" \"cent\" \"R-6\" \"R-3\" \"symm\" +}");
	print2file("dgsa.inp", "{===>} nmr.noe.ave.mode.hbnd=\"cent\";");
	print2file("dgsa.inp", "{======================= 3-bond J-coupling data =======================}");
	print2file("dgsa.inp", "{* the default setup is for the phi dihedral *}");
	print2file("dgsa.inp", "{* Class 1 *}");
	print2file("dgsa.inp", "{* 3-bond J-coupling non-glycine restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.file.1=\"\";");
	print2file("dgsa.inp", "{* 3-bond J-coupling non-glycine potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" \"multiple\" +}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.pot.1=\"harmonic\";");
	print2file("dgsa.inp", "{* 3-bond J-coupling non-glycine force value *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.force.1.1=1;");
	print2file("dgsa.inp", "{* 3-bond j-coupling multiple class force second value *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.force.2.1=0;");
	print2file("dgsa.inp", "{* 3-bond j-coupling Karplus coefficients *}");
	print2file("dgsa.inp", "{* the default values are for phi *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.1.1=6.98;");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.2.1=-1.38;");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.3.1=1.72;");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.4.1=-60.0;");
	print2file("dgsa.inp", "{* Class 2 *}");
	print2file("dgsa.inp", "{* 3-bond j-coupling glycine restraints files *}");
	print2file("dgsa.inp", "{* The potential for the glycine class must be multiple *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.file.2=\"\";");
	print2file("dgsa.inp", "{* 3-bond J-coupling non-glycine potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" \"multiple\" +}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.pot.2=\"multiple\";");
	print2file("dgsa.inp", "{* 3-bond J-coupling first force value *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.force.1.2=1;");
	print2file("dgsa.inp", "{* 3-bond j-coupling glycine or multiple force second value *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.force.2.2=0;");
	print2file("dgsa.inp", "{* 3-bond j-coupling Karplus coefficients *}");
	print2file("dgsa.inp", "{* the default values are for glycine phi *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.1.2=6.98;");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.2.2=-1.38;");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.3.2=1.72;");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.4.2=0.0;");
	print2file("dgsa.inp", "{================ 1-bond heteronuclear J-coupling data ================}");
	print2file("dgsa.inp", "{* Class 1 *}");
	print2file("dgsa.inp", "{* 1-bond heteronuclear j-coupling file *}");
	print2file("dgsa.inp", "{===>} nmr.oneb.file.1=\"\";");
	print2file("dgsa.inp", "{* 1-bond heteronuclear j-coupling potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" +}");
	print2file("dgsa.inp", "{===>} nmr.oneb.pot.1=\"harmonic\";");
	print2file("dgsa.inp", "{* 1-bond heteronuclear j-coupling force value *}");
	print2file("dgsa.inp", "{===>} nmr.oneb.force.1=1.0;");
	print2file("dgsa.inp", "{=============== alpha/beta carbon chemical shift data ================}");
	print2file("dgsa.inp", "{* Class 1 *}");
	print2file("dgsa.inp", "{* carbon, alpha and beta, chemical shift restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.carb.file.1=\"\";");
	print2file("dgsa.inp", "{* carbon, alpha and beta, chemical shift restraint potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" +}");
	print2file("dgsa.inp", "{===>} nmr.carb.pot.1=\"harmonic\";");
	print2file("dgsa.inp", "{* carbon, alpha and beta, chemical shift restraint force value *}");
	print2file("dgsa.inp", "{===>} nmr.carb.force.1=0.5;");
	print2file("dgsa.inp", "{===================== proton chemical shift data =====================}");
	print2file("dgsa.inp", "{* Class 1 *}");
	print2file("dgsa.inp", "{* class 1 proton chemical shift restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.prot.file.1=\"\";");
	print2file("dgsa.inp", "{* class 1 proton chemical shift potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" \"multiple\" +}");
	print2file("dgsa.inp", "{===>} nmr.prot.pot.1=\"harmonic\";");
	print2file("dgsa.inp", "{* class 1 proton chemical shift force value *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.1.1=7.5;");
	print2file("dgsa.inp", "{* 2nd class 1 proton chemical shift force value for multi *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.2.1=0;");
	print2file("dgsa.inp", "{* class 1 proton chemical shift violation cutoff threshold *}");
	print2file("dgsa.inp", "{===>} nmr.prot.thresh.1=0.3;");
	print2file("dgsa.inp", "{* Class 2 *}");
	print2file("dgsa.inp", "{* class 2 proton chemical shift restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.prot.file.2=\"\";");
	print2file("dgsa.inp", "{* class 2 proton chemical shift potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" \"multiple\" +}");
	print2file("dgsa.inp", "{===>} nmr.prot.pot.2=\"harmonic\";");
	print2file("dgsa.inp", "{* class 2 proton chemical shift force value *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.1.2=7.5;");
	print2file("dgsa.inp", "{* 2nd class 2 proton chemical shift force value for multi *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.2.2=0;");
	print2file("dgsa.inp", "{* class 2 proton chemical shift violation cutoff threshold *}");
	print2file("dgsa.inp", "{===>} nmr.prot.thresh.2=0.3;");
	print2file("dgsa.inp", "{* Class 3 *}");
	print2file("dgsa.inp", "{* class 3 proton chemical shift restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.prot.file.3=\"\";");
	print2file("dgsa.inp", "{* class 3 proton chemical shift potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" \"multiple\" +}");
	print2file("dgsa.inp", "{===>} nmr.prot.pot.3=\"harmonic\";");
	print2file("dgsa.inp", "{* class 3 proton chemical shift force value *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.1.3=7.5;");
	print2file("dgsa.inp", "{* 2nd class 3 proton chemical shift force value for multi *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.2.3=0;");
	print2file("dgsa.inp", "{* class 3 proton chemical shift violation cutoff threshold *}");
	print2file("dgsa.inp", "{===>} nmr.prot.thresh.3=0.3;");
	print2file("dgsa.inp", "{* Class 4 *}");
	print2file("dgsa.inp", "{* class 4 proton chemical shift restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.prot.file.4=\"\";");
	print2file("dgsa.inp", "{* class 4 proton chemical shift potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" \"multiple\" +}");
	print2file("dgsa.inp", "{===>} nmr.prot.pot.4=\"multiple\";");
	print2file("dgsa.inp", "{* class 4 proton chemical shift force value *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.1.4=7.5;");
	print2file("dgsa.inp", "{* 2nd class 4 proton chemical shift force value for multi *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.2.4=0;");
	print2file("dgsa.inp", "{* class 4 proton chemical shift violation cutoff threshold *}");
	print2file("dgsa.inp", "{===>} nmr.prot.thresh.4=0.3;");
	print2file("dgsa.inp", "{================ diffusion anisotropy restraint data =================}");
	print2file("dgsa.inp", "{* fixed or harmonically restrained external axis *}");
	print2file("dgsa.inp", "{+ choice: \"fixed\" \"harm\" +}");
	print2file("dgsa.inp", "{===>} nmr.dani.axis=\"harm\";");
	print2file("dgsa.inp", "{* Class 1 *}");
	print2file("dgsa.inp", "{* diffusion anisotropy restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.dani.file.1=\"\";");
	print2file("dgsa.inp", "{* diffusion anisotropy potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" +}");
	print2file("dgsa.inp", "{===>} nmr.dani.pot.1=\"harmonic\";");
	print2file("dgsa.inp", "{* diffusion anisotropy initial force value *}");
	print2file("dgsa.inp", "{===>} nmr.dani.force.init.1=0.01;");
	print2file("dgsa.inp", "{* diffusion anisotropy final force value *}");
	print2file("dgsa.inp", "{===>} nmr.dani.force.finl.1=1.0;");
	print2file("dgsa.inp", "{* diffusion anisotropy coefficients *}");
	print2file("dgsa.inp", "{* coef: <Tc> <anis> <rhombicity> <wh> <wn> *}");
	print2file("dgsa.inp", "{* Tc = 1/2(Dx+Dy+Dz) in <ns> *} ");
	print2file("dgsa.inp", "{===>} nmr.dani.coef.1.1=13.1;");
	print2file("dgsa.inp", "{* anis = Dz/0.5*(Dx+Dy) *} ");
	print2file("dgsa.inp", "{===>} nmr.dani.coef.2.1=2.1;");
	print2file("dgsa.inp", "{* rhombicity = 1.5*(Dy-Dx)/(Dz-0.5*(Dy+Dx)) *} ");
	print2file("dgsa.inp", "{===>} nmr.dani.coef.3.1=0.0;");
	print2file("dgsa.inp", "{* wH in <MHz> *} ");
	print2file("dgsa.inp", "{===>} nmr.dani.coef.4.1=600.13;");
	print2file("dgsa.inp", "{* wN in <MHz> *}");
	print2file("dgsa.inp", "{===>} nmr.dani.coef.5.1=60.82;");
	print2file("dgsa.inp", "{============= susceptability anisotropy restraint data ===============}");
	print2file("dgsa.inp", "{* fixed or harmonically restrained external axis *}");
	print2file("dgsa.inp", "{+ choice: \"fixed\" \"harm\" +}");
	print2file("dgsa.inp", "{===>} nmr.sani.axis=\"harm\";");
	print2file("dgsa.inp", "{* Class 1 *}");
	print2file("dgsa.inp", "{* susceptability anisotropy restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.sani.file.1=\"\";");
	print2file("dgsa.inp", "{* susceptability anisotropy potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" +}");
	print2file("dgsa.inp", "{===>} nmr.sani.pot.1=\"harmonic\";");
	print2file("dgsa.inp", "{* susceptability anisotropy initial force value *}");
	print2file("dgsa.inp", "{===>} nmr.sani.force.init.1=0.01;");
	print2file("dgsa.inp", "{* susceptability anisotropy final force value *}");
	print2file("dgsa.inp", "{===>} nmr.sani.force.finl.1=50.0;");
	print2file("dgsa.inp", "{* susceptability anisotropy coefficients *}");
	print2file("dgsa.inp", "{* coef: <DFS> <axial > <rhombicity>;");
	print2file("dgsa.inp", "   a0+a1*(3*cos(theta)^2-1)+a2*(3/2)*sin(theta)^2*cos(2*phi) *}");
	print2file("dgsa.inp", "{* DFS = a0 *}");
	print2file("dgsa.inp", "{===>} nmr.sani.coef.1.1=-0.0601;");
	print2file("dgsa.inp", "{* axial = a0-a1-3/2*a2 *}");
	print2file("dgsa.inp", "{===>} nmr.sani.coef.2.1=-8.02;");
	print2file("dgsa.inp", "{* rhombicity = a2/a1 *}");
	print2file("dgsa.inp", "{===>} nmr.sani.coef.3.1=0.4;");
	print2file("dgsa.inp", "{======================== other restraint data ========================}");
	print2file("dgsa.inp", "{* dihedral angle restraints file *}");
	print2file("dgsa.inp", "{* Note: the restraint file MUST NOT contain restraints ");
	print2file("dgsa.inp", "         dihedral or end *}");
	print2file("dgsa.inp", "{===>} nmr.cdih.file=\"\";");
	print2file("dgsa.inp", "{* DNA-RNA base planarity restraints file *}");
	print2file("dgsa.inp", "{* Note: include weights as \$pscale in the restraint file *}");
	print2file("dgsa.inp", "{===>} nmr.plan.file=\"\";");
	print2file("dgsa.inp", "{* input planarity scale factor - this will be written into \$pscale *}");
	print2file("dgsa.inp", "{===>} nmr.plan.scale=150;");
	print2file("dgsa.inp", "{* NCS-restraints file *}");
	print2file("dgsa.inp", "{* example is in inputs/xtal_data/eg1_ncs_restrain.dat *}");
	print2file("dgsa.inp", "{===>} nmr.ncs.file=\"\";");
	print2file("dgsa.inp", "{======================== input/output files ==========================}");
	print2file("dgsa.inp", "{* base name for input coordinate files *}");
	print2file("dgsa.inp", "{* used for simulated annealing when distance geometry is not used *}");
	print2file("dgsa.inp", "{===>} pdb.in.name=\"dg_sub_embed\";");
	print2file("dgsa.inp", "{* base name for output coordinate files *}");
	print2file("dgsa.inp", "{===>} pdb.out.name=\"$ID\";");
	print2file("dgsa.inp", "{===========================================================================}");
	print2file("dgsa.inp", "{         things below this line do not normally need to be changed         }");
	print2file("dgsa.inp", "{         except for the torsion angle topology setup if you have           }");
	print2file("dgsa.inp", "{         molecules other than protein or nucleic acid                      }");
	print2file("dgsa.inp", "{===========================================================================}");
	print2file("dgsa.inp", "flg.cv.flag=false;");
	print2file("dgsa.inp", "flg.cv.noe=false;");
	print2file("dgsa.inp", "flg.cv.coup=false;");
	print2file("dgsa.inp", "flg.cv.cdih=false;");
	print2file("dgsa.inp", "nmr.cv.numpart=10;");
	print2file("dgsa.inp", " ) {- end block parameter definition -}");
	print2file("dgsa.inp", "checkversion 1.3");
	print2file("dgsa.inp", "evaluate (\$log_level=quiet)");
	print2file("dgsa.inp", "structure ");
	print2file("dgsa.inp", "   if  (&struct.1 # \"\") then");
	print2file("dgsa.inp", "      \@\@&struct.1 ");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if  (&struct.2 # \"\") then");
	print2file("dgsa.inp", "      \@\@&struct.2 ");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if  (&struct.3 # \"\") then");
	print2file("dgsa.inp", "      \@\@&struct.3 ");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if  (&struct.4 # \"\") then");
	print2file("dgsa.inp", "      \@\@&struct.4 ");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if  (&struct.5 # \"\") then");
	print2file("dgsa.inp", "      \@\@&struct.5 ");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "end");
	print2file("dgsa.inp", "if ( &BLANK%pdb.in.file.1 = false ) then");
	print2file("dgsa.inp", "   coor \@\@&pdb.in.file.1");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "if ( &BLANK%pdb.in.file.2 = false ) then");
	print2file("dgsa.inp", "   coor \@\@&pdb.in.file.2");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "if ( &BLANK%pdb.in.file.3 = false ) then");
	print2file("dgsa.inp", "   coor \@\@&pdb.in.file.3");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "parameter");
	print2file("dgsa.inp", "   if (&par.1 # \"\") then");
	print2file("dgsa.inp", "      \@\@&par.1");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if (&par.2 # \"\") then");
	print2file("dgsa.inp", "      \@\@&par.2");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if (&par.3 # \"\") then");
	print2file("dgsa.inp", "      \@\@&par.3");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if (&par.4 # \"\") then");
	print2file("dgsa.inp", "      \@\@&par.4");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if (&par.5 # \"\") then");
	print2file("dgsa.inp", "      \@\@&par.5");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "end");
	print2file("dgsa.inp", "if ( \$log_level = verbose ) then");
	print2file("dgsa.inp", "  set message=normal echo=on end");
	print2file("dgsa.inp", "else");
	print2file("dgsa.inp", "  set message=off echo=off end");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "parameter");
	print2file("dgsa.inp", "   nbonds");
	print2file("dgsa.inp", "      repel=0.5");
	print2file("dgsa.inp", "      rexp=2 irexp=2 rcon=1.");
	print2file("dgsa.inp", "      nbxmod=-2");
	print2file("dgsa.inp", "      wmin=0.01");
	print2file("dgsa.inp", "      cutnb=4.5 ctonnb=2.99 ctofnb=3.");
	print2file("dgsa.inp", "      tolerance=0.5");
	print2file("dgsa.inp", "   end");
	print2file("dgsa.inp", "end");
	print2file("dgsa.inp", "set seed=&md.seed end");
	print2file("dgsa.inp", "{- Read experimental data -}");
	print2file("dgsa.inp", "   \@CNS_NMRMODULE:readdata ( nmr=&nmr;");
	print2file("dgsa.inp", "                             flag=&flg;");
	print2file("dgsa.inp", "                             output=\$nmr; )");
	print2file("dgsa.inp", "{- Read and store the number of NMR restraints -}");
	print2file("dgsa.inp", "   \@CNS_NMRMODULE:restraintnumber ( num=\$num; )");
	print2file("dgsa.inp", "   ");
	print2file("dgsa.inp", "{- Set mass and parameter values -}");
	print2file("dgsa.inp", "   ");
	print2file("dgsa.inp", "do (fbeta=10) (all)");
	print2file("dgsa.inp", "do (mass=100) (all)");
	print2file("dgsa.inp", "parameter                  ");
	print2file("dgsa.inp", "   nbonds  ");
	print2file("dgsa.inp", "      repel=0.80  ");
	print2file("dgsa.inp", "      rexp=2 irexp=2 rcon=1. ");
	print2file("dgsa.inp", "      nbxmod=3  ");
	print2file("dgsa.inp", "      wmin=0.01  ");
	print2file("dgsa.inp", "      cutnb=6.0 ctonnb=2.99 ctofnb=3.  ");
	print2file("dgsa.inp", "      tolerance=1.5  ");
	print2file("dgsa.inp", "   end  ");
	print2file("dgsa.inp", "end");
	print2file("dgsa.inp", "evaluate (\$nmr.trial.count = 0)    {- Initialize current structure number   -}");
	print2file("dgsa.inp", "evaluate (\$nmr.accept.count = 0)    {- Initialize number accepted            -}");
	print2file("dgsa.inp", "evaluate (\$nmr.counter 	= 0)");
	print2file("dgsa.inp", "evaluate (\$nmr.prev.counter = -1)");
	print2file("dgsa.inp", "\@CNS_NMRMODULE:initave  ( ave=\$ave;");
	print2file("dgsa.inp", "                          ave2=\$ave2;");
	print2file("dgsa.inp", "                          cv=\$cv;");
	print2file("dgsa.inp", "                          ener1=\$ener1;");
	print2file("dgsa.inp", "                          ener2=\$ener2;");
	print2file("dgsa.inp", "                          flag=&flg;");
	print2file("dgsa.inp", "                          nmr.prot=&nmr.prot; )");
	print2file("dgsa.inp", "        ");
	print2file("dgsa.inp", "{- Zero the force constant of disulfide bonds. -}");
	print2file("dgsa.inp", "parameter");
	print2file("dgsa.inp", "   bonds ( name SG ) ( name SG ) 0. TOKEN ");
	print2file("dgsa.inp", "end");
	print2file("dgsa.inp", "{- define a distance restraints for each disulfide bond, i.e., ");
	print2file("dgsa.inp", "   treat it as if it were an NOE and break the bond. -}");
	print2file("dgsa.inp", "for \$ss_rm_id_1 in id ( name SG ) loop STRM");
	print2file("dgsa.inp", "  for \$ss_rm_id_2 in id ( name SG and ");
	print2file("dgsa.inp", "			  bondedto ( id \$ss_rm_id_1 )  ) loop STR2");
	print2file("dgsa.inp", "    if (\$ss_rm_id_1 > \$ss_rm_id_2) then");
	print2file("dgsa.inp", "      pick bond ( id \$ss_rm_id_1 ) ( id \$ss_rm_id_2 ) equil");
	print2file("dgsa.inp", "      evaluate (\$ss_bond=\$result) ");
	print2file("dgsa.inp", "      noe ");
	print2file("dgsa.inp", "         assign ( id \$ss_rm_id_1 ) ( id \$ss_rm_id_2 ) \$ss_bond 0.1 0.1");
	print2file("dgsa.inp", "      end ");
	print2file("dgsa.inp", "    end if");
	print2file("dgsa.inp", "  end loop STR2");
	print2file("dgsa.inp", "end loop STRM");
	print2file("dgsa.inp", "{- Count the number of residues and determine molecule type -}");
	print2file("dgsa.inp", "identify (store9) (tag)");
	print2file("dgsa.inp", "evaluate (\$nmr.rsn.num = \$SELECT)");
	print2file("dgsa.inp", "identify (store9) ( tag and ( resn THY or resn CYT or resn GUA or");
	print2file("dgsa.inp", "                              resn ADE or resn URI ))");
	print2file("dgsa.inp", "evaluate (\$nmr.nucl.num = \$SELECT)    ");
	print2file("dgsa.inp", "if ( &md.dg.ref = \"coordinate\" ) then");
	print2file("dgsa.inp", "   flag exclude * include bond angl impr vdw end ");
	print2file("dgsa.inp", "   minimize lbfgs nstep=2000 drop=10.  nprint=100 end");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "do (refx=x) ( all )");
	print2file("dgsa.inp", "do (refy=y) ( all )");
	print2file("dgsa.inp", "do (refz=z) ( all )");
	print2file("dgsa.inp", "{- generate and store a bounds matrix -}");
	print2file("dgsa.inp", "if (&flg.dg.flag=true) then");
	print2file("dgsa.inp", "   flags exclude * include bond angle dihedral improper vdw noe cdih end");
	print2file("dgsa.inp", "   mmdg");
	print2file("dgsa.inp", "      shortest-path-algorithm=&&md.dg.algo");
	print2file("dgsa.inp", "      scale=&md.dg.scale");
	print2file("dgsa.inp", "      exponent=&md.dg.expo");
	print2file("dgsa.inp", "      baccuracy=&md.dg.bacc");
	print2file("dgsa.inp", "      taccuracy=&md.dg.tacc");
	print2file("dgsa.inp", "      iaccuracy=&md.dg.iacc");
	print2file("dgsa.inp", "      paccuracy=&md.dg.pacc");
	print2file("dgsa.inp", "      if (&md.dg.type=\"sub\") then");
	print2file("dgsa.inp", "   	 reference=&&md.dg.ref");
	print2file("dgsa.inp", "   	 storebounds");
	print2file("dgsa.inp", "      else");
	print2file("dgsa.inp", "   	 reference=&&md.dg.ref");
	print2file("dgsa.inp", "   	 group &md.dg.group.slct &md.dg.group.err");
	print2file("dgsa.inp", "   	 storebounds");
	print2file("dgsa.inp", "      end if");
	print2file("dgsa.inp", "   end");
	print2file("dgsa.inp", "      ");
	print2file("dgsa.inp", "   {- Begin protocol to generate structures distance geometry structures -}");
	print2file("dgsa.inp", "   while (&pdb.end.count > \$nmr.counter) loop dg");
	print2file("dgsa.inp", "      evaluate (\$nmr.counter=\$nmr.counter + 1)");
	print2file("dgsa.inp", "      evaluate (\$embedded=false)");
	print2file("dgsa.inp", "      flags exclude * include dg end");
	print2file("dgsa.inp", "      if (&md.dg.type=\"sub\") then");
	print2file("dgsa.inp", "   	 igroup interaction=(&md.dg.select) (&md.dg.select) end");
	print2file("dgsa.inp", "      end if");
	print2file("dgsa.inp", "      coor init end");
	print2file("dgsa.inp", "      while (\$embedded = false) loop embed");
	print2file("dgsa.inp", "   	 mmdg");
	print2file("dgsa.inp", "   	    if (&md.dg.type=\"sub\") then");
	print2file("dgsa.inp", "   	       recallbounds");
	print2file("dgsa.inp", "   	       substructure=(&md.dg.select)");
	print2file("dgsa.inp", "   	       selection=(&md.dg.select)");
	print2file("dgsa.inp", "   	    else");
	print2file("dgsa.inp", "   	       recallbounds");
	print2file("dgsa.inp", "   	       selection=(all)");
	print2file("dgsa.inp", "   	       if (&md.dg.metr.flag=true) then");
	print2file("dgsa.inp", "   		  &&md.dg.ord");
	print2file("dgsa.inp", "   		  metrization=(&md.dg.metr.atom)=&md.dg.metr.num");
	print2file("dgsa.inp", "   	       end if");
	print2file("dgsa.inp", "   	    end if");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "      end loop embed");
	print2file("dgsa.inp", "      do (x = x * \$dgscale) (known)");
	print2file("dgsa.inp", "      do (y = y * \$dgscale) (known)");
	print2file("dgsa.inp", "      do (z = z * \$dgscale) (known)");
	print2file("dgsa.inp", "      minimize lbfgs");
	print2file("dgsa.inp", "   	 nstep=&md.dg.step drop=1. nprint=25");
	print2file("dgsa.inp", "      end");
	print2file("dgsa.inp", "      \@CNS_NMRMODULE:printdg ( md=&md;");
	print2file("dgsa.inp", "                               output=\$nmr;");
	print2file("dgsa.inp", "                               pdb=&pdb; )");
	print2file("dgsa.inp", "   end loop dg");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "{- initialize and set scaling factors for simulated annealing -}");
	print2file("dgsa.inp", "set seed=&md.seed end");
	print2file("dgsa.inp", "evaluate (\$nmr.trial.count = 0)    {- Initialize current structure number   -}");
	print2file("dgsa.inp", "evaluate (\$nmr.dg.count = 0)");
	print2file("dgsa.inp", "evaluate (\$nmr.accept.count = 0)   {- Initialize number accepted            -}");
	print2file("dgsa.inp", "evaluate (\$nmr.counter = 0)");
	print2file("dgsa.inp", "evaluate (\$coor_count_init=0.)");
	print2file("dgsa.inp", "evaluate (\$coor_input_count=0.)");
	print2file("dgsa.inp", "if (&flg.dg.flag=true) then");
	print2file("dgsa.inp", "   evaluate (\$coor_input_count=&pdb.end.count)");
	print2file("dgsa.inp", "else");
	print2file("dgsa.inp", "   evaluate (\$coor_input_count=&pdb.dg.count)");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "\@CNS_NMRMODULE:initave  ( flag=&flg;");
	print2file("dgsa.inp", "                          ave=\$ave;");
	print2file("dgsa.inp", "                          ave2=\$ave2;");
	print2file("dgsa.inp", "                          cv=\$cv;");
	print2file("dgsa.inp", "                          ener1=\$ener1;");
	print2file("dgsa.inp", "                          ener2=\$ener2;");
	print2file("dgsa.inp", "                          nmr.prot=&nmr.prot; )");
	print2file("dgsa.inp", "        ");
	print2file("dgsa.inp", "{- scaling of nmr restraint data during regularization -}");
	print2file("dgsa.inp", "\@CNS_NMRMODULE:scalehot ( md=&md;");
	print2file("dgsa.inp", "                          nmr=&nmr;");
	print2file("dgsa.inp", "                          input.noe.scale=&md.cool.noe;");
	print2file("dgsa.inp", "                          input.cdih.scale=&md.hot.cdih; )");
	print2file("dgsa.inp", "if (&nmr.dani.axis = \"harm\") then");
	print2file("dgsa.inp", "   do (harmonic=20.0) (resid 500 and name OO)");
	print2file("dgsa.inp", "   do (harmonic=0.0) (resid 500 and name Z )");
	print2file("dgsa.inp", "   do (harmonic=0.0) (resid 500 and name X )");
	print2file("dgsa.inp", "   do (harmonic=0.0) (resid 500 and name Y )");
	print2file("dgsa.inp", "   do (harmonic=0.0) (not (resid 500))");
	print2file("dgsa.inp", "   restraints harmonic exponent=2 end");
	print2file("dgsa.inp", "elseif (&nmr.sani.axis = \"harm\") then");
	print2file("dgsa.inp", "   do (harmonic=20.0) (resid 500 and name OO)");
	print2file("dgsa.inp", "   do (harmonic=0.0) (resid 500 and name Z )");
	print2file("dgsa.inp", "   do (harmonic=0.0) (resid 500 and name X )");
	print2file("dgsa.inp", "   do (harmonic=0.0) (resid 500 and name Y )");
	print2file("dgsa.inp", "   do (harmonic=0.0) (not (resid 500))");
	print2file("dgsa.inp", "   restraints harmonic exponent=2 end");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "{- Increase the disulfide bond force constants to their full strength -}");
	print2file("dgsa.inp", "   parameter");
	print2file("dgsa.inp", "      bonds ( name SG ) ( name SG ) 1000. TOKEN ");
	print2file("dgsa.inp", "   end");
	print2file("dgsa.inp", "{- Regularize structures generated by distance geometry - loop until done -}");
	print2file("dgsa.inp", "if (&flg.dgsa.flag=true) then");
	print2file("dgsa.inp", "   while (&pdb.end.count > \$nmr.counter) loop dgsa");
	print2file("dgsa.inp", "      {- Set parameter values -}");
	print2file("dgsa.inp", "      parameter");
	print2file("dgsa.inp", "         nbonds");
	print2file("dgsa.inp", "            repel=0.5");
	print2file("dgsa.inp", "            rexp=2 irexp=2 rcon=1.");
	print2file("dgsa.inp", "            nbxmod=-2");
	print2file("dgsa.inp", "            wmin=0.01");
	print2file("dgsa.inp", "            cutnb=4.5 ctonnb=2.99 ctofnb=3.");
	print2file("dgsa.inp", "            tolerance=0.5");
	print2file("dgsa.inp", "         end");
	print2file("dgsa.inp", "      end");
	print2file("dgsa.inp", "      evaluate (\$nmr.trial.count = \$nmr.trial.count + 1)");
	print2file("dgsa.inp", "      if (\$nmr.trial.count <= \$coor_input_count) then");
	print2file("dgsa.inp", "         evaluate (\$nmr.dg.count=\$nmr.dg.count+1)");
	print2file("dgsa.inp", "         evaluate (\$coor_count_init=0.)");
	print2file("dgsa.inp", "      else");
	print2file("dgsa.inp", "         evaluate (\$coor_count_init=\$coor_count_init+1)");
	print2file("dgsa.inp", "         if (\$coor_count_init > \$coor_input_count ) then");
	print2file("dgsa.inp", "            evaluate (\$coor_count_init=1)");
	print2file("dgsa.inp", "         end if");
	print2file("dgsa.inp", "   	 evaluate (\$nmr.dg.count=\$coor_count_init)");
	print2file("dgsa.inp", "      end if");
	print2file("dgsa.inp", "   {- \$prefix is generated in the macro printdg -}");
	print2file("dgsa.inp", "      if (&flg.dg.flag=true) then");
	print2file("dgsa.inp", "         evaluate (\$filename=\$nmr.prefix+encode(\$nmr.dg.count)+\".pdb\")");
	print2file("dgsa.inp", "      else");
	print2file("dgsa.inp", "         evaluate (\$filename=&pdb.in.name+\"_\"+encode(\$nmr.dg.count)+\".pdb\")");
	print2file("dgsa.inp", "      end if");
	print2file("dgsa.inp", "         ");
	print2file("dgsa.inp", "      {- Test for correct enantiomer -}");
	print2file("dgsa.inp", "      for \$image in ( 1 -1 ) loop imag");
	print2file("dgsa.inp", "         set remarks=reset end ");
	print2file("dgsa.inp", "   	 coor initialize end");
	print2file("dgsa.inp", "   	 coor \@\@\$filename");
	print2file("dgsa.inp", "   	 do (x=x * \$image) ( known )");
	print2file("dgsa.inp", "   	 identity (store1) (not known)");
	print2file("dgsa.inp", "   	 coor copy end");
	print2file("dgsa.inp", "   	 do (x=refx) ( all )");
	print2file("dgsa.inp", "   	 do (y=refy) ( all )");
	print2file("dgsa.inp", "   	 do (z=refz) ( all )");
	print2file("dgsa.inp", "   	 for \$id in id ( tag ) loop fit");
	print2file("dgsa.inp", "   	    coordinates");
	print2file("dgsa.inp", "   	       fit select = ( byresidue (id \$id) and not store1 )");
	print2file("dgsa.inp", "   	    end");
	print2file("dgsa.inp", "   	   coor copy selection=( byresidue (id \$id) ) end");
	print2file("dgsa.inp", "   	 end loop fit");
	print2file("dgsa.inp", "   	 coor swap end");
	print2file("dgsa.inp", "         if (&nmr.dani.axis = \"fixed\" ) then");
	print2file("dgsa.inp", "            fix");
	print2file("dgsa.inp", "               select=(resname ANI)");
	print2file("dgsa.inp", "            end");
	print2file("dgsa.inp", "         elseif (&nmr.sani.axis = \"fixed\" ) then");
	print2file("dgsa.inp", "            fix");
	print2file("dgsa.inp", "               select=(resname ANI)");
	print2file("dgsa.inp", "            end");
	print2file("dgsa.inp", "         end if");
	print2file("dgsa.inp", "   	 parameter");
	print2file("dgsa.inp", "   	    nbonds");
	print2file("dgsa.inp", "   	       nbxmod=-2");
	print2file("dgsa.inp", "   	       repel=0.5");
	print2file("dgsa.inp", "   	    end");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 flags exclude * include bond vdw noe cdih coup oneb ");
	print2file("dgsa.inp", "   	                         carb ncs dani sani harm end");
	print2file("dgsa.inp", "   	 igroup interaction (all) (all) weights * 1.  vdw 20. end end");
	print2file("dgsa.inp", "   	 minimize lbfgs nstep=100 nprint=10 end");
	print2file("dgsa.inp", "   	 flags include angl end");
	print2file("dgsa.inp", "   	 minimize lbfgs nstep=100 nprint=10 end");
	print2file("dgsa.inp", "   	 flags include impr dihe end");
	print2file("dgsa.inp", "   	 evaluate (\$nstep1 = int(&md.hot.step/8))");
	print2file("dgsa.inp", "   	 evaluate (\$nstep2 = int(&md.hot.step/2))");
	print2file("dgsa.inp", "   	 do ( vx = maxwell(0.5) ) ( all )");
	print2file("dgsa.inp", "   	 do ( vy = maxwell(0.5) ) ( all )");
	print2file("dgsa.inp", "   	 do ( vz = maxwell(0.5) ) ( all )");
	print2file("dgsa.inp", "   	 igroup inter (all) (all) weights * 0.1 impr 0.05 vdw 20. end end");
	print2file("dgsa.inp", "   	 dynamics cartesian");
	print2file("dgsa.inp", "   	    cmremove=true");
	print2file("dgsa.inp", "   	    vscaling=false");
	print2file("dgsa.inp", "   	    tcoupling=true");
	print2file("dgsa.inp", "   	    timestep=&md.hot.ss");
	print2file("dgsa.inp", "   	    nstep=\$nstep1");
	print2file("dgsa.inp", "   	    nprint=\$nstep1");
	print2file("dgsa.inp", "   	    temperature=&md.hot.temp");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 igroup inter (all) (all) weights * 0.2 impr 0.1  vdw 20. end end");
	print2file("dgsa.inp", "   	 dynamics cartesian");
	print2file("dgsa.inp", "   	    cmremove=true");
	print2file("dgsa.inp", "   	    vscaling=false");
	print2file("dgsa.inp", "   	    tcoupling=true");
	print2file("dgsa.inp", "   	    timestep=&md.hot.ss");
	print2file("dgsa.inp", "   	    nstep=\$nstep1");
	print2file("dgsa.inp", "   	    nprint=\$nstep1");
	print2file("dgsa.inp", "   	    temperature=&md.hot.temp");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 parameter  nbonds repel=0.9   end  end");
	print2file("dgsa.inp", "   	 igroup inter (all) (all) weights * 0.2 impr 0.2 vdw 0.01 end end");
	print2file("dgsa.inp", "   	 dynamics cartesian");
	print2file("dgsa.inp", "   	    cmremove=true");
	print2file("dgsa.inp", "   	    vscaling=false");
	print2file("dgsa.inp", "   	    tcoupling=true");
	print2file("dgsa.inp", "   	    timestep=&md.hot.ss");
	print2file("dgsa.inp", "   	    nstep=\$nstep1");
	print2file("dgsa.inp", "   	    nprint=\$nstep1");
	print2file("dgsa.inp", "   	    temperature=&md.hot.temp");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 parameter nbonds nbxmod=-3  end  end");
	print2file("dgsa.inp", "   	 igroup inter (all) (all) weights * 0.4 impr 0.4 vdw 0.003 end end");
	print2file("dgsa.inp", "   	 dynamics cartesian");
	print2file("dgsa.inp", "   	    cmremove=true");
	print2file("dgsa.inp", "   	    vscaling=false");
	print2file("dgsa.inp", "   	    tcoupling=true");
	print2file("dgsa.inp", "   	    timestep=&md.hot.ss");
	print2file("dgsa.inp", "   	    nstep=\$nstep2");
	print2file("dgsa.inp", "   	    nprint=\$nstep2");
	print2file("dgsa.inp", "   	    temperature=&md.hot.temp");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 igroup inter (all) (all) weights * 1.0 impr 1.0 vdw 0.003 end end");
	print2file("dgsa.inp", "   	 dynamics cartesian");
	print2file("dgsa.inp", "   	    cmremove=true");
	print2file("dgsa.inp", "   	    vscaling=false");
	print2file("dgsa.inp", "   	    tcoupling=true");
	print2file("dgsa.inp", "   	    timestep=&md.hot.ss");
	print2file("dgsa.inp", "   	    nstep=\$nstep1");
	print2file("dgsa.inp", "   	    nprint=\$nstep1");
	print2file("dgsa.inp", "   	    temperature=&md.hot.temp");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 if (\$image = 1) then");
	print2file("dgsa.inp", "   	    do (store7=x) ( all )");
	print2file("dgsa.inp", "   	    do (store8=y) ( all )");
	print2file("dgsa.inp", "   	    do (store9=z) ( all )");
	print2file("dgsa.inp", "   	    do (store4=vx) ( all )");
	print2file("dgsa.inp", "   	    do (store5=vy) ( all )");
	print2file("dgsa.inp", "   	    do (store6=vz) ( all )");
	print2file("dgsa.inp", "   	 end if");
	print2file("dgsa.inp", "      end loop imag");
	print2file("dgsa.inp", "      {- Establish the correct handedness of the structure -}");
	print2file("dgsa.inp", "      energy end");
	print2file("dgsa.inp", "      evaluate (\$e_minus=\$ener)");
	print2file("dgsa.inp", "      coor copy end");
	print2file("dgsa.inp", "      do (x=store7) ( all )");
	print2file("dgsa.inp", "      do (y=store8) ( all )");
	print2file("dgsa.inp", "      do (z=store9) ( all )");
	print2file("dgsa.inp", "      energy end");
	print2file("dgsa.inp", "      evaluate (\$e_plus=\$ener)");
	print2file("dgsa.inp", "      if ( \$e_plus > \$e_minus ) then");
	print2file("dgsa.inp", "   	 evaluate (\$hand=-1 )");
	print2file("dgsa.inp", "   	 coor swap end");
	print2file("dgsa.inp", "      else");
	print2file("dgsa.inp", "   	 evaluate (\$hand= 1 )");
	print2file("dgsa.inp", "   	 do (vx=store4) ( all )");
	print2file("dgsa.inp", "   	 do (vy=store5) ( all )");
	print2file("dgsa.inp", "   	 do (vz=store6) ( all )");
	print2file("dgsa.inp", "      end if");
	print2file("dgsa.inp", "   {- Slow-cooling with cartesian dynamics -}");
	print2file("dgsa.inp", "      parameter");
	print2file("dgsa.inp", "   	 nbonds");
	print2file("dgsa.inp", "   	    repel=0.80");
	print2file("dgsa.inp", "   	    rexp=2 irexp=2 rcon=1.");
	print2file("dgsa.inp", "   	    nbxmod=3");
	print2file("dgsa.inp", "   	    wmin=0.01");
	print2file("dgsa.inp", "   	    cutnb=6.0 ctonnb=2.99 ctofnb=3.");
	print2file("dgsa.inp", "   	    tolerance=0.5");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "      end");
	print2file("dgsa.inp", "      flags include plan end");
	print2file("dgsa.inp", "      evaluate (\$final_t = 0)");
	print2file("dgsa.inp", "      evaluate (\$ncycle = int((&md.hot.temp-\$final_t)/&md.cool.tmpstp))");
	print2file("dgsa.inp", "      evaluate (\$nstep = int(&md.cool.step/\$ncycle))");
	print2file("dgsa.inp", "      evaluate (\$vdw_step=(&md.cool.vdw.finl/&md.cool.vdw.init)^(1/\$ncycle))");
	print2file("dgsa.inp", "      evaluate (\$rad_step=(&md.cool.init.rad-&md.cool.fina.rad)/\$ncycle)");
	print2file("dgsa.inp", "      evaluate (\$radius=&&md.cool.init.rad)");
	print2file("dgsa.inp", "      {- set up nmr restraint scaling -}");
	print2file("dgsa.inp", "      evaluate (\$kdani.inter.flag=false)");
	print2file("dgsa.inp", "      evaluate (\$ksani.inter.flag=false)");
	print2file("dgsa.inp", "      evaluate (\$kdani.cart.flag=false)");
	print2file("dgsa.inp", "      evaluate (\$ksani.cart.flag=false)");
	print2file("dgsa.inp", "      \@CNS_NMRMODULE:scalecoolsetup ( kdani=\$kdani;");
	print2file("dgsa.inp", "                                      ksani=\$ksani;");
	print2file("dgsa.inp", "                                      nmr=&nmr;");
	print2file("dgsa.inp", "                                      input.noe.scale=&md.cool.noe;");
	print2file("dgsa.inp", "                                      input.cdih.scale=&md.cool.cdih;");
	print2file("dgsa.inp", "                                      input.ncycle=\$ncycle; )");
	print2file("dgsa.inp", "      evaluate (\$bath=&md.hot.temp)");
	print2file("dgsa.inp", "      evaluate (\$k_vdw=&md.cool.vdw.init)");
	print2file("dgsa.inp", "      evaluate (\$i_cool = 0)");
	print2file("dgsa.inp", "      while (\$i_cool <= \$ncycle) loop cool");
	print2file("dgsa.inp", "   	 evaluate (\$i_cool = \$i_cool + 1)");
	print2file("dgsa.inp", "   	 igroup");
	print2file("dgsa.inp", "   	    interaction (chemical h*) (all) weights * 1 vdw 0. elec 0. end");
	print2file("dgsa.inp", "   	    interaction (not chemical h*) (not chemical h*) weights * 1 vdw \$k_vdw end");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 dynamics  cartesian");
	print2file("dgsa.inp", "   	    cmremove=true");
	print2file("dgsa.inp", "   	    vscaling = true");
	print2file("dgsa.inp", "   	    tcoup = false");
	print2file("dgsa.inp", "   	    timestep = &md.cool.ss");
	print2file("dgsa.inp", "   	    nstep = \$nstep");
	print2file("dgsa.inp", "   	    nprint = \$nstep");
	print2file("dgsa.inp", "   	    temperature = \$bath");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 evaluate (\$radius=max(&md.cool.fina.rad,\$radius-\$rad_step))");
	print2file("dgsa.inp", "   	 parameter  nbonds repel=\$radius   end end");
	print2file("dgsa.inp", "   	 evaluate (\$k_vdw=min(&md.cool.vdw.finl,\$k_vdw*\$vdw_step))");
	print2file("dgsa.inp", "   	 evaluate (\$bath=\$bath-&md.cool.tmpstp)");
	print2file("dgsa.inp", "         \@CNS_NMRMODULE:scalecool ( kdani=\$kdani;");
	print2file("dgsa.inp", "                                    ksani=\$ksani;");
	print2file("dgsa.inp", "                                    nmr=&nmr; )");
	print2file("dgsa.inp", "      end loop cool");
	print2file("dgsa.inp", "   {- Final minimization -}");
	print2file("dgsa.inp", "      {- turn on proton chemical shifts -}");
	print2file("dgsa.inp", "      flags include prot end");
	print2file("dgsa.inp", "      ");
	print2file("dgsa.inp", "      if (\$nmr.nucl.num > 0) then");
	print2file("dgsa.inp", "         flags include elec end");
	print2file("dgsa.inp", "      end if");
	print2file("dgsa.inp", "      noe             ");
	print2file("dgsa.inp", "         scale * &md.pow.noe ");
	print2file("dgsa.inp", "      end");
	print2file("dgsa.inp", "        ");
	print2file("dgsa.inp", "      restraints dihedral  ");
	print2file("dgsa.inp", "         scale = &md.pow.cdih  ");
	print2file("dgsa.inp", "      end");
	print2file("dgsa.inp", " 													");
	print2file("dgsa.inp", "      igroup interaction ( all ) ( all ) weights * 1 end end");
	print2file("dgsa.inp", "      evaluate (\$count=0 )");
	print2file("dgsa.inp", "      while (&md.pow.cycl > \$count) loop pmini");
	print2file("dgsa.inp", "         evaluate (\$count=\$count + 1)");
	print2file("dgsa.inp", "         minimize lbfgs nstep=&md.pow.step drop=10.0 nprint=25 end");
	print2file("dgsa.inp", "      end loop pmini");
	print2file("dgsa.inp", "      evaluate (\$nmr.min.num = \$count * &md.pow.step)");
	print2file("dgsa.inp", "      {- translate the geometric center of the structure to the origin -}");
	print2file("dgsa.inp", "      if (\$num.dani > 0. ) then");
	print2file("dgsa.inp", "      elseif (\$num.sani > 0. ) then");
	print2file("dgsa.inp", "      else");
	print2file("dgsa.inp", "         show ave ( x ) ( all )");
	print2file("dgsa.inp", "         evaluate (\$geom_x=-\$result)");
	print2file("dgsa.inp", "         show ave ( y ) ( all )");
	print2file("dgsa.inp", "         evaluate (\$geom_y=-\$result)");
	print2file("dgsa.inp", "         show ave ( z ) ( all )");
	print2file("dgsa.inp", "         evaluate (\$geom_z=-\$result)");
	print2file("dgsa.inp", "         coor translate vector=( \$geom_x \$geom_y \$geom_z ) selection=( all ) end");
	print2file("dgsa.inp", "      end if");
	print2file("dgsa.inp", "      ");
	print2file("dgsa.inp", "      \@CNS_NMRMODULE:printaccept ( ave=\$ave;");
	print2file("dgsa.inp", "                                   ave2=\$ave2;");
	print2file("dgsa.inp", "                                   cv=\$cv;");
	print2file("dgsa.inp", "                                   ener1=\$ener1;");
	print2file("dgsa.inp", "                                   ener2=\$ener2;");
	print2file("dgsa.inp", "                                   flag=&flg;");
	print2file("dgsa.inp", "                                   md=&md;");
	print2file("dgsa.inp", "                                   nmr=&nmr;");
	print2file("dgsa.inp", "                                   num=\$num;");
	print2file("dgsa.inp", "                                   output=\$nmr;");
	print2file("dgsa.inp", "                                   pdb=&pdb;  )");
	print2file("dgsa.inp", "   end loop dgsa");
	print2file("dgsa.inp", "   \@CNS_NMRMODULE:calcave ( ave=\$ave;                 ");
	print2file("dgsa.inp", "                            ave2=\$ave2;               ");
	print2file("dgsa.inp", "                            cv=\$cv;                   ");
	print2file("dgsa.inp", "                            ener1=\$ener1;               ");
	print2file("dgsa.inp", "                            ener2=\$ener2;             ");
	print2file("dgsa.inp", "                            flag=&flg;               ");
	print2file("dgsa.inp", "                            md=&md;");
	print2file("dgsa.inp", "                            nmr=&nmr;");
	print2file("dgsa.inp", "                            num=\$num;                 ");
	print2file("dgsa.inp", "                            output=\$nmr;           ");
	print2file("dgsa.inp", "                            pdb=&pdb;  )");
	print2file("dgsa.inp", "	");
	print2file("dgsa.inp", "      ");
	print2file("dgsa.inp", "      ");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "stop");
}

sub write_cns_generate_seq_file{
	print2file("gseq.inp", "{+ file: generate_seq.inp +}");
	print2file("gseq.inp", "{+ directory: general +}");
	print2file("gseq.inp", "{+ description: Generate structure file for protein, dna/rna, water, ");
	print2file("gseq.inp", "                ligands and/or carbohydrate from sequence information only +}");
	print2file("gseq.inp", "{+ comment: modified by Brian Smith (Edinburgh University) to allow protein");
	print2file("gseq.inp", "            residue renumbering +}");
	print2file("gseq.inp", "{+ authors: Paul Adams, and Axel Brunger +}");
	print2file("gseq.inp", "{+ copyright: Yale University +}");
	print2file("gseq.inp", "{- Guidelines for using this file:");
	print2file("gseq.inp", "   - all strings must be quoted by double-quotes");
	print2file("gseq.inp", "   - logical variables (true/false) are not quoted");
	print2file("gseq.inp", "   - do not remove any evaluate statements from the file -}");
	print2file("gseq.inp", "{- Special patches will have to be entered manually at the relevant points");
	print2file("gseq.inp", "   in the file - see comments throughout the file -}");
	print2file("gseq.inp", "{- begin block parameter definition -} define(");
	print2file("gseq.inp", "{============ protein topology, linkage, and parameter files =============}");
	print2file("gseq.inp", "{* topology files *}");
	print2file("gseq.inp", "{===>} topology_infile_1=\"CNS_TOPPAR:protein.top\";");
	print2file("gseq.inp", "{===>} topology_infile_2=\"CNS_TOPPAR:dna-rna.top\";");
	print2file("gseq.inp", "{===>} topology_infile_3=\"CNS_TOPPAR:water.top\";");
	print2file("gseq.inp", "{===>} topology_infile_4=\"CNS_TOPPAR:ion.top\";");
	print2file("gseq.inp", "{===>} topology_infile_5=\"CNS_TOPPAR:carbohydrate.top\";");
	print2file("gseq.inp", "{===>} topology_infile_6=\"\";");
	print2file("gseq.inp", "{===>} topology_infile_7=\"\";");
	print2file("gseq.inp", "{===>} topology_infile_8=\"\";");
	print2file("gseq.inp", "{* linkage files for linear, continuous polymers (protein, DNA, RNA) *}");
	print2file("gseq.inp", "{===>} link_infile_1=\"CNS_TOPPAR:protein.link\";");
	print2file("gseq.inp", "{===>} link_infile_2=\"CNS_TOPPAR:dna-rna-pho.link\";");
	print2file("gseq.inp", "{===>} link_infile_3=\"\";");
	print2file("gseq.inp", "{* parameter files *}");
	print2file("gseq.inp", "{===>} parameter_infile_1=\"CNS_TOPPAR:protein.param\";");
	print2file("gseq.inp", "{===>} parameter_infile_2=\"CNS_TOPPAR:dna-rna_rep.param\";");
	print2file("gseq.inp", "{===>} parameter_infile_3=\"CNS_TOPPAR:water_rep.param\";");
	print2file("gseq.inp", "{===>} parameter_infile_4=\"CNS_TOPPAR:ion.param\";");
	print2file("gseq.inp", "{===>} parameter_infile_5=\"CNS_TOPPAR:carbohydrate.param\";");
	print2file("gseq.inp", "{===>} parameter_infile_6=\"\";");
	print2file("gseq.inp", "{===>} parameter_infile_7=\"\";");
	print2file("gseq.inp", "{===>} parameter_infile_8=\"\";");
	print2file("gseq.inp", "{====================== other linkages and modifications  ==================}");
	print2file("gseq.inp", "{* extra linkages and modifications by custom patches *}");
	print2file("gseq.inp", "{===>} patch_infile=\"\";");
	print2file("gseq.inp", "{============================= sequence files ==============================}");
	print2file("gseq.inp", "{* multiple sequence files of the same type can be defined by duplicating");
	print2file("gseq.inp", "   the entries below and incrementing the file number *}");
	print2file("gseq.inp", "{* protein sequence file 1 *}");
	print2file("gseq.inp", "{===>} prot_sequence_infile_1=\"input.seq\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} prot_segid_1=\"\";");
	print2file("gseq.inp", "{* start residue numbering at *}");
	print2file("gseq.inp", "{===>} renumber_1=1;");
	print2file("gseq.inp", "{* protein sequence file 2 *}");
	print2file("gseq.inp", "{===>} prot_sequence_infile_2=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} prot_segid_2=\"\";");
	print2file("gseq.inp", "{* start residue numbering at *}");
	print2file("gseq.inp", "{===>} renumber_2=1;");
	print2file("gseq.inp", "{* protein sequence file 3 *}");
	print2file("gseq.inp", "{===>} prot_sequence_infile_3=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} prot_segid_3=\"\";");
	print2file("gseq.inp", "{* start residue numbering at *}");
	print2file("gseq.inp", "{===>} renumber_3=1;");
	print2file("gseq.inp", "{* protein sequence file 4 *}");
	print2file("gseq.inp", "{===>} prot_sequence_infile_4=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} prot_segid_4=\"\";");
	print2file("gseq.inp", "{* start residue numbering at *}");
	print2file("gseq.inp", "{===>} renumber_4=1;");
	print2file("gseq.inp", "{* nucleic acid sequence file 1 *}");
	print2file("gseq.inp", "{===>} nucl_sequence_infile_1=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} nucl_segid_1=\"\";");
	print2file("gseq.inp", "{* nucleic acid sequence file 2 *}");
	print2file("gseq.inp", "{===>} nucl_sequence_infile_2=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} nucl_segid_2=\"\";");
	print2file("gseq.inp", "{* nucleic acid sequence file 3 *}");
	print2file("gseq.inp", "{===>} nucl_sequence_infile_3=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} nucl_segid_3=\"\";");
	print2file("gseq.inp", "{* nucleic acid sequence file 4 *}");
	print2file("gseq.inp", "{===>} nucl_sequence_infile_4=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} nucl_segid_4=\"\";");
	print2file("gseq.inp", "{* water sequence file 1 *}");
	print2file("gseq.inp", "{===>} water_sequence_infile_1=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} water_segid_1=\"\";");
	print2file("gseq.inp", "{* water sequence file 2 *}");
	print2file("gseq.inp", "{===>} water_sequence_infile_2=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} water_segid_2=\"\";");
	print2file("gseq.inp", "{* carbohydrate sequence file 1 *}");
	print2file("gseq.inp", "{===>} carbo_sequence_infile_1=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} carbo_segid_1=\"\";");
	print2file("gseq.inp", "{* carbohydrate sequence file 2 *}");
	print2file("gseq.inp", "{===>} carbo_sequence_infile_2=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} carbo_segid_2=\"\";");
	print2file("gseq.inp", "{* prosthetic group sequence file 1 *}");
	print2file("gseq.inp", "{===>} prost_sequence_infile_1=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} prost_segid_1=\"\";");
	print2file("gseq.inp", "{* prosthetic group sequence file 2 *}");
	print2file("gseq.inp", "{===>} prost_sequence_infile_2=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} prost_segid_2=\"\";");
	print2file("gseq.inp", "{* ligand sequence file 1 *}");
	print2file("gseq.inp", "{===>} lig_sequence_infile_1=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} lig_segid_1=\"\";");
	print2file("gseq.inp", "{* ligand sequence file 2 *}");
	print2file("gseq.inp", "{===>} lig_sequence_infile_2=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} lig_segid_2=\"\";");
	print2file("gseq.inp", "{* ion sequence file 1 *}");
	print2file("gseq.inp", "{===>} ion_sequence_infile_1=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} ion_segid_1=\"\";");
	print2file("gseq.inp", "{* ion sequence file 2 *}");
	print2file("gseq.inp", "{===>} ion_sequence_infile_2=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} ion_segid_2=\"\";");
	print2file("gseq.inp", "{============================= output files ================================}");
	print2file("gseq.inp", "{* output structure file *}");
	print2file("gseq.inp", "{===>} structure_outfile=\"extended.mtf\";");
	print2file("gseq.inp", "{=========================== disulphide bonds ==============================}");
	print2file("gseq.inp", "{* Select pairs of cysteine residues that form disulphide bonds *}");
	print2file("gseq.inp", "{* First 2 entries are the segid and resid of the first cysteine (CYS A). *}");
	print2file("gseq.inp", "{* Second 2 entries are the segid and resid of the second cysteine (CYS B). *}");
	print2file("gseq.inp", "{+ table: rows=8 numbered");
	print2file("gseq.inp", "   cols=5 \"use\" \"segid CYS A\" \"resid CYS A\" \"segid CYS B\" \"resid CYS B\" +}");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_1=true;");
	print2file("gseq.inp", "{===>} ss_i_segid_1=\"\"; ss_i_resid_1=11;");
	print2file("gseq.inp", "{===>} ss_j_segid_1=\"\"; ss_j_resid_1=27;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_2=true;");
	print2file("gseq.inp", "{===>} ss_i_segid_2=\"\"; ss_i_resid_2=45;");
	print2file("gseq.inp", "{===>} ss_j_segid_2=\"\"; ss_j_resid_2=73;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_3=false;");
	print2file("gseq.inp", "{===>} ss_i_segid_3=\"\"; ss_i_resid_3=0;");
	print2file("gseq.inp", "{===>} ss_j_segid_3=\"\"; ss_j_resid_3=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_4=false;");
	print2file("gseq.inp", "{===>} ss_i_segid_4=\"\"; ss_i_resid_4=0;");
	print2file("gseq.inp", "{===>} ss_j_segid_4=\"\"; ss_j_resid_4=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_5=false;");
	print2file("gseq.inp", "{===>} ss_i_segid_5=\"\"; ss_i_resid_5=0;");
	print2file("gseq.inp", "{===>} ss_j_segid_5=\"\"; ss_j_resid_5=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_6=false;");
	print2file("gseq.inp", "{===>} ss_i_segid_6=\"\"; ss_i_resid_6=0;");
	print2file("gseq.inp", "{===>} ss_j_segid_6=\"\"; ss_j_resid_6=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_7=false;");
	print2file("gseq.inp", "{===>} ss_i_segid_7=\"\"; ss_i_resid_7=0;");
	print2file("gseq.inp", "{===>} ss_j_segid_7=\"\"; ss_j_resid_7=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_8=false;");
	print2file("gseq.inp", "{===>} ss_i_segid_8=\"\"; ss_i_resid_8=0;");
	print2file("gseq.inp", "{===>} ss_j_segid_8=\"\"; ss_j_resid_8=0;");
	print2file("gseq.inp", "{=========================== carbohydrate links  ===========================}");
	print2file("gseq.inp", "{* Select pairs of residues that are linked *}");
	print2file("gseq.inp", "{* First entry is the name of the patch residue. *}");
	print2file("gseq.inp", "{* Second and third entries are the resid and segid for the atoms");
	print2file("gseq.inp", "   referenced by \"-\" in the patch. *}");
	print2file("gseq.inp", "{* Fourth and fifth entries are the resid and segid for the atoms");
	print2file("gseq.inp", "   referenced by \"+\" in the patch *}");
	print2file("gseq.inp", "{+ table: rows=6 numbered");
	print2file("gseq.inp", "          cols=6 \"use\" \"patch name\" \"segid -\" \"resid -\" \"segid +\" \"resid +\" +}");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} carbo_use_1=false;");
	print2file("gseq.inp", "{===>} carbo_patch_1=\"B1N\";");
	print2file("gseq.inp", "{===>} carbo_i_segid_1=\"BBBB\"; carbo_i_resid_1=401;");
	print2file("gseq.inp", "{===>} carbo_j_segid_1=\"AAAA\"; carbo_j_resid_1=56;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} carbo_use_2=false;");
	print2file("gseq.inp", "{===>} carbo_patch_2=\"B1N\";");
	print2file("gseq.inp", "{===>} carbo_i_segid_2=\"BBBB\"; carbo_i_resid_2=402;");
	print2file("gseq.inp", "{===>} carbo_j_segid_2=\"AAAA\"; carbo_j_resid_2=182;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} carbo_use_3=false;");
	print2file("gseq.inp", "{===>} carbo_patch_3=\"\";");
	print2file("gseq.inp", "{===>} carbo_i_segid_3=\"\"; carbo_i_resid_3=0;");
	print2file("gseq.inp", "{===>} carbo_j_segid_3=\"\"; carbo_j_resid_3=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} carbo_use_4=false;");
	print2file("gseq.inp", "{===>} carbo_patch_4=\"\";");
	print2file("gseq.inp", "{===>} carbo_i_segid_4=\"\"; carbo_i_resid_4=0;");
	print2file("gseq.inp", "{===>} carbo_j_segid_4=\"\"; carbo_j_resid_4=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} carbo_use_5=false;");
	print2file("gseq.inp", "{===>} carbo_patch_5=\"\";");
	print2file("gseq.inp", "{===>} carbo_i_segid_5=\"\"; carbo_i_resid_5=0;");
	print2file("gseq.inp", "{===>} carbo_j_segid_5=\"\"; carbo_j_resid_5=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} carbo_use_6=false;");
	print2file("gseq.inp", "{===>} carbo_patch_6=\"\";");
	print2file("gseq.inp", "{===>} carbo_i_segid_6=\"\"; carbo_i_resid_6=0;");
	print2file("gseq.inp", "{===>} carbo_j_segid_6=\"\"; carbo_j_resid_6=0;");
	print2file("gseq.inp", "{========================= generate parameters =============================}");
	print2file("gseq.inp", "{* hydrogen flag - determines whether hydrogens will be retained *}");
	print2file("gseq.inp", "{* must be true for NMR, atomic resolution X-ray crystallography ");
	print2file("gseq.inp", "   or modelling.  Set to false for most X-ray crystallographic ");
	print2file("gseq.inp", "   applications at resolution > 1A *}");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} hydrogen_flag=false;");
	print2file("gseq.inp", "{* set bfactor flag *}");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} set_bfactor=true;");
	print2file("gseq.inp", "{* set bfactor value *}");
	print2file("gseq.inp", "{===>} bfactor=15.0;");
	print2file("gseq.inp", "{* set occupancy flag *}");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} set_occupancy=true;");
	print2file("gseq.inp", "{* set occupancy value *}");
	print2file("gseq.inp", "{===>} occupancy=1.0;");
	print2file("gseq.inp", "{===========================================================================}");
	print2file("gseq.inp", "{         things below this line do not need to be changed unless           }");
	print2file("gseq.inp", "{         you need to apply patches - at the appropriate places marked      }");
	print2file("gseq.inp", "{===========================================================================}");
	print2file("gseq.inp", " ) {- end block parameter definition -}");
	print2file("gseq.inp", " checkversion 1.3");
	print2file("gseq.inp", " evaluate (\$log_level=quiet)");
	print2file("gseq.inp", " {- read parameter files -}");
	print2file("gseq.inp", " parameter");
	print2file("gseq.inp", "  evaluate (\$counter=1)");
	print2file("gseq.inp", "  evaluate (\$done=false)");
	print2file("gseq.inp", "  while ( \$done = false ) loop read");
	print2file("gseq.inp", "   if ( &exist_parameter_infile_\$counter = true ) then");
	print2file("gseq.inp", "      if ( &BLANK%parameter_infile_\$counter = false ) then");
	print2file("gseq.inp", "         \@\@&parameter_infile_\$counter");
	print2file("gseq.inp", "      end if");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "    evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", "   evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "  end loop read");
	print2file("gseq.inp", " end");
	print2file("gseq.inp", " {- read topology files -}");
	print2file("gseq.inp", " topology");
	print2file("gseq.inp", "  evaluate (\$counter=1)");
	print2file("gseq.inp", "  evaluate (\$done=false)");
	print2file("gseq.inp", "  while ( \$done = false ) loop read");
	print2file("gseq.inp", "   if ( &exist_topology_infile_\$counter = true ) then");
	print2file("gseq.inp", "      if ( &BLANK%topology_infile_\$counter = false ) then");
	print2file("gseq.inp", "         \@\@&topology_infile_\$counter");
	print2file("gseq.inp", "      end if");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", "   evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "  end loop read");
	print2file("gseq.inp", " end");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop prot");
	print2file("gseq.inp", "   if ( &exist_prot_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%prot_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       do (refx=0) (all)");
	print2file("gseq.inp", "       segment");
	print2file("gseq.inp", "         chain");
	print2file("gseq.inp", "           evaluate (\$count=1)");
	print2file("gseq.inp", "           evaluate (\$done2=false)");
	print2file("gseq.inp", "           while ( \$done2 = false ) loop read");
	print2file("gseq.inp", "             if ( &exist_link_infile_\$count = true ) then");
	print2file("gseq.inp", "               if ( &BLANK%link_infile_\$count = false ) then");
	print2file("gseq.inp", "                  \@\@&link_infile_\$count");
	print2file("gseq.inp", "               end if");
	print2file("gseq.inp", "             else");
	print2file("gseq.inp", "               evaluate (\$done2=true)");
	print2file("gseq.inp", "             end if");
	print2file("gseq.inp", "             evaluate (\$count=\$count+1)");
	print2file("gseq.inp", "           end loop read");
	print2file("gseq.inp", "           sequence \@\@&prot_sequence_infile_\$counter end");
	print2file("gseq.inp", "         end");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "       do (segid=\"T^\" + encode(\$counter)) (attr refx=9999)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     if ( &exist_renumber_\$counter = true ) then");
	print2file("gseq.inp", "         if ( &BLANK%renumber_\$counter = false ) then");
	print2file("gseq.inp", "           evaluate (\$segid=\"T^\" + encode(\$counter))");
	print2file("gseq.inp", "           do ( resid = adjustl(format(\"I4\",decode(resid) + &renumber_\$counter - 1))) ");
	print2file("gseq.inp", "              ( (attr refx=9999) and segid \$segid )");
	print2file("gseq.inp", "         end if");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop prot");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop nseg");
	print2file("gseq.inp", "   if ( &exist_prot_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%prot_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       evaluate (\$segtmp=\"T^\" + encode(\$counter))");
	print2file("gseq.inp", "       do (segid=capitalize(&prot_segid_\$counter)) (segid \$segtmp)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop nseg");
	print2file("gseq.inp", " evaluate (\$ssc=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop ssbr");
	print2file("gseq.inp", "   if ( &exist_ss_use_\$ssc = true ) then");
	print2file("gseq.inp", "     if ( &ss_use_\$ssc = true ) then");
	print2file("gseq.inp", "       evaluate (\$segidtmp1=capitalize(&ss_i_segid_\$ssc))");
	print2file("gseq.inp", "       evaluate (\$segidtmp2=capitalize(&ss_j_segid_\$ssc))");
	print2file("gseq.inp", "       patch disu");
	print2file("gseq.inp", "         reference=1=(segid \$QUOTE%segidtmp1 and resid &ss_i_resid_\$ssc)");
	print2file("gseq.inp", "         reference=2=(segid \$QUOTE%segidtmp2 and resid &ss_j_resid_\$ssc)");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$ssc=\$ssc+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop ssbr");
	print2file("gseq.inp", " {* any special protein patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop nucl");
	print2file("gseq.inp", "   if ( &exist_nucl_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%nucl_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       do (refx=0) (all)");
	print2file("gseq.inp", "       segment");
	print2file("gseq.inp", "         chain");
	print2file("gseq.inp", "           evaluate (\$count=1)");
	print2file("gseq.inp", "           evaluate (\$done2=false)");
	print2file("gseq.inp", "           while ( \$done2 = false ) loop read");
	print2file("gseq.inp", "             if ( &exist_link_infile_\$count = true ) then");
	print2file("gseq.inp", "               if ( &BLANK%link_infile_\$count = false ) then");
	print2file("gseq.inp", "                  \@\@&link_infile_\$count");
	print2file("gseq.inp", "               end if");
	print2file("gseq.inp", "             else");
	print2file("gseq.inp", "               evaluate (\$done2=true)");
	print2file("gseq.inp", "             end if");
	print2file("gseq.inp", "             evaluate (\$count=\$count+1)");
	print2file("gseq.inp", "           end loop read");
	print2file("gseq.inp", "           sequence \@\@&nucl_sequence_infile_\$counter end");
	print2file("gseq.inp", "         end");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "       do (segid=capitalize(&nucl_segid_\$counter)) (attr refx=9999)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop nucl");
	print2file("gseq.inp", " {* patch rna sugars to dna here if needed - select the residues *}");
	print2file("gseq.inp", " {===>} ");
	print2file("gseq.inp", " for \$resid in () loop dna");
	print2file("gseq.inp", "   patch deox reference=nil=(resid \$resid) end");
	print2file("gseq.inp", " end loop dna");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " {* any special nucleic acid patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop carbo");
	print2file("gseq.inp", "   if ( &exist_carbo_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%carbo_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       do (refx=0) (all)");
	print2file("gseq.inp", "       segment");
	print2file("gseq.inp", "         chain");
	print2file("gseq.inp", "           sequence \@\@&carbo_sequence_infile_\$counter end");
	print2file("gseq.inp", "         end");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "       do (segid=capitalize(&carbo_segid_\$counter)) (attr refx=9999)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop carbo");
	print2file("gseq.inp", " evaluate (\$carc=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop cabr");
	print2file("gseq.inp", "   if ( &exist_carbo_use_\$carc = true ) then");
	print2file("gseq.inp", "     if ( &carbo_use_\$carc = true ) then");
	print2file("gseq.inp", "       evaluate (\$segidtmp1=capitalize(&carbo_i_segid_\$carc))");
	print2file("gseq.inp", "       evaluate (\$segidtmp2=capitalize(&carbo_j_segid_\$carc))");
	print2file("gseq.inp", "       patch &carbo_patch_\$carc");
	print2file("gseq.inp", "         reference=-=(segid \$QUOTE%segidtmp1 and");
	print2file("gseq.inp", "                      resid &carbo_i_resid_\$carc)");
	print2file("gseq.inp", "         reference=+=(segid \$QUOTE%segidtmp2 and");
	print2file("gseq.inp", "                      resid &carbo_j_resid_\$carc)");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$carc=\$carc+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop cabr");
	print2file("gseq.inp", " {* any special carbohydrate patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop prost");
	print2file("gseq.inp", "   if ( &exist_prost_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%prost_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       do (refx=0) (all)");
	print2file("gseq.inp", "       segment");
	print2file("gseq.inp", "         chain");
	print2file("gseq.inp", "           sequence \@\@&prost_sequence_infile_\$counter end");
	print2file("gseq.inp", "         end");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "       do (segid=capitalize(&prost_segid_\$counter)) (attr refx=9999)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop prost");
	print2file("gseq.inp", " {* any special prosthetic group patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop liga");
	print2file("gseq.inp", "   if ( &exist_lig_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%lig_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       do (refx=0) (all)");
	print2file("gseq.inp", "       segment");
	print2file("gseq.inp", "         chain");
	print2file("gseq.inp", "           sequence \@\@&lig_sequence_infile_\$counter end");
	print2file("gseq.inp", "         end");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "       do (segid=capitalize(&lig_segid_\$counter)) (attr refx=9999)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop liga");
	print2file("gseq.inp", " {* any special ligand patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop ion");
	print2file("gseq.inp", "   if ( &exist_ion_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%ion_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       do (refx=0) (all)");
	print2file("gseq.inp", "       segment");
	print2file("gseq.inp", "         chain");
	print2file("gseq.inp", "           sequence \@\@&ion_sequence_infile_\$counter end");
	print2file("gseq.inp", "         end");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "       do (segid=capitalize(&ion_segid_\$counter)) (attr refx=9999)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop ion");
	print2file("gseq.inp", " {* any special ion patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop water");
	print2file("gseq.inp", "   if ( &exist_water_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%water_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       do (refx=0) (all)");
	print2file("gseq.inp", "       segment");
	print2file("gseq.inp", "         chain");
	print2file("gseq.inp", "           sequence \@\@&water_sequence_infile_\$counter end");
	print2file("gseq.inp", "         end");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "       do (segid=capitalize(&water_segid_\$counter)) (attr refx=9999)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop water");
	print2file("gseq.inp", " {* any special water patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " {* any final patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " if (&hydrogen_flag=false) then");
	#	print2file("gseq.inp", "   delete selection=( hydrogen ) end");
	print2file("gseq.inp", " end if");
	print2file("gseq.inp", " if (&set_bfactor=true) then");
	print2file("gseq.inp", "   do (b=&bfactor) ( all )");
	print2file("gseq.inp", " end if");
	print2file("gseq.inp", " if (&set_occupancy=true) then");
	print2file("gseq.inp", "   do (q=&occupancy) ( all )");
	print2file("gseq.inp", " end if");
	print2file("gseq.inp", " write structure output=&structure_outfile end");
	print2file("gseq.inp", " stop");
}

sub write_cns_generate_extended_file{
	print2file("extn.inp", "{+ file: generate_extended.inp +}");
	print2file("extn.inp", "{+ directory: nmr_calc +}");
	print2file("extn.inp", "{+ description: Generates an extended strand with ideal geometry ");
	print2file("extn.inp", "                for each connected polymer.  ");
	print2file("extn.inp", "                The molecular structure file must not contain any ");
	print2file("extn.inp", "                closed loops except disulfide bonds which are automatically");
	print2file("extn.inp", "                excluded from the generation of the strand conformation.  +}");
	print2file("extn.inp", "{+ authors: Axel T. Brunger +}");
	print2file("extn.inp", "{+ copyright: Yale University +}");
	print2file("extn.inp", "{- begin block parameter definition -} define(");
	print2file("extn.inp", "{======================= molecular structure =========================}");
	print2file("extn.inp", "{* structure file(s) *}");
	print2file("extn.inp", "{===>} structure_file=\"extended.mtf\";");
	print2file("extn.inp", "{* parameter file(s) *}");
	# CNS_TOPPAR:protein-allhdg5-4.param
	print2file("extn.inp", "{===>} par_1=\"CNS_TOPPAR:protein.param\";");
	print2file("extn.inp", "{===>} par_2=\"\";");
	print2file("extn.inp", "{===>} par_3=\"\";");
	print2file("extn.inp", "{===>} par_4=\"\";");
	print2file("extn.inp", "{===>} par_5=\"\";");
	print2file("extn.inp", "{======================= input parameters ============================}");
	print2file("extn.inp", "{* maximum number of trials to generate an acceptable structure *}");
	print2file("extn.inp", "{===>} max_trial=10;");
	print2file("extn.inp", "{=========================== output files ============================}");
	print2file("extn.inp", "{* output coordinates *}");
	print2file("extn.inp", "{===>} output_coor=\"extended.pdb\";");
	print2file("extn.inp", "                                  ");
	print2file("extn.inp", "{===========================================================================}");
	print2file("extn.inp", "{        things below this line do not normally need to be changed          }");
	print2file("extn.inp", "{===========================================================================}");
	print2file("extn.inp", " ) {- end block parameter definition -}");
	print2file("extn.inp", " checkversion 1.3");
	print2file("extn.inp", " evaluate (\$log_level=quiet)");
	print2file("extn.inp", " structure \@&structure_file end");
	print2file("extn.inp", " parameter");
	print2file("extn.inp", "   if (&par_1 # \" \") then");
	print2file("extn.inp", "      \@\@&par_1");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "   if (&par_2 # \" \") then");
	print2file("extn.inp", "      \@\@&par_2");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "   if (&par_3 # \" \") then");
	print2file("extn.inp", "      \@\@&par_3");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "   if (&par_4 # \" \") then");
	print2file("extn.inp", "      \@\@&par_4");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "   if (&par_5 # \" \") then");
	print2file("extn.inp", "      \@\@&par_5");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", " end");
	print2file("extn.inp", "{ Set force constants for S-S bond lengths and angles to zero  }");
	print2file("extn.inp", "parameter");
	print2file("extn.inp", "   bonds ( name SG ) ( name SG ) 0. 1. ");
	print2file("extn.inp", "end");
	print2file("extn.inp", "igroup interaction=(all) (all) end");
	print2file("extn.inp", "ident (x) ( all )");
	print2file("extn.inp", "do (x=x/5.) ( all )");
	print2file("extn.inp", "do (y=random(0.5) ) ( all )");
	print2file("extn.inp", "do (z=random(0.5) ) ( all )");
	print2file("extn.inp", "flags exclude * include bond angle impr dihe vdw end");
	print2file("extn.inp", "parameter");
	print2file("extn.inp", "   nbonds");
	print2file("extn.inp", "      rcon=50. nbxmod=-3 repel=0.8 cutnb=6. ");
	print2file("extn.inp", "      rexp=2 irexp=2 inhibit=0.0 wmin=0.1 tolerance=0.5");
	print2file("extn.inp", "   end");
	print2file("extn.inp", "end");
	print2file("extn.inp", "evaluate (\$count=1) ");
	print2file("extn.inp", "while (\$count < 10 ) loop l1");
	print2file("extn.inp", "   do (x=x+gauss(0.1)) ( all ) ");
	print2file("extn.inp", "   do (y=y+gauss(0.1)) ( all ) ");
	print2file("extn.inp", "   do (z=z+gauss(0.1)) ( all ) ");
	print2file("extn.inp", "   minimize lbfgs nstep=200 nprint=10 end");
	print2file("extn.inp", "   evaluate (\$count=\$count+1)");
	print2file("extn.inp", "end loop l1");
	print2file("extn.inp", "evaluate (\$accept=false) ");
	print2file("extn.inp", "evaluate (\$trial=1) ");
	print2file("extn.inp", "while (\$accept=false) loop accp");
	print2file("extn.inp", "   for \$1 in id ( tag ) loop resi");
	print2file("extn.inp", "      igroup ");
	print2file("extn.inp", "         interaction=( byresidue (id \$1 ) and not name SG ) ");
	print2file("extn.inp", "                     ( not name SG ) ");
	print2file("extn.inp", "      end");
	print2file("extn.inp", "      evaluate (\$accept=true) ");
	print2file("extn.inp", "      print thres=0.1 bonds");
	print2file("extn.inp", "      if (\$violations > 0) then");
	print2file("extn.inp", "         evaluate (\$accept=false) ");
	print2file("extn.inp", "      end if");
	print2file("extn.inp", "      print thres=10. angles ");
	print2file("extn.inp", "      evaluate (\$angles=\$result)");
	print2file("extn.inp", "      if (\$violations > 0) then");
	print2file("extn.inp", "         evaluate (\$accept=false) ");
	print2file("extn.inp", "      end if");
	print2file("extn.inp", "      print thres=10. improper");
	print2file("extn.inp", "      if (\$violations > 0) then");
	print2file("extn.inp", "         evaluate (\$accept=false) ");
	print2file("extn.inp", "      end if");
	print2file("extn.inp", "      if (\$accept=false) then");
	print2file("extn.inp", "         do (x=x+gauss(0.3)) ( byresidue (id \$1 ) ) ");
	print2file("extn.inp", "         do (y=y+gauss(0.3)) ( byresidue (id \$1 ) ) ");
	print2file("extn.inp", "         do (z=z+gauss(0.3)) ( byresidue (id \$1 ) ) ");
	print2file("extn.inp", "      end if");
	print2file("extn.inp", "   end loop resi");
	print2file("extn.inp", "   igroup interaction=( all ) ( all ) end");
	print2file("extn.inp", "   parameter");
	print2file("extn.inp", "      nbonds");
	print2file("extn.inp", "         rcon=50. nbxmod=-3 repel=3. cutnb=10. ");
	print2file("extn.inp", "      end");
	print2file("extn.inp", "   end");
	print2file("extn.inp", "   flags exclude angle improper end");
	print2file("extn.inp", "   ");
	print2file("extn.inp", "   minimize lbfgs nstep=200 nprint=10 end");
	print2file("extn.inp", "   parameter");
	print2file("extn.inp", "      nbonds");
	print2file("extn.inp", "         rcon=50. nbxmod=-3 repel=0.8 cutnb=6. ");
	print2file("extn.inp", "      end");
	print2file("extn.inp", "   end");
	print2file("extn.inp", "   flags include angle improper end");
	print2file("extn.inp", "   ");
	print2file("extn.inp", "   evaluate (\$count=1) ");
	print2file("extn.inp", "   while (\$count < 5 ) loop l2");
	print2file("extn.inp", "      do (x=x+gauss(0.05)) ( all ) ");
	print2file("extn.inp", "      do (y=y+gauss(0.05)) ( all ) ");
	print2file("extn.inp", "      do (z=z+gauss(0.05)) ( all ) ");
	print2file("extn.inp", "      minimize lbfgs nstep=200 nprint=10 end");
	print2file("extn.inp", "      evaluate (\$count=\$count+1)");
	print2file("extn.inp", "   end loop l2");
	print2file("extn.inp", "   ");
	print2file("extn.inp", "   parameter");
	print2file("extn.inp", "      nbonds");
	print2file("extn.inp", "         rcon=50. nbxmod=3 repel=0.8 cutnb=6. ");
	print2file("extn.inp", "      end");
	print2file("extn.inp", "   end");
	print2file("extn.inp", "   ");
	print2file("extn.inp", "   minimize lbfgs nstep=300 nprint=10 end   ");
	print2file("extn.inp", "   minimize lbfgs nstep=300 nprint=10 end");
	print2file("extn.inp", "   igroup interaction=( not name SG ) ( not name SG ) end");
	print2file("extn.inp", "   energy end");
	print2file("extn.inp", "   evaluate (\$accept=true) ");
	print2file("extn.inp", "   print thres=0.05 bonds");
	print2file("extn.inp", "   evaluate (\$bonds=\$result)");
	print2file("extn.inp", "   if (\$violations > 0) then");
	print2file("extn.inp", "      evaluate (\$accept=false) ");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "   print thres=10. angles ");
	print2file("extn.inp", "   evaluate (\$angles=\$result)");
	print2file("extn.inp", "   if (\$violations > 0) then");
	print2file("extn.inp", "      evaluate (\$accept=false) ");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "   print thres=10. improper");
	print2file("extn.inp", "   evaluate (\$impr=\$result)");
	print2file("extn.inp", "   if (\$violations > 0) then");
	print2file("extn.inp", "      evaluate (\$accept=false) ");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "   print thres=180. dihedral ");
	print2file("extn.inp", "   evaluate (\$dihe=\$result)");
	print2file("extn.inp", "   evaluate (\$trial=\$trial + 1) ");
	print2file("extn.inp", "   if (\$trial > &max_trial ) then");
	print2file("extn.inp", "      exit loop accp");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "end loop accp");
	print2file("extn.inp", "remarks extended strand(s) generation");
	print2file("extn.inp", "remarks input molecular structure file=&structure_file ");
	print2file("extn.inp", "remarks final rms deviations (excluding disulfide bonds): ");
	print2file("extn.inp", "remarks    bonds=	 \$bonds[F8.4] A  ");
	print2file("extn.inp", "remarks    angles=	 \$angles[F8.4] degrees");
	print2file("extn.inp", "remarks    impropers= \$impr[F8.4] degrees");
	print2file("extn.inp", "remarks    dihedrals= \$dihe[F8.4] degrees (not used in some parameter sets!)");
	print2file("extn.inp", "remarks final van der Waals (repel) energy=\$vdw kcal/mole");
	print2file("extn.inp", "write coordinates output=&output_coor format=PDBO end  ");
	print2file("extn.inp", "stop");
}

sub print_usage{
	my $errormsg = shift;
	print "\nERROR!! $errormsg\n\n";
	my $param_info = <<EOF;
--------------------------------------------------------------------------------------
Chromosome3D v1.0 - build chromosome 3D models using HiC Interaction Frequency matrix
--------------------------------------------------------------------------------------

PARAMETERS:
i - Input Interaction Frequency matrix (mandatory)
o - Output directory (mandatory)
k - Scaling parameter (default = 11)
a - Alpha (default = 1/2)
m - Model Count (default = 20)

EXAMPLE:
perl ./chromosome3D.pl -i abc.matrix -o /tmp/abc -k 11 -a 0.5 -m 5

REFERENCE:
"Chromosome3D: a Distance Geometry Method for Reconstructing Three-Dimensional 
Chromosomal Structures from Hi-C Interaction Frequency Data" (submitted),
Badri Adhikari^, Tuan Trieu^, Jianlin Cheng
^These authors contributed equally to this work 
EOF
	print $param_info;
	print "\n";
	exit 1;
}
