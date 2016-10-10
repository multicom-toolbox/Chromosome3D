#!/usr/bin/perl -w
# Badri Adhikari, 6/16/2015

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;
use lib '/home/bap54/chr3D_projects/project01/scripts';
use confold;
use lib '/home/bap54/chr3D_projects/project01/scripts/Statistics-RankCorrelation-0.1204/lib';
use Statistics::RankCorrelation;

my $mat = shift;
my $pdb = shift;
my $range = shift; # threshold to define short range

confess "ERROR! Params: IF matrix, pdb\n" if not ($mat and $pdb);
confess "ERROR! IF matrix not found $mat\n" if not -f $mat;

$range = 3 if not defined $range;

my @pdb_list = ();
if (-f $pdb){
	push @pdb_list, $pdb;
}
else{
	@pdb_list = load_pdb($pdb);
}

my %srcc = ();
foreach my $p (@pdb_list){
	my @dist_mat = ();
	my @dist_pdb = ();
	my %ca_xyz1 = ca_xyz_pdb($p);
	my %ca_xyz2 = %ca_xyz1;
	if ($range >= (scalar keys %ca_xyz1)){
		print "Spearman Correlation coefficient = -";
		exit;
	}
	foreach my $r1 (sort {$a <=> $b} keys %ca_xyz1){
		foreach my $r2 (sort {$a <=> $b} keys %ca_xyz2){
			next if abs($r1 - $r2) < $range;
			my $d = calc_dist($ca_xyz1{$r1}, $ca_xyz2{$r2});
			#confess "\n\n".$pdb."\n\n";
			push @dist_pdb, (sprintf "%.3f", $d);
		}
	}
	my $h = 0;
	open MATRIX, $mat or confess $!;
	while(<MATRIX>){
		chomp $_;
		$_ =~ s/^\s+//;
		confess "??" if length($_) < 2;
		my @R = split /\s+/, $_;
		for(my $i = 0; $i <= $#R; $i++){
			next if abs($i - $h) < $range;
			push @dist_mat, $R[$i];
		}
		$h++
	}
	close MATRIX;
	confess "ERROR! mismatch in size!" if $#dist_mat != $#dist_pdb;
	my $c = Statistics::RankCorrelation->new( \@dist_mat, \@dist_pdb, sorted => 1 );
	my $n = $c -> spearman;
	$srcc{$p} = sprintf "%.3f", $n;
	my $c1 = Statistics::RankCorrelation->new( \@dist_pdb, \@dist_mat, sorted => 1 );
	my $n2 = $c1 -> spearman;
	confess "ERROR! correlation mismatch $n and $n2 !" if $n != $n2;
}

print "SRCC\tPDB\n"; 
foreach (sort {$srcc{$b} <=> $srcc{$a}} keys %srcc){
	printf $srcc{$_}."\t".$_."\n";
}
