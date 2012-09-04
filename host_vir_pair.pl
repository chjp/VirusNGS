#!/usr/bin/perl
=head1 Name

	host_vir_pair.pl
	version:1
	Date, Sep.04.2012
	Contact: kevinchjp@gmail.com

=head1 Description

	Find reads pair that span virus integration sites in host genome

=head1 Usage

	perl host_vir_pair.pl virSoapSe  hostSoapPe  hostSoapSe  pairSpan  unmapPair 
	Arguments explained
	virSoapSe:  file directory for virus SOAP single-end alignment output
	hostSoapPe: file dir for host SOAP paired-end output
	hostSoapSe: file dir for host SOAP single-end output
	pairSpan:   output containing paired reads than align to both host and virus genome
	unmapPair:  output containing potential PE with one mate map to virus genome but the other unmap to host genome. NA in current                        version

=cut

use strict;
use warnings;
use Getopt::Long;

my ($virSoapFile,$hostSoapPe, $hostSoapSe,  $pairSpan, $outUnmap) = @ARGV;
GetOptions(
	"virSoapSe" => \$virSoapFile,
	"hostSoapPe" => \$hostSoapPe,
	"hostSoapSe" => \$hostSoapSe,
	"pairSpan" => \$pairSpan,
	"unmapPair" => \$outUnmap
	);

die `pod2text $0` unless $virSoapFile && $hostSoapPe && $hostSoapSe && $pairSpan && $outUnmap;

my (%id);
my ($fhHBV, $fhHostPe, $fhHostSe, $fhUnmap);
#--------------------- open virus soap file-----------------
if ($virSoapFile =~ /\.gz$/){
	open $fhHBV, "gzip -dc $virSoapFile|" or die $!;
}else{
	open $fhHBV, $virSoapFile or die $!;
}

#---------------------read virus soap file--------------
while (<$fhHBV>){
	chomp;
	next if /^\s+$/;
	my @tmp = split;
	my ($id, $strand) = $tmp[0] =~ /^(\S+)([12])$/;

	if ($id{$id}){
		delete $id{$id};
	}else{
		$id{$id} = $_;
	}
}
close $fhHBV;

#--------------------read hostSoapPe, get read ids and check the %id for span--------
if ($hostSoapPe =~ /\.gz$/){
		open $fhHostPe, "gzip -dc $hostSoapPe|" or die $!;
}else{
		open $fhHostPe, $hostSoapPe or die $!;
}
my $fhout;
open $fhout, ">$pairSpan" or die $!;
while (<$fhHostPe>){
	chomp;
	next if /^\s+$/;
	my @tmp = split;
	my ($id, $strand) = $tmp[0] =~ /^(\S+)([12])$/;
	next unless $id && $strand;
	if($id{$id}){
			my ($hbvReadStr) = $id{$id} =~ /^\S+?([12])\t/; 
			print $fhout "$_\t$id{$id}\n" if ($strand!=$hbvReadStr);
	}
}
close $fhHostPe;
close $fhout;

#--------------------read hostSoapSe, get read ids and check the %id for span---------
if ($hostSoapSe =~ /\.gz$/){
    open $fhHostSe, "gzip -dc $hostSoapSe|" or die $!; 
}else{
     open $fhHostSe, $hostSoapSe or die $!; 
}

open $fhout, ">>$pairSpan" or die $!;
while (<$fhHostSe>){
	chomp;
	next if /^\s+$/;
	my @tmp = split;
	my ($id, $strand) = $tmp[0] =~ /^(\S+)([12])$/;
	next unless $id && $strand;
	print $fhout "$_\t$id{$id}\n"	if ($id{$id} && $strand ne $id{$id}); # if hostMapSe read pair with virusmap, print them. first column hostMap, second column virusmap
	delete $id{$id};
}
close  $fhHostSe;
close $fhout;

#---------------------open host soap unmap file-----------
#my $tmp = `file $hostUnmapFile`;
#if ($tmp =~ /gzip/){
#	open $fhUnmap, "gzip -dc $hostUnmapFile|" or die $!;
#}else{
#	open $fhUnmap, $hostUnmapFile or  die $!;
#}
#open $fhout, ">$outUnmap.virus" or die $!;
#open my $fhout2, ">$outUnmap.host" or die $!;
#
#if (keys %id){
#	while (my $line = <$fhUnmap>){ # fhUnmap reads hostUnmap
#		chomp $line;
#		my $seq = <$fhUnmap>;
#		chomp $seq;
#		my ($id, $strand) = $line =~ /^>(\S+)([12])$/;
#		die unless $id && $strand;
#		
#		next unless $id{$id}; # hash id keys are readID in HBV map file ; if the read from hostUnmap doesn't pair with HBVmap, next
## if the current read pairs with other reads in hostUnmap, do the remaining
#		my @tmp = split /\t/, $id{$id}; # split HBVmap soap line
#		next if $tmp[0] eq "$id$strand";
#		
#		print $fhout ">$tmp[0]\n$tmp[1]\n"; # HBVmap
#		print $fhout2 ">$id$strand\n$seq\n" ;# hostUnMap these two reads pair
#		
#		delete $id{$id};
#	}
#}
#close $fhUnmap;
#close $fhout;
