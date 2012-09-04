#!/usr/bin/perl
#contact: kevinchjp@gmail.com

use strict;
use warnings;
use Data::Dumper;


my ($pairSoapFile) = @ARGV;
die "perl readsCluster.pl pairSpan \n" unless $pairSoapFile;

my (%fEnd, %rEnd, %cluster);
parsePair($pairSoapFile, \%fEnd, \%rEnd);# write pe align info to %fEnd(host) %rEnd(virus) %exist;  fEnd is chr7, 84677378, +  rEnd is NC_009758, 8836, +    exist seem to no use
cluster(\%fEnd, \%rEnd, \%cluster);
display(\%fEnd, \%rEnd, \%cluster);

sub parsePair{ #function:we have pairSe, then exist{hostChr_start}{virChr_start}
	my ($file, $fEnd, $rEnd) = @_;
	my $fh;
	if ($file =~ /\.gz/){
		open $fh, "gzip -dc $file |" or die $!;
	}else{
		open $fh, $file or die $!;
	}
	
	my %exist;# store pair-info $exist{host_alig}{vir_alig}=1 or >1
	my $index = 0;
	my $i;
	while (<$fh>){#read pairSe : host align and HBV align
		my @tmp = split /\t/; #   pair with hostMapSe virus (hbv) map se;
		next if $tmp[3] > 1; # $tmp[3] is column 4, number of equal best hits;  unique hit
		for ($i = 9; $i < @tmp; $i++){
			last if $tmp[$i] eq "+" || $tmp[$i] eq "-"; 
		}

		my $next = $i - 6;#$next = vir readsID
		next if $tmp[$next + 3] > 1; # on HBV, hit > 1
		
		next if $exist{"$tmp[7]_$tmp[8]"}{"$tmp[$next + 7]_$tmp[$next + 8]"}; #exist{hostChr_start}{virChr_start}
		$fEnd->{$index} = [$tmp[7], $tmp[8], $tmp[6]]; #7--chr  8--start 6--strand
		$rEnd->{$index} = [$tmp[$next + 7], $tmp[$next + 8], $tmp[$next + 6]]; #notice $rEnd->{$index} is a ref to an array, so there actually 2 reference: 1 is $rEnd || 2 is the hash value is a ref to an array, see the next line test
#print "@{$rEnd->{$index}}\n";die;# the print is : NC_009758 8836 +
		$exist{"$tmp[7]_$tmp[8]"}{"$tmp[$next + 7]_$tmp[$next + 8]"}++;#exist{hostChr_start}{virChr_start}
		$index++;
	}
	close $fh;
}

sub cluster{ #info 1: 
	my ($fEnd, $rEnd, $clusterRef) = @_;
#fEnd--host rEnd--virus;  $fEnd->{index} index ~ line# in pairSe, $fEnd->{index} array ref--[chr,star,strand]	
	my @cluster;
	my ($fstart, $fend, $rstart, $rend) = (0) x 4;
	my ($fchr, $rchr) = ("") x 2;
	my $len = 90; #read length. If your sequencing lib is 75, changed! warning!
	
	my @index = (sort {$fEnd->{$a}[0] cmp $fEnd->{$b}[0] || $fEnd->{$a}[1] <=> $fEnd->{$b}[1]} keys %$fEnd);# pe align info   fEnd key is (1,2,3..) value is array: chr7, 84677378, +  rEnd is NC_009758, 8836, +
	for (my $i = 0; $i < @index; $i++){
#----------------------------------------
		unless ((defined $fchr) && ($fchr eq  $fEnd->{$index[$i]}[0] ) 
			&& 	 (defined $rchr)  && ($rchr eq $rEnd->{$index[$i]}[0])){ 	
			if (@cluster){
				push @{$clusterRef->{$fchr}{$rchr}}, [$fstart, $fend, $rstart, $rend, [@cluster]]; #push last pair info to @clusterRef
				@cluster = ();
			}

			$fchr = $fEnd->{$index[$i]}[0];   $rchr = $rEnd->{$index[$i]}[0]; # fchr rchr: f means host; r means virus
			$fstart = $fEnd->{$index[$i]}[1]; $rstart = $rEnd->{$index[$i]}[1];#fstart fend is the read align locus
			$fend = $fstart + $len - 1;	$rend = $rstart + $len - 1;
			push @cluster, $index[$i]; #@index[] is %$fEnd{}keys 1 2 3. sorted by chr and start
			next;
		}		
#----------------------------------------
		if (( ($fEnd->{$index[$i]}[1] < $fstart - $len + 1) || ($fEnd->{$index[$i]}[1] > $fend + $len - 1) ) # if current read doesn't overlap with the last read
		|| ( ($rEnd->{$index[$i]}[1] < $rstart - $len + 1) || ($rEnd->{$index[$i]}[1] > $rend + $len - 1) ) ){
			if (@cluster){
				push @{$clusterRef->{$fchr}{$rchr}}, [$fstart, $fend, $rstart, $rend, [@cluster]];
				@cluster = ();
			}
			
			$fchr = $fEnd->{$index[$i]}[0];   $rchr = $rEnd->{$index[$i]}[0];
			$fstart = $fEnd->{$index[$i]}[1]; $rstart = $rEnd->{$index[$i]}[1];
			$fend = $fstart + $len - 1;	$rend = $rstart + $len - 1;
			push @cluster, $index[$i];
			next;
		}

#-----------if current read overlaps with last read, then update current read.
		push @cluster, $index[$i];
		$fstart = $fEnd->{$index[$i]}[1] if $fEnd->{$index[$i]}[1] < $fstart; 
		$rstart = $rEnd->{$index[$i]}[1] if $rEnd->{$index[$i]}[1] < $rstart;
		$fend = $fEnd->{$index[$i]}[1] + $len -1 if $fEnd->{$index[$i]}[1] + $len -1 > $fend;
		$rend = $rEnd->{$index[$i]}[1] + $len - 1 if $rEnd->{$index[$i]}[1] + $len - 1 > $rend; 
	}
	
	if (@cluster){
		push @{$clusterRef->{$fchr}{$rchr}}, [$fstart, $fend, $rstart, $rend, [@cluster]];
		@cluster = ();
	}
}

sub display{
	my ($fEnd, $rEnd, $cluster) = @_;
	
	for my $fchr(keys %$cluster){
		for my $rchr (keys %{$cluster->{$fchr}}){
			for my $p ( @{$cluster->{$fchr}{$rchr}}){
				my ($fstart, $fend, $rstart, $rend, $indexRef) = @$p;
				my $size = @$indexRef;
				my ($fcluster, $rcluster, $fstrand, $rstrand);
				for (my $i = 0; $i < @$indexRef; $i++){
					$fcluster .= "$fEnd->{$$indexRef[$i]}[1],";
					$fstrand .=  "$fEnd->{$$indexRef[$i]}[2],";
					$rcluster .= "$rEnd->{$$indexRef[$i]}[1],";
					$rstrand .= "$rEnd->{$$indexRef[$i]}[2],";
				}
				print "$fchr\t$fstart\t$fend\t$size\t$fstrand\t$rchr\t$rstart\t$rend\t$rstrand\n";
			}
		}
	}
}
