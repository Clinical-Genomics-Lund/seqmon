#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use JSON;
use File::Slurp;

my @queue = `/data/bnf/scripts/cmd-job.pl --queue`;

my @ignore = read_file(".juneignore", chomp => 1); 


my %data;
my $which = 0;
foreach( @queue ) {
    chomp;
    next if /^\s*$/;

    if( /QUEUED\:/ ) {
	$which = "q";
    }
    elsif( /RUNNING\:/ ) {
	$which = "r";
    }
    elsif( $which ) {
	my @a = split /\t/;
	my $assay = which_assay( $a[1] );

	next if ($a[0] ~~ @ignore);
	
	if( $which eq "r" ) {
	    $a[7] =~ s/:\d\d//;
	    push( @{$data{running}}, { 'name'=>$a[0], 'start_time'=>$a[7], 'assay'=>$assay } );
	}
	else {
	    $a[2] =~ s/:\d\d//;	    
	    push( @{$data{queue}}, { 'name'=>$a[0], 'start_time'=>$a[2], 'assay'=>$assay } );	    
	}
    }
}

open( JSON, ">/home/bjorn/queue.json");
print JSON "var queue=";
print to_json(\%data);
close JSON;

system("scp /home/bjorn/queue.json pi\@10.0.224.47:/home/pi/seqmon");



sub which_assay {
    if( -e $_[0] ) {
	open( SH, $_[0] );
	while( <SH> ) {
	    if( /\/(myeloid|exome|NIPT|exome_tumor|rnaseq)\// ) {
		return "myeloisk panel" if $1 eq "myeloid";
		return "exom" if $1 =~ "exome";
		return "SWEA" if $1 =~ "swea";
		return "RNA-Seq" if $1 =~ "rnaseq";
		return "Somatiskt exom" if $1 =~ "exome_tumor";
	    
		return $1;
	    }
	}
    }
    return "okänd";
    
}

