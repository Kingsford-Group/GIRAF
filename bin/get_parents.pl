#!/usr/bin/perl

sub intersection_size {

    local($set1, $set2) = @_;
    local(%set1, $intersection, $size1, $size2) = ((), 0, 0, 0);
    map { $set1{$_} = 1; $size1++; } split(/ /, $set1);
    map { $intersection++ if($set1{$_}); $size2++; } (split(/ /, $set2));

    return ($intersection, $size1, $size2);
}
	
sub supersets {

    local($set, $prefix, *supersets) = @_;
    local($header, $total_trees);

    open(SPLIT, "$prefix\_splits");
    open(TREE, "$prefix\_trees");

    $header = <TREE>; $header =~ m/(\d+)/; $total_trees = $1;
    while(<SPLIT>) {

	$_ =~ m/\{(.*)\} \{(.*)\}/; local($left, $right) = ($1, $2);
	local($id, $count) = split(/\s+/, <TREE>); next if($count/$total_trees < 0.5);
	local($intersection, $size1) = &intersection_size($set, $left);

	if($intersection == $size1) {

	    $supersets{$left} = $count/$total_trees;
	}
	
	if($intersection == 0) {

	    $supersets{$right} = $count/$total_trees;
	}
    }
}

sub minimal_sets {

    local(*sets) = @_;
    local(@ordered_sets) = (sort {scalar(split(/\s+/, $a)) <=> scalar(split(/\s+/, $b))} keys(%sets));
    local($result, $conf) = ($ordered_sets[0], $sets{$ordered_sets[0]});

    for $i (0..$#ordered_sets) {

	next if($sets{$ordered_sets[$i]} eq "");
	($result, $conf) = ($ordered_sets[$i], $sets{$ordered_sets[$i]}) if($sets{$ordered_sets[$i]} > $conf &&
									    scalar(split(/\s+/, $ordered_sets[$i])) == scalar(split(/\s+/, $result)));

	for $j ($i+1..$#ordered_sets) {

	    next if($sets{$ordered_sets[$j]} eq "");
	    local($intersection, $size1) = &intersection_size($ordered_sets[$i], $ordered_sets[$j]);
	    delete $sets{$ordered_sets[$j]} if($intersection == $size1);
	}
    }

    return ($result, $conf);
}

($giraf_file) = @ARGV;

if($giraf_file eq "") {

    print "\nUsage: get_parents.pl <giraf_file>\n\n";
    exit;
}

open(IN, $giraf_file);
while(<IN>) {

    next if($_ =~ m/^\#/);
    ($prefix) = split;
    
    push @prefixes, $prefix;
}

if(@prefixes > 2) {

    print "ERROR: This version of the script can only analyze two segments\n";
    exit;
}

$report = `cat $prefixes[0]-$prefixes[1]_report`; print "<<Lineage Report>>\n\n";

while($report =~ m/Taxa: \{(.*)\}/g) {

    $candidate = $1;

    %supersets1 = %supersets2 = ();
    &supersets($candidate, $prefixes[0], *supersets1);
    ($result1, $conf1) = &minimal_sets(*supersets1);

    &supersets($candidate, $prefixes[1], *supersets2);
    ($result2, $conf2) = &minimal_sets(*supersets2);
    
    print "Reassortment = {$candidate}\n$prefixes[0] parent = {$result1}, Confidence = $conf1\n$prefixes[1] parent = {$result2}, Confidence = $conf2\n\n";
}
