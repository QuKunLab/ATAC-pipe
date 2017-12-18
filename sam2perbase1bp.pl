#!/usr/bin/perl

$half_window =0;
$window_size = 2*$half_window + 1;
$tempseq = "A";
$tempqual = "C";

$chrfile = $ARGV[0];
$infile = $ARGV[1];
$outfile = $ARGV[2];

if ($#ARGV ne 2) {
  print "commond line: perl bed2perbases.pl chrsize_file input_file output_file\n";
  exit;
}

open (in, "<$chrfile");
while ($line=<in>) {
  chomp $line;
  @data = split /\t/, $line;
  $chr = $data[0];
  $chrsize{$chr} = $data[1];
}
close in;

open (in, "<$infile");
open (out, ">$outfile");

while ($line=<in>) {

  $row_count ++;
  if ($row_count%1000000 eq 0) {
    print "$row_count\n";
  }

  chomp $line;
  @data = split /\t/, $line;
  $chr = $data[2];
  $posi = $data[3];
  $length = $data[8];


  if ($length>0) {
    $data[3] = $posi - $half_window;
    $data[7] = $data[3] + $length;
    $data[8] = $length + 2*$half_window;
  }
  else {
    $read_length = $data[5];
    $read_length =~s/M$//;
    $data[3] = $posi+$read_length-$half_window;
    $data[7] = $data[3] + $length;
    $data[8] = $length - 2*$half_window;
  }

  $data[5] = "$window_size"."M";
  $data[9] = $tempseq;
  $data[10] = $tempqual;

  if (($data[1] > 0) and ($data[2] <$chrsize{$chr})) {
    $pline = join "\t", @data;
    print out "$pline\n";
  }
}

close in;
close out;
