#!/usr/bin/perl

$plus_shift = 4;
$minus_shift = -5;


$chrfile = $ARGV[0];
$infile = $ARGV[1];
$outfile = $ARGV[2];

if ($#ARGV ne 2) {
  print "commond line: perl shift_sam_bases.pl chrsize_file input_file output_file\n";
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

  if (($line=~/^\@HD/) or ($line=~/^\@SQ/)) {
    print out $line;
  }
  
  else {
    $row_count ++;
    if ($row_count%1000000 eq 0) {
      print "$row_count\n";
    }

    chomp $line;
    @data = split /\t/, $line;
    $chr = $data[2];
    $plus_posi = $data[3];
    $minus_posi = $data[7];
    $length = $data[8];
    
    $flag = $data[1];
    $flagbin = dec2bin($flag);
    @flagbin = split "", $flagbin;
    
    if (($flagbin[-3] eq 0) and ($flagbin[-5] eq 0)) {
      $data[3] = $plus_posi + $plus_shift;
      $data[7] = $minus_posi + $minus_shift;
      $data[8] = abs($length) + $minus_shift - $plus_shift;
    }
    elsif (($flagbin[-3] eq 0) and ($flagbin[-5] eq 1)) {
      $data[3] = $plus_posi + $minus_shift;
      $data[7] = $minus_posi + $plus_shift;
      $data[8] = -(abs($length) + $minus_shift - $plus_shift);
    }
    
    if (($data[3] > 0) and ($data[3] <$chrsize{$chr})) {
      $pline = join "\t", @data;
      print out "$pline\n";
    }
  }

}

close in;
close out;


$decimalnum = bin2dec($binarynum);
# Convert a binary number to a decimal number
sub bin2dec {
  unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}


$binarynum = dec2bin($decimalnum);
# Convert a decimal to a binary
sub dec2bin {
  my $str = unpack("B32", pack("N", shift));
  $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
  return $str;
}

