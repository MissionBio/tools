#!/usr/bin/perl
use warnings;
use Getopt::Long;

my $input_1   = "read1.fastq.gz";
my $input_2   = "read2.fastq.gz";
my $output_1  = "read1.fastq";
my $output_2  = "read2.fastq";
my $discards = "discarded.seq";
my $verbose = undef;
my $help    = 0;
my $man     = 0;

GetOptions (
    "input-1=s"    => \$input_1,
    "input-2=s"    => \$input_2,
    "output-1=s"   => \$output_1,
    "output-2=s"   => \$output_2,
    "discards=s"   => \$discards,
    "verbose"      => \$verbose,
    'help|?'       => \$help, 
    'man'          => \$man
    );

my $fixed_part_1 = 'AGTACGTACGAGTC';
my $fixed_part_2 = 'GTACTCGCAGTAGTC';
my $len_bc1 = 9;
my $len_fixed1 = length($fixed_part_1);
my $len_bc2 = 9;
my $len_fixed2 = length($fixed_part_2);
my $min_read1_length = 18;
my $total_barcode_length = $len_fixed1 + $len_bc1 + $len_fixed2 + $len_bc2;
my $c1plusc2plusbc2_length = $len_fixed1 + $len_fixed2 + $len_bc2;
my $temp_file = "tag_RG.tmp";

open I1, "zcat $input_1 | " or die "$!\n";
open I2, "zcat $input_2 | " or die "$!\n";
open (O1, ">$output_1") or die "$!\n";
open (O2, ">$output_2") or die "$!\n";
open (OFF1, ">$discards") or die "$!\n";

my $n=0;
my $comments=0;
my %RG=();
$n = 0;
%item_1 = ();
%item_2 = ();
$offset0count=0;
$offset1count=0;
$offset2count=0;
$offset3count=0;
while ($line_1=<I1>) {
	$line_2=<I2>;
	chomp($line_1);
	chomp($line_2);
	$n++;
#	print $n,"\t";
	$mod = ($n % 4);
#	print $mod, "\n";
	$item_1{$mod} = $line_1;
	$item_2{$mod} = $line_2;
	if($mod == 0){
		$read1 = "";
		if(length($item_1{2}) >= $total_barcode_length + $min_read1_length){
			$offset = substr($item_1{2}, $len_bc1, 3);
			#print STDERR $offset."\n";
			if($offset eq 'AGT') {
			    $cell_barcode1='';
			    $cell_barcode2='';
                $cell_barcode = substr($item_1{2}, 0, $total_barcode_length+4);
                $read1 = substr($item_1{2}, $total_barcode_length+4);
                ($read_name_1) = split(/\s/, $item_1{1});
                ($read_name_2) = split(/\s/, $item_2{1});
                print O1 $read_name_1."_RG:Z:$cell_barcode\n";
                print O1 $read1, "\n";
                print O1 $item_1{3}, "\n";
                print O1 substr($item_1{0}, -length($read1)), "\n";

                print O2 $read_name_2."_RG:Z:$cell_barcode\n";
                print O2 $item_2{2}, "\n";
                print O2 $item_2{3}, "\n";
                print O2 $item_2{0}, "\n";
                $offset0count++;
			}
			elsif($offset eq 'CAG') {
			    $cell_barcode1='';
			    $cell_barcode2='';
                $cell_barcode1 = substr($item_1{2}, 0, $len_bc1);
                $cell_barcode2 = substr($item_1{2}, $len_bc1+1, $c1plusc2plusbc2_length+4);
                $cell_barcode=$cell_barcode1.$cell_barcode2;
                $read1 = substr($item_1{2}, $total_barcode_length+5);
                ($read_name_1) = split(/\s/, $item_1{1});
                ($read_name_2) = split(/\s/, $item_2{1});
                print O1 $read_name_1."_RG:Z:$cell_barcode\n";
                print O1 $read1, "\n";
                print O1 $item_1{3}, "\n";
                print O1 substr($item_1{0}, -length($read1)), "\n";

                print O2 $read_name_2."_RG:Z:$cell_barcode\n";
                print O2 $item_2{2}, "\n";
                print O2 $item_2{3}, "\n";
                print O2 $item_2{0}, "\n";
                $offset1count++;
			}
			elsif($offset eq 'TCA') {
			    $cell_barcode1='';
			    $cell_barcode2='';
                $cell_barcode1 = substr($item_1{2}, 0, $len_bc1);
                $cell_barcode2 = substr($item_1{2}, $len_bc1+2, $c1plusc2plusbc2_length+4);
                $cell_barcode=$cell_barcode1.$cell_barcode2;
                $read1 = substr($item_1{2}, $total_barcode_length+6);
                ($read_name_1) = split(/\s/, $item_1{1});
                ($read_name_2) = split(/\s/, $item_2{1});
                print O1 $read_name_1."_RG:Z:$cell_barcode\n";
                print O1 $read1, "\n";
                print O1 $item_1{3}, "\n";
                print O1 substr($item_1{0}, -length($read1)), "\n";

                print O2 $read_name_2."_RG:Z:$cell_barcode\n";
                print O2 $item_2{2}, "\n";
                print O2 $item_2{3}, "\n";
                print O2 $item_2{0}, "\n";
                $offset2count++;
			}
			elsif($offset eq 'GTC') {
			    $cell_barcode1='';
			    $cell_barcode2='';
                $cell_barcode1 = substr($item_1{2}, 0, $len_bc1);
                $cell_barcode2 = substr($item_1{2}, $len_bc1+3, $c1plusc2plusbc2_length+4);
                $cell_barcode=$cell_barcode1.$cell_barcode2;
                $read1 = substr($item_1{2}, $total_barcode_length+7);
                ($read_name_1) = split(/\s/, $item_1{1});
                ($read_name_2) = split(/\s/, $item_2{1});
                print O1 $read_name_1."_RG:Z:$cell_barcode\n";
                print O1 $read1, "\n";
                print O1 $item_1{3}, "\n";
                print O1 substr($item_1{0}, -length($read1)), "\n";

                print O2 $read_name_2."_RG:Z:$cell_barcode\n";
                print O2 $item_2{2}, "\n";
                print O2 $item_2{3}, "\n";
                print O2 $item_2{0}, "\n";
                $offset3count++;
			}
			else{
			print OFF1 $item_1{2}, "\n", $item_2{2}, "\n";
			}

		}else{
			print OFF1 $item_1{2}, "\n", $item_2{2}, "\n";
		}
	}
}
print STDERR $offset0count."\n";
print STDERR $offset1count."\n";
print STDERR $offset2count."\n";
print STDERR $offset3count."\n";
close I1;
close I2;
close O1;
close O2;
close OFF1;
