#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin);
use Cwd qw(abs_path);
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use lib "$Bin/../../lib/perl5";
use MyModule::GlobalVar qw($REF_H_DBPATH $TOOL_PATH);

my ($list, $result, $shell, $cpu, %flag, $watchDog, $reportdir,
  $line, $convert, $python3, $version, $size,
);
GetOptions(
  "i=s"   => \$list,
  "s=s"   => \$shell,
  "r=s"   => \$reportdir,
);
die "perl $0 -i input -s shell -r reportdir\n" unless defined $list && defined $shell && defined $reportdir;

$watchDog     = "$Bin/../watchDog.pl";
$convert      = "$TOOL_PATH/fqcheck/convert";
$python3      = "$TOOL_PATH/python3/python3";
$cpu        ||= 70;
$line         = 0;
$version      = "stLFR-reSeq V2.0.0.0";
$size       ||= "750x750";

#=============================================#
# build shell
#=============================================#
open S1,">$shell/run9.html.1.sh";

open LIST,$list;
while(<LIST>){
  chomp;
  next if /^#/;
  next if /^sample.*path/;
  my @info = split /\s+/;

  $flag{$info[0]}++;
  next if $flag{$info[0]} > 1;  # pass if duplicated samples

  $line = 0;

  # convert png
  foreach my $pdf(`ls $reportdir/$info[0]/$info[0].*.pdf`){
    chomp $pdf;
    my $png = $pdf;
    $png =~ s/.pdf/.png/;
    print S1 "$convert -resize $size $pdf $png\n";
    $line += 1;
  }
  print S1 "$convert -resize $size $Bin/legend_circos.pdf $reportdir/$info[0]/$info[0].legend_circos.png \n";
  print S1 "cp $Bin/legend_circos.pdf $reportdir/$info[0]/$info[0].legend_circos.pdf \n";
  $line += 2;

  # fix phasing table
  print S1 "perl $Bin/haplotype_fix.pl $reportdir/$info[0]/$info[0].haplotype.xls \n";
  $line += 1;

  # fix statistics table
  print S1 "perl $Bin/statistics_fix.pl $reportdir/Alignment.statistics.xls $info[0] \n";
  print S1 "perl $Bin/variant_fix.pl $reportdir/Variant.statistics.xls $info[0] \n";
  $line += 2;

  # build html
  print S1 "cd $reportdir && $python3 $Bin/stLFR/generate_DNApipe_report.py \'$version\' $info[0] ./ ./ \n";
  $line += 1;

}
close LIST;

close S1;
#=============================================#
# write main shell script
#=============================================#
open MAINSHELL,">>$shell/pipeline.sh";
print MAINSHELL "echo ========== 9.html start at : `date` ==========\n";
print MAINSHELL "perl $watchDog --mem 1g --num_paral $cpu --num_line $line $shell/run9.html.1.sh\n";
print MAINSHELL "echo ========== 9.html   end at : `date` ==========\n\n";
close MAINSHELL;

