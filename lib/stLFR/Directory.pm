package stLFR::Directory;

use strict;
use warnings;
use FindBin;

use base 'Exporter';
our @EXPORT = qw();
our %OPT;

sub run_d{
  %OPT = ref $_[0] eq "HASH" ? %{$_[0]} : @_;
  
  `rm -f $OPT{outputdir}/00.shell/pipeline.*`;
  `mkdir -p $OPT{outputdir}/00.shell/fragment`;
  `mkdir -p $OPT{outputdir}/01.filter` if $OPT{'analysis'} =~ /all|filter/;
  `mkdir -p $OPT{outputdir}/02.align`  if $OPT{'analysis'} =~ /all|align/;
  `mkdir -p $OPT{outputdir}/03.phase`  if $OPT{'analysis'} =~ /all|phase/;
  `mkdir -p $OPT{outputdir}/04.cnv`    if $OPT{'analysis'} =~ /all|cnvsv/;
  `mkdir -p $OPT{outputdir}/04.sv`     if $OPT{'analysis'} =~ /all|cnvsv/;
  `mkdir -p $OPT{outputdir}/report`;
  `mkdir -p $OPT{outputdir}/file`;

};

1;
