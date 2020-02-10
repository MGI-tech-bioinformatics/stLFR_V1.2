package stLFR::ParameterCheck;

use strict;
use warnings;
use FindBin;

use base 'Exporter';
our @EXPORT = qw();
our %OPT;

sub run_pc{
  %OPT = ref $_[0] eq "HASH" ? %{$_[0]} : @_;
  
  # check workflow
  if( ($OPT{'analysis'} =~ /align/ && !($OPT{'analysis'} =~ /all|base|filter/             || $OPT{'inputdir'})) ||
      ($OPT{'analysis'} =~ /phase/ && !($OPT{'analysis'} =~ /all|base|filter|align/       || $OPT{'inputdir'})) ||
      ($OPT{'analysis'} =~ /cnvsv/ && !($OPT{'analysis'} =~ /all|base|filter|align|phase/ || $OPT{'inputdir'}))
    ){
    die "ERROR: incomplete workflow!
    There is no INPUTDIR by -inputdir, please add it by -inputdir.
    Otherwise, use -analysis all or -analysis base to complete your workflow.

    For more information, please use -man.
    ";
  }

  # check task module
  if($OPT{'task'} == 2 && !$OPT{'name'}){
    die "ERROR: missing option!
    No project name for MONOTOR.PY since you chose monitor.
    Please add it by -name.

    For more information, please use -man.
    ";
  }  

};

1;
