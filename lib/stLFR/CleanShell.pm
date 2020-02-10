package stLFR::CleanShell;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../";
use lib "$FindBin::RealBin/../lib/perl5";
use lib "$FindBin::RealBin/../perl5";
use lib "$FindBin::RealBin/perl5";
use MyModule::GlobalVar qw($REF_H_DBPATH $TOOL_PATH $QUEUE_PM);
use stLFR::ResultLink;

use base 'Exporter';
our @EXPORT = qw();
my (
  $monitor, $python2,
);

$monitor = "$TOOL_PATH/monitor/monitor";

sub report_cs{
  my (
    $fqprefix,
    $alignprefix,
    $phaseprefix,
    $cnvsvprefix,
    $shelldir,
    $task,
    $project,
    $queue,
    $name,
  ) = (@_); 

  if(1){
    `cat $shelldir/t0141.filter_fastqdistribution.sh >> $shelldir/t0321.phase_phase.sh `;
    `cat $shelldir/t0142.filter_fastqdistribution.sh >> $shelldir/t0322.phase_phasestat.sh `;
    `cat $shelldir/t0231.align_alignstat.sh          >> $shelldir/t0321.phase_phase.sh `;
    `cat $shelldir/t0232.align_alignstat.sh          >> $shelldir/t0322.phase_phasestat.sh `;
  }

  if($task == 1){
    `echo -e "sh $shelldir/pipeline.sh 1>$shelldir/pipeline.sh.o 2>$shelldir/pipeline.sh.e &" > $shelldir/run.sh`;
    `rm $shelldir/pipeline.txt`;
  }else{
    `cat $shelldir/pipeline.txt | awk 'NF==2' > $shelldir/pipeline.txt.2`;
    `mv $shelldir/pipeline.txt.2 $shelldir/pipeline.txt`;
    `echo -e "
      source /opt/gridengine/default/common/settings.sh
      export LD_LIBRARY_PATH=/hwfssz1/ST_GCCNT/LIBRARYDEV/PROJECT/F16ZQSB1SY2523/lizhanqing/software/Python-2.7.10/lib:\\\${LD_LIBRARY_PATH}
      export PYTHONPATH=/ldfssz1/MGI_BIT/RUO/lizhanqing/Software/PyMonitor/V.1.5/lib:/hwfssz1/ST_GCCNT/LIBRARYDEV/PROJECT/F16ZQSB1SY2523/lizhanqing/software/Python-2.7.10/bin/python
      export DRMAA_LIBRARY_PATH=/opt/gridengine/lib/linux-x64/libdrmaa.so
      export PYMONITOR_PY_PATH=/ldfssz1/MGI_BIT/RUO/lizhanqing/Software/PyMonitor/V.1.6/pymonitor.py
      export PYMONITOR_LOG_PATH=~/.pymonitor.log
      export PYMONITOR_CONF_PATH=~/.pymonitor.conf
      $monitor taskmonitor -f 1 --q1 $queue --q2 $alignprefix --P1 $project --P2 $project -p $name -i $shelldir/pipeline.txt
    " > $shelldir/run.sh `;
    `rm $shelldir/pipeline.sh $shelldir/t0*sh`;
  }
  
};

1;
