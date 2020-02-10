package stLFR::ResultLink;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use lib "$FindBin::RealBin/../perl5";
use lib "$FindBin::RealBin/perl5";
use MyModule::GlobalVar qw($REF_H_DBPATH $TOOL_PATH $QUEUE_PM);

use base 'Exporter';
our @EXPORT = qw();

sub reportlink{
  my (
    $file,
    $predir,
    $postdir,
    $filehandle,
  ) = (@_);

  print $filehandle "
    cp $predir/$file $postdir
  ";
};

sub filelink{
  my (
    $file,
    $predir,
    $postdir,
    $filehandle,
  ) = (@_);

  print $filehandle "
    mv $predir/$file $postdir/
    ln -s $postdir/$file $predir/
  ";
  
};

1;