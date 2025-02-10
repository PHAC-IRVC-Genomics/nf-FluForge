my $usage = "perl rename-fasta.pl <fasta file name>";

if(scalar(@ARGV) != 1) { 
  die $usage;
}
my ($fa_file) = (@ARGV);

if($fa_file !~ m/\.fasta/) { 
  die "ERROR $fa_file expected to end with .fasta";
}

open(FA, $fa_file) || die "ERROR unable to open $fa_file";

my $fa_root = $fa_file;
$fa_root =~ s/\.fasta//;

while(my $line = <FA>) { 
  chomp $line;
  if($line =~ m/^\>/) { 
    if($line =~ /^\>(\S+)(.*)$/) { 
      my ($name, $desc) = ($1, $2); 
      printf(">%s.%s%s\n", $fa_root, $name, $desc);
    }
    else { 
      die "ERROR unable to parse line $line";
    }
  }
  else { 
    print $line . "\n";
  }
}
