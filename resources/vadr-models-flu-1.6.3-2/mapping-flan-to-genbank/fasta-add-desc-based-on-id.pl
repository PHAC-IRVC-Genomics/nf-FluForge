my $usage = "perl fasta-add-desc-based-on-id.pl <id file> <fasta file>";

if(scalar(@ARGV) != 2) { 
  die $usage;
}
my ($id_file, $fa_file) = (@ARGV);

my %flan_name2accn_H = ();
open(ID, $id_file) || die "ERROR unable to open $id_file";
my $line;
while($line = <ID>) { 
  chomp $line;
  #A-seg6_N8 gi|78069952|gb|CY004056.1| identical
  if($line =~ /^(\S+)\s++gi\|\d+\|\S+\|(\S+\.\d+)\|.*\s+identical$/) { 
    my ($flan_name, $accn) = ($1, $2);
    $flan_name2accn_H{$flan_name} = $accn;
  }
}
close(ID);

open(FA, $fa_file) || die "ERROR unable to open $fa_file for reading";

while($line = <FA>) { 
  #>A-seg1
  chomp $line;
  if($line =~ m/^\>(\S+)/) { 
    my $flan_name = $1;
    if(! defined $flan_name2accn_H{$flan_name}) { 
      die "ERROR flan name $flan_name did not exist in the id file $id_file";
    }
    printf(">$flan_name $flan_name2accn_H{$flan_name}\n");
  }
  else { 
    print $line . "\n";
  }
}
close(FA);
