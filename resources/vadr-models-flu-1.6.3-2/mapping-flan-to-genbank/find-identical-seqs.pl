#!/usr/bin/env perl
# 
# find-identical-seqs.pl: find identical sequences between two fasta files
#                        
# EPN, Fri Mar 17 11:49:58 2023
# 
use strict;
use warnings;
use Getopt::Long;

my $usage;
$usage  = "find-identical-seqs.pl\n\n";
$usage .= "Usage:\n";
$usage .= "\tfind-identical-seqs.pl [OPTIONS]\n";
$usage .= "\t<fasta file 1 (smaller of the two)\n";
$usage .= "\t<fasta file 2 (bigger of the two)\n";
$usage .= "\n\tOPTIONS:\n";
$usage .= "\t-s1: check if each seq in file 1 is a subseq of a seq in file 2\n";
$usage .= "\t-s2: check if each seq in file 2 is a subseq of a seq in file 1\n";
$usage .= "\n";

my $do_subseq1 = 0;    # set to '1' if -s1 used
my $do_subseq2 = 0;    # set to '1' if -s2 used

&GetOptions( "s1" => \$do_subseq1, 
             "s2" => \$do_subseq2);

if(scalar(@ARGV) != 2) { die $usage; }
my ($in_fa1, $in_fa2) = @ARGV;

if($do_subseq1 && $do_subseq2) { 
  die "ERROR, you can only use one of -s1 or -s2";
}
my $line;

# read in sequences in first input file into a hash
my %seq1_H = ();
my $seqname1 = undef;
my %seqname1_H = ();
open(IN1, $in_fa1) || die "ERROR unable to open $in_fa1 for reading";
while($line = <IN1>) {
  chomp $line;
  if($line =~ /^\>(\S+)/) { 
    $seqname1 = $1;
    if(defined $seqname1_H{$seqname1}) { 
      die "ERROR read $seqname1 twice in $in_fa1";
    }
    $seqname1_H{$seqname1} = 1;
  }
  else { 
    $line =~ s/\s//g;
    $line =~ tr/a-z/A-Z/;
    if(! defined $seqname1) { 
      die "ERROR didn't read sequence name line before sequence line in $in_fa1";
    }
    if(! defined $seq1_H{$seqname1}) { $seq1_H{$seqname1} = ""; }
    $seq1_H{$seqname1} .= $line;
  }
}
close(IN1);
foreach $seqname1 (sort keys %seq1_H) { 
  #printf(">$seqname1\n" . $seq1_H{$seqname1} . "\n");
}

# for each sequence in the second file, check if there's an identical seq from the first file
my %seqname2_H = ();
my $new_seqname2;
my $seqname2 = undef;
my $seq2 = "";
open(IN2, $in_fa2) || die "ERROR unable to open $in_fa2 for reading";
while($line = <IN2>) {
  chomp $line;
  if($line =~ /^\>(\S+)/) { 
    $new_seqname2 = $1;
    if(defined $seqname2) { # not the first sequence
      foreach $seqname1 (sort keys %seq1_H) { 
        if($do_subseq1) { 
          if($seq1_H{$seqname1} eq $seq2) {
            printf("$seqname1 $seqname2 identical\n");
          }
          elsif($seq2 =~ /$seq1_H{$seqname1}/) { 
            my $start = $-[0];
            my $stop  = $+[0];
            printf("$seqname1 is a subseq of $seqname2 %d..%d %d-missing-3'\n", $start+1, $stop, (length($seq2)-$stop));
          }
        }
        elsif($do_subseq2) { 
          if($seq1_H{$seqname1} eq $seq2) {
            printf("$seqname1 $seqname2 identical\n");
          }
          elsif($seq1_H{$seqname1} =~ /$seq2/) { 
            my $start = $-[0];
            my $stop  = $+[0];
            printf("$seqname2 is a subseq of $seqname1 %d..%d %d-missing-3'\n", $start+1, $stop, (length($seq1_H{$seqname1})-$stop));
          }
        }
        else { 
          if($seq1_H{$seqname1} eq $seq2) {
            printf("$seqname1 $seqname2 identical\n");
          }
          else { 
            ;#printf("$seqname1 $seqname2 not identical\n");
          }
        }
      }
    }
    $seqname2 = $new_seqname2;
    if(defined $seqname2_H{$seqname2}) { 
      die "ERROR read $seqname2 twice in $in_fa2";
    }
    $seqname2_H{$seqname2} = 1;
    $seq2 = "";
  }
  else { 
    $line =~ s/\s//g;
    $line =~ tr/a-z/A-Z/;
    if(! defined $seqname2) { 
      die "ERROR didn't read sequence name line before sequence line in $in_fa2";
    }
    $seq2 .= $line;
  }
}
close(IN2);

# and check final sequence
if(defined $seqname2) { # not the first sequence
  foreach $seqname1 (sort keys %seq1_H) { 
    if($do_subseq1) { 
      if($seq1_H{$seqname1} eq $seq2) {
        printf("$seqname1 $seqname2 identical\n");
      }
      elsif($seq2 =~ /$seq1_H{$seqname1}/) { 
        my $start = $-[0];
        my $stop  = $+[0];
        printf("$seqname1 is a subseq of $seqname2 %d..%d %d-missing-3'\n", $start+1, $stop, (length($seq2)-$stop));
      }
    }
    elsif($do_subseq2) { 
      if($seq1_H{$seqname1} eq $seq2) {
        printf("$seqname1 $seqname2 identical\n");
      }
      elsif($seq1_H{$seqname1} =~ /$seq2/) { 
        my $start = $-[0];
        my $stop  = $+[0];
        printf("$seqname2 is a subseq of $seqname1 %d..%d %d-missing-3'\n", $start+1, $stop, (length($seq1_H{$seqname1})-$stop));
      }
    }
    else { 
      if($seq1_H{$seqname1} eq $seq2) {
        printf("$seqname1 $seqname2 identical\n");
      }
      else { 
        ;#printf("$seqname1 $seqname2 not identical\n");
      }
    }
  }
}

