# chrbase chr base strnd coverage freqC freqT
open (IN,"$ARGV[0]");
open (OUT1,">${ARGV[0]}_CG_methykit.txt");
open (OUT2,">${ARGV[0]}_CHG_methykit.txt");
open (OUT3,">${ARGV[0]}_CHH_methykit.txt");
print OUT1 "chrbase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n";
print OUT2 "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n";
print OUT3 "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n";

while (<IN>) {
   chomp;
@a=split;
$chrbase=$a[0].".".$a[1];
$chr=$a[0];
$base=$a[1];
if ($a[2]=~/-/) { $strand= "R"; 
}
else { $strand ="F";}
$coverage=$a[3]+$a[4];
if ($coverage < 3) {next;}
$freqC= $a[3]/($a[3]+$a[4]) * 100;
$freqT= $a[4]/($a[3]+$a[4]) * 100;
$str=$chrbase."\t".$chr."\t".$base."\t".$strand."\t".$coverage."\t".$freqC."\t".$freqT."\n";
if ($a[5]=~ /CG/) { print OUT1 $str;
} 
elsif($a[5]=~/CHG/) {print OUT2 $str;}
elsif($a[5]=~/CHH/){ print OUT3 $str;}
else{print $_."\n";}
}
