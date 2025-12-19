mkdir("SNP");
mkdir("SNP/Homo_sapiens");
@chroms=();
%hashchrom=();
open AA,"SNPs_hg38.bed";
while(<AA>){
s/\r|\n//g;
@temp=split/\t/,$_;
	if($temp[0] ne "" && !exists $hashchrom{$temp[0]}){
	$hashchrom{$temp[0]}=$_."\n";
	push @chroms,$temp[0];
	}elsif($temp[0] ne "" && exists $hashchrom{$temp[0]}){
	$hashchrom{$temp[0]}.=$_."\n";
	}
}
close AA;

foreach $chrom (@chroms){
open AA,">SNP/Homo_sapiens/".$chrom;
print AA $hashchrom{$chrom};
close AA;
}
