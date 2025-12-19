$start = time;
mkdir("download/SNP");
@spes=("hg38");
@species=("Homo_sapiens"); #$

for($i=0;$i<scalar(@spes);$i++){
	#foreach $spe (@spes){
	mkdir("download/SNP/".$spes[$i]);
	my $dir="download/enhancer/".$species[$i];
	opendir(DIR,$dir) or "can't open the file";
	@dir=readdir DIR;
	@dir=grep{$_ ne "." && $_ ne ".."} @dir;
		foreach $cellfile (@dir){
		$cellname=$cellfile;
		$cellname=~ s/\.bed$//g;
		system("bedtools intersect -a ".$dir."/".$cellfile." -b SNPs_hg38.bed -wa -wb>download/SNP/".$spes[$i]."/".$cellfile);
		}
	}

$duration = time - $start;
print "All are done: $duration s\n";