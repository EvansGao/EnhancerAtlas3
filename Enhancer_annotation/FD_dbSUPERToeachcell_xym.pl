$start = time;
mkdir("download/dbSUPER");
@spes=("hg38","mm10");
@species=("Homo_sapiens","Mus_musculus"); #$

for($i=0;$i<scalar(@spes);$i++){
	#foreach $spe (@spes){
		mkdir("download/dbSUPER/".$spes[$i]);
		my $dir="download/enhancer/".$species[$i];
		opendir(DIR,$dir) or "can't open the file";
		@dir=readdir DIR;
		@dir=grep{$_ ne "." && $_ ne ".."} @dir;
		foreach $cellfile (@dir){
			$cellname=$cellfile;
			$cellname=~ s/\.bed$//g;
			system("bedtools intersect -a ".$dir."/".$cellfile." -b dbSUPER_".$spes[$i].".bed -f 0.1 -wa -wb>download/dbSUPER/".$spes[$i]."/".$cellfile);
		}
	}

$duration = time - $start;
print "All are done: $duration s\n";

#this pl is used to get the overlap between tipical enhancers and super enhancers.