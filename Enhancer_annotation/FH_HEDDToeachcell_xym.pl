$start = time;
mkdir("download/HEDD");
@spes=("hg38");
@species=("Homo_sapiens"); #$

for($i=0;$i<scalar(@spes);$i++){
	#foreach $spe (@spes){
	mkdir("download/HEDD/".$spes[$i]);
	my $dir="download/enhancer/".$species[$i];
	opendir(DIR,$dir) or "can't open the file";
	@dir=readdir DIR;
	@dir=grep{$_ ne "." && $_ ne ".."} @dir;
		foreach $cellfile (@dir){
		$cellname=$cellfile;
		$cellname=~ s/\.bed$//g;
		system("bedtools intersect -a ".$dir."/".$cellfile." -b HEDD_".$spes[$i].".bed -F 1.0 -wa -wb>download/HEDD/".$spes[$i]."/".$cellfile);
		}
	}

$duration = time - $start;
print "All are done: $duration s\n";

#HEDDraw数据从哪里来的？
#DataDownload-Alzheimer_EnhancerScore.txt
#DataDownload-Huntington_EnhancerScore.txt
#DataDownload-Obesity_EnhancerScore.txt
#DataDownload-Parkinson_EnhancerScore.txt
#DataDownload-Prostate_cancer_EnhancerScore.txt
#DataDownload-Schizophrenia_EnhancerScore.txt
#DataDownload-Sleep_disorder_EnhancerScore.txt
