$start = time;
#!note genome version
#@spes=("human","mouse","fly");
#@spes=("gg","mm","at","dr");
#@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
# # # @spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
@spes=("dm6","hg38","mm10");
#@spes=("Drosophila_melanogaster","Homo_sapiens","Mus_musculus");

# @spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
$cutoff=3;
for($i=0;$i<scalar(@spes);$i++){
#foreach $spe (@spes){
	my $dir="AllEPs/".$spes[$i];
	mkdir("download/TargetGene");
	mkdir("download/TargetGene/".$spes[$i]);
	opendir(DIR,$dir) or "can't open the file";
	@dirs=readdir DIR;
	@dirs=grep{$_ ne "." && $_ ne ".."} @dirs;
	foreach $file (@dirs){
	@enhancers=();
	%hashenhancerTogenes=();
	$cell=$file;
	$cell=~ s/\_EP\.txt//g;
	open BB,$dir."/".$file;
		while(<BB>){
		s/\r|\n//g;
		@temp=split/\_|\$/,$_;
		if(!exists $hashenhancerTogenes{$temp[0]}){
		$hashenhancerTogenes{$temp[0]}=$temp[1].":".$temp[2].":";
		push @enhancers,$temp[0];
		}elsif(exists $hashenhancerTogenes{$temp[0]} && !($hashenhancerTogenes{$temp[0]}=~ /\:$temp[2]\:/)){
		$hashenhancerTogenes{$temp[0]}.=";".$temp[1].":".$temp[2].":";
		}
		}
	close BB;
	open AA,">download/TargetGene/".$spes[$i]."/".$cell.".bed";
	foreach $enh (@enhancers){
	print AA $enh."\t".$hashenhancerTogenes{$enh}."\n";
	}
	close AA;
	}
}
$duration = time - $start;
print "All are done: $duration s\n";

#这是将enhancers，geneID和genename单独提取出来