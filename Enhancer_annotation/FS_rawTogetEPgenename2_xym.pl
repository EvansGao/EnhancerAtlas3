use List::MoreUtils qw(uniq);
#@spes=("human","mouse","fly");
#@spes=("mm","dr","at","gg");
#@spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
@spes=("dm6","hg38","mm10");

mkdir("EPgenename2");
for($i=0;$i<scalar(@spes);$i++){
#print $spes[$i]."\n";
mkdir("EPgenename2/".$spes[$i]);
my $dir="AllEPs/".$spes[$i];
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
foreach $filedir (@dir){
if($filedir ne "." && $filedir ne ".."){
	$cell=$filedir;
	$cell=~ s/\_EP\.txt//g;
	open CELL,$dir."/".$filedir;
	@enhs=();
	%hashenh=();
	while(<CELL>){
	s/\r|\n//g;
	@tmp=split/\t/,$_;
	@temp=split/\:|\-|\_|\$/,$tmp[0];
		if(!exists $hashenh{$temp[0]."\t".$temp[1]."\t".$temp[2]}){
		$hashenh{$temp[0]."\t".$temp[1]."\t".$temp[2]}=$temp[4];
		push @enhs,$temp[0]."\t".$temp[1]."\t".$temp[2];
		}else{
		$hashenh{$temp[0]."\t".$temp[1]."\t".$temp[2]}.=",".$temp[4];
		}
	
	}
	close CELL;
	open AA,">".$cell.".bed";
	foreach $enh (@enhs){
	@genenames=split/\,/,$hashenh{$enh};
	@genenames=uniq(@genenames);
	print AA $enh."\t".join(";",@genenames)."\n";
	}
	close AA;
	system("bedtools sort -i ".$cell.".bed>EPgenename2/".$spes[$i]."/".$cell.".bed");
	unlink($cell.".bed");
}
}


}
