#@spes=("hs","mm","dm");
#@speups=("HS","MM","DM");

# # # #@spes=("gg","mm","at","dr");
# # # @spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
# # # #@speups=("GG","MM","AT","DR");
# # # @speups=("CJ","MM","MF","CS","HS","GG","DR","PT","RM","DM","AT","OSJ","ZM");

@spes=("dm6","hg38","mm10");
@speups=("DM6","HG38","MM10");
# @speups=("TAIR10","CE11","CHLSAB1","DANRER11","DM6","ASM31002V1","GALGAL6","HG38","MM10","ASM251V1","RHEMAC10","IRGSP1",
#  "PANTRO6","ASM276V2","RN7","SACCER3","ASM294V2","SL3","SPUR5","SUSSCR11","TGA4","TRYBRU","XTRO10","ZM-B73");

mkdir("enhdetail");
for($i=0;$i<scalar(@spes);$i++){
mkdir("enhdetail/".$spes[$i]);
my $dir="enhs/".$spes[$i];    #dm6	hg38	mm10   HG3801   #pay attention to the counts of enhs/$spes[$i]
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
#@dir=grep{$_ ne "." && $_ ne ".."} @dir;
	foreach $file (@dir){
	if($file ne "." && $file ne ".."){
	%hashnum=();
	@arraynum=();
		open AA,$dir."/".$file;
#		print $file;
		mkdir("enhdetail/".$spes[$i]."/".$file);
		while(<AA>){
		chomp($_);
		@temp=split/\t/,$_;
		$temp[0]=~ s/$speups[$i]\d\d-//g;
		$start=int(($temp[0]-1)/100);
		$end=$start+1;
			if(!exists $hashnum{$start}){
			push @arraynum,$start;
			$hashnum{$start}=$start;
			open $start,">enhdetail/".$spes[$i]."/".$file."/".($start*100+1)."-".($start*100+100);
			print $start $_."\n";
			}else{
			print $start $_."\n";
			}
		}
		close AA;
		foreach $num (@arraynum){
		close $num;
		}
	}
}

}

