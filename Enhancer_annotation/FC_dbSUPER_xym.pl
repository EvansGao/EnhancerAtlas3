system("wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz");
system("wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz");
@spes=("hs","mm");
@spesliftover=("hg19ToHg38.over.chain.gz","mm9ToMm10.over.chain.gz");
@spever=("hg38","mm10");
@rawver=("hg19","mm9");

%hashrawverTOspe=();
%hashspeverTOspe=();
%hashspesliftoverTOspe=();
for($i=0;$i<scalar(@spes);$i++){
	$hashspesliftoverTOspe{$spes[$i]}=$spesliftover[$i];
	$hashspeverTOspe{$spes[$i]}=$spever[$i];
	$hashrawverTOspe{$spes[$i]}=$rawver[$i];
}

foreach $spe (@spes){
	my $dir="dbSUPERraw/".$spe;
	opendir(DIR,$dir) or "can't open the file";
	@dirs=readdir DIR;
	@dirs=grep{$_ ne "." && $_ ne ".."} @dirs;
	open AA,">dbSUPERpre_".$spe.".bed";
	foreach $file (@dirs){
	$cell=$file;
	$cell=~ s/^\d+\_|\.bed$//g;
	open BB,$dir."/".$file;
		while(<BB>){
		s/\r|\n//g;
		@temp=split/\t/,$_;
		if($temp[4]=~ /^SE\_\d+|^mSE\_\d+/){
		print AA $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[4].":".$cell."\n";
		}elsif($temp[3]=~ /^SE\_\d+|^mSE\_\d+/){
		print AA $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3].":".$cell."\n";
		}

		
		}
	close BB;
	}
	close AA;
	system("bedtools sort -i dbSUPERpre_".$spe.".bed>dbSUPER_".$spe.".".$hashrawverTOspe{$spe}.".bed");
	system("liftOver dbSUPER_".$spe.".".$hashrawverTOspe{$spe}.".bed ".$hashspesliftoverTOspe{$spe}." dbSUPER_".$hashspeverTOspe{$spe}.".bed unMapped");
	#system("liftOver dbSUPER_".$spe.".bed ".$hashspesliftoverTOspe{$spe}." dbSUPER_".$spe.$hashspeverTOspe{$spe}.".bed unMapped");
}

#dbSUPERraw/".$spe目录下的文件是利用单细胞测序技术鉴定的超级增强子，我的数据是增强子，因此不需要这个文件吧？
#注意文件中的数据是hg19和mm9的
#system("wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz")

#Command <- paste("liftOver",hg19file,"hg19ToHg38.over.chain.gz",file,"unMapped",sep=" ")
#system(Command)
#liftOver dbSUPERpre_mm.bed mm9ToMm10.over.chain.gz dbSUPERpre_mm10.bed unMapped
