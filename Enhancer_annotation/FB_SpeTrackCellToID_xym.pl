open AA,">cellsID.txt";
my $firstdir="download/enhancer";
opendir(DIR,$firstdir) or "can't open the file";
@firstdir=readdir DIR;
@firstdir=grep{$_ ne "." && $_ ne ".."} @firstdir;
foreach $spe (@firstdir){
$speup=uc($spe);   #uc() 函数用于将字符串转换为大写形式
	$seconddir=$firstdir."/".$spe;   #download/enhancer/ZEA_MAYS
	opendir(DIR,$seconddir) or "can't open the file";
	@seconddir=readdir DIR;
	@seconddir=grep{$_ ne "." && $_ ne ".."} @seconddir;
	$rankNum=1;
	foreach $cellfile (@seconddir){
	$cellname=$cellfile;
	$cellname=~ s/\.bed//g;
	print AA $cellname."\t".$spe."\t".$speup.sprintf("%03d",$rankNum)."\n";
	$rankNum++;
	}
}
close AA;

#这是对物种的各个细胞系进行编码统计