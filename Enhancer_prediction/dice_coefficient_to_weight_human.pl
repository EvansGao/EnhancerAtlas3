#system("wget https://hgdownload.soe.ucsc.edu/goldenPath/hs10/bigZips/latest/hg38.chrom.sizes");
#conda install -c bioconda perl-math-round -y
$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
use 5.010;
use Math::Round;
use List::MoreUtils qw(uniq);
use Cwd;
my $pwd = cwd();
$cellname=$pwd;
$cellname=~ s/.*\///g;

$dir="track";
opendir(DIR,$dir);
@track=readdir DIR;
@track=grep{$_ ne "." && $_ ne ".."} @track;#@track是TF-Binding POL2 Histone DHS

#get correlations between different tracks and weight for each track
%hashtracktrackCor=();
open COR,">cormatrix.txt";  #这是两两track计算相似性后的矩阵文件，可以用来画热图
print COR "\t".join("\t",@track)."\n";#这是向cor逐一添加变量@track如TF-Binding
%hashtrackweight=();
$weightsum=0;
for($ii=0;$ii<scalar(@track);$ii++){
	print COR $track[$ii];    #这是向cor循环添加变量TF-Binding POL2 Histone DHS，一排一个，形成两两比较的矩阵，cormatrix.txt可以直接画热图
	$hashtrackweight{$track[$ii]}=0; #赋初始值为0
	for($jj=0;$jj<scalar(@track);$jj++){
	print COR "\t".&Jaccard($track[$ii],$track[$jj]);#这是两两文件之间计算jaccard系数，并且添加到cor中，用\t分隔
	$hashtrackweight{$track[$ii]}+=&Jaccard($track[$ii],$track[$jj]);
	}
	print COR "\n";#每一行cor结尾用\n分隔提行，如此遍历完即5*5矩阵
	$hashtrackweight{$track[$ii]}+=-1; #减1是因为与自己的相似性是1
	$weightsum+=$hashtrackweight{$track[$ii]};#这是将每个track的权重相加
}
close COR;
open WEIGHT,">weight.txt";#这是单个track的权重
foreach $trk (@track){
$hashtrackweight{$trk}=$hashtrackweight{$trk}/$weightsum;
print WEIGHT $trk."\t".$hashtrackweight{$trk}."\n";
}
close WEIGHT;
#for i in  *asw.norm.bed; do for j in  *asw.norm.bed; do bedtools jaccard -a $i -b $j | cut -f3 | tail -n +2 | sed "s/^/$i\t$j\t/" ; done; done >ENH_hs_DHS_Histone_POL2_TF_Binding_jaccard_merge.txt
#从这个计算结果可知,Jaccard是自定义的函数,传入的参数是$trackA,$trackB,返回的值是return $numolap/$numall;

#get consensus peaks and consensus webfile
open CONS,">enh_consensusPre.bed";
$unionbd="";
foreach $trk (@track){#这是将所有的track文件合并
open TRACK,"track/".$trk."/ENH-hs-".$trk.".bed";
while(<TRACK>){
s/\r|\n//g;
@temp=split/\t/,$_;
print CONS $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$trk.":".($temp[2]-$temp[1]).":".($hashtrackweight{$trk}*$temp[3])."\n";#$temp[3]是Si''，即Wt*ScoreAt，$hashtrackweight{$trk}是weight.txt中单个track的权重
}
close TRACK;
$unionbd.=" track/".$trk."/ENH-hs-".$trk.".bed";
}
close CONS;
system("bedtools sort -i enh_consensusPre.bed>enh_consensusSort.bed");
system("bedtools merge -i enh_consensusSort.bed -c 4 -o collapse>enh_consensusMerge.bed");
open CONSMERGE,"enh_consensusMerge.bed";
open CONS,">ENH-hs-".$cellname."pre.bed";
$tracknumstd=round(scalar(@track)/2);
while(<CONSMERGE>){
s/\r|\n//g;
@temp=split/\t/,$_;
	if($temp[3] ne ""){
	@tmptracks=();
	@tracksiginfo=split/\,/,$temp[3];
	$sigarea=0;
		foreach $tracksig (@tracksiginfo){
		@peaktracksig=split/\:/,$tracksig;
		push @tmptracks,$peaktracksig[0];
		$sigarea+=$peaktracksig[1]*$peaktracksig[2];#$sigarea即是Lt*Wt*ScoreAt求和
		}
	@tmptracks=uniq(@tmptracks);
## exists in at least one half of all tracks
		if(scalar(@tmptracks)>=$tracknumstd){
		print CONS $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".($sigarea/($temp[2]-$temp[1]))."\n"; #ENH-hs-".$cellname."pre.bed";
		}#$sigarea/($temp[2]-$temp[1]是Lt*Wt*ScoreAt求和/Lt，$temp[0]."\t".$temp[1]."\t".$temp[2]是sort和merge以后的
			#目的是生成ENH-hs-".$cellname."pre.bed   $temp[2]-$temp[1]是merge以后的peak长度
	}
}
close CONSMERGE;
close CONS;
$hashtmpGSMTotargetfile{$cellname}="ENH-hs-".$cellname."pre.bed";
NORM($cellname,"3");  #由于调用的是$cellname,"3"，$hashtmpGSMTotargetfile{$GSM};因此"ENH-hs-".$cellname."pre.bed";不能删除
unlink("enh_consensusPre.bed"); #$#这里调用NORM生成的文件是ENH-hs-".$cellname.".bed
unlink("enh_consensusSort.bed");#$
unlink("enh_consensusMerge.bed");#$
#unlink("ENH-hs-".$cellname."pre.bed");#$不能删除ENH-hs-".$cellname."pre.bed文件，否则后面无法运行
system("bedtools unionbedg -i".$unionbd.">enhweb_consensusUnion.bed"); #这是把多个track的bed合并在一起  #$unionbd.=" track/".$trk."/ENH-hs-".$trk.".bed";
open WEB,"enhweb_consensusUnion.bed"; #$unionbd这里是明确指定了文件的，所以生成的enhweb_consensusUnion.bed不会有问题
open WEBPRE,">enhweb_consensusPre.bed";#enhweb_consensusUnion.bed跟ENH-hs-".$cellname."pre.bed没有关系
while(<WEB>){
s/\r|\n//g;
@temp=split/\t/,$_;
$sig=0;
	for($ii=0;$ii<scalar(@track);$ii++){
	$sig+=$temp[$ii+3]*$hashtrackweight{$track[$ii]};   #$sig是Wt*ScoreAt？应该是四个数值相加
	}
	if($temp[0] ne ""){
	print WEBPRE $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$sig."\n";
	}
}
close WEB;
close WEBPRE;
system("bedtools intersect -a enhweb_consensusPre.bed -b ENH-hs-".$cellname.".bed -wa>ENH-hs-".$cellname."Web.bed");
system("bedGraphToBigWig ENH-hs-".$cellname."Web.bed hg38.chrom.sizes ENH-hs-".$cellname."Web.bw");
#conda install -c bioconda ucsc-bedgraphtobigwig -y 

##get track enhancers  #@track是TF-Binding POL2 Histone DHS
foreach $trk (@track){
system("bedtools intersect -a track/".$trk."/ENH-hs-".$trk.".bed -b ENH-hs-".$cellname.".bed -wa>ENH-hs-".$cellname.$trk."pre.bed");
$hashtmpGSMTotargetfile{$cellname.$trk}="ENH-hs-".$cellname.$trk."pre.bed";   #与open $GSM,$hashtmpGSMTotargetfile{$GSM};一致才行
NORM($cellname.$trk,"3");  #NORM($cellname,"3");如果调用的$cellname，那么生成的文件是ENH-hs-".$cellname.".bed，所以会报错，这里生成的是ENH-hs-".$cellname.$trk.".bed
system("bedGraphToBigWig ENH-hs-".$cellname.$trk.".bed hg38.chrom.sizes ENH-hs-".$cellname.$trk."Web.bw");
#unlink("ENH-hs-".$cellname.$trk."pre.bed"); #$#由于调用的是$cellname,"3"，$hashtmpGSMTotargetfile{$GSM};因此"ENH-hs-".$cellname."pre.bed";不能删除
}#综合上下代码，没有发现ENH-hs-".$cellname.$trk.".bed出处，是调用NORM后生成的ENH-hs-".$GSM.".bed"，将$GSM替换为$cellname
#有了ENH-hs-".$cellname.$trk.".bed，才能进行第五个文件代码

$duration = time - $start;
print "All is done: $duration s\n";

sub Jaccard()  
{
	my ($trackA,$trackB)=@_;
	#system("cat track/".$trackA."/ENH-hs-".$trackA.".bed track/".$trackB."/ENH-hs-".$trackB.".bed>".$trackA.$trackB."pre.bed");   #这是将两个文件合并
	#system("bedtools sort -i ".$trackA.$trackB."pre.bed>".$trackA.$trackB."sort.bed");
	#system("bedtools merge -i ".$trackA.$trackB."sort.bed -c 4 -o collapse>".$trackA.$trackB."merge.bed");#这是两个文件取并集
	system("bedtools intersect -a track/".$trackA."/ENH-hs-".$trackA.".bed -b track/".$trackB."/ENH-hs-".$trackB.".bed -wao>".$trackA.$trackB."intersect_wao.bed");#这是两个文件取并集
	$fileflag=$trackA.$trackB;
	$fileflag=~ s/\-//g;
	open $fileflag,$trackA.$trackB."intersect_wao.bed";
	$numolap=0;
	$numall=0;
	$length1=0;
	$length2=0;
	while(<$fileflag>){
	s/\r|\n//g;
	@temp=split/\t/,$_;
		if($temp[8] ne ""){
			$numolap+=$temp[8];
	    }       
}
	close $fileflag;
	open $fileflag,"track/".$trackA."/ENH-hs-".$trackA.".bed";
	while(<$fileflag>){
		s/\r|\n//g;
		@temp=split/\t/,$_;
			$length1+=($temp[2]-$temp[1]);   #+1是因为end减去start本身多减了一个碱基
	}
	close $fileflag;

	open $fileflag,"track/".$trackB."/ENH-hs-".$trackB.".bed";
	while(<$fileflag>){
		s/\r|\n//g;
		@temp=split/\t/,$_;
			$length2+=($temp[2]-$temp[1]);
	}
	close $fileflag;
	$numolap=$numolap*2;
	$numall=$length1+$length2;
	return $numolap/$numall;
}
#相当于是bedtools intersect算出overlap，这个是总的相交的碱基数，然后在除以总长度
#子程序需要运行的时候需要调用，这里调用的应该就是jaccard（$trackA,$trackB），即Jaccard($track[$ii],$track[$jj])
#这是改良高老师的jaccard算法，相当于是计算了两个文件之间存在相交的peak对，每一对算1，当然其中有些与多个重叠相交也重叠计算了

sub   NORM()
{   
	my ($GSM,$numtype)=@_;
	@tmpRegions=();
	%hashtmpRegionSig=();
	open $GSM,$hashtmpGSMTotargetfile{$GSM};
	$lengsum=0;
	$sigsum=0;
	while(<$GSM>){
	s/\r|\n//g;
	my @tmp=split/\t/,$_;
		if($tmp[0]=~ /^chr/){
			push @tmpRegions,$tmp[0]."\t".$tmp[1]."\t".$tmp[2];
			$hashtmpRegionSig{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}=$tmp[$numtype];
			$lengsum+=$tmp[2]-$tmp[1];
			$sigsum+=($tmp[2]-$tmp[1])*$tmp[$numtype];
		}
	
	}
	close $GSM;
	$meansig=$sigsum/$lengsum;   #$meansig是加权平均数
	open $GSM,">ENH-hs-".$GSM.".bed";
	foreach $tmpRegion (@tmpRegions){
	print $GSM $tmpRegion."\t".(10*$hashtmpRegionSig{$tmpRegion}/$meansig)."\n";
	}
	close $GSM;
}


