#genomebin.bed  #设置的是单一的chr1的chrom.size文件
#standard_bin_promotor_exon.bed  #这是只是chr1的还是整个基因组的？
##################Combined.bed ————————> ENH-hs-".$name.".bed
#运行初始文件
#refGene.txt    #注意在网上下载的文件是refGene.txt，需要把文件名改为refgene.txt
#track
#ENH-hs-hs_male_adult_54_years_Peyers_patch_tissue.bed
#genomebin.bed
#standard_bin_promotor_exon.bed
##xym0_Cutoff_getpeaksnew950.pl
#运行程序后可得到最终的文件是hs_male_adult_54_years_Peyers_patch_tissue_combined.bed


#***************************************HSMMgood1 HUVECgood1,H9good1, HCT-116good1 HEK293good1
use Cwd;
my $pwd = cwd();
$name=$pwd;
$name=~ s/.*\///g;
print $name;   #这是文件目录最后一个目录名
$start = time;
##############open AA,"train_parameter.txt";   
open AA,"weight.txt"; 					#输入文件1
@array=();
%hash=();
while(<AA>){
chomp($_);
@temp=split/\t/,$_;
if(!exists $hash{$temp[0]}){  						#这是创建hash的惯用手法，将0和2作为键：值
push @array,$temp[0];        						#创建数组，temp[0]是trackname
##############$hash{$temp[0]}=$temp[2];   
$hash{$temp[0]}=$temp[1];     						#weight.txt文件中，DHS	0.444863419605257
																	 #ChIP_Histone	0.122071613736933
																	 #ChIP_POL2	0.43306496665781
}
}
#"train_parameter.txt"文件中的$temp[0]是track的名字如TF、Histone、pol2、P300、DHS等
#"train_parameter.txt"文件中的$temp[2]是每个track的权重，所以这个文件应该是weight.txt
close AA;
$meanscore=0;
open CUTOFF,">cutoffres950.txt";
$n=3;
for($i=0;$i<$n;$i++){
$unionbed="";
foreach $trackname (@array){
$flag=0;
###############open INPUT,$trackname."/".$trackname."best.bed";   	#输入文件2   
open INPUT,"track/".$trackname."/ENH-hs-".$trackname.".bed";   	#输入文件2    #$删除了best
	#这是/data1/xym/Species/human/group6/hsfemale_adult_51_years_gastroesophageal_sphincter_tissue/track/ChIP_Histone/ENH-hs-ChIP_Histone.bed
###############open OUTPUT,">".$trackname."/".$trackname."best950.bed";
open OUTPUT,">track/".$trackname."/".$trackname."best950.bed";
while(<INPUT>){
chomp($_);
@temp=split/\t/,$_;
if($temp[0] eq "chr1" && $flag==0){
$flag=1;
print OUTPUT $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3]."\n";
}elsif($temp[0] eq "chr1" && $flag==1){
print OUTPUT $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3]."\n";
}elsif($temp[0] ne "chr1" && $flag==1){
last;  												#last： 立即中止循环,就像C语言中的break；
													#next:并不要求立刻退出循环,但是需要立刻结束当前这次迭代,继续执行循环的下次迭代。
}
}
close INPUT;
close OUTPUT;
#以上的目的是只是提取其中chr1的数据，$trackname."/".$trackname."best.bed是指每个细胞类型或者组织中的鉴定的增强子数据吧？

#system("gawk -P -F '\t' '{if(\$1 == \"chr1\") print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4}' ".$name."//".$trackname."//".$trackname."best995.bed".">".$name."//chr1".$trackname."best995.bed");
system("bedtools shuffle -i track/".$trackname."/".$trackname."best950.bed -g genomebin.bed -excl standard_bin_promotor_exon.bed -chrom -noOverlapping>track/".$trackname."/".$trackname."best950random.bed");
#".$trackname."/".$trackname."best950.bed 文件中只包含了chr1的数据
#-i应在其中放置-i 中的要素的坐标 BED 文件。
#-excl 不应放置来自 -i 的特征（例如，基因组间隙）的坐标 BED 文件。
#-chrom 将 -i 中的特征保留在同一染色体上。仅仅改变它们在染色体上的位置。默认情况下，染色体和位置都是随机选择的。
#-noOverlapping 不允许打乱的间隔重叠。
#-g genomebin.bed应该就是hg38.chrom.size文件,这里我设置为单一的chr的长度文件
#这是随机提取文件数据
system("bedtools sort -i track/".$trackname."/".$trackname."best950random.bed>track/".$trackname."/".$trackname."best950randomsort.bed");
system("bedtools merge -i track/".$trackname."/".$trackname."best950randomsort.bed -c 4 -o mean>track/".$trackname."/".$trackname."best950randommerge.bed");
#-o mean取平均信号值
$duration = time - $start;       
#print $trackname."uniformed: $duration s\n";
$unionbed.=" track/".$trackname."/".$trackname."best950randommerge.bed";
}
system("bedtools unionbedg -i ".$unionbed.">randomall950.bed");
##########foreach $trackname (@array){
###########unlink("track/".$trackname."/".$trackname."best950.bed");
###########unlink("track/".$trackname."/".$trackname."best950random.bed");
###########unlink("track/".$trackname."/".$trackname."best950randomsort.bed");
###########unlink("track/".$trackname."/".$trackname."best950randommerge.bed");
##########}
#这里的trackname是每个细胞类型中的trackname？还是每个细胞类型最终的文件，如果是最终文件，只有一个，就不需要用到unionbedg。但如果不是，其他的文件按照track处理，没有进行算法处理？
#算法是计算每个peak的平均信号值。track是取超过一半track的位点。
#这个cutoff是将每个位点，随机抽取后，通过sort、merge计算平均信号后，对信号进行排序，选择前95%的值。然后取三次，作为最终的cutoff值。
#$trackname应该是每个细胞类型中的trackbed

open BB,"randomall950.bed";
@allScores=();
while(<BB>){
chomp($_);
@temp=split/\t/,$_;
$binscore=0;
for($k=0;$k<scalar(@array);$k++){
	$binscore+=$temp[$k+3]*$hash{$array[$k]};    #这是根据不同track的权重，计算相应的信号值
}
push @allScores,$binscore;

}
close BB;
##########@allScores=sort{$a<=>$b} @allScores;   #对数组 @allScores 中的元素进行排序，排序方式是按数字大小进行升序排列。
#以上代码是将数据从小到大排列，应该是从大到小排？如果是从小到大，得到的值是排在95的值，大值，如果是@allScores=sort{$b<=>$a} @allScores;从大到小，取得值就是排在小的95位值
@allScores=sort{$a<=>$b} @allScores;#$   #对数组 @allScores 中的元素进行排序，排序方式是按数字大小进行升序排列。
$meanscore+=$allScores[int(0.950*scalar(@allScores))];   #取排名int(0.950*scalar(@allScores))索引值
#目的是将数组 @allScores 中排在 95% 位置的元素值加到变量 $meanscore 中。
print CUTOFF $allScores[int(0.950*scalar(@allScores))]."\n";     #这是取信号值
#print $allScores[int(0.95*scalar(@allScores))]."\n";
}
#open CUTOFF,">cutoffres950.txt";文件将记录三个95%cutoff的值

$duration = time - $start;

close CUTOFF;
$cutoff=$meanscore/$n;   #$n=3   $cutoff是三次取值的95%的平均值
#unlink("randomall950.bed");

#以上的目的是将所有的track文件随机提取后，只是为了得到cutoff值？然后确定cutoff值后再在Combined.bed中提取满足条件的peak
##########################################

#################system("bedtools merge -i Combined.bed -c 4 -o max>Combinedmerge.bed");   #Combined.bed???是细胞类型的bed文件吗？
system("bedtools merge -i ENH-hs-".$name.".bed -c 4 -o max>Combinedmerge.bed");   #Combined.bed???是细胞类型的bed文件吗？ENH-hs-".$name.".bed这是前五个步骤得到的细胞类型bed
open CC,"Combinedmerge.bed";
open EE,">".$name."_Combined950.bed";
%hashcutoff=();
while(<CC>){
chomp($_);
@temp=split/\t/,$_;
if($temp[3]>=$cutoff){
	for($n=$temp[1];$n<$temp[2];$n++){
	##########$hashcutoff{$temp[0]."\t".$n}=$temp[0]."\t".($n*10)."\t".(($n+1)*10)."\t".$temp[3];
	$hashcutoff{$temp[0]."\t".$n}=$temp[0]."\t".$n."\t".($n+1)."\t".$temp[3];
	
	}
	###########print EE $temp[0]."\t".($temp[1]*10)."\t".($temp[2]*10)."\t".$temp[3]."\n";
	print EE $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3]."\n";

}#这里为什么start、end都扩大了10倍？
}
close CC;
#unlink("Combinedmerge.bed");
############open DD,"Combined.bed";
open DD,"ENH-hs-".$name.".bed";
open FF,">".$name."_combinedweb950.bed";
while(<DD>){
chomp($_);
@temp=split/\t/,$_;
if(exists $hashcutoff{$temp[0]."\t".$temp[1]}){
	
##########print FF $temp[0]."\t".($temp[1]*10)."\t".($temp[2]*10)."\t".$temp[3]."\n";
print FF $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3]."\n";
#"ENH-hs-".$name.".bed"这是扩展了10倍，与.$name."_Combined950.bed"不同在于信号值，一个时sort merge后的信号值，一个不是。
}

}
close DD;
close FF;
%hashcutoff=();
print $name."\t"."cutoff\t".($cutoff).":".$duration."\n";

#上述代码的目的是将确定好的cutoff值引入，将大于cutoff值得peak保存在$name."_combinedweb950.bed中。
#那么上面从小到大排列是对得？但这样过滤掉得peak会不会有点多？另外，如果是从大到小，过滤的peak应该会少些？
###########################################
use List::Util qw[min max];
open AA,"refgene.txt";
open BB,">refgenePro.bed";
open CC,">generegions.bed";
open DD,">exons.bed";
while(<AA>){
chomp($_);
@temp=split/\t/,$_;

if( !($temp[2]=~ /_/)){
print CC $temp[2]."\t".$temp[4]."\t".$temp[5]."\n";
@exonstarts=split/\,/,$temp[9];
@exonends=split/\,/,$temp[10];
for($i=0;$i<scalar(@exonstarts);$i++){
print DD $temp[2]."\t".$exonstarts[$i]."\t".$exonends[$i]."\n";
}
}

if($temp[3] eq "+" && !($temp[2]=~ /_/)){
print BB $temp[2]."\t".max(1,($temp[4]-10000))."\t".$temp[4]."\t"."ref\n";
}elsif($temp[3] eq "-" && !($temp[2]=~ /_/)){
print BB $temp[2]."\t".$temp[5]."\t".($temp[5]+10000)."\t"."ref\n";
}


#这里的启动子区域是上下游10000？？？？
}
close AA;
close BB;
close CC;
close DD;
system("bedtools sort -i refgenePro.bed>sortrefgenePro.bed");
#unlink("refgenePro.bed");
system("bedtools merge -i sortrefgenePro.bed>mergesortrefgenePro.bed");#启动子文件
#unlink("sortrefgenePro.bed");
system("bedtools sort -i generegions.bed>sortgeneregions.bed");
#unlink("generegions.bed");
system("bedtools merge -i sortgeneregions.bed>mergesortgeneregions.bed");#基因文件
#unlink("sortgeneregions.bed");
system("bedtools sort -i exons.bed>sortexons.bed");
##unlink("exons.bed");
system("bedtools merge -i sortexons.bed>mergesortexons.bed");   #外显子文件，这是部分外显子文件？
#unlink("sortexons.bed");
system("bedtools subtract -a mergesortgeneregions.bed -b mergesortexons.bed>introns.bed");   #内含子文件
#unlink("mergesortexons.bed");
#unlink("mergesortgeneregions.bed");
system("bedtools sort -i introns.bed>sortintrons.bed");
#unlink("introns.bed");
system("bedtools subtract -a mergesortrefgenePro.bed -b sortintrons.bed>refgeneProNointrons.bed");#这是将启动子文件中的内含子区域去除掉
#unlink("sortintrons.bed");
#unlink("mergesortrefgenePro.bed");
system("cat exons.bed refgeneProNointrons.bed>refgeneProExonNointrons.bed");  #这是将外显子和去掉内含子的启动子文件合并，然后sort，merge
system("bedtools sort -i refgeneProExonNointrons.bed>SortrefgeneProExonNointrons.bed");
system("bedtools merge -i SortrefgeneProExonNointrons.bed>refgene10kb.bed");   #这是最终pro+exon-introns的文件refgene10kb.bed
#unlink("refgeneProNointrons.bed");
#unlink("refgeneProExonNointrons.bed");
#unlink("SortrefgeneProExonNointrons.bed");

open CELLNUM,">celltracknum.txt";

system("bedtools subtract -a ".$name."_Combined950.bed -b refgene10kb.bed>".$name."_combinedpre.bed");  
#这是对cutoff后的cellline文件去除启动子、外显子、内含子
system("bedtools sort -i ".$name."_combinedpre.bed>".$name."_combined1.bed");
#unlink($name."_combinedpre.bed");
open CC,$name."_combined1.bed";
open DD,">".$name."_combinedbase.bed";

while(<CC>){
chomp($_);
@temp=split/\t/,$_;
#######print DD $temp[0]."\t".int($temp[1]/10)."\t".int($temp[2]/10)."\t".$temp[3]."\n";
print DD $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3]."\n";
}

#这里将之前start、end有缩小了十倍，但是已经在扩展十倍的基础上对文件启动子、外显子、内含子区域进行删减了，这是什么情况？
close CC;
close DD;
#unlink($name."_combined1.bed");

#以上代码的目的是生成去除了启动子、外显子、内含子等区域并且取值cutoff后的增强子文件$name."_combinedbase.bed，当然还没有去除CTCF
###########################################
#################@trackname=("MNase","DHS","EP300","CAGE","FAIRE","Histone","POL2","TF-Binding","STARR","GROseq","CHIA-PET");
@trackname=("MNase_seq","DHS","Cut_Run_POL2","Cut_Run_Histone","Cut_Run_TF_Binding","Cut_Run_P300","Cut_TAG_P300","Cut_TAG_POL2","Cut_TAG_Histone","Cut_TAG_TF_Binding","NET_seq","MPRA_seq","HiChIP_Seq","ChIP_P300","CAGE_seq","PRO_Seq","FAIRE_seq","ChIP_Histone","ChIP_POL2","ChIP_TF_Binding","STARR_Seq","GRO_Seq","CHIA-PET");
$celltracknum=0;
%hashtracknums=();
foreach $track (@trackname){
#############if (-e "track/".$track."/".$track."best.bed") {
if (-e "track/".$track."/ENH-hs-".$track.".bed") {

#############system("bedtools merge -i track/".$track."/".$track."best.bed -c 4 -o max>merged".$name.$track.".bed");
system("bedtools merge -i track/".$track."/ENH-hs-".$track.".bed -c 4 -o max>merged".$name.$track.".bed");
system("bedtools intersect -a ".$name."_combinedbase.bed -b merged".$name.$track.".bed -wa>intersect".$name.$track.".bed");
	#unlink("merged".$name.$track.".bed");   #这是得到每个track的增强子文件
	if(!(-z "intersect".$name.$track.".bed")){    #文件测试操作符 -z，用于检查文件是否为空
	system("bedtools merge -i intersect".$name.$track.".bed -c 4 -o mean>mergeintersect".$name.$track.".bed");
	#unlink("intersect".$name.$track.".bed");
	$eachtrackNum=&fileNum("mergeintersect".$name.$track.".bed");   #函数&fileNum得到的结果是文件的行数
	$celltracknum += $eachtrackNum;
	open AA,"mergeintersect".$name.$track.".bed";
	while(<AA>){
	chomp($_);
	if($_ ne ""){
	if(!exists $hashtracknums{$_}){
	$hashtracknums{$_}=$eachtrackNum;
	}else{
	$hashtracknums{$_} +=$eachtrackNum;
	}
	}

	}
	close AA;
	#unlink("mergeintersect".$name.$track.".bed");
	
	}

##unlink("intersect".$name.$track.".bed");


}
}
open TRACK,$name."_combinedbase.bed";
open TRACKZOOMIN,">track".$name."_combinedbase.bed";
open ZOOMIN,">".$name."_combined.bed";     #所以这个文件是最终的文件？？？？？？

$stanum=0;
print $celltracknum."\n";

while(<TRACK>){
chomp($_);
if(exists $hashtracknums{$_} && ($hashtracknums{$_}>=int($celltracknum*0.5))){
print CELLNUM $name."\t".$hashtracknums{$_}."\t".$celltracknum."\n";
@temp=split/\t/,$_;
$stanum++;
#####if(($temp[2]*10-$temp[1]*10)>50){     #意思是筛选出peak长度大于5的
if(($temp[2]-$temp[1])>5){     #意思是筛选出peak长度大于5的
print TRACKZOOMIN $_."\n";
###########print ZOOMIN $temp[0]."\t".int($temp[1]*10)."\t".int($temp[2]*10)."\t".$temp[3]."\n";
print ZOOMIN $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3]."\n";
}
#但为什么最终的增强子信息，start，end都扩展了十倍？
}

}
close TRACK;



#unlink($name."_combinedbase.bed");
close TRACKZOOMIN;
close ZOOMIN;
print $stanum;



############system("bedtools intersect -a Combined.bed -b track".$name."_combinedbase.bed>".$name."_combinedweb10.bed");
system("bedtools intersect -a ENH-hs-".$name.".bed -b track".$name."_combinedbase.bed>".$name."_combinedweb10.bed");
#unlink("track".$name."_combinedbase.bed");
open EE,$name."_combinedweb10.bed";
open FF,">".$name."_combinedweb.bed";

while(<EE>){
chomp($_);
@temp=split/\t/,$_;
###########print FF $temp[0]."\t".($temp[1]*10)."\t".($temp[2]*10)."\t".$temp[3]."\n";
print FF $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3]."\n";
}   #为什么上传网页的数据还要再次扩展十倍？
close EE;
#unlink($name."_combinedweb10.bed");
close FF;
close CELLNUM;
$duration = time - $start;
print $name."Getpeaksnew.pl is done: $duration s\n";


sub fileNum
{
  my($file)=@_;
  open NUM,$file;
  $linenum=0;
  while(<NUM>){
  chomp($_);
  if($_ ne ""){
  $linenum++;
  }
  }
  close NUM;
  return $linenum;
  
  }


#genomebin.bed  #设置的是单一的chr1的chrom.size文件
#standard_bin_promotor_exon.bed  #这是只是chr1的还是整个基因组的？
##################Combined.bed ————————> ENH-hs-".$name.".bed
#运行初始文件
#refGene.txt    #注意在网上下载的文件是refGene.txt，需要把文件名改为refgene.txt
#track
#ENH-hs-hs_male_adult_54_years_Peyers_patch_tissue.bed
#genomebin.bed
#standard_bin_promotor_exon.bed
##xym0_Cutoff_getpeaksnew950.pl
#运行程序后可得到最终的文件是hs_male_adult_54_years_Peyers_patch_tissue_combined.bed