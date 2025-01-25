#这个文件整合human code中的各个文件中的代码


$start = time;
use List::MoreUtils qw(uniq);
use List::Util qw(min max sum);
use List::MoreUtils qw(uniq);
# # # my $dir="gene_fpkm/";
# # # opendir(DIR,$dir) or "can't open the file";
# # # @dir=readdir DIR;
# # # @cells=();
# # # foreach $file (@dir){
# # # 	if($file=~ /(.*)\_gene_fpkm\.txt/){
# # # 	push @cells,$1;
# # # 	}
# # # }
#以上这一段是根据fpkm得到相应的文件，但由于我的fpkm和最终的细胞系不匹配，因此，以enh下的文件为准
my $dir = "enh/";   
opendir(DIR, $dir) or die "can't open the file";
my @dir = readdir DIR;
my @cells = grep { $_ ne "." && $_ ne ".." } @dir;
print "@cells\n";
#@cells=("Limb_E11.5","Cerebellum","Heart","Kidney","Brain_E14.5","Spleen","Thymus","3T3-L1","416B","BAT","Bone_marrow","C3H10Thalf","CD19+","CD4+CD8+","CMP","Cortex","ESC_Bruce4","ESC_NPC","Forelimb_E13","G1E-ER4","GMP","Heart_E14.5","HFSC","Intestine","Limb_E14.5","Large_intestine_epithelial","Liver_E14.5","Liver","Lung","MEL","NIH-3T3","NPC","Striatum","Testis","WAT");


#get the gene ensembl information of TSS
open TRANSCRIPT,"human_ensembl.txt";
@genes=();
%hashgenetogenename=();
%hashgenetogenestart=();
%hashgenetogeneend=();
%hashgenetogenestrand=();
%hashgenetochr=();
%hashgenetoTSS=();
while(<TRANSCRIPT>){
chomp($_);
@temp=split/\t/,$_;
if(!exists $hashgenetogenename{$temp[0]}){
$hashgenetogenename{$temp[0]}=$temp[6];
push @genes,$temp[0];
$hashgenetochr{$temp[0]}="chr".$temp[1];
$hashgenetogenestart{$temp[0]}=$temp[2];
$hashgenetogeneend{$temp[0]}=$temp[3];
if($temp[5] eq "1"){
$hashgenetoTSS{$temp[0]}=$temp[2];
$hashgenetogenestrand{$temp[0]}="+";
}elsif($temp[5] eq "-1"){
$hashgenetoTSS{$temp[0]}=$temp[3];
$hashgenetogenestrand{$temp[0]}="-";
}
}
}
close TRANSCRIPT;
#get the gene ensembl information of TSS

#get expression score of gene with special cell
%hashcellgenescore=();

foreach $cell (@cells){
open AA,"gene_fpkm/" . $cell . "_gene_fpkm.txt";
open GENESCORE,">genesigpre".$cell.".bed";
while(<AA>){
chomp($_);
@temp=split/\t/,$_;
$gene=$temp[0];
if(exists $hashgenetochr{$temp[0]}){
	if($temp[1]>0){
	print GENESCORE $hashgenetochr{$temp[0]}."\t".$hashgenetogenestart{$temp[0]}."\t".$hashgenetogeneend{$temp[0]}."\t".$temp[1]."\n";
	}
	if(!exists $hashcellgenescore{$cell."\t".$gene}){
	$hashcellgenescore{$cell."\t".$gene}=$temp[1];
	}elsif(exists $hashcellgenescore{$cell."\t".$gene} && $hashcellgenescore{$cell."\t".$gene}<$temp[1]){
	$hashcellgenescore{$cell."\t".$gene}=$temp[1];
	}	
}

}
close AA;
close GENESCORE;
system("bedtools sort -i genesigpre".$cell.".bed>genesig".$cell.".bed");
}
#get expression score of gene with special cell done
#####################################################<<<

#file format : chr2L start end fpkm

#set the distance range 1MB upstream/downstream the TSS 
open GENE,">genesinfo.bed";
foreach $gene (@genes){
print GENE $hashgenetochr{$gene}."\t".max(0,$hashgenetoTSS{$gene}-1000000)."\t".($hashgenetoTSS{$gene}+1000000)."\t".$gene."\n";
}#注意这里的起始位点是加减1000000
close GENE;
#system("bedtools sort -i genesinfo.bed>generange1MB.bed");
system("bedtools sort -i genesinfo.bed>genesinfosort.bed");
#set the distance range 1MB upstream/downstream the TSS done
#这是从dm的xym_ECb_GeneEnhEnhToCandidates.pl中移植过来的，下文中发现genesinfosort.bed，是generange1MB.bed？
#####################################################>>>

#####################################################<<<
#这一段要不要合并到这个文件，另说吧
$celllink="";
foreach $cell (@cells){
#$celllink.=" /data8/gts/human/allcells/".$cell."/".$cell."NOCTCF.bed";
$celllink.=" enh/".$cell."/".$cell."NOCTCF.bed";
}
system("bedtools unionbedg -i".$celllink.">union.bed");
open UNION,"union.bed";
open UNIONUNION,">unionintegrated.bed";
while(<UNION>){
chomp($_);
@temp=split/\t/,$_;
if($temp[0] ne ""){
	if(($temp[2]-$temp[1])>0){    #$temp[2]>$temp[1]
@newtemp = @temp[ 3 .. $#temp ];
print UNIONUNION $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".join(";",@newtemp)."\n";  #这是将所有有交集的peak的信号合并
}
}
}
close UNION;
close UNIONUNION;
system("cat ".$celllink.">allcellenhspre.bed");  #这是将所有的增强子peak合并，sort，merge
system("bedtools sort -i allcellenhspre.bed>allcellenhssort.bed");
system("bedtools merge -i allcellenhssort.bed -c 4 -o mean>allcellenhsmerge.bed");
$allenhNum=0;
@allenhs=();
open NEWF,"allcellenhsmerge.bed";
open NEWFFF,">allcellenhsmerge2500.bed";
while(<NEWF>){
chomp($_);
@temp=split/\t/,$_;
if(($temp[2]-$temp[1])<=2500  && ($temp[2]-$temp[1])>0){    #保证$temp[2]>$temp[1]
$allenhNum++;
print NEWFFF $_."\n";
push @allenhs,$temp[0]."\t".$temp[1]."\t".$temp[2];
}
}
close NEWF;
close NEWFFF;  #这是将所有的增强子peak合并，sort，merge，并删除end-start>2500的
print "All enhs are ".$allenhNum;
system("bedtools intersect -a allcellenhsmerge.bed -b unionintegrated.bed -wa -wb>allcellunion.bed");
@allenhs=();
%hashEnhScore=();
open BB,"allcellunion.bed";
while(<BB>){
chomp($_);
@temp=split/\t/,$_;
if(!exists $hashEnhScore{$temp[0]."\t".$temp[1]."\t".$temp[2]}){
$hashEnhScore{$temp[0]."\t".$temp[1]."\t".$temp[2]}=$temp[7];
push @allenhs,$temp[0]."\t".$temp[1]."\t".$temp[2];
}else{
$old=$hashEnhScore{$temp[0]."\t".$temp[1]."\t".$temp[2]};
$hashEnhScore{$temp[0]."\t".$temp[1]."\t".$temp[2]}=&MERGE($old,$temp[7]);
#选择大的信号值更大的作为最终的信号值
}
}
close BB;

#以上代码来自xym_1EnhancerInteraction.pl
#####################################################>>>


%hashcellenh=();
foreach $cell (@cells){
	#Integrate cell enhancers into standard format
	system("bedtools intersect -a allcellenhsmerge2500.bed -b enh/".$cell."/".$cell."NOCTCF.bed -wa -wb>allcellenhspre_".$cell.".bed");
	@tempenhs=();
	%hashtempenh=();
	open CELLENH,"allcellenhspre_".$cell.".bed";
	while(<CELLENH>){
	chomp($_);
	@temp=split/\t/,$_;
	if(!exists $hashtempenh{$temp[0]."\t".$temp[1]."\t".$temp[2]}){
	$hashtempenh{$temp[0]."\t".$temp[1]."\t".$temp[2]}=$temp[7];
	$hashcellenh{$cell."\t".$temp[0]."\t".$temp[1]."\t".$temp[2]}=$temp[7];
	push @tempenhs,$temp[0]."\t".$temp[1]."\t".$temp[2];
	}elsif(exists $hashtempenh{$temp[0]."\t".$temp[1]."\t".$temp[2]} && $hashtempenh{$temp[0]."\t".$temp[1]."\t".$temp[2]}<$temp[7]){
	$hashtempenh{$temp[0]."\t".$temp[1]."\t".$temp[2]}=$temp[7];
	$hashcellenh{$cell."\t".$temp[0]."\t".$temp[1]."\t".$temp[2]}=$temp[7];
	}
	}
	close CELLENH;
	unlink("allcellenhspre_".$cell.".bed");
	open CELLSTANDARD,">allenhs_".$cell.".bed";
	foreach $enh (@tempenhs){
	print CELLSTANDARD $enh."\t".$hashtempenh{$enh}."\n";
	}
	close CELLSTANDARD;
}
	#Integrate cell enhancers into standard format done

#get gene's enhancers
system("bedtools intersect -a genesinfosort.bed -b allcellenhsmerge2500.bed -wa -wb>AllgenesenhspairsPre.bed");
open ALLPAIR,"AllgenesenhspairsPre.bed";
open PAIRSCORE,">AllgenesenhspairsScore.bed";
%hashpairtocorrelation=();
@genes=();
%hashgeneEnhs=();
while(<ALLPAIR>){
chomp($_);
@temp=split/\t/,$_;
@signals=split/\;/,$temp[7];
$enh=$temp[4]."\t".$temp[5]."\t".$temp[6];
$gene=$temp[3];
if(!exists $hashgeneEnhs{$gene}){
$hashgeneEnhs{$gene}=$enh;
push @genes,$gene;
}else{
$hashgeneEnhs{$gene}.=";".$enh;
}
$correlation=&COR($enh,$gene);
$pair=$temp[4].":".$temp[5]."-".$temp[6]."_".$temp[3]."\$".$hashgenetogenename{$temp[3]}."\$".$hashgenetochr{$temp[3]}."\$".$hashgenetoTSS{$temp[3]}."\$".$hashgenetogenestrand{$temp[3]};
print PAIRSCORE $pair."\t".$correlation."\n";
$hashpairtocorrelation{$pair}=$correlation;

}
close ALLPAIR;
unlink("AllgenesenhspairsPre.bed");
close PAIRSCORE;


#get gene's enhancers done
########################################################<<<
#get gene's enhancers done
# # # # # # system("rm -R GeneEnhEnh");
# # # # # # mkdir("GeneEnhEnh");
# # # # # # my @processes=();
# # # # # # my $pronum=0;
# # # # # # for($i=0;$i<scalar(@genes);$i++){
# # # # # #    $processes[$pronum]=fork();
# # # # # # 	print $i."\t".$genes[$i]."\n";
# # # # # #    if($processes[$pronum]){   
# # # # # #     warn "launching child $child\n"; 
# # # # # #    }else{
# # # # # #    	sleep(4*int($i/40));
# # # # # #     GENEENHENH($genes[$i]);        # child handles 
# # # # # #     exit 0; 
# # # # # #   }
# # # # # # 	$pronum++;
# # # # # # }

# # # # # # for($k=0;$k<($pronum+1);$k++){
# # # # # # waitpid($processes[$k],0);
# # # # # # }
# # # # # # for($k=0;$k<($pronum+1);$k++){
# # # # # # undef($processes[$k]);
# # # # # # }


# # # # # # $duration = time - $start;
# # # # # # print "All is done: $duration s\n";
#这段注释掉的代码是human原本的代码，由于运行时占用内存空间太大，出现服务器崩溃情况
#因此用下面段代码来替代
##############################################<<<<<
#from mouse

mkdir("GeneEnhEnh"); #为什么这里会注释掉？而且代码没有生成相关的数据
my @processes=();
my $pronum=0;
$piece=120;
for($i=0;$i<scalar(@genes);$i++){
	print $i."\n";
   $processes[$pronum]=fork();
#	print $i."\t".$genes[$i]."\n";
   if($processes[$pronum]){   
#    warn "launching child $child\n"; 
   }else{
#   	sleep(3*int($i/40));
    GENEENHENH($genes[$i]);        # child handles 
    exit 0; 
  }
  
  if(($pronum+1)%$piece==0){
#  	print "pronum is:".($pronum+1)."\n";
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}

  }elsif(($pronum+1)%$piece!=0 && ($pronum+1)==$total_control){
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }
	$pronum++;
}

# # # system("rm -R GeneEnhEnh");
# # # mkdir("GeneEnhEnh");

# # # my @processes = ();
# # # my $pronum = 0;
# # # my $piece = 120;  # 控制并行进程数的变量
# # # my $total_control = scalar(@genes);

# # # for (my $i = 0; $i < $total_control; $i++) {
# # #     print $i . "\n";
# # #     $processes[$pronum] = fork();

# # #     if ($processes[$pronum]) {
# # #         warn "launching child $child\n";
# # #     } else {
# # #         sleep(4 * int($i / 40));
# # #         GENEENHENH($genes[$i]);  # 子进程处理
# # #         exit 0;
# # #     }

# # #     if (($pronum + 1) % $piece == 0 || ($pronum + 1) == $total_control) {
# # #         for (my $k = $piece * int(($pronum + 1) / $piece - 1); $k < ($pronum + 1); $k++) {
# # #             waitpid($processes[$k], 0);
# # #             undef($processes[$k]);
# # #         }
# # #     }

# # #     $pronum++;
# # # }

# # # $duration = time - $start;
# # # print "All is done: $duration s\n";
#这段代码是通过控制$piece大小，控制并行进程
##############################################>>>>>

#from mouse

sub   GENEENHENH()
{   
    my ($gene)=@_;
    %hashEnhEnhCor=();
	$num=int(rand(1000000));
	open $num,">GeneEnhEnh/".$gene.".bed";
	@tempenhs=split/\;/,$hashgeneEnhs{$gene};
if(scalar(@tempenhs)==1){
@temppos=split/\t/,$tempenhs[0];
$pair=$temppos[0].":".$temppos[1]."-".$temppos[2]."_".$gene."\$".$hashgenetogenename{$gene}."\$".$hashgenetochr{$gene}."\$".$hashgenetoTSS{$gene}."\$".$hashgenetogenestrand{$gene};
print $num $pair."\t"."1"."\n";
}else{
		%hashpairotherEnh=();
	for($ii=0;$ii<scalar(@tempenhs);$ii++){
		for($jj=$ii+1;$jj<scalar(@tempenhs);$jj++){
		$hashEnhEnhCor{$tempenhs[$ii]."\t".$tempenhs[$jj]}=&ENHCOR($tempenhs[$ii],$tempenhs[$jj]);  #返回的是增强子和增强子之间的pearson系数
		#求和$hashpairToEnhEnh{$pair}+=$hashEnhEnhCor{$tempenhs[$ii]."\t".$tempenhs[$jj]};
		$hashEnhEnhCor{$tempenhs[$jj]."\t".$tempenhs[$ii]}=$hashEnhEnhCor{$tempenhs[$ii]."\t".$tempenhs[$jj]};
		}
	}

	for($ii=0;$ii<scalar(@tempenhs);$ii++){
	@temppos=split/\t/,$tempenhs[$ii];
	$pair=$temppos[0].":".$temppos[1]."-".$temppos[2]."_".$gene."\$".$hashgenetogenename{$gene}."\$".$hashgenetochr{$gene}."\$".$hashgenetoTSS{$gene}."\$".$hashgenetogenestrand{$gene};
		$hashpairToEnhEnh{$pair}=0;
		for($jj=0;$jj<scalar(@tempenhs);$jj++){
			if($ii!=$jj){
			$hashpairToEnhEnh{$pair}+=$hashEnhEnhCor{$tempenhs[$ii]."\t".$tempenhs[$jj]};  #这是&ENHCOR函数求和的值
			}
		}
	print $num $pair."\t".$hashpairToEnhEnh{$pair}/(scalar(@tempenhs)-1)."\n";  #这是对ENHCOR函数求和的值取平均数
	}
}
    close $num;
}

# # # sub   GENEENHENH()
# # # {   
# # #     my ($gene)=@_;
# # #     %hashEnhEnhCor=();
# # # 	$num=int(rand(1000000));
# # # 	open $num,">GeneEnhEnh/".$gene.".bed";
# # # 	@tempenhs=split/\;/,$hashgeneEnhs{$gene};
# # # if(scalar(@tempenhs)==1){
# # # @temppos=split/\t/,$tempenhs[0];
# # # $pair=$temppos[0].":".$temppos[1]."-".$temppos[2]."_".$gene."\$".$hashgenetogenename{$gene}."\$".$hashgenetochr{$gene}."\$".$hashgenetoTSS{$gene}."\$".$hashgenetogenestrand{$gene};
# # # # # # print $num $pair."\t"."2"."\t"."OnlyOne"."\t"."OnlyOne"."\t"."OnlyOne"."\n";
# # # print $num $pair."\t"."1"."\n";
# # # }else{
# # # 		%hashpairotherEnh=();
# # # 	for($ii=0;$ii<scalar(@tempenhs);$ii++){
# # # 		for($jj=$ii+1;$jj<scalar(@tempenhs);$jj++){
# # # 		$hashEnhEnhCor{$tempenhs[$ii]."\t".$tempenhs[$jj]}=&ENHCOR($tempenhs[$ii],$tempenhs[$jj]);
# # # 		$hashEnhEnhCor{$tempenhs[$jj]."\t".$tempenhs[$ii]}=$hashEnhEnhCor{$tempenhs[$ii]."\t".$tempenhs[$jj]};
# # # 		}
# # # 	}

# # # 	for($ii=0;$ii<scalar(@tempenhs);$ii++){
# # # 		@temppos=split/\t/,$tempenhs[$ii];
# # # 	$pair=$temppos[0].":".$temppos[1]."-".$temppos[2]."_".$gene."\$".$hashgenetogenename{$gene}."\$".$hashgenetochr{$gene}."\$".$hashgenetoTSS{$gene}."\$".$hashgenetogenestrand{$gene};
# # # 		$hashpairToEnhEnh{$pair}=-100;
# # # 		$hashpairotherEnh{$pair}="";
# # # 		for($jj=0;$jj<scalar(@tempenhs);$jj++){
# # # 			if($ii!=$jj && $hashpairToEnhEnh{$pair}<$hashEnhEnhCor{$tempenhs[$ii]."\t".$tempenhs[$jj]}){
# # # 			$hashpairToEnhEnh{$pair}=$hashEnhEnhCor{$tempenhs[$ii]."\t".$tempenhs[$jj]};
# # # 			$hashpairotherEnh{$pair}=$tempenhs[$jj];
# # # 			}
# # # 		}
# # # 	print $num $pair."\t".$hashpairToEnhEnh{$pair}."\t".$hashpairotherEnh{$pair}."\n";
# # # 	}
# # # }
# # #     close $num;
# # # }



#######################################################>>>
%hashpairToEEScore=();
%hashpairToOtherEnh=();
for($k=0;$k<scalar(@genes);$k++){
	print $k."\n";
	open AA,"GeneEnhEnh/".$genes[$k].".bed";
	while(<AA>){
	chomp($_);
	@temp=split/\t/,$_;
	if($temp[0] ne ""){
	$hashpairToEEScore{$temp[0]}=$temp[1];
	$hashpairToOtherEnh{$temp[0]}=$temp[2].":".$temp[3]."-".$temp[4];
	}
	}
	close AA;
	print $k."\t".$genes[$k]."\n";
}

	#1MB upstream/downstream the TSS overlapped with cell enhancers with scores in all samples
foreach $cell (@cells){
#	if(-e "cellAllfeatures".$cell.".bed"){next;}
	system("bedtools intersect -a genesinfosort.bed -b allenhs_".$cell.".bed -wa -wb>cellEnhGeneSigPre".$cell.".bed");
	open ALL,"cellEnhGeneSigPre".$cell.".bed";
	open PAIR,">cellEnhGeneSig".$cell.".bed";
	open WINDOW,">cellEnhGeneWindowpre".$cell.".bed";
	@pairs=();
	%hashpairtoTwoFeatures=();
	%hashpairtoLength=();
	%hashpairtoGene=();
	%hashLocalGenetoEnhs=();
	while(<ALL>){
	chomp($_);
	@temp=split/\t/,$_;
	$pair=$temp[4].":".$temp[5]."-".$temp[6]."_".$temp[3]."\$".$hashgenetogenename{$temp[3]}."\$".$hashgenetochr{$temp[3]}."\$".$hashgenetoTSS{$temp[3]}."\$".$hashgenetogenestrand{$temp[3]};
		if($hashcellgenescore{$cell."\t".$temp[3]}>0){
		push @pairs,$pair;
		$hashpairtoGene{$pair}=$temp[3];
			if(!exists $hashLocalGenetoEnhs{$temp[3]}){
			$hashLocalGenetoEnhs{$temp[3]}=$temp[4].":".$temp[5]."-".$temp[6];
			}else{
			$hashLocalGenetoEnhs{$temp[3]}.=";".$temp[4].":".$temp[5]."-".$temp[6];
			}
		$hashpairtoTwoFeatures{$pair}=$hashpairtocorrelation{$pair}."\t".$hashcellgenescore{$cell."\t".$temp[3]};
			if($hashgenetoTSS{$temp[3]}<$temp[5]){
			print PAIR $pair."\t".$hashpairtocorrelation{$pair}."\t".$hashcellgenescore{$cell."\t".$temp[3]}."\t".$temp[7]."\n";
			print WINDOW $temp[4]."\t".$hashgenetoTSS{$temp[3]}."\t".$temp[5]."\t".$pair."\n";
			$hashpairtoLength{$pair}=$temp[5]-$hashgenetoTSS{$temp[3]};
			}elsif($hashgenetoTSS{$temp[3]}>$temp[6]){
			print PAIR $pair."\t".$hashpairtocorrelation{$pair}."\t".$hashcellgenescore{$cell."\t".$temp[3]}."\t".$temp[7]."\n";
			print WINDOW $temp[4]."\t".$temp[6]."\t".$hashgenetoTSS{$temp[3]}."\t".$pair."\n";
			$hashpairtoLength{$pair}=$hashgenetoTSS{$temp[3]}-$temp[6];
			}
		}
	}
	close ALL;
	unlink("cellEnhGeneSigPre".$cell.".bed");
	close PAIR;
	close WINDOW;
	system("bedtools sort -i cellEnhGeneWindowpre".$cell.".bed>cellEnhGeneWindow".$cell.".bed");
	system("bedtools intersect -a cellEnhGeneWindow".$cell.".bed -b allenhs_".$cell.".bed -wa -wb>cellEnhGeneWindowInEnhpre".$cell.".bed");
	open WINDOWENH,"cellEnhGeneWindowInEnhpre".$cell.".bed";
	%hashWindowInEnh=();
	while(<WINDOWENH>){
	chomp($_);
	@temp=split/\t/,$_;
	if(!exists $hashWindowInEnh{$temp[3]}){
	$hashWindowInEnh{$temp[3]}=($temp[6]-$temp[5])*$temp[7];
	}else{
	$hashWindowInEnh{$temp[3]}+=($temp[6]-$temp[5])*$temp[7];
	}
	
	}
	close WINDOWENH;
	unlink("cellEnhGeneWindowInEnhpre".$cell.".bed");
	system("bedtools intersect -a cellEnhGeneWindow".$cell.".bed -b genesig".$cell.".bed -wa -wb>cellEnhGeneWindowInGenepre".$cell.".bed");
	open WINDOWGENE,"cellEnhGeneWindowInGenepre".$cell.".bed";
	%hashWindowInGene=();
	while(<WINDOWGENE>){
	chomp($_);
	@temp=split/\t/,$_;
	if(!exists $hashWindowInGene{$temp[3]}){
	$hashWindowInGene{$temp[3]}=($temp[6]-$temp[5])*$temp[7];
	}else{
	$hashWindowInGene{$temp[3]}+=($temp[6]-$temp[5])*$temp[7];
	}
	
	}
	close WINDOWGENE;
	unlink("cellEnhGeneWindowInGenepre".$cell.".bed");
	open ALLPAIR,">cellAllfeatures".$cell.".bed";
	%hashpairToLocalEEScore=();
	foreach $pair (@pairs){
		if(!exists $hashWindowInEnh{$pair}){
		$hashWindowInEnh{$pair}=0;
		}
		if(!exists $hashWindowInGene{$pair}){
		$hashWindowInGene{$pair}=0;
		}
		$localgene=$hashpairtoGene{$pair};
		$enhstr=$hashLocalGenetoEnhs{$localgene};
		$otherEnh=$hashpairToOtherEnh{$temp[0]};
		if(exists $hashpairToEEScore{$pair} && ($enhstr=~ /$otherEnh/i || $otherEnh=~ /OnlyOne/i)){
		$hashpairToLocalEEScore{$pair}=$hashpairToEEScore{$pair};
		}
		
		if(exists $hashpairtoLength{$pair}){
		print ALLPAIR $pair."\t".$hashpairtoTwoFeatures{$pair}."\t".$hashpairtoLength{$pair}."\t".($hashWindowInEnh{$pair}/$hashpairtoLength{$pair})."\t".($hashWindowInGene{$pair}/$hashpairtoLength{$pair})."\t".$hashpairToLocalEEScore{$pair}."\n";
		}
	}
	close ALLPAIR;
	
#1MB upstream/downstream the TSS overlapped with general enhancers with scores in all samples done

}

#One purpose:get standard enhancer's score in each cell done
$duration = time - $start;
print "All is done: $duration s\n";


sub   MERGE()
{   
    my ($numstrold,$numstrnew)=@_;
    @olds=split/\;/,$numstrold;
    @news=split/\;/,$numstrnew;
    @recs=();
    $arrlen=scalar(@olds);
    for($i=0;$i<$arrlen;$i++){
    push @recs,max($olds[$i],$news[$i]);
    }
    $numstrrec=join(";",@recs);
    return $numstrrec;
    
}

sub   COR()
{
    my ($enh,$gene)=@_;
	$Ex=0;
	$ExA=0;
	$ExB=0;
	@Enhcellvals=();
	@Genecellvals=();
	foreach $eachcell (@cells){
		if(exists $hashcellenh{$eachcell."\t".$enh}){
		push @Enhcellvals,$hashcellenh{$eachcell."\t".$enh};
		}else{
		push @Enhcellvals,"0";
		}
		if(exists $hashcellgenescore{$eachcell."\t".$gene}){
		push @Genecellvals,$hashcellgenescore{$eachcell."\t".$gene};
		}else{
		push @Genecellvals,"0";
		}
		
	}
	$Enhcellsum = sum @Enhcellvals;
	$Enhcelllen=scalar(@Enhcellvals);
	$Enhcellmean=$Enhcellsum/$Enhcelllen;
	$Genecellsum = sum @Genecellvals;
	$Genecelllen=scalar(@Genecellvals);
	$Genecellmean=$Genecellsum/$Genecelllen;
	for($cori=0;$cori<scalar(@Enhcellvals);$cori++){
	$Ex+=($Enhcellvals[$cori]-$Enhcellmean)*($Genecellvals[$cori]-$Genecellmean);
	$ExA+=($Enhcellvals[$cori]-$Enhcellmean)**2;
	$ExB+=($Genecellvals[$cori]-$Genecellmean)**2;
	}
	if($ExA==0 || $ExB==0){
	return 0;
	}else{
	return $Ex/(sqrt($ExA)*sqrt($ExB));
	}

}


sub   ENHCOR()
{
    my ($enhA,$enhB)=@_;
	$Ex=0;
	$ExA=0;
	$ExB=0;
	@EnhAvals=split/\;/,$hashEnhScore{$enhA};
	@EnhBvals=split/\;/,$hashEnhScore{$enhB};
	$EnhAsum = sum @EnhAvals;
	$EnhAlen=scalar(@EnhAvals);
	$EnhAmean=$EnhAsum/$EnhAlen;
	$EnhBsum = sum @EnhBvals;
	$EnhBlen=scalar(@EnhBvals);
	$EnhBmean=$EnhBsum/$EnhBlen;
	for($cori=0;$cori<scalar(@EnhAvals);$cori++){
	$Ex+=($EnhAvals[$cori]-$EnhAmean)*($EnhBvals[$cori]-$EnhBmean);
	$ExA+=($EnhAvals[$cori]-$EnhAmean)**2;
	$ExB+=($EnhBvals[$cori]-$EnhBmean)**2;
	}
	if($ExA==0 || $ExB==0){
	return 0;
	}else{
	return $Ex/(sqrt($ExA)*sqrt($ExB));
	}

}
$duration = time - $start;
print "all is done: $duration s\n";