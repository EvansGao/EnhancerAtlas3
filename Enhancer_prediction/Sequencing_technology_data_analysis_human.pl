#运行环境是120.79.22.230的sudo su权限下的chipseq环境
$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
use 5.010;
mkdir("Reshs");
mkdir("ReshsCtr");
#sudo chown -R xianym:xianym /data/xianym/software
#use strict;
#use warnings;

##build hash in GSM files
%hashRESfile=(); #$
my $dir="Reshs";
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;   #读目录下所有文件，读完后返回null#$
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
foreach $file (@dir){		#$
if($file=~ /^ENH\-hs/i){
@temp=split/\-|\./,$file;   #$
$temp; #$
$hashRESfile{$temp[2]}=$file;     #原代码是$hashRESfile{$temp[2]}=$temp[3];但$temp[3]是bed文件
}#$temp[2]是GSM，$temp[3]是bed，%hashRESfile{GSM=>bed}
}#Reshs目录下的文件应该是后面的程序生成的文件ENH-hs-gsm??命名是？？$temp[2]}=$temp[3]分别是什么？？

%hashCTRfile=();#$
my $dir="ReshsCtr";    #这是在哪里创建的，能自动创建吗？没有报错，不纠结#$my
opendir(DIR,$dir) or "can't open the file";#将一个目录中的文件名（以及一些其它东西）读入
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
foreach $file (@dir){		#$
if($file=~ /^GSM/i){
@temp=split/\-|\./,$file;		#$
$temp;#$
$hashCTRfile{$temp[0]}=$file;
#print "hashCTRfile:".$file."\n";
}
}

##build the SRR information for GSM 
use Cwd;#$
my $pwd = cwd();#$
$cellname=$pwd; #$  #$pwd /data/xianym/software/enhdata/Kc167#$
$cellname=~ s/.*\///g;#$
open SRR,$cellname."GSMtoSRR.txt";#########################格式是gsm srx srr value  #open CELL,$cellname."St.bed"; #Kcl67St.bed#Kcl67GSMtoSRR.txt
%hashGSMtoSRRs=();#$  #$cellname =~ s/.*\///g;：使用正则表达式替换操作，将$cellname中最后一个斜杠（/）之前的所有内容替换为空。这样做是为了去除路径，只保留目录名。
while(<SRR>){
s/\r|\n//g;
@temp=split/\t/,$_;
$temp;#$
if($temp[0] ne ""){
$hashGSMtoSRRs{$temp[0]}=$temp[2];
print "$temp[0]\n";
print "$temp[2]\n";
}
}
#close SRR;

#use Cwd;  #我将这几行代码转移到了GSMtoSRR那里，因此注释掉了
#$pwd = cwd();#我将这几行代码转移到了GSMtoSRR那里，因此注释掉了
#$cellname=$pwd; #$  #$pwd /data/xianym/software/enhdata/Kc167#我将这几行代码转移到了GSMtoSRR那里，因此注释掉了
#$cellname=~ s/.*\///g;#我将这几行代码转移到了GSMtoSRR那里，因此注释掉了
print "###############";
print "$cellname\n\n\n\n"; #$cellname  Kc167 这里应该有很多以细胞系名称命名的文件 
print "###############";
##get data information
@AllGSM=();#$
@GSMforSRR=();#$
@GSMforSRRctr=();#$
@GSMforNonSRR=();#$
%hashGSMtoNonSRRfile=();#$
%hashExperCtr=();#$
%hashGSMtrack=();#$
%hashGSMdatatype=();#$
%hashGSMtoGSE=();#$
open CELL,$cellname."St.bed";#########################如Kcl67St.bed，格式有十几列
while(<CELL>){
s/\r|\n//g;
@temp=split/\t/,$_;#$
$temp;#$
$temp[0]=~ s/^\s+|\s+$|\,$//g;
$temp[2]=~ s/^\s+|\s+$|\,$//g;
$temp[3]=~ s/^\s+|\s+$|\,$//g;
$temp[9]=~ s/^\s+|\s+$|\,$//g;
$temp[10]=~ s/^\s+|\s+$|\,$//g;
$temp[14]=~ s/^\s+|\s+$|\,$//g;
$temp[15]=~ s/^\s+|\s+$|\,$//g;
push @AllGSM,$temp[0];#$  #push用于向数组末尾添加一个或多个元素
@tmpGSE=split/\,/,$temp[3];#$
$hashGSMtoGSE{$temp[0]}=$tmpGSE[0];
$hashGSMdatatype{$temp[0]}=$temp[14];
	if($temp[10]=~ /ATAC|DNase/i){
	$hashGSMtrack{$temp[0]}="DHS";
	}elsif($temp[10]=~ /MNase/i){
	$hashGSMtrack{$temp[0]}="MNase_seq";
	}elsif($temp[10]=~ /(Cut&Run)|(Cut_Run)|(Cut-Run)/i && $temp[9]=~ /(POLR2A)|(RNA polymerase II)|Pol II|PolII|Pol2|RNAPII/i){
	$hashGSMtrack{$temp[0]}="Cut_Run_POL2";
    }elsif($temp[10]=~ /(Cut&Run)|(Cut_Run)|(Cut-Run)/i && $temp[9]=~ /(p300)/i){
	$hashGSMtrack{$temp[0]}="Cut_Run_P300";######
	}elsif($temp[10]=~ /(Cut&Run)|(Cut_Run)|(Cut-Run)/i && $temp[9]=~ /(H3K4me1)|(H3K27ac)/i){
	$hashGSMtrack{$temp[0]}="Cut_Run_Histone";
	}elsif($temp[10]=~ /(Cut&Run)|(Cut_Run)|(Cut-Run)/i && !($temp[9]=~ /(H3K4me1)|(H3K27ac)|(POLR2A)|(RNA polymerase II)|Pol II|(P300)|PolII|Pol2|RNAPII/i)){
	$hashGSMtrack{$temp[0]}="Cut_Run_TF_Binding";
	}elsif($temp[10]=~ /(Cut&TAG)|(Cut_TAG)|(Cut-TAG)/i && $temp[9]=~ /(POLR2A)|(RNA polymerase II)|Pol II|PolII|Pol2|RNAPII/i){
	$hashGSMtrack{$temp[0]}="Cut_TAG_POL2";
    }elsif($temp[10]=~ /(Cut&TAG)|(Cut_TAG)|(Cut-TAG)/i && $temp[9]=~ /(p300)/i){
	$hashGSMtrack{$temp[0]}="Cut_TAG_P300";######
	}elsif($temp[10]=~ /(Cut&TAG)|(Cut_TAG)|(Cut-TAG)/i && $temp[9]=~ /(H3K4me1)|(H3K27ac)/i){
	$hashGSMtrack{$temp[0]}="Cut_TAG_Histone";
	}elsif($temp[10]=~ /(Cut&TAG)|(Cut_TAG)|(Cut-TAG)/i && !($temp[9]=~ /(H3K4me1)|(H3K27ac)|(POLR2A)|(RNA polymerase II)|Pol II|PolII|Pol2|(P300)|RNAPII/i)){
	$hashGSMtrack{$temp[0]}="Cut_TAG_TF_Binding";
	}elsif($temp[10]=~ /NET/i){
	$hashGSMtrack{$temp[0]}="NET_seq";
	}elsif($temp[10]=~ /MPRA/i){
	$hashGSMtrack{$temp[0]}="MPRA_seq";
	}elsif($temp[10]=~ /ChIPexo/i){
	$hashGSMtrack{$temp[0]}="ChIPexo_seq";
	}elsif($temp[10]=~ /CHIA_PET|ChIA-PET/i){
	$hashGSMtrack{$temp[0]}="CHIA_PET_seq";
	}elsif($temp[10]=~ /HiChIP/i){
	$hashGSMtrack{$temp[0]}="HiChIP_Seq";
	}elsif($temp[10]=~ /ChAR/i){
	$hashGSMtrack{$temp[0]}="ChAR_seq";
	}elsif($temp[10]=~ /CAGE|nAnTi-CAGE/i){
	$hashGSMtrack{$temp[0]}="CAGE_seq";
	}elsif($temp[10]=~ /FAIRE/i){
	$hashGSMtrack{$temp[0]}="FAIRE_seq";
	}elsif($temp[10]=~ /STARR/i){
	$hashGSMtrack{$temp[0]}="STARR_Seq";
	}elsif($temp[10]=~ /GRO/i){
	$hashGSMtrack{$temp[0]}="GRO_Seq";
	}elsif($temp[10]=~ /ChIP/i && $temp[9]=~ /(POLR2A)|(RNA polymerase II)|Pol II|PolII|Pol2|RNAPII/i){
	$hashGSMtrack{$temp[0]}="ChIP_POL2";
    }elsif($temp[10]=~ /ChIP/i && $temp[9]=~ /(p300)/i){
	$hashGSMtrack{$temp[0]}="ChIP_P300";######
	}elsif($temp[10]=~ /ChIP/i && $temp[9]=~ /(H3K4me1)|(H3K27ac)/i){
	$hashGSMtrack{$temp[0]}="ChIP_Histone";
	}elsif($temp[10]=~ /ChIP/i && !($temp[9]=~ /(H3K4me1)|(H3K27ac)|(POLR2A)|(RNA polymerase II)|Pol II|PolII|Pol2|(P300)|RNAPII/i)){
	$hashGSMtrack{$temp[0]}="ChIP_TF_Binding";
	}

	if($temp[15] eq "UsingSRR"){
		push @GSMforSRR,$temp[0];
		if($temp[2]=~ /^GSM/){
		push @GSMforSRRctr,$temp[2];
		$hashExperCtr{$temp[0]}=$temp[2];
		}
	}elsif($temp[15] ne "UsingSRR"){
		push @GSMforNonSRR,$temp[0];
		$hashGSMtoNonSRRfile{$temp[0]}=$temp[15];
	}

}
#close CELL;

#use main::Control;
$duration;#$
$start;#$
$duration = time - $start;
print "time before run SRR ctr files: $duration s\n";
#run SRR ctr files
@GSMforSRRctr=uniq(@GSMforSRRctr);
@GSMforSRRctrTRUE=();#$
foreach $GSMSRRctr (@GSMforSRRctr){#$
	if(!exists $hashCTRfile{$GSMSRRctr}){
	push @GSMforSRRctrTRUE,$GSMSRRctr;
	}
}
#$hashCTRfile{$GSMSRRctr}这是从temp[2]中分选出ReshsCtr文件中不存在的GSM到@GSMforSRRctrTRUE
$piece=6;#$
my @processes=();
my $pronum=0;
$GSMforSRRctrTRUE;#$
for($ii=0;$ii<scalar(@GSMforSRRctrTRUE);$ii++){ #$
   $processes[$pronum]=fork();  #从一个进程中创建两个进程
   if($processes[$pronum]){   
   print "control: ".$GSMforSRRctrTRUE[$ii]."\n";
   }else{ 
    Control($GSMforSRRctrTRUE[$ii]);         # child handles 如果不是零或空，则为 true，否则为 false。
    exit 0; 
  }
  if(($pronum+1)%$piece==0){
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }elsif(($pronum+1)%$piece!=0 && ($pronum+1)==scalar(@GSMforSRRctrTRUE)){
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }
	$pronum++;
}

sub   Control()
{   
	my ($ctrGSM)=@_;
	if(!exists $hashGSMtoSRRs{$ctrGSM}){
	open $ctrGSM,">".$ctrGSM.".wrongNoSRR";
	close $ctrGSM;
	system("mv ".$ctrGSM.".wrongNoSRR ReshsCtr");
	next;
	}
	@srrs=split/\;/,$hashGSMtoSRRs{$ctrGSM};#$  其实是一个GSM对应多个SRR，中间用；隔开，因此，这里是用；分隔后取$srrs[0]，也可以循环跑
	system("fastq-dump --split-3" .$srrs[0]);  #########################可以分成1，2fastq
		#if(-e "/home/tgao7/ncbi/public/sra/".$srrs[0].".sra.cache"){
		#unlink("/home/tgao7/ncbi/public/sra/".$srrs[0].".sra.cache");
		#}elsif(-e "/home/tgao7/ncbi/public/sra/".$srrs[0].".sra"){
		#unlink("/home/tgao7/ncbi/public/sra/".$srrs[0].".sra");
		#}
	if(-e $srrs[0].".fastq"){
		print "+++++++++++++++++++++++++\n";
		if ($hashGSMtrack{$GSM} eq "GRO_Seq" || $hashGSMtrack{$GSM} eq "PRO_seq" || $hashGSMtrack{$GSM} eq "CAGE_seq" || $hashGSMtrack{$GSM} eq "NET_seq") {
			system("bowtie2 --very-sensitive-local -I 1 -X 1000 -p 6 -x ../hg38/rRNAhg38 -U ".$srrs[0].".fastq --un ".$srrs[0]."_rmrRNA.fastq 2>sample_Map2rRNAStat.xls >".$srrs[0]."_rRNA.sam"); 
			system("bwa aln -t 20 ../hg38/hg38.fa ".$srrs[0]."_rmrRNA.fastq > ".$srrs[0]."_rmrRNA.sai ");
			system("bwa samse -n 1 ../hg38/hg38.fa ".$srrs[0]."_rmrRNA.sai ".$srrs[0]."_rmrRNA.fastq > ".$srrs[0]."_rmrRNA.sam");
			system("mv ".$srrs[0]."_rmrRNA.sam ".$srrs[0].".sam");
		}else{
			system("bowtie2 -p 8 -x ../hg38/hg38index -U ".$srrs[0].".fastq -S ".$srrs[0].".sam"); 
		}
	}else{
		print "-----------------------\n";
		if ($hashGSMtrack{$GSM} eq "GRO_Seq" || $hashGSMtrack{$GSM} eq "PRO_seq" || $hashGSMtrack{$GSM} eq "CAGE_seq" || $hashGSMtrack{$GSM} eq "NET_seq") {
			system("bowtie2 --very-sensitive-local -I 1 -X 1000 -p 6 -x ../hg38/rRNAhg38 -1 ".$srrs[0]."_1.fastq -2 ".$srrs[0]."_2.fastq --un ".$srrs[0]."_rmrRNA.fastq 2>".$srrs[0]."_Map2rRNAStat.xls >".$srrs[0]."_rRNA.sam"); 
			system("bwa aln -t 20 ../hg38/hg38.fa  ".$srrs[0]."_rmrRNA.1.fastq >  ".$srrs[0]."_rmrRNA.1.sai"); 
			system("bwa aln -t 20 ../hg38/hg38.fa  ".$srrs[0]."_rmrRNA.2.fastq >  ".$srrs[0]."_rmrRNA.2.sai"); 
			system("bwa samse ../hg38/hg38.fa  ".$srrs[0]."_rmrRNA.1.sai  ".$srrs[0]."_rmrRNA.2.sai  ".$srrs[0]."_rmrRNA.1.fastq  ".$srrs[0]."_rmrRNA.2.fastq >  ".$srrs[0]."_rmrRNA.sam"); 
			system("mv ".$srrs[0]."_rmrRNA.sam ".$srrs[0].".sam");
		}else{
			system("bowtie2 -p 8 -x ../hg38/hg38index -1 ".$srrs[0]."_1.fastq -2 ".$srrs[0]."_2.fastq -S ".$srrs[0].".sam");
		}	
	}
    #my @fastqs=glob("".$srrs[0].".*fastq");
    #my $sizec = @fastqs;
    #my $fastq1;
    #my $fastq2;
    #($fastq1,$fastq2)=@fastqs;
    #if($sizec == 1){
        #print "+++++++++++++++++++++++++\n";
        #system("bowtie2 -p 2 -x ../hg38/hg38index -U ".$srrs[0].".fastq -S ".$ctrGSM.".sam"); 
    #}
    #if($sizec == 2){
        #print "-----------------------\n";
        #system("bowtie2 -p 2 -x ../hg38/hg38index -1 $fastq1 -2 $fastq2 -S ".$ctrGSM.".sam");
    #}
	#system("bowtie2 -p 2 -x ../hg38/hg38index -U ".$srrs[0].".fastq -S ".$ctrGSM.".sam"); ######
	#unlink($srr[0].".fastq");   #bowtie2 -x ../hg38/hg38index -U $dir/$srr/${srr}.fastq -S $dir/$srr/$srr.sam
	if(-e $srrs[0].".fastq"){unlink($srrs[0].".fastq");}
	if(-e $srrs[0]."_1.fastq"){unlink($srrs[0]."_1.fastq");}
	if(-e $srrs[0]."_2.fastq"){unlink($srrs[0]."_2.fastq");}
	if(-e $ctrGSM.".sam"){
	$size = -s $ctrGSM.".sam";
		if($size>10485760){
		system("mv ".$ctrGSM.".sam ReshsCtr");
		}else{
		rename $ctrGSM.".sam", $ctrGSM.".wrongSmallSize";
		system("mv ".$ctrGSM.".wrongSmallSize ReshsCtr");
		}
	}else{
	open $ctrGSM,">".$ctrGSM.".wrongNoExist";
	close $ctrGSM;
	system("mv ".$ctrGSM.".wrongNoExist ReshsCtr");
	}
}#由于我后面收集的数据不存在control数据，因此这部分应该不会再运行

$duration = time - $start;
print "time processing run SRR ctr files: $duration s\n";

#use main::EXPER;
$duration = time - $start;
print "time before run SRR experiment files: $duration s\n";
#run SRR experiment files
@GSMforSRRTRUE=();#￥
foreach $GSMSRR (@GSMforSRR){#$
	if(!exists $hashRESfile{$GSMSRR}){
	push @GSMforSRRTRUE,$GSMSRR;#$
	}
}
$piece=6;  #20
my @processes=();#$
my $pronum=0;#$
for($ii=0;$ii<scalar(@GSMforSRRTRUE);$ii++){#$
   $processes[$pronum]=fork();
   if($processes[$pronum]){   
   print "experiment: ".$GSMforSRRTRUE[$ii]."\n";
   }else{ 
    EXPER($GSMforSRRTRUE[$ii]);         # child handles 
    exit 0; 
  }
  if(($pronum+1)%$piece==0){
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }elsif(($pronum+1)%$piece!=0 && ($pronum+1)==scalar(@GSMforSRRctrTRUE)){
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }
	$pronum++;
}


$duration = time - $start;
print "time processing run SRR experiment files: $duration s\n";

#$hashGSMtoSRRs{$temp[0]}=$temp[2];，$temp[0]是键，$temp[2]是值，!exists $hashGSMtoSRRs{$GSM}表示不存在对照的GSM数据
sub   EXPER()
{   
	my ($GSM)=@_;
	if(!exists $hashGSMtoSRRs{$GSM}){
	open $GSM,">".$GSM.".wrongNoSRR";
	close $GSM;
	system("mv ".$GSM.".wrongNoSRR Reshs");
	#next;
	}
	@srrs=split/\;/,$hashGSMtoSRRs{$GSM};   #将@GSMforSRR中GSM赋值给$GSM，因此@srrs存储的是Kc167GSMtoSRR.txt中的SRR
	$srrs;#$
	$ctrGSM=$hashExperCtr{$GSM};
	$ctrflag=0;
	if(exists $hashCTRfile{$ctrGSM} && !($hashCTRfile{$ctrGSM}=~ /wrong/i)){
	$ctrflag=1;
	}
	$srrfiles="";
	for($k=0;$k<scalar(@srrs);$k++){
	system("fastq-dump --split-3 ".$srrs[$k]);    #这部分才真正开始跑满足srr程序的GSM sample
		#if(-e "/home/tgao7/ncbi/public/sra/".$srrs[0].".sra.cache"){
		#unlink("/home/tgao7/ncbi/public/sra/".$srrs[0].".sra.cache");
		#}elsif(-e "/home/tgao7/ncbi/public/sra/".$srrs[0].".sra"){
		#unlink("/home/tgao7/ncbi/public/sra/".$srrs[0].".sra");
		#}
	if(-e $srrs[$k].".fastq"){
		print "+++++++++++++++++++++++++\n";
		if ($hashGSMtrack{$GSM} eq "GRO_Seq" || $hashGSMtrack{$GSM} eq "PRO_seq" || $hashGSMtrack{$GSM} eq "CAGE_seq") {
			system("bowtie2 --very-sensitive-local -I 1 -X 1000 -p 6 -x ../hg38/rRNAhg38 -U ".$srrs[$k].".fastq --un ".$srrs[$k]."_rmrRNA.fastq 2>sample_Map2rRNAStat.xls >".$srrs[$k]."_rRNA.sam"); 
			system("bwa aln -t 20 ../hg38/hg38.fa ".$srrs[$k]."_rmrRNA.fastq > ".$srrs[$k]."_rmrRNA.sai ");
			system("bwa samse -n 1 ../hg38/hg38.fa ".$srrs[$k]."_rmrRNA.sai ".$srrs[$k]."_rmrRNA.fastq > ".$srrs[$k]."_rmrRNA.sam");
			system("mv ".$srrs[$k]."_rmrRNA.sam ".$srrs[$k].".sam");
		}else{
			system("bowtie2 -p 8 -x ../hg38/hg38index -U ".$srrs[$k].".fastq -S ".$srrs[$k].".sam"); 
		}
	}else{
		print "-----------------------\n";
		if ($hashGSMtrack{$GSM} eq "GRO_Seq" || $hashGSMtrack{$GSM} eq "PRO_seq" || $hashGSMtrack{$GSM} eq "CAGE_seq") {
			system("bowtie2 --very-sensitive-local -I 1 -X 1000 -p 6 -x ../hg38/rRNAhg38 -1 ".$srrs[$k]."_1.fastq -2 ".$srrs[$k]."_2.fastq --un ".$srrs[$k]."_rmrRNA.fastq 2>".$srrs[$k]."_Map2rRNAStat.xls >".$srrs[$k]."_rRNA.sam"); 
			system("bwa aln -t 20 ../hg38/hg38.fa  ".$srrs[$k]."_rmrRNA.1.fastq >  ".$srrs[$k]."_rmrRNA.1.sai"); 
			system("bwa aln -t 20 ../hg38/hg38.fa  ".$srrs[$k]."_rmrRNA.2.fastq >  ".$srrs[$k]."_rmrRNA.2.sai"); 
			system("bwa samse ../hg38/hg38.fa  ".$srrs[$k]."_rmrRNA.1.sai  ".$srrs[$k]."_rmrRNA.2.sai  ".$srrs[$k]."_rmrRNA.1.fastq  ".$srrs[$k]."_rmrRNA.2.fastq >  ".$srrs[$k]."_rmrRNA.sam"); 
			system("mv ".$srrs[$k]."_rmrRNA.sam ".$srrs[$k].".sam");
		}else{
			system("bowtie2 -p 8 -x ../hg38/hg38index -1 ".$srrs[$k]."_1.fastq -2 ".$srrs[$k]."_2.fastq -S ".$srrs[$k].".sam");
		}	
	}
	#system("bowtie2  -p 8 -x hs3 -U ".$srrs[$k].".fastq -S ".$srrs[$k].".sam");
	#unlink($srrs[$k].".fastq");
	if(-e $srrs[$k].".fastq"){unlink($srrs[$k].".fastq");}
	if(-e $srrs[$k]."_1.fastq"){unlink($srrs[$k]."_1.fastq");}
	if(-e $srrs[$k]."_2.fastq"){unlink($srrs[$k]."_2.fastq");}
		if($ctrflag==1){
		system("macs2 callpeak -t ".$srrs[$k].".sam -c ".$hashCTRfile{$ctrGSM}." -f SAM -g hs -n ".$srrs." -q 0.05");####$$#####$srr->$srrs由于我的数据没有对照组GSM，所以不运行这一部分
		}else{
			if($hashGSMtrack{$GSM} eq "DHS" || $hashGSMtrack{$GSM} eq "MPRA_seq"){#只有DHS是例外，其他的track运行都一样
				system("macs2 callpeak -t ".$srrs[$k].".sam -f SAM -g hs -q 0.05 --nomodel --shift -100 --extsize 200 --keep-dup all -n ".$srrs[$k]);
			}else{
				system("macs2 callpeak -t ".$srrs[$k].".sam -f SAM -g hs -q 0.05 --nomodel --nolambda -n ".$srrs[$k]);
			}
		}
		
		if(-e $srrs[$k]."_peaks.narrowPeak"){
		$size = -s $srrs[$k]."_peaks.narrowPeak";
			if($size>=20480){
#			if($size>=0){
			#unlink($srrs[$k]."_peaks.xls");     #注释主要是为了查看循环了几次
			unlink($srrs[$k]."_summits.bed");
			$srrfiles.=";".$srrs[$k];#$
			rename $srrs[$k]."_peaks.narrowPeak", "ENH-hs-".$srrs[$k].".bed";
			unlink($srrs[$k].".sam");
			}else{
			#unlink($srrs[$k]."_peaks.xls");   #注释主要是为了查看循环了几次
			unlink($srrs[$k]."_summits.bed");
			unlink($srrs[$k]."_peaks.narrowPeak");
			}
		}
	unlink($srrs[$k].".sam");
	}
	#这里用的for循环，没有道理会有死循环啊
	
	%hashtmpGSMTotargetfile=();
	%hashtmpAlltofiles=();
	$srrfiles=~ s/\;$//g;
	$hashtmpAlltofiles{$GSM}=$srrfiles;#$
		if($srrfiles eq ""){
		open $GSM,">ENH-hs-".$GSM.".wrong";
		close $GSM;
		system("mv ENH-hs-".$GSM.".wrong Reshs");
		}else{
			if($srrfiles=~ /\;/){
			ASW($GSM,"6");
			$hashtmpGSMTotargetfile{$GSM}=$GSM.".bed";
			NORM($GSM,"3");
			system("mv ENH-hs-".$GSM.".bed Reshs");
			}else{
			$hashtmpGSMTotargetfile{$GSM}="ENH-hs-".$srrfiles.".bed";
			NORM($GSM,"6");
			system("mv ENH-hs-".$GSM.".bed Reshs");
			}
		}
}

#use main::NORM;  #这部分是开始运行非SRR的文件程序
%hashtmpGSMTotargetfile=();
system("wget https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz");
system("wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz");
system("wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chrom.sizes");
system("wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.chrom.sizes");
system("wget https://hgdownload.soe.ucsc.edu/goldenPath/hg18/bigZips/hg18.chrom.sizes");
$duration = time - $start;
print "time before run Non-SRR experiment files: $duration s\n";
#run Non-SRR experiment files
foreach $GSMNonSRR (@GSMforNonSRR){
	if(exists $hashRESfile{$GSMNonSRR} && !($hashRESfile{$GSMNonSRR}=~ /wrong/i)){next;}
	$datatype=$hashGSMdatatype{$GSMNonSRR};#$
	$filename=$hashGSMtoNonSRRfile{$GSMNonSRR};#$
	@datatypeKeys=split/\_/,$datatype;
	$indexchr=$datatypeKeys[1];
	$indexstart=$datatypeKeys[2];
	$indexend=$datatypeKeys[3];
	$indexsig=$datatypeKeys[4];
	$version=$datatypeKeys[5];
	if($datatypeKeys[0] eq "bd"){
		@str = $filename =~ /GSM\d+/ig;
		$str =join("",@str);
		print "$str\n";
		if(length($str)==9){
			print "--------------------------\n\n";
			@str3 =$str =~ /GSM\d{3}/ig;
			$str3 =join("",@str3);
			system("wget https://ftp.ncbi.nlm.nih.gov/geo/samples/${str3}nnn/$str/suppl/$filename");
		}elsif(length($str) ==10){
			print "+++++++++++++++++++++++++++++++++\n\n";
			@str4 =$str =~ /GSM\d{4}/ig;
			$str4 =join("",@str4);
			system("wget https://ftp.ncbi.nlm.nih.gov/geo/samples/${str4}nnn/$str/suppl/$filename");
		}
		if($filename=~ /^$GSMNonSRR.*\.gz$/i){
		 system("gzip -dc $filename > $GSMNonSRR.pre.bed");#删除了gsm/".$GSMNonSRR."/".
		}else{ system("cp $filename $GSMNonSRR.pre.bed");}
		open AA,$GSMNonSRR.".pre.bed";
		open BB,">".$GSMNonSRR.".bed";
		while(<AA>){
		s/\r|\n//g;
		@temp=split/\s+|\t/,$_;
		$temp;
		$sig = 1;
			#if($temp[$indexchr]=~ /1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|X|Y|M/i && !($temp[0]=~ /name|track|\#|browser/i)){
			if($temp[$indexchr]=~ /1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|X|Y/i && !($temp[0]=~ /name|track|\#|browser|GL|JH/i)){
				if(!($temp[$indexchr]=~ /^chr/i)){
				$temp[$indexchr]="chr".$temp[$indexchr];
				}
			##########
			#if(scalar @temp >= 5 && $fields[4] =~ /^\d+$/) {  # 第四列存在且为数值  # 在这里进行处理
			#	$sig=$temp[$indexsig];
			#}elsif(scalar @temp == 3 ){  ## 第四列不存在或者不是数值  # 在这里进行处理
			#	$sig=1;
			#	}
			#}
			###########
			if (scalar @temp == 3) {
                $sig = 1;
            } elsif (scalar @temp == 4) {
                #print "the number of temp is 4\n";
                if ($temp[3] =~ /^(\d+|\d+\.\d+)$/ && $temp[3] != 0) {
                    $sig = $temp[3];
                } else {
                    $sig = 1;
                }
            } elsif (scalar @temp == 5) {
                #print "the number of temp is 5\n";
                if ($temp[4] =~ /^(\d+|\d+\.\d+)$/ && $temp[4] != 0) {
                    $sig = $temp[4];
                } else {
                    $sig = 1;
                }
            } elsif (scalar @temp > 5) {
                #print "the number of temp is 5\n";
                if ($temp[4] =~ /^(\d+|\d+\.\d+)$/ && $temp[4] != 0) {
                    $sig = $temp[4];
                } else {
                    $sig = 1;
                }
            }
				#if($indexsig ne "non"){$sig=$temp[$indexsig];}
				if(($temp[$indexend]-$temp[$indexstart])<=2500){
				print BB $temp[$indexchr]."\t".$temp[$indexstart]."\t".$temp[$indexend]."\t".$sig."\n";
				}
			}
		}
		close AA;
		close BB;
			if($version=~ /hg18/i){system("liftOver ".$GSMNonSRR.".bed hg18ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bed unMapped");unlink($GSMNonSRR.".bed");rename $GSMNonSRR."hg38.bed", $GSMNonSRR.".bed";
			}elsif($version=~ /hg19/i){system("liftOver ".$GSMNonSRR.".bed hg19ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bed unMapped");unlink($GSMNonSRR.".bed");rename $GSMNonSRR."hg38.bed", $GSMNonSRR.".bed";
			}######
		$size = -s $GSMNonSRR.".bed";
		eval {if($size>10240){
		$hashtmpGSMTotargetfile{$GSMNonSRR}=$GSMNonSRR.".bed";
		NORM($GSMNonSRR,"3");
		system("mv ENH-hs-".$GSMNonSRR.".bed Reshs");
		}else{
		rename $GSMNonSRR.".bed", $GSMNonSRR.".wrongSmallSize";
		system("mv ".$GSMNonSRR.".wrongSmallSize Reshs");
		}
		};
		if($@){
		# 捕获到异常
		# 在这里处理异常或输出错误信息
		if($@ =~ /Illegal division by zero/){
			# 处理除以零的异常
			system("mv ".$GSMNonSRR.".Illegalzero Reshs");
			# 继续执行其他适当的操作或退出程序
		}else{
			# 处理其他异常
			print "其他错误：$@\n";
			# 继续执行其他适当的操作或退出程序
			system("mv ".$GSMNonSRR.".formatwrong Reshs");
		}
		};
		unlink($GSMNonSRR.".pre.bed");
		unlink($GSMNonSRR.".bed");
	}elsif($datatypeKeys[0] eq "bg"){
		@str = $filename =~ /GSM\d+/ig;
		$str =join("",@str);
		print "$str\n";
		if(length($str)==9){
			print "--------------------------\n\n";
			@str3 =$str =~ /GSM\d{3}/ig;
			$str3 =join("",@str3);
			system("wget https://ftp.ncbi.nlm.nih.gov/geo/samples/${str3}nnn/$str/suppl/$filename");
		}elsif(length($str) ==10){
			print "+++++++++++++++++++++++++++++++++\n\n";
			@str4 =$str =~ /GSM\d{4}/ig;
			$str4 =join("",@str4);
			system("wget https://ftp.ncbi.nlm.nih.gov/geo/samples/${str4}nnn/$str/suppl/$filename");
		}
		if($filename=~ /\.gz/i){ system("gzip -dc $filename > $GSMNonSRR.bedgraph");
		}else{ system("cp $filename $GSMNonSRR.bedgraph");}
		my $input_file = $GSMNonSRR.".bedgraph";
		my $output_file = $GSMNonSRR.".nohideandtrack.bedgraph";
		open my $input_fh, "<",$input_file;
		open my $output_fh, ">",$output_file;
		while (my $line = <$input_fh>) {
			chomp $line;
			my @fields = split("\t", $line);
			unless ($line =~ /hide/ || $line =~ /track/) {
				print $output_fh "$line\n";
			}
		}
		close($input_fh);
		close($output_fh);
		rename $output_file, "$GSMNonSRR.bedgraph";
			if($version=~ /hg18/i){
			system("liftOver ".$GSMNonSRR.".bedgraph hg18ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bedgraph unMapped");unlink($GSMNonSRR.".bedgraph");rename $GSMNonSRR."hg38.bedgraph", $GSMNonSRR.".bedgraph";
			}elsif($version=~ /hg19/i){
			system("liftOver ".$GSMNonSRR.".bedgraph hg19ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bedgraph unMapped");unlink($GSMNonSRR.".bedgraph");rename $GSMNonSRR."hg38.bedgraph", $GSMNonSRR.".bedgraph";
			}######
		system("macs2 bdgpeakcall -i ".$GSMNonSRR.".bedgraph -c 2 -o ".$GSMNonSRR."pre.bed");
		open AA,$GSMNonSRR."pre.bed";
		open BB,">".$GSMNonSRR.".bed";
		while(<AA>){
		s/\r|\n//g;
		@temp=split/\s+|\t/,$_;
		$temp;
			if(!($temp[0]=~ /track/i)){
				if(!($temp[0]=~ /^chr/i)){
				$temp[0]="chr".$temp[0];
				}
				if(($temp[2]-$temp[1])<=2500){
				print BB $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[4]."\n";
				}
			}
		}
		close AA;
		close BB;
		%hashtmpGSMTotargetfile;#$
		############
		$size = -s $GSMNonSRR.".bed";
		if($size>10240){
		$hashtmpGSMTotargetfile{$GSMNonSRR}=$GSMNonSRR.".bed";
		NORM($GSMNonSRR,"3");
		system("mv ENH-hs-".$GSMNonSRR.".bed Reshs");
		}else{
		rename $GSMNonSRR.".bed", $GSMNonSRR.".wrongSmallSize";
		system("mv ".$GSMNonSRR.".wrongSmallSize Reshs");
		}
		############
		#%hashtmpGSMTotargetfile;#$
		#$hashtmpGSMTotargetfile{$GSMNonSRR}=$GSMNonSRR.".bed";
		#NORM($GSMNonSRR,"3");
		#system("mv ENH-hs-".$GSMNonSRR.".bed Reshs");
		unlink($GSMNonSRR."pre.bed");
		unlink($GSMNonSRR.".bed");
		unlink($GSMNonSRR.".bedgraph");
	}elsif($datatypeKeys[0] eq "bw"){
		@str = $filename =~ /GSM\d+/ig;
		$str =join("",@str);
		print "$str\n";
		if(length($str)==9){
			print "--------------------------\n\n";
			@str3 =$str =~ /GSM\d{3}/ig;
			$str3 =join("",@str3);
			system("wget https://ftp.ncbi.nlm.nih.gov/geo/samples/${str3}nnn/$str/suppl/$filename");
		}elsif(length($str) ==10){
			print "+++++++++++++++++++++++++++++++++\n\n";
			@str4 =$str =~ /GSM\d{4}/ig;
			$str4 =join("",@str4);
			system("wget https://ftp.ncbi.nlm.nih.gov/geo/samples/${str4}nnn/$str/suppl/$filename");
		}
		if($filename=~ /\.gz/i){ system("gzip -dc $filename > $GSMNonSRR.bw");
		}else{ system("cp $filename $GSMNonSRR.bw");}
    	system("bigWigToBedGraph ".$GSMNonSRR.".bw ".$GSMNonSRR.".bedgraph");
		$size = -s $GSMNonSRR.".bedgraph";
		if($size==0){
			system("cp $GSMNonSRR.bedgraph Reshs");		##########################
			next;
			}else{
			if($version=~ /hg18/i){system("liftOver ".$GSMNonSRR.".bedgraph hg18ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bedgraph unMapped");unlink($GSMNonSRR.".bedgraph");rename $GSMNonSRR."hg38.bedgraph", $GSMNonSRR.".bedgraph";
			}elsif($version=~ /hg19/i){system("liftOver ".$GSMNonSRR.".bedgraph hg19ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bedgraph unMapped");unlink($GSMNonSRR.".bedgraph");rename $GSMNonSRR."hg38.bedgraph", $GSMNonSRR.".bedgraph";}
			system("macs2 bdgpeakcall -i ".$GSMNonSRR.".bedgraph -c 2 -o ".$GSMNonSRR."pre.bed");
			#system("macs2 bdgpeakcall -i ".$GSMNonSRR.".bedgraph -o ".$GSMNonSRR."pre.bed");
		} 
		# macs2 bdgbroadcall  -i GSM2948937.bedgraph -o GSM2948937.bdgbroadcall.bed
		open AA,$GSMNonSRR."pre.bed";
		open BB,">".$GSMNonSRR.".bed";
		while(<AA>){
		s/\r|\n//g;
		@temp=split/\s+|\t/,$_;
			if(!($temp[0]=~ /track/i)){
				if(!($temp[0]=~ /^chr/i)){
				$temp[0]="chr".$temp[0];
				}
				if(($temp[2]-$temp[1])<=2500){
				print BB $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[4]."\n";
				}
			}
		}
		close AA;
		close BB;
		$size = -s $GSMNonSRR.".bed";
		if($size>10240){   #if($size>1024){
		$hashtmpGSMTotargetfile{$GSMNonSRR}=$GSMNonSRR.".bed";
		NORM($GSMNonSRR,"3");
		system("mv ENH-hs-".$GSMNonSRR.".bed Reshs");
		}else{
		rename $GSMNonSRR.".bed", $GSMNonSRR.".wrongSmallSize";
		system("mv ".$GSMNonSRR.".wrongSmallSize Reshs");
		}
		#$hashtmpGSMTotargetfile{$GSMNonSRR}=$GSMNonSRR.".bed";
		#NORM($GSMNonSRR,"3");
		#system("mv ENH-hs-".$GSMNonSRR.".bed Reshs");
		unlink($GSMNonSRR."pre.bed");
		unlink($GSMNonSRR.".bed");
		unlink($GSMNonSRR.".bedgraph");
		unlink($GSMNonSRR.".bw");
	}elsif($datatypeKeys[0] eq "txt"){
		@str = $filename =~ /GSM\d+/ig;
		$str =join("",@str);
		print "$str\n";
		if(length($str)==9){
			print "--------------------------\n\n";
			@str3 =$str =~ /GSM\d{3}/ig;
			$str3 =join("",@str3);
			system("wget https://ftp.ncbi.nlm.nih.gov/geo/samples/${str3}nnn/$str/suppl/$filename");
		}elsif(length($str) ==10){
			print "+++++++++++++++++++++++++++++++++\n\n";
			@str4 =$str =~ /GSM\d{4}/ig;
			$str4 =join("",@str4);
			system("wget https://ftp.ncbi.nlm.nih.gov/geo/samples/${str4}nnn/$str/suppl/$filename");
		}
		if($filename=~ /^$GSMNonSRR.*\.gz$/i){
		 system("gzip -dc $filename > $GSMNonSRR.pre.txt");#删除了gsm/".$GSMNonSRR."/".
		}else{ system("cp $filename $GSMNonSRR.pre.txt");}
		my $input_file = $GSMNonSRR.".pre.txt";
		my $output_file = $GSMNonSRR.".pre.bed";
		open my $input_fh, "<",$input_file;
		open my $output_fh, ">",$output_file;
		while (my $line = <$input_fh>) {
			chomp $line;
			if(!($line =~ /name|track|\#|browser/i)){
				if($line =~ /^(\d+)\s+(chr\w+):(\d+)-(\d+)_(\w+)\s+(\d+)$/) {
					my $index = $1;
					my $chromosome = $2;
					my $start = $3;
					my $end = $4;
					my $annotation = $5;
					my $length = $6;
					print $output_fh  "$chromosome\t$start\t$end\t$length\n";
					}elsif($line =~ /^(chr\w+):(\d+)-(\d+)_(\w+)\s+(\d+)$/) {
					my $chromosome = $1;
					my $start = $2;
					my $end = $3;
					my $annotation = $4;
					my $length = $5;
					print $output_fh  "$chromosome\t$start\t$end\t$length\n";
					}elsif ($line =~ /^(\d+):(\d+)-(\d+)\s+(\d+)$/) {
					my $chromosome = "chr$1";
					my $start = $2;
					my $end = $3;
					my $length = $4;
    				print $output_fh "$chromosome\t$start\t$end\t$length\n";
					}elsif ($line =~ /^(chr\w+):(\d+)-(\d+)\s+(\d+)$/) {
					my $chromosome = $1;
					my $start = $2;
					my $end = $3;
					my $length = $4;
    				print $output_fh "$chromosome\t$start\t$end\t$length\n";
					}elsif ($line =~ /^(chr\w+).(\d+).(\d+)\s+(\d+)$/) {
					my $chromosome = $1;
					my $start = $2;
					my $end = $3;
					my $length = $4;
    				print $output_fh "$chromosome\t$start\t$end\t$length\n";
					}
					}
					}			
		close($input_fh);
		close($output_fh);
		open AA,$GSMNonSRR.".pre.bed";
		open BB,">".$GSMNonSRR.".bed";
		while(<AA>){
		s/\r|\n//g;
		@temp=split/\s+|\t/,$_;
		$temp;
			if(!($temp[0]=~ /track/i)){
				if(!($temp[0]=~ /^chr/i)){
				$temp[0]="chr".$temp[0];
				}
				if(($temp[2]-$temp[1])<=2500){
				print BB $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3]."\n";
				}
			}
		}
		close AA;
		close BB;
        if($version=~ /hg18/i){system("liftOver ".$GSMNonSRR.".bed hg18ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bed unMapped");unlink($GSMNonSRR.".bed");rename $GSMNonSRR."hg38.bed", $GSMNonSRR.".bed";
			}elsif($version=~ /hg19/i){system("liftOver ".$GSMNonSRR.".bed hg19ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bed unMapped");unlink($GSMNonSRR.".bed");rename $GSMNonSRR."hg38.bed", $GSMNonSRR.".bed";
			}######
		$size = -s $GSMNonSRR.".bed";
		if($size>10240){
		$hashtmpGSMTotargetfile{$GSMNonSRR}=$GSMNonSRR.".bed";
		NORM($GSMNonSRR,"3");
		system("mv ENH-hs-".$GSMNonSRR.".bed Reshs");
		}else{
		rename $GSMNonSRR.".bed", $GSMNonSRR.".wrongSmallSize";
		system("mv ".$GSMNonSRR.".wrongSmallSize Reshs");
		}
		#unlink($GSMNonSRR.".pre.bed");
		#unlink($GSMNonSRR.".bed");
		#unlink($GSMNonSRR.".pre.txt");
	}elsif($datatypeKeys[0] eq "np"){
		@str = $filename =~ /GSM\d+/ig;
		$str =join("",@str);
		print "$str\n";
		if(length($str)==9){
			print "--------------------------\n\n";
			@str3 =$str =~ /GSM\d{3}/ig;
			$str3 =join("",@str3);
			system("wget https://ftp.ncbi.nlm.nih.gov/geo/samples/${str3}nnn/$str/suppl/$filename");
		}elsif(length($str) ==10){
			print "+++++++++++++++++++++++++++++++++\n\n";
			@str4 =$str =~ /GSM\d{4}/ig;
			$str4 =join("",@str4);
			system("wget https://ftp.ncbi.nlm.nih.gov/geo/samples/${str4}nnn/$str/suppl/$filename");
		}
		if($filename=~ /^$GSMNonSRR.*\.gz$/i){
		 system("gzip -dc $filename > $GSMNonSRR.pre.bed");#删除了gsm/".$GSMNonSRR."/".
		}else{ system("cp $filename $GSMNonSRR.pre.bed");}
		open AA,$GSMNonSRR.".pre.bed";
		open BB,">".$GSMNonSRR.".bed";
		while(<AA>){
		s/\r|\n//g;
		@temp=split/\s+|\t/,$_;
		$temp;
			if(!($temp[0]=~ /name|track|\#|browser|chr\tstart/i)){
				if(!($temp[$indexchr]=~ /^chr/i)){
				$temp[$indexchr]="chr".$temp[$indexchr];
				}
				$sig=1;
				if($indexsig ne "non"){$sig=$temp[$indexsig];}
				if(($temp[$indexend]-$temp[$indexstart])<=2500){
				print BB $temp[$indexchr]."\t".$temp[$indexstart]."\t".$temp[$indexend]."\t".$sig."\n";
				}
			}
		}
		close AA;
		close BB;
			if($version=~ /hg18/i){system("liftOver ".$GSMNonSRR.".bed hg18ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bed unMapped");unlink($GSMNonSRR.".bed");rename $GSMNonSRR."hg38.bed", $GSMNonSRR.".bed";
			}elsif($version=~ /hg19/i){system("liftOver ".$GSMNonSRR.".bed hg19ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bed unMapped");unlink($GSMNonSRR.".bed");rename $GSMNonSRR."hg38.bed", $GSMNonSRR.".bed";
			}######
		$size = -s $GSMNonSRR.".bed";
		if($size>10240){
		$hashtmpGSMTotargetfile{$GSMNonSRR}=$GSMNonSRR.".bed";
		NORM($GSMNonSRR,"3");
		system("mv ENH-hs-".$GSMNonSRR.".bed Reshs");
		}else{
		rename $GSMNonSRR.".bed", $GSMNonSRR.".wrongSmallSize";
		system("mv ".$GSMNonSRR.".wrongSmallSize Reshs");
		}
		unlink($GSMNonSRR.".pre.bed");
		unlink($GSMNonSRR.".bed");
	}elsif($datatypeKeys[0] eq "bm"){
		@str = $filename =~ /GSM\d+/ig;
		$str =join("",@str);
		print "$str\n";
		if(length($str)==9){
			print "--------------------------\n\n";
			@str3 =$str =~ /GSM\d{3}/ig;
			$str3 =join("",@str3);
			system("wget https://ftp.ncbi.nlm.nih.gov/geo/samples/${str3}nnn/$str/suppl/$filename");
		}elsif(length($str) ==10){
			print "+++++++++++++++++++++++++++++++++\n\n";
			@str4 =$str =~ /GSM\d{4}/ig;
			$str4 =join("",@str4);
			system("wget https://ftp.ncbi.nlm.nih.gov/geo/samples/${str4}nnn/$str/suppl/$filename");
		}
		if($filename=~ /^$GSMNonSRR.*\.gz$/i){
		 system("gzip -dc $filename > $GSMNonSRR.bam");#删除了gsm/".$GSMNonSRR."/".
		}else{ system("cp $filename $GSMNonSRR.bam");}
		if($hashGSMtrack{$GSM} eq "DHS" || $hashGSMtrack{$GSM} eq "MPRA_seq"){#只有DHS是例外，其他的track运行都一样
			system("macs2 callpeak -t ".$GSMNonSRR.".bam -f BAM -g hs -q 0.05 --nomodel --shift -100 --extsize 200 --keep-dup all -n ".$GSMNonSRR);
		}else{
			system("macs2 callpeak -t ".$GSMNonSRR.".bam -f BAM -g hs -q 0.05 --nomodel --nolambda -n ".$GSMNonSRR);
		}
		open AA,$GSMNonSRR."_peaks.narrowPeak";
		open BB,">".$GSMNonSRR.".bed";
		while(<AA>){
		s/\r|\n//g;
		@temp=split/\s+|\t/,$_;
		$temp;
			if(!($temp[0]=~ /name|track|\#|browser|chr\tstart/i)){
				if(!($temp[$indexchr]=~ /^chr/i)){
				$temp[$indexchr]="chr".$temp[$indexchr];
				}
				$sig=1;
				if($indexsig ne "non"){$sig=$temp[$indexsig];}
				if(($temp[$indexend]-$temp[$indexstart])<=2500){
				print BB $temp[$indexchr]."\t".$temp[$indexstart]."\t".$temp[$indexend]."\t".$sig."\n";
				}
			}
		}
		close AA;
		close BB;
			if($version=~ /hg18/i){system("liftOver ".$GSMNonSRR.".bed hg18ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bed unMapped");unlink($GSMNonSRR.".bed");rename $GSMNonSRR."hg38.bed", $GSMNonSRR.".bed";
			}elsif($version=~ /hg19/i){system("liftOver ".$GSMNonSRR.".bed hg19ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bed unMapped");unlink($GSMNonSRR.".bed");rename $GSMNonSRR."hg38.bed", $GSMNonSRR.".bed";
			}######
		$size = -s $GSMNonSRR.".bed";
		if($size>10240){
		$hashtmpGSMTotargetfile{$GSMNonSRR}=$GSMNonSRR.".bed";
		NORM($GSMNonSRR,"3");
		system("mv ENH-hs-".$GSMNonSRR.".bed Reshs");
		}else{
		rename $GSMNonSRR.".bed", $GSMNonSRR.".wrongSmallSize";
		system("mv ".$GSMNonSRR.".wrongSmallSize Reshs");
		}
		unlink($GSMNonSRR."_peaks.narrowPeak");
		unlink($GSMNonSRR.".bed");
	}elsif($datatypeKeys[0] eq "wg"){
		@str = $filename =~ /GSM\d+/ig;
		$str =join("",@str);
		print "$str\n";
		if(length($str)==9){
			print "--------------------------\n\n";
			@str3 =$str =~ /GSM\d{3}/ig;
			$str3 =join("",@str3);
			system("wget https://ftp.ncbi.nlm.nih.gov/geo/samples/${str3}nnn/$str/suppl/$filename");
		}elsif(length($str) ==10){
			print "+++++++++++++++++++++++++++++++++\n\n";
			@str4 =$str =~ /GSM\d{4}/ig;
			$str4 =join("",@str4);
			system("wget https://ftp.ncbi.nlm.nih.gov/geo/samples/${str4}nnn/$str/suppl/$filename");
		}
		if($filename=~ /\.gz/i){ system("gzip -dc $filename > $GSMNonSRR.wig");
		}else{ system("cp $filename $GSMNonSRR.wig");}
        #system("wget https://hgdownload.soe.ucsc.edu/goldenPath/$version/bigZips/latest/$version.chrom.sizes"); ######
		system("wigToBigWig ".$GSMNonSRR.".wig $version.chrom.sizes ".$GSMNonSRR.".bw");     ######
		$size = -s $GSMNonSRR.".bw";
		if($size==0){
			system("cp $GSMNonSRR.bw Reshs");         ##########################
			next;
			}else{
    		system("bigWigToBedGraph ".$GSMNonSRR.".bw ".$GSMNonSRR.".bedgraph");
			if($version=~ /hg18/i){system("liftOver ".$GSMNonSRR.".bedgraph hg18ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bedgraph unMapped");unlink($GSMNonSRR.".bedgraph");rename $GSMNonSRR."hg38.bedgraph", $GSMNonSRR.".bedgraph";
			}elsif($version=~ /hg19/i){system("liftOver ".$GSMNonSRR.".bedgraph hg19ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bedgraph unMapped");unlink($GSMNonSRR.".bedgraph");rename $GSMNonSRR."hg38.bedgraph", $GSMNonSRR.".bedgraph";
			}
			system("macs2 bdgpeakcall -i ".$GSMNonSRR.".bedgraph -c 2 -o ".$GSMNonSRR."pre.bed");
		}
		open AA,$GSMNonSRR."pre.bed";
		open BB,">".$GSMNonSRR.".bed";
		while(<AA>){
		s/\r|\n//g;
		@temp=split/\s+|\t/,$_;
			if(!($temp[0]=~ /track/i)){
				if(!($temp[0]=~ /^chr/i)){
				$temp[0]="chr".$temp[0];
				}
				if(($temp[2]-$temp[1])<=2500){
				print BB $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[4]."\n";
				}
			}
		}
		close AA;
		close BB;
		$size = -s $GSMNonSRR.".bed";
		if($size>10240){   #if($size>1024){
		$hashtmpGSMTotargetfile{$GSMNonSRR}=$GSMNonSRR.".bed";
		NORM($GSMNonSRR,"3");
		system("mv ENH-hs-".$GSMNonSRR.".bed Reshs");
		}else{
		rename $GSMNonSRR.".bed", $GSMNonSRR.".wrongSmallSize";
		system("mv ".$GSMNonSRR.".wrongSmallSize Reshs");
		}
		#$hashtmpGSMTotargetfile{$GSMNonSRR}=$GSMNonSRR.".bed";
		#NORM($GSMNonSRR,"3");
		#system("mv ENH-hs-".$GSMNonSRR.".bed Reshs");
		unlink($GSMNonSRR."pre.bed");
		unlink($GSMNonSRR.".bed");
		unlink($GSMNonSRR.".bedgraph");
		unlink($GSMNonSRR.".bw");
	}
}
$duration = time - $start;
print "time processing run Non-SRR experiment files: $duration s\n";

#re-check the wrong files and classify correct files
%hashnewRESfile=();
$newdir="Reshs";
opendir(DIR,$newdir) or "can't open the file";
my @newdir=readdir DIR;
@newdir=grep{$_ ne "." && $_ ne ".."} @newdir;
foreach $file (@newdir){
if($file=~ /^ENH\-hs/i){
@temp=split/\-|\./,$file;
$hashnewRESfile{$temp[2]}=$file; #原本代码$hashnewRESfile{$temp[2]}=$temp[3];  $temp[2]是gsm名称
}
}
open WRONG,">GSMwrongList.bed";
open AVAIL,">GSMavail.bed";
open AVAILTRUE,">GSMavailTrue.bed";
open FINALTRACK,">GSMfinaltrack.bed";
open WRONGTRACK,">GSMwrongtrack.bed";
foreach $gsm (@AllGSM){
	if(exists $hashnewRESfile{$gsm} && $hashnewRESfile{$gsm} eq "bed"){
	print AVAIL $gsm."\t".$hashGSMtrack{$gsm}."\n";
	}else{
	print WRONG $gsm."\t".$hashnewRESfile{$gsm}."\n";
	}
	print AVAILTRUE $gsm."\t".$hashGSMtrack{$gsm}."\n";
	if(exists $hashnewRESfile{$gsm} && $hashnewRESfile{$gsm} =~ /bed/i){
	print FINALTRACK $gsm."\t".$hashGSMtrack{$gsm}."\n";
	}else{
	print WRONGTRACK $gsm."\t".$hashGSMtrack{$gsm}."\n";
	}
}
close WRONG;
close AVAIL;
close AVAILTRUE;
close FINALTRACK;
close WRONGTRACK;


sub   ASW()
{   
	my ($ALL,$numtype)=@_;
	%hashtmpAlltofiles;
	@allshortnames=split/\;/,$hashtmpAlltofiles{$ALL};
	open $ALL,">".$ALL."merge.bed";
	foreach $aswshort (@allshortnames){
	open $aswshort,"ENH-hs-".$aswshort.".bed";
	while(<$aswshort>){
	s/\r|\n//g;
	@tmp=split/\t/,$_;
		if($tmp[0]=~ /^chr/ && ($tmp[2]-$tmp[1])<=2500){
		print $ALL $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".($tmp[2]-$tmp[1])."-".$tmp[$numtype]."\n";
		}
	}
	close $aswshort;
	}
	close $ALL;
	if($aswshort=~ /^SRR/i){unlink("ENH-hs-".$aswshort.".bed");}
	system("bedtools sort -i ".$ALL."merge.bed>".$ALL."sort.bed");
	unlink($ALL."merge.bed");
	system("bedtools merge -i ".$ALL."sort.bed -c 4 -o collapse>".$ALL."pre.bed");
	unlink($ALL."sort.bed");
	open $ALL,$ALL."pre.bed";
	@ALLlines=();
	while(<$ALL>){
	s/\r|\n//g;
	@tmp=split/\t/,$_;
	$center=int(($tmp[1]+$tmp[2])/2);
	@groups=split/\,/,$tmp[3];
	$lengsum=0;
	$sigsum=0;
		foreach $group (@groups){
			@lensiginfo=split/\-/,$group;
			$lengsum+=$lensiginfo[0];
			$sigsum+=$lensiginfo[0]*$lensiginfo[1];
		}
	$meansig=$sigsum/$lengsum;
	$meanlen=$lengsum/scalar(@groups);
	$start=$center-int($meanlen/2);
	$end=$center+int($meanlen/2);
	push @ALLlines,$tmp[0]."\t".$start."\t".$end."\t".$meansig;
	
	}
	close $ALL;
	unlink($ALL."pre.bed");
	open $ALL,">".$ALL.".bed";
	print $ALL join("\n",@ALLlines);
	close $ALL;
}

sub   NORM()
{   
	my ($GSM,$numtype)=@_;  #比如NORM($GSMNonSRR,"3");
	@tmpRegions=();
	%hashtmpRegionSig=();
	open $GSM,$hashtmpGSMTotargetfile{$GSM}; #$hashtmpGSMTotargetfile{$GSMNonSRR}即$GSMNonSRR.bed，即open $GSMNonSRR，$GSMNonSRR.bed
	$lengsum=0;
	$sigsum=0;
	while(<$GSM>){
	s/\r|\n//g;
	@tmp=split/\t|\s+/,$_;
		if($tmp[0]=~ /^chr/ && ($tmp[2]-$tmp[1])<=2500){
			push @tmpRegions,$tmp[0]."\t".$tmp[1]."\t".$tmp[2];
			$hashtmpRegionSig{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}=$tmp[$numtype];
			$lengsum+=$tmp[2]-$tmp[1];
			$sigsum+=($tmp[2]-$tmp[1])*$tmp[$numtype];
		}
	
	}
	close $GSM;
	$meansig=$sigsum/$lengsum;
	open $GSM,">ENH-hs-".$GSM.".bed";
	foreach $tmpRegion (@tmpRegions){
	print $GSM $tmpRegion."\t".(10*$hashtmpRegionSig{$tmpRegion}/$meansig)."\n";
	}
	close $GSM;
}
$duration = time - $start;
print "All are done: $duration s\n";



