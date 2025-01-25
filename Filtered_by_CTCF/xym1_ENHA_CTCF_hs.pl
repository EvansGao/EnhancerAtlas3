#change1:danRer11/danRer11index
#change2:复制索引文件
#change3：-g 1.70e8  -g 1.70e8
#change4:$hashGSMtocellTypeName{$temp[0]}=$temp[1]

#运行环境是120.79.22.230的sudo su权限下的chipseq环境
$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
use 5.010;
mkdir("preCTCF"); #$#$#$#$#$#$#$
#sudo chown -R xianym:xianym /data/xianym/software
#use strict;
#use warnings;


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
%hashGSMtocellTypeName=();#$#$#$#$#$#$
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
$hashGSMtocellTypeName{$temp[0]}=$temp[1];
$hashGSMdatatype{$temp[0]}=$temp[14];

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
	if(-e $srrs[0].".fastq"){
		print "+++++++++++++++++++++++++\n";
			system("bowtie2 -p 8 -x hg38/hg38index -U ".$srrs[0].".fastq -S ".$srrs[0].".sam"); 
	}else{
		print "-----------------------\n";
			system("bowtie2 -p 8 -x hg38/hg38index -1 ".$srrs[0]."_1.fastq -2 ".$srrs[0]."_2.fastq -S ".$srrs[0].".sam");
			#bowtie2 -x /home/jiang/reference/hg19.index/bowtie2_index/hg19 -U ".$sourename."/".$filename.".fastq -S ".$sourename.".sam"
		}	
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
	if(-e $srrs[$k].".fastq"){
		print "+++++++++++++++++++++++++\n";
			system("bowtie2 -p 8 -x hg38/hg38index -U ".$srrs[$k].".fastq -S ".$srrs[$k].".sam"); 
			#bowtie2 -x /home/jiang/reference/hg19.index/bowtie2_index/hg19 -U ".$sourename."/".$filename.".fastq -S ".$sourename.".sam"
	}else{
		print "-----------------------\n";
			system("bowtie2 -p 8 -x hg38/hg38index -1 ".$srrs[$k]."_1.fastq -2 ".$srrs[$k]."_2.fastq -S ".$srrs[$k].".sam");
			#bowtie2 -x /home/jiang/reference/hg19.index/bowtie2_index/hg19 -U ".$sourename."/".$filename.".fastq -S ".$sourename.".sam"
	}
	if(-e $srrs[$k].".fastq"){unlink($srrs[$k].".fastq");}
	if(-e $srrs[$k]."_1.fastq"){unlink($srrs[$k]."_1.fastq");}
	if(-e $srrs[$k]."_2.fastq"){unlink($srrs[$k]."_2.fastq");}
		if($ctrflag==1){
		system("macs2 callpeak -t ".$srrs[$k].".sam -c ".$hashCTRfile{$ctrGSM}." -f SAM -g hs -n ".$srrs." -q 0.01");####$$#####$srr->$srrs由于我的数据没有对照组GSM，所以不运行这一部分
		#system("macs2 callpeak -t ".$sourename.".sam -c ".$hashwgettocotrolsource{$trackw}.".sam -f SAM -g hs -n ".$sourename." -B -q 0.01");
		}else{
				system("macs2 callpeak -t ".$srrs[$k].".sam -f SAM -g hs -q 0.01 --nomodel --nolambda -n ".$srrs[$k]);
				#system("macs2 callpeak -t ".$sourename.".sam -g hs --nomodel --nolambda -q 0.01 -n ".$sourename);
		}	
		if(-e $srrs[$k]."_peaks.narrowPeak"){
			unlink($srrs[$k]."_summits.bed");
			$srrfiles.=";".$srrs[$k];#$
			open AA,$srrs[$k]."_peaks.narrowPeak";
       		open BB,">preCTCF/".$hashGSMtocellTypeName{$GSM}.".bed";
			while(<AA>){
			chomp($_);
			@temp=split/\s+|\t|\r/,$_;
				if($temp[0]=~ /chr\d+|chrX|chrY/i && (length($temp[0])==4 || length($temp[0])==5) && $temp[6]>0){
					print BB $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[6]."\n";
				}
				}
			close AA;
			close BB;
			#rename $srrs[$k]."_peaks.narrowPeak", "ENH-hs-".$srrs[$k].".bed";
			unlink($srrs[$k].".sam");
			#cp $srrs[$k]."_peaks.narrowPeak", $hashGSMtocellTypeName{$GSM}.".bed";
			#system("cp $srrs[$k]_peaks.narrowPeak preCTCF/$hashGSMtocellTypeName{$GSM}.bed");  #$#$#$#$#$#$
			}else{
			unlink($srrs[$k]."_summits.bed");
			}
		
	unlink($srrs[$k].".sam");
}
}
	#这里用的for循环，没有道理会有死循环啊

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
	$cellTypeName=$hashGSMtocellTypeName{$GSMNonSRR};#$#$#$#$#$#$#$
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
				print BB $temp[$indexchr]."\t".$temp[$indexstart]."\t".$temp[$indexend]."\t".$sig."\n";
			}
		}
		close AA;
		close BB;
			if($version=~ /hg18/i){system("liftOver ".$GSMNonSRR.".bed hg18ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bed unMapped");unlink($GSMNonSRR.".bed");rename $GSMNonSRR."hg38.bed", $GSMNonSRR.".bed";
			}elsif($version=~ /hg19/i){system("liftOver ".$GSMNonSRR.".bed hg19ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bed unMapped");unlink($GSMNonSRR.".bed");rename $GSMNonSRR."hg38.bed", $GSMNonSRR.".bed";
			}######
		system("cp ".$GSMNonSRR.".bed  preCTCF/".$cellTypeName.".bed");
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
		open BB,">preCTCF/".$cellTypeName.".bed";
		while(<AA>){
		s/\r|\n//g;
		@temp=split/\s+|\t/,$_;
		$temp;
			if(!($temp[0]=~ /track/i)){
				if(!($temp[0]=~ /^chr/i)){
				$temp[0]="chr".$temp[0];
				}
				print BB $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[4]."\n";
			}
		}
		close AA;
		close BB;
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
		open BB,">preCTCF/".$cellTypeName.".bed";
		while(<AA>){
		s/\r|\n//g;
		@temp=split/\s+|\t/,$_;
			if(!($temp[0]=~ /track/i)){
				if(!($temp[0]=~ /^chr/i)){
				$temp[0]="chr".$temp[0];
				}
				print BB $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[4]."\n";
			}
		}
		close AA;
		close BB;
		#unlink($GSMNonSRR."pre.bed");
		#unlink($GSMNonSRR.".bed");
		unlink($GSMNonSRR.".bedgraph");
		unlink($GSMNonSRR.".bw");
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
				print BB $temp[$indexchr]."\t".$temp[$indexstart]."\t".$temp[$indexend]."\t".$sig."\n";
			}
		}
		close AA;
		close BB;
			if($version=~ /hg18/i){system("liftOver ".$GSMNonSRR.".bed hg18ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bed unMapped");unlink($GSMNonSRR.".bed");rename $GSMNonSRR."hg38.bed", $GSMNonSRR.".bed";
			}elsif($version=~ /hg19/i){system("liftOver ".$GSMNonSRR.".bed hg19ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bed unMapped");unlink($GSMNonSRR.".bed");rename $GSMNonSRR."hg38.bed", $GSMNonSRR.".bed";
			}######
		system("cp ".$GSMNonSRR.".bed  preCTCF/".$cellTypeName.".bed");
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
			system("macs2 callpeak -t ".$GSMNonSRR.".bam -f BAM -g hs -q 0.01 --nomodel --shift -100 --extsize 200 --keep-dup all -n ".$GSMNonSRR);
		}else{
			system("macs2 callpeak -t ".$GSMNonSRR.".bam -f BAM -g hs -q 0.01 --nomodel --nolambda -n ".$GSMNonSRR);
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
				print BB $temp[$indexchr]."\t".$temp[$indexstart]."\t".$temp[$indexend]."\t".$sig."\n";
			}
		}
		close AA;
		close BB;
			if($version=~ /hg18/i){system("liftOver ".$GSMNonSRR.".bed hg18ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bed unMapped");unlink($GSMNonSRR.".bed");rename $GSMNonSRR."hg38.bed", $GSMNonSRR.".bed";
			}elsif($version=~ /hg19/i){system("liftOver ".$GSMNonSRR.".bed hg19ToHg38.over.chain.gz ".$GSMNonSRR."hg38.bed unMapped");unlink($GSMNonSRR.".bed");rename $GSMNonSRR."hg38.bed", $GSMNonSRR.".bed";
			}######
		system("cp ".$GSMNonSRR.".bed  preCTCF/".$cellTypeName.".bed");
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
				print BB $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[4]."\n";
			}
		}
		close AA;
		close BB;
		system("cp ".$GSMNonSRR.".bed  preCTCF/".$cellTypeName.".bed");
		unlink($GSMNonSRR.".bedgraph");
		unlink($GSMNonSRR.".bw");
	}
}
$duration = time - $start;
print "time processing run Non-SRR experiment files: $duration s\n";

#re-check the wrong files and classify correct files
