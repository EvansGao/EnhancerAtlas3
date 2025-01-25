#/data1/xym/enhanceratlas3.0/CTCFfastq/test/preCTCF/".$cellname.".bed
#/data1/xym/enhanceratlas3.0/CTCFfastq/test/preCTCF/commonall.bed


$start = time;
$dir="/data1/xym/enhanceratlas3.0/CTCFfastq/test/preCTCF";
opendir(DIR,$dir) or "Can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;  
#这是CTCF的文件，需要把commonall.bed放进去

use Cwd;#$
my $pwd = cwd();#$
$cellname=$pwd; #$  #$pwd /data/xianym/software/enhdata/Kc167#$
$cellname=~ s/.*\///g;#$
#这是获取细胞类型的bed文件

#判断细胞类型的文件是否存在于@dir中
=pod
my $a = "some_value";  # 假设要检查的值为 "some_value"
my @d = ("value1", "value2", "some_value", "value3");  # 假设包含值的数组

if (grep { $_ eq $a } @d) {
    print "$a exists in the array \@d\n";
} else {
    print "$a does not exist in the array \@d\n";
}
=cut

open DD,">CTCTsta.txt";
#hs_male_adult_54_years_Peyers_patch_tissue_Combined950.bed
#hs_male_adult_54_years_Peyers_patch_tissue_combined.bed
if (grep { $_ =~ /$cellname/i } @d) {
    print "$cellname exists in the array \@d\n";
	system("bedtools intersect -a ".$cellname."_combined.bed -b /data1/xym/enhanceratlas3.0/CTCFfastq/test/preCTCF/".$cellname.".bed -wa>".$cellname."CTCFpre.bed");
	system("bedtools sort -i ".$cellname."CTCFpre.bed>".$cellname."CTCFpresort.bed");
	unlink($cellname."CTCFpre.bed");
	system("bedtools merge -i ".$cellname."CTCFpresort.bed>".$cellname."CTCF.bed");
	unlink($cellname."CTCFpresort.bed");
	system("bedtools subtract -a ".$cellname."_combined.bed -b ".$cellname."CTCF.bed>".$cellname."NOCTCF.bed");
	unlink($cellname."CTCF.bed");
	print DD $cellname."\t".&fileNum("/data1/xym/enhanceratlas3.0/CTCFfastq/test/preCTCF/".$cellname.".bed")."\t".&fileNum($cellname."_combined.bed")."\t".&fileNum($cellname."NOCTCF.bed")."\n";
} else {
    print "$cellname does not exist in the array \@d\n";
	system("bedtools intersect -a ".$cellname."_combined.bed -b /data1/xym/enhanceratlas3.0/CTCFfastq/test/preCTCF/commonall.bed -wa>".$cellname."CTCFpre.bed");
	system("bedtools sort -i ".$cellname."CTCFpre.bed>".$cellname."CTCFpresort.bed");
	unlink($cellname."CTCFpre.bed");
	system("bedtools merge -i ".$cellname."CTCFpresort.bed>".$cellname."CTCF.bed");
	unlink($cellname."CTCFpresort.bed");
	system("bedtools subtract -a ".$cellname."_combined.bed -b ".$cellname."CTCF.bed>".$cellname."NOCTCF.bed");
	unlink($cellname."CTCF.bed");
	print DD $cellname."\t"."NA"."\t".&fileNum($cellname."_combined.bed")."\t".&fileNum($cell."NOCTCF.bed")."\n";
#	unlink($cell."NOCTCF.bed");
}


close DD;



sub fileNum
{
  my($file)=@_;
  open NUM,$file;
  $num=0;
  while(<NUM>){
  chomp($_);
  if($_ ne ""){
  $num++;
  }
  }
  close NUM;
  return $num;
  
  }
#这个函数是为了计算行数

#思路
#生成单个细胞类型的CTCF文件直接就是fastq文件就行，这部分高老师说每个物种只要10+20种细胞类型就行了，不需要太多。
#让后将所有细胞类型CTCF利用bedtools unionbedg工具合并，然后循环，按照信号值，大于0的累计，总是超过一半的为共识CTCF，没有超过一半弃用
#最终得到的文件就是common文件