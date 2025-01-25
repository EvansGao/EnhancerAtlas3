$start = time;
use List::Util qw[min max sum];
use Statistics::Basic qw<median>;
use List::MoreUtils qw(uniq);
use List::Util qw(shuffle);
system("wget --no-check-certificate https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz");
system("gunzip refGene.txt.gz");
system("mv refGene.txt Refgene_hg38.txt");
open GENE,"Refgene_hg38.txt";
open PROEXON,">standard_promotor_exonpre.bed";
while(<GENE>){
s/\r|\n//g;
@temp=split/\t/,$_;
$temp[9]=~ s/\,$//g;
$temp[10]=~ s/\,$//g;
if($temp[3] eq "+"){
	if(!($temp[2] =~ /G|K|random/i)){
		if($temp[4]-5000 <= 0){
			$startNew =0;
		}else{
			$startNew = $temp[4]-5000;
		}
	print PROEXON $temp[2]."\t".$startNew."\t".($temp[4]+500)."\t1"."\n";#这是提取启动子区
	@exonstarts=split/\,/,$temp[9];
	@exonends=split/\,/,$temp[10];
	for($i=0;$i<scalar(@exonstarts);$i++){
	print PROEXON $temp[2]."\t".$exonstarts[$i]."\t".$exonends[$i]."\t1"."\n";
	}
	}#这是提取外显子区
}elsif($temp[3] eq "-"){
	if(!($temp[2] =~ /G|K|random/i)){
		if($temp[5]-500 <= 0){
			$startNew =0;
		}else{
			$startNew = $temp[5]-500;
		}
	print PROEXON $temp[2]."\t".$startNew."\t".($temp[5]+5000)."\t1"."\n";
	@exonstarts=split/\,/,$temp[9];
	@exonends=split/\,/,$temp[10];
	for($i=0;$i<scalar(@exonstarts);$i++){
	print PROEXON $temp[2]."\t".$exonstarts[$i]."\t".$exonends[$i]."\t1"."\n";
	}
	}
}
}
close GENE;
close PROEXON;
system("bedtools sort -i standard_promotor_exonpre.bed>standard_promotor_exonsort.bed");
system("bedtools merge -i standard_promotor_exonsort.bed -c 4 -o mean>standard_promotor_exon.bed");
unlink("standard_promotor_exonpre.bed");  #这是生成了只有启动子和外显子的文件
unlink("standard_promotor_exonsort.bed");


