open AA,"GeneHancer.txt";
open BB,">GHpre.bed";
@SNPs=();
%hashSNP=();
while(<AA>){
s/\r|\n//g;
@temp=split/\t/,$_;
@temGHs=split/\=|\;/,$temp[8];
	if($temp[2] eq "Enhancer"){
	print BB $temp[0]."\t".$temp[3]."\t".$temp[4]."\t".$temGHs[1]."\n";
	}
}
close AA;
close BB;
#system("liftOver GHpre.bed hg38Tohg19.over.chain GH_hg19pre.bed unMapped");
system("bedtools sort -i GHpre.bed>GH_hg38.bed");
unlink("./GHpre.bed");


#GeneHancer.txt 是怎么生成的？
#https://www.genecards.org/Guide/DatasetRequestConfirmation
#谷歌账号注册一个，然后contact us
#申请相关数据