#生成文件：EPEP".$spes[$i].".txt";
#生成文件：geneinfo".$spes[$i].".txt"
#it needs to change the format of per EP data before running this pl


$start = time;
use File::Copy;
%hashGene=();
#$spe="human";
#$speup="HS";
#$spelow="hs";

#@builds=("galGal5","danRer10","TAIR10");
#@spes=("chicken","zebrafish","arabidopsis");
#@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
# # # @builds=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
# # # @spes=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
# # # @speups=("CJ","MM","MF","CS","HS","GG","DR","PT","RM","DM","AT","OSJ","ZM");
# # # @spelows=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
#@fullnames=("Callithrix_jacchus","Mus_musculus","Macaca_fascicularis","Chlorocebus_sabaeus","Homo_sapiens","Gallus_gallus","Danio_rerio","Pan_troglodytes","Macaca_mulatta","Drosophila_melanogaster","Arabidopsis_thaliana","Oryza_sativa_Japonica","Zea_mays");

# @builds=("TAIR10","ce11","ChlSab1","danRer11","dm6","ASM31002v1","galGal6","hg38","mm10","ASM251v1","rheMac10","IRGSP1",
# "panTro6","ASM276v2","rn7","sacCer3","ASM294v2","SL3","Spur5","susScr11","TGA4","TryBru","Xtro10","Zm-B73");
# @spes=("Arabidopsis_thaliana","Caenorhabditis_elegans","Chlorocebus_sabaeus","Danio_rerio","Drosophila_melanogaster",
# "Ectocarpus_sp._Ec32","Gallus_gallus","Homo_sapiens","Mus_musculus","Kluyveromyces_lactis","Macaca_mulatta","Oryza_sativa",
#  "Pan_traglodytes","Plasmodium_falciparum","Rattus_norvegicus","Saccharomyces_cerevisiae","Schizosaccharomyces_pombe",
#  "Solanum_lycopersicum","Strongylocentrotus_purpuratus","sus_scrofa","Toxoplasma_gondii","Trypanosoma_brucei_brucei",
#  "Xenopus_tropicalis","Zea_mays");
# @speups=("TAIR10","CE11","CHLSAB1","DANRER11","DM6","ASM31002V1","GALGAL6","HG38","MM10","ASM251V1","RHEMAC10","IRGSP1",
# "PANTRO6","ASM276V2","RN7","SACCER3","ASM294V2","SL3","SPUR5","SUSSCR11","TGA4","TRYBRU","XTRO10","ZM-B73");

# @spelows=("tair10","ce11","chlsab1","danrer11","dm6","asm31002v1","galgal6","hg38","mm10","asm251v1","rhemac10","irgsp1",
# "pantro6","asm276v2","rn7","saccer3","asm294v2","sl3","spur5","susscr11","tga4","trybru","xtro10","zm-b73");

#@builds=("gg","mm","at","dr");
#@spes=("chicken","mouse","arabidopsis","zebrafish");
#@speups=("GG","MM","AT","DR");
#@spelows=("gg","mm","at","dr");
#我的EP数据分别在两个文件目录下，这两个文件目录我不想融合在一起，应为融合后会有点麻烦，不好区分哪些是2.0的数据，哪些是3.0的数据
#but my data from 3.0 and 2.0 had been intergated together before.


@builds=("dm6","hg38","mm10");
@spes=("Drosophila_melanogaster","Homo_sapiens","Mus_musculus");
@speups=("DM6","HG38","MM10");
@spelows=("dm6","hg38","mm10");

mkdir("AllEPs");
for($i=0;$i<scalar(@builds);$i++){
mkdir("AllEPs/".$builds[$i]);
if(-d "download/interaction/".$spes[$i]){
%hashUniprotTOensembl=();
open EE,"webgenes/uniprot".$spes[$i].".txt";
#Transcript stable ID	Gene stable ID	UniProtKB/Swiss-Prot ID	UniProtKB/TrEMBL ID	UniProtKB Gene Name ID
#ENST00000361390	ENSG00000198888	P03886	U5Z754	P03886

while(<EE>){
s/\r|\n//g;
@temp=split/\t/,$_;
$temp[3]=~ s/^\s+|\s+$//g;
$temp[4]=~ s/^\s+|\s+$//g;
	if(!exists $hashUniprotTOensembl{$temp[3]} && $temp[3] ne ""){
	$hashUniprotTOensembl{$temp[3]}=$temp[1];
	}
	if(!exists $hashUniprotTOensembl{$temp[4]} && $temp[4] ne ""){
	$hashUniprotTOensembl{$temp[4]}=$temp[1];
	}
}
close EE;
my $dir="download/interaction/".$spes[$i];
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
@cells=();
%hashnameRank=();
$rankNum=1;
%hashgeneEnhPositions=();
%hashgenetoScore=();
%hashegimap=();
%hashegicell=();

foreach $file (@dir){
$cell=$file;
$cell=~ s/\_interaction\.txt$//g;
open EP, ">AllEPs/".$builds[$i]."/".$cell."_EP.txt";
open CELL,$dir."/".$file;
#chr10:99975803-99977641|8.1293	chr10:100009086-100010746|DNMBP:A0A2X0U4N6|	0.207128214468315

#hs：chr1:264432-265023_ENSG00000237613$FAM138A$chr1$37595$-	0.723116
#mm：chr1:16050554-16051682_ENSMUSG00000067795$4930444P10Rik$chr1$16093325$-	25.946024
#dm：chr2L:1156643-1158046_FBgn0031287$CG4291$chr2L$852730$-1	6.913605
#格式不对啊
######################按照我的格式应该也是可以转换的
while(<CELL>){
s/\r|\n//g;
@tmp=split/\:|\-|\_|\||\t/,$_;
#@tmp=split/\:|\-|\_|\$|\t/,$_;
$tmpensembl="Unkown";
	if(exists $hashUniprotTOensembl{$tmp[8]}){
	$tmpensembl=$hashUniprotTOensembl{$tmp[8]};
	}elsif(!exists $hashUniprotTOensembl{$tmp[8]} && exists $hashUniprotTOensembl{$tmp[7]}){
	$tmpensembl=$hashUniprotTOensembl{$tmp[7]};
	}

print EP $tmp[0].":".$tmp[1]."-".$tmp[2]."_".$tmpensembl."\$".$tmp[7]."\$".$tmp[8]."\$".$tmp[4].":".$tmp[5]."-".$tmp[6]."\t".$tmp[10]."\n";
#chr10:99975803-99977641_Unkown$DNMBP$A0A2X0U4N6$chr10:100009086-100010746	0.207128214468315
#$enh geneID genename entryname $gene score
}
close CELL;
close EP;
}
#$tmp[10]哪里来的？

%hashGene=();
open AA,"geneMapdata".$spes[$i].".txt";  #EB程序生成的
%hash=();
while(<AA>){
chomp($_);
@temp=split/\t/,$_;
if($temp[0] ne ""){
$hashGene{$temp[0]}=$temp[1];
}
}
close AA;

mkdir("enhs");
mkdir("enhs/".$builds[$i]);
my $dir="AllEPs/".$builds[$i];
opendir(DIR,$dir) or "can't open the file";
@dir=readdir DIR;
@dir=grep{$_ ne "." && $_ ne ".."} @dir;
@cells=();
%hashnameRank=();
$rankNum=1;
%hashgeneEnhs=();
@genes=();
%hashgeneEnhPositions=();
%hashgenetoScore=();
open RANK,">cell_number".$builds[$i].".txt";
%hashegimap=();
%hashegicell=();
foreach $file (@dir){
$cell=$file;
$cell=~ s/\_EP\.txt//g;
$celluc=$cell;
$celluc=~s/\-|\_//g;
$celluc=uc($celluc);
$rankNumTwodigtal = sprintf("%02d",$rankNum);
$hashnameRank{$cell}=$speups[$i].$rankNumTwodigtal;
#print $hashnameRank{$cell}."\n";
print RANK $hashnameRank{$cell}."\t".$celluc."\t".$cell."\t"."Cell"."\n";
	open CELL,$dir."/".$file;
	open ENHPOS,">enhs/".$builds[$i]."/".$hashnameRank{$cell};
	$n=1;
	%hashexist=();
	while(<CELL>){
	s/\r|\n//g;
	@tmp=split/\t/,$_;
	@temp=split/\:|\-|\_|\$/,$tmp[0];
	if(!exists $hashGene{$temp[3]}){next;}
#chr3L:4467113-4467853_FBgn0035543$CG15020$chr3L$4421662$	20.445884
	if(!exists $hashexist{$temp[0].":".$temp[1]."-".$temp[2]}){
	$nn=sprintf("%05d",$n);
	$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}=$hashnameRank{$cell}."-".$nn;
	$noChr=$temp[0];
	$noChr=~ s/chr//g;
	print ENHPOS $hashnameRank{$cell}."-".sprintf("%05d",$n)."\t".$noChr.":".$temp[1]."-".$temp[2]."\n";
		if(!exists $hashgeneEnhs{$temp[3]}){
		$hashgeneEnhs{$temp[3]}=$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}.",";
		$hashgeneEnhPositions{$temp[3]}=$noChr.":".$temp[1]."-".$temp[2];
		$hashgenetoScore{$temp[3]}=$tmp[1];
		push @genes,$temp[3];
		}else{
			if(!($hashgeneEnhs{$temp[3]}=~ /$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}\,/)){
			$hashgeneEnhs{$temp[3]}.=$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}.",";
			$hashgeneEnhPositions{$temp[3]}.=",".$noChr.":".$temp[1]."-".$temp[2];
			$hashgenetoScore{$temp[3]}.=",".$tmp[1];
			}
		}
	$n++;
	}else{
	$nn=sprintf("%05d",$n);
	$noChr=$temp[0];
	$noChr=~ s/chr//g;
		if(!exists $hashgeneEnhs{$temp[3]}){
		$hashgeneEnhs{$temp[3]}=$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}.",";
		$hashgeneEnhPositions{$temp[3]}=$noChr.":".$temp[1]."-".$temp[2];
		$hashgenetoScore{$temp[3]}=$tmp[1];
		push @genes,$temp[3];
		}else{
			if(!($hashgeneEnhs{$temp[3]}=~ /$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}\,/)){
			$hashgeneEnhs{$temp[3]}.=$hashexist{$temp[0].":".$temp[1]."-".$temp[2]}.",";
			$hashgeneEnhPositions{$temp[3]}.=",".$noChr.":".$temp[1]."-".$temp[2];
			$hashgenetoScore{$temp[3]}.=",".$tmp[1];
			}
		}
	
	
	
	
	}
	
	if(!exists $hashegimap{$temp[3]} && !exists $hashegicell{$temp[3]}){
	$hashegimap{$temp[3]}=$hashGene{$temp[3]};
	$hashegicell{$temp[3]}=$hashnameRank{$cell};
	}elsif(exists $hashegimap{$temp[3]} && !($hashegicell{$temp[3]}=~ /$hashnameRank{$cell}/)){
	$hashegicell{$temp[3]}.=",".$hashnameRank{$cell};
	}
}
close CELL;
close ENHPOS;
$rankNum++;
}
close RANK;

open DATA,">EPEP".$spes[$i].".txt";
open EGIMAP,">geneinfo".$spes[$i].".txt";
foreach $gene (@genes){
$hashgeneEnhs{$gene}=~ s/\,$//g;
#print DATA $transcript."\t".$hashtranscriptGene{$transcript}."\t".$hashtranscriptEnhs{$transcript}."\t".$hashtranscriptEnhPositions{$transcript}."\n";	
print DATA $gene."\t".$hashgenetoScore{$gene}."\t".$hashgeneEnhs{$gene}."\n";
print EGIMAP $gene."\t".$hashegicell{$gene}."\t".$hashGene{$gene}."\n";
}
close DATA;
close EGIMAP;

}
}
$duration = time - $start;
print "All are done: $duration s\n";