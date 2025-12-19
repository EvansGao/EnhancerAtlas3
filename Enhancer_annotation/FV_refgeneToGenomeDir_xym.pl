if(-d "./Genome"){
system("rm -R Genome");
}
mkdir("Genome");
#@spes=("ss","sc","gg","dr","ce","rn","dm");
#@spevers=("susScr3","sacCer3","galGal4","danRer10","ce10","rn5","dm3");
#@spevers=("susScr3","sc3","gg4","dr10","ce10","rn5","dm3");

# # # @spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
# # # @spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");


@spes=("Arabidopsis_thaliana","Caenorhabditis_elegans","Chlorocebus_sabaeus","Danio_rerio","Drosophila_melanogaster",
	"Ectocarpus_sp._Ec32","Gallus_gallus","Homo_sapiens","Mus_musculus","Kluyveromyces_lactis","Macaca_mulatta","Oryza_sativa",
	"Pan_traglodytes","Plasmodium_falciparum","Rattus_norvegicus","Saccharomyces_cerevisiae","Schizosaccharomyces_pombe",
	"Solanum_lycopersicum","Strongylocentrotus_purpuratus","sus_scrofa","Toxoplasma_gondii","Trypanosoma_brucei_brucei",
	"Xenopus_tropicalis","Zea_mays");
@spevers=("TAIR10","ce11","ChlSab1","danRer11","dm6","ASM31002v1","galGal6","hg38","mm10","ASM251v1","rheMac10","IRGSP1",
	"panTro5","ASM276v2","rn7","sacCer3","ASM294v2","SL3","Spur5","susScr11","TGA4","TryBru","Xtro10","Zm-B73");

@speschromsize=("Arabidopsis_thaliana.TAIR10","ce11","Chlorocebus_sabaeus.ChlSab1","danRer11","dm6","Ectocarpus_siliculosus.ASM31002v1",
	"galGal6","hg38","mm10","Kluyveromyces_lactis.ASM251v1","rheMac10","Oryza_sativa.IRGSP1","panTro5","Plasmodium_falciparum.ASM276v2",
	"rn7","sacCer3","Schizosaccharomyces_pombe.ASM294v2","Solanum_lycopersicum.SL3","Strongylocentrotus_purpuratus.Spur5",
	"susScr11","Toxoplasma_gondii.TGA4","Trypanosoma_brucei.TryBru","Xenopus_tropicalis.Xtro10","Zea_mays.Zm-B73");


for($i=0;$i<scalar(@spes);$i++){
%hashchr=();
open AA,"../project/chromsizes/".$speschromsize[$i].".chrom.sizes";
while(<AA>){
s/\r|\n//g;
@tmp=split/\t/,$_;
	if($tmp[0] ne "" && !($tmp[0]=~ /random|chrUn/i)){
	$hashchr{$tmp[0]}="";
	}
}
close AA;

open CC,">".$spevers[$i].".bed";
if(-e "../project/Refgene/Refgene_".$spevers[$i].".txt"){
open BB,"../project/Refgene/Refgene_".$spevers[$i].".txt";
while(<BB>){
s/\r|\n//g;
@tmp=split/\t/,$_;
	if(exists $hashchr{$tmp[2]}){
	print CC $tmp[2]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[3]."\t".$tmp[12]."\t".$tmp[9]."\t".$tmp[10]."\n";
	#     	  chrom	 	genestart	   geneend	     strand	     genename	   exonstart	exonend
	}
}
#1089	NM_001112689	chr13	-	66142314	66148434	66142354	66148418	6	66142314,66144461,66144732,66147081,66148118,66148413,	66142517,66144649,66144891,66147235,66148196,66148434,	0	CIDEC	cmpl	cmpl	2,0,0,2,2,0,
#[编号]	    [1]		 	  [2]  [3]	  [4]		  [5]		[TSSend?]	[TSSend?]	[外显子区域数量？]	     		[9]		  											[10]							 [0] [geneName]	[unk]	[unk]	[-1,]
#[0]	    [1]		 	  [2]  [3]	  [4]		  [5]		[6]			   [7]		[8]	     		[9]		  											[10]							 			   [11] [12]	[13]	[14]	[15]

close BB;
}elsif(-e "../project/Refgene/genes_exon_".$spevers[$i].".txt"){
@genepositions=();
%hashgenepositionTOstrand=();
%hashgenepositionTOname=();
%hashgenepositionTOexonStart=();
%hashgenepositionTOexonEnd=();
open BB,"../project/Refgene/genes_exon_".$spevers[$i].".txt";
#Gene stable ID	Transcript stable ID	Exon region start (bp)	Exon region end (bp)	Chromosome/scaffold name	Gene start (bp)	Gene end (bp)	Gene name	Strand	Transcript start (bp)	Transcript end (bp)
#ENSRNA049474694	ENSRNA049474694-T1	12120216	12120334	11	12120216	12120334	5S_rRNA	1	12120216	12120334

while(<BB>){
s/\r|\n//g;
@tmp=split/\t/,$_;
$tmp[4]="chr".$tmp[4];
if($tmp[8] eq "1"){
$tmp[8]="+";
}elsif($tmp[8] eq "-1"){
$tmp[8]="-";
}
	if(!($tmp[0]=~ /Gene/i) && !exists $hashgenepositionTOstrand{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}){
	$hashgenepositionTOstrand{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}=$tmp[8];
	$hashgenepositionTOname{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}=$tmp[7];
	$hashgenepositionTOexonStart{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}=$tmp[2];
	$hashgenepositionTOexonEnd{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}=$tmp[3];
	push @genepositions,$tmp[4]."\t".$tmp[5]."\t".$tmp[6];
	}elsif(!($tmp[0]=~ /Gene/i) && exists $hashgenepositionTOstrand{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}){
	$hashgenepositionTOexonStart{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}.=",".$tmp[2];
	$hashgenepositionTOexonEnd{$tmp[4]."\t".$tmp[5]."\t".$tmp[6]}.=",".$tmp[3];
	}
}
close BB;
	foreach $geneposition (@genepositions){
	print CC $geneposition."\t".$hashgenepositionTOstrand{$geneposition}."\t".$hashgenepositionTOname{$geneposition}."\t".$hashgenepositionTOexonStart{$geneposition}."\t".$hashgenepositionTOexonEnd{$geneposition}."\n";
	}

}

close CC;
system("bedtools sort -i ".$spevers[$i].".bed>".$spevers[$i]."sort.bed");
mkdir("Genome/".$spes[$i]);
%hashchr=();
@chrs=();
open CHR,$spevers[$i]."sort.bed";
while(<CHR>){
s/\r|\n//g;
@tmp=split/\t/,$_;
	if($tmp[0] ne "" && !exists $hashchr{$tmp[0]}){
	$hashchr{$tmp[0]}="";
	push @chrs,$tmp[0];
	$tmpchr=$tmp[0];
	open $tmpchr,">"."Genome/".$spes[$i]."/".$tmpchr;
	print $tmpchr $_."\n";
	}elsif($tmp[0] ne "" && exists $hashchr{$tmp[0]}){
	$tmpchr=$tmp[0];
	print $tmpchr $_."\n";
	}
}
close CHR;
	foreach $chrom (@chrs){
	close $chrom;
	}
}




###############################
#以下是需要在ensembl中下载相应的gene_exon的文件
# @speschromsize=("Kluyveromyces_lactis","Schizosaccharomyces_pombe","Plasmodium_falciparum","Ectocarpus_siliculosus","Chlorocebus_sabaeus.ChlSab1",
# "Oryza_sativa.IRGSP1","Solanum_lycopersicum.SL3","Strongylocentrotus_purpuratus.Spur5","Arabidopsis_thaliana.TAIR10","Toxoplasma_gondii.TGA4",
# "Trypanosoma_brucei.TryBru","Xenopus_tropicalis.Xtro10","Zea_mays.Zm-B73");


#########################################
#genes_exon_TAIR10 	 https://plants.ensembl.org/biomart/martview/2b9ddec0bca26c936126e569b1a1295e   "TAIR10"
#genes_exon_Zm-B73	  60 (Zm-B73-REFERENCE-NAM-5.0)
#genes_exon_SL3    60  
#genes_exon_IRGSP	60 
#genes_exon_ChlSab1  Chlorocebus_sabaeus https://asia.ensembl.org/biomart/martview/ffc1b51b784cf958b678dfc4a4731a3e  113  ChlSab1  Vervet-AGM (ChlSab1.1) 
#genes_exon_ASM276v2   https://protists.ensembl.org/biomart/martview/a6bae05a1b023d6710e336592eb59d36 60   (ASM276v2)
#genes_exon_TGA4	60
#genes_exon_TryBru	60
#genes_exon_ASM294v2    genesSchizosaccharomyces_pombe https://fungi.ensembl.org/biomart/martview/39a165747bef091a9aac5a936c77c251 60 ASM294v2
#genes_exon_Spur5    genesStrongylocentrotus_purpuratus https://metazoa.ensembl.org/biomart/martview/6eea85f8b904fbc76b3a9d7ab3bb96d2 60  Strongylocentrotus purpuratus (Purple sea urchin, Spur 01) genes (Spur_5.0)
#genes_exon_Xtro10	genesXenopus_tropicalis  https://www.ensembl.org/biomart/martview/3fd9943d70094de6d7fc0ff5f57dae94  113 Tropical clawed frog (UCB_Xtro_10.0) ▼
############################
