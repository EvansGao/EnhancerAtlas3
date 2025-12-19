#这是将enhancer2.0和3.0的enh和EP数据复制到download目录下
#可行
# # # 表示源程序
# 表示我的注释
#D:\xianym\代码\爬虫\enhancer2.0数据分析\网页数据格式转换\dbscATAC代码
# @allspevers=("TAIR10","ce11","ChlSab1","danRer11","dm6","ASM31002v1","galGal6","hg38","ASM251v1","rheMac10","mm10","IRGSP1",
# "panTro6","ASM276v2","rn7","sacCer3","ASM294v2","SL3","Spur5","susScr11","TGA4","TryBru","Xtro10","Zm-B73");

#这个pl存在目的是创建目录和将enhancer和interaction文件拷贝到自己想要转化和分析的文件目录下
#enhancer和interaction将文件重命名为：cell_line.spever.enh.bed \cell_line.spever.interaction.txt
#需要将物种区分开来复制，其中有七个物种："hs38","mm10","galGal6","rheMac10","dm6","danRer11","susScr11",是celllineNOCTCF.bed  ——hs7250NOCTCF.bed
#"TAIR10","ce11","ChlSab1","ASM31002v1","ASM251v1","IRGSP1","panTro6","ASM276v2","rn7","sacCer3","ASM294v2","SL3","Spur5","TGA4","TryBru","Xtro10","Zm-B73",是cellline_combined.bed   ——hs7250_combined.bed
$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
use 5.010;
#@spevers=("mm10","TAIR10","danRer10","galGal5","hg38","rheMac10","panTro5");
#@spes=("mm","at","dr","gg","hs","mam","pt");
# # # @spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
# # # @spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");

# my $dir = "/data1/xym/enhanceratlas3.0/";   
# opendir(DIR, $dir) or die "can't open the file";
# my @dir = readdir DIR;
# my @cells = grep { $_ ne "." && $_ ne ".." && $_ ne "CTCFfastq" && $_ ne "OTHERS" && $_ ne "copy.txt" && $_ ne "nohup.out"} @dir;
# print "@cells\n";
# foreach $cell (@cells){
#     print "$cell\n";   #可以这里建一个hash，这样
	# $dirpath1 = $dir."/".$cell."/".$cell;   
	# opendir(DIR, $dirpath1) or die "can't open the file";
	# @dirpath1 = readdir DIR;
	# @projects = grep { $_ ne "." && $_ ne ".."} @dir;
	# print "@projects\n";   #这是所有物种的细胞系都会读取到
# }
# # # @ids=("GEO:GSE204761","GEO:GSE196794","GEO:GSE199739","GEO:GSE228270","GEO:GSE214082","GEO:GSE192772","GEO:GSE179705","GEO:GSE151230","GEO:GSE111586","GEO:GSE131688","GEO:GSE155178","GEO:GSE163697","GEO:GSE149683","GEO:GSE214132","10X:10X_MouseBrain","GEO:GSE173834");
#代码是读取所有物种的所有细胞系或组织

#!(chipseq)/data1/xym/DataProcess$
#this pl can run in both /data1/xym/ and /data1/xym/DataProcess

system("rm -R download");
mkdir("download");
mkdir("download/enhancer");
mkdir("download/interaction");
# @spevers=("TAIR10","ce","ChlSab1","danRer","dm","ASM31002v1","GalGal","hs","ASM251v1","rheMac","mm","IRGSP1",
# "panTro","ASM276v2","rn","sacCer","ASM294v2","SL3","Spur5","susScr","TGA4","TryBru","Xtro10","Zm-B73");

@spevers=("TAIR10","ce11","ChlSab1","danRer11","dm6","ASM31002v1","galGal6","hg38","mm10","ASM251v1","rheMac10","IRGSP1",
"panTro5","ASM276v2","rn7","sacCer3","ASM294v2","SL3","Spur5","susScr11","TGA4","TryBru","Xtro10","Zm-B73");
@spes=("Arabidopsis_thaliana","Caenorhabditis_elegans","Chlorocebus_sabaeus","Danio_rerio","Drosophila_melanogaster",
"Ectocarpus_sp._Ec32","Gallus_gallus","Homo_sapiens","Mus_musculus","Kluyveromyces_lactis","Macaca_mulatta","Oryza_sativa",
 "Pan_traglodytes","Plasmodium_falciparum","Rattus_norvegicus","Saccharomyces_cerevisiae","Schizosaccharomyces_pombe",
 "Solanum_lycopersicum","Strongylocentrotus_purpuratus","sus_scrofa","Toxoplasma_gondii","Trypanosoma_brucei_brucei",
 "Xenopus_tropicalis","Zea_mays");
 @speversinteraction=("dm6","hg38","mm10");
 @name=("Drosophila_melanogaster_","hs","mus");
# @speversNOCTCF =("hs38","mm10","galGal6","rheMac10","dm6","danRer11","susScr11");  #,是celllineNOCTCF.bed  ——hs7250NOCTCF.bed
# @spevers_combined =("TAIR10","ce11","ChlSab1","ASM31002v1","ASM251v1","IRGSP1","panTro6","ASM276v2","rn7","sacCer3","ASM294v2","SL3","Spur5","TGA4","TryBru","Xtro10","Zm-B73"); #是cellline_combined.bed   ——hs7250_combined.bed

#Oryza_sativa_Japonica_Group 把物种归到Oryza_sativa
#Saccharomyces_cerevisiae_BY4741 把物种归到Saccharomyces_cerevisiae

$dir ="/data1/xym/";
%hashspeverTOspe=();
for($i=0;$i<scalar(@spevers);$i++){
	$hashspeverTOspe{$spevers[$i]}=$spes[$i];
	$hashspeverTOspe{$spes[$i]}=$spevers[$i];
}

%hashspeinteractionTOname=();
for($j=0;$j<scalar(@name);$j++){
	$hashspeinteractionTOname{$speversinteraction[$j]}=$name[$j];
	$hashspeinteractionTOname{$name[$j]}=$speversinteraction[$j];
}
#创建hash:$hashspeverTOspe{"TAIR10"}="Arabidopsis_thaliana";
#创建hash:$hashspeverTOspe{"Arabidopsis_thaliana"}="TAIR10";


foreach $spevers (@spevers){
	mkdir("download/enhancer/".$hashspeverTOspe{$spevers});   #download/enhancer/Arabidopsis_thaliana
	#mkdir("download/interaction/".$hashspeverTOspe{$spevers});  

	$dirpath1 = $dir."enhanceratlas3.0/".$hashspeverTOspe{$spevers}."/".$hashspeverTOspe{$spevers};   
	opendir(DIR, $dirpath1) or die "can't open the file";
	@dirpath1 = readdir DIR;
	@projects1 = grep { $_ ne "." && $_ ne ".."} @dirpath1;  #$projects1="Drosophila_melanogaster_E4-8"
	print "@projects1\n";
	foreach $projects1 (@projects1){
		#if(exists $dirpath1."/".$projects1."/".$projects1."NOCTCF.bed"){
		if(-e "$dirpath1/$projects1/${projects1}NOCTCF.bed"){
		#copy($dirpath1."/".$projects1."/".$projects1."NOCTCF.bed", "download/enhancer/".$hashspeverTOspe{$spevers}."/".$projects1.".".$spevers.".enh.bed");
		copy($dirpath1."/".$projects1."/".$projects1."NOCTCF.bed", "download/enhancer/".$hashspeverTOspe{$spevers}."/".$projects1.".bed");
	}else{
		#copy($dirpath1."/".$projects1."/".$projects1."_combined.bed", "download/enhancer/".$hashspeverTOspe{$spevers}."/".$projects1.".".$spevers.".enh.bed");
		copy($dirpath1."/".$projects1."/".$projects1."_combined.bed", "download/enhancer/".$hashspeverTOspe{$spevers}."/".$projects1.".bed");
	}
}


#enhancer2.0 data copy
$dirpath2 =$dir."enhancer2.0Allcelline/Download_enhancers";
opendir(DIR, $dirpath2) or die "can't open the file";
@dirpath2 = readdir DIR;
@species = grep { $_ ne "." && $_ ne ".."} @dirpath2;   #Caenorhabditis_elegans
foreach $species (@species){	
		$dirfile=$dirpath2."/".$species;  
		opendir(DIR, $dirfile) or die "can't open the file";
		@dirfile = readdir DIR;
		@projects2 = grep { $_ ne "." && $_ ne ".." && $_ ne "other"} @dirfile;   #$projects2="Caenorhabditis_elegans_Larvae_L1"
		foreach $projects2 (@projects2){
			#if(exists $dirpath2."/".$species."/".$projects2."/".$projects2."NOCTCF.bed"){
			if(-e "$dirpath2/$species/$projects2/${projects2}NOCTCF.bed"){
		#copy($dirpath2."/".$species."/".$projects2."/".$projects2."NOCTCF.bed", "download/enhancer/".$species."/".$projects2.".".$hashspeverTOspe{$species}.".enh.bed");
		copy($dirpath2."/".$species."/".$projects2."/".$projects2."NOCTCF.bed", "download/enhancer/".$species."/".$projects2.".bed");
	}else{
		#copy($dirpath2."/".$species."/".$projects2."/".$projects2."_combined.bed", "download/enhancer/".$species."/".$projects2.".".$hashspeverTOspe{$species}.".enh.bed");
		copy($dirpath2."/".$species."/".$projects2."/".$projects2."_combined.bed", "download/enhancer/".$species."/".$projects2.".bed");
	}
	}
}
#如法炮制，将combined文件copy
}

#复制interactionsfiles
foreach $speversinteraction (@speversinteraction){
	#mkdir("download/enhancer/".$hashspeverTOspe{$spevers});   #download/enhancer/Arabidopsis_thaliana
	mkdir("download/interaction/".$hashspeverTOspe{$speversinteraction});
}

foreach $speversinteraction (@speversinteraction){
	#3.0
	$dirpath3 =$dir."enhancer3.0EP/".$hashspeverTOspe{$speversinteraction};
	opendir(DIR, $dirpath3) or die "can't open the file";
	@dirpath3 = readdir DIR;
	@interactions = grep { $_ ne "." && $_ ne ".."} @dirpath3;   #$interactions=Drosophila_melanogaster_Regenerating_Wing_Disc_L3.interaction.txt
	foreach $interactions (@interactions){
		@interactionsinfo=split/\.interaction\.txt/, $interactions;
		$interactionfile=$interactionsinfo[0];
		#copy($dirpath3."/".$interactions, "download/interaction/".$hashspeverTOspe{$speversinteraction}."/".$interactionfile.".".$speversinteraction.".interaction.txt");
		copy($dirpath3."/".$interactions, "download/interaction/".$hashspeverTOspe{$speversinteraction}."/".$interactionfile."_interaction.txt");
	}
	#2.0   #/data1/xym/enhancer2.0Allcelline/Download_enhancer-gene_interactions
	$dirpath4 =$dir."enhancer2.0Allcelline/Download_enhancer-gene_interactions/".$hashspeverTOspe{$speversinteraction};
	opendir(DIR, $dirpath4) or die "can't open the file";
	@dirpath4 = readdir DIR;
	@interactions = grep { $_ ne "." && $_ ne ".." && $_ ne "other"} @dirpath4;   #$interactions=BG3-c2_EP.txt
	foreach $interactions (@interactions){
		@interactionsinfo=split/\_EP\.txt/, $interactions;
		$interactionfile=$hashspeinteractionTOname{$speversinteraction}.$interactionsinfo[0];
		#copy($dirpath4."/".$interactions."/".$interactionfile."EP.txt", "download/interaction/".$hashspeverTOspe{$speversinteraction}."/".$interactionfile.".".$speversinteraction.".interaction.txt");
		copy($dirpath4."/".$interactions, "download/interaction/".$hashspeverTOspe{$speversinteraction}."/".$interactionfile."_interaction.txt");
	}
}

#需要将文件名重新命名为Prediction_cellline.txt
#Prediction_Kc167_Undamaged_Wing_Disc_L3.txt
#Prediction_K562_hsWI-38.txt
#Prediction_placenta_musMull_M_N-rCr263.4.txt




