#Motif/Motif_".$spes[$i].".bed  文件如何得到的？
#Motif_at.bed
#Motif_dm.bed
#Motif_dr.bed
#Motif_hs.bed
#Motif_mm.bed

$start = time;
mkdir("download/Motif");
@spevers=("hg38","mm10","TAIR10","danRer11","dm6");  
@spes=("hs","mm","at","dr","dm");
@species=("Homo_sapiens","Mus_musculus","Arabidopsis_thaliana","Danio_rerio","Drosophila_melanogaster");
for($i=0;$i<scalar(@spes);$i++){
	#foreach $spe (@spes){
	mkdir("download/Motif/".$spevers[$i]);
	my $dir="download/enhancer/".$species[$i];
	opendir(DIR,$dir) or "can't open the file";
	@dir=readdir DIR;
	@dir=grep{$_ ne "." && $_ ne ".."} @dir;
		foreach $cellfile (@dir){
		$cellname=$cellfile;
		$cellname=~ s/\.bed$//g;
			if(!(-e "download/Motif/".$spevers[$i]."/".$cellfile)){
			system("bedtools intersect -a ".$dir."/".$cellfile." -b Motif/Motif_".$spes[$i].".bed -F 1.0 -wa -wb>download/Motif/".$spevers[$i]."/".$cellfile);
			}
		}
	}

$duration = time - $start;
print "All are done: $duration s\n";

#scp -r /data/gts/dbscATAC/DataProcess/Motif datacenter@192.168.36.72:/data1/xym
#rsync -avz -e ssh /data/gts/dbscATAC/DataProcess/Motif datacenter@192.168.36.72:/data1/xym/motifxym

#以下操作可行
#sudo chown -R xianym:xianym /data/gts/dbscATAC/DataProcess   #需在120上获得相关权限
#rsync -avz xianym@120.79.22.230:/data/gts/dbscATAC/DataProcess/Motif  /data1/xym/motifxym
#运行rsync在本地/data1/xym/motifxym运行，输入120服务器的密码
############################
# @spevers=("TAIR10","ce11","ChlSab1","danRer11","dm6","ASM31002v1","galGal6","hg38","mm10","ASM251v1","rheMac10","IRGSP1",
# "panTro6","ASM276v2","rn7","sacCer3","ASM294v2","SL3","Spur5","susScr11","TGA4","TryBru","Xtro10","Zm-B73");
# @spes=("Arabidopsis_thaliana","Caenorhabditis_elegans","Chlorocebus_sabaeus","Danio_rerio","Drosophila_melanogaster",
# "Ectocarpus_sp._Ec32","Gallus_gallus","Homo_sapiens","Mus_musculus","Kluyveromyces_lactis","Macaca_mulatta","Oryza_sativa",
#  "Pan_traglodytes","Plasmodium_falciparum","Rattus_norvegicus","Saccharomyces_cerevisiae","Schizosaccharomyces_pombe",
#  "Solanum_lycopersicum","Strongylocentrotus_purpuratus","sus_scrofa","Toxoplasma_gondii","Trypanosoma_brucei_brucei",
#  "Xenopus_tropicalis","Zea_mays");

#@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
# # # @spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");

#liftOver Motif_dr10.bed danRer10ToDanRer11.over.chain.gz Motif_dr.bed unMapped