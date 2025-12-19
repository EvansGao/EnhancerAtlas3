use List::Util qw[min max];

# # # @spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
# # # @spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
# # # @upstreams=("5000","5000","5000","4000","5000","2000","2800","5000","5000","250","250","800","3000");
# # # @downstreams=("500","500","500","400","500","200","280","500","500","25","25","80","300");
#@spes=("mm","at","dr","gg","hs","mam","pt");
#@spevers=("mm10","TAIR10","danRer10","galGal6","hg38","rheMac10","panTro5");
#@upstreams=("5000","250","2800","2000","5000","5000","5000");
#@downstreams=("500","25","280","200","500","500","500");
@speverschroms=("Arabidopsis_thaliana.TAIR10","ce11","Chlorocebus_sabaeus.ChlSab1","danRer11","dm6",
 "galGal6","hg38","mm10","rheMac10","Oryza_sativa.IRGSP1",
 "panTro5","Plasmodium_falciparum.ASM276v2","rn7","sacCer3","Schizosaccharomyces_pombe.ASM294v2",
 "Solanum_lycopersicum.SL3","Strongylocentrotus_purpuratus.Spur5","susScr11",
 "Toxoplasma_gondii.TGA4","Trypanosoma_brucei.TryBru","Xenopus_tropicalis.Xtro10","Zea_mays.Zm-B73");

@spevers=("TAIR10","ce11","ChlSab1","danRer11","dm6","galGal6","hg38","mm10","rheMac10","IRGSP1",
  "panTro5","ASM276v2","rn7","sacCer3","ASM294v2","SL3","Spur5","susScr11","TGA4","TryBru","Xtro10","Zm-B73");


# 处理代码中的其他逻辑和修复问题

# 修复代码中的错误和完善逻辑
#"Ectocarpus_siliculosus.ASM31002v1", "Kluyveromyces_lactis.ASM251v1",
#,"ASM31002v1","ASM251v1"
#"350","35","25","0",
# @spes=("Arabidopsis_thaliana","Caenorhabditis_elegans","Chlorocebus_sabaeus","Danio_rerio","Drosophila_melanogaster",
#  "Ectocarpus_sp._Ec32","Gallus_gallus","Homo_sapiens","Mus_musculus","Kluyveromyces_lactis","Macaca_mulatta","Oryza_sativa",
#   "Pan_traglodytes","Plasmodium_falciparum","Rattus_norvegicus","Saccharomyces_cerevisiae","Schizosaccharomyces_pombe",
#   "Solanum_lycopersicum","Strongylocentrotus_purpuratus","sus_scrofa","Toxoplasma_gondii","Trypanosoma_brucei_brucei",
#   "Xenopus_tropicalis","Zea_mays");

@upstreams=("200","200","5000","2800","250","1800","5000","5000","5000","600","5000","50","4200","25","25","1300","1500","4000","150","180","2300","3500");
@downstreams=("20","20","500","280","25","180","500","500","500","60","500","0","420","0","0","130","150","400","0","0","230","350");#得补充#######################


#Oryza_sativa.IRGSP1.chrom.sizes
for($i=0;$i<scalar(@spevers);$i++){
open CHROM,"../project/chromsizes/".$speverschroms[$i].".chrom.sizes";
%hashchr=();
while(<CHROM>){
s/\r|\n//g;
@tmp=split/\t/,$_;
if(length($tmp[0])<8){
$hashchr{$tmp[0]}="";
}
}
close CHROM;

if(-e "../project/Refgene/Refgene_".$spevers[$i].".txt"){
open AA,"../project/Refgene/Refgene_".$spevers[$i].".txt";
%hashpro=();
open PRO,">promoterpre".$spevers[$i].".bed";
while(<AA>){
	s/\r|\n//g;
	@tmp=split/\t/,$_;
	$tmp[9]=~ s/\,$//g;
	$tmp[10]=~ s/\,$//g;
	$tmp[12]=~ s/^\s+|\s+$//g;
	if(exists $hashchr{$tmp[2]} && $tmp[3] eq "+" && ($tmp[4]-$upstreams[$i])>0 && $tmp[12] ne ""){
				if(!exists $hashpro{$tmp[2]."\t".($tmp[4]-$upstreams[$i])."\t".($tmp[4]+$downstreams[$i])."\t".$tmp[12]}){
				$hashpro{$tmp[2]."\t".($tmp[4]-$upstreams[$i])."\t".($tmp[4]+$downstreams[$i])."\t".$tmp[12]}="";
				print PRO $tmp[2]."\t".($tmp[4]-$upstreams[$i])."\t".($tmp[4]+$downstreams[$i])."\t".$tmp[12]."\n";
				}
		}elsif(exists $hashchr{$tmp[2]} && $tmp[3] eq "-" && ($tmp[5]-$downstreams[$i])>0 && $tmp[12] ne ""){
				if(!exists $hashpro{$tmp[2]."\t".($tmp[5]-$downstreams[$i])."\t".($tmp[5]+$upstreams[$i])."\t".$tmp[12]}){
				$hashpro{$tmp[2]."\t".($tmp[5]-$downstreams[$i])."\t".($tmp[5]+$upstreams[$i])."\t".$tmp[12]}="";
				print PRO $tmp[2]."\t".($tmp[5]-$downstreams[$i])."\t".($tmp[5]+$upstreams[$i])."\t".$tmp[12]."\n";
				}
		}
	}
	close AA;
	close PRO;
}elsif(-e "../project/Refgene/genes_exon_".$spevers[$i].".txt"){
	open AA,"../project/Refgene/genes_exon_".$spevers[$i].".txt";    #genes_exon_ASM294v2.txt
	%hashpro=();
	open PRO,">promoterpre".$spevers[$i].".bed";
	while(<AA>){
	s/\r|\n//g;
	s/\r|\n//g;
	@tmp=split/\t/,$_;
	$tmp[2]=~ s/^\s+|\s+$//g;
	$tmp[3]=~ s/^\s+|\s+$//g;
	$tmp[4]=~ s/^\s+|\s+$//g;
	#$tmp[4]="chr".$tmp[4];
	$tmp[5]=~ s/^\s+|\s+$//g;
	$tmp[6]=~ s/^\s+|\s+$//g;
	$tmp[7]=~ s/^\s+|\s+$//g;
	if (!($spevers[$i] =~ /ASM31002v1|TGA4|Spur5/i)) {
		$tmp[4]="chr".$tmp[4];
		}
	if(exists $hashchr{$tmp[4]} && $tmp[8] eq "1" && ($tmp[5]-$upstreams[$i])>0  && $tmp[7] ne ""){
			if(!exists $hashpro{$tmp[4]."\t".($tmp[5]-$upstreams[$i])."\t".($tmp[5]+$downstreams[$i])."\t".$standard}){
			$hashpro{$tmp[4]."\t".($tmp[5]-$upstreams[$i])."\t".($tmp[5]+$downstreams[$i])."\t".$tmp[7]}="";
			print PRO $tmp[4]."\t".($tmp[5]-$upstreams[$i])."\t".($tmp[5]+$downstreams[$i])."\t".$tmp[7]."\n";
			}
	}elsif(exists $hashchr{$tmp[4]} && $tmp[8] eq "-1" && ($tmp[6]-$downstreams[$i])>0 && $tmp[7] ne ""){
			if(!exists $hashpro{$tmp[4]."\t".($tmp[6]-$downstreams[$i])."\t".($tmp[6]+$upstreams[$i])."\t".$standard}){
			$hashpro{$tmp[4]."\t".($tmp[6]-$downstreams[$i])."\t".($tmp[6]+$upstreams[$i])."\t".$tmp[7]}="";
			print PRO $tmp[4]."\t".($tmp[6]-$downstreams[$i])."\t".($tmp[6]+$upstreams[$i])."\t".$tmp[7]."\n";
			}
	}
	}
	close AA;
	close PRO;
}

system("bedtools sort -i promoterpre".$spevers[$i].".bed>promotersort".$spevers[$i].".bed");
open CC,"promotersort".$spevers[$i].".bed";
open DD,">promoter".$spevers[$i].".bed";
@peaks=();
%hashpeak=();
while(<CC>){
chomp($_);
@temp=split/\t/,$_;
$temp[3]=~ s/\(.*//g;
if($temp[0] ne "" && !exists $hashpeak{$temp[0]."\t".$temp[1]."\t".$temp[2]}){
push @peaks,$temp[0]."\t".$temp[1]."\t".$temp[2];
$hashpeak{$temp[0]."\t".$temp[1]."\t".$temp[2]}=$temp[3];
}elsif(exists $hashpeak{$temp[0]."\t".$temp[1]."\t".$temp[2]} && !($hashpeak{$temp[0]."\t".$temp[1]."\t".$temp[2]}=~ /$temp[3]/)){
$hashpeak{$temp[0]."\t".$temp[1]."\t".$temp[2]}.=",".$temp[3];
}

}
close CC;
foreach $peak (@peaks){
@names=split/\,/,$hashpeak{$peak};
@names=sort { lc($a) cmp lc($b) } @names;
$hashpeak{$peak}=join(",",@names);
print DD $peak."\t".$hashpeak{$peak}."\n";
}
close DD;
unlink("promoterpre".$spevers[$i].".bed");
unlink("promotersort".$spevers[$i].".bed");

}


#############################

# Homo sapiens	3209286105	5000	5000/500
# Mus musculus	2730871774	4254.640572	5000/500
# Caenorhabditis_elegans	100286401	156.2440956	200/20
# Danio_rerio	1679203469	2616.163555	2800/280
# Drosophila_melanogaster	143726002	223.9220769	250/25
# Gallus_gallus	1065365425	1659.816841	1800/180
# Rattus_norvegicus	2647915728	4125.396804	4200/420
# Saccharomyces_cerevisiae	12157105	18.94051294	25
# Sus_scrofa	2501912388	3897.926682	4000/400
# Ambystoma_mexicanum	28206906136	43945.76428	44000/4400
# Arabidopsis_thaliana	119667750	186.4398282	200/20
# Bos_taurus	2715853792	4231.242873	4300/430
# Chlorocebus_sabaeus	2789656328	4346.225666	5000/500
# Cricetulus_griseus	2358167390	3673.97501	4000/400
# Ectocarpus_sp._Ec32	195810619	305.0688106	350/35
# Zea_mays	2182075994	3399.628333	3500/350
# Kluyveromyces_lactis	10689156	16.65347939	25
# Macaca_mulatta	2971331530	4629.271796	5000/500
# Myotis_lucifugus	2034575300	3169.825365	3500/350
# Oryza_sativa	375049285	584.3188683	600/60
# Oryza_sativa_Japonica_Group	375049285	584.3188683	600/60
# Pan_troglodytes	3050398082	4752.455814	5000/500
# Plasmodium_falciparum	23292622	36.28941334	50
# Saccharomyces_cerevisiae_BY4741	12157105	18.94051294	25
# Schizosaccharomyces_pombe	12631379	19.67942182	25
# Solanum_lycopersicum	827747456	1289.613062	1300/130
# Strongylocentrotus_purpuratus	921855793	1436.231864	1500/150
# Tetrahymena_thermophila	103014375	160.4942215	50
# Toxoplasma_gondii	65669794	102.3121527	150
# Trypanosoma_brucei_brucei	25789186	40.17900735	180
# Xenopus_tropicalis	1451301209	2261.096645	2300/230


##############################
# #######################
# Ambystoma_mexicanum.AmbMex60DD.chrom.sizes
# Arabidopsis_thaliana.TAIR10.chrom.sizes
# bosTau9.chrom.sizes
# ce11.chrom.sizes
# Chlorocebus_sabaeus.ChlSab1.chrom.sizes
# Cricetulus_griseus.HDv1.chrom.sizes
# danRer11.chrom.sizes
# dm6.chrom.sizes
# Ectocarpus_siliculosus.ASM31002v1.chrom.sizes
# galGal6.chrom.sizes
# hg38.chrom.sizes
# Kluyveromyces_lactis.ASM251v1.chrom.sizes
# mm10.chrom.sizes
# Oryza_sativa.IRGSP1.chrom.sizes
# panTro5.chrom.sizes
# Plasmodium_falciparum.ASM276v2.chrom.sizes
# rheMac10.chrom.sizes
# rn7.chrom.sizes
# sacCer3.chrom.sizes
# Schizosaccharomyces_pombe.ASM294v2.chrom.sizes
# Solanum_lycopersicum.SL3.chrom.sizes
# Strongylocentrotus_purpuratus.Spur5.chrom.sizes
# susScr11.chrom.sizes
# Tetrahymena_thermophila.JCVI-TTA1.chrom.sizes
# Toxoplasma_gondii.TGA4.chrom.sizes
# Trypanosoma_brucei.TryBru.chrom.sizes
# Xenopus_tropicalis.Xtro10.chrom.sizes
# Zea_mays.Zm-B73.chrom.sizes

# ##########################