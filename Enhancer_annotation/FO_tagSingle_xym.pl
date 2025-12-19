# @spevers=("TAIR10","ce11","ChlSab1","danRer11","dm6","ASM31002v1","galGal6","hg38","mm10","ASM251v1","rheMac10","IRGSP1",
# "panTro6","ASM276v2","rn7","sacCer3","ASM294v2","SL3","Spur5","susScr11","TGA4","TryBru","Xtro10","Zm-B73");
# 


$start = time;
open IDID,"cellsID.txt";
#Xenopus_tropicalis_blastula_stage_9	Xenopus_tropicalis	XENOPUS_TROPICALIS001
%hashcellspeTOid=();
while(<IDID>){
s/\r|\n//g;
@tmp=split/\t/,$_;
$hashcellspeTOid{$tmp[0]."\t".$tmp[1]}=$tmp[2];
}
close IDID;

mkdir("single");
#@spes=("hs","mm","dm");
#@spes=("hs","mm","dm");
#@spes=("gg","mm","at","dr");
# @spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
#@spes=("mm","at","dr","gg","hs","mam","pt");
@species=("Arabidopsis_thaliana","Caenorhabditis_elegans","Chlorocebus_sabaeus","Danio_rerio","Drosophila_melanogaster",
 "Ectocarpus_sp._Ec32","Gallus_gallus","Homo_sapiens","Mus_musculus","Kluyveromyces_lactis","Macaca_mulatta","Oryza_sativa",
  "Pan_traglodytes","Plasmodium_falciparum","Rattus_norvegicus","Saccharomyces_cerevisiae","Schizosaccharomyces_pombe",
  "Solanum_lycopersicum","Strongylocentrotus_purpuratus","sus_scrofa","Toxoplasma_gondii","Trypanosoma_brucei_brucei",
  "Xenopus_tropicalis","Zea_mays");
@spes=("TAIR10","ce11","ChlSab1","danRer11","dm6","ASM31002v1","galGal6","hg38","mm10","ASM251v1","rheMac10","IRGSP1",
"panTro5","ASM276v2","rn7","sacCer3","ASM294v2","SL3","Spur5","susScr11","TGA4","TryBru","Xtro10","Zm-B73");

for($i=0;$i<scalar(@species);$i++){
#foreach $spe (@spes){
	mkdir("single/".$spes[$i]);
	my $dir="download/enhancer/".$species[$i];
	opendir(DIR,$dir) or "can't open the file";
	@dirs=readdir DIR;
	@dirs=grep{$_ ne "." && $_ ne ".."} @dirs;
	foreach $file (@dirs){
		$cell=$file;
		$cell=~ s/\.bed//g;
#SNP hg38
		%hashenhSNP=();
		if(-e "download/SNP/".$spes[$i]."/".$file){
		open SNPSNP,"download/SNP/".$spes[$i]."/".$file;
		#chr1	1683817	1684411	7.133975987	chr1	1683980	1683982	rs188187201
			while(<SNPSNP>){
			s/\r|\n//g;
			@tmp=split/\t/,$_;
			if(!exists $hashenhSNP{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$hashenhSNP{$tmp[0].":".$tmp[1]."-".$tmp[2]}=$tmp[7].":".$tmp[4].":".$tmp[5];
			}else{
			$hashenhSNP{$tmp[0].":".$tmp[1]."-".$tmp[2]}.=";".$tmp[7].":".$tmp[4].":".$tmp[5];
			}

			}
		close SNPSNP;
		}
#GeneHancer hg38
		%hashenhGH=();
		if(-e "download/GeneHancer/".$spes[$i]."/".$file){
		open GH,"download/GeneHancer/".$spes[$i]."/".$file;
			while(<GH>){
			s/\r|\n//g;
			@tmp=split/\t/,$_;
			if(!exists $hashenhGH{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$hashenhGH{$tmp[0].":".$tmp[1]."-".$tmp[2]}=$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}else{
			$hashenhGH{$tmp[0].":".$tmp[1]."-".$tmp[2]}.=";".$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}

			}
		close GH;
		}
#dbSUPER  hg38
		%hashenhdbSUPER=();
		if(-e "download/dbSUPER/".$spes[$i]."/".$file){
		open DBSUPER,"download/dbSUPER/".$spes[$i]."/".$file;
			while(<DBSUPER>){
			s/\r|\n//g;
			@tmp=split/\t/,$_;
			if(!exists $hashenhdbSUPER{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$hashenhdbSUPER{$tmp[0].":".$tmp[1]."-".$tmp[2]}=$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}else{
			$hashenhdbSUPER{$tmp[0].":".$tmp[1]."-".$tmp[2]}.=";".$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}

			}
		close DBSUPER;
		}
#Motif 
		%hashenhMotif=();
		if(-e "download/Motif/".$spes[$i]."/".$file){
		open MOTIF,"download/Motif/".$spes[$i]."/".$file;
			while(<MOTIF>){
			s/\r|\n//g;
			@tmp=split/\t/,$_;
			if(!exists $hashenhMotif{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$hashenhMotif{$tmp[0].":".$tmp[1]."-".$tmp[2]}=$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}else{
			$hashenhMotif{$tmp[0].":".$tmp[1]."-".$tmp[2]}.=";".$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}

			}
		close MOTIF;
		}		

#HEDD  hg38
		%hashenhHEDD=();
		if(-e "download/HEDD/".$spes[$i]."/".$file){
		open HEDD,"download/HEDD/".$spes[$i]."/".$file;
			while(<HEDD>){
			s/\r|\n//g;
			@tmp=split/\t/,$_;
			if(!exists $hashenhHEDD{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$hashenhHEDD{$tmp[0].":".$tmp[1]."-".$tmp[2]}=$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}else{
			$hashenhHEDD{$tmp[0].":".$tmp[1]."-".$tmp[2]}.=";".$tmp[7].":".$tmp[4].":".$tmp[5].":".$tmp[6];
			}

			}
		close HEDD;
		}
		 
#TargetGene  hg38
		%hashenhTargetGene=();
		if(-e "download/TargetGene/".$spes[$i]."/".$file){
		open TargetGene,"download/TargetGene/".$spes[$i]."/".$file;
			while(<TargetGene>){
			s/\r|\n//g;
			@tmp=split/\t/,$_;
			if($tmp[0] ne ""){
			$hashenhTargetGene{$tmp[0]}=$tmp[1];
			}
			}
		close TargetGene;
		}
	
		open AA,"download/enhancer/".$species[$i]."/".$file;
		open BB,">"."single/".$spes[$i]."/".$file;
		$num=1;
		while(<AA>){
		s/\r|\n//g;
		@tmp=split/\t/,$_;
		$ID=$hashcellspeTOid{$cell."\t".$spes[$i]}."-".sprintf("%05d",$num);
		$str=$ID."\t".$tmp[0].":".$tmp[1]."-".$tmp[2];
			if(exists $hashenhSNP{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$str.="\t".$hashenhSNP{$tmp[0].":".$tmp[1]."-".$tmp[2]};
			}else{
			$str.="\t"."NA";
			}
			if(exists $hashenhGH{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$str.="\t".$hashenhGH{$tmp[0].":".$tmp[1]."-".$tmp[2]};
			}else{
			$str.="\t"."NA";
			}
			if(exists $hashenhdbSUPER{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$str.="\t".$hashenhdbSUPER{$tmp[0].":".$tmp[1]."-".$tmp[2]};
			}else{
			$str.="\t"."NA";
			}
			if(exists $hashenhMotif{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$str.="\t".$hashenhMotif{$tmp[0].":".$tmp[1]."-".$tmp[2]};
			}else{
			$str.="\t"."NA";
			}
			if(exists $hashenhHEDD{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$str.="\t".$hashenhHEDD{$tmp[0].":".$tmp[1]."-".$tmp[2]};
			}else{
			$str.="\t"."NA";
			}
			if(exists $hashenhTargetGene{$tmp[0].":".$tmp[1]."-".$tmp[2]}){
			$str.="\t".$hashenhTargetGene{$tmp[0].":".$tmp[1]."-".$tmp[2]};
			}else{
			$str.="\t"."NA";
			}
			
			
		print BB $str."\n";
		$num++;
		}
		close AA;
		close BB;
	
	}
}

$duration = time - $start;
print "All are done: $duration s\n";