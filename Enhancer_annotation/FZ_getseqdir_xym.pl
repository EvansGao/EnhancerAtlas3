$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );

system("rm -R seq");
mkdir("seq");
############
# # # @builds=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
# # # @spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
# # # @spechrs=("chr1:chr2:chr3:chr4:chr5:chr6:chr7:chrX:chr10:chr11:chr8:chr9:chr12:chr13:chr14:chr15:chr16:chr17:chr21:chr19:chr22:chr18:chr20:chrY","chr1:chr2:chrX:chr3:chr4:chr5:chr6:chr7:chr10:chr8:chr14:chr9:chr11:chr13:chr12:chr15:chr16:chr17:chrY:chr18:chr19","chrMT:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr1:chr20:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chrX","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr20:chr21:chr22:chr23:chr24:chr25:chr26:chr27:chr28:chr29:chrMT:chrX:chrY","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chrX:chr8:chr9:chr11:chr10:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr20:chr19:chrY:chr22:chr21:chrM","chr1:chr2:chr3:chr4:chrZ:chr5:chr7:chr6:chr8:chr9:chr10:chr12:chr11:chr13:chr14:chr20:chr15:chr18:chr17:chr19:chr27:chr33:chr21:chrW:chr24:chr31:chr23:chr26:chr22:chr28:chr25:chr16:chr30:chr32:chrM","chr4:chr7:chr5:chr3:chr6:chr2:chr1:chr9:chr16:chr20:chr8:chr17:chr14:chr13:chr18:chr12:chr19:chr15:chr23:chr21:chr10:chr11:chr24:chr22:chr25:chrM","chr1:chr3:chr4:chr5:chr6:chr7:chrX:chr8:chr12:chr10:chr11:chr2B:chr9:chr2A:chr13:chr14:chr15:chr17:chr16:chr18:chr20:chr19:chr22:chr21:chrY","chr1:chr2:chr5:chr3:chr6:chr4:chr7:chrX:chr8:chr9:chr11:chr12:chr14:chr15:chr13:chr10:chr17:chr16:chr20:chr18:chr19:chrY","chr3R:chr3L:chr2R:chrX:chr2L:chrY:chr4:chrM","chr1:chr2:chr3:chr4:chr5:chrMt:chrPt","chr10:chr11:chr12:chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chrMt:chrPt","chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chrMt:chrPt");
# # # #@spenames=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
# # # @linkkeys=("https://ftp.ensembl.org/pub/release-80/fasta/callithrix_jacchus/dna/Callithrix_jacchus.C_jacchus3.2.1.dna.chromosome",
# # # "https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome",
# # # "https://ftp.ensembl.org/pub/release-112/fasta/macaca_fascicularis/dna/Macaca_fascicularis.Macaca_fascicularis_6.0.dna.primary_assembly",
# # # "https://ftp.ensembl.org/pub/release-112/fasta/chlorocebus_sabaeus/dna/Chlorocebus_sabaeus.ChlSab1.1.dna.chromosome",
# # # "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome",
# # # "https://ftp.ensembl.org/pub/release-106/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.chromosome",
# # # "https://ftp.ensembl.org/pub/release-80/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome",
# # # "https://ftp.ensembl.org/pub/release-112/fasta/pan_troglodytes/dna/Pan_troglodytes.Pan_tro_3.0.dna.chromosome",
# # # "https://ftp.ensembl.org/pub/release-112/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly",
# # # "https://ftp.ensembl.org/pub/release-111/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly",
# # # "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome",
# # # "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.chromosome",
# # # "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-49/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.chromosome");
# # # ################
#https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-49/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.chromosome.1.fa.gz
#https://ftp.ensembl.org/pub/release-102/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.1.fa.gz
#https://ftp.ensembl.org/pub/release-111/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.2L.fa.gz
#https://ftp.ensembl.org/pub/release-112/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly.1.fa.gz
#chimphttp://useast.ensembl.org/biomart/martview/54c798e06daa89b37378af92404d7adf 112 Pan_tro_3.0
#zebrafishhttp://may2015.archive.ensembl.org/biomart/martview/84ea87f76305b6cae6ab8ca3f1c69826 GRCz10/danRer10
#chickenhttps://apr2022.archive.ensembl.org/biomart/martview/58809e781a3afb67a2db9bf606be1aad  GRCg6a/galGal6
#cynomolgushttp://useast.ensembl.org/biomart/martview/54c798e06daa89b37378af92404d7adf 112
#maizehttps://nov2020-plants.ensembl.org/biomart/martview/906720c0d1dfd493210d4f6750047291
#flyhttp://jan2024.archive.ensembl.org/biomart/martview/087e12dede9a9f9f681702b754afcae5
#marmosethttp://may2015.archive.ensembl.org/biomart/martview/84ea87f76305b6cae6ab8ca3f1c69826
#rhesushttps://useast.ensembl.org/biomart/martview/45a1148f4a7082997e275fa16019ece2  112   Mmul_10
#ricehttps://plants.ensembl.org/biomart/martview/8cfcdb4ff2f653e19c58e6ddbb8f97b4   59
#vervethttps://useast.ensembl.org/biomart/martview/d908531fed14c15ccf599d64c3dbd042 112
#@spes=("human","chicken","zebrafish","arabidopsis","mouse");
#@fullnames=("Callithrix_jacchus","Mus_musculus","Macaca_fascicularis","Chlorocebus_sabaeus","Homo_sapiens","Gallus_gallus","Danio_rerio","Pan_troglodytes","Macaca_mulatta","Drosophila_melanogaster","Arabidopsis_thaliana","Oryza_sativa_Japonica","Zea_mays");

#@builds=("TAIR10");
#@spes=("at");
#@spechrs=("chr1:chr2:chr3:chr4:chr5:chrMt:chrPt");
#@linkkeys=("https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome");
@builds=("TAIR10","ce11","ChlSab1","danRer11","dm6","ASM31002v1","galGal6","hg38","mm10","ASM251v1","rheMac10","IRGSP1",
 "panTro5","ASM276v2","rn7","sacCer3","ASM294v2","SL3","Spur5","susScr11","TGA4","TryBru","Xtro10","Zm-B73");

@spes=("Arabidopsis_thaliana","Caenorhabditis_elegans","Chlorocebus_sabaeus","Danio_rerio","Drosophila_melanogaster",
  "Ectocarpus_sp._Ec32","Gallus_gallus","Homo_sapiens","Mus_musculus","Kluyveromyces_lactis","Macaca_mulatta","Oryza_sativa",
   "Pan_traglodytes","Plasmodium_falciparum","Rattus_norvegicus","Saccharomyces_cerevisiae","Schizosaccharomyces_pombe",
   "Solanum_lycopersicum","Strongylocentrotus_purpuratus","sus_scrofa","Toxoplasma_gondii","Trypanosoma_brucei_brucei",
   "Xenopus_tropicalis","Zea_mays");

@spechrs=("chr1:chr2:chr3:chr4:chr5:chrMt:chrPt",  #
"chrV:chrX:chrIV:chrII:chrI:chrIII:chrMtDNA", #
"chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr20:chr21:chr22:chr23:chr24:chr25:chr26:chr27:chr28:chr29:chrX:chrY:chrMT",	
"chr4:chr7:chr5:chr3:chr6:chr2:chr1:chr9:chr16:chr20:chr8:chr17:chr14:chr13:chr18:chr12:chr19:chr15:chr23:chr21:chr11:chr10:chr24:chr22:chr25:chrMT",
"chr3R:chr3L:chr2R:chrX:chr2L:chrY:chr4",
"LG01:LG02:LG03:LG04:LG05:LG06:LG07:LG08:LG09:LG10:LG11:LG12:LG13:LG14:LG15:LG16:LG17:LG18:LG19:LG20:LG21:LG22:LG23:LG24:LG25:LG26:LG27:LG28:LG29:LG30:LG31:LG32:LG33:LG34",
"chr1:chr2:chr3:chr4:chrZ:chr5:chr7:chr6:chr8:chr9:chr10:chr12:chr11:chr13:chr14:chr20:chr15:chr18:chr17:chr19:chr27:chr33:chr21:chrW:chr24:chr31:chr23:chr26:chr22:chr28:chr25:chr16:chr30:chr32:chrZ:chrMT", #
"chr1:chr2:chr3:chr4:chr5:chr6:chr7:chrX:chr8:chr9:chr11:chr10:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr20:chr19:chrY:chr22:chr21:chrMT", #
"chr1:chr2:chrX:chr3:chr4:chr5:chr6:chr7:chr10:chr8:chr14:chr9:chr11:chr13:chr12:chr15:chr16:chr17:chrY:chr18:chr19:chrMT", #
"chrA:chrB:chrC:chrD:chrE:chrF",   #
"chr1:chr2:chr5:chr3:chr6:chr4:chr7:chrX:chr8:chr9:chr11:chr12:chr14:chr15:chr13:chr10:chr17:chr16:chr20:chr18:chr19:chrY:chrMT",  #:chrMT
"chr1:chrMt:chrPt:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12",   #chrMt:chrPt
"chr1:chr3:chr4:chr5:chr6:chr7:chrX:chr8:chr12:chr10:chr11:chr2B:chr9:chr2A:chr13:chr14:chr15:chr17:chr16:chr18:chr20:chr19:chr22:chr21:chrY:chrMT",  #:chrMT
"chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14",
"chr1:chr2:chr4:chr3:chr5:chrX:chr6:chr7:chr8:chr9:chr10:chr13:chr14:chr15:chr17:chr11:chr16:chr18:chr19:chr20:chr12:chrY:chrMT",  #:chrMT
"chrIV:chrXV:chrVII:chrXII:chrXVI:chrXIII:chrII:chrXIV:chrX:chrXI:chrV:chrVIII:chrIX:chrIII:chrVI:chrI:chrMito",  #chrMito
"chrI:chrII:chrIII:chrMT",  #
"chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12",
"",
"chr1:chr13:chr6:chr2:chr14:chr15:chr9:chr8:chr3:chr4:chrX:chr7:chr5:chr16:chr11:chr10:chr17:chr12:chr18:chrY:chrMT",
"TGME49_chrIa:TGME49_chrIb:TGME49_chrII:TGME49_chrIII:TGME49_chrIV:TGME49_chrV:TGME49_chrVI:TGME49_chrVIIa:TGME49_chrVIIb:TGME49_chrVIII:TGME49_chrIX:TGME49_chrX:TGME49_chrXI:TGME49_chrXII",
"chr4:chr1:chr11:chr9:chr3:chr5:chr10:chr2:chr7:chr6:chr8",
"chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chrMT",
"chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:Mt:Pt");#得核对下网上得chr和fa得chrom是否一致


@linkkeys=("https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome",
"https://ftp.ensembl.org/pub/release-102/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.chromosome",
"https://ftp.ensembl.org/pub/release-112/fasta/chlorocebus_sabaeus/dna/Chlorocebus_sabaeus.ChlSab1.1.dna.chromosome",
"https://ftp.ensembl.org/pub/release-113/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.chromosome",
"https://ftp.ensembl.org/pub/release-111/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.primary_assembly",
"http://ftp.ensemblgenomes.org/pub/protists/release-60/fasta/protists_stramenopiles1_collection/ectocarpus_siliculosus_gca_000310025/dna/Ectocarpus_siliculosus_gca_000310025.ASM31002v1.dna.chromosome",
"https://ftp.ensembl.org/pub/release-106/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.chromosome",
"https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome",
"https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome",
"http://ftp.ensemblgenomes.org/pub/fungi/release-60/fasta/fungi_ascomycota1_collection/kluyveromyces_lactis_gca_000002515/dna/Kluyveromyces_lactis_gca_000002515.ASM251v1.dna.chromosome",
"https://ftp.ensembl.org/pub/release-112/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly",
"https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.chromosome",
"https://ftp.ensembl.org/pub/release-112/fasta/pan_troglodytes/dna/Pan_troglodytes.Pan_tro_3.0.dna.chromosome",
"http://ftp.ensemblgenomes.org/pub/protists/release-60/fasta/plasmodium_falciparum/dna/Plasmodium_falciparum.ASM276v2.dna.chromosome",
"https://ftp.ensembl.org/pub/release-113/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.primary_assembly",
"https://ftp.ensembl.org/pub/release-113/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome",
"http://ftp.ensemblgenomes.org/pub/fungi/release-60/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.dna.chromosome",
"https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/solanum_lycopersicum/dna/Solanum_lycopersicum.SL3.0.dna.chromosome",
#"http://ftp.ensemblgenomes.org/pub/metazoa/release-60/fasta/strongylocentrotus_purpuratus/dna/Strongylocentrotus_purpuratus.Spur_5.0.dna.primary_assembly.NC_001453",
"http://ftp.ensemblgenomes.org/pub/metazoa/release-60/fasta/strongylocentrotus_purpuratus/dna/Strongylocentrotus_purpuratus.Spur_5.0.dna.toplevel",
"https://ftp.ensembl.org/pub/release-113/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.primary_assembly",
"http://ftp.ensemblgenomes.org/pub/protists/release-60/fasta/toxoplasma_gondii/dna/Toxoplasma_gondii.TGA4.dna.chromosome",
"http://ftp.ensemblgenomes.org/pub/protists/release-60/fasta/trypanosoma_brucei/dna/Trypanosoma_brucei.TryBru_Apr2005_chr11.dna.chromosome",
"https://ftp.ensembl.org/pub/release-113/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.UCB_Xtro_10.0.dna.primary_assembly",
"https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-49/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.chromosome");

#有大小写得怎么办？

for($i=0;$i<scalar(@builds);$i++){
mkdir("seq/".$spes[$i]);
@chrarray=split/\:/,$spechrs[$i];
$piece=8;
my @processes=();
my $pronum=0;
for($ii=0;$ii<scalar(@chrarray);$ii++){
   $processes[$pronum]=fork();
   if($processes[$pronum]){   
#   print $ii."\t".$celllines[$ii]."\n";
   }else{ 
   Seqindex($chrarray[$ii],$spes[$i],$linkkeys[$i]);         # child handles 
    exit 0; 
  }
  if(($pronum+1)%$piece==0){
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*(($pronum+1)/$piece-1);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }elsif(($pronum+1)%$piece!=0 && ($pronum+1)==scalar(@chrarray)){
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	waitpid($processes[$k],0);
  	}
  	for($k=$piece*int(($pronum+1)/$piece);$k<($pronum+1);$k++){
   	undef($processes[$k]);
  	}
  }
	$pronum++;
}


}


$duration = time - $start;
print "All are done: $duration s\n";


sub  Seqindex(){
my ($chrom, $spe, $linkkey)=@_;
$chromNum=$chrom;
if($chromNum =~ /chr/i){
	$chromNum=~s/^chr//g;
}
#$wholelink=$linkkey.$chrom.".fa.gz";
#$tmpname="";
#if($wholelink=~ /dna\/(.*)$/){
#$tmpname=$1;
#}
#https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
#ftp://ftp.ensembl.org/pub/release-80/fasta/callithrix_jacchus/dna/Callithrix_jacchus.C_jacchus3.2.1.dna.chromosome.1.fa.gz
#print "wget '".$linkkey.".".$chromNum.".fa.gz' -q -o /dev/null -O seq/".$spe."/".$chrom.".fa.gz\n";
if(!(-e "seq/".$spe."/".$chrom.".fa")){
system("wget '".$linkkey.".".$chromNum.".fa.gz' -q -o /dev/null -O seq/".$spe."/".$chrom.".fa.gz");
system("gunzip seq/".$spe."/".$chrom.".fa.gz");
system("sed -i '1d' seq/".$spe."/".$chrom.".fa");
system("sed -i '1"."s\/^\/>".$chrom."\\n\/' seq/".$spe."/".$chrom.".fa");
if(-e "seq/".$spe."/".$chrom.".fa"){
system("samtools faidx seq/".$spe."/".$chrom.".fa");
}
}

}

