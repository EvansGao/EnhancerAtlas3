#生成geneMapdata".$spes[$i].".txt
#it needs "webgenes/" and "/data1/xym/project/uniprot" to run this pl 


use List::MoreUtils qw(uniq);
#@spes=("human","chicken","zebrafish","arabidopsis");      #Mus musculus  Homo sapiens
#@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
#@builds=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
#marmoset (calJac3),mouse (mm10),cynomolgus (macFas6),vervet (ChlSab1_1),human (hg38),chicken (galGal6), zebrafish (danRer10), chimp (panTro5), rhesus (rheMac10), fly (dm6), arabidopsis (TAIR10), rice (IRGSP1), maize (B73v4)
#human (hg38), mouse (mm10), marmoset (calJac3), cynomolgus (macFas6),vervet (ChlSab1_1), chicken (galGal6), zebrafish (danRer10), chimp (panTro5), rhesus (rheMac10), fly (dm6), arabidopsis (TAIR10), rice (IRGSP1), maize (B73v4)
# # # @spes=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
@spes=("Arabidopsis_thaliana","Caenorhabditis_elegans","Chlorocebus_sabaeus","Danio_rerio","Drosophila_melanogaster",
"Ectocarpus_sp._Ec32","Gallus_gallus","Homo_sapiens","Mus_musculus","Kluyveromyces_lactis","Macaca_mulatta","Oryza_sativa",
 "Pan_traglodytes","Plasmodium_falciparum","Rattus_norvegicus","Saccharomyces_cerevisiae","Schizosaccharomyces_pombe",
 "Solanum_lycopersicum","Strongylocentrotus_purpuratus","sus_scrofa","Toxoplasma_gondii","Trypanosoma_brucei_brucei",
 "Xenopus_tropicalis","Zea_mays");
#@spevers=("TAIR10","ce11","ChlSab1","danRer11","dm6","ASM31002v1","galGal6","hg38","ASM251v1","rheMac10","mm10","IRGSP1",
#"panTro6","ASM276v2","rn7","sacCer3","ASM294v2","SL3","Spur5","susScr11","TGA4","TryBru","Xtro10","Zm-B73");
#@spes=("human","chicken","zebrafish","arabidopsis","mouse");
#@fullnames=("Callithrix_jacchus","Mus_musculus","Macaca_fascicularis","Chlorocebus_sabaeus","Homo_sapiens","Gallus_gallus","Danio_rerio","Pan_troglodytes","Macaca_mulatta","Drosophila_melanogaster","Arabidopsis_thaliana","Oryza_sativa_Japonica","Zea_mays");
#(chipseq)/data1/xym$

for($i=0;$i<scalar(@spes);$i++){
open AA,"webgenes/genes".$spes[$i].".txt";
#Transcript stable ID	Gene stable ID	Protein stable ID	Transcript start (bp)	Transcript end (bp)	Strand	Gene name	Transcript name	Chromosome/scaffold name
#ENSCSAT00000000001	ENSCSAG00000000001		583	656	1			MT
#Transcript stable ID	Gene stable ID	Protein stable ID	Transcript start (bp)	Transcript end (bp)	Strand	Gene name	Transcript name	Chromosome/scaffold name
#FBtr0343930	FBgn0029843	FBpp0310420	6091205	6103971	1	Nep1	Nep1-RD	X
%hashIDs=();
@IDs=();
%hashtogene=();
%hashtoprotein= ();
%hashtoPosition= ();
%hashtoTranscriptName= ();
%hashtoGeneName= ();
#%hashto
while(<AA>){
s/\r|\n//g;
@temp=split/\t/,$_;
if($spes[23] eq "Zea_mays"){    #$spes[23]
$temp[23]=$temp[1];
}
#print $temp[6]."\n";
if(!($temp[0]=~ /Transcript/i) && $temp[0] ne ""){
#if($temp[0]=~ /^EN|^FB|^AT|^Zm/i){	
	if(!exists $hashIDs{$temp[0]}){
	$hashIDs{$temp[0]}=$temp[0];  #Transcript stable ID
	push @IDs,$temp[0];
	}
if(!exists $hashtogene{$temp[0]}){
$hashtogene{$temp[0]}=$temp[1];  #Gene stable ID
}elsif(exists $hashtogene{$temp[0]} && !($hashtogene{$temp[0]}=~ /$temp[1]/)){
$hashtogene{$temp[0]}.=";".$temp[1];
}

if(!exists $hashtoprotein{$temp[0]}){
$hashtoprotein{$temp[0]}=$temp[2];   #Protein stable ID
}elsif(exists $hashtoprotein{$temp[0]} && !($hashtoprotein{$temp[0]}=~ /$temp[2]/)){
$hashtoprotein{$temp[0]}.=";".$temp[2];
}


$hashtoPosition{$temp[0]}="chr".$temp[8].":".$temp[3]."-".$temp[4].":".$temp[5];  #Chromosome/scaffold name	Transcript start (bp)	Transcript end (bp)	Strand		
$hashtoTranscriptName{$temp[0]}=$temp[7];  #Transcript name
$hashtoGeneName{$temp[0]}=$temp[6];  #Gene name

#print $hashtoGeneName[$temp[0]]."\n";
}
}
close AA;
%hashtoPDB= ();
%hashtoEntrezGene= ();
%hashtoHGNC= ();

open BB,"webgenes/pdb".$spes[$i].".txt";
#Transcript stable ID	Gene stable ID	HGNC ID	NCBI gene (formerly Entrezgene) ID	PDB ID
#ENSGALT00000080481	ENSGALG00000042750		39116926	


while(<BB>){
s/\r|\n//g;
@temp=split/\t/,$_;
if(!($temp[0]=~ /Transcript/i) && $temp[0] ne ""){
	if(!exists $hashIDs{$temp[0]}){
	$hashIDs{$temp[0]}=$temp[0];  #Transcript stable ID
	push @IDs,$temp[0];
	}
	
if(!exists $hashtoHGNC{$temp[0]}){
$hashtoHGNC{$temp[0]}=$temp[2]; #HGNC ID
}elsif(exists $hashtoHGNC{$temp[0]} && !($hashtoHGNC{$temp[0]}=~ /$temp[2]/)){
$hashtoHGNC{$temp[0]}.=";".$temp[2];
}

if(!exists $hashtoEntrezGene{$temp[0]}){
$hashtoEntrezGene{$temp[0]}=$temp[3];
}elsif(exists $hashtoEntrezGene{$temp[0]} && !($hashtoEntrezGene{$temp[0]}=~ /$temp[3]/)){
$hashtoEntrezGene{$temp[0]}.=";".$temp[3];
}

if(!exists $hashtoPDB{$temp[0]}){
$hashtoPDB{$temp[0]}=$temp[4];
}elsif(exists $hashtoPDB{$temp[0]} && !($hashtoPDB{$temp[0]}=~ /$temp[4]/)){
$hashtoPDB{$temp[0]}.=";".$temp[4];
}

}
}
close BB;

%hashRefSeqmRNA=();
%hashRefSeqncRNA=();
%hashChEMBL=();
open CC,"webgenes/ref".$spes[$i].".txt";
#Transcript stable ID	Gene stable ID	RefSeq mRNA ID	RefSeq ncRNA ID	HGNC ID
#ENSGALT00000078426	ENSGALG00000041922			

while(<CC>){
s/\r|\n//g;
@temp=split/\t/,$_;
if(!($temp[0]=~ /Transcript/i) && $temp[0] ne ""){
	if(!exists $hashIDs{$temp[0]}){
	$hashIDs{$temp[0]}=$temp[0];
	push @IDs,$temp[0];
	}
	
if(!exists $hashChEMBL{$temp[0]}){
$hashChEMBL{$temp[0]}=$temp[4];
}elsif(exists $hashChEMBL{$temp[0]} && !($hashChEMBL{$temp[0]}=~ /$temp[4]/)){
$hashChEMBL{$temp[0]}.=";".$temp[4];
}

if(!exists $hashRefSeqmRNA{$temp[0]}){
$hashRefSeqmRNA{$temp[0]}=$temp[2];
}elsif(exists $hashRefSeqmRNA{$temp[0]} && !($hashRefSeqmRNA{$temp[0]}=~ /$temp[2]/)){
$hashRefSeqmRNA{$temp[0]}.=";".$temp[2];
}

if(!exists $hashRefSeqncRNA{$temp[0]}){
$hashRefSeqncRNA{$temp[0]}=$temp[3];
}elsif(exists $hashRefSeqncRNA{$temp[0]} && !($hashRefSeqncRNA{$temp[0]}=~ /$temp[3]/)){
$hashRefSeqncRNA{$temp[0]}.=";".$temp[3];
}
}
}
close CC;

%hashRefSeqProteinID= ();
%hashUCSCID=();
%hashWikiGeneID=();
open DD,"webgenes/UCSC".$spes[$i].".txt";
#Transcript stable ID	Gene stable ID	RefSeq peptide ID	INSDC protein ID	WikiGene ID
#ENSGALT00000080481	ENSGALG00000042750	YP_009558652	QBG37922	39116926

while(<DD>){
s/\r|\n//g;
@temp=split/\t/,$_;
if(!($temp[0]=~ /Transcript/i) && $temp[0] ne ""){
	if(!exists $hashIDs{$temp[0]}){
	$hashIDs{$temp[0]}=$temp[0];
	push @IDs,$temp[0];
	}
	
if(!exists $hashRefSeqProteinID{$temp[0]}){
$hashRefSeqProteinID{$temp[0]}=$temp[2];
}elsif(exists $hashRefSeqProteinID{$temp[0]} && !($hashRefSeqProteinID{$temp[0]}=~ /$temp[2]/)){
$hashRefSeqProteinID{$temp[0]}.=";".$temp[2];
}

if(!exists $UCSCID{$temp[0]}){
$UCSCID{$temp[0]}=$temp[3];
}elsif(exists $UCSCID{$temp[0]} && !($UCSCID{$temp[0]}=~ /$temp[3]/)){
$UCSCID{$temp[0]}.=";".$temp[3];
}

if(!exists $hashWikiGeneID{$temp[0]}){
$hashWikiGeneID{$temp[0]}=$temp[4];
}elsif(exists $hashWikiGeneID{$temp[0]} && !($hashWikiGeneID{$temp[0]}=~ /$temp[4]/)){
$hashWikiGeneID{$temp[0]}.=";".$temp[4];
}
}
}
close DD;

open UNIPROT,"../project/uniprot/uniprot_genes_".$spes[$i].".txt";
#Entry	Entry name	Status	Protein names	Gene names	Organism	Length
#Q00266	METK1_HUMAN	reviewed	S-adenosylmethionine synthase isoform type-1 (AdoMet synthase 1) (EC 2.5.1.6) (Methionine adenosyltransferase 1) (MAT 1) (Methionine adenosyltransferase I/III) (MAT-I/III)	MAT1A AMS1 MATA1	Homo sapiens (Human)	395
#A0A1D6E0S8	WAK17_MAIZE	reviewed	Wall-associated receptor kinase 17 (ZmWAK17) (EC 2.7.11.1) (Receptor-like protein kinase 12)	WAK17 RLK12 ZEAMMB73_Zm00001d002447	Zea mays (Maize)	856
#A0A1D6EFT8	EDSB_MAIZE	reviewed	Eudesmanediol synthase (ZmEDS) (EC 4.2.3.197) (Terpene synthase 17) (Terpene synthase 7)	EDS TPS17 GRMZM2G010356 ZEAMMB73_Zm00001d004509	Zea mays (Maize)	557

%hashuniprotTOnames=();
while(<UNIPROT>){
s/\r|\n//g;
@tmp=split/\t/,$_;
@genenames=split/\s+/,$tmp[4];
for($m=0;$m < scalar(@genenames);$m++){
$genenames[$m]=~ s/ZEAMMB73\_//g;
}
	if($tmp[0] ne ""){
#	$hashuniprotTOnames{$tmp[0]}=join(";",@genenames[1..$#genenames]);
	$hashuniprotTOnames{$tmp[0]}=join(";",@genenames);
	}
}
close UNIPROT;



%hashSwissProtID= ();
%hashTrEMBLAccession=();
%hashSwissProtAccession=();
%hashAlias=();
open EE,"webgenes/uniprot".$spes[$i].".txt";
while(<EE>){
s/\r|\n//g;
@temp=split/\t/,$_;
#print $temp[2];

if(!($temp[0]=~ /Transcript/i) && $temp[0] ne ""){
	if(!exists $hashIDs{$temp[0]}){
	$hashIDs{$temp[0]}=$temp[1];
	push @IDs,$temp[0];
	}

	if(!exists $hashAlias{$temp[0]}){
	$hashAlias{$temp[0]}="";
		if(exists $hashuniprotTOnames{$temp[2]}){
		$hashAlias{$temp[0]}=$hashuniprotTOnames{$temp[2]};
		}
		if(exists $hashuniprotTOnames{$temp[3]}){
		$hashAlias{$temp[0]}.=";".$hashuniprotTOnames{$temp[3]};
		}
		if(exists $hashuniprotTOnames{$temp[4]}){
		$hashAlias{$temp[0]}.=";".$hashuniprotTOnames{$temp[4]};
		}
	}elsif(exists $hashAlias{$temp[0]}){
		if(exists $hashuniprotTOnames{$temp[2]}){
		$hashAlias{$temp[0]}.=";".$hashuniprotTOnames{$temp[2]};
		}
		if(exists $hashuniprotTOnames{$temp[3]}){
		$hashAlias{$temp[0]}.=";".$hashuniprotTOnames{$temp[3]};
		}
		if(exists $hashuniprotTOnames{$temp[4]}){
		$hashAlias{$temp[0]}.=";".$hashuniprotTOnames{$temp[4]};
		}
	}

if(!exists $hashSwissProtID{$temp[0]}){
$hashSwissProtID{$temp[0]}=$temp[2];
}elsif(exists $hashSwissProtID{$temp[0]} && !($hashSwissProtID{$temp[0]}=~ /$temp[2]/)){
$hashSwissProtID{$temp[0]}.=";".$temp[2];
}

if(!exists $hashTrEMBLAccession{$temp[0]}){
$hashTrEMBLAccession{$temp[0]}=$temp[3];
}elsif(exists $hashTrEMBLAccession{$temp[0]} && !($hashTrEMBLAccession{$temp[0]}=~ /$temp[3]/)){
$hashTrEMBLAccession{$temp[0]}.=";".$temp[3];
}

if(!exists $hashSwissProtAccession{$temp[0]}){
$hashSwissProtAccession{$temp[0]}=$temp[4];
}elsif(exists $hashSwissProtAccession{$temp[0]} && !($hashSwissProtAccession{$temp[0]}=~ /$temp[4]/)){
$hashSwissProtAccession{$temp[0]}.=";".$temp[4];
}

}
}
close EE;

@genes=();
%hashgeneinfo=();
foreach $ID (@IDs){
if($ID ne ""){
$hashtogene{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashtoprotein{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashtoPosition{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashtoTranscriptName{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashtoGeneName{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashtoHGNC{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashChEMBL{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashRefSeqmRNA{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashRefSeqncRNA{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashtoEntrezGene{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashtoPDB{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashRefSeqProteinID{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$UCSCID{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashFlyBasenameGene{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashSwissProtID{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashAlias{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashTrEMBLAccession{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$hashSwissProtAccession{$ID}=~ s/^(\;)(\1+)|^\s+(\;)(\1+)|^\;|^\s+\;|(\;)(\1+)$|(\;)(\1+)\s+$|\;$|\;\s+$//g;
$gene=$hashtogene{$ID};

	if(!exists $hashgeneinfo{$gene}){
	$hashgeneinfo{$gene}=$ID.";".$hashtoprotein{$ID}.";".$hashtoTranscriptName{$ID}.";".$hashtoGeneName{$ID}.";".$hashtoHGNC{$ID}.";".$hashChEMBL{$ID}.";".$hashRefSeqmRNA{$ID}.";".$hashRefSeqncRNA{$ID}.";".$hashtoEntrezGene{$ID}.";".$hashtoPDB{$ID}.";".$hashRefSeqProteinID{$ID}.";".$UCSCID{$ID}.";".$hashWikiGeneID{$ID}.";".$hashSwissProtID{$ID}.";".$hashAlias{$ID}.";".$hashTrEMBLAccession{$ID}.";".$hashSwissProtAccession{$ID};
	push @genes,$gene;
	}else{
	$hashgeneinfo{$gene}.=";".$ID.";".$hashtoprotein{$ID}.";".$hashtoTranscriptName{$ID}.";".$hashtoGeneName{$ID}.";".$hashtoHGNC{$ID}.";".$hashChEMBL{$ID}.";".$hashRefSeqmRNA{$ID}.";".$hashRefSeqncRNA{$ID}.";".$hashtoEntrezGene{$ID}.";".$hashtoPDB{$ID}.";".$hashRefSeqProteinID{$ID}.";".$UCSCID{$ID}.";".$hashWikiGeneID{$ID}.";".$hashSwissProtID{$ID}.";".$hashAlias{$ID}.";".$hashTrEMBLAccession{$ID}.";".$hashSwissProtAccession{$ID};
	}
#print $hashUniProtID{$ID}."\n";

}
}


open FF,">geneMapdata".$spes[$i].".txt";
foreach $gene (@genes){
@geneinfo=split/\;/,$hashgeneinfo{$gene};
@geneinfo=uniq(@geneinfo);
@geneinfo=grep {$_ ne ""} @geneinfo;
	if($gene ne ""){
	print FF $gene."\t".join(",",@geneinfo)."\n";
	}
}
close FF;

}





#genesDanio_rerio http://asia.ensembl.org/biomart/martview/d0b53f53fe7f2a87cc75d5b7ae833c4d  113 GRcz11 
#pdbDanio_rerio
#refDanio_rerio
#UCSCDanio_rerio
#uniprotDanio_rerio    
#uniprot_genes_Danio_rerio  
#genesCaenorhabditis_elegans https://asia.ensembl.org/biomart/martview/ffc1b51b784cf958b678dfc4a4731a3e  113  WBcel235/ce11
#pdbCaenorhabditis_elegans    #Expression Atlas ID
#refCaenorhabditis_elegans
#UCSCCaenorhabditis_elegans
#uniprotCaenorhabditis_elegans
#uniprot_genes_Caenorhabditis_elegans $
#genesPlasmodium_falciparum  https://protists.ensembl.org/biomart/martview/a6bae05a1b023d6710e336592eb59d36 60   (ASM276v2)
#pdbPlasmodium_falciparum
#refPlasmodium_falciparum
#UCSCPlasmodium_falciparum
#uniprotPlasmodium_falciparum 
#uniprot_genes_Plasmodium_falciparum  $
#genesRattus_norvegicus https://www.ensembl.org/biomart/martview/3fd9943d70094de6d7fc0ff5f57dae94 113 Rat (mRatBN7.2) 
#pdbRattus_norvegicus
#refRattus_norvegicus
#UCSCRattus_norvegicus
#uniprotRattus_norvegicus
#uniprot_genes_Rattus_norvegicus   $
#genesSaccharomyces_cerevisiae https://fungi.ensembl.org/biomart/martview/39a165747bef091a9aac5a936c77c251 60  R64 - sacCer3
#pdbSaccharomyces_cerevisiae
#refSaccharomyces_cerevisiae
#UCSCSaccharomyces_cerevisiae
#uniprotSaccharomyces_cerevisiae
#uniprot_genes_Saccharomyces_cerevisiae  $
#genesSchizosaccharomyces_pombe https://fungi.ensembl.org/biomart/martview/39a165747bef091a9aac5a936c77c251 60 ASM294v2
#pdbSchizosaccharomyces_pombe
#refSchizosaccharomyces_pombe
#UCSCSchizosaccharomyces_pombe
#uniprotSchizosaccharomyces_pombe
#uniprot_genes_Schizosaccharomyces_pombe $
#genesSolanum_lycopersicum https://plants.ensembl.org/biomart/martview/2b9ddec0bca26c936126e569b1a1295e  60 (SL3.0)
#pdbSolanum_lycopersicum
#refSolanum_lycopersicum
#UCSCSolanum_lycopersicum
#uniprotSolanum_lycopersicum
#uniprot_genes_Solanum_lycopersicum  $
#genesStrongylocentrotus_purpuratus https://metazoa.ensembl.org/biomart/martview/6eea85f8b904fbc76b3a9d7ab3bb96d2 60  Strongylocentrotus purpuratus (Purple sea urchin, Spur 01) genes (Spur_5.0)
#pdbStrongylocentrotus_purpuratus
#refStrongylocentrotus_purpuratus
#UCSCStrongylocentrotus_purpuratus
#uniprotStrongylocentrotus_purpuratus
#uniprot_genes_Strongylocentrotus_purpuratus  $
#genessus_scrofa https://www.ensembl.org/biomart/martview/3fd9943d70094de6d7fc0ff5f57dae94 113 Pig genes (Sscrofa11.1)
#pdbsus_scrofa
#refsus_scrofa
#UCSCsus_scrofa
#uniprotsus_scrofa
#uniprot_genes_sus_scrofa  $
#genesToxoplasma_gondii https://protists.ensembl.org/biomart/martview/a6bae05a1b023d6710e336592eb59d36 60 (TGA4)
#pdbToxoplasma_gondii
#refToxoplasma_gondii
#UCSCToxoplasma_gondii
#uniprotToxoplasma_gondii
#uniprot_genes_Toxoplasma_gondii  $
#genesTrypanosoma_brucei_brucei  https://protists.ensembl.org/biomart/martview/a6bae05a1b023d6710e336592eb59d36 60 (TryBru_Apr2005_chr11)
#pdbTrypanosoma_brucei_brucei
#refTrypanosoma_brucei_brucei
#UCSCTrypanosoma_brucei_brucei
#uniprotTrypanosoma_brucei_brucei
#uniprot_genes_Trypanosoma_brucei_brucei 
#genesXenopus_tropicalis  https://www.ensembl.org/biomart/martview/3fd9943d70094de6d7fc0ff5f57dae94  113 Tropical clawed frog (UCB_Xtro_10.0) ▼
#pdbXenopus_tropicalis
#refXenopus_tropicalis
#UCSCXenopus_tropicalis
#uniprotXenopus_tropicalis
#uniprot_genes_Xenopus_tropicalis  $文件的顺序是：Entry	Entry name	Reviewed	Protein names	Gene name	organism	length
#https://www.uniprot.org/uniprotkb?query=%28organism_name%3AToxoplasma_gondii%29
    # #genesEctocarpus_sp._Ec32    Ectocarpus siliculosus str. Ec 32 (CCAP 1310/04) (GCA_000310025.1) (ASM31002v1) ▼
#genesEctocarpus_sp._Ec32 https://www.ncbi.nlm.nih.gov/datasets/gene/GCA_000310025.1/
	#genesKluyveromyces_lactis https://fungi.ensembl.org/Kluyveromyces_lactis_gca_000002515/Info/Index  Kluyveromyces lactis str. NRRL Y-1140 (ASM251v1) ▼
	#https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_000002515.2/



#genesArabidopsis_thaliana.txt https://plants.ensembl.org/biomart/martview/2b9ddec0bca26c936126e569b1a1295e   "TAIR10"
#genesDrosophila_melanogaster.txt http://asia.ensembl.org/biomart/martview/ffc1b51b784cf958b678dfc4a4731a3e  113 dm6
#genesGallus_gallus https://apr2022.archive.ensembl.org/biomart/martview/58809e781a3afb67a2db9bf606be1aad  106  Chicken genes (GRCg6a)
#genesChlorocebus_sabaeus   113  ChlSab1  Vervet-AGM (ChlSab1.1) 
#geneshuman https://www.ensembl.org/biomart/martview/3fd9943d70094de6d7fc0ff5f57dae94   113 (GRCh38.p14)
#genesMacaca_mulatta  https://www.ensembl.org/biomart/martview/3fd9943d70094de6d7fc0ff5f57dae94 113 Mmul_10 - rheMac10
#genesmouse https://nov2020.archive.ensembl.org/biomart/martview/11f032c236fa7fc326a1817dee60736e 102 GRCm38 - mm10
#genesOryza_sativa   https://plants.ensembl.org/biomart/martview/2b9ddec0bca26c936126e569b1a1295e 60  (IRGSP-1.0)
#genesPan_traglodytes https://www.ensembl.org/biomart/martview/3fd9943d70094de6d7fc0ff5f57dae94 113 Chimpanzee genes (Pan_tro_3.0)   #########
#genesZea_mays https://plants.ensembl.org/biomart/martview/2b9ddec0bca26c936126e569b1a1295e 60 (Zm-B73-REFERENCE-NAM-5.0)



# # #chimphttp://useast.ensembl.org/biomart/martview/54c798e06daa89b37378af92404d7adf 112 Pan_tro_3.0
# # #zebrafishhttp://may2015.archive.ensembl.org/biomart/martview/84ea87f76305b6cae6ab8ca3f1c69826 GRCz10/danRer10
# # #chickenhttps://apr2022.archive.ensembl.org/biomart/martview/58809e781a3afb67a2db9bf606be1aad  GRCg6a/galGal6
# # #cynomolgushttp://useast.ensembl.org/biomart/martview/54c798e06daa89b37378af92404d7adf 112
# # #maizehttps://nov2020-plants.ensembl.org/biomart/martview/906720c0d1dfd493210d4f6750047291
# # #flyhttp://jan2024.archive.ensembl.org/biomart/martview/087e12dede9a9f9f681702b754afcae5
# # #marmosethttp://may2015.archive.ensembl.org/biomart/martview/84ea87f76305b6cae6ab8ca3f1c69826
# # #rhesushttps://useast.ensembl.org/biomart/martview/45a1148f4a7082997e275fa16019ece2  112   Mmul_10
# # #ricehttps://plants.ensembl.org/biomart/martview/8cfcdb4ff2f653e19c58e6ddbb8f97b4   59
# # #arabidopsishttps://plants.ensembl.org/biomart/martview/8cfcdb4ff2f653e19c58e6ddbb8f97b4 59
# # #vervethttps://useast.ensembl.org/biomart/martview/d908531fed14c15ccf599d64c3dbd042 112
