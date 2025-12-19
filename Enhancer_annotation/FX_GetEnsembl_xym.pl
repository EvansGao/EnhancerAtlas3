system("rm -R Ensembl");
#@spes=("hs","mm","at","dr");
#@spevers=("hg38","mm10","TAIR10","danRer10");
# # # @spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
# # # @spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");

@spevers=("TAIR10","ce11","ChlSab1","danRer11","dm6","ASM31002v1","galGal6","hg38","mm10","ASM251v1","rheMac10","IRGSP1",
 "panTro5","ASM276v2","rn7","sacCer3","ASM294v2","SL3","Spur5","susScr11","TGA4","TryBru","Xtro10","Zm-B73");
@spes=("Arabidopsis_thaliana","Caenorhabditis_elegans","Chlorocebus_sabaeus","Danio_rerio","Drosophila_melanogaster",
 "Ectocarpus_sp._Ec32","Gallus_gallus","Homo_sapiens","Mus_musculus","Kluyveromyces_lactis","Macaca_mulatta","Oryza_sativa",
  "Pan_traglodytes","Plasmodium_falciparum","Rattus_norvegicus","Saccharomyces_cerevisiae","Schizosaccharomyces_pombe",
  "Solanum_lycopersicum","Strongylocentrotus_purpuratus","sus_scrofa","Toxoplasma_gondii","Trypanosoma_brucei_brucei",
  "Xenopus_tropicalis","Zea_mays");


mkdir("Ensembl");
for ($i = 0; $i < scalar(@spes); $i++) {
    open AA, "../project/Refgene/ensembl_" . $spevers[$i] . ".txt";
    %chrgene = ();
    %hashgene = ();
    %hashgeneExonstart = ();
    %hashgeneName = ();
    %hashgeneStrand = ();
    %hashgeneStartEnd = ();
    @chrarray = ();

    while (<AA>) {
        chomp($_);
        @tmp = split /\t/, $_;

        if (!($spevers[$i] =~ /ASM31002v1|TGA4|Spur5/i)) {
            if (!($tmp[6] =~ /Name/i)) {
                if (!exists $chrgene{"chr" . $tmp[6]}) {
                    $chrgene{"chr" . $tmp[6]} = $tmp[0];
                    push @chrarray, "chr" . $tmp[6];
                } else {
                    if (!($chrgene{"chr" . $tmp[6]} =~ /$tmp[0]/)) {
                        $chrgene{"chr" . $tmp[6]} .= "," . $tmp[0];
                    }
                }
            }
        } else {
            if (!($tmp[6] =~ /Name/i)) {
                if (!exists $chrgene{$tmp[6]}) {
                    $chrgene{$tmp[6]} = $tmp[0];
                    push @chrarray, $tmp[6];
                } else {
                    if (!($chrgene{$tmp[6]} =~ /$tmp[0]/)) {
                        $chrgene{$tmp[6]} .= "," . $tmp[0];
                    }
                }
            }
        }

        $hashgeneExonstart{$tmp[0] . "\t" . $tmp[1]} = $tmp[2];

        if (!exists $hashgene{$tmp[0]}) {
            $hashgene{$tmp[0]} = $tmp[1];
            $hashgeneName{$tmp[0]} = $tmp[7];
            if ($tmp[5] eq "1") {
                $hashgeneStrand{$tmp[0]} = "+";
            } else {
                $hashgeneStrand{$tmp[0]} = "-";
            }

            $hashgeneStartEnd{$tmp[0]} = $tmp[3] . "\t" . $tmp[4];
        } else {
            $hashgene{$tmp[0]} .= "," . $tmp[1];
        }
    }

    close AA;
    @chrarray = sort { $a <=> $b } @chrarray;

    mkdir("Ensembl/" . $spes[$i]);

    foreach $chrom (@chrarray) {
        if (length($chrom) < 8) {
            open my $chrom_fh, ">", "Ensembl/" . $spes[$i] . "/" . $chrom . ".bed";

            @gene = split /,/, $chrgene{$chrom};

            foreach $singlegene (@gene) {
                @exonstarts = split /,/, $hashgene{$singlegene};
                @exonstarts = sort { $a <=> $b } @exonstarts;
                $exonS = "";
                $exonE = "";

                for ($j = 0; $j < scalar(@exonstarts); $j++) {
                    $exonS .= "," . $exonstarts[$j];
                    $exonE .= "," . $hashgeneExonstart{$singlegene . "\t" . $exonstarts[$j]};
                }

                $exonS =~ s/^\,//g;
                $exonE =~ s/^\,//g;

                print $chrom_fh $chrom . "\t" . $hashgeneStartEnd{$singlegene} . "\t" . $hashgeneStrand{$singlegene} . "\t" . $hashgeneName{$singlegene} . "\t" . $exonS . "\t" . $exonE . "\t" . $singlegene . "\n";
            }

            close $chrom_fh;

            system("bedtools sort -i Ensembl/" . $spes[$i] . "/" . $chrom . ".bed > Ensembl/" . $spes[$i] . "/" . $chrom);
            unlink("Ensembl/" . $spes[$i] . "/" . $chrom . ".bed");
        }
    }
}




#########################

# # # for($i=0;$i<scalar(@spes);$i++){
# # # open AA,"../project/Refgene/ensembl_".$spevers[$i].".txt";
# # # %chrgene=();
# # # %hashgene=();
# # # %hashgeneExonstart=();
# # # %hashgeneName=();
# # # %hashgeneStrand=();
# # # %hashgeneStartEnd=();
# # # @chrarray=();
# # # while(<AA>){
# # # chomp($_);
# # # @tmp=split/\t/,$_;
# # # if(&& !($spevers[$i] =~ /ASM31002v1|TGA4|Spur5/i)){
# # # 	if(!($tmp[6]=~ /Name/i)){
# # # 		if(!exists $chrgene{"chr".$tmp[6]}){   #注意这里的挑选出不是以chr开头的物种"ASM31002v1","TGA4","Spur5"
# # # 				$chrgene{"chr".$tmp[6]}=$tmp[0];
# # # 			push @chrarray,"chr".$tmp[6];
# # # 		}else{
# # # 			if(!($chrgene{"chr".$tmp[6]}=~ /$tmp[0]/)){
# # # 			$chrgene{"chr".$tmp[6]}.=",".$tmp[0];
# # # 			}
# # # 		}
# # # 	}
# # # }else{
# # # 	if(!($tmp[6]=~ /Name/i)){
# # # 		if(!exists $chrgene{$tmp[6]}){   #注意这里的挑选出不是以chr开头的物种"ASM31002v1","TGA4","Spur5"
# # # 				$chrgene{$tmp[6]}=$tmp[0];
# # # 			push @chrarray,$tmp[6];
# # # 		}else{
# # # 			if(!($chrgene{$tmp[6]}=~ /$tmp[0]/)){
# # # 			$chrgene{$tmp[6]}.=",".$tmp[0];
# # # 			}
# # # 		}
# # # 	}

# # # }

# # # $hashgeneExonstart{$tmp[0]."\t".$tmp[1]}=$tmp[2];

# # # if(!exists $hashgene{$tmp[0]}){
# # # $hashgene{$tmp[0]}=$tmp[1];
# # # $hashgeneName{$tmp[0]}=$tmp[7];
# # # 		if($tmp[5] eq "1"){
# # # 		$hashgeneStrand{$tmp[0]}="+";
# # # 		}else{
# # # 		$hashgeneStrand{$tmp[0]}="-";
# # # 		}

# # # $hashgeneStartEnd{$tmp[0]}=$tmp[3]."\t".$tmp[4];
# # # }else{
# # # $hashgene{$tmp[0]}.=",".$tmp[1];
# # # }

# # # }
# # # close AA;
# # # @chrarray=sort{$a<=>$b} @chrarray;

# # # mkdir("Ensembl/".$spes[$i]);
# # # foreach $chrom (@chrarray){
# # # 	if(length($chrom)<8){
# # # 	open $chrom,">Ensembl/".$spes[$i]."/".$chrom.".bed";
# # # 	@gene=split/\,/,$chrgene{$chrom};
# # # 	foreach $singlegene (@gene){
# # # 		@exonstarts=split/\,/,$hashgene{$singlegene};
# # # 		@exonstarts=sort{$a<=>$b} @exonstarts;
# # # 		$exonS="";
# # # 		$exonE="";
# # # 		for($j=0;$j<scalar(@exonstarts);$j++){
# # # 		$exonS.=",".$exonstarts[$j];
# # # 		$exonE.=",".$hashgeneExonstart{$singlegene."\t".$exonstarts[$j]};
# # # 		}
# # # 		$exonS=~ s/^\,//g;
# # # 		$exonE=~ s/^\,//g;
# # # 		print $chrom $chrom."\t".$hashgeneStartEnd{$singlegene}."\t".$hashgeneStrand{$singlegene}."\t".$hashgeneName{$singlegene}."\t".$exonS."\t".$exonE."\t".$singlegene."\n";
# # # 	}
# # # 	close $chrom;
# # # 	system("bedtools sort -i Ensembl/".$spes[$i]."/".$chrom.".bed >Ensembl/".$spes[$i]."/".$chrom);
# # # 	unlink("Ensembl/".$spes[$i]."/".$chrom.".bed");
# # # 	}
# # # }

# # # }

#ensembl_Spur5   https://metazoa.ensembl.org/biomart/martview/706058286df8c0b413b7b95cc61122d4
#genesDanio_rerio http://asia.ensembl.org/biomart/martview/d0b53f53fe7f2a87cc75d5b7ae833c4d  113 GRcz11 
#ensembl_danRer11
#genesCaenorhabditis_elegans https://asia.ensembl.org/biomart/martview/ffc1b51b784cf958b678dfc4a4731a3e  113  WBcel235/ce11
#ensembl_ce11   
#genesPlasmodium_falciparum  https://protists.ensembl.org/biomart/martview/a6bae05a1b023d6710e336592eb59d36 60   (ASM276v2)
#ensembl_ASM276v2
#genesRattus_norvegicus https://www.ensembl.org/biomart/martview/3fd9943d70094de6d7fc0ff5f57dae94 113 Rat (mRatBN7.2) 
#ensembl_rn7
#genesSaccharomyces_cerevisiae https://fungi.ensembl.org/biomart/martview/39a165747bef091a9aac5a936c77c251 60  R64 - sacCer3
#ensembl_sacCer3
#genesSchizosaccharomyces_pombe https://fungi.ensembl.org/biomart/martview/39a165747bef091a9aac5a936c77c251 60 ASM294v2
#ensembl_ASM294v2
#genesSolanum_lycopersicum https://plants.ensembl.org/biomart/martview/2b9ddec0bca26c936126e569b1a1295e  60 (SL3.0)
#ensembl_SL3   
#genesStrongylocentrotus_purpuratus https://metazoa.ensembl.org/biomart/martview/6eea85f8b904fbc76b3a9d7ab3bb96d2 60  Strongylocentrotus purpuratus (Purple sea urchin, Spur 01) genes (Spur_5.0)
#ensembl_  $
#genessus_scrofa https://www.ensembl.org/biomart/martview/3fd9943d70094de6d7fc0ff5f57dae94 113 Pig genes (Sscrofa11.1)
#ensembl_susScr11
#genesToxoplasma_gondii https://protists.ensembl.org/biomart/martview/a6bae05a1b023d6710e336592eb59d36 60 (TGA4)
#ensembl_TGA4
#genesTrypanosoma_brucei_brucei  https://protists.ensembl.org/biomart/martview/a6bae05a1b023d6710e336592eb59d36 60 (TryBru_Apr2005_chr11)
#ensembl_TryBru 
#genesXenopus_tropicalis  https://www.ensembl.org/biomart/martview/3fd9943d70094de6d7fc0ff5f57dae94  113 Tropical clawed frog (UCB_Xtro_10.0) ▼
#ensembl_Xtro10
#uniprot_genes_Xenopus_tropicalis  $文件的顺序是：Entry	Entry name	Reviewed	Protein names	Gene name	organism	length
#https://www.uniprot.org/uniprotkb?query=%28organism_name%3AToxoplasma_gondii%29
    # #genesEctocarpus_sp._Ec32    Ectocarpus siliculosus str. Ec 32 (CCAP 1310/04) (GCA_000310025.1) (ASM31002v1) ▼
#genesEctocarpus_sp._Ec32 https://www.ncbi.nlm.nih.gov/datasets/gene/GCA_000310025.1/
	#genesKluyveromyces_lactis https://fungi.ensembl.org/Kluyveromyces_lactis_gca_000002515/Info/Index  Kluyveromyces lactis str. NRRL Y-1140 (ASM251v1) ▼
	#https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_000002515.2/