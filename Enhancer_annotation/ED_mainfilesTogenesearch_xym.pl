#生成文件：geneforsearch".$spes[$i].".txt";

# # # #@spes=("chicken","human","mouse","arabidopsis");
# # # #@spetags=("ENSDARG","ENSG","ENSMUSG","AT");
# # # @spes=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
# # # #@spetags=("ENSGALG","ENSDARG","AT","ENSMUSG");
# # # #$spe="fly";
# # # #$spetag="FBgn";

#(chipseq)/data1/xym/DataProcess$
@spes=("Arabidopsis_thaliana","Caenorhabditis_elegans","Chlorocebus_sabaeus","Danio_rerio","Drosophila_melanogaster",
"Ectocarpus_sp._Ec32","Gallus_gallus","Homo_sapiens","Mus_musculus","Kluyveromyces_lactis","Macaca_mulatta","Oryza_sativa",
 "Pan_traglodytes","Plasmodium_falciparum","Rattus_norvegicus","Saccharomyces_cerevisiae","Schizosaccharomyces_pombe",
 "Solanum_lycopersicum","Strongylocentrotus_purpuratus","sus_scrofa","Toxoplasma_gondii","Trypanosoma_brucei_brucei",
 "Xenopus_tropicalis","Zea_mays");

for($i=0;$i<scalar(@spes);$i++){
open AA,"geneMapdata".$spes[$i].".txt";
#ENSGALG00000042750	ENSGALT00000080481,ENSGALP00000056694,ND1-201,ND1,39116926,YP_009558652,QBG37922,BAE16026,QAV56922,BAD11052,BAD11114,BAD11039,AVM38890,AYU71493,BAC57575,ARN58715,ARN58754,ASF62418,ARJ60438,APD80705,APD80718,ARB52444,ALN98097,AMD61849,AMP87941,AOR53637,ALN98088,ALN98078,ALN98025,ALN98042,ALN98051,ALN98060,ALN98069,ALN98007,AJK91538,AJR19309,AJR19322,AKF00068,AKJ83408,AJK29994,AJD80461,AJD80409,AJD80422,AJD80435,AJD80448,AJD80370,AJD80383,AJD80396,AJD22570,AIJ02084,AIU94502,AJD22557,AHG53940,AHI95689,AHL24729,AIE44485,ADW41565,ADW41526,ADW41539,ADW41552,ADB07155,ADB07090,ADB07103,ADB07129,ADB07142,ADB07064,ADB07077,ADB07051,ADB06973,ADB07038,ADB06908,ADB06921,ADB06765,ADB06791,ADB06830,ADB06882,ADB06895,ADB06726,ADB06752,ADB06648,ADB06661,ADB06687,ADB06700,ADB06635,ADB06596,ADB06609,ADB06622,AAO44997,ABF70963,ADB06583,AAO44985,CAA36625,P18936,MT-ND1,MTND1,NADH1,Q7GTV3
#the first one is gene ID,then the second are protein ID and gene names?
%hashIDtoInfo=();
@genes=();
while(<AA>){
s/\r|\n//g;
@temp=split/\t/,$_;
#if(($temp[0]=~ /$spetags[$i]/i) && $temp[0] ne ""){
if($temp[0] ne ""){
$hashIDtoInfo{$temp[0]}=$temp[1];
push @genes,$temp[0];
}
}
close AA;


%hashGeneIDtoCell=();
open CELL,"geneinfo".$spes[$i].".txt";
#ENSG00000154263	HG3802,HG3863,HG38121,HG38138,HG38146,HG38180,HG38187,HG38195,HG38232,HG38240,HG38278,HG38323	ENST00000522787,ENSP00000468351,ABCA10-207,ABCA10,HGNC:30,10349,uc060jgz.1,K7ERP5,Q8WWZ4,E5RFN6,E5RFP5,ENST00000521538,ENSP00000466506,ABCA10-205,uc060jha.1,K7EMH5,ENST00000690296,ENSP00000509702,ABCA10-213,NM_001377321,NP_001364250,uc315rob.1,ENST00000269081,ENSP00000269081,ABCA10-201,NM_080282,NP_525021,uc010dfa.2,ENST00000522406,ENSP00000429853,ABCA10-206,uc010dfb.2,ENST00000524231,ABCA10-210,ENST00000519732,ABCA10-203,ENST00000518929,ENSP00000430341,ABCA10-202,uc060jhc.1,ENST00000523419,ENSP00000428032,ABCA10-208,uc060jhd.1,ENST00000588514,ABCA10-212,ENST00000523512,ENSP00000429945,ABCA10-209,uc060jhf.1,ENST00000524273,ENSP00000429795,ABCA10-211,uc060jhg.1,H0YBL8,ENST00000521526,ABCA10-204
while(<CELL>){
s/\r|\n//g;
@temp=split/\t/,$_;
if(!exists $hashGeneIDtoCell{$temp[0]}){
$hashGeneIDtoCell{$temp[0]}=$temp[1];
}
}
close CELL;



open II,">geneforsearch".$spes[$i].".txt";
foreach $gene (@genes){
if(!exists $hashGeneIDtoCell{$gene}){
print II $gene."\t".$hashIDtoInfo{$gene}."\t"."NO"."\n";
}else{
print II $gene."\t".$hashIDtoInfo{$gene}."\t".$hashGeneIDtoCell{$gene}."\n";
}
}
close II;

}

