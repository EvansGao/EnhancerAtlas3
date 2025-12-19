open AA,"gwas_catalog_v1.0.2-associations_e111_r2024-03-28.tsv";
@SNPs=();
%hashSNP=();
while(<AA>){
s/\r|\n//g;
@temp=split/\t/,$_;
@temSNPs=split/\;\s+|\sx\s/,$temp[21];
	foreach $temSNP (@temSNPs){
		if(!exists $hashSNP{$temSNP} && $temp[11] ne "" && !($temp[11]=~ /\sx\s/) && $temp[3] ne "DATE" && !($temp[21]=~ /\;\s+|\sx\s|\,\s+/)){
		$hashSNP{$temSNP}="chr".$temp[11]."\t".$temp[12]."\t".($temp[12]+2);
		push @SNPs,$temSNP;
		}
	}
}
close AA;
print scalar(@SNPs)."\n";
open BB,">SNPs_hg38pre.bed";
foreach $snp (@SNPs){
print BB $hashSNP{$snp}."\t".$snp."\n";
}
close BB;
#system("liftOver SNPs_hg38.bed hg38Tohg19.over.chain SNPs_hg19pre.bed unMapped");
system("bedtools sort -i SNPs_hg38pre.bed>SNPs_hg38.bed");
#bedtools sort -i Motif_dm3.bed>Motif_dm.bed


#https://www.ebi.ac.uk/gwas/docs/file-downloads
	#All associations v1.0.2 - with added ontology annotations, GWAS Catalog study accession numbers and genotyping technology
	#Click to download
#https://www.ebi.ac.uk/gwas/api/search/downloads/alternative