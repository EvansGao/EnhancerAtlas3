##the pl is used to caculate the enhancer counts of per enhancer bed

$start = time;
#@species=("Gallus_gallus","Mus_musculus","Arabidopsis_thaliana","Danio_rerio");
#@spevers=("calJac3","mm10","macFas6","ChlSab1_1","hg38","galGal6","danRer10","panTro5","rheMac10","dm6","TAIR10","IRGSP1","B73v4");
#@spenames=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
#@spes=("gg","mm","at","dr");
#@tracks=("combined","P300","POLR2A","Histone","TF-binding","DHS","FAIRE","MNase-seq","GRO-seq","CAGE","MPRA","STARR-seq","CHIA-PET");
#@tracks=("Consensus");
#@firstchrs=("chr1","chr1","chr1","chr1");

# # # @species=("marmoset","mouse","cynomolgus","vervet","human","chicken","zebrafish","chimp","rhesus","fly","arabidopsis","rice","maize");
# # # @spes=("cj","mm","mf","cs","hs","gg","dr","pt","rm","dm","at","osj","zm");
# # # @tracks=("Consensus");
# # # @firstchrs=("chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr2R","chr1","chr1","chr1");



@species=("Arabidopsis_thaliana","Caenorhabditis_elegans","Chlorocebus_sabaeus","Danio_rerio","Drosophila_melanogaster",
"Ectocarpus_sp._Ec32","Gallus_gallus","Homo_sapiens","Mus_musculus","Kluyveromyces_lactis","Macaca_mulatta","Oryza_sativa",
"Pan_traglodytes","Plasmodium_falciparum","Rattus_norvegicus","Saccharomyces_cerevisiae","Schizosaccharomyces_pombe",
"Solanum_lycopersicum","Strongylocentrotus_purpuratus","sus_scrofa","Toxoplasma_gondii","Trypanosoma_brucei_brucei",
"Xenopus_tropicalis","Zea_mays");
@spes=("TAIR10","ce11","ChlSab1","danRer11","dm6","ASM31002v1","galGal6","hg38","mm10","ASM251v1","rheMac10","IRGSP1",
"panTro5","ASM276v2","rn7","sacCer3","ASM294v2","SL3","Spur5","susScr11","TGA4","TryBru","Xtro10","Zm-B73");
@tracks=("Consensus");   #what is this?
@firstchrs=("chr1","chrI","chr1","chr1","chr2L","LG01","chr1","chr1","chr1","chra","chr1","chr1","chr1","chr1","chr1","chrI","chri","chr1","AAGJ06000001.1","chr1","TGME49_chrIa","chr1","chr1","chr1");


%hashspetrackTocells=();
%hashspetrackTocellenhs=();
for($i=0;$i<scalar(@species);$i++){
	my $dir="cells/".$spes[$i];
	opendir(DIR,$dir) or "can't open the file";
	@dir=readdir DIR;
	@dir=grep{$_ ne "." && $_ ne ".."} @dir;
	foreach $file (@dir){
		print $file."\n";
		my $seconddir="cells/".$spes[$i]."/".$file."/".$firstchrs[$i];
		opendir(DIR,$seconddir) or "can't open the file";
		@seconddir=readdir DIR;
		@seconddir=grep{$_ ne "." && $_ ne ".."} @seconddir;
		foreach $track (@seconddir){
			if($track eq "combined"){
			$track="Consensus";
			}
			
			if(!exists $hashspetrackTocells{$spes[$i]."\t".$track}){    #$hashspetrackTocells{TAIR10	Consensus}=Ectocarpus_sp._Ec32_gametophyte?
			$hashspetrackTocells{$spes[$i]."\t".$track}=$file;
			$hashspetrackTocellenhs{$spes[$i]."\t".$track}=&LINES($species[$i],$file.".bed");
			}else{
			$hashspetrackTocells{$spes[$i]."\t".$track}.=";".$file;
			$hashspetrackTocellenhs{$spes[$i]."\t".$track}.=";".&LINES($species[$i],$file.".bed");   #the function is used to caculate the enhancer counts of per enhancer bed
			}
		}
	}
}

open SPETRACK,">speciestrackinfo.txt";
for($i=0;$i<scalar(@species);$i++){
	for($j=0;$j<scalar(@tracks);$j++){
		if(exists $hashspetrackTocells{$spes[$i]."\t".$tracks[$j]}){
		print SPETRACK $species[$i]."\t".$tracks[$j]."\t".$hashspetrackTocells{$spes[$i]."\t".$tracks[$j]}."\t".$hashspetrackTocellenhs{$spes[$i]."\t".$tracks[$j]}."\n";
		
		}

	}
}
close SPETRACK;

$duration = time - $start;
print "All are done: $duration s\n";

sub  LINES()
{   
    my ($species,$filename)=@_;
$linenum=0;
open FILE,"download/enhancer/".$species."/".$filename;
while (<FILE>) { $linenum++; }
close FILE;
return $linenum;
}
#这是计算增强子文件有多少行，相当于是统计增强子个数