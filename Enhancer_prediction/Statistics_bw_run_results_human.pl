use Cwd;
my $pwd = cwd();
$cellname=$pwd; #$  #$pwd /data/xianym/software/enhdata/Kc167
$cellname=~ s/.*\///g;
print "###############";
print "$cellname\n\n\n\n"; #$cellname  Kc167 这里应该有很多以细胞系名称命名的文件
print "###############";
##get data information
@AllGSM=();#$
@GSMforSRR=();#$
@GSMforSRRctr=();#$
@GSMforNonSRR=();#$
%hashGSMtoNonSRRfile=();#$
%hashExperCtr=();#$
%hashGSMtrack=();#$
%hashGSMdatatype=();#$
%hashGSMtoGSE=();#$
open CELL,$cellname."St.bed";#########################如Kcal67St.bed，格式有十几列
while(<CELL>){
s/\r|\n//g;
@temp=split/\t/,$_;#$
$temp;#$
$temp[0]=~ s/^\s+|\s+$|\,$//g;
$temp[2]=~ s/^\s+|\s+$|\,$//g;
$temp[3]=~ s/^\s+|\s+$|\,$//g;
$temp[9]=~ s/^\s+|\s+$|\,$//g;
$temp[10]=~ s/^\s+|\s+$|\,$//g;
$temp[14]=~ s/^\s+|\s+$|\,$//g;
$temp[15]=~ s/^\s+|\s+$|\,$//g;
push @AllGSM,$temp[0];#$
@tmpGSE=split/\,/,$temp[3];#$
$hashGSMtoGSE{$temp[0]}=$tmpGSE[0];
$hashGSMdatatype{$temp[0]}=$temp[14];
	if($temp[10]=~ /ATAC|DNase/i){
	$hashGSMtrack{$temp[0]}="DHS";
	}elsif($temp[10]=~ /MNase/i){
	$hashGSMtrack{$temp[0]}="MNase_seq";
	}elsif($temp[10]=~ /(Cut&Run)|(Cut_Run)|(Cut-Run)/i && $temp[9]=~ /(POLR2A)|(RNA polymerase II)|Pol II|PolII|Pol2|RNAPII/i){
	$hashGSMtrack{$temp[0]}="Cut_Run_POL2";
    }elsif($temp[10]=~ /(Cut&Run)|(Cut_Run)|(Cut-Run)/i && $temp[9]=~ /(p300)/i){
	$hashGSMtrack{$temp[0]}="Cut_Run_P300";######
	}elsif($temp[10]=~ /(Cut&Run)|(Cut_Run)|(Cut-Run)/i && $temp[9]=~ /(H3K4me1)|(H3K27ac)/i){
	$hashGSMtrack{$temp[0]}="Cut_Run_Histone";
	}elsif($temp[10]=~ /(Cut&Run)|(Cut_Run)|(Cut-Run)/i && !($temp[9]=~ /(H3K4me1)|(H3K27ac)|(POLR2A)|(RNA polymerase II)|Pol II|(P300)|PolII|Pol2|RNAPII/i)){
	$hashGSMtrack{$temp[0]}="Cut_Run_TF_Binding";
	}elsif($temp[10]=~ /(Cut&TAG)|(Cut_TAG)|(Cut-TAG)/i && $temp[9]=~ /(POLR2A)|(RNA polymerase II)|Pol II|PolII|Pol2|RNAPII/i){
	$hashGSMtrack{$temp[0]}="Cut_TAG_POL2";
    }elsif($temp[10]=~ /(Cut&TAG)|(Cut_TAG)|(Cut-TAG)/i && $temp[9]=~ /(p300)/i){
	$hashGSMtrack{$temp[0]}="Cut_TAG_P300";######
	}elsif($temp[10]=~ /(Cut&TAG)|(Cut_TAG)|(Cut-TAG)/i && $temp[9]=~ /(H3K4me1)|(H3K27ac)/i){
	$hashGSMtrack{$temp[0]}="Cut_TAG_Histone";
	}elsif($temp[10]=~ /(Cut&TAG)|(Cut_TAG)|(Cut-TAG)/i && !($temp[9]=~ /(H3K4me1)|(H3K27ac)|(POLR2A)|(RNA polymerase II)|Pol II|PolII|Pol2|(P300)|RNAPII/i)){
	$hashGSMtrack{$temp[0]}="Cut_TAG_TF_Binding";
	}elsif($temp[10]=~ /NET/i){
	$hashGSMtrack{$temp[0]}="NET_seq";
	}elsif($temp[10]=~ /MPRA/i){
	$hashGSMtrack{$temp[0]}="MPRA_seq";
	}elsif($temp[10]=~ /ChIPexo/i){
	$hashGSMtrack{$temp[0]}="ChIPexo_seq";
	}elsif($temp[10]=~ /CHIA_PET|ChIA-PET/i){
	$hashGSMtrack{$temp[0]}="CHIA_PET_seq";
	}elsif($temp[10]=~ /HiChIP/i){
	$hashGSMtrack{$temp[0]}="HiChIP_Seq";
	}elsif($temp[10]=~ /ChAR/i){
	$hashGSMtrack{$temp[0]}="ChAR_seq";
	}elsif($temp[10]=~ /CAGE|nAnTi-CAGE/i){
	$hashGSMtrack{$temp[0]}="CAGE_seq";
	}elsif($temp[10]=~ /FAIRE/i){
	$hashGSMtrack{$temp[0]}="FAIRE_seq";
	}elsif($temp[10]=~ /STARR/i){
	$hashGSMtrack{$temp[0]}="STARR_Seq";
	}elsif($temp[10]=~ /GRO/i){
	$hashGSMtrack{$temp[0]}="GRO_Seq";
	}elsif($temp[10]=~ /PRO/i){
	$hashGSMtrack{$temp[0]}="PRO_Seq";
	}elsif($temp[10]=~ /ChIP/i && $temp[9]=~ /(POLR2A)|(RNA polymerase II)|Pol II|PolII|Pol2|RNAPII/i){
	$hashGSMtrack{$temp[0]}="ChIP_POL2";
    }elsif($temp[10]=~ /ChIP/i && $temp[9]=~ /(p300)/i){
	$hashGSMtrack{$temp[0]}="ChIP_P300";######
	}elsif($temp[10]=~ /ChIP/i && $temp[9]=~ /(H3K4me1)|(H3K27ac)/i){
	$hashGSMtrack{$temp[0]}="ChIP_Histone";
	}elsif($temp[10]=~ /ChIP/i && !($temp[9]=~ /(H3K4me1)|(H3K27ac)|(POLR2A)|(RNA polymerase II)|Pol II|PolII|Pol2|(P300)|RNAPII/i)){
	$hashGSMtrack{$temp[0]}="ChIP_TF_Binding";
	}

	if($temp[15] eq "UsingSRR"){
		push @GSMforSRR,$temp[0];
		if($temp[2]=~ /^GSM/){
		push @GSMforSRRctr,$temp[2];
		$hashExperCtr{$temp[0]}=$temp[2];
		}
	}elsif($temp[15] ne "UsingSRR"){
		push @GSMforNonSRR,$temp[0];
		$hashGSMtoNonSRRfile{$temp[0]}=$temp[15];
	}

}
#close CELL;

#re-check the wrong files and classify correct files
%hashnewRESfile=();
$newdir="bwhs";
opendir(DIR,$newdir) or "can't open the file";
my @newdir=readdir DIR;
@newdir=grep{$_ ne "." && $_ ne ".."} @newdir;
foreach $file (@newdir){
if($file=~ /^ENH\-hs/i){
@temp=split/\-|\./,$file;
$hashnewRESfile{$temp[2]}=$file; #原本代码$hashnewRESfile{$temp[2]}=$temp[3];
}
}
open WRONG,">bwwrongList.bed";
open AVAIL,">bwavail.bed";
open AVAILTRUE,">bwavailTrue.bed";
open FINALTRACK,">bwfinaltrack.bed";
open WRONGTRACK,">bwwrongtrack.bed";
foreach $gsm (@AllGSM){
	if(exists $hashnewRESfile{$gsm} && $hashnewRESfile{$gsm} eq "bw"){
	print AVAIL $gsm."\t".$hashGSMtrack{$gsm}."\n";
	}else{
	print WRONG $gsm."\t".$hashnewRESfile{$gsm}."\n";
	}
	print AVAILTRUE $gsm."\t".$hashGSMtrack{$gsm}."\n";
	if(exists $hashnewRESfile{$gsm} && $hashnewRESfile{$gsm} =~ /bw/i){
	print FINALTRACK $gsm."\t".$hashGSMtrack{$gsm}."\n";
	}else{
	print WRONGTRACK $gsm."\t".$hashGSMtrack{$gsm}."\n";
	}
}
close WRONG;
close AVAIL;
close AVAILTRUE;
close FINALTRACK;
close WRONGTRACK;
#$filetype=$temp[3]; 
#if(exists $hashnewRESfile{$gsm} && $filetype eq "bed"){}#
#if(exists $hashnewRESfile{$gsm} && $hashnewRESfile{$gsm} eq "bed"){