$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
use 5.010;
open AVAILTRUE,"GSMfinaltrack.bed";
@tracks=();
%hashTrackToSingle=();
while(<AVAILTRUE>){
chomp($_);
s/\r|\n//g;
@temp=split/\t/,$_;
	if(!exists $hashTrackToSingle{$temp[1]}){
	$hashTrackToSingle{$temp[1]}=$temp[0];
	push @tracks,$temp[1];  #@tracks TF-Binding POL2 Histone
	}else{
	$hashTrackToSingle{$temp[1]}.=";".$temp[0];   #TF-BindingGSM1572289;GSM1572290
	}
}
close AVAILTRUE;
if(scalar(@tracks)<3){
print "Tracks are not enough\n";
}else{
mkdir("track");
	foreach $track (@tracks){
	mkdir("track"."/".$track);
	ASWTRACK($track,"3");   #这是调用了ASWTRACK子程序，即 ($ALL,$numtype)=@_;赋值为$track,"3"
	system("bedtools subtract -a ".$track.".bed -b standard_promotor_exon.bed -A>".$track."NoPX.bed");  #bedtools subtract将".$track.".bed中的外显子和启动子区去除
	$hashtmpGSMTotargetfile{$track}=$track."NoPX.bed";
	NORMTRACK($track,"3");	#这是调用了NORMTRACK子程序，即 ($GSM,$numtype)=@_;赋值为$track,"3"
	system("mv ENH-hs-".$track.".bed track/".$track);
#	unlink($track.".bed");
#	unlink($track."NoPX.bed");
	}

}
$duration = time - $start;
print "All are done: $duration s\n";


sub   ASWTRACK()
{   
	my ($ALL,$numtype)=@_;
	@allshortnames=split/\;/,$hashTrackToSingle{$ALL};  #TF-BindingGSM1572289;GSM1572290
	open $ALL,">".$ALL."merge.bed";
	foreach $aswshort (@allshortnames){
	open $aswshort,"Reshs/ENH-hs-".$aswshort.".bed";     #这是将Reshs所有的文件逐个打开后，一行行读取，然后把所有的文件读取到同一个文件中。等同于cat一个track的所有bed进行合并
	while(<$aswshort>){				
	s/\r|\n//g;
	@tmp=split/\t/,$_;
		if($tmp[0]=~ /^chr/ && ($tmp[2]-$tmp[1])<=2500){
		print $ALL $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".($tmp[2]-$tmp[1])."-".$tmp[$numtype]."\n";
		}
	}
	close $aswshort;
	}
	close $ALL;
	system("bedtools sort -i ".$ALL."merge.bed>".$ALL."sort.bed");   #先进行sort和merge，然后再进行ASW，以防染色体之间发生混乱ASW
	unlink($ALL."merge.bed");
	system("bedtools merge -i ".$ALL."sort.bed -c 4 -o collapse>".$ALL."pre.bed");
	unlink($ALL."sort.bed");
	open $ALL,$ALL."pre.bed";
	@ALLlines=();
	while(<$ALL>){
	s/\r|\n//g;
	@tmp=split/\t/,$_;
	$center=int(($tmp[1]+$tmp[2])/2);
	@groups=split/\,/,$tmp[3];
	$lengsum=0;
	$sigsum=0;
		foreach $group (@groups){
			@lensiginfo=split/\-/,$group;
			$lengsum+=$lensiginfo[0];
			$sigsum+=$lensiginfo[0]*$lensiginfo[1];
		}
	$meansig=$sigsum/$lengsum;
	$meanlen=$lengsum/scalar(@groups);
	$start=$center-int($meanlen/2);
	$end=$center+int($meanlen/2);
	push @ALLlines,$tmp[0]."\t".$start."\t".$end."\t".$meansig;
	
	}
	close $ALL;
	unlink($ALL."pre.bed");
	open $ALL,">".$ALL.".bed";
	print $ALL join("\n",@ALLlines);
	close $ALL;
}

sub   NORMTRACK()
{   
	my ($GSM,$numtype)=@_;
	@tmpRegions=();
	%hashtmpRegionSig=();
	open $GSM,$hashtmpGSMTotargetfile{$GSM};
	$lengsum=0;
	$sigsum=0;
	while(<$GSM>){
	s/\r|\n//g;
	@tmp=split/\t/,$_;
		if($tmp[0]=~ /^chr/ && ($tmp[2]-$tmp[1])<=2500){
			push @tmpRegions,$tmp[0]."\t".$tmp[1]."\t".$tmp[2];
			$hashtmpRegionSig{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}=$tmp[$numtype];
			$lengsum+=$tmp[2]-$tmp[1];
			$sigsum+=($tmp[2]-$tmp[1])*$tmp[$numtype];
		}
	
	}
	close $GSM;
	$meansig=$sigsum/$lengsum;
	open $GSM,">ENH-hs-".$GSM.".bed";
	foreach $tmpRegion (@tmpRegions){
	print $GSM $tmpRegion."\t".(10*$hashtmpRegionSig{$tmpRegion}/$meansig)."\n";
	}
	close $GSM;
}

#需要更改的是命名部分
