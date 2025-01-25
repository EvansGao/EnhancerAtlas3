$start = time;
use File::Copy;
use List::MoreUtils qw(uniq);
$dir="track";
opendir(DIR,$dir);
@track=readdir DIR;
@track=grep{$_ ne "." && $_ ne ".."} @track;
use Cwd;
my $pwd = cwd();
$cellname=$pwd;
$cellname=~ s/.*\///g;
foreach $trk (@track){
open AA,"ENH-hs-".$cellname.$trk.".bed";
open BB,">ENH-hs-".$cellname.$trk."C.bed";
	while(<AA>){
	s/\r|\n//g;
	@temp=split/\t/,$_;
			for($i=$temp[1];$i<$temp[2];$i++){
			print BB $temp[0]."\t".$i."\t".($i+1)."\t".&ND($i,$temp[1],$temp[2],$temp[3])."\n";
			}
	}
close AA;
close BB;
system("bedGraphToBigWig ENH-hs-".$cellname.$trk."C.bed hg38.chrom.sizes ENH-hs-".$cellname.$trk."WebC.bw");
}

open AA,"ENH-hs-".$cellname.".bed";
open BB,">ENH-hs-".$cellname."C.bed";
while(<AA>){
s/\r|\n//g;
@temp=split/\t/,$_;
	for($i=$temp[1];$i<$temp[2];$i++){
	print BB $temp[0]."\t".$i."\t".($i+1)."\t".&ND($i,$temp[1],$temp[2],$temp[3])."\n";
	}
}
close AA;
close BB;
system("bedGraphToBigWig ENH-hs-".$cellname."C.bed hg38.chrom.sizes ENH-hs-".$cellname."WebC.bw");

mkdir("bwhs");
$dir="Reshs";
opendir(DIR,$dir);
@singlefiles=readdir DIR;
@singlefiles=grep{$_ ne "." && $_ ne ".."} @singlefiles;
foreach $singlefile (@singlefiles){
open AA,"Reshs/".$singlefile;
$singlebw=$singlefile;
$singlebw=~ s/bed$/bw/g;
open BB,">C".$singlefile;
	while(<AA>){
	s/\r|\n//g;
	@temp=split/\t/,$_;
			for($i=$temp[1];$i<$temp[2];$i++){
			print BB $temp[0]."\t".$i."\t".($i+1)."\t".&ND($i,$temp[1],$temp[2],$temp[3])."\n";
			}
	}
close AA;
close BB;
system("sort -k 1,1 -k2,2n C".$singlefile."> C".$singlefile."Sort");
system("bedGraphToBigWig C".$singlefile."Sort hg38.chrom.sizes bwhs/".$singlebw);
if(-e "ex_rerun_bw/ENH-hs-'$temp[0]'.bw"){unlink("CENH-hs-'$temp[0]'.bed");}
if(-e "ex_rerun_bw/ENH-hs-'$temp[0]'.bw"){unlink("CENH-hs-'$temp[0]'.bedSort");}
#unlink("C".$singlefile);
}
#这是将reshs目录下所有的bed文件都转化为bw文件

$duration = time - $start;
print "All is done: $duration s\n";

sub ND
{
  my($position,$start,$end,$score)=@_;
  $center=int(($end-$start)/2);
  if(($position-$start)<$center){
  $xx=$position-$start-$center;
  $xi=log($score/0.00001);
  $sigma=sqrt(($center**2)/(2*$xi));
  return $score*exp(-$xx**2/(2*$sigma**2));
  }elsif(($position-$start)>=$center){
  $xx=$position-$end+$center;
  $xi=log($score/0.00001);
  $sigma=sqrt((($end-$start-$center)**2)/(2*$xi));
  if($sigma==0){
  return 0;
  }else{
  return $score*exp(-$xx**2/(2*$sigma**2));
  }
  }
}
#ND是EAGLE算法