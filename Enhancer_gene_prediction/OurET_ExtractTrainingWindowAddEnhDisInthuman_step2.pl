#这个与mouse、dm转的运行训练集的数据一样，只是文件命名不同


$start = time;
use List::Util qw[min max];
use List::MoreUtils qw(uniq);
use List::Util qw(shuffle);
#get the promoters from ensembl gene information
open GENE,"human_ensembl.txt";
%hashgenepromoter=();
%hashpromotergene=();
while(<GENE>){
chomp($_);
@temp=split/\t/,$_;
if($temp[0] ne "" && !exists $hashgenepromoter{$temp[0]} && $temp[5] eq "1"){
$hashgenepromoter{$temp[0]}="chr".$temp[1]."\t".max(1,($temp[2]-5000))."\t".($temp[2]+500);
$hashpromotergene{"chr".$temp[1].":".max(1,($temp[2]-5000))."-".($temp[2]+500)}=$temp[0];
}elsif($temp[0] ne "" && !exists $hashgenepromoter{$temp[0]} && $temp[5] eq "-1"){
$hashgenepromoter{$temp[0]}="chr".$temp[1]."\t".($temp[3]-500)."\t".($temp[3]+5000);
$hashpromotergene{"chr".$temp[1].":".($temp[3]-500)."-".($temp[3]+5000)}=$temp[0];
}

}
close GENE;
#get the promoters from ensembl gene information done

$cell="hsK562";
$standard="hsK562_CHIA-PET_Rep0.bed";
#$standard="GM12878_HiC.bed";
#integrating the pairs with errors into the standard pair format
open AA,"cellAllfeatures".$cell.".bed";
%hashpairgene=();
@pairgenes=();

#hash pair to the gene name and feature scores
%hashpairToNameScore=();
#hash pair to the gene name and feature scores done
%hashpairtoCor=();
%hashpairtoGeneScore=();
%hashpairtoEnhScore=();
%hashpairtoEnhWindowScore=();
%hashpairtoGeneWindowScore=();
open BB,">Pairs_".$cell.".bed";
while(<AA>){
chomp($_);
@tempinfo=split/\$/,$_;
@tempones=split/\:|\-|_/,$tempinfo[0];
@tempfeatures=split/\t/,$tempinfo[4];
$featurescores=join("\t",@tempfeatures);
$hashpairToNameScore{$tempinfo[0]}=$tempinfo[0]."\$".$tempinfo[1]."\$".$tempinfo[2]."\$".$tempinfo[3]."\$".$featurescores;
$hashpairtoCor{$tempones[3]."\t".$tempones[0].":".$tempones[1]."-".$tempones[2]}=$tempfeatures[1];
$hashpairtoGeneScore{$tempones[3]."\t".$tempones[0].":".$tempones[1]."-".$tempones[2]}=$tempfeatures[2];
$hashpairtoEnhScore{$tempones[3]."\t".$tempones[0].":".$tempones[1]."-".$tempones[2]}=$tempfeatures[3];
$hashpairtoEnhWindowScore{$tempones[3]."\t".$tempones[0].":".$tempones[1]."-".$tempones[2]}=$tempfeatures[4];
$hashpairtoGeneWindowScore{$tempones[3]."\t".$tempones[0].":".$tempones[1]."-".$tempones[2]}=$tempfeatures[5];
if(!exists $hashpairgene{$tempones[3]}){
$hashpairgene{$tempones[3]}=$tempones[0].":".$tempones[1]."-".$tempones[2];
push @pairgenes,$tempones[3];
}else{
$hashpairgene{$tempones[3]}.=";".$tempones[0].":".$tempones[1]."-".$tempones[2];
}
print BB $tempones[0]."\t".$tempones[1]."\t".$tempones[2]."\t".$hashgenepromoter{$tempones[3]}."\n";
}
close AA;
close BB;
#integrating the JEME pairs with errors into the standard pair format done

#calculating their own overlapped pairs for JEME and HiC files
$pairsname="Pairs_".$cell.".bed";
&OVERLAP($pairsname,$standard);
#calculating their own overlapped pairs for JEME and HiC files done

#convert overlapped pair format into JEME and hash genes with enhancer regions from overlapped pairs
open AA,"Int_".$pairsname;
%hashoverlapgenes=();
@overlapgenes=();
open BB,">Intpairs_".$pairsname;
while(<AA>){
chomp($_);
@temp=split/\:|\-|\t/,$_;
if(exists $hashpromotergene{$temp[3].":".$temp[4]."-".$temp[5]}){
print BB $temp[0].":".$temp[1]."-".$temp[2]."_".$hashpromotergene{$temp[3].":".$temp[4]."-".$temp[5]}."\n";
	$gene=$hashpromotergene{$temp[3].":".$temp[4]."-".$temp[5]};
	if(!exists $hashoverlapgenes{$gene}){
	$hashoverlapgenes{$gene}=$temp[0].":".$temp[1]."-".$temp[2];
	push @overlapgenes,$gene;
	}elsif(exists $hashoverlapgenes{$gene}){
	$hashoverlapgenes{$gene}.=";".$temp[0].":".$temp[1]."-".$temp[2];
	}
}
if(exists $hashpromotergene{$temp[0].":".$temp[1]."-".$temp[2]}){
print BB $temp[3].":".$temp[4]."-".$temp[5]."_".$hashpromotergene{$temp[0].":".$temp[1]."-".$temp[2]}."\n";
	$gene=$hashpromotergene{$temp[0].":".$temp[1]."-".$temp[2]};
	if(!exists $hashoverlapgenes{$gene}){
	$hashoverlapgenes{$gene}=$temp[3].":".$temp[4]."-".$temp[5];
	push @overlapgenes,$gene;
	}elsif(exists $hashoverlapgenes{$gene}){
	$hashoverlapgenes{$gene}.=";".$temp[3].":".$temp[4]."-".$temp[5];
	}
}

}
close AA;
close BB;
open CC,">trainingWindowIntdisStForR_".$cell.".bed";
open DD,">trainingIntrandomdisStForR_".$cell.".bed";
#positive training data
$positivenum=0;
%hashpositivegene=();
%hashpositiveenhancer=();
print scalar(@overlapgenes)."\n";
foreach $overlapgene (@overlapgenes){
$hashpositivegene{$overlapgene}="";
@enhs=split/\;/,$hashoverlapgenes{$overlapgene};
	foreach $enh (@enhs){
#		if($hashpairtoCor{$overlapgene."\t".$enh}>0 && $hashpairtoGeneScore{$overlapgene."\t".$enh}>0 && $hashpairtoEnhScore{$overlapgene."\t".$enh}>0){
#		if($hashpairtoGeneScore{$overlapgene."\t".$enh}>0 && $hashpairtoEnhScore{$overlapgene."\t".$enh}>0){
#		print CC $hashpairToNameScore{$enh."_".$overlapgene}."\t"."1"."\n";
			if($hashpairToNameScore{$enh."_".$overlapgene} ne ""){
			print CC $hashpairToNameScore{$enh."_".$overlapgene}."\t"."1"."\n";
			$hashpositiveenhancer{$enh}="";
			$positivenum++;
			}

#		}
	}
#$randenh=$enhs[rand @enhs];
}
#positive training data done

#negative training data
@pairgenes=shuffle(@pairgenes);
$num=0;
$randomnum=0;
print scalar(@pairgenes)."\n";
foreach $pairgene (@pairgenes){
	@enhs=split/\;/,$hashpairgene{$pairgene};
	foreach $enh (@enhs){
		if($randomnum<$positivenum){
		print DD $hashpairToNameScore{$enh."_".$pairgene}."\t"."1"."\n";
		$randomnum++;
		}elsif($randomnum>=$positivenum && $randomnum<2*$positivenum){
		print DD $hashpairToNameScore{$enh."_".$pairgene}."\t"."2"."\n";
		$randomnum++;
		}
	}
	
	if(!exists $hashpositivegene{$pairgene}){
	foreach $enh (@enhs){
		if(!exists $hashpositiveenhancer{$enh} && $num<=$positivenum){
		print CC $hashpairToNameScore{$enh."_".$pairgene}."\t"."2"."\n";
		$num++;
		}
	}
#	$randenh=$enhs[rand @enhs];
	}
#	if($num>(5*$positivenum)){last;}
#	if($num>$positivenum){last;}
}
close CC;
close DD;
print "aa".$positivenum."\t".$num."\n";
#negative training data done


#convert overlapped pair format into JEME and extract genes with enhancer regions from overlapped pairs done



sub   OVERLAP()
{  
    my ($pairsname,$standard)=@_;
    $filenum1=int(rand(1000000));
    $filenum11=int(rand(1000000));
    %hashenh1=();
    %hashpair1=();
    open $filenum1,$pairsname;
    unlink("enh_".$pairsname);
    open $filenum11,">enh_".$pairsname;
    while(<$filenum1>){
    chomp($_);
    @temp=split/\t/,$_;
    if(!exists $hashenh1{$temp[0]."\t".$temp[1]."\t".$temp[2]}){
    $hashenh1{$temp[0]."\t".$temp[1]."\t".$temp[2]}="";
    print $filenum11 $temp[0]."\t".$temp[1]."\t".$temp[2]."\n";
    }
    
    if(!exists $hashenh1{$temp[3]."\t".$temp[4]."\t".$temp[5]}){
    $hashenh1{$temp[3]."\t".$temp[4]."\t".$temp[5]}="";
    print $filenum11 $temp[3]."\t".$temp[4]."\t".$temp[5]."\n";
    }
    $hashpair1{$temp[0].":".$temp[1]."-".$temp[2]."\t".$temp[3].":".$temp[4]."-".$temp[5]}="";
    $hashpair1{$temp[3].":".$temp[4]."-".$temp[5]."\t".$temp[0].":".$temp[1]."-".$temp[2]}="";
    }
    close $filenum1;
    close $filenum11;
    
    $filenum2=int(rand(1000000));
    $filenum22=int(rand(1000000));
    %hashenh2=();
	@pairs=();
    open $filenum2,$standard;
    unlink("enh_".$standard);
    open $filenum22,">enh_".$standard;
    while(<$filenum2>){
    chomp($_);
    @temp=split/\t/,$_;
    push @pairs,$temp[0].":".$temp[1]."-".$temp[2]."\t".$temp[3].":".$temp[4]."-".$temp[5];
    if(!exists $hashenh2{$temp[0]."\t".$temp[1]."\t".$temp[2]}){
    $hashenh2{$temp[0]."\t".$temp[1]."\t".$temp[2]}="";
    print $filenum22 $temp[0]."\t".$temp[1]."\t".$temp[2]."\n";
    }
    
    if(!exists $hashenh2{$temp[3]."\t".$temp[4]."\t".$temp[5]}){
    $hashenh2{$temp[3]."\t".$temp[4]."\t".$temp[5]}="";
    print $filenum22 $temp[3]."\t".$temp[4]."\t".$temp[5]."\n";
    }
    }
    close $filenum2;
    close $filenum22;

    system("bedtools sort -i enh_".$pairsname.">enhsort_".$pairsname);
    system("bedtools sort -i enh_".$standard.">enhsort_".$standard);
    $randnum=int(rand(1000000));
    system("bedtools intersect -a enhsort_".$standard." -b enhsort_".$pairsname." -wa -wb>intersect".$randnum.".bed");
    open $randnum,"intersect".$randnum.".bed";
    %hashenh2toenh1=();
    while(<$randnum>){
    chomp($_);
    @temp=split/\t/,$_;
    if(!exists $hashenh2toenh1{$temp[0].":".$temp[1]."-".$temp[2]}){
    $hashenh2toenh1{$temp[0].":".$temp[1]."-".$temp[2]}=$temp[3].":".$temp[4]."-".$temp[5];
    }else{
    $hashenh2toenh1{$temp[0].":".$temp[1]."-".$temp[2]}.=";".$temp[3].":".$temp[4]."-".$temp[5];
    }
	}
    close $randnum;
    open $filenum1,">Int_".$pairsname;
    open $filenum2,">Int_".$standard;
    open $randnum,">Int_combined".$randnum.".bed";
    @fileones=();
    @filetwos=();
    foreach $pair (@pairs){
    @temp=split/\t/,$pair;
    @formers=split/\;/,$hashenh2toenh1{$temp[0]};
    @latters=split/\;/,$hashenh2toenh1{$temp[1]};
    if($temp[0] eq "chr1:712491-714111"){
    print scalar(@formers)."\t".scalar(@latters)."\n";
    }
    	foreach $former (@formers){
    		foreach $latter (@latters){
    			if(exists $hashpair1{$former."\t".$latter}){
    			push @fileones,$former."\t".$latter;
    			push @filetwos,$pair;
    			print $randnum $former."\t".$latter."|".$pair."\n";
    			}
    		}
    	}
    }
    close $randnum;
    @fileones=uniq(@fileones);
    @filetwos=uniq(@filetwos);
    foreach $fileone (@fileones){
    print $filenum1 $fileone."\n";
    }
    foreach $filetwo (@filetwos){
    print $filenum2 $filetwo."\n";
    }
    close $filenum1;
    close $filenum2;

    
    }


$duration = time - $start;
print "all is done: $duration s\n";

