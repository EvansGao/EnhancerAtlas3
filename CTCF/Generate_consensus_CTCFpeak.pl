	#不用change
	#这个文件是用来得到共识CTCF

	$start = time;
	$dir="preCTCF";
	opendir(DIR,$dir) or "Can't open the file";
	@dir=readdir DIR;
    @dir=grep{$_ ne "." && $_ ne ".."} @dir;  #读取文件时，通常会有.及..两个文件，因此，为了防止读取这两个文件，选择过滤掉这两个文件
    #print $dir;
    print "\n----------\n";
    $unionbedg="";
    $cellnum=0;
    foreach $file (@dir){
        print "$file\n";
        print "###############\n";
        $unionbedg.=" preCTCF/".$file;  #前面加空格是为让文件名称与文件名称之间有空格隔开。相当于，作用
        $cellnum++;
    }
    
    print "$cellnum\n";
    print "$unionbedg\n";
    system("bedtools unionbedg -i ".$unionbedg." >allCTCF.bed");

    open AA,"allCTCF.bed";
	open BB,">common.bed";
	while(<AA>){
	chomp($_);
	$num=0;
	@temp=split/\t/,$_;
	for($i=3;$i<scalar(@temp);$i++){
	if($temp[$i]>0){
	$num++;
	}
	}
	if($num>=int($cellnum/2)){
	print BB $temp[0]."\t".$temp[1]."\t".$temp[2]."\n";
	}
	}
	close AA;
	close AA;
	close BB;
	system("bedtools sort -i common.bed>sortcommon.bed");
	unlink("common.bed");
	system("bedtools merge -i sortcommon.bed>commonall.bed");
	$duration = time - $start;
	print "Done: $duration s\n";