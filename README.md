# EnhancerAtlas 3.0: an updated enhancer resource for 1428 tissue/cell types across 24 species
![headnew](https://github.com/user-attachments/assets/7255d18c-0f95-4936-bab9-30551184373f)
EnhancerAtlas 3.0 provides enhancer annotation in 26 species, including arabidopsis (TAIR10), worm (ce11), cow(bosTau9), Chinese hamster(CHOK1GS_HDv1), vervet(ChlSab1), fly (dm6), zebrafish (danRer11),E. sp. Ec32(ASM31002v1), human (hg38), chicken (galGal6), mouse (mm10), K. lactis(ASM251v1), rhesus (rheMac10), chimp (panTro5), rice (IRGSP1), rat(rn7), P. falciparum(ASM276v2), S. pombe(ASM294v2), yeast(sacCer3), Tomato(SL3), S. purpuratus(Spur5), pig(susScr11), T. gondii(TGA4), T. brucei brucei(TryBru), frog(Xtro10), and maize (B73v4). The consensus enhancers were predicted based on multiple high throughput experimental datasets. Currently, EnhancerAtlas 3.0 comprises a total of 57,264,271 annotated consensus enhancers derived from 38,504 datasets that encompass 26 species and 1,428 distinct tissues and cell types. The datasets incorporate information from 22 high-throughput sequencing technologies (e.g. CUT&TAG, CUT&RUN, HiChIP, NET-seq, PRO-seq, ChIP-seq, DHS, STARR-seq, MPRA). Furthermore, we predicted enhancer-target interactions for human, mouse, and fly, with 1,642,747, 3,114,890, and 132,119 interactions identified, respectively. 

# Dependencies
To conduct prediction of typical enhancers, and enhancer-gene interactions, following softwares are required:<br />
✯ ROSE<br />
✯ R version > 4.2<br />
✯ Signac 1.12.0<br />
✯ Seurat 5.1.0<br />
✯ cicero 1.3.9<br />
✯ Perl >5.16.3<br />
✯ bedtools > 2<br />
✯ macs2 2.2.6<br />
✯ sinto 0.10.1<br />
✯ deeptools > 2<br />
✯ python3<br />

# Tutorial
**[Prediction-of-enhancers](https://github.com/EvansGao/dbscATAC/wiki/Prediction-of-enhancers)**   : Predict consensus enhancers by the GSE149683 dataset for mouse<br />
**[Prediction-of-enhancer-gene-interactions](https://github.com/EvansGao/dbscATAC/wiki/Identification-of-markers)**   : Identify enhancer-gene interactions using the GSE163697 dataset for fly<br />

  
# Download of results
All the single-cell chromatin accessibility annotation data are available at:<br />
http://singlecelldb.com/EnhancerAtlas/downloadv2.php
