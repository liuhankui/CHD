# CHD

# required software
```
fastp, https://github.com/OpenGene/fastp
BWA, https://github.com/lh3/bwa
samblaster, https://github.com/GregoryFaust/samblaster
sambamba, https://github.com/biod/sambamba
Samtools, https://github.com/samtools/samtools
GATK4, https://gatk.Broadinst-itute.org
vawk, https://github.com/cc2qe/vawk
bcftools, https://github.com/samtools/bcftools
tabix,bgzip, https://github.com/tabixio/tabix
mosdepth, https://github.com/brentp/mosdepth
vep, https://github.com/Ensembl/ensembl-vep
```

# fq2gvcf
```
fastp -q 20 -u 10 -n 5 --in1 fq1.gz --in2 fq2.gz --out1 clean.fq1.gz --out2 clean.fq2.gz
bwa mem -t 4 -R '@RG\tID:id\tPL:illumina\tPU:id\tLB:sample\tSM:id\tCN:BGI' GRCh38.fa clean.fq1.gz clean.fq2.gz|samblaster --excludeDups --ignoreUnmated --maxSplitCount 2 --minNonOverlap 20 |samtools view -Sb -|samtools sort - -O CRAM -o sort.cram --reference GRCh38.fa
mosdepth qc sort.cram -f GRCh38.fa --fast-mode --no-per-base --by V5.bed --thresholds 1,4,10,20
java -Xmx3g -jar gatk-package-4.1.2.0-local.jar HaplotypeCaller --QUIET true -R GRCh38.fa -L V5.bed -I sort.cram -O gatk.gvcf.gz -ERC GVCF -A ClippingRankSumTest -A LikelihoodRankSumTest -A MappingQualityZero
```

# gvcf2trio.vcf
```
java -Xmx3g -jar gatk-package-4.1.2.0-local.jar CombineGVCFs -R GRCh38.fa -L V5.bed -V child.gvcf.gz -V father.gvcf.gz -V mother.gvcf.gz  -O trio.gvcf.gz
java -Xmx3g -jar gatk-package-4.1.2.0-local.jar GenotypeGVCFs --QUIET true -R GRCh38.fa -L V5.bed -V trio.gvcf.gz -O trio.vcf.gz
java -Xmx3g -jar gatk-package-4.1.2.0-local.jar VariantRecalibrator -RGRCh38.fa -V trio.vcf.gz --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_146.hg38.vcf.gz -an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -O SNP.recal --tranches-file SNP.tranches --max-gaussians 4
java -Xmx3g -jar gatk-package-4.1.2.0-local.jar ApplyVQSR -R GRCh38.fa -V trio.vcf.gz -O trio.VQSR.SNP.vcf.gz --truth-sensitivity-filter-level 99.0 --tranches-file SNP.tranches --recal-file SNP.recal -mode SNP
java -Xmx3g -jar gatk-package-4.1.2.0-local.jar VariantRecalibrator -R GRCh38.fa -V trio.SNP.vcf.gz --resource:mills,known=true,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode INDEL -O INDEL.recal --tranches-file INDEL.tranches --max-gaussians 2
java -Xmx3g -jar gatk-package-4.1.2.0-local.jar ApplyVQSR -R GRCh38.fa -V trio.VQSR.SNP.vcf.gz -O trio.VQSR.SNP.INDEL.vcf.gz --truth-sensitivity-filter-level 99.0 --tranches-file INDEL.tranches --recal-file INDEL.recal -mode INDEL
```

# vcf2dnv
```
python3.6 pydnm -v trio.VQSR.SNP.INDEL.vcf.gz -f trio.fam -o trio.DNV
cat trio.DNV.txt|awk '{if(NR>1){if(length($4)==length($5)){if($NF>0.5)print $1,$2,$2,$4,$5,$6,$NF}else{if($NF>0.5 && length($4)+length($5)<22)print$1,$2,$2,$4,$5,$6,$NF}}}'|tr ' ' '\t' > s.txt
bcftools view trio.VQSR.SNP.INDEL.vcf.gz -R s.txt|awk '{if(NR==1){print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC"};print $0}'|vawk --header '{split(S$A$AD,a,",");if(S$B$DP>9 && S$C$DP>9 && S$B$AD~/,0/ && S$C$AD~/,0/ && a[1] >=8 && a[2] >= 8 && (I$Prob > 0.9 || I$MQ==60))print}'|gzip -f > DNV.vcf.gz
vep --assembly GRCh38 --fork 4 -i DNV.vcf.gz -o DNV.vep.vcf.gz --vcf --compress_output bgzip --no_stats  --merged --force_overwrite --offline --use_given_ref \
--total_length --numbers --ccds --hgvs --symbol --canonical --protein --biotype --tsl --nearest symbol \
--fasta GRCh38.fa \
--dir_cache vep_path/cache \
--dir_plugins vep_path/VEP_plugins \
--custom gnomAD.vcf.gz,gnomADg,vcf,exact,0,AF,AF_eas,AF_popmax,popmax,AF_XX,AF_XY,AF_eas_XX,AF_eas_XY,nhomalt,nhomalt_XX,nhomalt_XY,nhomalt_eas,nhomalt_eas_XX,nhomalt_eas_XY,nhomalt_popmax \
--custom mbiobank_ChinaMAP.phase1.vcf.gz,ChinaMAP,vcf,exact,0,AF \
--plugin CADD,vep_path/GRCh38/whole_genome_SNVs.tsv.gz,vep_path/InDels.tsv.gz \
--plugin REVEL,vep_path/new_tabbed_revel_grch38.tsv.gz \
--plugin PrimateAI,vep_path/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz \
--plugin SpliceAI,snv=vep_path/spliceai_scores.raw.snv.hg38.vcf.gz,indel=vep_path/spliceai_scores.raw.indel.hg38.vcf.gz \
--plugin dbNSFP,vep_path/dbNSFP4.3a.gz,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,M-CAP_score,M-CAP_pred,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id \
--plugin DisGeNET,file=vep_path/DisGeNET/all_variant_disease_pmid_associations_final.tsv.gz,disease=1
bcftools +split-vep DNV.vep.vcf.gz -l|awk '{print $2}'|tr '\n' ' '|awk '{print "CHR POS REF ALT FAM Prob MQ GT AD GT AD GT AD "$0}'|sed 's/ $//'|tr ' ' '\t' > DNV.xls
bcftools +split-vep DNV.vep.vcf.gz -f '%CHROM %POS %REF %ALT %INFO/ID %INFO/Prob %INFO/MQ [%GT %AD ]%CSQ\n' -d -A '\t'|tr '/' '|'|tr ' ' '\t' >> DNV.xls
cat DNV.xls|awk 'NR==1 || (($73=="." || $73<=0.001) && ($74=="." || $74<=0.001) && ($89=="." || $89<=0.001) && $17!="." && $37=="YES"){split($15,a,"&");$15=a[1];print}'|tr ' ' '\t' > DNV.rare.txt
```

# recc burden
```
cumAF<-function(x){
  af<-0
  for(i in x){
    if(is.na(i)){i<-0}
    af<-af+(1-af)*i
  }
  return(af)
}

cadd<-25.3
maf<-0.001
df0<-read.table(gzfile('gnomad.all.gz'))
names(df0)<-c('V5','V6','V7','V9','V10','V11','AF')
df0$V5<-apply(df0,1,function(x){unlist(strsplit(x[1],'&'))[1]})
df0<-df0[!(df0$V6=='MODERATE' & df0$V5!='missense_variant'),]
df0<-df0[df0$V7!='.' & df0$AF!='.' & df0$AF>0,]
df0<-df0[df0$AF<maf,]
df0$V9[df0$V9=='.']<-0
df0$V11[df0$V11=='.']<-0
df0$V9<-as.numeric(df0$V9)
df0$V11<-as.numeric(df0$V11)

df0<-df0[df0$V6=='HIGH' | df0$V5=='synonymous_variant' | df0$V5=='missense_variant',]
df0$V6[df0$V6=='MODERATE' & df0$V5=='missense_variant' & df0$V9>=cadd]<-'HIGH'

df0$AF<-as.numeric(df0$AF)
adf<-ddply(df0,.(V6,V7,V10),summarise,af=cumAF(AF))
adf$af<-adf$af*adf$af
adf$tag<-paste(adf$V6,adf$V7,adf$V10,sep='_')
ratedf<-adf[,1:4]
names(ratedf)<-c('Model','Gene','Transcript','E_Number')

df<-read.table(gzfile('CHD.rare.gz'),sep='\t')
df$ll<-nchar(df[,3])+nchar(df[,4])
df<-df[df$ll<12,]
df<-df[df$V7!='.',]
df<-df[df[,34]/df[,35]<0.005 & df[,35]/9066>0.8,]#
df$V5<-apply(df,1,function(x){unlist(strsplit(x[5],'&'))[1]})
df<-df[!(df$V6=='MODERATE' & df$V5!='missense_variant'),]

df$V9[df$V9=='.']<-0
df$V11[df$V11=='.']<-0
df$V9<-as.numeric(df$V9)
df$V11<-as.numeric(df$V11)
df<-df[df$V6=='HIGH' | df$V5=='synonymous_variant' | df$V5=='missense_variant',]
df$V6[df$V6=='MODERATE' & df$V5=='missense_variant' & df$V9>=cadd]<-'HIGH'

df1<-df[df$V44=='1/1'& df$V45=='0/1' & df$V46=='0/1',]
sdf1<-ddply(df1,.(V6,V7,V10,V42,V43),summarise,num=length(unique(V42)))

df2<-df[df$V44!='1/1',]
df2$source<-ifelse(df2$V45=='0/1','f','m')
sdf2<-ddply(df2,.(V6,V7,V10,V42,V43),summarise,num=length(unique(source)))
sdf2<-sdf2[sdf2$num==2,]

sdf<-rbind(sdf1,sdf2)

sdf<-ddply(sdf,.(V6,V7,V10,V43),summarise,num=length(unique(V42)))
odf<-dcast(sdf[sdf$V7!='.',],V6+V7+V10~V43)
odf[is.na(odf)]<-0

odf$tag<-paste(odf$V6,odf$V7,odf$V10,sep='_')
ndf<-merge(ndf,adf[,c(4,5)],by='tag',all.x=T)
ndf$af[is.na(ndf$af)]<-1/100000
ndf$p<-1- ppois(ndf$Case-1,ndf$af*820)
ndf$cutoff<-qpois(1-0.05/20000,ndf$af*820)
ndf<-ndf[order(ndf$af),]
ndf$V6<-factor(ndf$V6,levels = c('HIGH','MODERATE','LOW'),labels=c('LoF+D-Missense','T-Missense','Synonymous'),order=T)

tmpdf<-ndf[,c(2,3,5,8,9,10)]
names(tmpdf)<-c("Model","Gene","O_high_case","E_Number","Pvalue","cutoff")
tmpdf$symbol<-tmpdf$Gene
tmpdf$E_num<-tmpdf$E_Number*820

tmpdf$sign<-ifelse(tmpdf$Pvalue<0.05/20000,'Genome-wide','No evidence')
tmpdf$sign<-factor(tmpdf$sign,levels=c('Genome-wide','No evidence'),order=T)
tmpldf<-tmpdf[tmpdf$Pvalue<0.05/20000, ]
set.seed(2023)
ggplot(tmpdf[tmpdf$O_high_case>0 & tmpdf$Model!='T-Missense' & tmpdf$Pvalue>0.05/20000,],aes(E_num/487,O_high_case))+
  geom_path(aes(E_num/820,cutoff,group=Model),colour='black')+
  geom_point(shape=21,colour='black',aes(fill=sign,size=-log10(Pvalue)))+
  theme_classic()+
  facet_wrap(~Model,nrow=1,scale='free')+
  scale_size_continuous(expression(paste('-',log[10],"P-value")),breaks = c(1,2,3,4,5,6,7),range = c(0.1,3)) +
  scale_fill_manual("Significance",values=c('steelblue'))+
  scale_x_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^(x)),
    labels = trans_format("log10", math_format(10^.x))
  )+
  #scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  xlab('Expected RG rate')+ylab('Observed RG')+
  theme(#legend.position='none',
    axis.text=element_text(size=12,colour="black"),
    axis.title = element_text(size=15,colour="black"),
    strip.text = element_text(size=12,colour="black"),
    plot.title = element_text(size = 15))
```
