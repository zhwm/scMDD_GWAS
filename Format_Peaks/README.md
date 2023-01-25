# Format scATAC-seq peaks for calculating LD scores

## Format differential expression peaks

### Get cell types list
```
sed "s/\"//g" top5kDAR_combined_subcluster.csv | cut -f 3 -d ',' | sed 1d | sort | uniq > top5kDAR.celltypes.lst
```

### Extract peak information for each cell type in hg38
```
while read a; do for i in top5kDAR_combined_subcluster.csv top5kDAR_Female.csv top5kDAR_Male.csv top5kDAR_interaction.csv; do sed "s/\"//g" $i | cut -f 2,3 -d ',' | awk -F ',' -v var=$a '$2==var{print}' | cut -f 1 -d ',' | awk -F '-' '{print $1"\t"$2"\t"$3"\t"$0}' > hg38/$i\_$a\.bed ; done; done < top5kDAR.celltypes.lst
```

### Extract all peaks as control set
```
for i in top5kDAR_combined_subcluster.csv top5kDAR_Female.csv top5kDAR_Male.csv top5kDAR_interaction.csv; do cut -f 2 $i -d ',' | sed 1d | sed "s/\"//g" | sort | uniq | sed "s/-/\t/g" > hg38/$i\.bed; done
```

## Format marker peaks
```
for i in MarkerSubPeaks/*.tsv; do cut -f 1,3,4 $i | sed 1d | awk '{print $0"\t"$1"-"$2"-"$3}' >  hg38/*.bed ; done
```

## Format all peaks
```
for i in Ast End ExN InN Mic Oli OPC; do Rscript get_peaks.R BroadPeaks_without_constraint/$i\-reproduciblePeaks.gr.rds 250 hg38/$i\.bed; done
while read a; do Rscript get_peaks.R SubClusterPeaks_without_constraint/$a\-reproduciblePeaks.gr.rds 250 hg38/$a\.bed; done < sub.lst 
```

## Convert peak coordinates to hg19
```
for i in hg38/*.bed; do ~/utils/liftOver $i ~/utils/hg38ToHg19.over.chain.gz $(echo $i | sed "s/38/19/g") unmapped ; done
cat hg19/*_markerpeaks.tsv.bed | sort | uniq > hg19/sub_markerpeaks.tsv.bed
for i in Ast End ExN InN Mic Oli OPC; do cat hg19/$i.bed ; done > hg19/broad.bed
while read a; do cat hg19/$a.bed ; done < sub.lst > hg19/sub.bed
```
