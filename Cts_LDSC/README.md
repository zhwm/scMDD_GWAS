# Cell type specific analysis

## Cell type specific analysis for differential peaks
```
awk '{print$1"\t../Calculate_LD_Score/top5kDAR_combined_subcluster.csv_"$1".,../Calculate_LD_Score/top5kDAR_combined_subcluster.csv."}' ../Format_Peaks/top5kDAR.celltypes.lst > top5kDAR_combined_subcluster.csv.ldcts
for i in ../Format_GWAS/*.sumstats.gz; do python ~/utils/ldsc/ldsc.py --h2-cts $i --ref-ld-chr ../1KG/1000G_EUR_Phase3_baseline/baseline. --out $(echo $i | cut -d '/' -f 3 | sed "s/.sumstats.gz//g")_top5kDAR_combined_subcluster.csv  --ref-ld-chr-cts top5kDAR_combined_subcluster.csv.ldcts --w-ld-chr ../1KG/weights_hm3_no_hla/weights.; done
```

## Cell type specific analysis for female peaks  
```
awk '{print$1"\t../Calculate_LD_Score/top5kDAR_Female.csv_"$1".,../Calculate_LD_Score/top5kDAR_Female.csv."}' ../Format_Peaks/top5kDAR.celltypes.lst > top5kDAR_Female.csv.ldcts
for i in ../Format_GWAS/*.sumstats.gz; do python ~/utils/ldsc/ldsc.py --h2-cts $i --ref-ld-chr ../1KG/1000G_EUR_Phase3_baseline/baseline. --out $(echo $i | cut -d '/' -f 3 | sed "s/.sumstats.gz//g")_top5kDAR_Female.csv  --ref-ld-chr-cts top5kDAR_Female.csv.ldcts --w-ld-chr ../1KG/weights_hm3_no_hla/weights.; done
```

## Cell type specific analysis for male peaks
```
awk '{print$1"\t../Calculate_LD_Score/top5kDAR_Male.csv_"$1".,../Calculate_LD_Score/top5kDAR_Male.csv."}' ../Format_Peaks/top5kDAR.celltypes.lst > top5kDAR_Male.csv.ldcts
for i in ../Format_GWAS/*.sumstats.gz; do python ~/utils/ldsc/ldsc.py --h2-cts $i --ref-ld-chr ../1KG/1000G_EUR_Phase3_baseline/baseline. --out $(echo $i | cut -d '/' -f 3 | sed "s/.sumstats.gz//g")_top5kDAR_Male.csv  --ref-ld-chr-cts top5kDAR_Male.csv.ldcts --w-ld-chr ../1KG/weights_hm3_no_hla/weights.; done
```

## Cell type specific analysis for interaction peaks
```
awk '{print$1"\t../Calculate_LD_Score/top5kDAR_interaction.csv_"$1".,../Calculate_LD_Score/top5kDAR_interaction.csv."}' ../Format_Peaks/top5kDAR.celltypes.lst > top5kDAR_interaction.csv.ldcts
for i in ../Format_GWAS/*.sumstats.gz; do python ~/utils/ldsc/ldsc.py --h2-cts $i --ref-ld-chr ../1KG/1000G_EUR_Phase3_baseline/baseline. --out $(echo $i | cut -d '/' -f 3 | sed "s/.sumstats.gz//g")_top5kDAR_interaction.csv  --ref-ld-chr-cts top5kDAR_interaction.csv.ldcts --w-ld-chr ../1KG/weights_hm3_no_hla/weights.; done
```

## Cell type specific analysis for marker peaks
```
awk '{print$1"\t../Old_LD_Score/"$1"_markerpeaks.tsv.,../Old_LD_Score/sub_markerpeaks.tsv."}' ../Format_Peaks/sub_markerpeaks.tsv > sub_markerpeaks.tsv.ldcts
for i in ../Format_GWAS/*.sumstats.gz; do python ~/utils/ldsc/ldsc.py --h2-cts $i --ref-ld-chr ../1KG/1000G_EUR_Phase3_baseline/baseline. --out $(echo $i | cut -d '/' -f 3 | sed "s/.sumstats.gz//g")_sub_markerpeaks.tsv  --ref-ld-chr-cts sub_markerpeaks.tsv.ldcts --w-ld-chr ../1KG/weights_hm3_no_hla/weights.; done
```

## Cell type specific analysis for broad cell groups
```
awk '{print$1"\t../Old_LD_Score/"$1".,../Old_LD_Score/sub."}' ../Format_Peaks/sub.lst > sub_allpeaks.tsv.ldcts
for i in ../Format_GWAS/*.sumstats.gz; do python ~/utils/ldsc/ldsc.py --h2-cts $i --ref-ld-chr ../1KG/1000G_EUR_Phase3_baseline/baseline. --out $(echo $i | cut -d '/' -f 3 | sed "s/.sumstats.gz//g")_sub_allpeaks.tsv  --ref-ld-chr-cts sub_allpeaks.tsv.ldcts --w-ld-chr ../1KG/weights_hm3_no_hla/weights.; done
```

## Cell type specific analysis for subtypes
```  
awk '{print$1"\t../Old_LD_Score/"$1".,../Old_LD_Score/broad."}' ../Format_Peaks/broad.lst > broad_allpeaks.tsv.ldcts
for i in ../Format_GWAS/*.sumstats.gz; do python ~/utils/ldsc/ldsc.py --h2-cts $i --ref-ld-chr ../1KG/1000G_EUR_Phase3_baseline/baseline. --out $(echo $i | cut -d '/' -f 3 | sed "s/.sumstats.gz//g")_broad_allpeaks.tsv  --ref-ld-chr-cts broad_allpeaks.tsv.ldcts --w-ld-chr ../1KG/weights_hm3_no_hla/weights.; done
```
