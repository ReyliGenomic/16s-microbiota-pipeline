library(vegan)

dist_bc <- vegdist(taxa_corrected2, method = "bray")

adonis2(
  dist_bc ~ batchid + agey + Sx + BMI + glu + SBP + chol + Class..SBP_CTRL,
  data = meta_df,
  permutations = 999,
  by = "margin"
)

library(vegan)
library(compositions)

# CLR transform
otu_clr <- clr(taxa_corrected2 + 1)

dist_ait <- dist(otu_clr)

adonis2(
  dist_ait ~ batchid + agey + Sx + BMI + glu + SBP + chol + Class..SBP_CTRL,
  data = meta_df,
  permutations = 999,
  by = "margin"
)

anova(betadisper(dist_ait, meta_df$Lot))

bd <- betadisper(dist_ait, meta_df$Lot)
plot(bd)

#Permanova sin batch 
adonis2(
  dist_ait ~ agey + Sx + BMI + glu + SBP + chol + Class..SBP_CTRL,
  data = meta_df,
  permutations = 999,
  by = "margin"
)


