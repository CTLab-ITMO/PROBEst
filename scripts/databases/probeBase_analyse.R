library(tidyverse)
library(Biostrings)
#BiocManager::install("pwalign")

df1 <- read.csv("./data/databases/open/probeBase_formatted.csv") %>%
  mutate(id = as.numeric(id))
colnames(df1) <- colnames(df1) %>% str_remove_all("\\.")

# corrections
df1$Formamide <- df1$Formamide %>% as.numeric()
df1$Formamide[is.na(df1$Formamide)] <- 0
df1$Formamide %>% lattice::densityplot()
table(df1$Formamide) %>% sapply(function(x) x/2640) %>% round(4) %>% as.data.frame()

df1$GCcontent %>% lattice::densityplot()
df1$Lengthnt %>% lattice::densityplot()

df2 <- read.csv("./data/databases/open/probeBase_genome_results.csv")

df <- full_join(df1, df2, by = "id") %>%
  mutate(sseq = Sequence %>% str_remove_all("[^ATGCUI]")) %>%
  mutate(qseq = paste0(left_flank, target, right_flank)) %>%
  select(id, sseq, qseq,left_flank,right_flank,target,
         `Formamide`,`GCcontent`,`Hybridizationefficiency`,
         `Lengthnt`,`Modifiedversions`,
         evalue, mismatches, length, bitscore,identity) %>%
  mutate(
    qseq_aln = NA,
    sseq_aln = NA,
    score = NA,
    tmp = NA,
    G = NA
  )

for (i in 1:nrow(df)){
  if(qseq != "NANANA") {
    align = pwalign::pairwiseAlignment(df$sseq[i], df$qseq[i], type = "overlap")
    df$qseq_aln[i] <- substr(
      df$qseq[i],
      align@subject@range@start,
      align@subject@range@start+align@pattern@range@width-1)

    df$sseq_aln[i] <- substr(
      df$sseq[i],
      align@pattern@range@start,
      align@pattern@range@start+align@pattern@range@width-1)

    df$score[i] <- align@score
  }
}


