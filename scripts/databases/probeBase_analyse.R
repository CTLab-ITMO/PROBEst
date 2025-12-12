library(tidyverse)
library(Biostrings)
#BiocManager::install("pwalign")
library(reticulate)

# Initialize Python environment and import the module
use_python("~/miniconda3/bin/python")
rna_structure <- import("PROBESt.rna_structure")

# Function to calculate hairpin probability using Python
calculate_hairpin_prob <- function(sequence) {
  return(rna_structure$calculate_hairpin_prob(sequence))
}

# Function to calculate dimer Gibbs energy using Python
calculate_dimer_G <- function(string1, string2, type1 = "RNA", type2 = "DNA") {
  return(rna_structure$calculate_dimer_G(string1, string2, type1, type2))
}

df_combine <- function(p1, p2) {
  df1 <- read.csv(p1) %>%
    mutate(id = as.numeric(id))
  colnames(df1) <- colnames(df1) %>% str_remove_all("\\.")

  # corrections
  df1$Formamide <- df1$Formamide %>% as.numeric()
  df1$Formamide[is.na(df1$Formamide)] <- 0
  df1$Formamide %>% lattice::densityplot()
  table(df1$Formamide) %>% sapply(function(x) x/2640) %>% round(4) %>% as.data.frame()

  df1$GCcontent %>% lattice::densityplot()
  df1$Lengthnt %>% lattice::densityplot()

  df2 <- read.csv(p2)

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
      hairpin_prob = NA,
      dimer_DNA = NA,
      dimer_DNA_flank = NA,
      dimer_probe = NA,
      dimer_probe_DNA = NA,
      tmp = NA,
      G = NA
    ) %>%
    subset(qseq != "NANANA")

  for (i in 1:nrow(df)){
    if(df$qseq[i] != "NANANA") {
      align = pwalign::pairwiseAlignment(df$sseq[i], df$qseq[i], type = "overlap")
      df$qseq_aln[i] <- substr(
        df$qseq[i],
        align@subject@range@start,
        align@subject@range@start+align@pattern@range@width-1)

      df$sseq_aln[i] <- substr(
        df$sseq[i],
        align@pattern@range@start,
        align@pattern@range@start+align@pattern@range@width-1)

      # Calculate hairpin probability using Python function
      df$hairpin_prob[i] <- calculate_hairpin_prob(df$sseq[i])

      # Calculate dimer energies using Python functions
      df$dimer_DNA[i] <- calculate_dimer_G(df$sseq[i], string2 = NULL, type1 = "DNA", type2 = "DNA")
      df$dimer_DNA_flank[i] <- calculate_dimer_G(paste0(df$left_flank[i], df$sseq[i], df$right_flank[i]),
                                                 string2 = NULL, type1 = "DNA", type2 = "DNA")
      df$dimer_probe[i] <- calculate_dimer_G(df$sseq[i], df$sseq[i], type1 = "RNA", type2 = "RNA")
      df$dimer_probe_DNA[i] <- calculate_dimer_G(df$sseq[i], df$qseq[i], type1 = "RNA", type2 = "DNA")

      df$score[i] <- align@score
    }
  }

  return(df)
}

p1 = "./data/databases/open/probeBase_formatted.csv"
p2 = "./data/databases/open/probeBase_genome_results.csv"

df_true <- df_combine(p1,p2) %>% mutate(type = TRUE)

p1 = "./data/databases/open/probeBase_false_formatted.csv"
p2 = "./data/databases/open/probeBase_false_genome_results.csv"
df_false <- df_combine(p1,p2) %>% mutate(type = FALSE)

df <- rbind(df_true, df_false) %>%
  select(!c(id,target,tmp,G,Hybridizationefficiency))

df %>%
  write_csv("./data/databases/open/test_ML_database.csv")

df %>%
  ggplot(aes(color = type, group = type, fill = type)) +
  geom_point(aes(dimer_DNA_flank, dimer_probe_DNA), alpha = .5) +
  theme_minimal()

