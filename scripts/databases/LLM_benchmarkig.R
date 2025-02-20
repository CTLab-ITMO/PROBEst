library(jsonlite)
library(tidyverse)

indir <- "../../data/databases/articles/benchmark_json"
postf <- "pitz2021.json"

fnames <- dir(indir, pattern = postf, full.names = T)


#read
res <- list()
for (fname in fnames) {
  sname <- basename(fname) %>% str_remove_all("_.*")
  res[[sname]] <- read_json(fname)
}

#parse
pfun <- function(json, sname) {
  r1 <- length(unlist(json))
  r2 <- length(json[["article_data"]][["probe_samples"]][["sample_group"]][[1]][["probe_group"]][["sequences"]][["probe_sequences"]][["probes"]][["probe"]] %>% 
                 map("probe_sequence") %>% 
                 unlist)
  
  df <- data.frame(
    name = sname,
    length = r1,
    probes = r2
  )
}


res_df <- tibble()
for (sname in names(res)) {
  res_df <- rbind(res_df, pfun(res[[sname]], sname))
}

res_df %>% 
  arrange(-probes) %>% 
  mutate(name = fct_inorder(name)) %>% 
  pivot_longer(cols = -1, names_to = "param") %>% 
  ggplot(aes(name, value, fill = name)) +
  geom_col(show.legend = F) +
  facet_grid(param~., scales = "free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")
