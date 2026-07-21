library(tidyverse)
library(plyranges)

inGFF <- read_gff(file = "../../0_data/reference/plasmoDB52/PlasmoDB-52_Pfalciparum3D7.gff")

inGFF %>% 
  filter(type %in% c("protein_coding_gene", "pseudogene", "ncRNA_gene")) %>% 
  mutate(tss = ifelse(strand == "+", start,end),
         tts = ifelse(strand == "+", end,start),
		 name = ID) %>% 
  write_gff3(file = "../../0_data/reference/pf_nucDyn.gff")
