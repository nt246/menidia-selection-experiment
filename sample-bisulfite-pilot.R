library(tidyverse)


dat <- read_tsv("data/D_PolygenicDosage_beagle_Top0.99DiffStat_2.txt")
dat_0_5 <- dat %>% 
  filter(!str_detect(Population, "Gen10"))

dim(dat)
dim(dat_0_5)

ggplot(dat_0_5) +
  geom_histogram(aes(IncrSums, fill = Population), color = "black") +
  theme_classic()


ggplot(dat_0_5) +
  geom_histogram(aes(STD_Length, fill = Population), color = "black")


ggplot(dat_0_5) +
  geom_point(aes(IncrSums, STD_Length, color = Population)) +
  theme_classic()



# Look at how selected individuals fall on distribution

bisulfite <- read_tsv("data/bisulfite_inds.txt")
bisulfite$SampleID

ggplot(dat_0_5) +
  geom_point(aes(IncrSums, STD_Length, color = Population)) +
  geom_point(data = filter(dat_0_5, SampleID %in% bisulfite$SampleID), aes(IncrSums, STD_Length, color = Population), pch = 8, size = 8) +
  theme_classic()

ggsave(file = "M_menidia_for_bisulfite.png")


# Select individuals with intermediate values

dat_0_5 %>% 
  group_by(Population) %>% 
  summarize(mean_length = mean(STD_Length), mean_polyg = mean(IncrSums))


write_csv(dat_0_5, file = "PolygenicDosage_and_Length.csv")



dat_0_5 %>% 
  filter(str_detect(Population, "D1"), STD_Length > 48, STD_Length < 52, IncrSums > 28200, IncrSums < 28800)

dat_0_5 %>% 
  filter(str_detect(Population, "D2"), STD_Length > 48, STD_Length < 52, IncrSums > 28100, IncrSums < 28900)

dat_0_5 %>% 
  filter(SampleID %in% SS, STD_Length > 48, STD_Length < 52, IncrSums > 28100, IncrSums < 28900)

dat_0_5 %>% 
  filter(SampleID %in% SS)


dat_0_5 %>% 
  filter(str_detect(Population, "RGen0"), STD_Length > 70, STD_Length < 75, IncrSums > 23000, IncrSums < 24000)


dat_0_5 %>% 
  filter(str_detect(Population, "U1"), STD_Length > 68, STD_Length < 72, IncrSums > 17500, IncrSums < 18500)


dat_0_5 %>% 
  filter(str_detect(Population, "U2"), STD_Length > 68, STD_Length < 72, IncrSums > 17700, IncrSums < 18700)



# inport inversion data

inv_dat <- read_tsv("data/D2Gen5_MinusMixUps_Chr24_signifFDRcorrPfst_Genotype.txt")
table(inv_dat$SNP)

dim(inv_dat)

inv_dat %>% 
  select(starts_with("D2Gen5")) %>% 
  summarize(across(everything(), ~ sum(.x, na.rm = TRUE))) %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)  %>% 
  rename("allele_count" = `1`) %>% 
  ggplot() +
  geom_histogram(aes(allele_count))



SS <- inv_dat %>% 
  select(starts_with("D2Gen5")) %>% 
  summarize(across(everything(), ~ sum(.x, na.rm = TRUE))) %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)  %>% 
  rename("allele_count" = `1`) %>% 
  filter(str_detect(name, "D2")) %>% 
  arrange(allele_count) %>% 
  filter(allele_count > 9000) %>% 
  .$name



  View()
