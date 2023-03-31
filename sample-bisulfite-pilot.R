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
  theme_classic() +
  labs(x = "Polygenic score (count of slow-growing alleles)", y = "Length (mm)")

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



# import inversion data
# I got this mixed up. It's the northern homozygote that's rare. So SS with high counts should be NN

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

ggsave(file = "Inversion_allele_count_histogram.png")


SS <- inv_dat %>% 
  select(starts_with("D2Gen5")) %>% 
  summarize(across(everything(), ~ sum(.x, na.rm = TRUE))) %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)  %>% 
  rename("allele_count" = `1`) %>% 
  filter(str_detect(name, "D2")) %>% 
  arrange(allele_count) %>% 
  filter(allele_count > 5000) %>% 
  .$name


inv_allelecounts <- inv_dat %>% 
  select(starts_with("D2Gen5")) %>% 
  summarize(across(everything(), ~ sum(.x, na.rm = TRUE))) %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)  %>% 
  rename("inv_allele_count" = `1`) %>% 
  filter(str_detect(name, "D2")) %>% 
  arrange(name) 





# Sample selection round 2 (March 2023) -----------------------------------

# Need the five largest and five smallest fish from each of the groups D1, D2, U1, U2

sample_selection <- vector("character", length = 40)

pops <- c("D1Gen5", "D2Gen5", "U1Gen5", "U2Gen5")

p=3 #didn't write out the for loop, just did it manually
for (p in seq_along(pops)) 
  
  dat_sort <- dat_0_5 %>% 
    filter(Population == pops[p]) %>% 
    arrange(STD_Length) 
  
  sample_selection[31:35] <- dat_sort %>% 
    head(5) %>% 
    pull(SampleID)
  
  sample_selection[36:40] <- dat_sort %>% 
    tail(5) %>% 
    pull(SampleID)
    
write_lines(sample_selection, file = "data/bisulfite_inds_v2_temp.txt")



# Visualize samples on the plot
bisulfite$SampleID

ggplot(dat_0_5) +
  geom_point(aes(IncrSums, STD_Length, color = Population)) +
  geom_point(data = filter(dat_0_5, SampleID %in% sample_selection), aes(IncrSums, STD_Length, color = Population), pch = 8, size = 8) +
  theme_classic() +
  labs(x = "Polygenic score (count of slow-growing alleles)", y = "Length (mm)")

ggsave(file = "M_menidia_for_bisulfite_size_extremes.png")



ggplot(dat_0_5) +
  geom_point(aes(IncrSums, STD_Length, color = Population)) +
  geom_point(data = filter(dat_0_5, SampleID %in% c(sample_selection, bisulfite$SampleID)), aes(IncrSums, STD_Length, color = Population), pch = 8, size = 8) +
  theme_classic() +
  labs(x = "Polygenic score (count of slow-growing alleles)", y = "Length (mm)")




## Examine the inversion type of D2 inds

write_csv(filter(dat_0_5, SampleID %in% sample_selection), file = "data/M_menidia_for_bisulfite_size_extremes.csv")

dat_0_5 %>% 
  left_join(inv_allelecounts, by = c("SampleID" = "name")) %>% 
  filter(Population == "D2Gen5") %>% 
  mutate(inv_type = ifelse(inv_allele_count < 5000, "SS", ifelse(inv_allele_count >9000, "NN", "NS")),
         IncrSums = round(IncrSums, 1),
         inv_allele_count = round(inv_allele_count, 1)) %>% 
  arrange(STD_Length) %>% 
  write_csv(file = "data/D2Gen5_data.csv")


# Visualize
manual_round2 <- read_lines("data/Adjusted_round2_selection.txt")
length(manual_round2)

bisulfite$SampleID

ggplot(dat_0_5) +
  geom_point(aes(IncrSums, STD_Length, color = Population)) +
  geom_point(data = filter(dat_0_5, SampleID %in% c(bisulfite$SampleID, manual_round2)), aes(IncrSums, STD_Length, color = Population), pch = 8, size = 8) +
  theme_classic() +
  labs(x = "Polygenic score (count of slow-growing alleles)", y = "Length (mm)")

ggsave(file = "Menidia_bisulfite_round2_size_extremes.png")

