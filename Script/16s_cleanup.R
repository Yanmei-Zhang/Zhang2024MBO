# This is code to replicate the analyses and pictures from my DNA extraction method paper 
# *Code developed by Yanmei zhang*

# 16s amplicon cleanup in Cloquet DNA optimization trial
  
## Load packages and data
ps.16s = readRDS("../03_tophyloseq/results/ps.16s.rds")
sample_data(ps.16s)$Protocol <-factor(sample_data(ps.16s)$Protocol, level = c("QIAGEN", "CTAB"))

## Statistic of raw reads befpre processing
depth_track <- data.frame(sample_data(ps.16s)) %>%
  mutate(raw = sample_sums(ps.16s))  %>%
  mutate(SampleID = row.names(.))

depth_track %>%
  group_by(Sample_or_Control) %>%
  dplyr::summarise(range(raw), mean(raw), sd(raw))

depth_track %>%
  filter(Sample_or_Control == "Sample") %>%
  mutate(status = ifelse(Year_since_decay_start %in% c(0), "non_decayed", "decayed")) %>%
  group_by(Tree_host, status) %>%
  dplyr::summarise(mean(raw), sd(raw))

## Taxonomic based filtering
###The non-bacterial taxa need to be removed.
ntaxa(subset_taxa(ps.16s, Order == "Chloroplast")) #163
ntaxa(subset_taxa(ps.16s, Family == "Mitochondria")) #430

### Remove non-bacteria taxa
ps.bacteria = subset_taxa(ps.16s, (Family !="Mitochondria"|is.na(Family)) & (Order !=
                                                                               "Chloroplast"|is.na(Order))) # 593 ASVs (163 +430) were removed

ps.bacteria  = subset_taxa(ps.bacteria, ! is.na(Kingdom)) # there was 9 asv NOT in Kingdom


### Pick up the NA phylum (132)
ps.na = subset_taxa(ps.bacteria, is.na(Phylum)) 
### Recover DNA sequence and save in fasta format
rep.seq.na <- refseq(ps.na)
### Blast the NA sequence to pick put the host sequence
### The NA sequence are bacteria sequence

### Make a data frame with a column for the read counts of each sample
depth_track1 <- data.frame(sample_data(ps.bacteria)) %>%
  mutate(nonbacterial_filtered = sample_sums(ps.bacteria)) %>%
  right_join(depth_track)

### Statistical summary of ASV classified as choloroplast (Order = "Chloroplast" )
length(get_taxa_unique(ps.16s,taxonomic.rank="Order")) #222

Taxonomies_order <- ps.16s %>%
  tax_glom_wt(ranks = "Order") %>% 
  transform_sample_counts(function(x) {x/sum(x)} )%>% 
  psmelt() %>%
  #filter(Abundance > 0.05) %>%
  arrange(Order)

Order_stat <- Taxonomies_order %>%
  filter(Sample_or_Control == "Sample") %>%
  mutate(status = ifelse(Year_since_decay_start %in% c(0), "non_decayed", "decayed")) %>%
  group_by(Tree_host, status, Order) %>%
  dplyr::summarise(mean=mean(Abundance)*100, sd=sd(Abundance)) %>%
  arrange(dplyr::desc(mean)) %>%
  as.data.frame()

# Statistical summary of ASV classified as mitochondria (Family = "Mitochondria")
length(get_taxa_unique(ps.16s,taxonomic.rank="Family")) #332

Taxonomies_family <- ps.16s %>%
  tax_glom_wt(ranks = "Family") %>% 
  transform_sample_counts(function(x) {x/sum(x)} )%>% 
  psmelt() %>%
  #filter(Abundance > 0.05) %>%
  arrange(Family)

Family_stat <- Taxonomies_family %>%
  filter(Sample_or_Control == "Sample") %>%
  mutate(status = ifelse(Year_since_decay_start %in% c(0), "non_decayed", "decayed")) %>%
  group_by(Tree_host, status, Family) %>%
  dplyr::summarise(mean=mean(Abundance)*100, sd=sd(Abundance)) %>%
  arrange(dplyr::desc(mean)) %>%
  as.data.frame()


## Identifying and remove contaminants in marker-gene data
### Here is a link <https://benjjneb.github.io/decontam/vignettes/decontam_intro.html> to decotam instruction.
### Inspect Library Sizes
depth_df <- depth_track1[order(depth_track1$nonbacterial_filtered),]
depth_df$Index <- seq(nrow(depth_df))
ggplot(data=depth_df, aes(x=Index, y=nonbacterial_filtered, color=Sample_or_Control)) + 
  geom_point()

### Identify Contaminants - Prevalence
sample_data(ps.bacteria)$is.neg <- sample_data(ps.bacteria)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(ps.bacteria, method="prevalence", neg="is.neg")
head(contamdf.prev)

table(contamdf.prev$contaminant)

head(which(contamdf.prev$contaminant))

### More aggressive classification threshold 
contamdf.prev05 <- isContaminant(ps.bacteria, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
head(which(contamdf.prev05$contaminant))

### Make phyloseq object of presence-absence in negative controls and true samples
ps.bacteria.pa <- transform_sample_counts(ps.bacteria, function(abund) 1*(abund>0))
ps.bacteria.pa.neg <- prune_samples(sample_data(ps.bacteria.pa)$Sample_or_Control == "Control", ps.bacteria.pa)
ps.bacteria.pa.pos <- prune_samples(sample_data(ps.bacteria.pa)$Sample_or_Control == "Sample", ps.bacteria.pa)
### Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.bacteria.pa.pos), pa.neg=taxa_sums(ps.bacteria.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

### Remove contaminants
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps.bacteria)

### Remove negative control
ps.noncontam <- subset_samples(ps.noncontam, ! sample_names(ps.noncontam) %in% c("C0", "Q0"))
ps.noncontam <- prune_taxa(taxa_sums(ps.noncontam)>0, ps.noncontam)

### Make a data frame with a column for the read counts of each sample
depth_track2 <- data.frame(sample_data(ps.noncontam)) %>%
  mutate(decontam = sample_sums(ps.noncontam)) %>%
  right_join(depth_track1)


## Quality control 
### I will remove microbes that do not show up regularly. These are noisy OTUs, and while they may be real, they most likely popped up through contamination, sequencing and PCR error. They make the downstream stats noisy, and so I will remove them. In this dataset, I will remove OTUs that represented by less than ten reads in a given sample.
ps1 <- subset_samples(ps.noncontam, sample_sums(ps.noncontam)>=100) # 5 samples removed
ps1 = prune_taxa(taxa_sums(ps1)>0, ps1)

ps2 = filter_taxa(ps1, function(x) sum(x >= 10) >= 1, TRUE)

### Make a data frame with a column for the read counts of each sample
### Make a data frame with a column for the read counts of each sample
depth_track3 <- data.frame(sample_data(ps1)) %>%
  mutate(quality_filtered1 = sample_sums(ps1)) %>%
  right_join(depth_track2)

depth_track3 <- data.frame(sample_data(ps2)) %>%
  mutate(quality_filtered2 = sample_sums(ps2)) %>%
  right_join(depth_track3)