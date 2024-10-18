# This is code to replicate the analyses and pictures from my DNA extraction method paper 
# *Code developed by Yanmei zhang*

# ITS2 amplicon cleanup in Cloquet DNA optimization trial"

## Load packages and data
ps.its = readRDS("ps.its.rds")
sample_data(ps.its)$Protocol <-factor(sample_data(ps.its)$Protocol, level = c("QIAGEN",
                                                                              "CTAB"))
## Statistic of raw reads befpre processing
depth_track <- data.frame(sample_data(ps.its)) %>%
  mutate(raw = sample_sums(ps.its)) %>%
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
### The non-fungal taxa need to be removed.
unique(tax_table(ps.its)[, "Kingdom"])
length(unique(tax_table(ps.its)[, "Kingdom"])) # 11

### Pick up the NA kingdom (111)
ps.k.na = subset_taxa(ps.its, is.na(Kingdom)) # 360
### Recover DNA sequence and save in fasta format
rep.seq.k.na <- refseq(ps.k.na)
### Blast the NA sequence and I found some of them were fungi and the question is how to filter them based on a criterion and add them to the taxa table? I will remove those NA sequences as They're not so much

### Remove non-fungal taxa
ps.fungi = subset_taxa(ps.its, Kingdom == "k__Fungi") # there was 2456 fungal taxa

ps.nonfungi = subset_taxa(ps.its, ! Kingdom == "k__Fungi"| is.na(Kingdom)) # there was 2399 nonfungal taxa

### Pick up the NA phylum (51)
ps.na = subset_taxa(ps.fungi, is.na(Phylum)) 
### Recover DNA sequence and save in fasta format
rep.seq.na <- refseq(ps.na)
### Blast the NA sequence to pick put the host sequence. The NA sequence are fungal sequence

### Make a data frame with a column for the read counts of each sample
depth_track1 <- data.frame(sample_data(ps.fungi)) %>%
  mutate(nonfungal_filtered = sample_sums(ps.fungi)) %>%
  right_join(depth_track)

### Statistical summary of ASV classified as non fungi (kingdom = "k__Fungi" )
length(get_taxa_unique(ps.its,taxonomic.rank="Kingdom"))

Kingdoms <- ps.its %>%
  aggregate_taxa("Kingdom") %>% 
  transform_sample_counts(function(x) {x/sum(x)}) %>% 
  psmelt() %>%
  arrange(Kingdom)
King_stat <- Kingdoms %>%
  dplyr::group_by(Tree_host,  Year_since_decay_start, Kingdom) %>%
  dplyr::summarise(mean=mean(Abundance)*100, sd=sd(Abundance)) %>%
  arrange(dplyr::desc(mean)) %>%
  as.data.frame()

## Identifying and remove contaminants in marker-gene data
Here is a link <https://benjjneb.github.io/decontam/vignettes/decontam_intro.html> to decotam instruction.
### Inspect Library Sizes
depth_df <- depth_track1[order(depth_track1$nonfungal_filtered),]
depth_df$Index <- seq(nrow(depth_df))
ggplot(data=depth_df, aes(x=Index, y=nonfungal_filtered, color=Sample_or_Control)) + 
  geom_point()

### Identify Contaminants - Prevalence
sample_data(ps.fungi)$is.neg <- sample_data(ps.fungi)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(ps.fungi, method="prevalence", neg="is.neg")
head(contamdf.prev)

table(contamdf.prev$contaminant)

head(which(contamdf.prev$contaminant))

### More aggressive classification threshold 
contamdf.prev05 <- isContaminant(ps.fungi, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
head(which(contamdf.prev05$contaminant))

### Remove contaminants
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps.fungi)

### Remove negative control
ps.noncontam <- subset_samples(ps.noncontam, ! sample_names(ps.noncontam) %in% c("C0", "Q0"))
ps.noncontam <- prune_taxa(taxa_sums(ps.noncontam)>0, ps.noncontam)

### Make a data frame with a column for the read counts of each sample
depth_track2 <- data.frame(sample_data(ps.noncontam)) %>%
  mutate(decontam = sample_sums(ps.noncontam)) %>%
  right_join(depth_track1)

## Quality control 
### I will remove microbes that do not show up regularly. These are noisy OTUs, and while they may be real, they most likely popped up through contamination, sequencing and PCR error. They make the downstream stats noisy, and so I will remove them. In this dataset, I will remove OTUs that represented by less than ten reads in a given sample.
ps1 <- subset_samples(ps.noncontam, sample_sums(ps.noncontam)>=200) # 2 samples removed
ps1 = prune_taxa(taxa_sums(ps1)>0, ps1)

ps2 = filter_taxa(ps1, function(x) sum(x >= 10) >= 1, TRUE)

### Make a data frame with a column for the read counts of each sample
depth_track3 <- data.frame(sample_data(ps2)) %>%
  mutate(quality_filtered = sample_sums(ps2)) %>%
  right_join(depth_track2)
