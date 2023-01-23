# Script categoriezes samples. Check global variables to adjust rules used to define the different categories.

# Load libraries:
library(tidyverse)
# library(janitor)
library(epitools) # for statistical test
library(nnet)  # Multinomial model

# Categorize samples

# Global variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Categories are as further, samples will be categorized in the first category conditions are met
# Category 1: samples ready to go without further processing 
# (i.e. top clone >= 95% of LC)
cat1_thresh <- 0.95
# Category 2: top2 clones push over same threshold
# (i.e. top 2 clones same chain, 100% identical over overlap, >= 200 bp alignment length,
# fraction higher than cat1_thresh)
# Note: identity is already expressed as a percent
cat2_identity_thresh <- 100
cat2_length_thresh <- 200
# Category 2b:
# Top 2 clones not the same type
# Category 3: Definite major clone, but doesn't pass threshold
# (i.e. top clones > 50%, second clone <= 5%)
cat3_minfreqtopclone <- 0.5
cat3_maxfreqsecondclone <- 0.05
# Category 4: Still unsatisfied

# Scientific notation options:
options(scipen=10000)

# Get command line arguments:
# (mixcrstats_topclones_blast_combined.txt) (output directory)
args = commandArgs(trailingOnly=TRUE)
top2 = args[1]
out_dir = args[2]

# # Files (if not command line):
# top2 <- "sra/800MiXCR/SUMMARY/mixcrstats_topclones_blast_combined.txt"
# # Set out directory:
# out_dir = "sra/800MiXCR/SUMMARY/"


# Read in df:
combined_df <- as_tibble(read_tsv(top2))

# Recalculate frequencies:
combined_df$firstclone_CorrectedFraction <- combined_df$first_cloneCount / combined_df$light
combined_df$secondclone_CorrectedFraction <- combined_df$second_cloneCount / combined_df$light

# Count of top two chains together:
combined_df$top2count <- combined_df$first_cloneCount + combined_df$second_cloneCount
# Make NA if chains don't match:
combined_df$top2count[combined_df$first_Chain != combined_df$second_Chain] <- NA
# Combined fraction of top two chains:
combined_df$top2fraction <- combined_df$top2count / combined_df$light

# Check if there are any NA's for CDR3:
sum(is.na(combined_df$first_nSeqCDR3))
sum(is.na(combined_df$second_nSeqCDR3))

# Does top two CDR3 match?
combined_df <- combined_df %>% 
  mutate(top2CDR3 = ifelse(first_nSeqCDR3==second_nSeqCDR3, "CDR3match", "CDR3mismatch"))
combined_df$top2CDR3 <- as_factor(combined_df$top2CDR3)
summary(combined_df$top2CDR3)


# Write df:
# write.csv(combined_df, file=paste(out_dir, "mixcrstats_topclones_blast_combined_correctedFractions.csv", sep=""), row.names=FALSE)
write.table(combined_df, file=paste(out_dir, "mixcrstats_topclones_blast_combined_correctedFractions.txt", sep=""), 
            sep="\t", row.names=FALSE)

# Total count:
total_n <- length(combined_df$sample)
print(sprintf("Total sample count: %d", total_n))

# Category 1: Total where first clone >= 0.95%
cat1 <- combined_df[combined_df$firstclone_CorrectedFraction >= cat1_thresh,]
cat1_n <- length(cat1$sample)
print(sprintf("Count in category 1: %d", cat1_n))
not_cat1 <- combined_df[combined_df$firstclone_CorrectedFraction < cat1_thresh,]
not_cat1_n <- length(not_cat1$sample)
print(sprintf("Count not in category 1: %d", not_cat1_n))
summary(cat1$top2CDR3)

# Category 2:
# Category 2: top2 clones push over same threshold
# (i.e. top 2 clones same chain, 100% identical over overlap, >= 200 bp alignment length,
# fraction higher than cat1_thresh)
# No CDR3 contraint cat2 <- not_cat1 %>% 
# No CDR3 contraint   filter(chain12_match & 
# No CDR3 contraint            blast_identity_percent==cat2_identity_thresh &
# No CDR3 contraint            blast_alignment_length>=cat2_length_thresh & 
# No CDR3 contraint            top2fraction>=cat1_thresh)
# Add CDR3 contraint:

# Look at summary for relevant information for defining category 2:
summary(not_cat1 %>% select(chain12_match, 
                            blast_identity_percent, 
                            blast_alignment_length,
                            top2fraction,
                            first_nSeqCDR3, second_nSeqCDR3,
                            blast_assembled_string))

cat2 <- not_cat1 %>% 
  filter(chain12_match & 
           blast_identity_percent>=cat2_identity_thresh &
           blast_alignment_length>=cat2_length_thresh & 
           top2fraction>=cat1_thresh & 
           first_nSeqCDR3==second_nSeqCDR3 &
           blast_assembled_string!="NotAvailable" &
           !is.na(blast_assembled_string))

cat2_n <- length(cat2$sample)
print(sprintf("Count in category 2: %d", cat2_n))
not_cat2 <- not_cat1 %>% 
  filter(!sample %in% cat2$sample)
not_cat2_n <- length(not_cat2$sample)
print(sprintf("Count not in category 2: %d", not_cat2_n))
summary(cat2$top2CDR3)

# Category 3:  Definite major clone, but doesn't pash threshold:
cat3 <- not_cat2 %>% 
  filter(firstclone_CorrectedFraction>cat3_minfreqtopclone & 
           secondclone_CorrectedFraction<=cat3_maxfreqsecondclone)
cat3_n <- length(cat3$sample)
print(sprintf("Count in category 3: %d", cat3_n))
not_cat3 <- not_cat2 %>% 
  filter(!sample %in% cat3$sample)
not_cat3_n <- length(not_cat3$sample)
print(sprintf("Count not in category 3: %d", not_cat3_n))  
summary(cat3$top2CDR3)
summary(not_cat3$top2CDR3)

# Save lists of samples in each category
write.table(pull(cat1, sample), paste(out_dir, "cat1list_ReadyToGo.txt", sep=""), 
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pull(cat2, sample), paste(out_dir, "cat2list_Top2Together.txt", sep=""), 
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pull(cat3, sample), paste(out_dir, "cat3list_MajorClonepresentbelowthreshold.txt", sep=""), 
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(pull(not_cat3, sample), paste(out_dir, "cat4list_uncategorizedsamples.txt", sep=""), 
            row.names=FALSE, col.names=FALSE, quote=FALSE)

# Add categorization column to main table:
# First do piece wise categorization:
cat1_samples <- tibble(
  sample = cat1$sample, 
  CATEGORY = "Category1_ReadyToGo"
)
cat2_samples <- tibble(
  sample = cat2$sample, 
  CATEGORY = "Category2_Top2Together"  
)
cat3_samples <- tibble(
  sample = cat3$sample, 
  CATEGORY = "Category3_MajorClonepresentbelowthreshold"  
)
cat4_samples <- tibble(
  sample = not_cat3$sample, 
  CATEGORY = "Category4_uncategorizedsamples"  
)
# Bind together the pieces:
sample_categories <- bind_rows(cat1_samples, cat2_samples, cat3_samples, cat4_samples)
# Combine with main df:
combined_df2 <- full_join(combined_df, sample_categories, by="sample")

# Check out category counts:
combined_df2$CATEGORY <- factor(combined_df2$CATEGORY, levels=c("Category1_ReadyToGo", 
                                                                "Category2_Top2Together", 
                                                                "Category3_MajorClonepresentbelowthreshold",
                                                                "Category4_uncategorizedsamples"))
summary(combined_df2$CATEGORY)

# Save new dataframe:
write.table(combined_df2, file=paste(out_dir, "mixcrstats_topclones_blast_combined_correctedFractions_withCategorizations.txt", sep=""), 
            sep="\t", row.names=FALSE)

# #############################################################################
# Potential polyclonality
# Number where the top two clones sum to less than 50
sum(not_cat3$top2fraction < 0.50, na.rm=TRUE)
# Inverse:
sum(not_cat3$top2fraction >= 0.50, na.rm=TRUE)
# IGK - IGL mismatch:
sum(not_cat3$chain12_matchf!="match")
# top 2 < 50% + clone 1>=5% + clone 2>=5%
sum((not_cat3$top2fraction < 0.50) & 
      (not_cat3$firstclone_CorrectedFraction>=0.05) & 
      (not_cat3$secondclone_CorrectedFraction>=0.05), na.rm=TRUE)
# top 2 < 50% + clone 1>=10% + clone 2>=10%
sum((not_cat3$top2fraction < 0.50) & 
      (not_cat3$firstclone_CorrectedFraction>=0.1) & 
      (not_cat3$secondclone_CorrectedFraction>=0.1), na.rm=TRUE)
# top 2 < 50% + clone 1>=5% + clone 2>=5% + CDR3 mismatch
sum((not_cat3$top2fraction < 0.50) & 
      (not_cat3$firstclone_CorrectedFraction>=0.05) & 
      (not_cat3$secondclone_CorrectedFraction>=0.05) &
      (not_cat3$top2CDR3=="CDR3mismatch"), na.rm=TRUE)

# Look at fraction compared with match status:
print(not_cat3[, c("top2fraction", "chain12_matchf", "blast_alignment_count", 
                   "firstclone_CorrectedFraction", "secondclone_CorrectedFraction")] 
      %>% arrange(desc(chain12_matchf)), n=100)

# Top LC clone <95%, first and second >=10%
sum((not_cat3$firstclone_CorrectedFraction<0.95) &
  (not_cat3$firstclone_CorrectedFraction>=0.1) & 
      (not_cat3$secondclone_CorrectedFraction>=0.1), na.rm=TRUE)
# Top LC clone <95%, first and second >=10% + CDR3 mismatch:
sum((not_cat3$firstclone_CorrectedFraction<0.95) &
      (not_cat3$firstclone_CorrectedFraction>=0.10) & 
      (not_cat3$secondclone_CorrectedFraction>=0.10) &
      (not_cat3$top2CDR3=="CDR3mismatch"), na.rm=TRUE)
# Top LC clone <95%, first and second >=10% + no BLAST alignment
sum((not_cat3$firstclone_CorrectedFraction<0.95) &
      (not_cat3$firstclone_CorrectedFraction>=0.1) & 
      (not_cat3$secondclone_CorrectedFraction>=0.1) &
      (not_cat3$blast_alignment_count=="notAligned"), na.rm=TRUE)
# Top LC clone <95%, first and second >=10% + no BLAST alignment + CDR3 mismatch
sum((not_cat3$firstclone_CorrectedFraction<0.95) &
      (not_cat3$firstclone_CorrectedFraction>=0.1) & 
      (not_cat3$secondclone_CorrectedFraction>=0.1) &
      (not_cat3$blast_alignment_count=="notAligned") &
      (not_cat3$top2CDR3=="CDR3mismatch"), na.rm=TRUE)

# Top LC clone <95%, first and second each >=10%, % identity <=98%
sum((not_cat3$firstclone_CorrectedFraction<0.95) &
      (not_cat3$firstclone_CorrectedFraction>=0.10) & 
      (not_cat3$secondclone_CorrectedFraction>=0.10) &
      (not_cat3$blast_identity_percent<=98), na.rm=TRUE)
# Top LC clone <95%, first and second each >=10%, % identity <=98%
sum((not_cat3$firstclone_CorrectedFraction<0.95) &
      (not_cat3$firstclone_CorrectedFraction>=0.10) & 
      (not_cat3$secondclone_CorrectedFraction>=0.10) &
      (not_cat3$blast_identity_percent<=98) &
      (not_cat3$top2CDR3=="CDR3mismatch"), na.rm=TRUE)
# Top LC clone <95% second >=10% + % identity <=98% + alignment length >= 500
sum((not_cat3$firstclone_CorrectedFraction<0.95) &
      (not_cat3$firstclone_CorrectedFraction>=0.10) & 
      (not_cat3$secondclone_CorrectedFraction>=0.10) &
      (not_cat3$blast_identity_percent<=98) &
      (not_cat3$first_maxSubLength>=500) & 
      (not_cat3$second_maxSubLength>=500), na.rm=TRUE)
# Top LC clone <95%, first and second >=10% + % identity <=98% + alignment length >= 500 + CDR3 mismatch
sum((not_cat3$firstclone_CorrectedFraction<0.95) &
      (not_cat3$firstclone_CorrectedFraction>=0.10) & 
      (not_cat3$secondclone_CorrectedFraction>=0.10) &
      (not_cat3$blast_identity_percent<=98) &
      (not_cat3$first_maxSubLength>=500) & 
      (not_cat3$second_maxSubLength>=500) &
      (not_cat3$top2CDR3=="CDR3mismatch"), na.rm=TRUE)

# Top LC clone <95%, first and second >=10%, + IGK IGL chain mismatch
sum((not_cat3$firstclone_CorrectedFraction<0.95) &
      (not_cat3$firstclone_CorrectedFraction>=0.10) & 
      (not_cat3$secondclone_CorrectedFraction>=0.10) &
      (not_cat3$chain12_matchf!="match"), na.rm=TRUE)
# Top LC clone <95%, first and second >=10%, + IGK IGL chain mismatch + CDR3 mismatch
sum((not_cat3$firstclone_CorrectedFraction<0.95) &
      (not_cat3$firstclone_CorrectedFraction>=0.10) & 
      (not_cat3$secondclone_CorrectedFraction>=0.10) &
      (not_cat3$chain12_matchf!="match") &
      (not_cat3$top2CDR3=="CDR3mismatch"), na.rm=TRUE)


# #############################################################################
# IGK vs IGL for each category

# Add category column to each:
cat1$category <- "Cat1"
cat2$category <- "Cat2"
cat3$category <- "Cat3"
not_cat3$category <- "Cat4"
recombined_df <- bind_rows(cat1, cat2, cat3, not_cat3)
recombined_df$category <- as_factor(recombined_df$category)
recombined_df$domClone_Chain <- as_factor(recombined_df$domClone_Chain)

# Tabulate just IGK and IGL
chain_table <- table(recombined_df$domClone_Chain)
chain_table
# Proportion:
prop.table(chain_table)

# Tabulate:
cat_chain_table <- table(recombined_df$category, recombined_df$domClone_Chain)
cat_chain_table
# Proportion table:
prop.table(cat_chain_table)
# Proportion each category out of 1:
prop.table(cat_chain_table, margin=1)
# Proportion table IGK IGL out of 1:
prop.table(cat_chain_table, margin=2)

# Tabulate category 4:
table(not_cat3$first_Chain, not_cat3$second_Chain)

# Relabel match mismatch:
# not_cat3$chain12_matchf <- factor(not_cat3$chain12_matchf, 
#                                   levels=c("match (n=686)", "MISMATCH (n=81)"), 
#                                   labels=c("match", "MISMATCH"))
# Tabulate category 4 match & failed blast:
table(not_cat3$chain12_matchf, not_cat3$blast_alignment_count)
# #########################################################################

# Tabulate in other direction:
cat_chain_table_T <- table(recombined_df$domClone_Chain, recombined_df$category)
cat_chain_table_T
# Proportion table IGK IGL out of 1:
prop.table(cat_chain_table_T, margin=1)

# Chi-squared check:
cat_chain_table_T
cat_res_T <- chisq.test(cat_chain_table_T, correct=FALSE)
# Chi-squared check if expected more than 5:
cat_res_T$expected
# Chi-squared
cat_res_T
# Other way:
cat_chain_table
cat_res <- chisq.test(cat_chain_table, correct=FALSE)
# Chi-squared check if expected more than 5:
cat_res$expected
# Chi-squared
cat_res

# Odds ratio:  
or <- oddsratio(cat_chain_table)
or

# Do it by logistic regression:
glm1 <- glm(recombined_df$domClone_Chain ~ recombined_df$category, family=binomial(link="logit"))
sum1 <- summary(glm1)
sum1
# Odds ratio:
exp(cbind(OR=coef(glm1), confint(glm1)))

# Do it the other way: 
# https://www.r-bloggers.com/2020/05/multinomial-logistic-regression-with-r/
glm2 <- multinom(recombined_df$category ~ recombined_df$domClone_Chain)
sum2 <- summary(glm2)
sum2
# Odds ratio:
exp(coef(glm2))
exp(confint(glm2))

# #########################################################################
# Who had blast alignment and yet had a mismatch?
cat4_blast_mismatch  <- not_cat3 %>% 
  filter(chain12_matchf=="MISMATCH", blast_alignment_count=="1")

# Samples in category 4 with low identity
cat4_lowidentity <- not_cat3 %>% 
  filter(blast_identity_percent<=98) %>% 
  arrange(blast_identity_percent) %>% 
  pull(sample)
cat4_lowidentity
# Save to text file:
write.table(cat4_lowidentity, paste(out_dir, "cat4sublist_lowidentiy.txt", sep=""), 
            row.names=FALSE, col.names=FALSE, quote=FALSE)
# Look at table with applicable results:
cat4_lowidentity_df <- not_cat3 %>% 
  filter(blast_identity_percent<=98) %>% 
  select(sample, blast_identity_percent, blast_alignment_length, top2CDR3) %>% 
  arrange(blast_identity_percent)
cat4_lowidentity_df

# Get minimum value of identity for CDR3 match:
not_cat3 %>% 
  filter(top2CDR3 == "CDR3match") %>% 
  top_n(-1, wt=blast_identity_percent) %>% 
  select(sample, blast_identity_percent, blast_alignment_length, top2CDR3)
# Get maximum value of identity for CDR3 mismatch:
not_cat3 %>% 
  filter(top2CDR3 == "CDR3mismatch") %>% 
  top_n(1, wt=blast_identity_percent) %>% 
  select(sample, blast_identity_percent, blast_alignment_length, top2CDR3)

# Figure category 4 for clone 1 & 2:
cat4_fraction_plot <- not_cat3 %>% 
  ggplot(aes(x=firstclone_CorrectedFraction, y=secondclone_CorrectedFraction, color=first_Chain, shape=chain12_matchf)) +
  geom_point(size=1.5) +
  ggtitle("Category 4 Frequency of Top 2 Clones") +
  xlab("First Clone Fraction") + 
  ylab("Second Clone Fraction")
show(cat4_fraction_plot)
ggsave(paste(out_dir, "cat4_fraction_plot.png", sep=""), cat4_fraction_plot, device=png, width=6, height=4, units="in")
# Label low identity:
cat4_fraction_plot_v2 <- cat4_fraction_plot +
  geom_text(data=subset(not_cat3, sample %in% cat4_lowidentity), 
            aes(label=sample), 
            hjust=-0.2, size=1) 
show(cat4_fraction_plot_v2)
ggsave(paste(out_dir, "cat4_fraction_plot_v2.png", sep=""), cat4_fraction_plot_v2, device=png, width=9, height=6, units="in", dpi=600)
# Color for CDR3 mismatch:
cat4_fraction_plot_v3 <- not_cat3 %>% 
  ggplot(aes(x=firstclone_CorrectedFraction, y=secondclone_CorrectedFraction, shape=first_Chain, color=top2CDR3)) +
  geom_point(size=1.5) +
  ggtitle("Category 4 Frequency of Top 2 Clones") +
  xlab("First Clone Fraction") + 
  ylab("Second Clone Fraction") +
  geom_text(data=subset(not_cat3, sample %in% cat4_lowidentity), 
            aes(label=sample), 
            hjust=-0.2, size=1) 
show(cat4_fraction_plot_v3)
ggsave(paste(out_dir, "cat4_fraction_plot_v3.png", sep=""), cat4_fraction_plot_v3, device=png, width=9, height=6, units="in", dpi=600)


# Mismatches don't get graphed since there is no sum value:
cat4_fraction_plot2 <- not_cat3 %>% 
  ggplot(aes(x=firstclone_CorrectedFraction, y=top2fraction, color=first_Chain, shape=blast_alignment_count)) +
  geom_point(size=1.5) +
  ggtitle("Category 4 Frequency of Top 2 Clones") +
  xlab("First Clone Fraction") + 
  ylab("First + Second Clone Fraction\n(IGK IGL mismatches excluded)")
show(cat4_fraction_plot2)
ggsave(paste(out_dir, "cat4_fraction_plot2.png", sep=""), cat4_fraction_plot2, device=png, width=6, height=4, units="in")
# Label low identity:
cat4_fraction_plot2_v2 <- cat4_fraction_plot2 +
  geom_text(data=subset(not_cat3, sample %in% cat4_lowidentity), 
            aes(label=sample), 
            hjust=-0.2, size=1) 
show(cat4_fraction_plot_v2)
ggsave(paste(out_dir, "cat4_fraction_plot2_v2.png", sep=""), cat4_fraction_plot2_v2, device=png, width=9, height=6, units="in", dpi=600)
# Color for CDR3 mismatch:
cat4_fraction_plot2_v3 <- not_cat3 %>% 
  ggplot(aes(x=firstclone_CorrectedFraction, y=top2fraction, color=top2CDR3, shape=blast_alignment_count)) +
  geom_point(size=1.5) +
  ggtitle("Category 4 Frequency of Top 2 Clones") +
  xlab("First Clone Fraction") + 
  ylab("First + Second Clone Fraction\n(IGK IGL mismatches excluded)") +
  geom_text(data=subset(not_cat3, sample %in% cat4_lowidentity), 
            aes(label=sample), 
            hjust=-0.2, size=1) 
show(cat4_fraction_plot2_v3)
ggsave(paste(out_dir, "cat4_fraction_plot2_v3.png", sep=""), cat4_fraction_plot2_v3, width=9, height=6, units="in", dpi=600)

# Category 4 identity and alignment length:
cat4_blast_plot <- not_cat3 %>% 
  ggplot(aes(x=blast_identity_percent, y=blast_alignment_length, color=top2fraction, shape=top2CDR3)) +
  geom_point(size=1.5) +
  ggtitle("Category 4 Identity and Alignment") +
  xlab("% Identity between Aligned Region Clone 1 & 2") + 
  ylab("Alignment Length")
show(cat4_blast_plot)
ggsave(paste(out_dir, "cat4_blast_plot.png", sep=""), cat4_blast_plot, device=png, width=6, height=4, units="in")
# Add sample name text to figure:
cat4_blast_plot_v2 <- cat4_blast_plot +
  geom_text(data=subset(not_cat3, sample %in% cat4_lowidentity), 
            aes(label=sample), 
            hjust=-0.2, size=1) 
show(cat4_blast_plot_v2)
ggsave(paste(out_dir, "cat4_blast_plot_v2.png", sep=""), cat4_blast_plot_v2, device=png, width=9, height=6, units="in", dpi=600)
# Shape for CDR3 mismatch:
cat4_blast_plot_v3 <- not_cat3 %>% 
  ggplot(aes(x=blast_identity_percent, y=blast_alignment_length, color=top2fraction, shape=top2CDR3)) +
  geom_point(size=1.5) +
  ggtitle("Category 4 Identity and Alignment") +
  xlab("% Identity between Aligned Region Clone 1 & 2") + 
  ylab("Alignment Length") +
  geom_text(data=subset(not_cat3, sample %in% cat4_lowidentity), 
            aes(label=sample), 
            hjust=-0.2, size=1) 
show(cat4_blast_plot_v3)
ggsave(paste(out_dir, "cat4_blast_plot_v3.png", sep=""), cat4_blast_plot_v3, device=png, width=9, height=6, units="in", dpi=600)

# #############################################################################
# Barplot of subseq counts:
combined_df2 %>% 
  ggplot(aes(x=first_subseq_count, fill=CATEGORY, color=CATEGORY)) +
  geom_bar(position=position_dodge()) +
  ggtitle("Number of Subsequences for the Top Clone") +
  scale_y_continuous(n.breaks=8) +
  xlab("Number of MiXCR segments") +
  geom_text(stat="count", aes(label=..count..), vjust=-0.5, position=position_dodge(width=0.9))
combined_df2 %>% 
  ggplot(aes(x=second_subseq_count, fill=CATEGORY, color=CATEGORY)) +
  geom_bar(position=position_dodge()) +
  ggtitle("Number of Subsequences for the Second Clone") +
  scale_y_continuous(n.breaks=8) +
  xlab("Number of MiXCR segments") +
  geom_text(stat="count", aes(label=..count..), vjust=-0.5, position=position_dodge(width=0.9)) 
# Do stacked with chain
# First convert chain to factor
combined_df2$domClone_Chain <- factor(combined_df2$domClone_Chain)
combined_df2$first_Chain <- factor(combined_df2$first_Chain)
combined_df2$second_Chain <- factor(combined_df2$second_Chain)
# Stacked for First Clone:
combined_df2 %>% 
  ggplot(aes(x=first_subseq_count, fill=first_Chain)) +
  geom_bar() +
  facet_wrap(.~CATEGORY, scales="free", nrow=1) +
  ggtitle("Number of Subsequences for the Top Clone") +
  xlab("Number of MiXCR segments") +
  geom_text(stat="count", aes(label=..count..), position=position_stack(vjust=0.5))
  # scale_x_continuous(n.breaks=4)
# Stacked for Second Clone:
combined_df2 %>% 
  ggplot(aes(x=second_subseq_count, fill=second_Chain)) +
  geom_bar() +
  facet_wrap(.~CATEGORY, scales="free", nrow=1) +
  ggtitle("Number of Subsequences for the Second Clone") +
  xlab("Number of MiXCR segments") +
  geom_text(stat="count", aes(label=..count..), position=position_stack(vjust=0.5))

# Histogram of subseq lengths:
# Pivot_longer:
subseq_length_columns <- c(grep(glob2rx("first_subseq*len"), names(combined_df2), value=TRUE), 
                           grep(glob2rx("second_subseq*len"), names(combined_df2), value=TRUE))
length_long <- combined_df2 %>% 
  pivot_longer(cols=one_of(subseq_length_columns), names_to="graph_var", values_to="graph_value")
# Create a variable to separate out first and second and subseq:
length_long$graph_var2 <- length_long$graph_var
length_long <- length_long %>% 
  separate("graph_var", c("graph_clone", "graph_subseq"), "_")
# Faceted histogram:
length_long %>% 
  filter(graph_subseq %in% c("subseq1len", "subseq2len", "subseq3len")) %>% 
  # filter(graph_value>0) %>% 
  ggplot(aes(x=graph_value, color=CATEGORY)) +
  geom_histogram(position="identity", fill="transparent", bins=10) +
  facet_wrap(CATEGORY ~ graph_var2, nrow=4, scales="free") +
  # facet_wrap(graph_var2 ~ CATEGORY, scales="free") +
  xlab("Subsequence Length") +
  stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-0.5, bins=10, size=3) +
  scale_y_continuous(expand=expansion(mult=c(0,0.2)))

# #############################################################################
# List samples with multiple segments from category 1:
cat1 %>% 
  filter(first_subseq_count > 1) %>% 
  select(sample, first_subseq_count, first_subseq1len, first_subseq2len) %>% 
  arrange(desc(first_subseq_count), first_subseq1len) %>% 
  print(n=Inf)

# ############################################################################## 
# Subseq1 vs Subseq2 length of top clone
combined_df2 %>% 
  filter(first_subseq_count > 1) %>% 
  ggplot(aes(x=first_subseq1len, y=first_subseq2len, color=CATEGORY)) +
  geom_point() +
  ggtitle("Top Clone Subseq1 and Subseq2 Length, when multiple Subsequences are present") +
  xlab("Subsequence 1 Length") +
  ylab("Subsequence 2 Length") +
  facet_wrap(~ CATEGORY, nrow=2)
combined_df2 %>% 
  ggplot(aes(x=first_subseq1len, y=first_subseq2len, color=CATEGORY)) +
  geom_point() +
  ggtitle("Top Clone Subseq1 and Subseq2 Length") +
  xlab("Subsequence 1 Length") +
  ylab("Subsequence 2 Length") +
  facet_wrap(~ CATEGORY, nrow=2)


# #############################################################################
# Figures colored by category

# LC count vs dom clone frequency
lcN_domcloneN <- combined_df2  %>% 
  ggplot(aes(x=light, y=firstclone_CorrectedFraction)) +
  geom_point(aes(color=CATEGORY), size=2, alpha=0.75) +
  ggtitle("LC Count vs\nDominant Light Chain Clone Frequency") +
  xlab("LC Count") +
  ylab("Dominant LC Clone Frequency")
# show(lcN_domcloneN)
ggsave(paste(out_dir, "lcN_domcloneN_category.png", sep=""), lcN_domcloneN, device=png, width=8, height=6, units="in")

# First clone vs second clone by category
lcN_domcloneN <- combined_df2  %>% 
  ggplot(aes(y=secondclone_CorrectedFraction, x=firstclone_CorrectedFraction)) +
  geom_point(aes(color=CATEGORY), size=2, alpha=0.75) +
  ggtitle("Dominant vs Second Light Chain Clone Frequency") +
  ylab("Second clone Frequency") +
  xlab("Dominant LC Clone Frequency")


# LC count vs number of clones in log scale:
lcN_domfreq <- combined_df2   %>% 
  ggplot(aes(y=number_of_clones, x=light)) +
  geom_point(aes(color=CATEGORY), size=2, alpha=0.75) +
  ggtitle("Number of LC Clones\nvs LC Count ") +
  ylab("Number of Distinct LC Clones") +
  xlab("LC Count")
lcN_domfreqlog <- lcN_domfreq +
  scale_y_continuous(trans="log10") +
  ylab("Number of Distinct LC Clones (log10)")
# show(lcN_domfreqlog)
ggsave(paste(out_dir, "lcN_domfreqlog_category.png", sep=""), lcN_domfreqlog, device=png, width=8, height=6, units="in")

# Identity vs Alignment length
modified_df <- combined_df2 %>% 
  replace_na(list(blast_alignment_length = 0, blast_identity_percent = 0))
identity_alignmentlength <- modified_df %>% 
  ggplot(aes(x=blast_alignment_length, y=blast_identity_percent)) +
  geom_point(aes(color=CATEGORY), size=2, alpha=0.5) +
  ggtitle("Identity and alignment length") +
  xlab("Alignment length (nt)") +
  ylab("Identity (%)") # +
  # scale_color_gradient(low="red", high="blue", na.value="gray50")
show(identity_alignmentlength)

# Identity by category:
modified_df %>% 
  ggplot(aes(x=blast_identity_percent)) +
  geom_histogram() +
  facet_wrap(~CATEGORY, nrow=2)
# Drop 0s:
modified_df %>% 
  filter(blast_identity_percent != 0) %>% 
  ggplot(aes(x=blast_identity_percent)) +
  geom_histogram(bins=50) +
  facet_wrap(~CATEGORY, nrow=2)


# #############################################################################
