# Takes a look at Blast Results

# Load libraries:
library(tidyverse)
library(janitor)
library(naniar)  # Tools for NA

# Scientific notation options:
options(scipen=10)

# Get command line arguments:
# (CSV with MiXCR ratio and stats made using R) (full mixcr results) (blast results) (output directory)
args = commandArgs(trailingOnly=TRUE)
mixcr_ratio = args[1]
fullmixcr = args[2]
blastfile = args[3]
out_dir = args[4]


# Set working directory: 
# setwd("pipeline/")

# Read in files:
blast_results <- as_tibble(read_csv(blastfile))
# Strip "\n" from blast alignment string:
blast_results$alignment_string <-  str_replace_all(blast_results$alignment_string, 
                                                   pattern="\n", 
                                                   replacement="  ")
# Drop "X1" column:
blast_results <- blast_results %>% 
  select(-X1)

mixcr_ratio_results <- as_tibble(read_csv(mixcr_ratio))
fullmixcr_results <- as_tibble(read_tsv(fullmixcr, na=c(""),
                                        col_types=cols(.default="?",
                                                       subseq1="c",
                                                       subseq2="c",
                                                       subseq3="c",
                                                       subseq4="c",
                                                       subseq5="c",
                                                       subseq6="c",
                                                       subseq7="c",
                                                       subseq8="c",
                                                       subseq9="c")))

# #############################################################################
# Filter BLAST results so only top clone is subject and 
# second clone is BLAST
# For only subseq1
blast_filter <- blast_results %>% 
  filter(meta_subject_subseq=="subseq1" & meta_query_subseq =="subseq1") %>% 
  filter(meta_subject_chain != "IGH") %>% 
  filter(meta_subject_cRank==0, meta_query_cRank==1) %>% 
  group_by(meta_subject_name) %>% 
  slice_min(n=1, meta_subject_cID) %>% 
  ungroup() 
# The above will not catch errors produced by the top clone not matching the 
# most frequent light chain

# Check length:
print("Dimensions of BLAST results:")
dim(blast_results)
print("Dimensions of filtered blast results:")
dim(blast_filter)

# Rename blast columns:
colnames(blast_filter) <- paste("blast", colnames(blast_filter), sep="_")
# Save blast column names so can "NA" appropriate 
blast_colnames <- colnames(blast_filter)

# #############################################################################
# Combine filtered blast results with MiXCR:
# (Join by dominant chain and subject name)
# (Note: right_join keep=TRUE did not return the correct number of rowss)
# blast_combined <- right_join(x=blast_filter, y=mixcr_ratio_results, 
#                         by=c("blast_meta_subject_chain"="lc_with_most", "blast_meta_subject_name"="sample"))
blast_combined <- left_join(x=mixcr_ratio_results, y=blast_filter, 
                             by=c("lc_with_most"="blast_meta_subject_chain", "sample"="blast_meta_subject_name"), keep=TRUE)


print("Dimensions of MiXCR ratio results:")
dim(mixcr_ratio_results)
print("Dimensions of blast combined results:")
dim(blast_combined)


# #############################################################################
# Look at top two light chains for each:
top_2 <- fullmixcr_results %>% 
  filter(Chain != "IGH") %>% 
  group_by(sample) %>% 
  top_n(-cloneId, n=2)
#  # Keep only some columns:
#  top_2 <- top_2 %>% 
#    select(Chain, cloneRank, sample, cloneId, cloneCount, cloneFraction, 
#           allVHitsWithScore, allJHitsWithScore, nSeqCDR3, 
#           maxSubLength, subseq_count,
#           subseq1len, subseq2len, subseq3len, subseq4len, subseq5len,
#           subseq6len, subseq7len, subseq8len, subseq9len)
# Separate top 1 and 2:
firstClone <- top_2 %>% 
  top_n(-cloneId, n=1)
secondClone <- top_2 %>% 
  top_n(cloneId, n=1)
# Add suffix to variable names:
colnames(firstClone) <- paste("first", colnames(firstClone), sep="_")
colnames(secondClone) <- paste("second", colnames(secondClone), sep="_")
firstClone <- firstClone %>% 
  rename(sample=first_sample)
secondClone <- secondClone %>% 
  rename(sample=second_sample)

# Combine two dataframes:
top_2_wide <- full_join(firstClone, secondClone, by="sample")

# Check if the chain of the first equals the second:
top_2_wide$chain12_match <- top_2_wide$first_Chain == top_2_wide$second_Chain

# Save results:
# write.csv(top_2_wide, file=paste(out_dir, "top2_clones.csv", sep=""), row.names=FALSE)
write.table(top_2_wide, file=paste(out_dir, "top2_clones.txt", sep=""), 
            sep="\t", row.names=FALSE)

# Get counts:
sample_count <- length(top_2_wide$chain12_match)
print(sprintf("Sample count: %d", sample_count))
top2_match_count <- length(top_2_wide$chain12_match[top_2_wide$chain12_match])
print(sprintf("Number of match samples (where top 2 light chains are either IGK==IGK or IGL==IGL): %d", top2_match_count))
top2_mismatch_count <- length(top_2_wide$chain12_match[!top_2_wide$chain12_match])
print(sprintf("Number of mismatch samples (where top 2 light chains are either IGK&IGL or IGL&IGK): %d", top2_mismatch_count))
top_2_wide$chain12_matchf <- factor(top_2_wide$chain12_match, 
                                    levels=c(TRUE, FALSE),
                                    labels=c(paste("match", sep=""), 
                                             paste("MISMATCH", sep="")))
# top_2_wide$chain12_matchf <- factor(top_2_wide$chain12_match, 
#                                     levels=c(TRUE, FALSE),
#                                     labels=c(paste("match (n=", top2_match_count, ")", sep=""), 
#                                              paste("MISMATCH (n=", top2_mismatch_count, ")", sep="")))

# #############################################################################
# Combine top two clones with MiXCR and blast results:
# (NOTE: when keep=FALSE, not all sample names were kept.)
giant_combined <- full_join(blast_combined, top_2_wide, by=c("sample"), keep=TRUE)
# Check if sample.x equals sample.y:
print(sprintf("Dimensions: %d; Number of correctly matched sample names: %d", 
              length(giant_combined$sample.x), sum(giant_combined$sample.x==giant_combined$sample.y)))
if (length(giant_combined$sample.x) == sum(giant_combined$sample.x==giant_combined$sample.y)){
  print("Correct number of matching sample names")
} else {
  print("ERROR: After full_join in postblast.R, there is the wrong number of sample names")
}
# Rename sample column names
giant_combined <- giant_combined %>% 
  rename(sample=sample.x) %>% 
  select(-sample.y)

# New dimensions:
dim(giant_combined)

# Check if cloneID of top 2 clones matches appropriate alignment
blast_mismatched_cid <- giant_combined %>% 
  filter(blast_meta_subject_cID!=first_cloneId | 
           blast_meta_query_cID!=second_cloneId) %>% 
  select(sample, blast_meta_subject_cID, blast_meta_query_cID,
       first_cloneId, second_cloneId)
print("Samples where BLAST cloneId does not match top two clones:")
blast_mismatched_cid

# For inappropriate use of BLAST alignment, zero back down to "NA"
for (col in blast_colnames){
  # Make copy of original vector:
  temp_vector <- giant_combined[[col]]
  # Create temp vector with inapprorpriate BLAST replaced with NA:
  giant_combined <- giant_combined %>% 
    mutate(temp = ifelse(((blast_meta_subject_cID!=first_cloneId | blast_meta_query_cID!=second_cloneId) |
                             is.na(blast_meta_subject_cID)), NA, temp_vector))
  # Put temp vector back in right place:
  giant_combined[[col]] <- giant_combined$temp
  # Remove temp vector:
  giant_combined <- giant_combined %>% select(-temp)
}


# Make it so not found alignment is labeled:
giant_combined$blast_alignment_count <- as.character(giant_combined$blast_alignment_count)
giant_combined$blast_alignment_count <- replace_na(giant_combined$blast_alignment_count, "notAligned")
giant_combined$blast_alignment_count <- as_factor(giant_combined$blast_alignment_count)
summary(giant_combined$blast_alignment_count)


# Save df:
# write.csv(giant_combined, file=paste(out_dir, "mixcrstats_topclones_blast_combined.csv", sep=""), row.names=FALSE)
write.table(giant_combined, file=paste(out_dir, "mixcrstats_topclones_blast_combined.txt", sep=""), 
            sep="\t", row.names=FALSE)
# #############################################################################
# List those with NA alignment and low fraction LC top clone
print("Those missing an Alignment with a low fraction LC top clone AND Second Most LC clone is same chain:")
samples_of_interest <- giant_combined %>% 
  filter(blast_alignment_count=="notAligned" & first_cloneFraction<0.9 & chain12_match) %>% 
  arrange(first_cloneFraction) %>% 
  select(sample)
samples_of_interest[[1]]


# List those where most frequent light chain don't match top clone chain:
print("Samples where the most frequent light chain (IGL vs IGK) doesn't match top clone chain:")
chainMost_domChain_mismatch <- mixcr_ratio_results %>% 
  filter(lc_with_most != domClone_Chain)
chainMost_domChain_mismatch$sample

# #############################################################################
# Update blast_combined dataframe in order to handle fixing inappropriate alignments:
blast_combined_colnames <- names(blast_combined)
blast_combined <- giant_combined %>% 
  select(all_of(blast_combined_colnames))

# #############################################################################
# Length of dominant clone with and without MiXCR results:
# Violin plot for dominant clone subsequence length, maximum length, for light and heavy:
# Make longer version of dataframe:
length_longer <- blast_combined %>% 
  pivot_longer(cols=c(domClone_maxSubLength, longest_maxSubLength, 
                      domCloneH_maxSubLength, longestHeavy_maxSubLength), 
               names_to="graph_var", values_to="graph_value")
length_longer$graph_var <- factor(length_longer$graph_var, 
                                             levels=c("domClone_maxSubLength", 
                                                      "longest_maxSubLength", 
                                                      "domCloneH_maxSubLength", 
                                                      "longestHeavy_maxSubLength"))
# Plot Violin of length considering alignment count:
violin_length <- length_longer  %>% 
  ggplot(aes(x=graph_var, y=graph_value, fill=blast_alignment_count)) +
  geom_violin(position=position_dodge(0.8), width=0.8) +
  # geom_boxplot(position=position_dodge(0.8), width=0.05, color="black") +
  ggtitle("MiXCR Sequence Length vs BLAST Alignment Status") +
  ylab("Sequence Length") +
  xlab("")
# show(violin_length)
ggsave(paste(out_dir, "blast_violin_length.png", sep=""), violin_length, device=png, width=12, height=8, units="in")


# #############################################################################
# Add qualitative column:
blast_combined$blast_identity_percent_100f <- NA
blast_combined$blast_identity_percent_100f[blast_combined$blast_identity_percent == 100] <- "100%"
blast_combined$blast_identity_percent_100f[blast_combined$blast_identity_percent < 100] <- "<100%"
blast_combined$blast_identity_percent_100f <- factor(blast_combined$blast_identity_percent_100f)
summary(blast_combined$blast_identity_percent_100f)

# Identity, Coverage, other BLAST stats:
# Round e_value:
blast_combined$blast_E_value <- signif(blast_combined$blast_e_value, digits=3)
temp_test <- blast_combined[c("blast_e_value", "blast_E_value")]

# Make longer version of dataset:
blast_longer <- blast_combined %>% 
  pivot_longer(cols=c(blast_query_length, 
                      blast_subject_length, 
                      blast_score, 
                      blast_E_value, 
                      blast_alignment_length, 
                      blast_query_start, 
                      blast_subject_start,
                      blast_query_end, 
                      blast_subject_end, 
                      blast_identities, 
                      blast_identity_percent, 
                      blast_extension_available_from_query, 
                      blast_extension_available_from_subject, 
                      blast_query_coverage_percent, 
                      blast_subject_coverage_percent,
                      blast_gaps),
               names_to="graph_var", values_to="graph_value")
blast_longer$graph_var <- factor(blast_longer$graph_var, 
                                 levels=c("blast_query_length", 
                                          "blast_subject_length", 
                                          "blast_alignment_length", 
                                          "blast_score", 
                                          "blast_query_start", 
                                          "blast_subject_start",
                                          "blast_query_end", 
                                          "blast_subject_end", 
                                          "blast_identities", 
                                          "blast_identity_percent", 
                                          "blast_extension_available_from_query", 
                                          "blast_extension_available_from_subject", 
                                          "blast_query_coverage_percent", 
                                          "blast_subject_coverage_percent",
                                          "blast_gaps", 
                                          "blast_E_value"),
                                 labels=c("query_length", 
                                          "subject_length", 
                                          "alignment_length", 
                                          "score", 
                                          "query_start", 
                                          "subject_start",
                                          "query_end", 
                                          "subject_end", 
                                          "identities", 
                                          "identity_percent", 
                                          "extension_available_from_query", 
                                          "extension_available_from_subject", 
                                          "query_coverage_percent", 
                                          "subject_coverage_percent",
                                          "gaps", 
                                          "E_value"))

# Violin plot of blast stats:
blast_stats_violin <- blast_longer %>% 
  ggplot(aes(x=graph_var, y=graph_value, color=graph_var)) +
  geom_violin(position=position_dodge(0.8), width=0.8) +
  # geom_boxplot(position=position_dodge(0.8), width=0.05, color="black") +
  ggtitle("BLAST statistics") +
  xlab("") +
  ylab("") +
  facet_wrap(~graph_var, scales="free") +
  theme(legend.position = "none")
# show(blast_stats_violin)
ggsave(paste(out_dir, "blast_violin_stats.png", sep=""), blast_stats_violin, device=png, width=12, height=9, units="in")


# Separate out by identity percent:
blast_stats_violin2 <- blast_longer %>% 
  ggplot(aes(x=graph_var, y=graph_value, color=blast_identity_percent_100f)) +
  geom_violin(position=position_dodge(0.8), width=0.8) +
  ggtitle("BLAST statistics") +
  xlab("") +
  ylab("") +
  theme(axis.ticks.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(size=0.1)) +
  facet_wrap(~graph_var, scales="free") 
  # theme(legend.position = "none")
# show(blast_stats_violin2)
ggsave(paste(out_dir, "blast_violin_stats2.png", sep=""), blast_stats_violin2, device=png, width=16, height=9, units="in")


# #############################################################################
# Does chain for top two light chain clones match?
top2_boxplot <- top_2_wide %>% 
  ggplot(aes(x=chain12_matchf, y=first_cloneFraction, color=chain12_matchf)) +
  # geom_violin(width=0.8) +
  geom_boxplot(varwidth=TRUE) +
  ggtitle("Top 2 Light Chain Clones both IGK or IGL vs Top Clone Fraction") +
  xlab("Do the top 2 LC Clones use the same LC locus? (i.e. both IGK or IGL)") +
  ylab("Top Clone Clone Fraction") +
  theme(legend.position = "none")
# show(top2_boxplot)
ggsave(paste(out_dir, "top2_boxplot.png", sep=""), top2_boxplot, device=png, width=8, height=8, units="in")

# For the mismatch, look at top clone fraction vs light chain count:
mismatch_fractioncount <- giant_combined %>% 
  ggplot(aes(x=first_cloneFraction, y=light, color=blast_identity_percent)) +
  geom_point(size=1) +
  ggtitle("Top Clone Fraction, Total Light chain,\nTop 2 clone Chain Status") +
  xlab("Fraction of Light Chain Sequences are Top Clone") +
  ylab("Count Total Light Chain") +
  scale_color_gradient(low="red", high="blue", na.value="gray50") +
  facet_wrap(~chain12_matchf)
# show(mismatch_fractioncount)
ggsave(paste(out_dir, "mismatch_fractioncount.png", sep=""), mismatch_fractioncount, device=png, width=12, height=8, units="in")
# For the mismatch, look at top clone fraction vs light chain count:
mismatch_fractioncount2 <- giant_combined %>% 
  ggplot(aes(x=first_cloneFraction, y=light, color=blast_identity_percent)) +
  geom_point(size=1) +
  ggtitle("Top Clone Fraction, Total Light chain,\nTop 2 clone Chain Status\n(NA Identity in Red)") +
  xlab("Fraction of Light Chain Sequences are Top Clone") +
  ylab("Count Total Light Chain") +
  scale_color_gradient(low="black", high="gray75", na.value="red") +
  facet_wrap(~chain12_matchf)
# show(mismatch_fractioncount2)
ggsave(paste(out_dir, "mismatch_fractioncount2.png", sep=""), mismatch_fractioncount2, device=png, width=12, height=8, units="in")
# Add identity NA labels:
mismatch_fractioncount3 <- giant_combined %>% 
  ggplot(aes(x=first_cloneFraction, y=light, color=blast_identity_percent)) +
  geom_point(size=1) +
  geom_text(data=subset(giant_combined, blast_alignment_count=="notAligned" & first_cloneFraction<0.9 & chain12_match), 
            aes(label=blast_meta_subject_name), hjust=-0.1, size=2) + 
  ggtitle("Top Clone Fraction, Total Light chain,\nTop 2 clone Chain Status\n(NA Identity in Red)") +
  xlab("Fraction of Light Chain Sequences are Top Clone") +
  ylab("Count Total Light Chain") +
  scale_color_gradient(low="black", high="gray75", na.value="red") +
  facet_wrap(~chain12_matchf)
# show(mismatch_fractioncount3)
ggsave(paste(out_dir, "mismatch_fractioncount3.png", sep=""), mismatch_fractioncount3, device=png, width=14, height=8, units="in", dpi=600)


# #############################################################################
# Identity vs Coverage
identity_querycoverage <- giant_combined %>% 
  ggplot(aes(x=blast_query_coverage_percent, y=blast_identity_percent, 
             color=blast_subject_coverage_percent)) +
  geom_point(size=1) +
  ggtitle("Identity & Query Coverage (Second Clone)") +
  xlab("Query Coverage (%)") +
  ylab("Identity (%)") +
  scale_color_gradient(low="red", high="blue", na.value="gray50")
# show(identity_querycoverage)
ggsave(paste(out_dir, "identity_querycoverage.png", sep=""), identity_querycoverage, device=png, width=12, height=8, units="in")
identity_subjectcoverage <- giant_combined %>% 
  ggplot(aes(x=blast_subject_coverage_percent, y=blast_identity_percent, 
             color=blast_query_coverage_percent)) +
  geom_point(size=1) +
  ggtitle("Identity & Subject Coverage (First Clone)") +
  xlab("Subject Coverage (%)") +
  ylab("Identity (%)") +
  scale_color_gradient(low="red", high="blue", na.value="gray50")
# show(identity_subjectcoverage)
ggsave(paste(out_dir, "identity_subjectcoverage.png", sep=""), identity_subjectcoverage, device=png, width=12, height=8, units="in")

# Extension available:
extension_plot <- giant_combined %>% 
  ggplot(aes(x=blast_extension_available_from_subject, 
             y=blast_extension_available_from_query, 
             color=blast_identity_percent)) +
  geom_point(size=1) +
  ggtitle("Length not present in other clone") +
  xlab("Subject Length Only (First clone)") +
  ylab("Query Length Only (Second clone)") +
  scale_color_gradient(low="red", high="blue", na.value="gray50")
# show(extension_plot)
ggsave(paste(out_dir, "extension_plot.png", sep=""), extension_plot, device=png, width=12, height=8, units="in")


# Query length only vs length of dominant clone:
extension_domlength <- giant_combined %>% 
  ggplot(aes(x=domClone_maxSubLength, 
             y=blast_extension_available_from_query, 
             color=blast_identity_percent)) +
  geom_point(size=1) +
  ggtitle("Length not present in top clone vs Top Clone Length") +
  xlab("Dominant Clone Length") +
  ylab("Exclusive Query Length (Second clone)") +
  scale_color_gradient(low="red", high="blue", na.value="gray50")
# show(extension_domlength)
ggsave(paste(out_dir, "extension_domlength.png", sep=""), extension_domlength, device=png, width=12, height=8, units="in")
# #############################################################################
# Identity percent:
png(filename=paste(out_dir, "hist_percent_identity.png", sep=""), width=8, height=6, units="in", res=300)
hist(giant_combined$blast_identity_percent, breaks=100, main="Percent Identity",
     xlab="Percent Identity")
dev.off()

identity_percents <- giant_combined$blast_identity_percent[!is.na(giant_combined$blast_identity_percent)]
length(identity_percents)
length(identity_percents[identity_percents >= 98])
length(identity_percents[identity_percents >= 99])
length(identity_percents[identity_percents >= 99.5])
length(identity_percents[identity_percents == 100])


# #############################################################################
# Look at second clone longest alignment onto first clone smaller segments
blast_filter_reflectback <- blast_results %>% 
  filter(meta_subject_subseq =="subseq1" & meta_query_subseq != "subseq1") %>% 
  filter(meta_subject_chain != "IGH") %>% 
  filter(meta_subject_cRank==1, meta_query_cRank==0) %>% 
  group_by(meta_subject_name) %>% 
  slice_min(n=1, meta_query_cID) %>% 
  ungroup() %>% 
  distinct()
# The above will not catch errors produced by the top clone not matching the 
# most frequent light chain

# Check length:
print("Dimensions of BLAST results:")
dim(blast_results)
print("Dimensions of filtered blast results, comparing longest segment clone 2 to small pieces clone 1:")
dim(blast_filter_reflectback)

# Save the results:
write.csv(blast_filter_reflectback, file=paste(out_dir, "blastClone1ShorterSegmentsAgainstSecondCloneLongest.csv", sep=""), row.names=FALSE)


# #############################################################################

