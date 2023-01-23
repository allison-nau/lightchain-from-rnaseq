# Script to look at sequence coverage using modified IMGT numbering based on High V-Quest results

# Load libraries:
library(tidyverse)
library(janitor)
library(stringr)
library(stringi)
library(forcats)

# Files to load in:
# HighV-Quest results for top2 clones:
top2_file <- "sra/800MiXCR/SUMMARY/imgt_highvquest_results/20220515_top2clones/vquest_airr.tsv"
# HighV-Quest results for CAP3 assembled contigs:
cap3_file <- "sra/800MiXCR/SUMMARY/imgt_highvquest_results/20220515_cap3contigs/vquest_airr.tsv"
# Main results:
results_file <- "sra/800MiXCR/SUMMARY/mixcrstats_topclones_blast_combined_correctedFractions_withCategorizations.txt"
# Original MiXCR ratios:
mixcr_ratio = "sra/800MiXCR/SUMMARY/mixcr_ratio.csv"
# Output directory:
out_dir = "sra/800MiXCR/SUMMARY/"

# Treat as constant start position:
CONSTANT_START <- 390

# Read in files:
top2_df <- as_tibble(read_tsv(top2_file))
cap3_df <- as_tibble(read_tsv(cap3_file))
results_df <- as_tibble(read_tsv(results_file , na=c("", "NotAvailable", "NA")))
mixcr_df <- as_tibble(read_csv(mixcr_ratio))

# Get main light chain ID:
samples_light <- mixcr_df %>% 
  select(sample, lc_with_most)


# Split out columns:
# Top2 dataframe clean up:
top2_df <- top2_df %>% 
  separate(sequence_id, c("sample", "chain", "cloneID", "cloneRank", "subseq"), sep="_")
top2_df <- top2_df %>% 
  mutate(cloneID = gsub("cID-", "", cloneID),
         cloneRank = gsub("cRank-", "", cloneRank)) %>% 
  mutate_at(c("cloneID", "cloneRank"), as.numeric)
# CAP3 df clean up:
cap3_df <- cap3_df %>% 
  separate(sequence_id, c("sample", "chain", "subseq"))
cap3_df$cloneID <- -999
cap3_df$cloneRank <- -999
# Combine together:
imgt_combined <- bind_rows(top2_df, cap3_df)

# Convert subsequence to factor:
imgt_combined$subseq <- as_factor(imgt_combined$subseq)




# Get part of sequence that is 5'UTR+Leader, V+J, Constant:
# Start of V: "v_sequence_start"
# End of J: "j_sequence_end"
imgt_combined$fiveprime <- str_sub(imgt_combined$sequence, start=1, end=imgt_combined$v_sequence_start-1)
imgt_combined$threeprime <- str_sub(imgt_combined$sequence, start=imgt_combined$j_sequence_end+1, end=nchar(imgt_combined$sequence))
imgt_combined$sequence_length <- nchar(imgt_combined$sequence)
imgt_combined$fiveprime_length <- nchar(imgt_combined$fiveprime)
imgt_combined$threeprime_length <- nchar(imgt_combined$threeprime)
imgt_combined$V_length <- nchar(imgt_combined$v_sequence_alignment)
imgt_combined$J_length <- nchar(imgt_combined$j_sequence_alignment)
# Fill NA's with O
imgt_combined <- imgt_combined %>% 
  mutate(sequence_length = ifelse(is.na(sequence_length), 0, sequence_length),
         fiveprime_length = ifelse(is.na(fiveprime_length), 0, fiveprime_length),
         threeprime_length = ifelse(is.na(threeprime_length), 0, threeprime_length),
         V_length = ifelse(is.na(V_length), 0, V_length),
         J_length = ifelse(is.na(J_length), 0, J_length))

# Get leading and trailing gap:
# TODO: could have probably skipped swapping it to spaces
imgt_combined$sequence_alignment_withspaces <- str_replace_all(imgt_combined$sequence_alignment, "\\.", " ")
imgt_combined$sequence_alignment_withspaces_reversed <- stringi::stri_reverse(imgt_combined$sequence_alignment_withspaces)
imgt_combined$leading_gap <- stringr::str_count(imgt_combined$sequence_alignment_withspaces, "\\G ") 
imgt_combined$ending_gap <- stringr::str_count(imgt_combined$sequence_alignment_withspaces_reversed, "\\G ")
# Do same with V:
imgt_combined$v_sequence_alignment_withspaces <- str_replace_all(imgt_combined$v_sequence_alignment, "\\.", " ")
imgt_combined$v_sequence_alignment_withspaces_reversed <- stringi::stri_reverse(imgt_combined$v_sequence_alignment_withspaces)
imgt_combined$v_leading_gap <- stringr::str_count(imgt_combined$v_sequence_alignment_withspaces, "\\G ") 
imgt_combined$v_ending_gap <- stringr::str_count(imgt_combined$v_sequence_alignment_withspaces_reversed, "\\G ")
# Do same with J:
imgt_combined$j_sequence_alignment_withspaces <- str_replace_all(imgt_combined$j_sequence_alignment, "\\.", " ")
imgt_combined$j_sequence_alignment_withspaces_reversed <- stringi::stri_reverse(imgt_combined$j_sequence_alignment_withspaces)
imgt_combined$j_leading_gap <- stringr::str_count(imgt_combined$j_sequence_alignment_withspaces, "\\G ") 
imgt_combined$j_ending_gap <- stringr::str_count(imgt_combined$j_sequence_alignment_withspaces_reversed, "\\G ")

# Get maximum V-J length by chain:
imgt_combined$VJ_length <- nchar(imgt_combined$sequence_alignment)
length_summary <- imgt_combined %>% 
  group_by(chain) %>% 
  summarise(min = min(VJ_length, na.rm=TRUE),
            mean = mean(VJ_length, na.rm=TRUE),
            max = max(VJ_length, na.rm=TRUE), 
            n = n(),
            minV = min(V_length, na.rm=TRUE),
            meanV = mean(V_length, na.rm=TRUE),
            maxV = max(V_length, na.rm=TRUE), 
            minJ = min(J_length, na.rm=TRUE),
            meanJ = mean(J_length, na.rm=TRUE),
            maxJ = max(J_length, na.rm=TRUE))


# TODO except this math doesn't work because J numbering starts again at 1...
# Start of segment is 1 - length(upstream) if length(upstream) > 0, else Start of segment is v_germline_start if not NA, else start of segment is j_germline_start
# End of segment is CONSTANT_START + length(constant) if lenght(constant) > 0, else End of segment is j_germline_end if not NA, else End of segment is v_germline_end
imgt_combined$master_start <- NA_real_
imgt_combined$master_end <- NA_real_
# TODO: do these if else statements make sense?
# TODO: test ending gap...
# TODO: am I Handling the J leading gap stuff correctly?
imgt_combined$master_start <- if_else(!is.na(imgt_combined$fiveprime_length) & imgt_combined$fiveprime_length > 0, 1-imgt_combined$fiveprime_length, 
                                      if_else(!is.na(imgt_combined$v_germline_start), imgt_combined$v_germline_start + imgt_combined$v_leading_gap, 
                                              if_else(!is.na(imgt_combined$j_germline_start), CONSTANT_START-imgt_combined$J_length+imgt_combined$j_leading_gap, NA_real_)))
imgt_combined$master_end <- if_else(!is.na(imgt_combined$threeprime_length) & imgt_combined$threeprime_length>0, CONSTANT_START+imgt_combined$threeprime_length, 
                                    ifelse(!is.na(imgt_combined$j_germline_end), CONSTANT_START- + (-39 + imgt_combined$j_germline_end) - imgt_combined$j_ending_gap,
                                           ifelse(!is.na(imgt_combined$v_germline_end), imgt_combined$v_germline_end - imgt_combined$v_ending_gap, NA_real_)))


# Just look at the most relevant pieces of information:
imgt_smalldf <- imgt_combined %>% 
  dplyr::select(sample, chain, subseq, sequence, sequence_alignment, sequence_alignment_withspaces, fiveprime, threeprime, 
         v_sequence_start, v_sequence_end, j_sequence_start, j_sequence_end, 
         v_germline_start, v_germline_end, j_germline_start, j_germline_end, 
         sequence_length, fiveprime_length, threeprime_length,
         VJ_length, V_length, J_length, 
         master_start, master_end, leading_gap, ending_gap,
         v_leading_gap, v_ending_gap, j_leading_gap, j_ending_gap) %>% 
  filter(chain != "IGH")

summary(imgt_smalldf)

# ------------------------------------------------------------------------------
# Select on applicable columns:
results_df2 <- results_df %>% 
  select(sample, CATEGORY, lc_with_most, first_subseq_count, blast_assembled_string)
# Combine with BLAST results:
results_df2 <- full_join(results_df2, imgt_combined, by=c("sample"="sample"))
# Create factors where needed:
results_df2$CATEGORY <- factor(results_df2$CATEGORY, levels=c("Category1_ReadyToGo", 
                                                              "Category2_Top2Together", 
                                                              "Category3_MajorClonepresentbelowthreshold", 
                                                              "Category4_uncategorizedsamples"))
results_df3 <- results_df2 %>% 
  mutate(across(c(sample, lc_with_most, subseq), as.factor))

# Keep only where dominant LC equals the light chain we are working with from IMGT:
results_df3 <- results_df3 %>% 
  filter(lc_with_most == chain)

print(sprintf("Number of samples from MiXCR stats: %d", length(mixcr_df$sample)))
print(sprintf("Number of unique samples: %d", length(unique(results_df3$sample))))
# NOTE: MMRF121622 IGK and IGL clones are all too low COUNT to be considered
print("Missing samples (if any):")
print(setdiff(mixcr_df$sample, unique(results_df3$sample)))
print(sprintf("Number of unique samples+segments: %d", nrow(unique(results_df3[, c("sample", "subseq")]))))

# Order:
results_df3 <- results_df3 %>% 
  group_by(sample) %>% 
  mutate(firstStart = min(master_start)) %>% 
  ungroup()

# TODO what am I doing (see Prof. Morgan's script)

# For coloring purpose:
# TODO do I need this? results_df3$sample <- fct_reorder(results_df3$sample, results_df3$firstStart)

# For coloring purpose, cRank and subseq combined column:
results_df3$cRank_and_segment <- paste(results_df3$cloneRank, results_df3$subseq, sep="-")
results_df3$cRank_and_segment <- factor(results_df3$cRank_and_segment)


# ------------------------------------------------------------------------------
# IMGT length distribution:
imgt_combined %>% 
  ggplot(aes(x=VJ_length)) +
  geom_histogram(aes(fill=chain), binwidth=10) +
  xlab("Length of IMGT alignment") +
  ylab("Count") +
  stat_bin(geom="text", colour="black", size=2, binwidth=10,
           aes(label=..count.., group=chain, y=1.1*(..count..))) +
  # stat_bin(binwidth=25, geom="text", colour="black", size=3.5,
  #          aes(label=..count.., group=chain, y=1.2*(..count..))) +
  ggtitle("Length of alignment string produced by IMGT") +
  facet_grid(rows=vars(chain))
# Look across subseq:
imgt_combined %>% 
  ggplot(aes(x=VJ_length)) +
  geom_histogram(aes(fill=chain), binwidth=10) +
  xlab("Length of IMGT alignment") +
  ylab("Count") +
  # stat_bin(binwidth=25, geom="text", colour="black", size=3.5,
  #          aes(label=..count.., group=chain, y=1.2*(..count..))) +
  ggtitle("Length of alignment string produced by IMGT") +
  facet_grid(rows=vars(subseq), cols=vars(chain), scales="free")

# IMGT length:
imgt_combined %>% 
  ggplot(aes(x=VJ_length, y=sequence_length, color=chain)) +
  geom_point(size=1) +
  xlab("Length of IMGT alignment") +
  ylab("Sequence length") +
  # stat_bin(binwidth=25, geom="text", colour="black", size=3.5,
  #          aes(label=..count.., group=chain, y=1.2*(..count..))) +
  ggtitle("Length of alignment string produced by IMGT") # +
  # facet_grid(rows=vars(subseq), cols=vars(chain), scales="free")
imgt_combined %>% 
  ggplot(aes(x=VJ_length, y=sequence_length, color=chain)) +
  geom_point(size=0.5) +
  xlab("Length of IMGT alignment") +
  ylab("Sequence length") +
  # stat_bin(binwidth=25, geom="text", colour="black", size=3.5,
  #          aes(label=..count.., group=chain, y=1.2*(..count..))) +
  ggtitle("Length of alignment string produced by IMGT")  +
  facet_wrap(~subseq, scales="free")

# ------------------------------------------------------------------------------
# Look at alignment position
results_df3 %>% 
  # slice_min(sample, n=100) %>%   # TODO remove
  ggplot(aes(sample, color=subseq, group=sample)) +
  scale_x_discrete(breaks=NULL) +
  scale_y_discrete(drop=TRUE) +  # TODO does this do what I want (don't want samples not present taking space?
  labs(y="Sequence position")  +
  coord_flip() +  # TODO is this needed?
  geom_linerange(aes(ymin=master_start, ymax=master_end), size=0.2) +
  facet_grid(lc_with_most ~ CATEGORY, scales="free")

# Lets skip the Grid:
for (cat in levels(results_df3$CATEGORY)){
  for (chain in levels(results_df3$lc_with_most)){
    temp_fig <- results_df3 %>% 
      filter(CATEGORY==cat & lc_with_most==chain & !(!CATEGORY=="Category2_Top2Together" & first_subseq_count<2)) %>% 
      # slice_min(sample, n=100) %>%   # TODO remove
      ggplot(aes(sample, color=subseq, group=sample)) +
      # scale_x_discrete(breaks=NULL) +
      # scale_y_discrete(drop=TRUE) +  # TODO does this do what I want (don't want samples not present taking space?
      labs(y="Sequence position")  +
      coord_flip() +  # TODO is this needed?
      geom_linerange(aes(ymin=master_start, ymax=master_end), 
                     size=ifelse(cat=="Category1_ReadyToGo", 3, 3),
                     alpha=ifelse(cat=="Category2_Top2Together", 0.5, 0.5)) +
      ggtitle(sprintf("%s %s", cat, chain))
    # show(temp_fig)
  }
}


# Lets just do category 2 and chains separate:
for (ch in c("IGK", "IGL")){
  temp <- results_df3 %>%  
    filter(CATEGORY=="Category2_Top2Together" & chain==ch) %>% 
    ggplot(aes(sample, color=cRank_and_segment, group=sample)) +
    # scale_x_discrete(breaks=NULL) +
    # scale_y_discrete(drop=TRUE) +  # TODO does this do what I want (don't want samples not present taking space?
    labs(y="Sequence position")  +
    coord_flip() +  # TODO is this needed?
    geom_linerange(aes(ymin=master_start, ymax=master_end), 
                   size=2,
                   alpha=0.5) +
    ggtitle(paste("Category 2 ", ch, sep="")) +
    # theme(axis.text.y=element_text(size=4)) +
    facet_grid(~ cloneRank, scales="free_y") +
    geom_hline(yintercept = c(1, CONSTANT_START), colour="grey60")
  # show(temp)
}

table(results_df3 %>%  
          filter(CATEGORY=="Category2_Top2Together") %>% 
          select(cRank_and_segment))


# -----------------------------------------------------------------------------
# Filter on contigs for CATEGORY 2, and top clone for everything else:
filter_top <- results_df3 %>% 
  filter((cloneRank==0 & CATEGORY!="Category2_Top2Together") | 
           (cloneRank==-999 & CATEGORY=="Category2_Top2Together"))

for (cat in levels(filter_top$CATEGORY)){
  for (ch in c("IGK", "IGL")){
    temp <- filter_top %>% 
      filter(CATEGORY==cat & chain==ch) %>% 
      ggplot(aes(sample, color=cRank_and_segment, group=sample)) +
      # scale_x_discrete(breaks=NULL) +
      # scale_y_discrete(drop=TRUE) +  # TODO does this do what I want (don't want samples not present taking space?
      labs(y="Sequence position")  +
      coord_flip() +  # TODO is this needed?
      geom_linerange(aes(ymin=master_start, ymax=master_end), 
                     size=ifelse(cat=="Category1_ReadyToGo", 0.5, 2),
                     alpha=0.5) +
      ggtitle(paste(cat, ch, sep=" ")) +
      geom_hline(yintercept = c(1, CONSTANT_START), colour="grey60")
    # show(temp)
  }
}

# Focus on those that have a V_start > 1 OR multiple segments:
for (cat in levels(filter_top$CATEGORY)){
  for (ch in c("IGK", "IGL")){
    temp <- filter_top %>% 
      filter(CATEGORY==cat & chain==ch) %>% 
      filter(master_start>1 | first_subseq_count > 1 | master_end < CONSTANT_START) %>% 
      ggplot(aes(sample, color=cRank_and_segment, group=sample)) +
      # scale_x_discrete(breaks=NULL) +
      # scale_y_discrete(drop=TRUE) +  # TODO does this do what I want (don't want samples not present taking space?
      labs(y="Sequence position")  +
      coord_flip() +  # TODO is this needed?
      geom_linerange(aes(ymin=master_start, ymax=master_end), 
                     size=3,
                     alpha=0.5) +
      ggtitle(paste(cat, ch, sep=" ")) +
      geom_hline(yintercept = c(1, CONSTANT_START), colour="grey60")
    # show(temp)
  }
}

# -----------------------------------------------------------------------------
# Create a modified dataframe, with sample IDs indexed by group:
filter_top2_IGK <- filter_top %>% 
  filter(chain == "IGK") %>% 
  arrange(CATEGORY, master_start) %>% 
  mutate(sample_ordered=factor(sample, levels=unique(sample[order(dplyr::desc(CATEGORY), master_start, dplyr::desc(master_end), sample)]))) %>% 
  mutate(sample_graphing=as.numeric(sample_ordered))
filter_top2_IGL <- filter_top %>% 
  filter(chain == "IGL") %>% 
  arrange(CATEGORY, master_start) %>% 
  mutate(sample_ordered=factor(sample, levels=unique(sample[order(dplyr::desc(CATEGORY), master_start, dplyr::desc(master_end), sample)]))) %>% 
  mutate(sample_graphing=as.numeric(sample_ordered))
filter_top2 <- bind_rows(filter_top2_IGK, filter_top2_IGL) %>% 
  mutate(sample_graphing = factor(sample_graphing))

# Save for graphing:
write.table(filter_top2, file=paste(out_dir, "imgt_results_for_graphing.txt", sep=""), 
            sep="\t", row.names=FALSE)


# Counts by group:
temp_counts <- filter_top2 %>% 
  select(sample, chain, CATEGORY)
temp_counts <- distinct(temp_counts)
table(temp_counts$chain, temp_counts$CATEGORY)

# Summary coverage figure
filter_top %>% 
  arrange(CATEGORY, master_start) %>% 
  mutate(sample_ordered=factor(sample, levels=unique(sample[order(CATEGORY, master_start, dplyr::desc(master_end), sample)]))) %>% 
  ggplot(aes(sample_ordered, color=CATEGORY)) +
  # scale_x_discrete(breaks=NULL) +
  # scale_y_discrete(drop=TRUE) +  # TODO does this do what I want (don't want samples not present taking space?
  labs(x="Sample", y="Sequence position")  +
  coord_flip() +  # TODO is this needed?
  geom_linerange(aes(ymin=master_start, ymax=master_end), 
                 size=0.25,
                 alpha=0.5) +
  ggtitle("Coverage") +
  geom_hline(yintercept = c(1, CONSTANT_START), colour="grey60") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  # facet_grid(CATEGORY~chain, scales="free") +
  geom_text(aes(x=0, y=0, label="Start V-J"), colour="black", hjust=-0.08, vjust=-0.8, size=3) +
  facet_wrap(~chain, ncol=2, scales="free")
# Do labeleing differently:
df_line <- data.frame(
  y_intercept = c(1, CONSTANT_START),
  label = c("Vstart", "Cstart"),
  linetype = c("dashed", "solid")
)
summary_coverage_fig <- filter_top2 %>% 
  # arrange(CATEGORY, master_start) %>% 
  # mutate(sample_ordered=factor(sample, levels=unique(sample[order(dplyr::desc(CATEGORY), master_start, dplyr::desc(master_end), sample)]))) %>% 
  ggplot(aes(sample_graphing, color=CATEGORY)) +
  # scale_x_discrete(breaks=NULL) +
  # scale_y_discrete(drop=TRUE) +  # TODO does this do what I want (don't want samples not present taking space?
  # labs(x="Sample", y="Sequence position")  +
  xlab("Sample") +
  ylab("Sequence position") +
  coord_flip() +  # TODO is this needed?
  geom_linerange(aes(ymin=master_start, ymax=master_end), 
                 size=1,  # Was 0.2, ifelse(chain=="IGK", 0.001, 1)
                 alpha=1) +
  ggtitle("Coverage") +
  geom_hline(data=df_line, aes(yintercept=y_intercept, labels=label), colour="gray20", linetype="dashed") +
  scale_y_continuous(sec.axis = sec_axis(~., breaks = df_line$y_intercept, labels = df_line$label)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.placement="outside", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="black", colour="black")) +  # TODO switch back to element blank?
  # facet_grid(CATEGORY~chain, scales="free") +
  facet_wrap(~chain, ncol=2, scales="free_x")
# show(summary_coverage_fig)
ggsave(paste(out_dir, "summary_coverage_fig.svg", sep=""), summary_coverage_fig,  width=16, height=9, units="in", dpi=600)
# Larger figure:
summary_coverage_fig <- filter_top %>% 
  arrange(CATEGORY, master_start) %>% 
  mutate(sample_ordered=factor(sample, levels=unique(sample[order(dplyr::desc(CATEGORY), master_start, dplyr::desc(master_end), sample)]))) %>% 
  ggplot(aes(sample_ordered, color=CATEGORY)) +
  # scale_x_discrete(breaks=NULL) +
  # scale_y_discrete(drop=TRUE) +  # TODO does this do what I want (don't want samples not present taking space?
  # labs(x="Sample", y="Sequence position")  +
  xlab("Sample") +
  ylab("Sequence position") +
  coord_flip() +  # TODO is this needed?
  geom_linerange(aes(ymin=master_start, ymax=master_end), 
                 size=1,  # Was 0.2, ifelse(chain=="IGK", 0.001, 1)
                 alpha=1) +
  ggtitle("Coverage") +
  geom_hline(data=df_line, aes(yintercept=y_intercept, labels=label), colour="gray20", linetype="dashed") +
  scale_y_continuous(sec.axis = sec_axis(~., breaks = df_line$y_intercept, labels = df_line$label)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.placement="outside", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="black", colour="black")) +  # TODO switch back to element blank?
  # facet_grid(CATEGORY~chain, scales="free") +
  facet_wrap(~chain, ncol=2, scales="free")
# show(summary_coverage_fig)
ggsave(paste(out_dir, "summary_coverage_fig_big.svg", sep=""), summary_coverage_fig,  width=24, height=13.5, units="in", dpi=1200)
# Just part of the figure:
summary_coverage_fig <- filter_top %>% 
  arrange(CATEGORY, master_start) %>% 
  mutate(sample_ordered=factor(sample, levels=unique(sample[order(dplyr::desc(CATEGORY), master_start, dplyr::desc(master_end), sample)]))) %>% 
  filter(chain=="IGK" & CATEGORY=="Category1_ReadyToGo") %>% 
  ggplot(aes(sample_ordered, color=CATEGORY)) +
  # scale_x_discrete(breaks=NULL) +
  # scale_y_discrete(drop=TRUE) +  # TODO does this do what I want (don't want samples not present taking space?
  # labs(x="Sample", y="Sequence position")  +
  xlab("Sample") +
  ylab("Sequence position") +
  coord_flip() +  # TODO is this needed?
  geom_linerange(aes(ymin=master_start, ymax=master_end), 
                 size=2,  # Was 0.2, ifelse(chain=="IGK", 0.001, 1)
                 alpha=1) +
  ggtitle("Coverage") +
  geom_hline(data=df_line, aes(yintercept=y_intercept, labels=label), colour="gray20", linetype="dashed") +
  scale_y_continuous(sec.axis = sec_axis(~., breaks = df_line$y_intercept, labels = df_line$label)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.placement="outside", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="black", colour="black")) +  # TODO switch back to element blank?
  # facet_grid(CATEGORY~chain, scales="free") +
  facet_wrap(~chain, ncol=2, scales="free")
# show(summary_coverage_fig)
ggsave(paste(out_dir, "summary_coverage_fig_partial_IGK.svg", sep=""), summary_coverage_fig,  width=15, height=16, units="in", dpi=1200)

# TODO: add length of longest segment to order?

# -----------------------------------------------------------------------------
# 3' length:
constant_length_fig <- filter_top %>% 
  group_by(sample) %>% 
  filter(master_end == max(master_end, na.rm=TRUE)) %>% 
  # ggplot(aes(x=threeprime_length, fill=CATEGORY)) +
  ggplot(aes(x=threeprime_length, fill=chain)) +
  geom_histogram(binwidth=10) +
  ggtitle("Constant Region Length") +
  xlab("Length (from 3' segment)") +
  # stat_bin(binwidth=10, geom="text", colour="black", size=3.5,
  #          aes(label=..count.., group=chain, y=(..count..))) +
  # facet_wrap(CATEGORY~chain, ncol=2, scales="free")
  facet_wrap(~chain, ncol=2, scales="free") 
# show(constant_length_fig)
ggsave(paste(out_dir, "constant_length_fig.svg", sep=""), constant_length_fig,  width=8, height=6, units="in", dpi=600)

# Get average constant region:
filter_top %>% 
  group_by(sample) %>% 
  filter(master_end == max(master_end, na.rm=TRUE)) %>% 
  ungroup() %>% 
  summarise(Mean = mean(threeprime_length, na.rm=TRUE))
# By group:
filter_top %>% 
  group_by(sample) %>% 
  filter(master_end == max(master_end, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(chain) %>% 
  summarise(Mean = mean(threeprime_length, na.rm=TRUE))

# -----------------------------------------------------------------------------
# Determine number with complete VDJ (based on longest segment)
longest_filter_top <- filter_top %>% 
  group_by(sample) %>% 
  top_n(1, sequence_length) %>% 
  ungroup()
table(longest_filter_top$CATEGORY, longest_filter_top$complete_vdj)
table(longest_filter_top$lc_with_most, longest_filter_top$complete_vdj)
table(longest_filter_top$chain, longest_filter_top$complete_vdj, longest_filter_top$CATEGORY)
longest_filter_top %>% 
  # filter(complete_vdj == TRUE) %>% 
  count(complete_vdj, CATEGORY, chain, sort=FALSE)

# -----------------------------------------------------------------------------
# Take a lot at Category 2 and number of segments per group:
filter_top_cat2 <- filter_top %>% 
  filter(CATEGORY=="Category2_Top2Together")
count_cat2_pieces <- filter_top_cat2 %>% count(sample)
cat2_multiplepieces <- count_cat2_pieces %>% 
  filter(n>1) %>% 
  pull(sample)
# Look at the multiple pieces samples:
cat2_multiplepieces_df <- filter_top_cat2 %>% 
  filter(sample %in% cat2_multiplepieces)
