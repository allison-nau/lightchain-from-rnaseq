# Calculate MiXCR ratios

# Import libraries:
library(tidyverse)
library(ggpubr)

# Get command line arguments:
# (giant mixcr combined tsv file) (location of accounting.txt file made using 
# qacct) (directory to output results to)
args = commandArgs(trailingOnly=TRUE)
read_in = args[1]
accounting_in = args[2]
out_dir = args[3]
# For running code without command args:
# read_in = "sra/800MiXCR/SUMMARY/mixcr_all_clonotypes_v0.txt"
# accounting_in = "sra/800MiXCR/SUMMARY/job_account.csv"
# stdout_in = "sra/800MiXCR/SUMMARY/processed_stdout_wide_NOTrestricted.csv"
# out_dir = "sra/800MiXCR/SUMMARY/"


# Flag thresholds:
# Minimum light chain length:
MINLENGTH <- 500
# Minimum dominant clone freq:
MINDOMCLONE <- 0.75
# Minimum dominant light chain freq:
MINDOMLC <- 0.9

# Scientific notation options:
options(scipen=10000)

# Read in results:
results <- as_tibble(read_tsv(read_in, na=c(""),
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

# Group by sample and chain
my_summary <- results %>% group_by(sample, Chain) %>% 
  summarize(total_count=sum(cloneCount, na.rm=TRUE))
# my_summary

# Convert from long to wide:
my_summary <- my_summary %>% spread(Chain, total_count)

# Light Chain total:
my_summary <- my_summary %>%  rowwise() %>% 
  mutate(light = sum(IGK, IGL, na.rm=TRUE))

# Swap NA's with 0:
my_summary <- my_summary %>% replace_na(list(IGH=0, IGK=0, IGL=0, light=0))

# Get count of light chain that is dominant:
my_summary <- my_summary %>% 
  rowwise() %>% 
  mutate(lc_dom=max(IGK, IGL))

# IGK:IGL ratio & Light to heavy ratio:
my_summary <- my_summary %>% rowwise() %>% 
  mutate("lc_dom/light" = lc_dom/light, 
         "light/IGH" = light/IGH,
         "IGH/light" = IGH/light)

# Add all IGK, IGL, IGH:
my_summary <- my_summary %>%  rowwise() %>% 
  mutate(total = sum(IGK, IGL, IGH, na.rm=TRUE)) 

# Current summary:
my_summary

# ----------------------------------------------------------------------------
# What fraction of reads are in the dominant light chain
# Keep only certain columns
keep_col <- c("Chain", "cloneRank", "sample", "cloneId", "cloneCount", 
              "cloneFraction", "mixcr_length", "subseq1len", "subseq2len", 
              "subseq3len", "subseq4len", "subseq5len", "subseq6len", 
              "subseq7len", "subseq8len", "subseq9len", "maxSubLength")
dominant <- results[, keep_col]

# Remove heavy chain
dominant <- dominant %>% 
  filter(Chain != "IGH")
# Make chain factor:
dominant$Chain <- factor(dominant$Chain)
# Will finish getting the dominant after getting the max_sublength

# Get total number of clones:
number_of_clones_df <- dominant %>% count(sample)
number_of_clones_df <- number_of_clones_df %>% 
  rename(number_of_clones = n)

# Group by sample to get max sublength:
max_sublength <- dominant %>% 
  group_by(sample) %>% 
  arrange(desc(maxSubLength)) %>% 
  dplyr::slice(1)
# Keep only certain columns and rename:
max_sublength <- max_sublength[, c("Chain", "cloneRank", "sample", "cloneId", "cloneCount", 
                                   "cloneFraction", "maxSubLength")]
# Add prefix to column names:
colnames(max_sublength) <- paste("longest", colnames(max_sublength), sep="_")
# Remove prefix from sample name:
max_sublength <- max_sublength %>% 
  dplyr::rename("sample" = "longest_sample")

# Get the dominant light chain row:
dominant <- dominant %>% 
  group_by(sample) %>% 
  arrange(desc(cloneCount)) %>% 
  dplyr::slice(1)
# Add prefix to column names:
colnames(dominant) <- paste("domClone", colnames(dominant), sep="_")
# Remove prefix from sample name:
dominant <- dominant %>% 
  dplyr::rename("sample" = "domClone_sample")


# ----------------------------------------------------------------------------
# Combine summary with dominant and max_sublength:
comb_summary <- full_join(my_summary, dominant, by="sample")
comb_summary <- full_join(comb_summary, max_sublength, by="sample")
comb_summary <- full_join(comb_summary, number_of_clones_df, by="sample")

# Calculate what fraction of light chain is the dominant clone:
comb_summary <- comb_summary %>% rowwise() %>% 
  mutate("dom_clone/light" = domClone_cloneCount/light)

# Determine if IGK or IGL should have been dominant based on total clone count:
comb_summary$lc_with_most = NA
comb_summary$lc_with_most[comb_summary$IGK > comb_summary$IGL] <- "IGK"
comb_summary$lc_with_most[comb_summary$IGK < comb_summary$IGL] <- "IGL"

# Does dominant chain by most frequent clone match by total number from that chain?
comb_summary$domClone_domChain <- comb_summary$lc_with_most == comb_summary$domClone_Chain

# Rearrange columns:
comb_summary <- comb_summary %>% relocate(lc_with_most, .after=light)
comb_summary <- comb_summary %>% relocate(`dom_clone/light`, .before=domClone_Chain)
comb_summary <- comb_summary %>% relocate(`domClone_cloneCount`, .before=`dom_clone/light`)

# Across samples, what proportion have IGK as dominant light chain?
t <- table(comb_summary$lc_with_most)
print("Number IGK or IGL dominant:")
t
print("Proportion IGK or IGL dominant:")
prop.table(t)

# Across samples, how many have the top clonotype as IGK vs IGL?
t2 <- table(comb_summary$domClone_Chain)
print("Number IGK or IGL top clone:")
t2
print("Proportion IGK or IGL tope clone:")
prop.table(t2)

# Comparing IGK vs IGL total count vs top clone
table(comb_summary$lc_with_most, comb_summary$domClone_Chain)


# Make lc with most factor:
comb_summary$lc_with_most <- factor(comb_summary$lc_with_most)

# get ones that Longest light chain subsequence length != dominant LC Clone Longest subsequent length:
print("Samples where the longest light chain is not the dominant LC clone:")
samples_mainclone_notlongest <- comb_summary$sample[comb_summary$domClone_maxSubLength != comb_summary$longest_maxSubLength]
samples_mainclone_notlongest
print(sprintf("Number of samples dominant clone is not longest: %d", length(samples_mainclone_notlongest)))
# get ones that dom_clone/light < 90%
print(paste("Samples where the dominant light chain clone is not at least ", MINDOMCLONE*100, "% of light chain sequences:", sep=""))
# samples_mainclone_infrequent <- comb_summary$sample[comb_summary$`dom_clone/light` <= 0.85]
# samples_mainclone_infrequent <- comb_summary$sample[comb_summary$`dom_clone/light` <= 0.95]
samples_mainclone_infrequent <- comb_summary$sample[comb_summary$`dom_clone/light` < MINDOMCLONE]
samples_mainclone_infrequent
print(sprintf("Number of samples dominant light chain clone is not at least %d%% of light chain sequences: %d", MINDOMCLONE*100, length(samples_mainclone_infrequent)))
# Get ones where dominant LC clone < 500 nt:
print(sprintf("Samples where the longest segment in the dominant light chain clone is <%dnt", MINLENGTH))
samples_mainclone_short <- comb_summary$sample[comb_summary$domClone_maxSubLength < MINLENGTH]
samples_mainclone_short
print(sprintf("Number of samples dominant light chain is not at least %d nt: %d", MINLENGTH, length(samples_mainclone_short)))
# Where dominant light chain isn't at least 95% of light chain 
print(sprintf("Samples where the dominant light chain is not at least %d%% of light chain sequences:", MINDOMLC*100))
print(sprintf("i.e. samples where the non-dominant light chain type (IGK or IGL) is more than %d%% of light chain sequences", round((1-MINDOMLC)*100)))
# samples_mainclone_infrequent <- comb_summary$sample[comb_summary$`dom_clone/light` <= 0.85]
# samples_mainclone_infrequent <- comb_summary$sample[comb_summary$`dom_clone/light` <= 0.95]
samples_mainchain_infrequent <- comb_summary$sample[comb_summary$`lc_dom/light` < MINDOMLC]
samples_mainchain_infrequent
print(sprintf("Number of samples dominant light chain is not at least %d%% of light chain sequences: %d", MINDOMLC*100, length(samples_mainchain_infrequent)))

# Samples to watch:
print("Samples to watch:")
samples_watch <- unique(c(samples_mainclone_notlongest, samples_mainclone_infrequent, samples_mainclone_short, samples_mainchain_infrequent))
samples_watch
print(sprintf("Number of samples to watch: %d", length(samples_watch)))
# Add to comb_summary:
comb_summary$mainclone_notlongest = "False"
comb_summary$mainclone_notlongest[comb_summary$sample %in% samples_mainclone_notlongest] = "True"
comb_summary$mainclone_lowfreq = "False"
comb_summary$mainclone_lowfreq[comb_summary$sample %in% samples_mainclone_infrequent] = "True"
comb_summary$mainchain_lowfreq = "False"
comb_summary$mainclone_lowfreq[comb_summary$sample %in% samples_mainchain_infrequent] = "True"
comb_summary$mainclone_short = "False"
comb_summary$mainclone_short[comb_summary$sample %in% samples_mainclone_short] = "True"
comb_summary$lc_flag = "OK"
comb_summary$lc_flag[comb_summary$mainclone_notlongest == "True" & comb_summary$mainclone_lowfreq ==  "False"] = "NotLongest"
comb_summary$lc_flag[comb_summary$mainclone_lowfreq == "True" & comb_summary$mainclone_notlongest == "False"] = "LowFreq"
comb_summary$lc_flag[comb_summary$mainclone_notlongest == "True" & comb_summary$mainclone_lowfreq == "True"] = "Both"
comb_summary$mainclone_notlongest <- factor(comb_summary$mainclone_notlongest)
comb_summary$mainclone_lowfreq <- factor(comb_summary$mainclone_lowfreq)
comb_summary$lc_flag = factor(comb_summary$lc_flag, levels=c("OK", "NotLongest", "LowFreq", "Both"))


# Get dominant IGH:
dom_igh <- results %>% 
  filter(Chain=="IGH", cloneRank==0) %>% 
  select(sample, cloneCount)
dom_igh <- dom_igh %>% 
  dplyr::rename(IGHdomCount = cloneCount)
# Add to comb_summary:
comb_summary <- full_join(comb_summary, dom_igh, by="sample")
# Fill NA with 0:
comb_summary <- comb_summary %>% replace_na(list(IGHdomCount=0))
# What fraction is dominant count out of total IGH?
comb_summary$domIghFrac <- comb_summary$IGHdomCount / comb_summary$IGH



# ----------------------------------------------------------------------------
# Get the dominant heavy chain row:
dominantH <- results[, keep_col]
# Remove heavy chain
dominantH <- dominantH %>% 
  filter(Chain == "IGH")
# Make chain factor:
dominantH$Chain <- factor(dominantH$Chain)
# Will finish getting the dominant after getting the max_sublength
# Group by sample to get max sublength:
maxH_sublength <- dominantH %>% 
  group_by(sample) %>% 
  arrange(desc(maxSubLength)) %>% 
  dplyr::slice(1)
# Keep only certain columns and rename:
maxH_sublength <- maxH_sublength[, c("Chain", "cloneRank", "sample", "cloneId", "cloneCount", 
                                     "cloneFraction", "maxSubLength")]
# Add prefix to column names:
colnames(maxH_sublength) <- paste("longestHeavy", colnames(maxH_sublength), sep="_")
# Remove prefix from sample name:
maxH_sublength <- maxH_sublength %>% 
  dplyr::rename("sample" = "longestHeavy_sample")
# Get the dominant light chain row:
dominantH <- dominantH %>% 
  group_by(sample) %>% 
  arrange(desc(cloneCount)) %>% 
  dplyr::slice(1)
# Add prefix to column names:
colnames(dominantH) <- paste("domCloneH", colnames(dominantH), sep="_")
# Remove prefix from sample name:
dominantH <- dominantH %>% 
  dplyr::rename("sample" = "domCloneH_sample")
# Combine summary with dominant and max_sublength:
comb_summary <- full_join(comb_summary, dominantH, by="sample")
comb_summary <- full_join(comb_summary, maxH_sublength, by="sample")

# ----------------------------------------------------------------------------

# Save summary:
write.csv(comb_summary, file=paste(out_dir, "mixcr_ratio.csv", sep=""), row.names=FALSE)

# Samples to double checK:
write.table(samples_watch, 
            paste(out_dir, "samples_to_manually_check.txt", sep=""), 
            row.names=FALSE, col.names=FALSE, quote=FALSE)
# ----------------------------------------------------------------------------


# plots: ---------------------------------------------------------------------
# Create some nicer names:
comb_summary$`Most Frequent LC` <- comb_summary$lc_with_most

# Violin plot of light/IGH ratio:
violin_light_igh <- comb_summary %>% 
  ggplot(aes(x=lc_with_most, y=`light/IGH`, color=lc_with_most)) +
  geom_violin() +
  geom_boxplot(width=0.1, outlier.shape=NA) +
  geom_jitter(data=subset(comb_summary, `light/IGH`<0.2 | `light/IGH` > 64), inherit.aes=TRUE, position=position_jitter(0.2), size=0.5) +
  ggtitle("light/IGH") +
  xlab("Dominant Light Chain") +
  scale_y_continuous(trans="log2", n.breaks=10, minor_breaks=NULL) +
  # scale_x_continuous(breaks=NULL) +
  theme(axis.ticks.x = element_blank())
# show(violin_light_igh)
ggsave(paste(out_dir, "violin_light_igh.svg", sep=""), violin_light_igh,  width=5, height=6, units="in")

# Violin plot lc_dom/light ratio
violin_lcdom_light <- comb_summary %>% 
  ggplot(aes(x=`Most Frequent LC`, y=`lc_dom/light`, color=`Most Frequent LC`)) +
  geom_violin() +
  geom_jitter(data=subset(comb_summary, `lc_dom/light`<0.95), inherit.aes=TRUE, position=position_jitter(0.2), size=0.5) +
  ggtitle("Dominant Light Chain Count / Total Light Chain Count") +
  xlab("Most Frequent Light Chain") +
  ylab("(Dominant Light Chain site count) / (Total Light Chain Count)") +
  # scale_x_continuous(breaks=NULL) +
  theme(axis.ticks.x = element_blank(),
        plot.title=element_text(size=12, hjust=0.4)) 
# show(violin_lcdom_light)
ggsave(paste(out_dir, "violin_lcdom_light.svg", sep=""), violin_lcdom_light,  width=5, height=6, units="in")

# Violin plot lc_dom/light ratio
violin_domclone_light <- comb_summary %>% 
  ggplot(aes(x=`Most Frequent LC`, y=`dom_clone/light`, color=`Most Frequent LC`)) +
  geom_violin() +
  geom_boxplot(width=0.1, outlier.shape=NA) +
  geom_jitter(data=subset(comb_summary, `dom_clone/light`<0.75), inherit.aes=TRUE, position=position_jitter(0.2), size=0.5) +
  # geom_text(data=subset(comb_summary, sample %in% samples_watch), 
  #           aes(label=sample), 
  #           hjust=-1.2, size=2.5) +
  ggtitle("Top Light Chain clone count / Total Light Chain Count") +
  xlab("Most Frequent Light Chain") +
  ylab("(Top Light Chain Clone Count) / (Total Light Chain Count)") +
  # scale_x_continuous(breaks=NULL) +
  theme(axis.ticks.x = element_blank(),
        plot.title=element_text(size=12, hjust=0.4))
# show(violin_domclone_light)
ggsave(paste(out_dir, "violin_domclone_light.svg", sep=""), violin_domclone_light,  width=5, height=6, units="in")

# Combined violin plots:
combined_violin <- ggarrange(violin_light_igh, violin_lcdom_light, violin_domclone_light, nrow=1)
# show(combined_violin)
ggsave(paste(out_dir, "combined_violin.svg", sep=""), combined_violin,  width=15, height=6, units="in")


# Look at `dom_clone/light` vs `lc_dom/light`
plot_clone_flags <- comb_summary %>% 
  ggplot(aes(x=`lc_dom/light`, y=`dom_clone/light`, color=domClone_maxSubLength, shape=mainclone_notlongest)) +
  geom_point() +
  ggtitle("Flagging Parameters") +
  xlab("Fraction dominant light chain out of all LC") +
  ylab("Fraction dominant light chain clone out of all LC") +
  scale_color_gradient(low="red", high="blue") +
  scale_shape_manual(values=c(20, 3))
# show(plot_clone_flags)
ggsave(paste(out_dir, "plot_clone_flags.svg", sep=""), plot_clone_flags,  width=9, height=6, units="in", dpi=600)
# plot_clone_flags_log <- comb_summary %>% 
#   ggplot(aes(x=`lc_dom/light`, y=`dom_clone/light`, color=domClone_maxSubLength, shape=mainclone_notlongest)) +
#   geom_point() +
#   ggtitle("Flagging Parameters") +
#   xlab("Fraction dominant light chain out of all LC") +
#   ylab("Fraction dominant light chain clone out of all LC") +
#   scale_color_gradient(low="red", high="blue") +
#   scale_shape_manual(values=c(20, 3)) +
#   scale_x_continuous(trans="log2") +
#   scale_y_continuous(trans="log2")
# show(plot_clone_flags_log)


# Total IGH count vs IGH dom ratio:
plot_ighRatio <- comb_summary %>% 
  replace_na(list(domIghFrac=0)) %>% 
  ggplot(aes(x=IGH, y=domIghFrac, color=lc_with_most, shape=lc_flag)) +
  geom_point() +
  # geom_text(data=subset(comb_summary, sample %in% samples_watch), 
  #           aes(x=IGH, y=domIghFrac, label=sample), 
  #           hjust=-0.1, size=2.5) +
  ggtitle("Total Heavy Chain Sequences\nvs\nDominant IGH Clone Ratio") +
  xlab("Total Heavy Chain Sequences") +
  ylab("Dominant Heavy Chain Sequence Ratio")
# show(plot_ighRatio)
ggsave(paste(out_dir, "plot_ighRatio.svg", sep=""), plot_ighRatio,  width=8, height=6, units="in")

# Graph maxSubLength for top clone vs maxsublength across all light:
plot_lightmaxsublength <- comb_summary %>% 
  ggplot(aes(x=longest_maxSubLength, y=domClone_maxSubLength)) +
  geom_point(aes(color=lc_with_most, shape=lc_flag)) +
  geom_text(data=subset(comb_summary, sample %in% samples_mainclone_notlongest), 
            aes(x=longest_maxSubLength, y=domClone_maxSubLength, label=sample), 
            hjust=-0.1, size=1) +
  geom_abline(intercept=0, slope=1) +
  ggtitle("Longest Light Chain Subsequence Length\nvs\nDominant Clone Longest Subsequence Length") +
  xlab("Longest Light Chain Subsequence Length (bp)") +
  ylab("Dominant LC Clone Longest Subsequence Length")
# show(plot_lightmaxsublength)
ggsave(paste(out_dir, "plot_lightmaxsublength.svg", sep=""), plot_lightmaxsublength,  width=8, height=6, units="in", dpi=600)
# Do against IG sequence count
plot_domlength <- comb_summary %>% 
  ggplot(aes(x=total, y=domClone_maxSubLength)) +
  geom_point(aes(color=lc_with_most, shape=lc_flag), size=2) +
  geom_text(data=subset(comb_summary, domClone_maxSubLength < 500), 
            aes(x=total, y=domClone_maxSubLength, label=sample), 
            hjust=-0.1, size=1) +
  geom_smooth(method="gam", se=FALSE, color="black") + # , label="GAM"
  geom_smooth(method="lm", se=FALSE, color="green") + # , label="lm"
  ggtitle("Total Ig vs\nDominant Clone Longest Subsequence Length") +
  xlab("Total Ig") +
  ylab("Dominant LC Clone Longest Subsequence Length") # +
  # legend()
# show(plot_domlength)
ggsave(paste(out_dir, "plot_domCloneSublength.svg", sep=""), plot_domlength,  width=8, height=6, units="in", dpi=600)
# Do against IG sequence count, no labels
plot_domlength_nolabels <- comb_summary %>% 
  ggplot(aes(x=total, y=domClone_maxSubLength)) +
  geom_point(aes(color=lc_with_most, shape=lc_flag), size=2) +
  geom_smooth(method="gam", se=FALSE, color="black") + # , label="GAM"
  geom_smooth(method="lm", se=FALSE, color="green") + # , label="lm"
  ggtitle("Total Ig vs\nDominant Clone Longest Subsequence Length") +
  xlab("Total Ig") +
  ylab("Dominant LC Clone Longest Subsequence Length") # +
# legend()
# show(plot_domlength_nolabels)
ggsave(paste(out_dir, "plot_domCloneSublength_nolabels.svg", sep=""), plot_domlength_nolabels,  width=8, height=6, units="in")
# Do against IG sequence count, all labels
plot_domlength_alllabels <- comb_summary %>% 
  ggplot(aes(x=total, y=domClone_maxSubLength)) +
  geom_point(aes(color=lc_with_most, shape=lc_flag), size=2) +
  geom_text(aes(x=total, y=domClone_maxSubLength, label=sample), 
            hjust=-0.1, size=1) +
  geom_smooth(method="gam", se=FALSE, color="black") + # , label="GAM"
  geom_smooth(method="lm", se=FALSE, color="green") + # , label="lm"
  ggtitle("Total Ig vs\nDominant Clone Longest Subsequence Length") +
  xlab("Total Ig") +
  ylab("Dominant LC Clone Longest Subsequence Length") # +
# show(plot_domlength_alllabels)
ggsave(paste(out_dir, "plot_domCloneSublength_alllabels.svg", sep=""), plot_domlength_alllabels,  width=8, height=6, units="in", dpi=600)
# Do against light count, no labels
plot_domlength_againstlight_nolabels <- comb_summary %>% 
  ggplot(aes(x=light, y=domClone_maxSubLength)) +
  geom_point(aes(color=lc_with_most, shape=lc_flag), size=2) +
  geom_smooth(method="gam", se=FALSE, color="black") + # , label="GAM"
  geom_smooth(method="lm", se=FALSE, color="green") + # , label="lm"
  ggtitle("Total Light Chain vs\nDominant Clone Longest Subsequence Length") +
  xlab("Total Light Chain") +
  ylab("Dominant LC Clone Longest Subsequence Length") # +
# legend()
# show(plot_domlength_againstlight_nolabels)
ggsave(paste(out_dir, "plot_domlength_againstlight_nolabels.svg", sep=""), plot_domlength_againstlight_nolabels,  width=8, height=6, units="in")
# Do against light count, no labels, simplified
plot_domlength_againstlight_nolabels_simplified <- comb_summary %>% 
  mutate(`Most Frequent LC` = lc_with_most) %>% 
  ggplot(aes(x=light, y=domClone_maxSubLength)) +
  geom_point(aes(color=`Most Frequent LC`, shape=`Most Frequent LC`), size=2) +
  # geom_smooth(method="gam", se=FALSE, color="black") + # , label="GAM"
  geom_smooth(method="lm", se=FALSE, color="black") + # , label="lm"
  ggtitle("Total Light Chain vs\nTop Clone Longest Subsequence Length") +
  xlab("Total Light Chain") +
  ylab("Top LC Clone Longest Subsequence Length") # +
# legend()
# show(plot_domlength_againstlight_nolabels)
ggsave(paste(out_dir, "plot_domlength_againstlight_nolabels_simplified.svg", sep=""), plot_domlength_againstlight_nolabels_simplified,  width=8, height=6, units="in")
# Do against IG sequence count, no labels, simplified
plot_domlength_againstlight_nolabels_simplified <- comb_summary %>% 
  mutate(`Most Frequent LC` = lc_with_most) %>% 
  ggplot(aes(x=total, y=domClone_maxSubLength)) +
  geom_point(aes(color=`Most Frequent LC`, shape=`Most Frequent LC`), size=2) +
  # geom_smooth(method="gam", se=FALSE, color="black") + # , label="GAM"
  geom_smooth(method="lm", se=FALSE, color="black") + # , label="lm"
  ggtitle("Total Ig vs\nTop Clone Longest Subsequence Length") +
  xlab("Total Ig") +
  ylab("Top LC Clone Longest Subsequence Length") # +
# legend()
# show(plot_domlength_nolabels)
ggsave(paste(out_dir, "plot_domlength_againstig_nolabels_simplified.svg", sep=""), plot_domlength_againstlight_nolabels_simplified,  width=8, height=6, units="in")


# Max light length, even if not dominant:
plot_lightmaxsublength2 <- comb_summary %>% 
  ggplot(aes(x=total, y=longest_maxSubLength)) +
  geom_point(aes(color=lc_with_most, shape=lc_flag), size=2) +
  geom_text(data=subset(comb_summary, longest_maxSubLength < MINLENGTH), 
            aes(x=total, y=longest_maxSubLength, label=sample), 
            hjust=-0.1, size=1) +
  geom_smooth(method="gam", se=FALSE, color="black") + # , label="GAM"
  geom_smooth(method="lm", se=FALSE, color="green") + # , label="lm"
  ggtitle("Total Ig vs\nLongest LC Subsequence Length") +
  xlab("Total Ig") +
  ylab("Longest LC Subsequence Length")
# show(plot_lightmaxsublength2)
ggsave(paste(out_dir, "plot_lightmaxsublength2.svg", sep=""), plot_lightmaxsublength2,  width=8, height=6, units="in", dpi=600)
# Max light length, even if not dominant, no labels::
plot_lightmaxsublength2_nolabels <- comb_summary %>% 
  ggplot(aes(x=total, y=longest_maxSubLength)) +
  geom_point(aes(color=lc_with_most, shape=lc_flag), size=2) +
  geom_smooth(method="gam", se=FALSE, color="black") + # , label="GAM"
  geom_smooth(method="lm", se=FALSE, color="green") + # , label="lm"
  ggtitle("Total Ig vs\nLongest LC Subsequence Length") +
  xlab("Total Ig") +
  ylab("Longest LC Subsequence Length")
# show(plot_lightmaxsublength2_nolabels)
ggsave(paste(out_dir, "plot_lightmaxsublength2_nolabels.svg", sep=""), plot_lightmaxsublength2_nolabels,  width=8, height=6, units="in")
# Max light length, even if not dominant, all labels:
plot_lightmaxsublength_alllabels2 <- comb_summary %>% 
  ggplot(aes(x=total, y=longest_maxSubLength)) +
  geom_point(aes(color=lc_with_most, shape=lc_flag), size=2) +
  geom_text(aes(x=total, y=longest_maxSubLength, label=sample), 
            hjust=-0.1, size=1) +
  geom_smooth(method="gam", se=FALSE, color="black") + # , label="GAM"
  geom_smooth(method="lm", se=FALSE, color="green") + # , label="lm"
  ggtitle("Total Ig vs\nLongest LC Subsequence Length") +
  xlab("Total Ig") +
  ylab("Longest LC Subsequence Length")
# show(plot_lightmaxsublength_alllabels2)
ggsave(paste(out_dir, "plot_lightmaxsublength_alllabels2.svg", sep=""), plot_lightmaxsublength_alllabels2,  width=8, height=6, units="in", dpi=600)

# Graph maxSubLength for top clone HEAVY vs maxsublength across all heavy:
plotH_maxsublength <- comb_summary %>% 
  ggplot(aes(x=longestHeavy_maxSubLength, y=domCloneH_maxSubLength)) +
  geom_point(aes(color=lc_with_most, shape=lc_flag), size=2) +
  geom_abline(intercept=0, slope=1) +
  ggtitle("Longest Heavy Chain Subsequence Length\nvs\nDominant Clone Longest Subsequence Length") +
  xlab("Longest Heavy Chain Subsequence Length (bp)") +
  ylab("Dominant IGH Clone Longest Subsequence Length")
# show(plotH_maxsublength)
ggsave(paste(out_dir, "plotH_maxsublength.svg", sep=""), plotH_maxsublength,  width=8, height=6, units="in")
# Graph maxSubLength for top clone HEAVY vs maxsublength across all heavy, all labels:
plotH_maxsublength_alllabels <- comb_summary %>% 
  ggplot(aes(x=longestHeavy_maxSubLength, y=domCloneH_maxSubLength)) +
  geom_point(aes(color=lc_with_most, shape=lc_flag), size=2) +
  geom_text(aes(x=longestHeavy_maxSubLength, y=domCloneH_maxSubLength, label=sample), 
            hjust=-0.1, size=1) +
  geom_abline(intercept=0, slope=1) +
  ggtitle("Longest Heavy Chain Subsequence Length\nvs\nDominant Clone Longest Subsequence Length") +
  xlab("Longest Heavy Chain Subsequence Length (bp)") +
  ylab("Dominant IGH Clone Longest Subsequence Length")
# show(plotH_maxsublength_alllabels)
ggsave(paste(out_dir, "plotH_maxsublength_alllabels.svg", sep=""), plotH_maxsublength_alllabels,  width=8, height=6, units="in", dpi=600)

# Do against total IG for heavy:
plotH_domlength <- comb_summary %>% 
  ggplot(aes(x=total, y=domCloneH_maxSubLength)) +
  geom_point(aes(color=lc_with_most, shape=lc_flag), size=2) +
  geom_smooth(method="gam", se=FALSE, color="black") + # , label="GAM"
  geom_smooth(method="lm", se=FALSE, color="green") + # , label="lm"
  ggtitle("Total Ig vs\nDominant Heavy Clone Longest Subsequence Length") +
  xlab("Total Ig") +
  ylab("Dominant Heavy Clone Longest Subsequence Length")
# show(plotH_domlength)
ggsave(paste(out_dir, "plotH_domCloneSublength.svg", sep=""), plotH_domlength,  width=8, height=6, units="in")
# Max heavy length, even if not dominant:
plotH_maxsublength2 <- comb_summary %>% 
  ggplot(aes(x=total, y=longestHeavy_maxSubLength)) +
  geom_point(aes(color=lc_with_most, shape=lc_flag), size=2) +
  geom_smooth(method="gam", se=FALSE, color="black") + # , label="GAM"
  geom_smooth(method="lm", se=FALSE, color="green") + # , label="lm"
  ggtitle("Total Ig vs\nLongest Heavy Subsequence Length") +
  xlab("Total Ig") +
  ylab("Longest Heavy Subsequence Length")
# show(plotH_maxsublength2)
ggsave(paste(out_dir, "plotH_maxsublength2.svg", sep=""), plotH_maxsublength2,  width=8, height=6, units="in")


# Graph total light and total heavy
plot_total_counts <- comb_summary %>% 
  ggplot(aes(x=IGH, y=light)) +
  geom_point(aes(color=mainclone_lowfreq, shape=mainclone_notlongest)) +
  geom_text(data=subset(comb_summary, sample %in% samples_watch), 
            aes(x=IGH, y=light, label=sample), 
            hjust=-0.1, size=1) +
  scale_shape_manual(values=c(16, 3)) +
  scale_color_manual(values=c("blue", "red")) +
  ggtitle("IGH vs Light Count") +
  xlab("Number of IGH sequences") +
  ylab("Number of light chain sequences") +
  scale_x_continuous(trans="log2", oob=scales::squish_infinite) +
  scale_y_continuous(trans="log2")
# show(plot_total_counts)
ggsave(paste(out_dir, "plot_total_counts.svg", sep=""), plot_total_counts,  width=8, height=6, units="in", dpi=600)


# Number of IGK, IGL, IGH, light chaing vs IG count
comb_summary_long <- comb_summary %>% 
  pivot_longer(cols=c(IGK, IGL, IGH, light), names_to="graph_var", values_to="graph_value")
comb_summary_long$graph_var <- factor(comb_summary_long$graph_var, levels=c("IGK", "IGL", "light", "IGH"))


# Plot chain count vs IG count:
plot_total_ig <- comb_summary_long %>% 
  ggplot(aes(x=total, y=graph_value)) +
  geom_point(aes(color=domClone_Chain, shape=lc_flag), size=2) +
  geom_smooth(method="gam", se=FALSE, aes(color=domClone_Chain)) +
  ggtitle("Count vs IG Count") +
  xlab("IG Count") +
  ylab("Chain Count") +
  facet_wrap(vars(graph_var))
# show(plot_total_ig)
ggsave(paste(out_dir, "plot_total_ig.svg", sep=""), plot_total_ig,  width=12, height=8, units="in")


# Combined scattered plots:
combined_scattered <- ggarrange(plot_lightmaxsublength, plotH_maxsublength, nrow=1)
# show(combined_scattered)
ggsave(paste(out_dir, "combined_lengthvslength.svg", sep=""), combined_scattered,  width=14, height=6, units="in", dpi=600)
# Combined LC length plots:
combined_lc <- ggarrange(plot_domlength, plot_lightmaxsublength2, nrow=1)
# show(combined_lc)
ggsave(paste(out_dir, "combined_lc.svg", sep=""), combined_lc,  width=14, height=6, units="in", dpi=600)
# Combined heavy length plots:
combined_H <- ggarrange(plotH_domlength, plotH_maxsublength2, nrow=1)
# show(combined_H)
ggsave(paste(out_dir, "combined_h.svg", sep=""), combined_H,  width=14, height=6, units="in", dpi=600)


# #############################################################################
# Violin plot for dominant clone subsequence length, maximum length, for light and heavy:
# Make longer version of dataframe:
comb_summary_long_length <- comb_summary %>% 
  pivot_longer(cols=c(domClone_maxSubLength, longest_maxSubLength, 
                      domCloneH_maxSubLength, longestHeavy_maxSubLength), 
               names_to="graph_var", values_to="graph_value")
comb_summary_long_length$graph_var <- factor(comb_summary_long_length$graph_var, 
                                             levels=c("domClone_maxSubLength", 
                                                      "longest_maxSubLength", 
                                                      "domCloneH_maxSubLength", 
                                                      "longestHeavy_maxSubLength"))
# Plot Violin of comb_summary_long_length:
violin_length <- comb_summary_long_length %>% 
  ggplot(aes(x=graph_var, y=graph_value)) +
  geom_violin(width=0.8) +
  geom_boxplot(width=0.1, outlier.shape=NA) +
  geom_jitter(inherit.aes=FALSE, shape=16, position=position_jitter(0.2), 
              size=0.5, aes(x=graph_var, y=graph_value, color=graph_var)) +
  ggtitle("MiXCR Sequence Length") +
  ylab("Sequence Length") +
  theme(legend.position = "none")
# show(violin_length)
ggsave(paste(out_dir, "violin_length.svg", sep=""), violin_length,  width=12, height=8, units="in")

# Number of clones per sample vs Fraction of dominant clones
cloneN_domfreq <- comb_summary  %>% 
  ggplot(aes(x=number_of_clones, y=`dom_clone/light`)) +
  geom_point(aes(color=lc_with_most), size=1) +
  ggtitle("Number of LC Clones vs\nDominant Light Chain Clone Frequency") +
  xlab("Number of Distinct LC Clones") +
  ylab("Dominant Light Chain Clone Frequency")
# show(cloneN_domfreq)
ggsave(paste(out_dir, "cloneN_domfreq.svg", sep=""), cloneN_domfreq,  width=8, height=6, units="in")
# Do same in log scale:
cloneN_domfreqlog <- cloneN_domfreq +
  scale_x_continuous(trans="log10") +
  xlab("Number of Distinct LC Clones (log10)")
# show(cloneN_domfreqlog)
ggsave(paste(out_dir, "cloneN_domfreqlog.svg", sep=""), cloneN_domfreqlog,  width=8, height=6, units="in")

# LC count vs number of clones
lcN_domfreq <- comb_summary  %>% 
  ggplot(aes(y=number_of_clones, x=light)) +
  geom_point(aes(color=lc_with_most), size=1) +
  ggtitle("Number of LC Clones\nvs LC Count ") +
  ylab("Number of Distinct LC Clones") +
  xlab("LC Count")
# show(lcN_domfreq)
ggsave(paste(out_dir, "lcN_domfreq.svg", sep=""), lcN_domfreq,  width=8, height=6, units="in")
# Do same in log scale:
lcN_domfreqlog <- lcN_domfreq +
  scale_y_continuous(trans="log10") +
  ylab("Number of Distinct LC Clones (log10)")
# show(lcN_domfreqlog)
ggsave(paste(out_dir, "lcN_domfreqlog.svg", sep=""), lcN_domfreqlog,  width=8, height=6, units="in")


# LC count vs dom clone frequency
lcN_domcloneN <- comb_summary  %>% 
  ggplot(aes(x=light, y=`lc_dom/light`)) +
  geom_point(aes(color=lc_with_most), size=1) +
  ggtitle("LC Count vs\nDominant Light Chain Clone Frequency") +
  xlab("LC Count") +
  ylab("Dominant LC Clone Frequency")
# show(lcN_domcloneN)
ggsave(paste(out_dir, "lcN_domcloneN.svg", sep=""), lcN_domcloneN,  width=8, height=6, units="in")

# IG count vs dom clone frequency
igN_domcloneN <- comb_summary  %>% 
  ggplot(aes(x=total, y=`lc_dom/light`)) +
  geom_point(aes(color=lc_with_most), size=1) +
  ggtitle("Total Ig vs\nDominant Light Chain Clone Frequency") +
  xlab("Total Ig") +
  ylab("Dominant LC Clone Frequency")
# show(igN_domcloneN)
ggsave(paste(out_dir, "igN_domcloneN.svg", sep=""), igN_domcloneN,  width=8, height=6, units="in")