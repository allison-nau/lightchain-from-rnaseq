# Look at MiXCR timings & MultiQC/FastQC related statistics


# Import libraries:
library(tidyverse)
library(ggplot2)
library(GGally)
library(janitor)
# library(lubridate)
# library(tibbletime)
# library(chron)

# Get command line arguments:
# (CSV with MiXCR ratio and stats made using R) (location of standardout 
# produced using grep and processed using python, 
# leaving off column with original filenames) (directory to output files to)
args = commandArgs(trailingOnly=TRUE)
read_in = args[1]
accounting_in = args[2]
stdout_in = args[3]
out_dir = args[4]
# For running code without command args:
# read_in = "sra/800MiXCR/SUMMARY/mixcr_ratio.csv"
# accounting_in = "sra/800MiXCR/SUMMARY/job_account.csv"
# stdout_in = "sra/800MiXCR/SUMMARY/processed_stdout_qc_wide_NOTrestricted.csv"
# out_dir = "sra/800MiXCR/SUMMARY/"

# Scientific notation options:
options(scipen=10000)


# TODO set working directory

# ##############################################################################
# Read in timings:
timings <- as_tibble(read_csv(accounting_in, 
                              col_types=cols(.default="?",
                                             taskid="i", 
                                             jobnumber="i", 
                                             ru_wallclock="d",
                                             ru_utime="d",
                                             ru_stime="d",
                                             cpu="d",
                                             maxvmem="c"),
                              na = c("", "NA", "NaN", "undefined")))
# Create job_task column:
timings$job_task <- paste(timings$jobname, ".o", timings$jobnumber, ".", 
                           timings$taskid, sep="")

# Remove duplicate rows:
timings <- timings %>% distinct(jobname, jobnumber, taskid, qsub_time, 
                                .keep_all=TRUE)

# ONLY KEEP pipeline rows
timings <- timings %>% 
  filter(jobname == "pipeline.qsub")
timings <- timings %>% 
  group_by(taskid) %>% 
  top_n(1, jobnumber) %>% 
  ungroup()


#   # Convert columns to time:
#   timings <- timings %>%
#     rowwise() %>% 
#     mutate(User_min = as.numeric(as.difftime(User.Time, format="%H:%M:%S", units="mins")),
#            System_min = as.numeric(as.difftime(System.Time, format="%H:%M:%S", units="mins")),
#            Wallclock_min = as.numeric(as.difftime(Wallclock.Time, format="%H:%M:%S", units="mins")),
#            CPU_min = as.numeric(as.difftime(CPU, format="%H:%M:%S", units="mins")))
# Convert VMEM to number:
timings <- timings %>% 
  rowwise() %>% 
  mutate(maxvmem = as.numeric(gsub("G", "", maxvmem)))
#   # Separate sample index:
#   timings <- timings %>% 
#     separate(file, c("Job", "Index"), sep="\\.")


# Convert timings to MINUTES (/60^2)
time_cols <- c("ru_wallclock", "ru_utime", "ru_stime", "cpu")
new_names <- c("Wallclock_min", "User_min", "System_min", "CPU_min")
for (i in 1:length(time_cols)){
  timings[,new_names[i]] <- timings[,time_cols[i]] / (60)
}

# Read in mixer summary results:
summary <- as_tibble(read.csv(read_in))

# Read in standard out:
stdoutdf <- as_tibble(read_csv(stdout_in, na=c("")))
# Change column name:
stdoutdf <- stdoutdf %>% 
  rename(sample = `Before MiXCR, file name will be changed to`,
         toomanypartial = `WARNING: too many partial alignments detected, consider skipping assemblePartial (enriched library?).`)
# Fill NA in too many partials with "NA"
stdoutdf <- stdoutdf %>%
  replace_na(list(toomanypartial =  "NA"))
# Convert to factor:
stdoutdf$toomanypartial <- factor(stdoutdf$toomanypartial)


# Combine MiXCR summary and standard out:
summary2 <- merge(summary, stdoutdf, by="sample", all=TRUE)
# Check which rows are missing MiXCR values:
print("Samples missing MiXCR values:")
summary2$sample[is.na(summary2$total)]
summary2[is.na(summary2$total), c("sample", "job_task")]
# Add in timings:
timings <- merge(timings, summary2, by="job_task", all=TRUE)

# Drop duplicated rows, keeping latest version:
timings <- timings %>% 
  group_by(taskid) %>% 
  top_n(1, jobnumber) %>% 
  ungroup()


timings$samplegroup <- "middle"
timings <- timings %>% relocate(samplegroup)  # TODO remove
timings$samplegroup[timings$taskid<=50] <- "smallest"
timings$samplegroup[timings$taskid>700] <- "largest"
timings$samplegroup <- factor(timings$samplegroup, 
                              levels=c("largest", "middle", "smallest"))

# Save:
write.csv(timings, file=paste(out_dir, "mixcr_ratio_timings_stdout.csv", sep=""), row.names=FALSE)


# #####################################################################
# Boxplot timing & vmem:

# Columns to keep
keep_time <- c("taskid", "samplegroup", "Wallclock_min", "User_min", "System_min", "CPU_min", "maxvmem")
keep_time2 <- c("Wallclock_min", "User_min", "System_min", "CPU_min", "maxvmem")
# Convert timings to long:
timings_long <- timings %>% 
  select(all_of(keep_time)) %>% 
  pivot_longer(cols=keep_time2)
timings_long$name = factor(timings_long$name, levels=keep_time2)

# Make violin plot:
timing_violin <- timings_long %>% 
  ggplot(aes(x=name, y=value)) +
  geom_violin(width=0.8) +
  geom_boxplot(width=0.1, outlier.shape=NA) +
  geom_jitter(inherit.aes=FALSE, shape=16, position=position_jitter(0.2), size=0.5, aes(x=name, y=value, color=samplegroup)) +
  ggtitle("Timings (min) + vmem (GB)") +
  xlab("") +
  ylab("") +
  facet_wrap(~name, scales="free")
# show(timing_violin)
ggsave(paste(out_dir, "timing_violin.png", sep=""), timing_violin, device=png, width=12, height=8, units="in")


# #############################################################################
# MultiQC/FastQC 

# Make Mean quality score a factor:
timings$fastqc_per_base_sequence_quality <- factor(timings$fastqc_per_base_sequence_quality, 
                                                   levels=c("pass", "warn", "fail"))
# Make longest sequence length a factor:
timings$fastqc_length_longest_f <- "other"
timings$fastqc_length_longest_f[timings$fastqc_1_length_longest == 82 | timings$fastqc_1_length_longest == 83] <- "82-83"
timings$fastqc_length_longest_f[timings$fastqc_1_length_longest == 87] <- "87"
timings$fastqc_length_longest_f[timings$fastqc_1_length_longest == 100 | timings$fastqc_1_length_longest == 101] <- "100-101"
timings$fastqc_length_longest_f[timings$fastqc_1_length_longest == 109] <- "109"
timings$fastqc_length_longest_f <- factor(timings$fastqc_length_longest_f, levels=c("82-83", "87", "100-101", "109"))

# MultiQC Mean Quality Scores vs Light Chain Length
mean_quality <- timings %>% 
  ggplot(aes(x=fastqc_per_base_sequence_quality, y=domClone_maxSubLength)) +
  geom_violin(width=0.8) +
  geom_boxplot(width=0.1, outlier.shape=NA) +
  geom_jitter(inherit.aes=FALSE, shape=16, position=position_jitter(0.2), 
              size=0.5, aes(x=fastqc_per_base_sequence_quality, 
                            y=domClone_maxSubLength, color=samplegroup)) +
  ggtitle("Mean Quality Scores and Dominant Clone Maximum Sublength") +
  ylab("Dominant Clone Max Sublength")
# show(mean_quality)
ggsave(paste(out_dir, "qualityandlength.png", sep=""), mean_quality, device=png, width=8, height=8, units="in")

# MultiQC Longest Sequence length VS Light Chain Length
sequencing_length <- timings %>% 
  ggplot(aes(x=fastqc_length_longest_f, y=domClone_maxSubLength)) +
  geom_violin(width=0.8) +
  geom_boxplot(width=0.1, outlier.shape=NA) +
  geom_jitter(inherit.aes=FALSE, shape=16, position=position_jitter(0.2), 
              size=0.5, aes(x=fastqc_length_longest_f, 
                            y=domClone_maxSubLength, color=fastqc_per_base_sequence_quality)) +
  scale_color_manual(values = c("pass" = "blue",
                                "warn" =  "orange",
                                "fail" = "red")) +
  ggtitle("Sequencing Length and Dominant Clone Maximum Sublength") +
  xlab("Sequencing Length") +
  ylab("Dominant Clone Max Sublength")
# show(sequencing_length)
ggsave(paste(out_dir, "sequencinglengthandmixcrlength.png", sep=""), sequencing_length, device=png, width=12, height=8, units="in")





# #####################################################################







