# Create FASTA's for germline

# Get command line arguments:
# (CSV with MIXCR ratio) (fasta file location) (directory to output results to) (IMGT reference file, csv) (IMGT reference leader sequences, csv)
args = commandArgs(trailingOnly=TRUE)
mixcr_ratio = args[1]
fasta_located_at = args[2]
outdir = args[3]
imgt_file = args[4]
imgt_leader_file = args[5]


# If not using Command Line:
#  mixcr_ratio = "sra/800MiXCR/SUMMARY/mixcr_ratio.csv"
#  fasta_located_at = "sra/800MiXCR/SUMMARY/fasta/"
#  outdir = "sra/800MiXCR/SUMMARY/fasta_germline/"
#  # IMGT reference:
#  imgt_file = "IMGTref.csv"
#  # IMGT leader sequences:
#  imgt_leader_file = "IMGTleadersequences.csv"


# Load libraries:
library(tidyverse)
library("Biostrings")  # Not sure which package I want to use
library("msa")
library("stringr")
# library("DECIPHER")

# Read in IMGT dataframe:
imgt_df <- as_tibble(read_csv(imgt_file))
imgt_leader_df <- as_tibble(read_csv(imgt_leader_file))
# IMGT leader sequences:
# >M64856|IGKV1-33*01|Homo sapiens|F|L-PART1+L-PART2|131..185+310..320|66 nt|1| | | | |66+0=66| | |
igk_leader <- "atggacatgagggtccctgctcagctcctggggctcctgcagctctggctctcaggtgccagatgt" 
# >X71966|IGLV3-21*01|Homo sapiens|F|L-PART1+L-PART2|194..239+677..687|57 nt|1| | | | |57+0=57| | |
igl_leader <- "atggcctggaccgttctcctcctcggcctcctctctcactgcacaggctctgtgacc"

# Keep only functional genes:
imgt_df <- imgt_df %>% 
  filter(Functional=="F")
imgt_leader_df <- imgt_leader_df %>% 
  filter(Functional=="F")

print("Saving fasta to:")
print(outdir)
# Remove trailing slash if present of output directory:
if (stringr::str_sub(outdir, -1, -1) == "/"){
  outdir_name <- str_sub(outdir, 1, -2)
} else {
  outdir_name <- outdir
}

# Set working directory:
setwd(outdir_name)
print("Working directory:")
print(getwd())

# Get samples:
mixcr_df <- as_tibble(read_csv(mixcr_ratio))
samples <- mixcr_df$sample

# Get main light chain ID:
samples_light <- mixcr_df %>% 
  filter(sample %in% samples) %>% 
  select(sample, lc_with_most)


# Function for fasta reference filenames: --------------------------------------
get_fasta_fn <- function(sample_name, chain){
  filename <- paste(fasta_located_at, sample_name, "_", chain, ".fasta", sep="")
  return(filename)
}

# Function for fasta output filenames: -----------------------------------------
get_fasta_out <- function(sample_name, chain){
  filename <- paste(outdir, sample_name, "_", chain, "_germline.fasta", sep="")
  return(filename)
}

# Function for looking up germline: -------------------------------------------
# (IMGT reference dataframe, IG chain, gene to lookup)
get_germline_seq <- function(imgt_table1, chain, gene){
  # Clean up gene name:
  pieces <- str_split(gene, "-")[[1]]
  if (length(pieces) > 1){
    gene <- paste(pieces[1], toString(as.integer(pieces[2])), sep="-")
    gene_org <- gene
  } else {
    gene_org <- gene
  }
  # Look up seq:
  seq <- imgt_table1$SeqNoGaps[imgt_table1$Locus==chain & imgt_table1$AlleleGroupNoChain==gene &
                             imgt_table1$SpecificAllele=="01"]
  # If lookup didn't work, do something different for constant region:
  if ((length(seq)<1) & (substr(gene, 4, 4) == "C")){
    seq <- imgt_table1$SeqNoGaps[imgt_table1$Locus==chain & imgt_table1$AlleleGroup==gene &
                               imgt_table1$SpecificAllele=="01"]
  }
  # Try with adding D to the name:
  if((length(seq)<1) & (length(pieces) > 1)){
    gene <- paste(pieces[1], "D-", toString(as.integer(pieces[2])), sep="")
    # Look up seq:
    seq <- imgt_table1$SeqNoGaps[imgt_table1$Locus==chain & imgt_table1$AlleleGroupNoChain==gene &
                               imgt_table1$SpecificAllele=="01"]
  }
  # Try looking for *02 instead:
  if (length(seq)<1){
    gene <- gene_org
    # Look up seq:
    seq <- imgt_table1$SeqNoGaps[imgt_table1$Locus==chain & imgt_table1$AlleleGroupNoChain==gene &
                                   imgt_table1$SpecificAllele=="02"]
  }
  # Try with adding D to the name:
  if((length(seq)<1) & (length(pieces) > 1)){
    gene <- paste(pieces[1], "D-", toString(as.integer(pieces[2])), sep="")
    # Look up seq:
    seq <- imgt_table1$SeqNoGaps[imgt_table1$Locus==chain & imgt_table1$AlleleGroupNoChain==gene &
                                   imgt_table1$SpecificAllele=="02"]
  }
  # Try looking for *03 instead:
  if (length(seq)<1){
    gene <- gene_org
    # Look up seq:
    seq <- imgt_table1$SeqNoGaps[imgt_table1$Locus==chain & imgt_table1$AlleleGroupNoChain==gene &
                                   imgt_table1$SpecificAllele=="03"]
  }
  # Try with adding D to the name:
  if((length(seq)<1) & (length(pieces) > 1)){
    gene <- paste(pieces[1], "D-", toString(as.integer(pieces[2])), sep="")
    # Look up seq:
    seq <- imgt_table1$SeqNoGaps[imgt_table1$Locus==chain & imgt_table1$AlleleGroupNoChain==gene &
                                   imgt_table1$SpecificAllele=="03"]
  }
  return(seq)
}


# Initiate an empty tibble to store germline information:
germline_df <- tibble(
  sample = character(),
  chain = character(),
  germline_reference_name = character(),
  germline_reference = character(),
  germline_reference_translated = character()
)


# ---------------------------------------------------------------------------
# Add fasta filenames for column that 
samples_light <- samples_light %>% 
  mutate(fasta = get_fasta_fn(sample, lc_with_most),
         fasta_out = get_fasta_out(sample, lc_with_most))

# Produce FASTA's:
for (i in c(1:length(samples_light$sample))){
  sample_name <- samples_light$sample[i]
  chain <- samples_light$lc_with_most[i]
  print(sprintf("Working on sample %s", sample_name))
  # Read in file:
  fasta_file <- samples_light$fasta[i]
  # Fasta out:
  fasta_out_fn <- samples_light$fasta_out[i]
  # Check if file exists:
  if (file.exists(fasta_file)){
    the_sequences <- readDNAStringSet(fasta_file)
      # Get top 2 sequences from DNA string set:  
      stringset_df <- tibble(names=names(the_sequences), seq=as.character(the_sequences))
      # Get sample name:
      stringset_df <- stringset_df %>% 
        separate(col=names, into=c("sample", "Chain", "cloneID", "cloneRank", "cloneCount", "V", "J", "C", "subseq"), sep="_")
      # Clean up column entries:
      stringset_df <- stringset_df %>% 
        mutate_at("cloneID", str_replace, "cID-", "") %>% 
        mutate_at("cloneCount", str_replace, "cCount-", "") %>% 
        mutate_at("cloneRank", str_replace, "-", "") %>% 
        mutate_at("V", str_replace, "V-", "") %>% 
        mutate_at("J", str_replace, "J-", "") %>% 
        mutate_at("C", str_replace, "C-", "")
      # Look up VJ information: 
      germline_V <- stringset_df$V[stringset_df$cloneRank=="cRank0" & (stringset_df$subseq=="subseq1" | stringset_df$subseq=="allTogether")][1]
      germline_V_seq <- get_germline_seq(imgt_table1=imgt_df, chain=chain, gene=germline_V)
      germline_leader_seq <- get_germline_seq(imgt_table1=imgt_leader_df, chain=chain, gene=germline_V)
      # metadata(germline_V_seq)$name <- paste("germline", germline_V, sep="_")
      # print(sprintf("Germline V: %s", germline_V))
      # print(germline_V_seq)
      # print("Leader sequence:")
      # print(germline_leader_seq)
      # Check for missing V gene:
      if (length(germline_V_seq) == 0){
        print(sprintf("Missing V germline gene %s %s", chain, germline_V))
      } else if (is.na(germline_V_seq)){
        print(sprintf("Missing V germline gene %s %s", chain, germline_V))
      } else if (nchar(germline_V_seq, allowNA=TRUE) < 10){
        print(sprintf("Missing V germline gene %s %s", chain, germline_V))
      }
      # Check for missing leader gene:
      if (length(germline_leader_seq) == 0){
        print(sprintf("Missing Leader gene %s %s", chain, germline_V))
      } else if (is.na(germline_leader_seq)){
        print(sprintf("Missing Leader gene %s %s", chain, germline_V))
      } else if (nchar(germline_leader_seq, allowNA=TRUE) < 10){
        print(sprintf("Missing Leader gene %s %s", chain, germline_V))
      }
      #
      germline_J <- stringset_df$J[stringset_df$cloneRank=="cRank0" & (stringset_df$subseq=="subseq1" | stringset_df$subseq=="allTogether")][1]
      germline_J_seq <- get_germline_seq(imgt_table1=imgt_df, chain=chain, gene=germline_J)
      germline_C <- stringset_df$C[stringset_df$cloneRank=="cRank0" & (stringset_df$subseq=="subseq1" | stringset_df$subseq=="allTogether")][1]
      germline_C_seq <- get_germline_seq(imgt_table1=imgt_df, chain=chain, gene=germline_C)
      # metadata(germline_J_seq)$name <- paste("germline", germline_J, sep="_")
      germline_seq_to_include <- c()
      names_to_include <- c()
      combined_seq <- ""
      combined_seq_name <- ""
      # If V gene is available (only include leader if V is available):
      if (!identical(germline_V_seq, character(0))){
        # If leader sequence is available:
        if (!identical(germline_leader_seq, character(0))){
          temp_name <- "germLead"
        } else {
          # If germline squence wasn't available
          if (chain == "IGK"){
            germline_leader_seq <- igk_leader
            temp_name <- "substituteGermLead"
          } else if (chain == "IGL"){
            germline_leader_seq <- igl_leader
            temp_name <- "substituteGermLead"
          }
        }
        germline_seq_to_include <- append(germline_seq_to_include, germline_leader_seq)
        # Get first 3 nucleotides of leader:
        leader_start <- substr(germline_leader_seq, 1, 3)
        leader_start <- toupper(leader_start)
        # Pad string if string not divisible by 3:
        while (nchar(germline_leader_seq) %% 3 != 0){
          # print(sprintf("In while loop 1, length: %d", nchar(germline_leader_seq)))
          if (leader_start == "ATG"){
            germline_leader_seq <- paste(germline_leader_seq, "n", sep="")
          } else {
            germline_leader_seq <- paste("n", germline_leader_seq, sep="")
          }
        }
        # Pad string if < 51:
        while (nchar(germline_leader_seq) < 51){
          # print(sprintf("In while loop 2, length: %d", nchar(germline_leader_seq)))
          if (leader_start == "ATG"){
            germline_leader_seq <- paste(germline_leader_seq, "n", sep="")
          } else {
            germline_leader_seq <- paste("n", germline_leader_seq, sep="")
          }
        }
        # Add leader sequence:
        names_to_include <- append(names_to_include, temp_name)
        combined_seq <- paste(combined_seq, germline_leader_seq, sep="")
        combined_seq_name <- paste(combined_seq_name, temp_name, sep="")
        # Add V gene sequence:
        germline_seq_to_include <- append(germline_seq_to_include, germline_V_seq)
        temp_name <- paste("germ_", chain, germline_V, sep="")
        names_to_include <- append(names_to_include, temp_name)
        combined_seq <- paste(combined_seq, germline_V_seq, sep="")
        combined_seq_name <- paste(combined_seq_name, temp_name, sep="__")
      }
      # If J gene is available:
      if (!identical(germline_J_seq, character(0))){
        germline_seq_to_include <- append(germline_seq_to_include, germline_J_seq)
        temp_name <- paste("germ_", chain, germline_J, sep="")
        names_to_include <- append(names_to_include, temp_name)
        combined_seq <- paste(combined_seq, germline_J_seq, sep="")
        combined_seq_name <- paste(combined_seq_name, temp_name, sep="__")
        # Pad string if string not divisible by 3:
        while (nchar(combined_seq) %% 3 != 0){
          # print(sprintf("In while loop 2, length: %d", nchar(combined_seq))) 
          combined_seq <- paste(combined_seq, "n", sep="")
        }
        # If constant gene is available:
        if (!identical(germline_C_seq, character(0))){
          germline_seq_to_include <- append(germline_seq_to_include, germline_C_seq)
          temp_name <- paste("germ", germline_C, sep="")
          names_to_include <- append(names_to_include, temp_name)
          combined_seq <- paste(combined_seq, germline_C_seq, sep="")
          combined_seq_name <- paste(combined_seq_name, temp_name, sep="__")
        }
      }
      # Get translated sequence:
      temp_sequence <- DNAString(combined_seq)
      temp_translated <- translate(temp_sequence, if.fuzzy.codon="solve")
      temp_translated <- toString(temp_translated)
      # Add row to tibble:
      germline_df <- germline_df %>% 
        add_row(sample = sample_name,
                chain = chain,
                germline_reference_name = combined_seq_name,
                germline_reference = combined_seq,
                germline_reference_translated = temp_translated)
      # Add in the combined scenario:
      germline_seq_to_include <- append(germline_seq_to_include, combined_seq)
      names_to_include <- append(names_to_include, combined_seq_name)
      # Create the Stringset:
      germline_seq <- DNAStringSet(germline_seq_to_include)
      names(germline_seq) <- names_to_include
      # Write germline DNAStringSet to fasta:
      writeXStringSet(germline_seq, filepath=fasta_out_fn, append=FALSE, 
                      compress=FALSE, format="fasta")
  }
}


# Save dataframe with germline:
write.csv(germline_df, file=paste(outdir, "reference_germline.csv", sep=""), row.names=FALSE)


# #############################################################################
# Lengths of leader sequences:
imgt_leader_df$length <- nchar(imgt_leader_df$SeqNoGaps)
imgt_leader_df$length_f <- factor(imgt_leader_df$length)

imgt_leader_df %>% 
  ggplot(aes(x=length_f, group=Locus, fill=Locus)) +
  geom_bar(position="dodge") +
  ggtitle("Length of Leader Sequence")

table(imgt_leader_df$length_f, imgt_leader_df$Locus)
