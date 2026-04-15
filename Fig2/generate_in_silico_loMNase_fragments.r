# Load required packages
if (!require("dplyr")) install.packages("dplyr")
if (!require("data.table")) install.packages("data.table")
library(dplyr)
library(data.table)

# Function: Read CTCF BED file and fragment length frequency file
read_input_files <- function(ctcf_file, freq_file) {
  # Read CTCF BED file
  ctcf_data <- fread(ctcf_file, sep = "\t", header = FALSE,
                     col.names = c("chrom", "start", "end", "name", "score", "strand"))
  
  # Calculate length of each region
  ctcf_data$region_length <- ctcf_data$end - ctcf_data$start
  
  # Read fragment length frequency file
  freq_data <- fread(freq_file, sep = "\t", header = FALSE,
                     col.names = c("length", "frequency"))
  
  # Normalize frequencies to sum to 1
  freq_data$frequency <- freq_data$frequency / sum(freq_data$frequency)
  
  return(list(ctcf_data = ctcf_data, freq_data = freq_data))
}

# Function: Generate random fragments based on frequency distribution
generate_random_fragments <- function(ctcf_data, freq_data, fragments_per_region = 1, seed = 123) {
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Check if fragment lengths are feasible within regions
  max_region_length <- max(ctcf_data$region_length)
  available_lengths <- freq_data$length
  min_required_length <- min(available_lengths)
  
  # Check for regions shorter than the minimum available fragment length
  short_regions <- ctcf_data$region_length[ctcf_data$region_length < min_required_length]
  if (length(short_regions) > 0) {
    warning(paste(length(short_regions), 
                  "regions are shorter than the minimum available fragment length (", 
                  min_required_length, "bp) and cannot generate fragments"))
  }
  
  # Generate random fragments for each region
  all_fragments <- list()
  fragment_counter <- 1
  
  for (i in 1:nrow(ctcf_data)) {
    chrom <- ctcf_data$chrom[i]
    region_start <- ctcf_data$start[i]
    region_end <- ctcf_data$end[i]
    region_length <- ctcf_data$region_length[i]
    strand <- ctcf_data$strand[i]
    original_name <- ctcf_data$name[i]
    
    # Generate specified number of fragments for this region
    for (j in 1:fragments_per_region) {
      # Sample fragment length from frequency distribution
      # Only consider fragment lengths <= region length
      valid_lengths <- freq_data$length[freq_data$length <= region_length]
      valid_freqs <- freq_data$frequency[freq_data$length <= region_length]
      
      # If valid lengths exist
      if (length(valid_lengths) > 0 && sum(valid_freqs) > 0) {
        # Re-normalize probabilities
        valid_freqs <- valid_freqs / sum(valid_freqs)
        
        # Randomly select a length from valid lengths
        fragment_length <- sample(valid_lengths, 1, prob = valid_freqs)
        
        # Calculate available start position range
        max_start <- region_end - fragment_length
        
        # Randomly select start position if region length is sufficient
        if (max_start >= region_start) {
          fragment_start <- sample(region_start:max_start, 1)
        } else {
          fragment_start <- region_start
        }
        
        fragment_end <- fragment_start + fragment_length
        
        # Create fragment record
        fragment <- data.frame(
          chrom = chrom,
          start = fragment_start,
          end = fragment_end,
          name = paste0(original_name, "_frag", j),
          score = fragment_length,  # Use fragment length as score
          strand = strand,
          original_name = original_name,
          fragment_length = fragment_length,
          stringsAsFactors = FALSE
        )
        
        all_fragments[[fragment_counter]] <- fragment
        fragment_counter <- fragment_counter + 1
      } else {
        warning(paste("Region", i, "(", chrom, ":", region_start, "-", region_end, 
                      ") cannot generate fragments because region length (", region_length, 
                      "bp) is shorter than all available fragment lengths"))
      }
    }
  }
  
  # Combine all fragments
  if (length(all_fragments) > 0) {
    fragments_df <- bind_rows(all_fragments)
    
    # Reorder columns
    fragments_df <- fragments_df[, c("chrom", "start", "end", "name", "score", "strand", 
                                     "original_name", "fragment_length")]
    
    return(fragments_df)
  } else {
    warning("No fragments were generated")
    return(NULL)
  }
}

# Function: Save results to files
save_results <- function(fragments_df, output_prefix = "random_fragments") {
  if (is.null(fragments_df)) {
    warning("No fragments to save")
    return(NULL)
  }
  
  # Save BED file
  bed_file <- paste0(output_prefix, ".bed")
  fwrite(fragments_df[, 1:6], bed_file, sep = "\t", 
         row.names = FALSE, col.names = FALSE, quote = FALSE)
  message(paste("BED file saved to:", bed_file))
  
  # Save complete results (with additional information)
  full_file <- paste0(output_prefix, "_full.tsv")
  fwrite(fragments_df, full_file, sep = "\t", 
         row.names = FALSE, col.names = TRUE, quote = FALSE)
  message(paste("Complete results saved to:", full_file))
  
  return(list(bed_file = bed_file, full_file = full_file))
}

# Function: Validate generated fragments
validate_fragments <- function(fragments_df, ctcf_data, freq_data) {
  if (is.null(fragments_df)) {
    return(NULL)
  }
  
  message("=== Fragment Validation Results ===")
  
  # 1. Verify fragments are within original regions
  validation_results <- data.frame(
    fragment = fragments_df$name,
    in_original_region = FALSE,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(fragments_df)) {
    chrom <- fragments_df$chrom[i]
    start <- fragments_df$start[i]
    end <- fragments_df$end[i]
    original_name <- fragments_df$original_name[i]
    
    # Find corresponding original region
    original_row <- ctcf_data[ctcf_data$name == original_name, ]
    
    if (nrow(original_row) == 1) {
      original_start <- original_row$start
      original_end <- original_row$end
      
      validation_results$in_original_region[i] <- (
        start >= original_start & 
          end <= original_end
      )
    }
  }
  
  in_region_percent <- mean(validation_results$in_original_region) * 100
  message(paste("Percentage of fragments within original regions:", round(in_region_percent, 2), "%"))
  
  # 2. Verify fragment length distribution
  fragment_length_counts <- table(fragments_df$fragment_length)
  fragment_length_freq <- fragment_length_counts / sum(fragment_length_counts)
  
  # Create length distribution comparison
  length_comparison <- data.frame(
    length = as.numeric(names(fragment_length_freq)),
    generated_frequency = as.numeric(fragment_length_freq)
  )
  
  # Merge with target frequencies
  length_comparison <- merge(length_comparison, 
                             freq_data[, c("length", "frequency")], 
                             by = "length", all.x = TRUE)
  colnames(length_comparison) <- c("length", "generated_frequency", "target_frequency")
  
  # Calculate differences
  length_comparison$difference <- length_comparison$generated_frequency - 
    length_comparison$target_frequency
  
  message("\nFragment Length Distribution Comparison:")
  print(length_comparison)
  
  # 3. Basic statistics
  message("\nBasic Statistics:")
  message(paste("Total fragments:", nrow(fragments_df)))
  message(paste("Average fragment length:", round(mean(fragments_df$fragment_length), 2), "bp"))
  message(paste("Minimum fragment length:", min(fragments_df$fragment_length), "bp"))
  message(paste("Maximum fragment length:", max(fragments_df$fragment_length), "bp"))
  message(paste("Number of distinct lengths used:", length(unique(fragments_df$fragment_length))))
  
  return(list(
    validation = validation_results,
    length_comparison = length_comparison
  ))
}

# Main function
generate_ctcf_fragments <- function(ctcf_file, freq_file, 
                                    fragments_per_region = 1, 
                                    output_prefix = "random_fragments",
                                    seed = 123) {
  message("Starting random fragment generation...")
  message(paste("CTCF file:", ctcf_file))
  message(paste("Length frequency file:", freq_file))
  message(paste("Fragments per region:", fragments_per_region))
  message(paste("Random seed:", seed))
  
  # 1. Read input files
  input_data <- read_input_files(ctcf_file, freq_file)
  ctcf_data <- input_data$ctcf_data
  freq_data <- input_data$freq_data
  
  message(paste("Read", nrow(ctcf_data), "CTCF regions"))
  message(paste("Read", nrow(freq_data), "fragment length types"))
  
  # 2. Generate random fragments
  fragments_df <- generate_random_fragments(ctcf_data, freq_data, 
                                            fragments_per_region, seed)
  
  if (is.null(fragments_df)) {
    message("Failed to generate any fragments")
    return(NULL)
  }
  
  message(paste("Successfully generated", nrow(fragments_df), "fragments"))
  
  # 3. Save results
  output_files <- save_results(fragments_df, output_prefix)
  
  # 4. Validate results
  validation_results <- validate_fragments(fragments_df, ctcf_data, freq_data)
  
  # 5. Return results
  result <- list(
    fragments = fragments_df,
    output_files = output_files,
    validation = validation_results,
    ctcf_data = ctcf_data,
    freq_data = freq_data
  )
  
  message("\n=== Generation Complete ===")
  return(result)
}