library(MALDIquant)
library(MALDIquantForeign)
library(ggplot2)
library(stringr)

# Update this with the path to your working directory containing data.
setwd('path_to_data')

# Load Data
# Update this with the path to your raw data.
data_dir <- file.path('path_to_MALDI_data_from_MSV000091443')
spectra_list_raw <- import(data_dir)

# QC
any(sapply(spectra_list_raw, isEmpty))  # FALSE
table(sapply(spectra_list_raw, length))
all(sapply(spectra_list_raw, isRegular))  # FALSE

# Preprocessing
# Parameters can be updated for different instruments.
spectra_list <- trim(spectra_list_raw, c(4000, 20000))  # Update mass range as needed.
spectra_list <- transformIntensity(spectra_list, method='sqrt')
spectra_list <- smoothIntensity(spectra_list, method='SavitzkyGolay', halfWindowSize=10)
spectra_list <- removeBaseline(spectra_list, method='SNIP', iterations=100)
spectra_list <- calibrateIntensity(spectra_list, method='TIC')

# Spectra Averaging - Average technical replicate spectra for each sample.
avg_spectra_list <- list()
count <- 1
# Replicate numbers should be updated as needed for other datasets.
# Current numbers are for 54 biological samples with 48 technical replicates each.
for (i in seq(1, 2592, 48)) {
  avg_spectra_list[[count]] <- averageMassSpectra(spectra_list[i:i+47], method='mean')
  count <- count + 1
}

# Peak Detection
peaks <- detectPeaks(avg_spectra_list, method='MAD', halfWindowSize=20, SNR=3)
peaks <- binPeaks(peaks, tolerance=0.2)

# Create Feature Matrix
feature_matrix <- intensityMatrix(peaks, avg_spectra_list)
filenames <- basename(unlist(lapply(peaks, function(x) attributes(x)$metaData$file)))
attributes(feature_matrix)$dimnames[[1]] <- str_split_i(filenames, pattern='_', 1)
feature_df <- as.data.frame(feature_matrix)

# Load Clinical Annotations for Primary Cancer
annotations <- read.csv('Tampon_Annotations.csv')
feature_df$Class <- annotations$Primary_Cancer

# Subset data to only include "Ovarian Cancer" and "Benign" samples.
ovarian <- feature_df[!is.na(feature_df$Class),]
ovarian <- ovarian[(ovarian$Class == 'Ovarian' | ovarian$Class == 'Benign'),]

# Write out feature matrix to CSV file.
write.csv(ovarian, file='TP_Ovarian.csv')
