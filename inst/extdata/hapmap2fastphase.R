## Convert gziped HapMap file into simplified fastPHASE inp file. The output is
## filename.inp, where `filename` is a basename extracted from input
## filename.txt.gz file. Nonbiallelic sites are filtered out.
## Usage: Rscript hapmap2fastphase.R filename.txt.gz unrelated.txt
## Author: Gennady Khvorykh, inzilico.com

isBi <- function(x) {
  tmp <- unlist(strsplit(x, split = ""))
  tmp <- tmp[tmp != "?"]
  length(unique(tmp)) <= 2
}

args <- commandArgs(TRUE)

# Check input
if(is.na(args[1]) | is.na(args[2]))
  stop("Missing input.\nUsage: Rscript hapmap2fastphase.R filename.txt.gz unrelated.txt", call. = F )

if(!file.exists(args[1])) stop(args[1], " doesn't exist!")
if(!file.exists(args[2])) stop(args[2], " doesn't exist!")

## TODO: Remove after debug
# args <- c("~/data/hapmap/phaseI/genotypes_chr22_CEU.txt.gz",
#           "~/data/hapmap/phaseI/unrelated.txt")

# Get basename
fn <- basename(args[1])
fn <- sub(".txt.gz", "", fn)
fn <- paste0(fn, ".inp")

# Show input
message("HapMap file: ", args[1])
message("Converting...")

# Load data
cmd <- sprintf('zcat -cq %s | cut -d" " -f-5,12-', args[1])
data <- system(cmd, intern = T)

# Transform into matrix
data <- strsplit(data, " ")
data <- do.call(rbind, data)

# Get sample ids
ids <- data[1, -c(1:5)]

# Remove first line containing column names
data <- data[-1, ]

# Subset observations with strand +
ind <- data[, 5] == "+"
data <- data[ind, ]

# Subset genotypes
data <- data[, -c(1:5)]

# Subset unrelated individuals
unrelated <- read.table(args[2], stringsAsFactors = F)[, 2]
data <- data[, ids %in% unrelated]

# Arrange haplotypes
hap <- apply(data, 2, strsplit, split = "")
hap <- lapply(hap, do.call, what = "rbind")
hap <- do.call(cbind, hap)

# Replace any character except A, T, C, G by 0 (missing allele)
hap <- apply(hap, 1, gsub, pattern = "[^ATCG]", replacement = "?")

# Subset biallelic sites
bi <- apply(hap, 2, isBi)
hap <- hap[, bi]

if (sum(!bi) > 0) message("Non biallelic: ", sum(!bi))

nhap <- nrow(hap)
nind <- nhap/2
nsnp <- ncol(hap)

message("Haplotypes: ", nhap)
message("SNPs: ", nsnp)

# Save as fastPHASE inp simplified file
write(nind, fn)
write(nsnp, fn, append = T)
write(t(hap), fn, ncolumns = nsnp, append = T, sep = "")

message("All done!")




