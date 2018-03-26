## Convert gziped HapMap file into standard text Plink ped file accompanied with
## map file. The output are filename.ped and filename.map, where `filename` is a
## basename extracted from input filename.txt.gz file.
## Usage: Rscript hapmap2plink.R filename.txt.gz
## Author: Gennady Khvorykh, inzilico.com

isBiallelic <- function(x) {
  tmp <- unlist(strsplit(x, split = ""))
  tmp <- tmp[tmp != 0]
  length(unique(tmp)) <= 2
}

args <- commandArgs(TRUE)

# Check input
if(is.na(args[1]))
  stop("Missing input.\nUsage: Rscript hapmap2plink.R filename.txt.gz", call. = F )

if(!file.exists(args[1])) stop(args[1], " doesn't exist!")

# Get basename
fn <- basename(args[1])
fn <- sub(".txt.gz", "", fn)

# Show input
message("HapMap file: ", args[1])
message("Converting...")

# Load data
cmd <- sprintf('zcat -cq %s | cut -d" " -f-5,12-', args[1])
data <- system(cmd, intern = T)

# Transform into matrix
data <- strsplit(data, " ")
data <- do.call(rbind, data)

# Create map data
map <- data[-1, c(3, 1, 4, 5)]

# Leave only observations with strand +
ind <- map[, 4] == "+"
map <- map[ind, ]

# Accomplish map data creation
map[, 1] <- tolower(map[, 1])
map <- cbind(map, 0)

chrom <- unique(map[, 1])
message("Chromosome: ", chrom)

# Create ped data
ped <- data[, -c(1:5)]
ids <- ped[1, ]
message("Individuals: ", length(ids))

ped <- ped[-1, ]
ped <- ped[ind, ]
ped <- t(ped)

# Replace any letter except A, T, C, G by 0 (missing allele)
ped <- apply(ped, 2, gsub, pattern = "[^ATCG]", replacement = "0")

# Leave only biallelic sites
biallelic <- apply(ped, 2, isBiallelic)
map <- map[biallelic, ]
ped <- ped[, biallelic]

if (sum(!biallelic) > 0) message("Non biallelic: ", sum(!biallelic))
message("SNPs: ", nrow(map))

#  Accomplish the creation of ped
ped <- cbind(0, ids, 0, 0, 0, -9, ped)

# Save ped/map
write.table(map[, c(1, 2, 5, 3)], paste0(fn, ".map"), col.names = F, row.names = F, quote = F)
write.table(ped, paste0(fn, ".ped"), col.names = F, row.names = F, quote = F)

message("All done!")




