#!/usr/bin/env Rscript
##
## Build the `human-dups` fixture: opg-exampledata/human-symbol plus 40 duplicated
## gene symbols. Duplicated feature names are the case that made the app's cleanX
## silently drop rows (999 -> 969 on the shipped salmon-ENSOKI), and no
## Human/RNA-seq fixture exercised it -- the app rejects salmon's organism
## ("Oncorhynchus kisutch" is not in its namespace list), so the browser path needs
## a Human one.
##
##   Rscript fixtures/make_human_dups.R [src_dir] [out_dir]
##
## The duplicated rows carry slightly different values (x1.37) so a collapse shows
## up in the data, not just in the row count.
##

argv <- commandArgs(trailingOnly = TRUE)
## Resolve our own directory from --file=, the same way the other harness scripts
## do, so the fixture lands next to this script regardless of the caller's cwd.
HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))

src <- if (length(argv) >= 1) argv[1] else
  "/home/massagno/bigomics/GitHub/opg-exampledata/human-symbol"
out <- if (length(argv) >= 2) argv[2] else file.path(HERE, "human-dups")
if (!dir.exists(src)) stop("source fixture not found: ", src)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

## NB: read as a MATRIX. A data.frame cannot hold duplicated row names -- rbind()
## silently uniquifies them and the fixture ends up with zero duplicates.
cnt <- as.matrix(read.csv(file.path(src, "counts.csv"), row.names = 1, check.names = FALSE))

set.seed(7)
pick <- sort(sample(nrow(cnt), 40))
dup <- cnt[pick, , drop = FALSE]
dup[] <- pmax(0, round(dup * 1.37))
cnt2 <- rbind(cnt, dup)
cnt2 <- cnt2[order(rownames(cnt2)), , drop = FALSE]

## Write the feature name as an explicit first column, so duplicates survive the
## round-trip (write.csv(row.names=TRUE) would be fine too, but this is explicit).
df <- data.frame(feature = rownames(cnt2), cnt2, check.names = FALSE)
write.csv(df, file.path(out, "counts.csv"), row.names = FALSE)

file.copy(file.path(src, "samples.csv"), file.path(out, "samples.csv"), overwrite = TRUE)
file.copy(file.path(src, "contrasts.csv"), file.path(out, "contrasts.csv"), overwrite = TRUE)
writeLines(c("organism: Human", "data_type: RNA-seq"), file.path(out, "info.yml"))

chk <- read.csv(file.path(out, "counts.csv"), check.names = FALSE)
message(sprintf("[make_human_dups] %d features, %d duplicated names -> %s",
  nrow(chk), sum(duplicated(chk[[1]])), out))
stopifnot(sum(duplicated(chk[[1]])) == 40)
