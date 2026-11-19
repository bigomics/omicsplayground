## test-module-heatmap-labels.R
##
## Regression guard for issue #1642: feature names of the multi-omics heatmap
## (Multi-Omics -> Factor heatmap / Multi-omics heatmap) are long semicolon
## joined gene symbol lists and used to take over the whole plot.
## playbase::mofa.plot_heatmap shortened the short feature IDs before mapping
## them onto those long names, so the labels drawn were never shortened.

.board_dir <- if (dir.exists("components/board.mofa")) {
  "components/board.mofa"
} else {
  "../../components/board.mofa"
}
source(file.path(.board_dir, "R", "plot_moduleheatmap.R"), encoding = "UTF-8", local = TRUE)

long_symbol <- "H3C1;H3C4;H3C3;H3C6;H3C11;H3C8;H3C12;H3C10;H3C2;H3C7"

shorten <- function(maxchar) {
  features <- c("gene:f1", "gene:f2", "protein:p1", "protein:p2")
  x <- matrix(0, 4, 2, dimnames = list(features, c("A", "B")))
  genes <- data.frame(
    feature = features,
    symbol = c("H3C1", "H2AC15", long_symbol, "CLIC4"),
    gene_name = c(long_symbol, "H2AC15;H2AC14", long_symbol, "P045;CLIC4"),
    data_type = c("gene", "gene", "protein", "protein"),
    row.names = features,
    stringsAsFactors = FALSE
  )
  rownames(moduleheatmap_shorten_labels(x, maxchar = maxchar, gene_table = genes))
}

test_that("long feature names are shortened with an ellipsis", {
  lab <- shorten(maxchar = 30)

  expect_true(all(nchar(sub("^[A-Za-z]+:", "", lab)) <= 30))
  expect_true(any(grepl("H3C1;H3C4;H3C3;H3C6;H3", lab, fixed = TRUE)))
  expect_true(all(grepl("\\.\\.\\.$", lab[grepl("H3C1;H3C4", lab)])))
})

test_that("short feature names are left untouched", {
  lab <- shorten(maxchar = 80)

  expect_setequal(
    lab,
    c(
      "gene:H3C1;H3C4;H3C3;H3C6;H3C11;H3C8;H3C12;H3C10;H3C2;H3C7", "gene:H2AC15;H2AC14",
      "protein:H3C1;H3C4;H3C3;H3C6;H3C11;H3C8;H3C12;H3C10;H3C2;H3C7", "protein:P045;CLIC4"
    )
  )
  expect_false(any(grepl("\\.\\.\\.", lab)))
})

test_that("data-type prefixes and unique row names are kept", {
  lab <- shorten(maxchar = 10)

  expect_true(all(grepl("^gene:|^protein:", lab)))
  expect_false(any(duplicated(lab)))

  ## identical gene names in different data types stay distinct labels
  x <- matrix(0, 2, 2, dimnames = list(c("gene:f1", "protein:p1"), c("A", "B")))
  genes <- data.frame(
    feature = c("gene:f1", "protein:p1"),
    symbol = long_symbol, gene_name = long_symbol,
    data_type = c("gene", "protein"),
    row.names = c("gene:f1", "protein:p1"), stringsAsFactors = FALSE
  )
  lab <- rownames(moduleheatmap_shorten_labels(x, maxchar = 10, gene_table = genes))
  expect_setequal(lab, c("gene:H3...", "protein:H3"))
  ## the limit covers the whole label, prefix included, for every label
  expect_true(all(nchar(lab) <= 10))
})

test_that("disambiguation keeps the data-type prefix off the chopping block", {
  ## CodeRabbit regression: one duplicate pair used to truncate the whole
  ## label, prefix included, yielding "ge..." and breaking the type split.
  features <- c("gene:f1", "gene:f2")
  x <- matrix(0, 2, 2, dimnames = list(features, c("A", "B")))
  genes <- data.frame(
    feature = features,
    symbol = c("H3C1", "H3C1"),
    gene_name = c(long_symbol, long_symbol),
    data_type = c("gene", "gene"),
    row.names = features, stringsAsFactors = FALSE
  )
  for (mc in c(10, 12, 15, 30)) {
    lab <- rownames(moduleheatmap_shorten_labels(x, maxchar = mc, gene_table = genes))
    expect_true(all(grepl("^gene:", lab)), paste("maxchar", mc))
    expect_equal(anyDuplicated(lab), 0L, info = paste("maxchar", mc))
    expect_true(all(nchar(lab) <= mc), paste("maxchar", mc))
  }
})

## Two collisions inside one data type: the same gene name twice, and two
## different names that become identical once truncated. Heatmaps reject
## duplicate row names, so both must be disambiguated.
collide <- function(maxchar = 20) {
  features <- c("gene:f1", "gene:f2", "gene:f3", "gene:f4", "protein:p1")
  x <- matrix(0, 5, 2, dimnames = list(features, c("A", "B")))
  ## two gene rows share a symbol, and two more gene rows differ only after
  ## the 20-character limit, so both pairs collide once shortened
  genes <- data.frame(
    feature = features,
    symbol = c("H3C1", "H3C1", "ABCDEFGHIJKLMNOPQRST-one", "ABCDEFGHIJKLMNOPQRST-two", "CLIC4"),
    gene_name = c(
      long_symbol, long_symbol,
      "ABCDEFGHIJKLMNOPQRST-one", "ABCDEFGHIJKLMNOPQRST-two", "CLIC4"
    ),
    data_type = c("gene", "gene", "gene", "gene", "protein"),
    row.names = features, stringsAsFactors = FALSE
  )
  rownames(moduleheatmap_shorten_labels(x, maxchar = maxchar, gene_table = genes))
}

test_that("labels that collide within one data type are made unique", {
  lab <- collide(maxchar = 20)

  ## every gene row is shortened to the limit and the collisions disambiguated
  gene_lab <- sub("^gene:", "", lab[grepl("^gene:", lab)])
  expect_length(gene_lab, 4)
  expect_false(any(duplicated(gene_lab)))
  expect_true(all(grepl("^gene:|^protein:", lab)))

  ## Collisions are resolved by shortening further, so make the assertion
  ## about the property, not a particular suffix shape.
  h3 <- gene_lab[grepl("^H3C1", gene_lab)]
  expect_length(h3, 2)
  expect_false(h3[1] == h3[2])

  ## the two names that only differ past the limit stay distinguishable
  ab <- gene_lab[grepl("^ABCDEF", gene_lab)]
  expect_length(ab, 2)
  expect_false(ab[1] == ab[2])
})

test_that("a very small maxchar still keeps prefixes and unique names", {
  lab <- collide(maxchar = 4)

  expect_false(any(duplicated(lab)))
  ## the data-type prefix survives the floor, so rows still split by type
  expect_true(all(grepl("^gene:|^protein:", lab)))
})

test_that("a prefix longer than the limit keeps the data type", {
  ## "metabolite:" is 11 characters, so nothing can fit at maxchar = 10. The
  ## prefix is kept rather than silently dropping the data type.
  features <- c("metabolite:c1")
  x <- matrix(0, 1, 2, dimnames = list(features, c("A", "B")))
  genes <- data.frame(
    symbol = "c1", gene_name = "Cholesterol", data_type = "metabolite",
    stringsAsFactors = FALSE
  )
  genes$feature <- features
  rownames(genes) <- features
  expect_identical(
    rownames(moduleheatmap_shorten_labels(x, maxchar = 10, gene_table = genes)),
    "metabolite:"
  )
  ## one character of room is enough to shorten again
  expect_identical(
    rownames(moduleheatmap_shorten_labels(x, maxchar = 12, gene_table = genes)),
    "metabolite:C"
  )
})

test_that("many duplicates keep the suffix within the limit", {
  ## The compact " (n)" suffix cannot fit at maxchar = 10, and neither can
  ## " 10"; the labels must still stay unique and within the limit.
  features <- paste0("protein:p", 1:12)
  x <- matrix(0, 12, 2, dimnames = list(features, c("A", "B")))
  genes <- data.frame(
    feature = features, symbol = "ER",
    gene_name = rep("ER-alpha_pS118", 12),
    data_type = "protein",
    row.names = features, stringsAsFactors = FALSE
  )
  lab <- rownames(moduleheatmap_shorten_labels(x, maxchar = 10, gene_table = genes))

  expect_false(anyDuplicated(lab) > 0)
  expect_true(all(grepl("^protein:", lab)))
  expect_true(all(nchar(lab) <= 10), paste(lab, collapse = " | "))
})

test_that("numbered labels stay within the limit while a number fits", {
  ## "protein:" leaves 2 characters at maxchar = 10, so single- and
  ## double-digit numbers fit. Up to 99 duplicates stay unique and in limit.
  dup_labels <- function(n, maxchar = 10) {
    features <- paste0("protein:p", seq_len(n))
    x <- matrix(0, n, 2, dimnames = list(features, c("A", "B")))
    genes <- data.frame(
      symbol = "ER", gene_name = rep("ER-alpha", n), data_type = "protein",
      stringsAsFactors = FALSE
    )
    genes$feature <- features
    rownames(genes) <- features
    rownames(moduleheatmap_shorten_labels(x, maxchar = maxchar, gene_table = genes))
  }

  lab <- dup_labels(99)
  expect_false(anyDuplicated(lab) > 0)
  expect_true(all(nchar(lab) <= 10))

  ## A three-digit number cannot fit in the remaining budget. Uniqueness is
  ## kept, which is what a heatmap requires; only that row exceeds the limit.
  lab <- dup_labels(101)
  expect_false(anyDuplicated(lab) > 0)
  expect_equal(sum(nchar(lab) > 10), 1L)
})

test_that("long-prefix labels survive the truncation playbase applies", {
  ## playbase::mofa.plot_heatmap() truncates rownames to maxchar again and then
  ## derives the data-type split from the truncated text. At maxchar = 10 an
  ## 11-character "metabolite:" prefix would become "metabol...", losing the
  ## type and making the two rows identical. The server therefore asks playbase
  ## for no further truncation once the labels are prepared.
  features <- c("metabolite:m1", "metabolite:m2")
  x <- matrix(0, 2, 2, dimnames = list(features, c("A", "B")))
  genes <- data.frame(
    symbol = "m", gene_name = rep("Cholesterol", 2), data_type = "metabolite",
    stringsAsFactors = FALSE
  )
  genes$feature <- features
  rownames(genes) <- features
  prepared <- rownames(moduleheatmap_shorten_labels(x, maxchar = 10, gene_table = genes))

  ## what the server passes to playbase
  limit <- max(10, max(nchar(prepared)))
  through <- stringr::str_trunc(prepared, limit)

  expect_identical(through, prepared)
  expect_true(all(startsWith(through, "metabolite:")))
  expect_false(anyDuplicated(through) > 0)
  expect_identical(unique(sub(":.*", "", through)), "metabolite")
})

test_that("invalid maxchar or missing gene names leave the matrix alone", {
  features <- c("gene:f1", "gene:f2")
  x <- matrix(0, 2, 2, dimnames = list(features, c("A", "B")))

  expect_identical(moduleheatmap_shorten_labels(x, maxchar = NA), x)
  expect_identical(moduleheatmap_shorten_labels(x, gene_table = NULL), x)
  ## no gene table: feature IDs are shown as they are
  expect_identical(rownames(moduleheatmap_shorten_labels(x, gene_table = NULL)), features)
})
