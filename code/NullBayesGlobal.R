# script for generating null distributions of conditional probabilities
# of having isoniazid resistance given a particular random mutation
# from a set of completely random mutations across all isolates

# load all libraries
{ 
  library(tidyverse)
  library(dplyr)
  library(reshape2)
}
#################################################################################
# setting command line arguments for adding job and array ids to file name
args <- commandArgs(TRUE)
job.id <- args[1]
array.id <- args[2]
#################################################################################
# writing a function to calculate bayes probability
calculateBayes <- function(total.n, total.n.with.a, total.n.with.b, total.n.with.a.and.b) {
  # Calculate the conditional probability P(A|B)
  bayes.prob <- total.n.with.a.and.b / total.n.with.b
  return(bayes.prob)
}
# function for using apply to get bayes probs for each locus
calculateMetrics <- function(count.table, phenotypes) {
  # Get the sample IDs with mutations (non-zero counts)
  isomut <- colnames(count.table)[apply(count.table, 2, function(x) any(x > 0))]
  isolate.with.mutation.n <- length(isomut)
  # Get the count of resistant isolates with mutations
  R.isolate.n.with.mutation <- sum(resistant_isolates %in% isomut)
  # Calculate the Bayes probability
  bayes.prob <- calculateBayes(isolate.n, R.isolate.n, isolate.with.mutation.n, R.isolate.n.with.mutation)
  return(bayes.prob)
}
# function for converting matrix of isolate IDs and mutations to wide form
flattenCounts <- function(isolate.locus.matrix) {
  # transforming table to wide form with SNP occurrences
  count.table.wide <- isolate.locus.matrix %>%
    group_by(Sample.ID, Gene.ID) %>%
    summarise(n = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Sample.ID, values_from = n, values_fill = list(n = 0))
  # convert to dataframe
  count.table.wide <- as.data.frame(count.table.wide)
  count.table.wide <- count.table.wide[complete.cases(count.table.wide), ]
  # set rownames and remove locus column
  rownames(count.table.wide) <- count.table.wide$Gene.ID
  count.table.wide <- count.table.wide %>% dplyr::select(-Gene.ID)
  return(count.table.wide)
}
#################################################################################
# reading in snp table, clinical metadata, and list of genes in Mtb
genotypes <- read_csv(gzfile('/proj/omics4tb2/epepper/AMR/manuscript-data/LDF/all-isolate-locus-matrix.csv.gz'))
phenotypes <- read_csv(file = '/proj/omics4tb2/epepper/AMR/manuscript-data/tb-profiler-calls.csv')
#################################################################################
# setting random sample to size of 5000
isolate.n <- 5000
#################################################################################
# performing randomization of mutation counts and calculating bayes on each
n.iterations <- 500
for (x in 1:n.iterations) {
  # start timer
  start.time <- proc.time()
  # save iteration number
  iteration.number <- x
  # randomly sampling data to collect 5000 isolates and perform bayes calculations on their mutations
  randomly.sampled.isolates <- sample(x = unique(phenotypes$Sample.ID), size = isolate.n, replace = F)
  # sampling genotypes and phenotypes
  genotypes.sample <- genotypes %>% filter(Sample.ID %in% randomly.sampled.isolates)
  phenotypes.sample <- phenotypes %>% filter(Sample.ID %in% randomly.sampled.isolates)
  # getting list of all locus tags we have information on
  locus.tags <- data.frame('locus' = unique(genotypes.sample$Gene.ID))
  # getting table of counts for each isolate locus tag pair
  snp.data.wide <- flattenCounts(genotypes.sample)
  # Get the number of isolates that are resistant to Isoniazid
  resistant_isolates <- phenotypes.sample %>%
    filter(Sample.ID %in% colnames(snp.data.wide) & Isoniazid == 1) %>%
    pull(Sample.ID)
  R.isolate.n <- length(resistant_isolates)
  # randomization step
  snp.data.wide.random <- data.frame(sapply(1:ncol(snp.data.wide),
                                            function(x){snp.data.wide[sample(1:nrow(snp.data.wide), nrow(snp.data.wide)), x]}))
  # reset rownames and colnames on randomized table
  rownames(snp.data.wide.random) <- locus.tags[[1]]
  colnames(snp.data.wide.random) <- colnames(snp.data.wide)
  #################################################################################
  # now calculating null conditional probabilities for each gene
  system.time(bayes.null <- data.frame(sapply(1:nrow(snp.data.wide.random),
                                              function(x){ calculateMetrics(snp.data.wide.random[x, ], phenotypes.sample) })))
  # resetting columns names
  bayes.null$locus.tag <- rownames(snp.data.wide.random)
  rownames(bayes.null) <- NULL
  names(bayes.null) <- c('bayes', 'locus')
  # create a table that sampled loci, and calculated bayes probability and add table from one iteration to list
  compute.time <- proc.time() - start.time
  print(paste('done with iteration ', iteration.number, ' taking ',
              round((compute.time[[3]] / 60), digits = 2), ' minutes', sep = ''))
  write_csv(bayes.null, file = paste('/proj/omics4tb2/epepper/AMR/manuscript-data/LDF/bayes-null-distributions//INH-global/bayes-null-permutation-', x, '-', job.id, '.', array.id, '.csv', sep = ''))
  bayes.null <- NULL
}