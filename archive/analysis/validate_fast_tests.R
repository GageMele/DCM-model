# Meaningful checks for the new subset logic, manifests, and sampler CSV reader.
package <- "DCM_Ibex_fast_tests"
for (path in list.files(package, pattern = "\\.R$", full.names = TRUE)) parse(path)
source(file.path(package, "common.R"))
cases <- read.delim(file.path(package, "cases.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(cases) == 24L, !anyDuplicated(cases$case_id),
  sum(cases$mode == "posterior") == 22L, sum(cases$mode == "prior") == 2L,
  all(cases$chains >= 2), all(cases$iter > cases$warmup),
  all(cases$cpus == cases$chains * pmax(1, cases$threads)))
stopifnot(sum(cases$cpus) == 58L)
baseline <- cases[cases$case_id == 3, ]
thread_cases <- cases[cases$case_id %in% c(3, 21, 22), ]
for (field in c("subset", "prior", "mode", "chains", "iter", "warmup", "adapt_delta", "treedepth"))
  stopifnot(length(unique(thread_cases[[field]])) == 1L)

# Duplicate cast-level metadata is represented correctly; support is retained.
metadata <- data.frame(cast_key = sprintf("cast_%02d", 1:30),
  log10_distance = rep(1:10, 3), MLD_temp_MIX = rep(c(4, 8, 13, 20, 31, 45), 5))
order_one <- supporting_cast_order(metadata, 123)
order_two <- supporting_cast_order(metadata, 123)
stopifnot(identical(order_one, order_two), !anyDuplicated(order_one$ids),
  setequal(order_one$ids, metadata$cast_key))
support <- metadata[match(head(order_one$ids, order_one$minimum), metadata$cast_key), ]
stopifnot(length(unique(support$log10_distance)) >= 5,
          length(unique(support$MLD_temp_MIX)) >= 5)
unsupported <- metadata
unsupported$MLD_temp_MIX <- 1
stopifnot(inherits(try(supporting_cast_order(unsupported, 123), silent = TRUE), "try-error"))

# Allocation is nested and sums to the requested target when feasible.
sizes <- c(A = 50L, B = 250L, C = 700L)
lower <- c(A = 5L, B = 6L, C = 5L)
small <- allocate_cast_counts(sizes, lower, .05)
middle <- allocate_cast_counts(sizes, small, .10)
large <- allocate_cast_counts(sizes, middle, .20)
stopifnot(sum(small) == 50, sum(middle) == 100, sum(large) == 200,
  all(small >= lower), all(middle >= small), all(large >= middle), all(large <= sizes))
too_small <- allocate_cast_counts(sizes, lower, .01)
stopifnot(sum(too_small) == sum(lower))

# Equal cast weighting prevents a long cast from dominating the prior location.
toy <- data.frame(cast_key = c(rep("long", 100), "short"),
                  chl_norm = c(rep(0, 100), 100))
stats <- full_response_statistics(toy)
stopifnot(stats$casts == 2, stats$empirical_scale >= .5,
          stats$cast_balanced_q75 == 100)
stopifnot(is.na(finite_min(c(400, NA))), is.na(finite_max(numeric())))

# Read warmup and sampling diagnostics while skipping parameter columns.
dir.create("validation_fast_tests", showWarnings = FALSE)
csv <- file.path("validation_fast_tests", "mock_chain.csv")
writeLines(c("# num_warmup = 2", "# save_warmup = 1",
  "lp__,accept_stat__,stepsize__,treedepth__,n_leapfrog__,divergent__,energy__,b_example",
  "-2,0.8,0.1,4,15,1,2,20", "-1,0.9,0.1,3,7,0,3,21",
  "# Adaptation terminated", "-1,0.95,0.05,2,3,0,2.5,22",
  "-1,0.96,0.05,2,3,0,3.5,23", "# Elapsed Time: 15.4 seconds (Warm-up)",
  "#               2.6 seconds (Sampling)"), csv)
sampler <- read_sampler_csv(csv, 4L, 2L)
times <- chain_elapsed(csv)
stopifnot(nrow(sampler) == 4L, !"b_example" %in% names(sampler),
  identical(sampler$phase, c("warmup", "warmup", "sampling", "sampling")),
  times$warmup_seconds == 15.4, times$sampling_seconds == 2.6)

# Confirm the original scientific specification can be identified unambiguously.
original <- new.env(parent = globalenv())
sys.source("DCM_Ibex_workflow/02_fit_model.R", envir = original)
expressions <- as.list(body(original$fit_model))[-1]
for (name in c("andrew_form", "priors")) {
  matches <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
    identical(x[[2]], as.name(name)), expressions)
  stopifnot(length(matches) == 1L)
}
cat("R syntax, 24-case design, support selection, nested allocation, weighting,\n",
    "diagnostic CSV parsing, and original-specification checks passed.\n")
