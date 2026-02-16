# Shared functions for sRange analysis
# This file contains reusable functions for both canids and penguins analysis

#' Parse XML sampling dates and extract species information
#'
#' @param text XML text containing samplingDates elements
#' @return A list containing:
#'   - priors: data.frame with species, start, and end columns
#'   - ranges_single: vector of species with single samples
#'   - ranges: vector of species with paired samples (first/last)
#'   - ranges_: vector of all range names without suffixes
parse_sampling_dates <- function(text) {
  # Extract taxon names without the "@", and the lower and upper values
  sp <- gsub('@', '', regmatches(text, gregexpr('(?<=taxon="@)[^"]+', text, perl=TRUE))[[1]])
  start <- as.numeric(regmatches(text, gregexpr('(?<=upper=")[^"]+', text, perl=TRUE))[[1]])
  end <- as.numeric(regmatches(text, gregexpr('(?<=lower=")[^"]+', text, perl=TRUE))[[1]])

  # Remove '_first' or '_last' from taxon names
  ranges_ <- gsub('_(first|last)$', '', sp)

  priors <- data.frame(species=sp, start=start, end=end)

  all_string_listwise <- unlist(lapply(ranges_, unique))
  ranges_single <- names(which(table(all_string_listwise)==1))
  ranges <- names(which(table(all_string_listwise)>1))

  return(list(
    priors = priors,
    ranges_single = ranges_single,
    ranges = ranges,
    ranges_ = ranges_
  ))
}


#' Create species data frame for CSV output
#'
#' @param ranges_single Vector of single-sampled species
#' @param ranges Vector of paired-sampled species
#' @param sp Full species names with suffixes
#' @param start Start dates
#' @param end End dates
#' @return data.frame with species range information
create_species_dataframe <- function(ranges_single, ranges, sp, start, end, ranges_) {
  id_first <- which(ranges_ %in% ranges_single)
  sp_first <- data.frame(
    "Species" = ranges_[id_first],
    "First start" = start[id_first],
    "First end" = end[id_first],
    "Last start" = rep("-", length(id_first)),
    "Last end" = rep("-", length(id_first))
  )

  sp_both <- data.frame(
    "Species" = c(),
    "First start" = c(),
    "First end" = c(),
    "Last start" = c(),
    "Last end" = c()
  )

  for (r in ranges) {
    id1 <- which(sp == paste0(r, "_first"))
    id2 <- which(sp == paste0(r, "_last"))
    sp_both <- rbind(
      sp_both,
      data.frame(
        "Species" = r,
        "First start" = start[id1],
        "First end" = end[id1],
        "Last start" = start[id2],
        "Last end" = end[id2]
      )
    )
  }

  return(rbind(sp_first, sp_both))
}


#' Calculate IQR statistics for species ranges
#'
#' @param log Log data frame with posterior samples
#' @param priors Prior data frame with species, start, end
#' @param ranges_single Vector of single-sampled species
#' @param ranges Vector of paired-sampled species
#' @return List containing ranges_full_s, ranges_full_e, ranges_full, ranges_stat
get_sRange_full <- function(log, priors, ranges_single, ranges) {
  start_r <- c()
  start_r_low <- c()
  start_r_high <- c()
  start_r_min <- c()
  start_r_max <- c()
  start_r_iqr <- c()
  start_r_iqr_prior <- c()
  end_r <- c()
  end_r_low <- c()
  end_r_high <- c()
  end_r_min <- c()
  end_r_max <- c()
  end_r_iqr <- c()
  end_r_iqr_prior <- c()

  distr <- c()
  endpoint <- c()
  range <- c()
  range_s <- c()
  range_e <- c()
  median_s <- c()
  median_s_prior <- c()
  median_e <- c()
  median_e_prior <- c()

  # Process single-sampled species
  for (r in ranges_single) {
    if ((max(log[, paste0(r, "_first")]) - min(log[, paste0(r, "_first")])) > 0.1) {
      start_l <- length(log[, paste0(r, "_first")])
      distr <- c(distr, log[, paste0(r, "_first")])
      endpoint <- c(endpoint, rep("start", start_l))
      range <- c(range, rep(r, start_l))

      id <- which(priors$species == paste0(r, "_first"))
      un <- runif(5000, min=priors$end[id], max=priors$start[id])
      start_r_iqr <- c(start_r_iqr, IQR(log[, paste0(r, "_first")]))
      start_r_iqr_prior <- c(start_r_iqr_prior, IQR(un))

      range_s <- c(range_s, r)
      median_s <- c(median_s, median(log[, paste0(r, "_first")]))
      median_s_prior <- c(median_s_prior, median(un))
    }

    start_r <- c(start_r, median(log[, paste0(r, "_first")]))
    end_r <- c(end_r, NA)
    hpd_f <- HPDinterval(as.mcmc(log[, paste0(r, "_first")]))
    start_r_low <- c(start_r_low, hpd_f[1])
    start_r_high <- c(start_r_high, hpd_f[2])
    start_r_min <- c(start_r_min, min(log[, paste0(r, "_first")]))
    start_r_max <- c(start_r_max, max(log[, paste0(r, "_first")]))

    end_r_low <- c(end_r_low, NA)
    end_r_high <- c(end_r_high, NA)
    end_r_min <- c(end_r_min, NA)
    end_r_max <- c(end_r_max, NA)
  }

  # Process paired-sampled species
  for (r in ranges) {
    if ((max(log[, paste0(r, "_first")]) - min(log[, paste0(r, "_first")])) > 0.1) {
      start_l <- length(log[, paste0(r, "_first")])
      distr <- c(distr, log[, paste0(r, "_first")])
      endpoint <- c(endpoint, rep("start", start_l))
      range <- c(range, rep(r, start_l))

      id <- which(priors$species == paste0(r, "_first"))
      un <- runif(5000, min=priors$end[id], max=priors$start[id])
      start_r_iqr <- c(start_r_iqr, IQR(log[, paste0(r, "_first")]))
      start_r_iqr_prior <- c(start_r_iqr_prior, IQR(un))
      range_s <- c(range_s, r)
      median_s <- c(median_s, median(log[, paste0(r, "_first")]))
      median_s_prior <- c(median_s_prior, median(un))
    }

    if ((max(log[, paste0(r, "_last")]) - min(log[, paste0(r, "_last")])) > 0.1) {
      end_l <- length(log[, paste0(r, "_last")])
      distr <- c(distr, log[, paste0(r, "_last")])
      endpoint <- c(endpoint, rep("end", end_l))
      range <- c(range, rep(r, end_l))

      id <- which(priors$species == paste0(r, "_last"))
      un <- runif(5000, min=priors$end[id], max=priors$start[id])
      end_r_iqr <- c(end_r_iqr, IQR(log[, paste0(r, "_last")]))
      end_r_iqr_prior <- c(end_r_iqr_prior, IQR(un))
      range_e <- c(range_e, r)
      median_e <- c(median_e, median(log[, paste0(r, "_last")]))
      median_e_prior <- c(median_e_prior, median(un))
    }

    start_r <- c(start_r, median(log[, paste0(r, "_first")]))
    end_r <- c(end_r, median(log[, paste0(r, "_last")]))

    hpd_f <- HPDinterval(as.mcmc(log[, paste0(r, "_first")]))
    start_r_low <- c(start_r_low, hpd_f[1])
    start_r_high <- c(start_r_high, hpd_f[2])
    start_r_min <- c(start_r_min, min(log[, paste0(r, "_first")]))
    start_r_max <- c(start_r_max, max(log[, paste0(r, "_first")]))

    hpd_l <- HPDinterval(as.mcmc(log[, paste0(r, "_last")]))
    end_r_low <- c(end_r_low, hpd_l[1])
    end_r_high <- c(end_r_high, hpd_l[2])
    end_r_min <- c(end_r_min, min(log[, paste0(r, "_last")]))
    end_r_max <- c(end_r_max, max(log[, paste0(r, "_last")]))
  }

  ranges_full_s <- data.frame(
    range = gsub("_", " ", range_s),
    start_r_iqr = start_r_iqr,
    start_r_iqr_prior = start_r_iqr_prior,
    median_s_prior = median_s_prior,
    median_s = median_s
  )

  ranges_full_e <- data.frame(
    range = gsub("_", " ", range_e),
    end_r_iqr = end_r_iqr,
    end_r_iqr_prior = end_r_iqr_prior,
    median_e_prior = median_e_prior,
    median_e = median_e
  )

  ranges_full <- data.frame(distr=distr, endpoint=endpoint, range=range)

  length(end_r) <- length(start_r)
  ranges_stat <- data.frame(
    range = c(ranges_single, ranges),
    median_start = start_r,
    median_end = end_r,
    start_r_min = start_r_min,
    start_r_low = start_r_low,
    start_r_max = start_r_max,
    start_r_high = start_r_high,
    end_r_min = end_r_min,
    end_r_low = end_r_low,
    end_r_max = end_r_max,
    end_r_high = end_r_high
  )

  return(list(ranges_full_s, ranges_full_e, ranges_full, ranges_stat))
}


#' Create IQR ratio plot
#'
#' @param ranges_full_s_comb Combined start ranges dataframe
#' @param ranges_full_e_comb Combined end ranges dataframe
#' @param y_scale_params Optional list with breaks, limits, labels for y-axis
#' @param show_legend Logical, whether to show legend
#' @param text_size Text size for theme
#' @return ggplot object
create_iqr_plot <- function(ranges_full_s_comb, ranges_full_e_comb,
                            y_scale_params = NULL,
                            show_legend = FALSE,
                            text_size = 21) {
  p <- ggplot(ranges_full_s_comb, aes(x = range, y = start_r_iqr_prior/start_r_iqr)) +
    geom_point(aes(color = "First", shape=model), size = 3, alpha=0.7) +
    geom_point(data=ranges_full_e_comb,
               aes(x = range, y = end_r_iqr_prior/end_r_iqr, color = "Last", shape=model),
               size = 3, alpha=0.7) +
    geom_hline(yintercept = 1) +
    ylab("IQR ratio (prior/posterior)") +
    xlab("Species") +
    scale_color_manual(values = c("First" = "#FF7F50", "Last" = "#43AA8B")) +
    theme_minimal() +
    labs(color="Sample", shape="Morph. data attached to") +
    theme(
      legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid"),
      legend.justification = c("right", "bottom"),
      legend.box = "horizontal",
      text = element_text(size=text_size),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    ) +
    guides(color = guide_legend(override.aes = list(shape=15, size=5)))

  if (!show_legend) {
    p <- p + theme(legend.position = "None")
  }

  if (!is.null(y_scale_params)) {
    p <- p + scale_y_continuous(
      breaks = y_scale_params$breaks,
      limits = y_scale_params$limits,
      labels = y_scale_params$labels
    )
  }

  return(p)
}


#' Process sRange analysis workflow
#'
#' @param log_both Log data from "both" model
#' @param log_first Log data from "first" model
#' @param priors Prior dataframe
#' @param ranges_single Single-sampled species
#' @param ranges Paired-sampled species
#' @return List with ranges_full_s_comb and ranges_full_e_comb
process_srange_workflow <- function(log_both, log_first, priors, ranges_single, ranges) {
  all_both <- get_sRange_full(log_both, priors, ranges_single, ranges)
  ranges_full_s <- all_both[[1]]
  ranges_full_e <- all_both[[2]]

  all_first <- get_sRange_full(log_first, priors, ranges_single, ranges)
  ranges_full_s1 <- all_first[[1]]
  ranges_full_e1 <- all_first[[2]]

  ranges_full_s1$model <- "First"
  ranges_full_e1$model <- "First"

  ranges_full_s$model <- "Both"
  ranges_full_e$model <- "Both"

  ranges_full_s_comb <- rbind(ranges_full_s1, ranges_full_s)
  ranges_full_e_comb <- rbind(ranges_full_e1, ranges_full_e)

  return(list(
    ranges_full_s_comb = ranges_full_s_comb,
    ranges_full_e_comb = ranges_full_e_comb
  ))
}
