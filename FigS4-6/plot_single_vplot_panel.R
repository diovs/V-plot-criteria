#!/usr/bin/env Rscript
## Plot one V-plot panel for Fig. S4-S6.
##
## Input format is the pipeline *_fragL_dist.txt table:
##   site_id    fragment_length    signed_distance
##
## Fig. S4: loMNase, Fig. S5: DNase, Fig. S6: ATAC.

suppressMessages({
  library(data.table)
  library(ggplot2)
})

if (!requireNamespace("ggpointdensity", quietly = TRUE)) {
  stop("Package 'ggpointdensity' is required. Install it with install.packages('ggpointdensity').")
}

ASSAY_PARS <- list(
  loMNase = list(xlim = c(-60, 60),   xbreaks = seq(-60, 60, 30),
                 ylim = c(0, 100),    ybreaks = seq(0, 100, 25),
                 width = 4.8, height = 4.3),
  DNase   = list(xlim = c(-60, 60),   xbreaks = seq(-60, 60, 30),
                 ylim = c(0, 100),    ybreaks = seq(0, 100, 25),
                 width = 4.8, height = 4.3),
  ATAC    = list(xlim = c(-100, 100), xbreaks = seq(-100, 100, 50),
                 ylim = c(0, 150),    ybreaks = seq(0, 150, 50),
                 width = 5.0, height = 4.1)
)

defaults <- list(
  input = NA_character_,
  output_prefix = NA_character_,
  title = NA_character_,
  assay = "loMNase",
  apex_x = NA_real_,
  apex_y = NA_real_,
  channel = FALSE,
  channel_mode = "full",
  axis_titles = FALSE,
  cap = 200000L,
  seed = 11L,
  point_size = 0.2,
  dpi = 200L,
  formats = "png,pdf",
  header = FALSE
)

usage <- function(status = 0) {
  cat(
    "Usage:\n",
    "  Rscript plot_single_vplot_panel.R --input FILE --assay loMNase --title CTCF --output-prefix OUT\n\n",
    "Required:\n",
    "  --input FILE             *_fragL_dist.txt table: site_id, fragment_length, signed_distance\n\n",
    "Common options:\n",
    "  --assay NAME             loMNase, DNase, ATAC, or custom [loMNase]\n",
    "  --title TEXT             panel title; default = input file stem\n",
    "  --output-prefix OUT      output path without extension; default = input file stem\n",
    "  --apex-x VALUE           fitted apex x-position for V-channel overlay\n",
    "  --apex-y VALUE           fitted V-channel width for V-channel overlay\n",
    "  --channel true|false     draw red dashed V-channel if apex is supplied [false]\n",
    "  --channel-mode full|inner full draws the four guide lines used for TF panels [full]\n",
    "  --axis-titles true|false add per-panel axis titles [false]\n",
    "  --cap N                  subsample rows after window filtering [200000]\n",
    "  --formats png,pdf        comma-separated output formats [png,pdf]\n",
    "  --dpi N                  PNG DPI [200]\n",
    "  --header true|false      set true only if the input table has column names [false]\n\n",
    "Custom axis options, used when --assay custom:\n",
    "  --xlim MIN,MAX --ylim MIN,MAX --xbreaks A,B,C --ybreaks A,B,C --width W --height H\n\n",
    "Examples:\n",
    "  Rscript plot_single_vplot_panel.R --input CTCF_loMNase_fragL_dist.txt --assay loMNase \\\n",
    "      --title CTCF --apex-x -1.69 --apex-y 31.13 --channel true --output-prefix S4_CTCF\n",
    "  Rscript plot_single_vplot_panel.R --input MNase_Bias_AAG_loMNase_fragL_dist.txt --assay loMNase \\\n",
    "      --title AAG --channel false --output-prefix S4_AAG\n",
    sep = ""
  )
  quit(save = "no", status = status)
}

parse_bool <- function(x) {
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

parse_num_vec <- function(x) {
  as.numeric(strsplit(as.character(x), ",", fixed = TRUE)[[1]])
}

parse_args <- function(argv) {
  opts <- defaults
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (key %in% c("-h", "--help")) usage(0)
    if (!startsWith(key, "--")) stop("Unknown argument: ", key)
    name <- sub("^--", "", key)
    name <- gsub("-", "_", name)
    if (i == length(argv)) stop("Missing value after ", key)
    value <- argv[[i + 1]]
    if (name %in% c("input", "output_prefix", "title", "assay", "channel_mode", "formats")) {
      opts[[name]] <- value
    } else if (name %in% c("apex_x", "apex_y", "point_size", "width", "height")) {
      opts[[name]] <- as.numeric(value)
    } else if (name %in% c("cap", "seed", "dpi")) {
      opts[[name]] <- as.integer(value)
    } else if (name %in% c("channel", "axis_titles", "header")) {
      opts[[name]] <- parse_bool(value)
    } else if (name %in% c("xlim", "ylim", "xbreaks", "ybreaks")) {
      opts[[name]] <- parse_num_vec(value)
    } else {
      stop("Unknown option: ", key)
    }
    i <- i + 2
  }
  opts
}

read_fragl_dist <- function(path, header = FALSE) {
  d <- fread(path, header = header)
  if (ncol(d) < 3) stop("Input must have at least 3 columns: site_id, fragment_length, signed_distance")
  if (header) {
    nms <- tolower(names(d))
    frag_col <- which(nms %in% c("fragment", "fragment_length", "fragl", "frag_len", "length"))[1]
    dist_col <- which(nms %in% c("distance", "signed_distance", "dist"))[1]
    if (is.na(frag_col) || is.na(dist_col)) {
      stop("Header input must contain fragment and distance columns.")
    }
    out <- data.table(fragment = as.numeric(d[[frag_col]]), distance = as.numeric(d[[dist_col]]))
  } else {
    out <- data.table(fragment = as.numeric(d[[2]]), distance = as.numeric(d[[3]]))
  }
  out[is.finite(fragment) & is.finite(distance)]
}

apex_segments <- function(ax, ay, ymax, ymin = 0, mode = "full") {
  if (!is.finite(ax) || !is.finite(ay)) return(NULL)
  if (mode == "inner") {
    return(data.frame(
      x = c(ax, ax),
      y = c(ay, ay),
      xend = c(ax - (ymax - ay) / 2, ax + (ymax - ay) / 2),
      yend = c(ymax, ymax)
    ))
  }
  data.frame(
    x    = c(ax,             ax,             (ax - ay) - (ymin - ay) / 2, (ax + ay) + (ymin - ay) / 2),
    y    = c(ay,             ay,             ymin,                         ymin),
    xend = c(ax - (ymax - ay) / 2, ax + (ymax - ay) / 2, (ax - ay) - (ymax - ay) / 2, (ax + ay) + (ymax - ay) / 2),
    yend = c(ymax,           ymax,           ymax,                         ymax)
  )
}

choose_pars <- function(opts) {
  if (opts$assay %in% names(ASSAY_PARS)) {
    pars <- ASSAY_PARS[[opts$assay]]
  } else if (opts$assay == "custom") {
    pars <- list()
  } else {
    stop("--assay must be loMNase, DNase, ATAC, or custom")
  }
  for (name in c("xlim", "ylim", "xbreaks", "ybreaks", "width", "height")) {
    if (!is.null(opts[[name]])) pars[[name]] <- opts[[name]]
  }
  required <- c("xlim", "ylim", "xbreaks", "ybreaks", "width", "height")
  missing <- required[!vapply(required, function(x) !is.null(pars[[x]]), logical(1))]
  if (length(missing) > 0) stop("Missing custom plotting parameter(s): ", paste(missing, collapse = ", "))
  pars
}

plot_one <- function(opts) {
  if (is.na(opts$input)) usage(1)
  input <- normalizePath(opts$input, mustWork = TRUE)
  if (is.na(opts$title)) opts$title <- tools::file_path_sans_ext(basename(input))
  if (is.na(opts$output_prefix)) {
    opts$output_prefix <- file.path(dirname(input), tools::file_path_sans_ext(basename(input)))
  }
  pars <- choose_pars(opts)
  set.seed(opts$seed)

  d <- read_fragl_dist(input, header = opts$header)
  original_rows <- nrow(d)
  d <- d[
    distance >= pars$xlim[1] & distance <= pars$xlim[2] &
      fragment >= pars$ylim[1] & fragment <= pars$ylim[2]
  ]
  filtered_rows <- nrow(d)
  if (nrow(d) > opts$cap) d <- d[sample(.N, opts$cap)]
  plotted_rows <- nrow(d)

  p <- ggplot(d, aes(distance, fragment)) +
    ggpointdensity::geom_pointdensity(size = opts$point_size, show.legend = FALSE) +
    scale_color_gradient2(low = "#88a7c5", high = "#104e8b") +
    scale_x_continuous(breaks = pars$xbreaks, expand = c(0, 0)) +
    scale_y_continuous(breaks = pars$ybreaks, expand = c(0, 0)) +
    theme_classic() +
    ggtitle(opts$title) +
    coord_fixed(ratio = 1, xlim = pars$xlim, ylim = pars$ylim, expand = FALSE) +
    theme(
      legend.position = "none",
      panel.border = element_rect(fill = NA, linewidth = 1, color = "black"),
      axis.text = element_text(size = 15, colour = "black"),
      plot.title = element_text(hjust = 0.5, size = 20, colour = "black"),
      axis.ticks = element_line(color = "black"),
      text = element_text(size = 15, family = "sans", colour = "black"),
      plot.margin = unit(c(0.3, 0.5, 0.2, 0.3), "cm")
    )

  if (opts$axis_titles) {
    p <- p + labs(x = "Distance from motif (bp)", y = "Fragment length (bp)")
  } else {
    p <- p + theme(axis.title = element_blank())
  }

  if (opts$channel) {
    segs <- apex_segments(opts$apex_x, opts$apex_y, pars$ylim[2], pars$ylim[1], opts$channel_mode)
    if (is.null(segs)) stop("--channel true requires finite --apex-x and --apex-y")
    p <- p + geom_segment(
      data = segs,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "red",
      linewidth = 0.6
    )
  }

  formats <- trimws(strsplit(opts$formats, ",", fixed = TRUE)[[1]])
  outputs <- character()
  for (fmt in formats) {
    out <- paste0(opts$output_prefix, ".", fmt)
    ggsave(out, plot = p, width = pars$width, height = pars$height, dpi = opts$dpi)
    outputs <- c(outputs, out)
  }

  cat("input:", input, "\n")
  cat("original rows:", original_rows, "\n")
  cat("rows in window:", filtered_rows, "\n")
  cat("plotted rows:", plotted_rows, "\n")
  cat("outputs:\n")
  cat(paste0("  ", outputs, collapse = "\n"), "\n")
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
plot_one(opts)
