set.seed(42)

library(ggplot2)
library(ggrepel)
library(latex2exp)
library(patchwork)

theme_set(theme_bw())
theme_update(
  axis.title.y = element_text(size = 12, angle = 0, vjust = 0.5, margin = margin(t = 0, r = 7.5, b = 0, l = 0)),
  axis.title.x = element_text(size = 8))

width <- 0.3528 * 410
height <- 0.3528 * 630
width.prsnt <- 0.3528 * 307
height.prsnt <- 0.3528 * 0.8 * 269

# https://stackoverflow.com/questions/36801488/how-to-distinguish-jumps-in-a-ggplot-line-plot
geom_conline <- function(
                         mapping = NULL,
                         data = NULL,
                         stat = "identity",
                         position = "identity",
                         ...,
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE
)
{
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomConLine,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}

GeomConLine <- ggproto(
  "GeomConLine",
  Geom,
  default_aes = aes(colour = "black", size = 0.5, linetype = 1, alpha = NA),
  draw_key = draw_key_vpath,
  # required_aes = c("x", "ymin", "ymax"),
  required_aes = c("x", "y", "con"),
  draw_panel = function(data, panel_scales, coord) {
    # data <- transform(data, xend = x, y = ymin, yend = ymax)
    data <- subset(
      transform(
        data,
        'xend' = c(x[2:nrow(data)], NA),
        'yend' = c(y[2:nrow(data)], NA),
        'con' = c(con[1:nrow(data) - 1], FALSE)
      ),
      'subset' = con,
      'select' = -con
    )
    ggplot2:::ggname("geom_stem", GeomSegment$draw_panel(data, panel_scales, coord))
  }
)

# simulation setup
t.start <- 0
t.end <- 8
t <- seq(t.start, t.end, length.out = 60 * 8) # sample every minute

arrivals <- function(lambda) {
  intervals <- c()
  done <- FALSE
  while (!done) {
    interval <- rexp(1, lambda)
    intervals <- c(intervals, interval)
    done <- sum(intervals) >= t.end
  }
  arvls <- cumsum(intervals)
  counts <- sapply(t, function(time) {
    sum(arvls <= time)
  })
  out <- list(arrivals = arvls, counts = counts)
  return(out)
}

serviceTimes <- function(n, alpha, beta) rgamma(n, alpha, beta)

workload <- function(t, w0, input) {
  W_t <- rep(0, length(t))
  idle.periods <- data.frame(begin = numeric(), end = numeric())

  idle.started <- FALSE
  for (i in seq_along(t)) {
    W_t[i] <- max(max(input[1:i]), w0) - input[i]
  }
  for (i in seq_along(t)) {
    if (W_t[i] == 0) {
      if (idle.started) {
        if (i < length(t) && W_t[i + 1] > 0) {
          idle.end <- t[i]
          idle.started <- FALSE
          idle.periods <- rbind.data.frame(idle.periods, data.frame(begin = idle.begin, end = idle.end))
        }
      } else {
        idle.begin <- t[i]
        idle.started <- TRUE
      }
    }
  }
  out <- list(workload = W_t, idle.periods = idle.periods)
  return(out)
}

# illustration
lambda <- 10
mu <- 1 / 10
sigma <- 2.5 / 60
alpha <- mu**2 / sigma**2
beta <- mu / sigma**2

res <- arrivals(lambda)
arvls <- res$arrivals
counts <- res$counts
service.times <- serviceTimes(length(arvls), alpha, beta)

X_t <- rep(0, length(t))
X_t <- sapply(seq_along(t), function(i) { t[i] - ifelse(counts[i] > 0, sum(service.times[1:counts[i]]), 0) })

w0 <- 1 / 6 # 10 minutes backlog
res <- workload(t, w0, X_t)
W_t <- res$workload
idle.periods <- res$idle.periods

connections <- sapply(arvls, function(arrival) { which.min(abs(arrival - t)) })
connections <- ifelse(seq_along(t) %in% connections, FALSE, TRUE)

p <- ggplot(data.frame(t = t, X_t = X_t)) +
  geom_conline(aes(t, X_t, "con" = connections)) +
  geom_conline(aes(t, X_t, "con" = !connections), linetype = 2) +
  geom_point(aes(x = arrivals, y = -Inf), data.frame(arrivals = arvls), size = 2) +
  coord_cartesian(clip = "off") +
  xlab("Time") +
  ylab(TeX("$X_t(\\omega)$", italic = TRUE))
# project
ggsave(file.path("plots", "input.eps"), p, units = "mm", width = height / 2.5, height = width / 1.25, device = cairo_ps, symbolfamily = cairoSymbolFont("OpenSymbol"))
# presentation
ggsave(file.path("plots", "input_prsnt.eps"), p, units = "mm", width = width.prsnt, height = height.prsnt, device = cairo_ps, symbolfamily = cairoSymbolFont("OpenSymbol"))

p <- ggplot(data.frame(t = t, W_t = W_t)) +
  geom_conline(aes(t, W_t, "con" = connections)) +
  geom_conline(aes(t, W_t, "con" = !connections), linetype = 2) +
  geom_point(aes(x = arrivals, y = -Inf), data = data.frame(arrivals = arvls), size = 2) +
  coord_cartesian(clip = "off") +
  xlab("Time") +
  ylab(TeX("$W_t(\\omega)$", italic = TRUE))
# project
ggsave(file.path("plots", "workload.eps"), p, units = "mm", width = height / 2.5, height = width / 1.25, device = cairo_ps, symbolfamily = cairoSymbolFont("OpenSymbol"))
# presentation
ggsave(file.path("plots", "workload_prsnt.eps"), p, units = "mm", width = width.prsnt, height = height.prsnt, device = cairo_ps, symbolfamily = cairoSymbolFont("OpenSymbol"))

# lambda * mu <= 1
lambda <- 12
mu <- 4 / 60
sigma <- 2.5 / 60
alpha <- mu**2 / sigma**2
beta <- mu / sigma**2

res <- arrivals(lambda)
arvls <- res$arrivals
counts <- res$counts
service.times <- serviceTimes(length(arvls), alpha, beta)

X_t <- rep(0, length(t))
X_t <- sapply(seq_along(t), function(i) { t[i] - ifelse(counts[i] > 0, sum(service.times[1:counts[i]]), 0) })

w0 <- 1 / 6 # 10 minutes backlog
res <- workload(t, w0, X_t)
W_t <- res$workload
idle.periods <- res$idle.periods

connections <- sapply(arvls, function(arrival) { which.min(abs(arrival - t)) })
connections <- ifelse(seq_along(t) %in% connections, FALSE, TRUE)

p <- ggplot(data.frame(t = t, X_t = X_t)) +
  geom_conline(aes(t, X_t, "con" = connections)) +
  geom_conline(aes(t, X_t, "con" = !connections), linetype = 2) +
  geom_point(aes(x = arrivals, y = -Inf), data.frame(arrivals = arvls), size = 2) +
  coord_cartesian(clip = "off") +
  xlab("Time") +
  ylab(TeX("$X_t(\\omega)$", italic = TRUE))
# project
ggsave(file.path("plots", "input_leq_1.eps"), p, units = "mm", width = height / 2.5, height = width / 1.25, device = cairo_ps, symbolfamily = cairoSymbolFont("OpenSymbol"))
# presentation
ggsave(file.path("plots", "input_leq_1_prsnt.eps"), p, units = "mm", width = width.prsnt, height = height.prsnt, device = cairo_ps, symbolfamily = cairoSymbolFont("OpenSymbol"))

p <- ggplot(data.frame(t = t, W_t = W_t)) +
  geom_conline(aes(t, W_t, "con" = connections)) +
  geom_conline(aes(t, W_t, "con" = !connections), linetype = 2) +
  geom_point(aes(x = arrivals, y = -Inf), data = data.frame(arrivals = arvls), size = 2) +
  geom_rect(aes(xmin = begin, xmax = end, ymin = -Inf, ymax = Inf), idle.periods, fill = "red", alpha = 0.3) +
  coord_cartesian(clip = "off") +
  xlab("Time") +
  ylab(TeX("$W_t(\\omega)$", italic = TRUE))
# project
ggsave(file.path("plots", "workload_leq_1.eps"), p, units = "mm", width = height / 2.5, height = width / 1.25, device = cairo_ps, symbolfamily = cairoSymbolFont("OpenSymbol"))
# presentation
ggsave(file.path("plots", "workload_leq_1_prsnt.eps"), p, units = "mm", width = width.prsnt, height = height.prsnt, device = cairo_ps, symbolfamily = cairoSymbolFont("OpenSymbol"))

# lambda * mu > 1
lambda <- 12
mu <- 1 / 10
sigma <- 2.5 / 60
alpha <- mu**2 / sigma**2
beta <- mu / sigma**2

res <- arrivals(lambda)
arvls <- res$arrivals
counts <- res$counts
service.times <- serviceTimes(length(arvls), alpha, beta)

X_t <- rep(0, length(t))
X_t <- sapply(seq_along(t), function(i) { t[i] - ifelse(counts[i] > 0, sum(service.times[1:counts[i]]), 0) })

w0 <- 1 / 6 # 10 minutes backlog
res <- workload(t, w0, X_t)
W_t <- res$workload
idle.periods <- res$idle.periods

connections <- sapply(arvls, function(arrival) { which.min(abs(arrival - t)) })
connections <- ifelse(seq_along(t) %in% connections, FALSE, TRUE)

p <- ggplot(data.frame(t = t, X_t = X_t)) +
  geom_conline(aes(t, X_t, "con" = connections)) +
  geom_conline(aes(t, X_t, "con" = !connections), linetype = 2) +
  geom_point(aes(x = arrivals, y = -Inf), data.frame(arrivals = arvls), size = 2) +
  coord_cartesian(clip = "off") +
  xlab("Time") +
  ylab(TeX("$X_t(\\omega)$", italic = TRUE))
# project
ggsave(file.path("plots", "input_ge_1.eps"), p, units = "mm", width = height / 2.5, height = width / 1.25, device = cairo_ps, symbolfamily = cairoSymbolFont("OpenSymbol"))
# presentation
ggsave(file.path("plots", "input_ge_1_prsnt.eps"), p, units = "mm", width = width.prsnt, height = height.prsnt, device = cairo_ps, symbolfamily = cairoSymbolFont("OpenSymbol"))

p <- ggplot(data.frame(t = t, W_t = W_t)) +
  geom_conline(aes(t, W_t, "con" = connections)) +
  geom_conline(aes(t, W_t, "con" = !connections), linetype = 2) +
  geom_point(aes(x = arrivals, y = -Inf), data = data.frame(arrivals = arvls), size = 2) +
  geom_rect(aes(xmin = begin, xmax = end, ymin = -Inf, ymax = Inf), idle.periods, fill = "red", alpha = 0.3) +
  coord_cartesian(clip = "off") +
  xlab("Time") +
  ylab(TeX("$W_t(\\omega)$", italic = TRUE))
# project
ggsave(file.path("plots", "workload_ge_1.eps"), p, units = "mm", width = height / 2.5, height = width / 1.25, device = cairo_ps, symbolfamily = cairoSymbolFont("OpenSymbol"))
# presentation
ggsave(file.path("plots", "workload_ge_1_prsnt.eps"), p, units = "mm", width = width.prsnt, height = height.prsnt, device = cairo_ps, symbolfamily = cairoSymbolFont("OpenSymbol"))

