library(R.matlab)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ragg)
library(akima)


data <- readMat("mainfoldheatmap.mat")

xypoint1 <- data.frame(x = data$x1, y = data$y1)
xypoint2 <- data.frame(x = data$x2, y = data$y2)


df_hex <- data.frame(
  x = rep(data$xx, times = data$ndensity),
  y = rep(data$yy, times = data$ndensity)
)

df_density <- df_hex %>%
  count(x, y, name = "z")


interp_data <- with(df_density, interp(x, y, z = z, duplicate = "mean", nx = 1000, ny = 1000))
df_interp <- expand.grid(x = interp_data$x, y = interp_data$y)
df_interp$z <- as.vector(interp_data$z)
df_interp$z[is.na(df_interp$z)] <- 0 

df_interp$z <- floor(pmax(0, df_interp$z))

df_hex <- data.frame(
  x = rep(df_interp$x, times = df_interp$z),
  y = rep(df_interp$y, times = df_interp$z)
)

ggplot(df_hex, aes(x = x, y = y)) +
  geom_hex() +
  scale_fill_gradientn(
    colors = c("white", "red", "yellow"),
    values = scales::rescale(c(0, 0.06, 1)),
    na.value = "white"
  )+ 
  geom_point(data = xypoint1, aes(x = x, y = y),
             shape = 16, fill = "blue", color = "blue", size = 6, alpha = 0.5) +
  geom_point(data = xypoint2, aes(x = x, y = y),
             shape = 16, fill = "blue", color = "blue", size = 6, alpha = 0.5) +
  coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) -> p

p_temp_file <- tempfile(tmpdir = "./tmp", fileext = '.png')
agg_png(p_temp_file, width = 5.33, height = 5.33, units = "in", res = 300)
suppressWarnings(print(p))
invisible(dev.off())

