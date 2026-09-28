library(spatstat.geom)
library(spatstat.explore)
library(ggplot2)
library(ggpubr)

dna.vis.list = list()
dna.plot.list = list()
snap = 1

# output manually from TrackMate
df = read.csv("AOX-sample9-240frames-3DDriftCorrect-cell1_spots.csv")
df = df[4:nrow(df),]
df$POSITION_X = as.numeric(df$POSITION_X)
df$POSITION_Y = as.numeric(df$POSITION_Y)
df$FRAME = as.numeric(df$FRAME)
df$TRACK_ID = as.numeric(df$TRACK_ID)
df <- df[order(df$TRACK_ID, df$FRAME), ]

ggplot(df, aes(x=POSITION_X, y=POSITION_Y, colour=factor(TRACK_ID))) + geom_path() + theme(legend.position = "none")

subdf = df[,c("TRACK_ID", "FRAME", "POSITION_X", "POSITION_Y")]
colnames(subdf) = c("traj", "t", "x", "y")

snap = 0
for(i in c(1, 100, 200)) {
  snap = snap+1
  sub.traj.t = subdf[subdf$t == i,]
  pts = sub.traj.t[,c("x", "y")]
  
  # convex hull from your points
  h <- chull(pts$x, pts$y)
  h <- rev(h)
  
  # get domain of points
  win <- owin(poly = list(x = pts$x[h], y = pts$y[h]))
  
  # point pattern within our cell
  pp <- ppp(x = pts$x, y = pts$y, window = win)
  
  win_df <- data.frame(
    x = pts$x[h],
    y = pts$y[h]
  )
  
  dna.vis.list[[snap]] = ggplot() + 
    geom_polygon(data = win_df, aes(x = x, y = y),
                 fill = NA, colour = "black", linewidth = 1) +
    geom_point(data = pts, aes(x = x, y = y),
               size = 1) +
    coord_equal() + theme_classic() + labs(x = "x", y = "y")
  
  L <- Lest(pp, correction = "translation")  # robust default
  plot(L)
  abline(h = 0, lty = 2)
  env <- envelope(
    pp,
    fun = Lest,
    nsim = 99,
    correction = "translation",
    simulate = expression(runifpoint(npoints(pp), win = win)) # null model in our cell
  )
  
  #      plot(env)
  #     abline(h = 0, lty = 2)
  df <- data.frame(
    r = env$r,
    obs = env$obs,
    mean = env$mmean,
    lo = env$lo,
    hi = env$hi)
  
  dna.plot.list[[snap]] = ggplot(df, aes(x = r)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") +
    geom_line(aes(y = obs), linewidth = 1) +
    geom_line(aes(y = mean), linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme_classic() +
    labs(y = "L(r)", x = "Distance (r)")
}

ggarrange(plotlist = c(dna.vis.list, dna.plot.list), nrow=2, ncol = length(dna.vis.list))

#####

library(dplyr)
library(tidyr)

# example inputs:
# df with columns: traj, t, x, y
df = subdf  
D = 10
d = 5

# 1. compute total displacement per trajectory
disp <- df %>%
  arrange(traj, t) %>%
  group_by(traj) %>%
  summarise(
    disp = sqrt((last(x) - first(x))^2 + (last(y) - first(y))^2),
    .groups = "drop"
  )

# keep only trajectories that move enough
active_traj <- disp %>%
  filter(disp >= D) %>%
  pull(traj)

df2 <- df %>%
  filter(traj %in% active_traj)

# 2. self-join by time to compare trajectories at same t
pairs <- df2 %>%
  inner_join(df2, by = "t", suffix = c("_1", "_2")) %>%
  filter(traj_1 < traj_2) %>%
  mutate(
    dist = sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2)
  )

# 3. find pairs that are ever within distance d
close_pairs <- pairs %>%
  group_by(traj_1, traj_2) %>%
  summarise(
#    ever_close = any(dist <= d),
    ever_close = mean(dist <= d) >= 0.5,
    min_dist = min(dist),
    .groups = "drop"
  ) %>%
  filter(ever_close)

close_pairs
pair.set = unlist(close_pairs[,1:2], use.names = FALSE)
ggplot(data=subdf,
       aes(x=x, y=y, color=factor(traj))) +
  geom_path(alpha=0.1) +
  geom_path(data=subdf[subdf$traj %in% pair.set,], 
            aes(x=x, y=y, color=factor(traj))) +
  theme_minimal() + 
  theme(legend.position = "none")

