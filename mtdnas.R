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
