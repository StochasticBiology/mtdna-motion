library(spatstat.geom)
library(spatstat.explore)
library(ggplot2)
library(ggpubr)

src.traj.files = paste("../Mitosim/optimal-mitos-main/plant-mito-dynamics/mtgfp-rawtrajectories/mtGFP-", 1:5, ".xml-rawtrajs.csv", sep="")

snap = 0
vis.list = plot.list = list()

for(expt in 1:length(src.traj.files)) {
  # initialise storage
  # extract data for this experiment
  sub.traj = read.csv(src.traj.files[expt])
  if(colnames(sub.traj)[1] == "frame") {
    colnames(sub.traj) = c("t", "traj", "active", "x", "y")
  }
  for(i in c(1)) { #}, max(sub.traj$t)/4, 2*max(sub.traj$t)/4, 3*max(sub.traj$t)/4, max(sub.traj$t))) {
    snap = snap + 1
    # pull trajectories up to this time point, and plot
    sub.traj.t = sub.traj[sub.traj$t == i,]
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
    
    vis.list[[snap]] = ggplot() + 
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
    
    plot.list[[snap]] = ggplot(df, aes(x = r)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") +
      geom_line(aes(y = obs), linewidth = 1) +
      geom_line(aes(y = mean), linetype = "dashed") +
      geom_hline(yintercept = 0, linetype = "dotted") +
      theme_classic() +
      labs(y = "L(r)", x = "Distance (r)")
  }
}
ggarrange(plotlist = c(vis.list, plot.list), nrow=2, ncol = length(vis.list))
