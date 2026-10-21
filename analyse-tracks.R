library(ggplot2)
library(ggpubr)
library(dplyr)
library(spatstat.geom)
library(spatstat.explore)

# pull in spot tracking data from FIJI
# extract a dataframe of mito-time-mindist-ID of closest partner
# spatial unit is microns
df = result = data.frame()
for (f in list.files(pattern = "\\.csv$")) {
  tdf <- read.csv(f)
  tdf = tdf[4:nrow(tdf),]
  tdf1 = data.frame(id = as.numeric(tdf$TRACK_ID), t = as.numeric(tdf$POSITION_T),
                  x = as.numeric(tdf$POSITION_X), y = as.numeric(tdf$POSITION_Y))
  tdf1 <- tdf1[order(tdf1$id, tdf1$t), ]
  tdf1$ref = f
  
  tresult <- tdf1 |>
    group_by(t) |>
    group_modify(\(d, ...) {
      distmat <- as.matrix(dist(d[, c("x", "y")]))
      diag(distmat) <- Inf
      j <- max.col(-distmat, ties.method = "first")
      tibble(
        id = d$id,
        mindist = distmat[cbind(seq_len(nrow(d)), j)],
        partner = d$id[j]
      )
    }) |>
    ungroup() |>
    select(id, t, mindist, partner)
  tresult$ref = f
  df = rbind(df, tdf1)
  result = rbind(result, tresult)
}

refs = unique(df$ref)

# plot trajectories by sample to check
ggplot(df, aes(x=x, y=y, color=factor(id))) + 
  geom_path() + theme(legend.position = "none") +
  facet_wrap(~ ref)

# distribution of minimum distances across frames
# peak around 2-3 microns
ggplot(result, aes(x = mindist)) + scale_x_log10() + geom_histogram() + facet_wrap(~ref)

# structured "partnering" plot
ggplot(result, aes(x = id, y = partner)) + geom_point(alpha = 0.05)

# attempt at trajectories
ggplot(result[result$id == 1,], aes(x = t, y = mindist, color = factor(partner))) + 
  geom_line() + scale_y_log10() + theme(legend.position = "none")

library(dplyr)

# loop through experiments
pp.plot = list()
pp.prop = c()
for(expt.ref in 1:length(refs)) {
  # pull examples where two mitos are within dist.thresh of each other for nrow.thresh consecutive rows
  
  nrow.thresh = 15
dist.thresh = 2

tracts <- result[result$ref == refs[expt.ref] & result$mindist < dist.thresh,] |>
  group_by(id) |>
  mutate(run = cumsum(partner != lag(partner, default = first(partner)))) |>
  group_by(id, run) |>
  filter(n() > nrow.thresh) |>
  ungroup()

# add the x,y coordinates of both partners involved in these pairs for plotting 
tracts.xy <- tracts[tracts$id > tracts$partner,] |>
  dplyr::left_join(
    df[df$ref==refs[expt.ref],] |> dplyr::rename(x1 = x, y1 = y),
    by = c("id", "t")
  ) |>
  dplyr::left_join(
    df[df$ref==refs[expt.ref],] |> dplyr::rename(id = id, x2 = x, y2 = y),
    by = c("partner" = "id", "t")
  )

# plot persistent pairs on the background of all trajectories
pp.plot[[expt.ref]] = ggplot() +
  geom_path(data=df[df$ref==refs[expt.ref],], aes(x=x, y=y, group=id), color="#EEEEEE") +
  geom_path(data=tracts.xy, aes(x=x1,y=y1,color=factor(id))) +
  geom_path(data=tracts.xy, aes(x=x2,y=y2,color=factor(id))) +
  theme_minimal() + theme(legend.position = "none")

# what proportion of mito-frames are in the persistent pair set?
pp.prop = c(pp.prop, nrow(tracts.xy)/nrow(df[df$ref==refs[expt.ref],]))
}

ggarrange(plotlist = pp.plot)
pp.prop
c(length(pp.prop), mean(pp.prop), sd(pp.prop))

sf = 2
png("pp-example.png", width=200*sf, height=500*sf, res=72*sf)
print(pp.plot[[1]])
dev.off()


ex.set = df[df$ref == refs[4] & df$id %in% c(tracts$id, tracts$partner),]
ggplot(ex.set, aes(x = x, y = y, color = factor(id))) + geom_path()

###### separation analysis


snap = 0
dna.vis.list = dna.plot.list = list()
for(i in 1:length(refs)) {
  subdf = df[df$ref == refs[i],]
  snap = snap+1
  sub.traj.t = subdf[subdf$t == max(subdf$t),]
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
  dat.df <- data.frame(
    r = env$r,
    obs = env$obs,
    mean = env$mmean,
    lo = env$lo,
    hi = env$hi)
  
  dna.plot.list[[snap]] = ggplot(dat.df, aes(x = r)) +
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

disp <- subdf %>%
  arrange(id, t) %>%
  group_by(id) %>%
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

