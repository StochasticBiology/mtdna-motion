library(ggplot2)
library(ggpubr)
library(dplyr)
library(spatstat.geom)
library(spatstat.explore)
library(viridis)

# pull in nucleoid spot tracking data from FIJI
# [
# these are found using TrackMate on the nucleoid channel of the Chustecki et al. datasets
# with a 0.4um diameter, 2-3-2 settings for traj connections, and variable quality thresholds (see methods spreadsheet)
# ]
# extract a dataframe of mito-time-mindist-ID of closest partner
# spatial unit is microns
df = result = data.frame()
for (f in list.files("WheatTracks", pattern = "dnas\\.csv$", full.names = TRUE)) {
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

### speed distribution
dna.speed.df <- df |>
  dplyr::group_by(id, ref) |>
  dplyr::filter(dplyr::n() > 2) |>
  dplyr::slice(2:(dplyr::n() - 1)) |>
  dplyr::mutate(
    speed = sqrt((x - dplyr::lag(x))^2 +
                   (y - dplyr::lag(y))^2) /
      (t - dplyr::lag(t))
  ) |>
  dplyr::ungroup() |>
  dplyr::select(id, t, speed, ref)

ggplot(dna.speed.df, aes(x = log10(speed), fill = ref)) + 
  geom_histogram() 
# weird y axis?

### counts
dna.count.df <- df |>
  dplyr::group_by(ref, t) |>
  dplyr::summarise(n_ids = dplyr::n_distinct(id), .groups = "drop")

ggplot(dna.count.df, aes(x = n_ids, fill=ref)) + geom_histogram()

### areas
dna.cell.area.df <- df |>
  dplyr::filter(t == 0) |>
  dplyr::group_by(ref) |>
  dplyr::summarise(
    area = {
      h <- grDevices::chull(x, y)
      xh <- x[h]
      yh <- y[h]
      abs(sum(xh * dplyr::lead(yh, default = yh[1]) -
                yh * dplyr::lead(xh, default = xh[1]))) / 2
    },
    .groups = "drop"
  )

#### densities

dna.density.df <- dna.count.df |>
  dplyr::left_join(dna.cell.area.df, by = "ref") |>
  dplyr::mutate(density = n_ids / area)

ggplot(dna.density.df, aes(x = density, fill=ref)) + geom_histogram()

#### end physical stats


refs = unique(df$ref)

# plot trajectories by sample to check
ggplot(df, aes(x=x, y=y, color=factor(id))) + 
  geom_path() + theme(legend.position = "none") +
  facet_wrap(~ ref)

# distribution of minimum distances across frames
# peak around 2-3 microns
ggplot(result, aes(x = mindist)) + scale_x_log10() + geom_histogram() + facet_wrap(~ref)
sep.plot = ggplot(result, aes(x = mindist, fill = ref)) +
  scale_fill_viridis_d() +
  scale_x_log10() + geom_histogram() +
  theme_minimal() + theme(legend.position = "none")
sep.plot

# structured "partnering" plot
ggplot(result, aes(x = id, y = partner)) + geom_point(alpha = 0.05)

# attempt at trajectories
ggplot(result[result$id == 1,], aes(x = t, y = mindist, color = factor(partner))) + 
  geom_line() + scale_y_log10() + theme(legend.position = "none")

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
  geom_path(data=tracts.xy, aes(x=x2,y=y2,color=factor(id)), alpha = 0.5) +
  theme_minimal() + theme(legend.position = "none")

# what proportion of mito-frames are in the persistent pair set?
pp.prop = c(pp.prop, nrow(tracts.xy)/nrow(df[df$ref==refs[expt.ref],]))
}

ggarrange(plotlist = pp.plot)
pp.prop
c(length(pp.prop), mean(pp.prop), sd(pp.prop))

sf = 2
png("wheat-pp-example.png", width=200*sf, height=500*sf, res=72*sf)
print(pp.plot[[1]])
dev.off()

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

ggarrange(dna.vis.list[[1]], dna.plot.list[[1]])

#############


mito.df = data.frame()
for (f in list.files("WheatTracks", pattern = "mitos\\.csv$", full.names = TRUE)) {
  tdf <- read.csv(f)
  tdf = tdf[4:nrow(tdf),]
  tdf1 = data.frame(id = as.numeric(tdf$TRACK_ID), t = as.numeric(tdf$POSITION_T),
                    x = as.numeric(tdf$POSITION_X), y = as.numeric(tdf$POSITION_Y))
  tdf1 <- tdf1[order(tdf1$id, tdf1$t), ]
  tdf1$ref = f
  mito.df = rbind(mito.df, tdf1)
}

### speed distribution
mito.speed.df <- mito.df |>
  dplyr::group_by(id, ref) |>
  dplyr::filter(dplyr::n() > 2) |>
  dplyr::slice(2:(dplyr::n() - 1)) |>
  dplyr::mutate(
    speed = sqrt((x - dplyr::lag(x))^2 +
                   (y - dplyr::lag(y))^2) /
      (t - dplyr::lag(t))
  ) |>
  dplyr::ungroup() |>
  dplyr::select(id, t, speed, ref)

ggplot(mito.speed.df, aes(x = log10(speed), fill = ref)) + 
  geom_histogram() 
# weird y axis?

### counts
mito.count.df <- mito.df |>
  dplyr::group_by(ref, t) |>
  dplyr::summarise(n_ids = dplyr::n_distinct(id), .groups = "drop")

ggplot(mito.count.df, aes(x = n_ids, fill=ref)) + geom_histogram()

### areas
mito.cell.area.df <- mito.df |>
  dplyr::filter(t == 0) |>
  dplyr::group_by(ref) |>
  dplyr::summarise(
    area = {
      h <- grDevices::chull(x, y)
      xh <- x[h]
      yh <- y[h]
      abs(sum(xh * dplyr::lead(yh, default = yh[1]) -
                yh * dplyr::lead(xh, default = xh[1]))) / 2
    },
    .groups = "drop"
  )

#### densities

mito.density.df <- mito.count.df |>
  dplyr::left_join(mito.cell.area.df, by = "ref") |>
  dplyr::mutate(density = n_ids / area)

ggplot(mito.density.df, aes(x = density, fill=ref)) + geom_histogram()

#### end physical stats

p.style = theme_minimal() + 
  theme(legend.position = "none") 
phys.plot = ggarrange(
  ggplot(dna.speed.df, aes(x=log10(speed), fill=ref)) + 
    geom_histogram() + scale_fill_viridis_d() + p.style,# + xlim(-4,0),
  ggplot(dna.count.df, aes(x=n_ids, fill=ref)) + 
    geom_histogram() + scale_fill_viridis_d() + p.style,# + xlim(0,120),
  ggplot(dna.density.df, aes(x=density, fill=ref)) + 
    geom_histogram() + scale_fill_viridis_d() + p.style,# + xlim(0.05,0.3),
  ggplot(mito.speed.df, aes(x=log10(speed), fill=ref)) + 
    geom_histogram() + scale_fill_viridis_d() + p.style,# + xlim(-4,0),
  ggplot(mito.count.df, aes(x=n_ids, fill=ref)) + 
    geom_histogram() + scale_fill_viridis_d() + p.style,# + xlim(0,120),
  ggplot(mito.density.df, aes(x=density, fill=ref)) + 
    geom_histogram() + scale_fill_viridis_d() + p.style,# + xlim(0.05,0.3),
  labels = c("Ai", "ii", "iii", "Bi", "ii", "iii")
)
phys.plot.wheat = ggarrange(
  ggplot(dna.speed.df, aes(x=log10(speed), fill=ref)) + 
    geom_histogram() + scale_fill_viridis_d() + p.style,# + xlim(-4,0),
  ggplot(dna.density.df, aes(x=density, fill=ref)) + 
    geom_histogram() + scale_fill_viridis_d() + p.style,# + xlim(0.05,0.3),
  ggplot(mito.speed.df, aes(x=log10(speed), fill=ref)) + 
    geom_histogram() + scale_fill_viridis_d() + p.style,# + xlim(-4,0),
  ggplot(mito.density.df, aes(x=density, fill=ref)) + 
    geom_histogram() + scale_fill_viridis_d() + p.style,# + xlim(0.05,0.3),
  labels = c("Ci", "ii", "Di", "ii")
)

phys.plot.wheat

####

png("wheat-collective-2.png", width=500*sf, height=300*sf, res=72*sf)
ggarrange(phys.plot.wheat, pp.plot[[2]],
         
          nrow=2)
dev.off()

png("wheat-collective-3.png", width=500*sf, height=200*sf, res=72*sf)
ggarrange(phys.plot.wheat)
dev.off()


