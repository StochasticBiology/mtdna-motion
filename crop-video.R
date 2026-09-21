crop_video_following_traj_png <- function(
    input_video,
    output_video,
    traj,
    crop_width,
    crop_height,
    fps = 30
) {
  
  # ---- ensure ordering ----
  traj <- traj[order(traj$frame), ]
  
  n <- nrow(traj)
  
  # convert to top-left crop coordinates
  x0 <- round(traj$x - crop_width / 2)
  y0 <- round(traj$y - crop_height / 2)
  
  # ---- temp directories ----
  tmpdir <- tempfile("frame_extract_")
  dir.create(tmpdir)
  
  cropdir <- tempfile("frame_crop_")
  dir.create(cropdir)
  
  # ---- step 1: extract frames ----
  message("Extracting frames...")
  system(sprintf(
    "ffmpeg -y -i \"%s\" \"%s/frame_%%05d.png\"",
    input_video,
    tmpdir
  ))
  
  # ---- step 2: crop selected frames ----
  message("Cropping frames...")
  
  for (i in seq_len(n)) {
    
    frame_file <- sprintf("%s/frame_%05d.png", tmpdir, traj$frame[i] + 1)
    out_file   <- sprintf("%s/frame_%05d.png", cropdir, i)
    
    # skip if frame doesn't exist
    if (!file.exists(frame_file)) next
    
    system(sprintf(
      "ffmpeg -y -i \"%s\" -vf \"crop=%d:%d:%d:%d\" \"%s\"",
      frame_file,
      crop_width,
      crop_height,
      x0[i],
      y0[i],
      out_file
    ))
  }
  
  # ---- step 3: reassemble video ----
  message("Reassembling video...")
  
  system(sprintf(
    "ffmpeg -y -framerate %d -i \"%s/frame_%%05d.png\" -c:v libx264 -pix_fmt yuv420p \"%s\"",
    fps,
    cropdir,
    output_video
  ))
  
  message("Done.")
}

traj = list()

traj[[1]] <- data.frame(
  file = "AnalysedLiveSYBRVideos_seedling1-45minsybr001-cropped-cell1-igj.avi",
  frame = 0:47,
  x = rep(265, 48),
  y = rep(388, 48)
)
traj[[2]] <- data.frame(
  file = "AnalysedLiveSYBRVideos_seedling1-45minsybr001-cropped-cell1-igj.avi",
  frame = 44:70,
  x = rep(358, 27),
  y = rep(350, 27)
)

for(i in 1:length(traj)) {
  fileout = sprintf("output-%i.mp4", i)
  

crop_video_following_traj_png(
  input_video = traj[[i]]$file,
  output_video = fileout,
  traj = traj[[i]],
  crop_width = 30,
  crop_height = 30,
  fps = 6
)
}
