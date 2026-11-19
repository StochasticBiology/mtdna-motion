# mtdna-motion

Data and analysis for the manuscript "Physical constraints on the social-genetic dynamics of mtDNA inside plant mitochondria".

<img width="813" height="130" alt="image" src="https://github.com/user-attachments/assets/d1128299-cb3c-42e8-a2cc-4c893a4bf89e" />

Contents
----
This repository makes use of the data from http://github.com/jo-c-bio/physical-genetic-dynamics .

* `fig-case-studies.csv` -- description of source files and coordinates in videos of the case study behaviours in the manuscript

Code
---
* `analyse-tracks[-wheat].R` -- R code to analyse the TrackMate trajectories from video data
* `jie-model.R` -- graphics for the mathematical model of nucleoid dynamics in the manuscript
* `mito-spcacing.R` -- XXX
* `mtdnas.R` -- XXX
* `crop-video.R` -- XXX

Data
---
* `Videos/` -- raw video data (`.oir` and `.czi` format; use FIJI) from wheat and Arabidopsis mitochondrial and mtDNA live imaging. Several are compressed due to Github file size limits; unzip before use.
* `Tracks2/` -- TrackMate trajectories (`.csv` format) for mitochondria and nucleoids from Arabidopsis video data
* `WheatTracks/` -- TrackMate trajectories (`.csv` format) for mitochondria and nucleoids from wheat video data
