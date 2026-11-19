git clone https://github.com/jo-c-bio/physical-genetic-dynamics
cd Videos
for f in *.zip; do
    [ -f "$f" ] && unzip -o "$f"
done
cp ../physical-genetic-dynamics/AnalysedLiveAOX-mt-HI-NESSVideos/SourceCZIs/*.czi .
