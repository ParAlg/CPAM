#!/bin/bash
# Download and build graph inputs

db_source="https://snap.stanford.edu/data/bigdata/communities/com-dblp.ungraph.txt.gz"
yt_source="https://snap.stanford.edu/data/bigdata/communities/com-youtube.ungraph.txt.gz"
lj_source="https://snap.stanford.edu/data/soc-LiveJournal1.txt.gz"
co_source="https://snap.stanford.edu/data/bigdata/communities/com-orkut.ungraph.txt.gz"

# fs_source="https://snap.stanford.edu/data/bigdata/communities/com-friendster.ungraph.txt.gz"
# ru_source="http://nrvis.com/download/data/misc/road_usa.zip"

sources=($db_source $yt_source $lj_source $co_source)

cd ../other_systems/gbbs/utils
make snap_converter
cd ../../../inputs

snap_converter="../other_systems/gbbs/utils/snap_converter"

for f in ${sources[@]}; do
  wget $f
  compressed_fname="${f##*/}"
  gunzip $compressed_fname
  txt_filename=${compressed_fname%.gz}
  filename=${txt_filename%.txt}
  $snap_converter -s -i $txt_filename -o $filename.adj
  rm $txt_filename
done


