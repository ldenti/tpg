# Targeted Pangenome Graphs

Utility to extract subgraphs (or subhaplotypes) from a pangenome graph (in GBZ format). Regions of interest are provided as BED file.

**Note 1:** BED entries must have names (4th column).

**Note 2:** to locate the regions on the graph, we will use the reference paths included in the graph. By default, we use "CHM13" to select the reference paths. If that's not the case, please use the `-r` option (e.g., `-r GRCh38`).

### HowTo

``` sh
git clone --recursive https://github.com/ldenti/tpg
cd tpg
make -j4

# build R-index
./tpg build example/chrM.d9.gbz

# get subgraphs (in GFA format)
./tpg extract example/chrM.d9.gbz example/regions.bed > regions.gfa
# get subhaplotypes (in FASTA format)
./tpg extract -f example/chrM.d9.gbz example/regions.bed > regions.fa
```