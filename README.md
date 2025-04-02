# Targeted Pangenome Graphs

``` sh
git clone --recursive https://github.com/ldenti/tpg
cd tpg
make -j4

# build sampled haplotype orders
./tpg build example/chrM.d9.gbz
# get subgraphs
./tpg extract example/chrM.d9.gbz example/regions.bed > regions.gfa
```

#### TODO
- [ ] (re)implement gfa-based extraction
