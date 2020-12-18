# concorR (development version)

* Added a check for vertex names to `concor_igraph_apply()`, add them if 
missing.
* Fixed `make_blk()` and `concor_make_igraph()` to carry vertex names through if 
they had to be added.
* Added a `stat` option to `make_reduced()`, so it can use density to threshold 
reduced network links (the original behavior) or average degree.
* Updated `plot_reduced()` color palette so that it auto-generates the number of needed colors. Using viridis now instead of rainbow.
* Added a function `make_reduced_from_partition()` which makes a reduced network matrix based on a given partition.

# concorR 0.2.1

* `concor_make_igraph()` and `concor_igraph_apply()` fixed to use edge weights 
if present in the input adjacency matrix (`concor_make_igraph()`) or igraph 
object (`concor_igraph_apply()`) (#6).
