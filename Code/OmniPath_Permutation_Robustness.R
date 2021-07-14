# 0. Overview:
{
# The Idea with this script is to:
# .	run all methods on a single resource (OP)
# .	take the topmost ranked interactions for each method -> R-zero
# .	create a modified OP resource in which the topmost ranked interactions
#     remain, but of the remainder x% of the interactions have been removed and 
#     replaced with entirely random pairs of genes derived from the test data 
#     that do not exist in the resource, x =10,20,40% etc.
# .	Rerun methods on modified omnipath resource, get top ranks -> R-modified
# .	plot percentage of R-zero in R-modified over x and investigate result
}

# 1. 