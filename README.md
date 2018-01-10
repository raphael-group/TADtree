# TADtree

TADtree is an algorithm the identification of hierarchical topological domains in Hi-C data.

It was developed by the [Raphael research group](http://compbio.cs.brown.edu) 
in the [Department of Computer Science](http://cs.brown.edu) and [Center for Computational Molecular Biology](http://brown.edu/ccmb) 
at [Brown University](http://brown.edu).

## Usage ##
To run, use command:
python <path to TADtree.py> <path to control file>

See below for guidance on choice of parameters. Note that the parameter "N" specifies the
maximum number of TADs to compute for each chromosome, but TADtree also computes TAD sets 
for all values of N less than the one given. The resulting TAD sets are outputted to:

<output_directory>/<contact_map_name>/N#.txt

where # is the number of TADs in the given TAD set, before removal of duplicates. The 
proportion of TADs in each TAD set removed as duplicates is recorded in 

<output_directory>/<contact_map_name>/proportion_duplicates.txt

## Parameters ##
S, the maximum allowable size for a TAD, is measured in bins. Thus, S must be increased to 
for higher-resolution data to detect TADs of the same physical size. Note, however, that
the run-time of TADtree is ~O(S^5) so high values of S will have a high performance
penalty.

Gamma determines the trade-off between sensitivity to deviations from background the 
contact map, and specificity for those deviations that are delineated by clear boundaries.
The default value of gamma is likely to work for a variety of data.  However, if TAD 
boundaries are poorly localized, it might be good to increase gamma. On the other hand, 
if multiple TADs accumulate at strong boundary points, gamma is probably too high.

M = 10 should always work.

p and q reflect the minimum scale over which changes in interaction preference can be 
robustly detected. My intuition is that these can stay the same, since this scale would 
naturally decrease for high resolution data, and p and q are given in units of bins.

As described above, a maximum value of N (the total number of TADs allowed on a given 
chromosome) is specified in the control file, but all smaller values of N are actually 
computed. This means that one can do manual model selection (i.e. choosing the final value 
of N) without re-running the algorithm. One principled way to choose the value of N is to 
look at the proportion of duplicate TADs, which is also outputted (see above). For example, 
it might be appropriate to choose the value of N for which 1-2% of all outputted TADs are
duplicates, but this is just a rule of thumb. Also note that duplicate TAD calls are 
filtered from the outputted TAD sets, as described in the paper.

## Reference
Caleb Weinreb and Benjamin J. Raphael.[Identification of hierarchical chromatin domains.](http://bioinformatics.oxfordjournals.org/content/early/2015/09/23/bioinformatics.btv485.abstract)
	Bioinformatics, (2015) pii: btv485.

[Website](http://compbio.cs.brown.edu/projects/tadtree)

    
