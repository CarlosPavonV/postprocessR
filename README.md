# postprocessR
Random R scripts to process the output of some software that is commonly used in evolutionary biology:<br/><br/>
biogeobears2simmap creates a simmap object from the output of stochastic character mapping with BioGeoBears (https://github.com/nmatzke/BioGeoBEARS).<br/><br/>
get_rrphylo_rtt estimates a summary statistic for the absolute rates through time from the RRphylo (https://github.com/cran/RRphylo) output (e.g., you can calculate the mean of the rates of trait evolution for n points in time). Right now it can only handle univariate data and extant-only trees.<br/><br/>
sampleMCMCTree generates a sample of tress based on the 95% CI of divergence times estimated by MCMCTree (http://abacus.gene.ucl.ac.uk/software/paml.html).<br/><br/>
sampleBPPtrees takes the FigTree-readable output of BPP analysis A00 (https://github.com/bpp/bpp) and generates a sample of trees (in a multiPhylo object) based on the 95% HPD of divergence times.<br/><br/>
