#######################################################
#           mandatory variables
#######################################################

# SAMPLE should be a data frame, row names must be sample id and there must be a
#   'subgroup' column which corresponds to classes of samples. There can also
#   be other annotation columns.
SAMPLE = data.frame(subgroup = ..., ...)
rownames(SAMPLE) = ...

# COLOR: a list of color settings corresponding to annotation column in 
#   'SAMPLE'. The value of each list element must be either a named vector
#   or a color mapping function. 'COLOR$subgroup' must be defined or random color
#   will be assigned. Names of other color settings should be same as
#   corresponding columns in 'SAMPLE'.
COLOR = list(subgroup = c(...), ...)

# CHROMOSOME: a vector of chromosome names.
CHROMOSOME = paste0("chr", 1:22)

# GENOME: abbreviation of species.
GENOME = "hg19"

# PROJECT_DIR: path of output directory. Several sub directories will be created.
PROJECT_DIR = ...

# GENOMIC_FEATURE_LIST: a list of genomic features as GRanges objects. There
#   must be a element named 'cgi' and 'cgi_shore' will be added accordingly. 
#   If `TXDB` is provided, gene, exon, intron, tss, intergenic
#   will be added automatically.
GENOMIC_FEATURE_LIST = list(
    cgi = ...
)

# how to get the object which stores data by chromosome
methylation_hooks$get_by_chr = function(chr) {
	...
    obj = list(gr = ...,
               raw = ...,
               cov = ...,
               meth = ...)
	return(obj)
})


#######################################################
#           optional variables
#######################################################


# TXDB (optional): a `GenomicFeatures::TxDb` object.
TXDB = loadDb(...)

# GTF file which is used to build 'TXDB'. If it is null, `metadata(TXDB)[3, "value"]` will be used
GTF_FILE = ...

# EXPR (optional): a matrix which contains expression values. Column names 
#   should be sample id and row names should be gene ids. Note gene ids in the 
#   expression matrix should be same type as genes in `GenomicFeatures::TxDb`.
EXPR = ...

# MARKS (optional): a vector of ChIP-Seq markers.
MAKRS = c(...)

# chipseq_hooks() is optional unless you want to do integrative analysis.
chipseq_hooks$sample_id = function(mark) {

})

chipseq_hooks$peak = function(mark, sid, ...) {

})

chipseq_hooks$chromHMM = function(sid, ...) {

})

CGI_SHORE_EXTEND = 2000
