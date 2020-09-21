get.processors <- function(processors){
  
  if(Sys.info()["sysname"] == 'Windows'){
    # mclapply is not supported on windows
    # so we give a single processor,
    # in which case mclapply calls fall back
    # on lapply
    return(1)
  }
  
  if(is.null(processors)){ 
    processors = detectCores(all.tests = FALSE, logical = FALSE)
  }
  
  return(processors)
  
}
############################
#' Create detailed summaries of all ABI sequencing reads in a folder (recursively searched)
#' 
#'
#' @keywords chromatogram, peak, mismatch
#'
#' @export summarise.abi.folder



summarise.abi.folder <- function(input.folder, trim.cutoff = 0.0001, secondary.peak.ratio = 0.33, write.secondary.peak.files = FALSE, processors = NULL){
  
  processors = get.processors(processors)
  
  print("Looking for .ab1 files...")
  abi.fnames = list.files(input.folder, pattern = "\\.ab1$", full.names = T, recursive = T)
  
  print(sprintf(("Found %d .ab1 files..."), length(abi.fnames)))
  
  
  print("Loading reads...")
  abi.seqs = mclapply(abi.fnames, read.abif, mc.cores = processors)
  
  print("Calculating read summaries...")
  # now make a data.frame of summaries of all the files
  summaries.dat = mclapply(abi.seqs, 
                           summarise.abi.file,
                           trim.cutoff = trim.cutoff,
                           secondary.peak.ratio = secondary.peak.ratio,
                           processors = 1,
                           mc.cores = processors  
  )
  
  print("Cleaning up")
  summaries = mclapply(summaries.dat, function(x) x[["summary"]], mc.cores = processors)
  summaries = do.call(rbind, summaries)
  
  reads = mclapply(summaries.dat, function(x) x[["read"]], mc.cores = processors)
  
  folder.names = basename(dirname(abi.fnames))
  file.names = basename(abi.fnames)
  
  summaries = cbind.data.frame("file.path" = as.character(abi.fnames), "folder.name" = as.character(folder.names), "file.name" = file.names, summaries, stringsAsFactors = FALSE)
  
  return(list("summaries" = summaries, "reads" = reads))
  
}
################
#' Trim an alignment to the limits of a supplied reference sequence
#'
#' This function profile aligns the alignment to the
#' reference sequence. It then trims the alignment to the limits of the 
#' reference sequence, and returns the trimmed alignment.
#'
#' @param alignment a DNAStringSet object
#' @param reference a DNA sequence as a DNAString sequence
#'
#' @export trim.to.reference

trim.to.reference <- function(alignment, reference){
  
  if(class(alignment)!='DNAStringSet'){ stop("alignment must be a DNAStringSet object")}
  if(class(reference)!='DNAString'){ stop("reference must be a DNAString object")}
  
  ref = DNAStringSet(reference)
  names(ref) = 'ref'
  aln = AlignProfiles(alignment, ref)
  
  ref.aligned = as.character(aln['ref'])
  not.gaps = str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
  ref.start = min(not.gaps)
  ref.finish = max(not.gaps)
  
  aln.trimmed = subseq(alignment, ref.start, ref.finish)
  
  return(aln.trimmed)    
  
}
#################
#' Produce a one-row summary of a merged read
#' 
#' @param merged.read a merged.read object produced by the merge.reads() function
#'
#' @return A numeric vector including (note that insertions, deletions, and stop codons will only be called if the read was merged with a reference amino acid sequence):
#'          \enumerate{
#'              \item {consensus.len}: the length of the consensus sequence\cr
#'              \item {n.reads}: the number of reads used to build the consensus\cr
#'              \item {distance.min}: the minimum distance between two reads \cr
#'              \item {distance.max}: the maximum distance between two reads \cr
#'              \item {distance.med}: the median distance between pairs of reads \cr
#'              \item {readlen.min}: the length of the smallest read \cr
#'              \item {readlen.max}: the length of the longest read \cr
#'              \item {readlen.med}: the median of all read lengths \cr
#'              \item {insertions.min}: the smallest number of insertions inferred in a read \cr
#'              \item {insertions.max}: the larges number of insertions inferred in a read \cr
#'              \item {insertions.med}: the median number of insertions inferred in a read \cr
#'              \item {deletions.min}: the smallest number of deletions inferred in a read \cr
#'              \item {deletions.max}: the larges number of deletions inferred in a read \cr
#'              \item {deletions.med}: the median number of deletions inferred in a read \cr
#'              \item {stops.min}: the smallest number of stop codons inferred in a read \cr
#'              \item {stops.max}: the larges number of stop codons inferred in a read \cr
#'              \item {stops.med}: the median number of stop codons inferred in a read \cr
#'          }  
#'
#' @export summarise.merged.read


summarise.merged.read <- function(merged.read){
  
  if(class(merged.read) != 'merged.read') { stop("merged.read must be a merged.read object")}
  
  m = merged.read
  reads = m$alignment[1:(length(m$alignment)-1)]
  read.lens = unlist(lapply(reads, function(x) length(DNAString(paste(as.matrix(del.gaps(x)), collapse = '')))))
  
  # the NAs allow us to take min/max/med and get NA back
  # TODO: reduce these to refer to only the reads that made it...
  insertions = m$indels$insertions
  if(is.null(insertions)){ insertions = NA }
  deletions = m$indels$deletions
  if(is.null(deletions)){ deletions = NA }
  stops = m$stop.codons$stop.codons
  if(is.null(stops)){ stops = NA }
  
  spc = nrow(m$secondary.peak.columns)
  if(is.null(spc)){ spc = 0 }
  
  summary = list(
    "consensus.len"     = length(m$consensus),
    "n.reads"           = length(reads),
    "distance.min"      = min(m$distance.matrix),
    "distance.max"      = max(m$distance.matrix),
    "distance.med"      = median(m$distance.matrix),
    "readlen.min"       = min(read.lens),
    "readlen.max"       = max(read.lens),
    "readlen.med"       = median(read.lens),
    "insertions.min"    = min(insertions),
    "insertions.max"    = max(insertions),
    "insertions.med"    = median(insertions),
    "deletions.min"     = min(deletions),
    "deletions.max"     = max(deletions),
    "deletions.med"     = median(deletions),
    "stops.min"         = min(stops),
    "stops.max"         = max(stops),
    "stops.med"         = median(stops),
    "secondary.peak.cols" = spc
  )
  
  return(summary)
  
}
####################
#' Merge sequence reads into a single consensus sequence
#' 
#' This function attempts to merge any number of unaligned reads into a single consensus sequence. It calculates a number of statistics that should help you decide whether a sensible consensus sequence exists for your data.\cr
#'
#' @param readset a DNAStringSet object with at least 2 reads
#' @param ref.aa.seq an amino acid reference sequence supplied as a string or an AAString object. If your sequences are protein-coding DNA seuqences, and you want to have frameshifts automatically detected and corrected, supply a reference amino acid sequence via this argument. If this argument is supplied, the sequences are then kept in frame for the alignment step. Fwd sequences are assumed to come from the sense (i.e. coding, or "+") strand.
#' @param minInformation minimum fraction of the sequences required to call a consensus sequence at any given position (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.5, implying that 3/4 of all reads must be present in order to call a consensus.
#' @param threshold Numeric giving the maximum fraction of sequence information that can be lost in the consensus sequence (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.5, implying that each consensus base must use at least half of the infomration at a given position. 
#' @param processors The number of processors to use, or NULL (the default) for all available processors
#' @param genetic.code Named character vector in the same format as GENETIC_CODE (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
#' @param accept.stop.codons TRUE/FALSE. TRUE (the defualt): keep all reads, regardless of whether they have stop codons; FALSE: reject reads with stop codons. If FALSE is selected, then the number of stop codons is calculated after attempting to correct frameshift mutations (if applicable).
#' @param reading.frame 1, 2, or 3. Only used if accept.stop.codons == FALSE. This specifies the reading frame that is used to determine stop codons. If you use a ref.aa.seq, then the frame should always be 1, since all reads will be shifted to frame 1 during frameshift correction. Otherwise, you should select the appropriate reading frame. 
#' 
#' @return A merged.read object, which is a list with the following components:
#'          \enumerate{
#'              \item {consensus}: the consensus sequence from the merged reads.\cr
#'              \item {alignment}: the alignment of all the reads, with the called consensus sequence. E.g. if the list was called 'merged.reads', you could use BrowseSeqs(merged.reads$alignment) to view the alignment.
#'              \item {differences}: a data frame of the number of pairwise differences between each read and the consensus sequence, as well as the number of bases in each input read that did not contribute to the consensus sequence. Can assist in detecting incorrect reads, or reads with a lot of errors.
#'              \item {distance.matrix}: a distance matrix of genetic distances (corrected with the JC model) between all of the input reads.
#'              \item {dendrogram}: a dendrogram depicting the distance.matrix. E.g. if the list was called 'merged.reads', you could use plot(merged.reads$dendrogram) to see the dendrogram.
#'              \item {indels}: if you specified a reference sequence via ref.aa.seq, then this will be a data frame describing the number of indels and deletions that were made to each of the input reads in order to correct frameshift mutations.
#'              \item {secondary.peak.columns}: this is a data frame with one row for each column in the alignment that contained more than one secondary peak. The data frame has three columns: the column number of the alignment; the number of secondary peaks in that column; and the bases (with IUPAC ambiguity codes representing secondary peak calls) in that column represented as a string.
#'          }
#'
#' @keywords merge, reads, alignment, sequence
#'
#' @export merge.reads
#'

merge.reads <- function(readset, ref.aa.seq = NULL, minInformation = 0.5, threshold = 0.5, processors = NULL, genetic.code = GENETIC_CODE, accept.stop.codons = TRUE, reading.frame = 1){
  
  # check input options
  processors = get.processors(processors)
  
  # this sometimes happens when we automate things
  if(length(readset) < 2) {return(NULL)}
  
  # Try to correct frameshifts in the input sequences 
  if(!is.null(ref.aa.seq)) {
    
    print("Correcting frameshifts in reads using amino acid reference sequence")
    corrected = CorrectFrameshifts(myXStringSet = readset, myAAStringSet = AAStringSet(ref.aa.seq), geneticCode = genetic.code, type = 'both', processors = processors, verbose = FALSE)
    readset = corrected$sequences
    indels = get.indel.df(corrected$indels)
    stops = as.numeric(unlist(mclapply(readset, count.stop.codons, reading.frame, genetic.code, mc.cores = processors)))
    stops.df = data.frame("read" = names(readset), "stop.codons" = stops)
    readset.lengths = unlist(lapply(readset, function(x) length(x)))
    readset = readset[which(readset.lengths>0)]
  }else{
    indels = NULL
    stops.df = NULL
  }
  
  # we might end up removing sequences during frameshift correction
  if(length(readset) < 2) {return(NULL)}
  
  
  # Remove reads with stop codons
  if(accept.stop.codons == FALSE){
    print("Removing reads with stop codons")
    if(is.null(ref.aa.seq)){ # otherwise we already did it above
      stops = as.numeric(unlist(mclapply(readset, count.stop.codons, reading.frame, genetic.code, mc.cores = processors)))
      stops.df = data.frame("read" = names(readset), "stop.codons" = stops)
    }
    old_length = length(readset)
    readset = readset[which(stops==0)]
    print(sprintf("%d reads with stop codons removed", old_length - length(readset)))
  }
  
  if(length(readset) < 2) {return(NULL)}
  
  print("Aligning reads")
  if(!is.null(ref.aa.seq)){
    aln = AlignTranslation(readset, geneticCode = genetic.code, processors = processors, verbose = FALSE)
  }else{
    aln = AlignSeqs(readset, processors = processors, verbose = FALSE)
  }
  
  if(is.null(names(aln))){
    names(aln) = paste("read", 1:length(aln), sep="_")
  }
  
  # call consensus
  print("Calling consensus sequence")
  consensus = ConsensusSequence(aln,
                                minInformation = minInformation,
                                includeTerminalGaps = TRUE,
                                ignoreNonBases = TRUE,
                                threshold = threshold,
                                noConsensusChar = "-",
                                ambiguity = TRUE
  )[[1]]
  
  print("Calculating differences between reads and consensus")
  diffs = mclapply(aln, n.pairwise.diffs, subject = consensus, mc.cores = processors)
  diffs = do.call(rbind, diffs)
  diffs.df = data.frame("name" = names(aln), "pairwise.diffs.to.consensus" = diffs[,1], "unused.chars" = diffs[,2])
  rownames(diffs.df) = NULL
  
  # get a dendrogram
  dist = DistanceMatrix(aln, correction = "Jukes-Cantor", penalizeGapLetterMatches = FALSE, processors = processors, verbose = FALSE)
  dend = IdClusters(dist, type = "dendrogram", processors = processors, verbose = FALSE)
  
  # add consensus to alignment
  aln2 = c(aln, DNAStringSet(consensus))
  names(aln2)[length(aln2)] = "consensus"
  # strip gaps from consensus (must be an easier way!!)
  consensus.gapfree = DNAString(paste(as.matrix(del.gaps(consensus)), collapse = ''))
  
  # count columns in the alignment with >1 coincident secondary peaks
  sp.df = count.coincident.sp(aln, processors = processors)
  
  merged.read = list("consensus" = consensus.gapfree, 
                     "alignment" = aln2, 
                     "differences" = diffs.df, 
                     "distance.matrix" = dist,
                     "dendrogram" = dend,
                     "indels" = indels,
                     "stop.codons" = stops.df,
                     "secondary.peak.columns" = sp.df)
  
  class(merged.read) = "merged.read"
  
  return(merged.read)
}


get.indel.df <- function(indel.list){
  
  r = lapply(indel.list, indel.row)
  indel.df = data.frame(matrix(unlist(r), byrow = T, nrow = length(indel.list)))
  indel.df = cbind(names(indel.list), indel.df)
  names(indel.df) = c('read', 'insertions', 'deletions', 'distance')    
  return(indel.df)
  
}

indel.row <- function(row){
  
  n.ins = length(row$insertions)
  n.del = length(row$deletions)
  dist = row$distance
  return(c(n.ins, n.del, dist))
}

n.pairwise.diffs <- function(pattern, subject){
  # pairwise differences assuming pattern and subject are aligned
  comp = compareStrings(pattern, subject)
  qs = str_count(comp, '\\?')
  ps = str_count(comp, '\\+')
  return(c(qs, ps))
}

count.coincident.sp <- function(aln, processors){
  # make a data frame of columns in the alignment that have
  # more than one secondary peak
  
  is = 1:aln@ranges@width[1]
  r = mclapply(is, one_ambiguous_column, aln=aln, mc.cores = processors)
  r = Filter(Negate(is.null), r)
  
  if(length(r)>0){
    r = as.data.frame(matrix(unlist(r), nrow=length(r), byrow=T))
    names(r) = c('column.number', 'ambiguities', 'column')
    return(r)
  }else{
    return(NULL)
  }
  
}

ambiguous = names(IUPAC_CODE_MAP)[5:length(names(IUPAC_CODE_MAP))]

one_ambiguous_column <- function(i, aln){
  col = as.character(subseq(aln, i, i))
  str = paste(col, sep="", collapse="")
  amb = sum(col %in% ambiguous)
  if(amb>1){
    return(c(i, amb, str))
  }
}
####' Make a readset
#' 
#' @param fwd.fnames a list of full file paths to forward reads from ab1 files (i.e. those that do not need to be reverse-complemented). 
#' @param rev.fnames a list of full file paths to reverse reads from ab1 files (i.e. those that *do* need to be reverse-complemented). 
#' @param trim TRUE/FALSE trim sequences based on quality while creating readgroup? If TRUE, the trim.mott function is applied to each sequence before inclusion in the readgroup. Note, trimming only works if the raw data are stored in ab1 files with appropriate information.
#' @param trim.cutoff value passed to trim.mott as quality cutoff for sequencing trimming, only used if 'trim' == TRUE
#' @param max.secondary.peaks reads with more secondary peaks than this will not be included in the readset. The default (NULL) is to include all reads regardless of secondary peaks 
#' @param secondary.peak.ratio Only applies if max.secondary.peaks is not NULL. The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are counted. Those below the ratio are not. 
#' @param min.length reads shorter than this will not be included in the readset. The default (1) means that all reads with length of 1 or more will be included.
#' @param processors The number of processors to use, or NULL (the default) for all available processors
#'
#' @return A set of unaligned reads as a DNAstringset object, names are the input file paths. Note that secondary peaks are encoded in these reads as ambiguous bases using IUPAC ambiguity codes.
#'
#' @export make.readset


make.readset <- function(fwd.fnames, rev.fnames, trim = TRUE, trim.cutoff = 0.0001, max.secondary.peaks = NULL, secondary.peak.ratio = 0.33, min.length = 1, processors = NULL){
  
  processors = get.processors(processors)
  
  fwd.dat = mclapply(fwd.fnames, loadread, trim, trim.cutoff, revcomp = FALSE, max.secondary.peaks = max.secondary.peaks, secondary.peak.ratio = secondary.peak.ratio, min.length = min.length, processors = 1, mc.cores = processors)
  rev.dat = mclapply(rev.fnames, loadread, trim, trim.cutoff, revcomp = TRUE,  max.secondary.peaks = max.secondary.peaks, secondary.peak.ratio = secondary.peak.ratio, min.length = min.length, processors = 1, mc.cores = processors)
  
  fwd.reads = lapply(fwd.dat, function(x) x[["read"]])
  rev.reads = lapply(rev.dat, function(x) x[["read"]])
  names(fwd.reads) = fwd.fnames
  names(rev.reads) = rev.fnames
  
  # remove the NULL reads
  all.reads = c(fwd.reads, rev.reads)
  all.reads = Filter(Negate(is.null), all.reads)
  readset = DNAStringSet(all.reads)
  
  # build the summary data frame just as in summarise.abi.folder
  fwd.summaries = lapply(fwd.dat, function(x) x[["summary"]])
  rev.summaries = lapply(rev.dat, function(x) x[["summary"]])
  all.summaries = c(fwd.summaries, rev.summaries)
  all.summaries = do.call(rbind, all.summaries)
  abi.fnames = unlist(c(fwd.fnames, rev.fnames))
  folder.names = basename(dirname(abi.fnames))
  file.names = basename(abi.fnames)
  all.summaries = cbind.data.frame("file.path" = as.character(abi.fnames), "folder.name" = as.character(folder.names), "file.name" = file.names, all.summaries, stringsAsFactors = FALSE)
  
  
  used.reads = names(readset)
  all.summaries$read.included.in.readset = all.summaries$file.path %in% used.reads
  
  
  return(list("readset" = readset, "read.summaries" = all.summaries))
  
}

loadread <- function(fname, trim, trim.cutoff, revcomp, max.secondary.peaks, secondary.peak.ratio, min.length, processors){
  
  read.abi = read.abif(fname)
  
  s = summarise.abi.file(read.abi, trim.cutoff, secondary.peak.ratio, processors = processors)
  
  summary = s$summary
  
  # here we store the secondary peaks by storing them in a single sequence
  # as ambiguity codes. Note, this is version of a read that we use.
  # So in this package, a read has an ambiguit whereever there's a 
  # secondary peak
  d = c(DNAStringSet(s$read@primarySeq), DNAStringSet(s$read@secondarySeq))
  read = ConsensusSequence(d)[[1]]
  
  if(trim == TRUE){
    trim.start = summary["trim.start"]
    trim.finish = summary["trim.finish"]
    sp = summary["trimmed.secondary.peaks"]
    
  }else if(trim == FALSE){
    trim.start = 1
    trim.finish = length(read)
    sp = summary["raw.secondary.peaks"]
  }
  
  # FILTER read based on user specified limits
  read = read[trim.start:trim.finish]
  
  if(!is.null(max.secondary.peaks)){
    if(sp > max.secondary.peaks){
      read = NULL
    }
  }
  
  if(length(read) < min.length){
    read = NULL
  }    
  
  if(!is.null(read)) {
    if(revcomp == TRUE){
      read = reverseComplement(read)
    }
  }
  return(list('read' = read, summary = summary))
  
}


##########

#' Automatically make consensus sequences by grouping .ab1 files by name.
#' 
#' @param input.folder The parent folder of all of the reads contained in ab1 files you wish to analyse. Subfolders will be scanned recursively.
#' @param forward.suffix the suffix of the filenames for forward reads, i.e. reads that do not need to be reverse-complemented. Include the full suffix, e.g. "forward.ab1".
#' @param reverse.suffix the suffix of the filenames for reverse reads, i.e. reads that *do* need to be reverse-complemented. Include the full suffix, e.g. "reverse.ab1".
#' @param min.reads The minimum number of reads required to make a consensus sequence, must be 2 or more (default 2). 
#' @param trim TRUE/FALSE trim sequences based on quality while creating readgroup? If TRUE, the trim.mott function is applied to each sequence before inclusion in the readgroup. Note, trimming only works if the raw data are stored in ab1 files with appropriate information.
#' @param trim.cutoff value passed to trim.mott as quality cutoff for sequencing trimming, only used if 'trim' == TRUE
#' @param max.secondary.peaks reads with more secondary peaks than this will not be included in the readset. The default (NULL) is to include all reads regardless of secondary peaks 
#' @param secondary.peak.ratio Only applies if max.secondary.peaks is not NULL. The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are counted. Those below the ratio are not. 
#' @param min.length reads shorter than this will not be included in the readset. The default (20) means that all reads with length of 20 or more will be included. Note that this is the length of a read after it has been trimmed.
#' @param ref.aa.seq an amino acid reference sequence supplied as a string or an AAString object. If your sequences are protein-coding DNA seuqences, and you want to have frameshifts automatically detected and corrected, supply a reference amino acid sequence via this argument. If this argument is supplied, the sequences are then kept in frame for the alignment step. Fwd sequences are assumed to come from the sense (i.e. coding, or "+") strand.
#' @param minInformation minimum fraction of the sequences required to call a consensus sequence at any given position (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.5 implying that 3/4 of all reads must be present in order to call a consensus.
#' @param threshold Numeric giving the maximum fraction of sequence information that can be lost in the consensus sequence (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.5, implying that each consensus base can ignore at most 50 percent of the information at a given position. 
#' @param genetic.code Named character vector in the same format as GENETIC_CODE (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
#' @param accept.stop.codons TRUE/FALSE. TRUE (the defualt): keep all reads, regardless of whether they have stop codons; FALSE: reject reads with stop codons. If FALSE is selected, then the number of stop codons is calculated after attempting to correct frameshift mutations (if applicable).
#' @param reading.frame 1, 2, or 3. Only used if accept.stop.codons == FALSE. This specifies the reading frame that is used to determine stop codons. If you use a ref.aa.seq, then the frame should always be 1, since all reads will be shifted to frame 1 during frameshift correction. Otherwise, you should select the appropriate reading frame. 
#' @param processors The number of processors to use, or NULL (the default) for all available processors
#'
#' @return Note that the consensus.tree should not be used for inference, but only as a guide to help judge problems in the consensus sequences. This tree is built using ape's njs() function, and negative branch lengths are then converted to their absolute values. The latter process aids viewing clarity (so should help when looking for problems in the consensus seuqences), but has little biological validity. 

#' @export make.consensus.seqs

make.consensus.seqs <- function(input.folder, forward.suffix, reverse.suffix, min.reads = 2, trim = TRUE, trim.cutoff = 0.0001, min.length = 20, max.secondary.peaks = NULL, secondary.peak.ratio = 0.33, ref.aa.seq = NULL, minInformation = 0.5, threshold = 0.5, genetic.code = GENETIC_CODE, accept.stop.codons = TRUE, reading.frame = 1,  processors = NULL){
  
  processors = get.processors(processors)
  
  rs = make.readsets(input.folder = input.folder, 
                     forward.suffix = forward.suffix, 
                     reverse.suffix = reverse.suffix, 
                     trim = trim, 
                     trim.cutoff = trim.cutoff, 
                     min.length = min.length, 
                     max.secondary.peaks = max.secondary.peaks, 
                     secondary.peak.ratio = secondary.peak.ratio,
                     processors = processors
  )
  
  # Process readset output, and filter based on number of reads
  print(sprintf("Filtering readsets with <%d reads...", min.reads))
  readsets = rs$readsets
  read.summaries = rs$read.summaries
  readset.lengths = unlist(lapply(readsets, function(x) length(x)))
  valid.readsets = readsets[which(readset.lengths >= min.reads)]
  valid.readset.lengths = unlist(lapply(valid.readsets, function(x) length(x)))
  
  print(sprintf("After filtering, %d of %d readsets remain", length(valid.readsets), length(readsets)))
  
  
  if(median(valid.readset.lengths) > length(valid.readsets)){
    # better to do readgroups sequentially, but parallelise each
    mc.cores = 1
    c.processors = processors
  }else{
    # better to do readgroups in parallel, but sequentially within each
    mc.cores = processors
    c.processors = 1
  }
  
  print("Building consensus sequences...")
  
  consensi = mclapply(valid.readsets,
                      merge.reads,
                      ref.aa.seq = ref.aa.seq, 
                      minInformation = minInformation, 
                      threshold = threshold, 
                      processors = c.processors, 
                      genetic.code = genetic.code, 
                      accept.stop.codons = accept.stop.codons, 
                      reading.frame = reading.frame,
                      mc.cores = mc.cores
  )
  
  # make the set of consensus sequences
  consensus.seqs = lapply(consensi, function(x) x$consensus)
  # Some are null, becuase reads can be removed for e.g. stop codons
  consensus.seqs = Filter(Negate(is.null), consensus.seqs)
  consensus.set  = DNAStringSet(consensus.seqs)
  
  print(sprintf("Successfully built %d consensus sequences", length(consensus.set)))
  
  # Group the summaries
  print("Summarising consensus sequences...")
  consensus.summaries = mclapply(consensi, summarise.merged.read, mc.cores = processors)
  consensus.summaries = do.call(rbind, consensus.summaries)
  consensus.summaries = cbind("consensus.name" = row.names(consensus.summaries), consensus.summaries)
  row.names(consensus.summaries) = NULL
  consensus.summaries = data.frame(consensus.summaries)
  
  # which reads made it to the consensus sequence
  # careful, this list also has 'consensus' in it a lot
  used.reads = unlist(lapply(consensi, function(x) as.character(names(x$alignment))))
  read.summaries$read.included.in.consensus = read.summaries$file.path %in% used.reads
  
  # a column for successful consensus sequence
  success = names(consensi)
  success.indices = which(read.summaries$readset.name %in% success)
  read.summaries$consensus.name = NA
  read.summaries$consensus.name[success.indices] = as.character(read.summaries$readset.name[success.indices])
  
  # Now we add more data from the read summaries
  used.read.summaries = subset(read.summaries, read.included.in.consensus==TRUE)    
  rsm = melt(used.read.summaries, id.vars = c("consensus.name", "folder.name", "file.name", "readset.name", "file.path", "read.included.in.readset", "read.included.in.consensus"))
  meds = dcast(rsm, consensus.name ~ variable, median)
  maxs = dcast(rsm, consensus.name ~ variable, max)
  mins = dcast(rsm, consensus.name ~ variable, min)
  more.summaries = data.frame("consensus.name" = as.character(meds$consensus.name),
                              "raw.secondary.peaks.min" = mins$raw.secondary.peaks, 
                              "raw.secondary.peaks.max" = maxs$raw.secondary.peaks,
                              "raw.secondary.peaks.med" = meds$raw.secondary.peaks,
                              "trimmed.secondary.peaks.min" = mins$trimmed.secondary.peaks, 
                              "trimmed.secondary.peaks.max" = maxs$trimmed.secondary.peaks,
                              "trimmed.secondary.peaks.med" = meds$trimmed.secondary.peaks,
                              "raw.mean.quality.min" = mins$raw.mean.quality, 
                              "raw.mean.quality.max" = maxs$raw.mean.quality,
                              "raw.mean.quality.med" = meds$raw.mean.quality,
                              "trimmed.mean.quality.min" = mins$trimmed.mean.quality, 
                              "trimmed.mean.quality.max" = maxs$trimmed.mean.quality,
                              "trimmed.mean.quality.med" = meds$trimmed.mean.quality
  )
  
  consensus.summaries = merge(consensus.summaries, more.summaries, by = "consensus.name", sort = FALSE)
  
  consensus.summaries$consensus.name = as.character(consensus.summaries$consensus.name)
  
  # align the consensus sequences
  if(length(consensus.set)>1){
    print("Aligning consensus sequences...")
    if(!is.null(ref.aa.seq)){
      aln = AlignTranslation(consensus.set, geneticCode = genetic.code, processors = processors, verbose = FALSE)
    }else{
      aln = AlignSeqs(consensus.set, processors = processors, verbose = FALSE)
    }
    
    # make a rough NJ tree. Labels are rows in the summary df
    print("Building tree of consensus sequences...")
    neat.labels = match(names(aln), 
                        as.character(consensus.summaries$consensus.name)
    )
    aln2 = aln
    names(aln2) = neat.labels
    aln.bin = as.DNAbin(aln2)
    aln.dist = dist.dna(aln.bin, pairwise.deletion = TRUE)
    
    # Sometimes it's impossible to make a tree...
    aln.tree = NULL
    try({
      aln.tree = bionjs(aln.dist)
      
      # deal with -ve branches
      # This is not necessarily accurate, but it is good enough to judge seuqences using the tree
      aln.tree$edge.length[which(aln.tree$edge.length<0)] = abs(aln.tree$edge.length[which(aln.tree$edge.length<0)])            },
      silent = TRUE
    )
  }else{
    aln = NA
    aln.tree = NA
  }
  return(list("read.summaries" = read.summaries, 
              "merged.reads" = consensi, 
              "consensus.summaries" = consensus.summaries, 
              "consensus.sequences" = consensus.set, 
              "consensus.alignment" = aln,
              "consensus.tree" = aln.tree
  )
  )
  
}
###############
#' Automatically load readsets by grouping .ab1 files by name.
#'
#' @param input.folder The parent folder of all of the reads contained in ab1 files you wish to analyse. Subfolders will be scanned recursively.
#' @param forward.suffix the suffix of the filenames for forward reads, i.e. reads that do not need to be reverse-complemented. Include the full suffix, e.g. "forward.ab1".
#' @param reverse.suffix the suffix of the filenames for reverse reads, i.e. reads that *do* need to be reverse-complemented. Include the full suffix, e.g. "reverse.ab1".
#' @param trim TRUE/FALSE trim sequences based on quality while creating readgroup? If TRUE, the trim.mott function is applied to each sequence before inclusion in the readgroup. Note, trimming only works if the raw data are stored in ab1 files with appropriate information.
#' @param trim.cutoff value passed to trim.mott as quality cutoff for sequencing trimming, only used if 'trim' == TRUE
#' @param max.secondary.peaks reads with more secondary peaks than this will not be included in the readset. The default (NULL) is to include all reads regardless of secondary peaks
#' @param secondary.peak.ratio Only applies if max.secondary.peaks is not NULL. The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are counted. Those below the ratio are not.
#' @param min.length reads shorter than this will not be included in the readset. The default (1) means that all reads with length of 1 or more will be included.
#' @param processors The number of processors to use, or NULL (the default) for all available processors
#'
#' @export make.readsets

make.readsets <- function(input.folder, forward.suffix, reverse.suffix, trim = TRUE, trim.cutoff = 0.0001, processors = NULL, min.length = 1, max.secondary.peaks = NULL, secondary.peak.ratio = 0.33){
  
  processors = get.processors(processors)
  
  print("Looking for .ab1 files...")
  abi.files = list.files(input.folder, pattern = "\\.ab1$", full.names = T, recursive = T)
  
  # get a set of unique filenames after removing suffixes
  print(sprintf("Grouping %d files into sets by filename...", length(abi.files)))
  group.dataframe = get.group.dataframe(abi.files, forward.suffix, reverse.suffix)
  groups = unique(group.dataframe$group)
  
  print(sprintf("%d of the files matched either the forward or reverse suffixes", nrow(group.dataframe)))
  
  print(sprintf("Grouped these into %d sets of files", length(groups)))
  
  # load full file paths for readgroups based on unique filenames
  print("Loading .ab1 files...")
  readset.fnames = mclapply(groups,
                            get.readgroup.fnames,
                            group.dataframe = group.dataframe,
                            forward.suffix = forward.suffix,
                            reverse.suffix = reverse.suffix,
                            mc.cores = processors)
  
  # how we parallelise depends on how many readgroups there are
  if(length(readset.fnames[[1]]) > length(readset.fnames)){
    # better to do readgroups sequentially, but parallelise each
    mc.cores = 1
    c.processors = processors
  }else{
    # better to do readgroups in parallel, but sequentially within each
    mc.cores = processors
    c.processors = 1
  }
  
  print("Constructing readsets...")
  rs = mclapply(readset.fnames, make.readset.from.list,
                trim = trim,
                trim.cutoff = trim.cutoff,
                max.secondary.peaks = max.secondary.peaks,
                secondary.peak.ratio = secondary.peak.ratio,
                min.length = min.length,
                processors = c.processors,
                mc.cores = mc.cores
  )
  
  names(rs) = groups
  
  print("Building read summaries...")
  readsets = mclapply(rs, function(x) x[["readset"]], mc.cores = processors)
  summaries = mclapply(rs, function(x) x[["read.summaries"]], mc.cores = processors)
  summaries = do.call(rbind, summaries)
  rownames(summaries) = NULL
  
  names(group.dataframe) = c('file.path', 'readset.name')
  
  summaries = merge(summaries, group.dataframe, by = 'file.path')
  
  return(list("readsets" = readsets, "read.summaries" = summaries))
  
}


make.readset.from.list <- function(fnames, trim, trim.cutoff, max.secondary.peaks, secondary.peak.ratio, min.length, processors){
  # this just unpacks the fnames and sends them on, so I can use mclapply
  fwd.reads = fnames$forward.reads
  rev.reads = fnames$reverse.reads
  
  readset = make.readset(fwd.reads, rev.reads,
                         trim = trim,
                         trim.cutoff = trim.cutoff,
                         max.secondary.peaks = max.secondary.peaks,
                         secondary.peak.ratio = secondary.peak.ratio,
                         min.length = min.length,
                         processors = processors)
  
  return(readset)
  
}

get.group.dataframe <- function(fname.list, forward.suffix, reverse.suffix){
  
  files.cleaned = fname.list
  f.matches = str_match(files.cleaned, forward.suffix)
  f.indices = which(!is.na(f.matches))
  r.matches = str_match(files.cleaned, reverse.suffix)
  r.indices = which(!is.na(r.matches))
  
  #ONLY keep files that do match the suffixes:
  keep = c(f.indices, r.indices)
  files.cleaned = files.cleaned[keep]
  
  # try with and without the .ab1 after the suffixes, just in case
  files.cleaned = str_replace(files.cleaned, forward.suffix, replacement = "")
  files.cleaned = str_replace(files.cleaned, reverse.suffix, replacement = "")
  
  return(data.frame("file.path" = fname.list[keep], "group" = files.cleaned))
  
}


load.sangerseqs <- function(filenames){
  seqs = lapply(filenames, readsangerseq)
  return(seqs)
}


get.readgroup.fnames <- function(group, group.dataframe, forward.suffix, reverse.suffix){
  
  readgroup.fnames = as.character(group.dataframe$file.path[which(group.dataframe$group == group)])
  
  f.matches = str_match(readgroup.fnames, forward.suffix)
  f.indices = which(!is.na(f.matches))
  fwd.fnames = readgroup.fnames[f.indices]
  
  r.matches = str_match(readgroup.fnames, reverse.suffix)
  r.indices = which(!is.na(r.matches))
  rev.fnames = readgroup.fnames[r.indices]
  
  readgroup.fnames = list("forward.reads" = fwd.fnames, "reverse.reads" = rev.fnames)
  
  return(readgroup.fnames)
  
}
#############
#' Count stop codons in a DNA sequence
#'
#' @param sequence a DNAString object
#' @param reading.frame a number from 1 to 3 denoting the reading frame of the sequence
#' @param genetic.code Named character vector in the same format as GENETIC_CODE (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
#'
#' @export count.stop.codons
#'

count.stop.codons <- function(sequence, reading.frame = 1, genetic.code = GENETIC_CODE){
  
  if(!reading.frame %in% c(1,2,3)){ stop("reading.frame must be 1, 2, or 3")}
  if(class(sequence)!='DNAString'){ stop("sequence must be a DNAString object")}
  if(!("*" %in% genetic.code)) { stop("Your genetic code does not specify any stop codons")}
  
  l = length(sequence) + 1 - reading.frame
  
  if(l < 3){
    warning(sprintf("Cannot calculate stop codons on sequence of length %d in reading frame %d", 
                    length(sequence), reading.frame))
    return(NULL)
  }
  
  # this comes almost straight from the BioStrings manual
  tri = trinucleotideFrequency(sequence[reading.frame:length(sequence)], step=3)
  
  names(tri) <- genetic.code[names(tri)]
  
  freqs = sapply(split(tri, names(tri)), sum)
  
  stops = freqs["*"]
  
  return(as.numeric(stops))
}
##############
#' Create a detailed summary of a single ABI sequencing file
#' 
#' @param seq.abif an abif.seq s4 object from the sangerseqR package
#' @param trim.cutoff the cutoff at which you consider a base to be bad. This works on a logarithmic scale, such that if you want to consider a score of 10 as bad, you set cutoff to 0.1; for 20 set it at 0.01; for 30 set it at 0.001; for 40 set it at 0.0001; and so on. Contiguous runs of bases below this quality will be removed from the start and end of the sequence. Given the high quality reads expected of most modern ABI sequencers, the defualt is 0.0001.
#' @param secondary.peak.ratio the ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are not. 
#' @param output.folder If output.folder is NA (the default) no files are written. If a valid folder is provided, two files are written to that folder: a .csv file of the secondary peaks (see description below) and a .pdf file of the chromatogram.
#' @param file.prefix If output.folder is specified, this is the prefix which will be appended to the .csv and the .pdf file. The default is "seq".
#' @param processors The number of processors to use, or NULL (the default) for all available processors
#'
#' @return A numeric vector including:
#'          \enumerate{
#'              \item {raw.length}: the length of the untrimmed sequence, note that this is the sequence after conversion to a sangerseq object, and then the recalling the bases with MakeBaseCalls from the sangerseqR package\cr
#'              \item {trimmed.length}: the length of the trimmed sequence, after trimming using trim.mott from this package and the parameter supplied to this function \cr
#'              \item {trim.start}: the start position of the good sequence, see trim.mott for more details\cr
#'              \item {trim.finish}: the finish position of the good sequence, see trim.mott for more details\cr
#'              \item {raw.secondary.peaks}: the number of secondary peaks in the raw sequence, called with the secondary.peaks function from this package and the parameters supplied to this function \cr
#'              \item {trimmed.secondary.peaks}: the number of secondary peaks in the trimmed sequence, called with the secondary.peaks function from this package and the parameters supplied to this function \cr
#'              \item {raw.mean.quality}: the mean quality score of the raw sequence \cr
#'              \item {trimmed.mean.quality}: the mean quality score of the trimmed sequence \cr
#'              \item {raw.min.quality}: the minimum quality score of the raw sequence \cr
#'              \item {trimmed.min.quality}: the minimum quality score of the trimmed sequence \cr
#'          }  
#'
#' @export summarise.abi.file

summarise.abi.file <- function(seq.abif, trim.cutoff = 0.0001, secondary.peak.ratio = 0.33, output.folder = NA, prefix = "seq", processors = NULL){
  
  seq.sanger = sangerseq(seq.abif)
  
  # first we get the secondary peaks
  # note that the secondary peaks correspond to the seq.sanger object AFTER we
  # have called makeBaseCalls. And that this means the trim locations and the 
  # secondary peak locations do not match, since makeBaseCalls usually calls
  # fewer bases than the standard ABI calls.
  secondary.peaks.data = secondary.peaks(seq.sanger, secondary.peak.ratio, output.folder, prefix, processors = processors)
  secondary.peaks = secondary.peaks.data$secondary.peaks
  seq.sanger = secondary.peaks.data$read
  
  # now we trim the sequence
  trims = trim.mott(seq.abif, cutoff = trim.cutoff)
  qual = seq.abif@data$PCON.2
  qual.trimmed = qual[trims$start:trims$finish]
  
  if(length(qual.trimmed)==0){qual.trimmed = c(NA)} # so we can summarise later
  
  # now we fix up the trim locations to correspond to the sangerseq primary seq object
  if(trims$start==1 && trims$finish==nchar(as.character(seq.abif@data$PBAS.2))){
    # there's nothing to do, because we didn't trim anything
    trim.start = 1
    trim.finish = length(primarySeq(seq.sanger))
    
  }else{
    trims.fixed = fix.trims(trims, seq.sanger, seq.abif, processors)
    trim.start = trims.fixed$start
    trim.finish = trims.fixed$finish
  }
  
  # get trimmed and untrimmed version of raw data
  seq.trimmed = seq.sanger@primarySeq[trim.start:trim.finish]
  secondary.peaks.trimmed = subset(secondary.peaks, position >= trim.start & position <= trim.finish)
  
  print(qual.trimmed)
  read.summary = c("raw.length"                       = length(seq.sanger@primarySeq), 
                   "trimmed.length"                   = length(seq.trimmed),
                   "trim.start"                       = trim.start,
                   "trim.finish"                      = trim.finish,
                   "raw.secondary.peaks"              = nrow(secondary.peaks),
                   "trimmed.secondary.peaks"          = nrow(secondary.peaks.trimmed),
                   "raw.mean.quality"                 = mean(qual),
                   "trimmed.mean.quality"             = mean(qual.trimmed),
                   "raw.min.quality"                  = min(qual),
                   "trimmed.min.quality"              = min(qual.trimmed)                     
  )
  
  return(list("summary" = read.summary, "read" = seq.sanger))
  
}


fix.trims <- function(trims, seq.sanger, seq.abif, processors){
  
  # transfer trim locations from one sequence (denoted in the trims list, and which 
  # correspond to the seq.abif object to another 
  # the primarySeq(seq.sanger) from the seq.sanger object
  
  if(trims$start == 0 & trims$finish == 0){
    # no need to do anything fancy here...
    return(trims)
  }
  
  # 1. First we trim the original sequence
  original.seq = seq.abif@data$PBAS.2
  
  original.trimmed = substring(original.seq, trims$start, trims$finish)
  
  # 2. Align the original and recalled sequences
  recalled = primarySeq(seq.sanger, string = TRUE)
  seqs = DNAStringSet(c(original.trimmed, recalled))
  pa = AlignSeqs(seqs, iterations = 0, refinements = 0, verbose = FALSE, processors = processors)
  
  # 3. Get the sequence out, and find the first and last gaps.
  aligned.trimmed = as.character(pa[[1]])
  not.gaps = str_locate_all(aligned.trimmed, pattern = "[^-]")[[1]][,1]
  
  start = min(not.gaps)
  finish = max(not.gaps)
  
  if(start < 1){start = 1}
  if(finish > nchar(recalled)){finish = nchar(recalled)}
  
  return(list("start" = start, "finish" = finish))
}
##########3
trim.mott <- function(abif.seq, cutoff = 0.0001){
  
  if(class(cutoff)!='numeric' | cutoff < 0){
    stop("cutoff must be a number of at least 0")
  }
  
  if(class(abif.seq)!='abif'){
    stop("abif.seq must be an 'abif' object from the sangerseqR package")
  }
  
  
  
  abif.seq = abif.seq@data
  start = FALSE # flag for starting position of trimmed sequence
  trim_start = 0 # init start index
  
  seqlen = nchar(abif.seq$PBAS.2)
  qual = abif.seq$PCON.2
  
  # calculate base score 
  score_list = cutoff - (10 ** (qual / -10.0))
  
  # calculate cummulative score 
  # if cumulative value < 0, set it to 0 
  # the BioPython implementation always trims the first base, 
  # this implementation does not. 
  score = score_list[1]
  if(score < 0){ 
    score = 0 
  }else{
    trim_start = 1
    start = TRUE
  }
  
  cummul_score = c(score)
  
  for(i in 2:length(score_list)){
    score = cummul_score[length(cummul_score)] + score_list[i]
    if(score <= 0){
      cummul_score = c(cummul_score, 0)
    }else{
      cummul_score = c(cummul_score, score)
      if(start == FALSE){
        # trim_start = value when cummulative score is first > 0 
        trim_start = i
        start = TRUE
      }
    }
    
    # trim_finish = index of highest cummulative score, 
    # marking the end of sequence segment with highest cummulative score 
    trim_finish = which.max(cummul_score)
    
  }
  
  # fix an edge case, where all scores are worse than the cutoff
  # in this case you wouldn't want to keep any bases at all
  if(sum(cummul_score)==0){trim_finish = 0}
  
  return(list("start" = trim_start, "finish" = trim_finish))
  
}
############
#' Check for secondary peaks in a sangerseq object
#' 
#' This function finds and reports secondary peaks in a sangerseq object. It returns a table of secondary peaks, and optionally saves an annotated chromatogram and a csv file of the peak locations.
#' 
#' @param s a sangerseq s4 object from the sangerseqR package
#' @param ratio the ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are not. 
#' @param output.folder If output.folder is NA (the default) no files are written. If a valid folder is provided, two files are written to that folder: a .csv file of the secondary peaks (see description below) and a .pdf file of the chromatogram.
#' @param file.prefix If output.folder is specified, this is the prefix which will be appended to the .csv and the .pdf file. The default is "seq".
#' @param processors The number of processors to use, or NULL (the default) for all available processors
#' 
#' @return A list with two elements:
#'          \enumerate{
#'              \item {secondary.peaks}: a data frame with one row per secondary peak above the ratio, and three columns: "position" is the position of the secondary peak relative to the primary sequence; "primary.basecall" is the primary base call; "secondary.basecall" is the secondary basecall. \cr
#'              \item {read}: the input sangerseq s4 object after having the makeBaseCalls() function from sangerseqR applied to it. This re-calls the primary and secondary bases in the sequence, and resets a lot of the internal data.
#'          }
#'
#' @keywords chromatogram, peak, mismatch
#'
#' @export secondary.peaks
#'

secondary.peaks <- function(s, ratio = 0.33, output.folder = NA, file.prefix = "seq", processors = NULL){
  
  # make secondary basecalls, and align them to the original sequence
  basecalls = makeBaseCalls(s, ratio = ratio)
  
  primary = primarySeq(basecalls, string = TRUE)
  secondary = secondarySeq(basecalls, string = TRUE)
  
  # perhaps we don't need to align...
  #seqs = DNAStringSet(c(primary, secondary))
  # these seuqences should be VERY similar...
  #pa = AlignSeqs(seqs, iterations = 0, refinements = 0, verbose = FALSE, processors = processors)
  
  
  # NB: it would seem to make more sense to use mismatchTable here, 
  # but I recoded it this way because mismatchTable had a bug.
  comp = compareStrings(primary, secondary)
  diffs = str_locate_all(pattern ='\\?',comp)[[1]][,1]
  primary.vector = strsplit(primary, split="")[[1]]
  secondary.vector = strsplit(secondary, split="")[[1]]
  
  primary.basecall    = primary.vector[diffs]
  secondary.basecall  = secondary.vector[diffs]
  
  r = data.frame("position" = diffs, "primary.basecall" = primary.basecall, "secondary.basecall" = secondary.basecall)
  
  if(!is.na(output.folder)){
    if(dir.exists(output.folder)){
      chromname = paste(file.prefix, "_", "chromatogram.pdf", sep='')
      tablename = paste(file.prefix, "_", "secondary_peaks.csv", sep='')
      chrom = chromatogram(basecalls, height = 2, showcalls = 'both', filename = file.path(output.folder, chromname))
      write.csv(r, file = file.path(output.folder, tablename))
    }else{
      warning(sprintf("Couldn't find directory '%s', no files saved", output.folder))
    }
  }
  
  return(list("secondary.peaks" = r, "read" = basecalls))
  
}
#2018 GitHub, Inc.
###########