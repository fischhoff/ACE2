ACE2
================
Ilya Fischhoff
4/29/2020

\#\#\#\#\#install packages

    ## 
    ## Attaching package: 'seqinr'

    ## The following object is masked from 'package:plyr':
    ## 
    ##     count

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:plyr':
    ## 
    ##     rename

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:plyr':
    ## 
    ##     desc

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'XVector'

    ## The following object is masked from 'package:plyr':
    ## 
    ##     compact

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:seqinr':
    ## 
    ##     translate

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     complement

    ## The following objects are masked from 'package:seqinr':
    ## 
    ##     as.alignment, consensus

    ## 
    ## Attaching package: 'phylotools'

    ## The following object is masked from 'package:seqinr':
    ## 
    ##     read.fasta

    ## 
    ## Attaching package: 'data.table'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     shift

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, second

\#\#check that set of mammals is contained in set of vertebrates

``` r
A <- read.csv("ACE2GenBank - GenBank.csv")
names(A)[names(A)=="AccessionNum"]="AccessionNumGene"
save(A, file= "A.Rdata")
```

Get protein based on accession ID

``` r
refIDs <- A$AccessionNum
 SEQS <- foreach(i = refIDs, .errorhandling = "pass") %do% {
  SRCH <- entrez_search(db = "nuccore", term = i)
  entrez_fetch(db = "protein", id = entrez_link(dbfrom = "nuccore", id = SRCH[[1]], db = "protein")$links[[1]], rettype = "fasta")
 }
 
inds_missing = NULL
for (a in 1:length(refIDs)){
  if( (SEQS[[a]][1]=="subscript out of bounds")==TRUE){
    inds_missing = c(inds_missing, a)
  }
} 

inds_empty = NULL
for (a in 1:length(refIDs)){
  if( (str_detect(SEQS[[a]][1],"empty"))==TRUE){
    inds_empty = c(inds_empty, a)
  }
} 

inds_none = c(inds_missing, inds_empty)


inds_okay = setdiff(seq(1, length(refIDs)),inds_none)

length(inds_okay)+length(inds_none)== length(refIDs)
```

    ## [1] TRUE

``` r
SEQS_ok = SEQS[inds_okay]

out = NULL
for (a in 1:length(inds_okay)){
  output= str_split(SEQS_ok[a],"\n")
  output = data.frame(output)
  lines = dim(output)[1]
  name = output[1,]
  acc = str_split(name, " ")
  accession = acc[[1]][1]
  accession = str_replace(accession, ">", "")
  sequence = str_c(output[c(2:lines),], collapse = "")
  tmp = data.frame(AccessionNumProtein = accession,
                   name = name,
                   sequence = sequence)
  out = rbind(out, tmp)
}

A_ok = A[inds_okay, ]

out$AccessionNumGene = refIDs[inds_okay]
out = merge(out, A_ok)
keep = c("AccessionNumGene", "AccessionNumProtein", "name", "sequence" ,  "Species" , "Order", "Class")
out = out[, keep]
# A_ok = cbind(A_ok, out$AccessionNumProtein)
# A_ok = cbind(A_ok, out$name)
# A_ok = cbind(A_ok, out$sequence)

refIDs[inds_none]
```

    ## [1] "XM_002930611.3" "KAB0345583.1"   "VFV30336.1"     "KAF0878287.1"  
    ## [5] "AAY57872.1"     "NP_001358344.1"

``` r
none_list = c("sequence_KAB0345583.1.fasta",
              "sequence_VFV30336.1.fasta",
              "sequence_KAF0878287.1.fasta",
              "sequence_AAY57872.1.fasta",
              "sequence_NP_001358344.1.fasta")

tmp = NULL
for (a in 1:length(none_list)){
  fastaFile = readAAStringSet(none_list[a])
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df1 <- data.frame(AccessionNumProtein = refIDs[inds_none][a],name = seq_name, sequence)
  A_tmp = subset(A, AccessionNumGene == df1$AccessionNumProtein)
  df1$Species = A_tmp$Species
  df1$Order = A_tmp$Order
  df1$Class = A_tmp$Class
  df1$AccessionNumGene = "NA"
  tmp = rbind(tmp, df1)
}

out = rbind(out, tmp)
save(out, file = "out.Rdata")
```

\#\#read in a fasta file from MEROPS

``` r
rm(list = ls())
load("out.Rdata")
out$database = "GenBank"

M<- read.csv("ACE2GenBank - MEROPS.csv")
```

    ## Warning in read.table(file = file, header = header, sep = sep, quote = quote, :
    ## incomplete final line found by readTableHeader on 'ACE2GenBank - MEROPS.csv'

``` r
M_list = paste0(M$AccessionNum, "_original.fasta")

tmp = NULL
for (a in 1:length(M_list)){
  fastaFile = readAAStringSet(M_list[a])
  name = names(fastaFile)
  sequence = paste(fastaFile)
  sequence = gsub('[[:digit:]]+', '', sequence)
  sequence = gsub("[[:space:]]", "", sequence)
  df1 <- data.frame(name, sequence)
  df1$AccessionNumProtein = M$AccessionNum[a]
  df1$AccessionNumGene = "NA"
  df1$Species = M$Species[a]
  df1$Order = M$Order[a]
  df1$Class = M$Class[a]
  df1$database = "MEROPS"
  tmp = rbind(tmp, df1)
}

out = rbind(out, tmp)

fish = c("Actinopterygii", "Sarcopterygii", "Holocephali", "Cephalaspidomorphi", "Elasmobranchii")
is_fish = rep(0, dim(out)[1])

fish_inds = which(out$Class %in% fish)
is_fish[fish_inds]=1
out$is_fish = is_fish

ACE2_sequences = out
ACE2_sequences$nchar = nchar(ACE2_sequences$sequence)
#remove any sequences that are shorter than 600 amino acids; nope, leaving them all in
# ACE2_sequences = subset(ACE2_sequences, nchar >= 600)

write.csv(ACE2_sequences, file = "ACE2_sequences.csv")
DF = ACE2_sequences[, c("name", "sequence")]
names(DF)[names(DF)=="name"]="seq.name"
names(DF)[names(DF)=="sequence"]="seq.text"

dat2fasta(DF, outfile = "ACE2.fasta")
```

    ## ACE2.fasta has been saved to  /Users/fischhoff/ilya documents/R/COVID19/ACE2
