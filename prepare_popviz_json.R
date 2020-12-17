# prepare viz json files

# clean environment
rm(list=ls())

# create viz_output_scaled.txt for Vassilios
library(stringr)
library(jsonlite)
library(plyr)

samples <- c('11116', '10960', '10962', '10961', '10643', '10590', '10591', '8967', '9052', '9053', '8968', '9508', '9062', '9063', '9064', '8987', '9507')
sample.paths <- list()
sample.paths['11116'] <- "Conde_11116/"
sample.paths['10960'] <- "Conde_10960/"
sample.paths['10961'] <- "Conde_10961/"
sample.paths['10962'] <- "Conde_10962/"
sample.paths['10643'] <- "Conde_10643/"
sample.paths['10590'] <- "Conde_10590/"
sample.paths['10591'] <- "Conde_10591/"
sample.paths['8967'] <- "Conde_8967-68/8967/"
sample.paths['9052'] <- "Conde_9052/"
sample.paths['9053'] <- "Conde_9053/"
sample.paths['9062'] <- "Conde_9062/"
sample.paths['9063'] <- "Conde_9063/"
sample.paths['9064'] <- "Conde_9064/"
sample.paths['8968'] <- "Conde_8967-68/8968/"
sample.paths['8987'] <- "Conde_8987/"
sample.paths['9507'] <- "Conde_9507/"
sample.paths['9508'] <- "Conde_9508/"

# data paths
my.date <- "200707"
json.filename <- paste0(my.date, "_slice_silac_viz")
json.path <- "output/"
main.outputs.path <- "norm_res/"
main.data.path <- "data/"
uniprot.file <- "/Users/admin/Work/PAF/projects/pumba/data/fasta/20170707_UP000005640_9606.csv"
pg.filename <-  "output_scaled.txt" 
cleavage.path <- "Human_Protein_Cleavage_data_190521.csv"


##########################################################
# load data
##########################################################

# load all data
pg.list <- list()
proteinAC.vec <- c()
pep.list <- list()
fit.list <- list()

for(sample.name in samples){
  # load proteinGroup table
  i <- which(sample.name == samples)
  my.sample.path <- paste0(main.outputs.path, sample.paths[sample.name])
  
  my.pg.path <- paste0(my.sample.path, pg.filename)
  my.pg.missing.path <- paste0(my.sample.path, pg.missing.filename)
  
  if(!file.exists(my.pg.path)){
    warning(paste0("File [", my.pg.path, "] does not exist."))
    next
  }
  
  pg <- read.table(my.pg.path, stringsAsFactors=FALSE,quote="\"", row.names=NULL,
                   header=TRUE, sep="\t", fill=TRUE, na.strings=c("Non Num\303\251rique"))
  
  
  pg.headers <- c("Protein.IDs", "Majority.protein.IDs", "Peptide.counts..all.", "Peptide.counts..razor.unique.", "Peptide.counts..unique.", "Protein.names", "Gene.names", "Fasta.headers", "Number.of.proteins", "Peptides", "Sequence.coverage....", "Mol..weight..kDa.", "Sequence.length", "Peptide.IDs")
  pg <- pg[,pg.headers]
  pg$proteinAC <- unlist(lapply(pg[,"Majority.protein.IDs"], function(x) strsplit(x, ";")[[1]][1]))
  
  proteinAC.vec <- unique(c(proteinAC.vec, pg$proteinAC))
  
  pg.list[[sample.name]] <- pg
  
  # load weight-slice fitting
  my.fit.path <- paste0(my.sample.path, "outputFittingFunction.RData")
  load(my.fit.path)
  my.fit.m <- as.matrix(predicted_masses)
  my.fit <- data.frame(slice=rownames(my.fit.m), theoMW=my.fit.m)
  fit.list[[sample.name]] <- my.fit
  
  # load peptides table
  my.pep.path <- paste0(main.data.path, sample.paths[sample.name], "peptides.txt")
  pep <- read.table(my.pep.path, stringsAsFactors=FALSE,quote="\"", row.names=NULL,
                    header=TRUE, sep="\t", fill=TRUE, na.strings=c("Non Num\303\251rique"))
  
  ratio.cols <- grep("Ratio.H.L.normalized.", colnames(pep))
  count.cols <- grep("Ratio.H.L.count.", colnames(pep))
  
  #new.ratio.cols <- paste0(colnames(pep)[ratio.cols], "_", my.fit$theoMW)
  
  pep.headers <- c("id", "Sequence", "Start.position", "End.position", "Amino.acid.before", "Amino.acid.after", colnames(pep)[ratio.cols], colnames(pep)[count.cols])
  pep.list[[sample.name]] <- pep[,pep.headers]
  
  rm(list="pep")
  rm(list="my.fit")
}

# load uniprot sequences
#uniprot <- read.table(uniprot.file, header=TRUE)
uniprot <- read.csv(uniprot.file, header=TRUE)

# load the cleavage info
cleavage <- read.csv(file=cleavage.path)
# replace yes/no by TRUE/FALSE
cleavage$Prediction <- mapvalues(cleavage$Prediction, from=c("yes", "no"), to=c(TRUE, FALSE))
cleavage$Signal.peptide <- mapvalues(cleavage$Signal.peptide, from=c("yes", "no"), to=c(TRUE, FALSE))
cleavage$Confirmed <- mapvalues(cleavage$Confirmed, from=c("yes", "no"), to=c(TRUE, FALSE))

##########################################################
# create JSON
##########################################################

# search term
search.terms <- data.frame(proteinAC=character(), searchTerm=character(), stringsAsFactors = FALSE)

# group them by first proteinAC
protein.list <- list()
protein.list.idx <- 1

# create a new json every "flush.size" proteins
flush.size <- 100
flush.idx <- 1

# proteinAC <- "P02786"
# proteinAC <- "P07339"
# proteinAC <- "A0A140TA64"
# proteinAC <- "P16104"
# proteinAC <- "Q9ULP9"
# proteinAC <- "P62937"
# proteinAC <- "P62987"

tot.loops <- length(proteinAC.vec)
current.loop <- 1

for(proteinAC in proteinAC.vec){
  print(paste0(current.loop, " - ", tot.loops))
  current.loop <- current.loop + 1
  
  protein.info <- list()
  protein.info[['proteinAC']] <- unbox(proteinAC)
  protein.info[['samples']] <- list()
  
  # add the cleavage table
  cleavage.ids <- which(cleavage$Uniprot == proteinAC)
  if(length(cleavage.ids > 0)){
    cleavage.list <- list()
    cleavage.list.id <- 1
    
    for(cleavage.id in cleavage.ids){
      cleavage.info <- list()
      cleavage.info[['pos']] <- unbox(as.numeric(cleavage$P1[cleavage.id]))
      cleavage.info[['seq']] <- unbox(as.character(cleavage$P5.P5.[cleavage.id]))
      cleavage.info[['probScore']] <- unbox(as.numeric(cleavage$Prob_score[cleavage.id]))
      cleavage.info[['prediction']] <- unbox(as.logical(cleavage$Prediction[cleavage.id]))
      cleavage.info[['signalPep']] <- unbox(as.logical(cleavage$Signal.peptide[cleavage.id]))
      cleavage.info[['confirmed']] <- unbox(as.logical(cleavage$Confirmed[cleavage.id]))
      cleavage.list[[cleavage.list.id]] <- cleavage.info
      cleavage.list.id <- cleavage.list.id + 1
    }

    protein.info[['cleavages']] <- cleavage.list
  }
  
  # min and max molweight
  min.molweight <- NA
  max.molweight <- NA  
  
  # sequence info
  seq <- as.character(uniprot$seq[uniprot$ac == proteinAC])
  seq.length <- nchar(seq)
  
  # ignore the entry if we didnt find it
  if(length(seq.length) < 1) next
  
  #protein.info[['sequence']] <- unbox(seq)
  #protein.info[['sequenceLength']] <- unbox(seq.length)
  
  all.protein.names <- c()
  all.gene.name.list <- list()
  all.description.list <- list()
  pep.info.list <- list()
  pep.info.list.idx <- 1
  
  # look for the protein in all samples
  for(sample.name in samples){
    pg <- pg.list[[sample.name]]
    #pg.idx <- which(pg$proteinAC == proteinAC)
    pg.idx <- grep(proteinAC, pg$Majority.protein.IDs)
    
    # skip the sample if the protein is not found (or to many)
    if(length(pg.idx) < 1){
      next
    } 
    if(length(pg.idx) > 1) stop(paste0("there are more than one entry for this proteinAC [" , proteinAC, "] in [", sample.name, "]"))
    
    sel.pg <- pg[pg.idx,]
    
    # load peptides
    pep.ids <- strsplit(sel.pg$Peptide.IDs, ";")[[1]]
    peps <- pep.list[[sample.name]]
    sel.peps <- peps[peps$id %in% pep.ids,]
    
    # we will need the slice-mw fit
    slice.mw.fit <- fit.list[[sample.name]]
    
    pep.seq.list <- list()
    pep.seq.list.idx <- 1
    
    for(sel.pep.idx in 1:nrow(sel.peps)){  
      p <- sel.peps[sel.pep.idx,]
      start.pos <- as.numeric(p$Start.position)
      end.pos <- as.numeric(p$End.position)
      pep.seq <- substr(seq, start.pos, end.pos)
      
      # we ignore wrong sequences
      if(p$Sequence != pep.seq){
        #warning(paste0("[", proteinAC, "] peptideID [", p$id, "] from sample [", sample.name, "] is [", pep.seq, "] instead of expected [", p$Sequence, "]"))
        start.pos <- NA
        end.pos <- NA
      }
      
      ratio.cols <- grep("Ratio.H.L.normalized.", colnames(p))
      count.cols <- grep("Ratio.H.L.count.", colnames(p))
      
      nr.of.valid.peps <- 0
      
      # create a peptide per slice
      for(col.idx in 1:length(ratio.cols)){
        pep.ratio <- p[,ratio.cols[col.idx]]
        
        # ignore NA
        if(is.na(pep.ratio)) next
        nr.of.valid.peps <- nr.of.valid.peps + 1
        
        mol.weigth <- slice.mw.fit$theoMW[slice.mw.fit$slice == col.idx]
        
        pep.info <- list()
        pep.count <- p[,count.cols[col.idx]]
        pep.info[['sampleName']] <- unbox(sample.name)
        pep.info[['log2ratio']] <- unbox(log2(pep.ratio))
        pep.info[['ratioCount']] <- unbox(pep.count)
        pep.info[['molWeight']] <- unbox(mol.weigth)
        pep.info[['sliceNr']] <- unbox(col.idx)
        pep.info[['sequence']] <- unbox(p$Sequence)
        pep.info[['id']] <- unbox(paste(sample.name, p$id, mol.weigth, sep=":"))
        
        if(! is.na(start.pos)){
          pep.info[['aminoAcidBefore']] <- unbox(p$Amino.acid.before)
          pep.info[['aminoAcidAfter']] <- unbox(p$Amino.acid.after)
          pep.info[['endPos']] <- unbox(end.pos)
          pep.info[['startPos']] <- unbox(start.pos)
        }
        
        pep.info[['ratio']] <- unbox(pep.ratio)
        min.molweight <- min(min.molweight, mol.weigth, na.rm=TRUE)
        max.molweight <- max(max.molweight, mol.weigth, na.rm=TRUE)
        
        pep.info.list[[pep.info.list.idx]] <- pep.info
        pep.info.list.idx <- pep.info.list.idx + 1
      }
      
      # skip this peptide if it does not have a single valid psm
      if(nr.of.valid.peps == 0) next
      
      pep.seq <- list()
      pep.seq[['sequence']] <- unbox(p$Sequence)
      if(! is.na(start.pos)){
        pep.seq[['startPos']] <- unbox(start.pos)
        pep.seq[['endPos']] <- unbox(end.pos)
      }
      
      # only add them if it does not exist yet
      if(! any(unlist(lapply(pep.seq.list, function(x) x$sequence == pep.seq[['sequence']])))){
        pep.seq.list[[pep.seq.list.idx]] <- pep.seq
        pep.seq.list.idx <- pep.seq.list.idx + 1
      }
    }
    
    sample.info <- list()
    sample.info[['peptideCountsAll']] <- unbox(sel.pg$Peptide.counts..all.)
    sample.info[['peptideCountsRazorUnique']] <- unbox(sel.pg$Peptide.counts..razor.unique.)
    sample.info[['peptideCountsUnique']] <- unbox(sel.pg$Peptide.counts..unique.)
    sample.info[['nrPeptides']] <- unbox(sel.pg$Peptides)
    sample.info[['seqCoverage']] <- unbox(sel.pg$Sequence.coverage....)
    sample.info[['proteinIDs']] <- strsplit(sel.pg$Majority.protein.IDs, ";")[[1]]
    sample.info[['peptideSequences']] <- pep.seq.list
    protein.info[['samples']][[sample.name]] <- sample.info
    
    # get the list of alternative proteins and gene names
    protein.names <- strsplit(sel.pg$Majority.protein.IDs, ";")[[1]]
    gene.names <- strsplit(sel.pg$Gene.names, ";")[[1]]
    protein.desc <- strsplit(sel.pg$Protein.names, ";")[[1]]
    
    # put entries to all gene.names.list
    for(pn.idx in 1:length(protein.names)){
      all.gene.name.list[[protein.names[pn.idx]]] <- gene.names[pn.idx]
      all.description.list[[protein.names[pn.idx]]] <- protein.desc[pn.idx]
    }
    
    # add to search.terms
    for(protein.idx in 1:length(protein.names)){
      protein.name <- protein.names[protein.idx]
      # check if it already exists
      if(length(which(search.terms$proteinAC == protein.name)) == 0){
        new.search.term <- c(protein.name)
        if(! is.na(gene.names[protein.idx])) new.search.term <- c(new.search.term, gene.names[protein.idx])
        if(! is.na(protein.desc[protein.idx])) new.search.term <- c(new.search.term, protein.desc[protein.idx])
        new.search.entry <- data.frame(proteinAC = protein.name, searchTerm = paste(new.search.term, collapse = " | "), stringsAsFactors = FALSE)
        search.terms <- rbind(search.terms, new.search.entry)
      }else{
        sel.search.idx <- which(search.terms$proteinAC == protein.name)
        sel.search.term <- search.terms[sel.search.idx,]
        if(sel.search.term$proteinAC == sel.search.term$searchTerm){
          new.search.term <- c(protein.name)
          if(! is.na(gene.names[protein.idx])) new.search.term <- c(new.search.term, gene.names[protein.idx])
          if(! is.na(protein.desc[protein.idx])) new.search.term <- c(new.search.term, protein.desc[protein.idx])
          search.terms$searchTerm[sel.search.idx] <- paste(new.search.term, collapse = " | ")
        }
      }
    }
    
    all.protein.names <- c(all.protein.names, protein.names)
  }
  
  all.protein.names <- unique(all.protein.names)
  all.gene.names <- unlist(lapply(all.protein.names, function(apn){ all.gene.name.list[[apn]] }))
  all.descriptions <- unlist(lapply(all.protein.names, function(apn){ 
      my.desc <- all.description.list[[apn]] 
      if(is.na(my.desc)){
        "-"
      }else{
        my.desc
      }
    }))
  
  protein.info[['alternativeProteinACs']] <- all.protein.names
  protein.info[['alternativeGeneNames']] <- all.gene.names
  protein.info[['geneName']] <- unbox(all.gene.names[1])
  protein.info[['theoMolWeight']] <- unbox(sel.pg$Mol..weight..kDa.)
  protein.info[['theoMolWeightLog10']] <- unbox(log10(sel.pg$Mol..weight..kDa.))
  protein.info[['fastaHeaders']] <- unbox(sel.pg$Fasta.headers)
  protein.info[['description']] <- unbox(paste(all.descriptions, collapse="; "))
  protein.info[['minMolWeight']] <- unbox(min.molweight)
  protein.info[['maxMolWeight']] <- unbox(max.molweight)
  protein.info[['peptides']] <- pep.info.list
  
  protein.list[[protein.list.idx]] <- protein.info
  protein.list.idx <- protein.list.idx + 1
  
  # write a json and create reinitialize protein.list
  if(protein.list.idx > flush.size){
    cat(toJSON(protein.list), file=paste0(json.path, json.filename, "_", flush.idx, ".json"))
    protein.list <- list()
    protein.list.idx <- 1
    flush.idx <- flush.idx + 1
  }
}
 # write the last json
cat(toJSON(protein.list), file=paste0(json.path, json.filename, "_", flush.idx, ".json"))

# write the search term json
cat(toJSON(search.terms), file=paste0(json.path, my.date, "_search_terms.json"))

# write sequences to json
#cat(toJSON(uniprot), file=paste0(json.path, json.filename, "_uniprot.json"))
