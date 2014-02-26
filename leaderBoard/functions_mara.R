
## The Broad Institute of MIT and Harvard / Cancer program.
# General utility functions
#

suppressPackageStartupMessages(library(multicore))
suppressPackageStartupMessages(library(gplots))

interval_subset <- function(af, int_file) {
#INTERVAL_SUBSET: Function to subset a File by specific interval. Similar to bedtools intersect which would be implemented via command line
#af: the original file you want to subset. It must contain Chromosome, Start_position, End_position but may contain anything else as well
#int_file: interval file to intersect with. Must contain Chromosome, Start_position, End_position.

print("loading in intervals package")
library(intervals)

chrom = append(seq(1,22), c("X","Y"))
af_subset = data.frame()

for (ii in seq(1,length(chrom))) {
	chrom_num = chrom[ii]
	maf.split = subset(af, Chromosome==chrom_num) #af is the maf file to reduce to the specific interval. 
	if (dim(maf.split)[1]==0){
		print(paste("No intervals for this chromosome in the original maf:", chrom_num))
	} else {
		maf.split.int = Intervals(maf.split[,c("Start_position", "End_position")], close = c(TRUE, TRUE), type = "Z")
		bed.subset = subset(int_file, Chromosome==chrom_num)
	} 
	if (dim(bed.subset)[1]==0){
		print(paste("No intervals for this chromosome:", chrom_num))
	} else {
		bed.subset.int = Intervals(bed.subset[,c("Start_position", "End_position")], close = c(TRUE,TRUE), type = "Z")
		overlap = unique(unlist(interval_overlap(bed.subset.int, maf.split.int)))
		overlap_complete = overlap[overlap > 0]
		maf.split.final = maf.split[overlap_complete,]
		af_subset = rbind(af_subset, maf.split.final) #final maf that only contains overlapping intervals
	}
}
reduced_file = af_subset
return(reduced_file)
}


#Jaccard clustering, comes from package prabclus
library(prabclus)
jaccard_missval <- function(regmat) 
{
    nart <- ncol(regmat)
    jdist <- rep(0, nart * nart)
    dim(jdist) <- c(nart, nart)
    reg.col.sum <- apply(regmat, 2, sum, na.rm=TRUE)
    reg.aggrement <- t(regmat) %*% regmat
    jdist <- 1 - reg.aggrement/(reg.col.sum - t(t(reg.aggrement) - 
        reg.col.sum))
    jdist
}



# used for genome plot
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

## shorthand listing largest objects in the workspace
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

## lapply(list, dim) shortcut
ldim = function(l)
  {
    return(lapply(l, function(x) {if (!is.null(dim(x))) dim(x) else length(x)}))
  }

############################
# qq_chisq 
#
# plots qq plot for observed chi squared vector obs against expected with df = df
#############################
qq_chisq = function(obs, df=1, highlight = c(), hexbins = NULL, samp = NULL)
  {
    if (!is.null(hexbins))
      require(hexbin)
    
    obs = obs[!is.na(obs)]

    if (!is.null(samp))
      if (samp<length(obs))
        obs = sample(obs, samp)
    
    exp = rchisq(length(obs), df = df)
    ord = order(obs)
    smax = max(c(exp, obs))
    colors = vector(mode = "character", length = length(obs)); colors[] = "black"; colors[highlight] = "red";
    colors = factor(colors[ord]);

    if (!is.null(hexbins))
      plot(hexbin(sort(exp), obs[ord], xbins = hexbins,  xlab = "Expected Chi^2", ylab = "Observed Chi^2"))
    else
      {
        plot(sort(exp), obs[ord], xlab = "Expected Chi^2", ylab = "Observed Chi^2", xlim = c(0, smax), ylim = c(0, smax), col = colors, pch = 18);
        lines(x=c(0, smax), y = c(0, smax), col = "gray");
      }
  }

#########
# get_filename
#
# grabs filenames from list of paths
########
get_filename = function(paths)
  {
    return(gsub('(^|(.*\\/))?([^\\/]+)', '\\3', paths))
  }

########
# splits a single string according to fixed widths contained in fw (ie each components i of fw denotes the width of field i in string str
########
strsplit.fwf = function(str, fw)
  {
    if (length(str)>1)
      {
        warning('String should be of length 1, only taking first element')
        str = str[1];
      }
    
    cs = cumsum(fw);
    return(substr(rep(str, length(fw)), cs-fw+1,c(cs[1:(length(cs)-1)], nchar(str))))
  }

########
# Returns vector of line counts for each file in path 
########
line.counts = function(paths)
  {
    out = rep(NA, length(paths))
    ix = which(file.exists(paths))
    out[ix] = sapply(paths, function(x) { p = pipe(paste('cat ', x, ' | wc -l ')); as.numeric(readLines(p)); close(p)});
    return(out)
  }


##########
# pk
#
# "pk"'s at vector or data frame ie samples n rows from it with replacement (or max rows / items if less than n)#
##########
pk = function(x, n = 0)
  {
    if (inherits(x, 'data.frame'))
      {
        n = min(nrow(x), n)
        return(x[sample(1:nrow(x), n), ])        
      }
    else
      {
        n = min(length(x), n)
        return(sample(x, n))
      }      
  }


##############
# border - orders rows of a logical / binary matrix treating each row as binary number with digits encoded as TRUE / FALSE values of entries
#
##############
border = function(B, na.rm = TRUE)
  {
    B = array(as.logical(B), dim = dim(B))
    tmp = vector(mode = "numeric", length = nrow(B));
    if (na.rm)
      B[is.na(B)] = FALSE;
    for (i in 1:ncol(B))
        tmp = tmp + 2^(ncol(B)-i)*as.numeric(B[,i]==1);
    return(order(tmp))
  }

#########
# Returns vector of column names for list of files, using sep as delimiter
#########
column.names = function(paths, sep = "\t")
  {
    out = list();
    out[paths] = NA;
    ix = file.exists(paths);
    out[ix] = strsplit(sapply(paths[ix], function(x) readLines(x, 1)), sep)
    return(out)
  }

############################
# qq_pval
#
# plots qq plot for observed pval vs uniform
#############################
qq_pval = function(obs, highlight = c(), hexbins = NULL, samp = NULL, lwd = 1, bestfit=T, color='black', input.pch=18, input.cex=1, conf.lines=T, input.MAX=NULL, qvalues=NULL, genes=NULL, ...)
{	
	obs = -log10(obs[!is.na(obs)])
	obs = obs[!is.infinite(obs)]
	
	if (!is.null(samp))
		if (samp<length(obs))
			obs = sample(obs, samp)
	
	N <- length(obs)
	## create the null distribution 
	## (-log10 of the uniform)
	exp <- -log(1:N/N,10)
	
	if (is.null(input.MAX))
		MAX <- max(obs,exp) + 0.5
	else
		MAX <- input.MAX
	
	c95 <- rep(0,N)
	c05 <- rep(0,N)
	
	for(i in 1:N){
		c95[i] <- qbeta(0.95,i,N-i+1)
		c05[i] <- qbeta(0.05,i,N-i+1)
	}
	
	if (conf.lines){
		## plot the two confidence lines
		plot(exp, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l", axes=FALSE, xlab="", ylab="")
		par(new=T)
		plot(exp, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l", axes=FALSE, xlab="", ylab="")
		par(new=T)
		p1 <- c(exp[1], exp[1])
		p2 <- c(-log(c95,10)[1], -log(c05,10)[1])
		lines(x=p1, y=p2)
		x.coords <- c(exp,rev(exp))
		y.coords <- c(-log(c95,10),rev(-log(c05,10)))
		polygon(x.coords, y.coords, col='light gray', border=NA)
		par(new=T)
	}
	
	ord = order(obs)
	
	colors = vector(mode = "character", length = length(obs)); colors[] = "black"; colors[highlight] = "red";
	colors = factor(colors[ord]);
	
        dat = data.frame(x = sort(exp), y = obs[ord]);
        plot(dat$x, dat$y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, MAX), ylim = c(0, MAX), pch=input.pch, cex=input.cex, bg=color, ...);
        if (!is.null(qvalues)){
          genes.to.label <- which(qvalues[ord] <= 0.2);
          genes <- genes[ord];
          if (length(genes.to.label > 0)){
            text(dat$x[genes.to.label], dat$y[genes.to.label], labels=genes[genes.to.label], pos=3);
          }
        }
        lines(x=c(0, MAX), y = c(0, MAX), col = "black", lwd = lwd);
        lambda = lm(y ~ x-1, dat)$coefficients;
        if (bestfit)
          {
            lines(x=c(0, MAX), y = c(0, lambda*MAX), col = "red", lty = 2, lwd = lwd);
            legend('bottomright',sprintf('lambda = %.2f', lambda), text.col='red', bty='n')
          }        
      }


##############################
# wfplot
#
# Quick waterfall plot
#
# data is a numeric vector
# labels are text labels of the same length as data
# col is either (1) an unamed list of unique colors (2) a named list mapping unique labels to colors
##############################
wfplot = function(data, labels = NULL, col = NULL, leg.pos = NULL, ...)
  {
    ix = order(data);
    ulab = unique(labels)

    if (!is.null(col))
      {
        if (is.null(names(col)))
          {
            if (length(ulab)>2)
              col = brewer.pal(length(ulab), 'Set3')
            else
              col = c('gray', 'red')

            col = col[match(labels, ulab)]
            names(col) = labels;
          }
        else
          {
            og.col = col;
            col = col[labels];
            ulab = intersect(names(og.col), names(col))
          }
      }
    
    barplot(data[ix], col = col[ix], border = FALSE, ...)

    if (is.null(leg.pos))
      leg.pos = c(mean(par('usr')[1:2])/4, 3*data[which.max(abs(data))]/4)
        
    legend(leg.pos[1], leg.pos[2], legend = ulab, fill = col[ulab])
  }
               
               


############
# Takes a character vector and makes an expression for creating that character vector
############
list.expr = function(x)
  {
    paste("c('", paste(x, sep = "", collapse = "', '"), "')", sep = "")
  }

########
# toggles options error recover / NULL
########
fuckr = function()
  {
    if (!is.null(options()$error))
      {
        options(error = NULL);
        print('Options error set to NULL');
      }
    else
      {
        options(error = recover);
        print('Options error set to recover');
      }
  }

#################
## flatten
##
## flattens 3rd dim of 3D array along cdim 1 (ie rows) or cdim 2 (ie cols) pasting together the appropriate combinations of dimnames with sep "sep"
## or if sep = NULL, then just dropping the 3rd dimension names
##
#################
flatten = function(A, cdim = 2, sep = "_")
{
  if (!(cdim==1 | cdim ==2))
    stop('cdim must be 1 or 2')

  ind = order(rep(c(1:dim(A)[cdim]), dim(A)[3]));
  
  out = A[,,1];

  if (cdim == 2)
    {
      if (dim(A)[3]>1)
        for (i in 2:dim(A)[3])
          out = cbind(out, A[,,i]);            
      dimnames(A)[[1]] = dimnames(A)[[1]]
    }
      
  if (cdim == 1)
    {
      if (dim(A)[3]>1)
        for (i in 2:dim(A)[3])
          out = rbind(out, A[,,i]);            
      dimnames(A)[[2]] = dimnames(A)[[2]]
    }

  out = out[,ind]; #reshuffle to get desired ordering  
  newdimnames = rep(dimnames(A)[[cdim]], each = dim(A)[3]);
  if (!is.null(sep))
    newdimnames = paste(newdimnames, dimnames(A)[[3]], sep = sep);  
  dimnames(out)[[cdim]] = newdimnames;

  return(out)
}

##################
# Makes bsub command that wraps shell command "cmd" to send to queue "queue"
# redirecting output / error etc streams to path prefixed by "jname",
# optional_args: maximum memory requirements "mem", "jlabel" job label
##################
bsub_cmd = function(cmd, queue, jname, jlabel=NULL, jgroup = NULL, mem=NULL, group = "cgafolk", cwd = NULL)
  {
    qjname = paste( "\"", jname, "\"", sep="" )
    qjout = paste( "\"", jname, ".bsub.out", "\" ", sep="" )
    qjerr = paste( "\"", jname, ".bsub.err", "\" ", sep="" )
    qjrout = paste( "\"", jname, ".R.out", "\" ", sep="" )
    out_cmd = paste( "bsub -q ", queue, " -o ", qjout, " -e ",  qjerr, " -P ", group);
    if (!is.null(mem)) out_cmd = paste(out_cmd, " -R \"rusage[mem=", mem, "]\" ", sep = "");
    if (!is.null(jlabel)) out_cmd = paste(out_cmd, " -J ", jlabel )
    if (!is.null(jgroup)) out_cmd = paste(out_cmd, " -g ", sub('^\\/*', '/', jgroup))
    if (!is.null(cwd)) out_cmd = paste(out_cmd, " -cwd ", cwd )
    out_cmd = paste(out_cmd," \"",  cmd, "\"", sep = "")
  }


##################
# chunk_and_dump
#
# Takes data frame and chunks into either "k" chunks or so that each chunk has <= "nrow" rows
# then dumps each chunk to file, returning a vector of resulting file paths
#
# If nrow is specified then "k" is ignored. 
##################
chunk_and_dump = function(df, file.prefix = 'chunk', k = 10, nrow = NULL, sep = "\t", quote = F, row.names = F,
  no.dump = F, # no.dump if you want just filenames
  ...) ## ... passed to write.table
  {
    if (!no.dump)
      system(paste('mkdir -p', gsub('\\/([^\\/]+)$', '\\/', file.prefix)))

    if (!is.null(nrow))
      ix = seq(1, nrow(df), nrow)
    else
      ix = round(seq(1, nrow(df), nrow(df)/k))

    if (ix[length(ix)] != nrow(df))
      ix = c(ix, nrow(df)+1)

    fn = paste(file.prefix, ix[1:(length(ix)-1)], sep = ".")

    if (!no.dump)
      {
        if (length(ix)==1)
          write.table(df[1, ], file = fn, sep = sep, quote = quote, row.names = F, ...)
        else
          sapply(1:(length(ix)-1), function(x) write.table(df[ix[x]:(ix[x+1]-1), ], file = fn[x], sep = sep, quote = quote, row.names = F, ...))
      }
    
    return(fn)
  }
      


##################################
# write.bed
#
# dumps segs df with fields $chr, $start, $end, $label to a bed file with filename fn
# also handles segs with chromStart, chromEnd, startbp / endbp, pos1 / pos2 nomenclature
################################## 
write.bed = function(segs, fn, chr = FALSE, header = T, uncollapse = F)
  {

    if (header)
      writeLines(paste("track name = ", fn, "color = 0,0,255"), con = fn)
    else
      system(paste('rm -f', fn))

    segs = standardize_segs(segs);

    if (is.null(segs$label))
      segs$label= paste('seg', 1:nrow(segs), sep = "");

    if (chr)
      segs$chr = paste('chr', gsub('chr', '', segs$chr), sep = "");

    if (uncollapse)
      segs$pos2[segs$pos1==segs$pos2] = segs$pos2[segs$pos1==segs$pos2] + 1;
    
    write.table(segs[, c('chr', 'pos1', 'pos2', 'label')], file = fn, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE);        
  }

read.bed = function(bed_fn)
{
  bed = read.delim(bed_fn, skip = 1, strings = F, header = FALSE);  
  colnames(bed)[1:4] = c('chr', 'start', 'end', 'name');  
  bed$chr = gsub('chr', '', bed$chr);
  return(bed)
}

read.gff = function(gff_fn)
{
  gff = read.delim(gff_fn, skip = 2, strings = F, header = FALSE);  
  colnames(gff) = c('chr', 'name', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'group');  
  gff$chr = gsub('chr', '', gff$chr);
  return(gff)
}



##################
## Produces (simple) R code calling function named "func" with args in list "argv", prepending with
## source() call to directories in the vector "sources" if specified.
##
## NOTE: args in argv can be lists or vectors consisting of numerical values or characters.  Lists can have named fields.  
## These will be assigned in a "hard coded" way in the Rcode, so these should be ideally scalars or
## pretty short vectors / lists. 
##
## For code to run properly, the names of "argv" must correspond to argument names of "func", or
## if the list has unnamed fields then they must be ordered in the order of the function args.  
##
## Useful for dumping tmp code files for farming when there are many arguments being passed
##################
func_code = function(func, argv, sources = c())
  {
    out = "";
    if (length(sources)>0)
      {
        out = sprintf('%s%s\n', out, paste("source(\"", sources, "\")", sep = "", collapse = "\n"));
      }

    argv_strings = vector(mode="character", length = length(argv));

    for (i in 1:length(argv))
      {
        this_arg = eval(argv[[i]]); # need to eval if data frame slice passed down as vector (i.e. as "call")
        
        if (is.list(this_arg) & is.null(dim(this_arg))) # checks we have a bona fide list ie not a data frame
          {
            if (max(unlist(lapply(this_arg, length)))>1)
            {
              print("Error: nested list arguments not allowed in argv");
              return( NA );
            }

            list_strings = as.vector(this_arg)
            
            chars = unlist(lapply(this_arg, is.character));  # put quotes around char list items
            list_strings[chars] = paste('\"', list_strings[chars], '\"', sep = "");
            
            if (!is.null(names(this_arg))) # take care of named list items if exists
              {
                named = names(this_arg) != "";        
                list_strings[named] = paste(names(this_arg)[named], " = ", list_strings[named],  sep = "");  # prepend "name=" to named items of list
              }
            
            argv_strings[[i]] = sprintf("list(%s)", paste(list_strings, collapse = ", "));  # pre-pend list constructor and comma concat list items            
          }
        else if (is.vector(this_arg) & is.null(dim(this_arg))) # make sure we have vector and not an array
          {
            vec_strings = this_arg;
            if (is.character(this_arg))
              vec_strings = paste('\"', vec_strings, '\"', sep = "");
            
            if (length(vec_strings)>1) # use c() if we have a vector greater than length 1
              argv_strings[i] = sprintf("c(%s)", paste(vec_strings, collapse = ","))
            else
              argv_strings[i] = vec_strings;
          }
        else if (is.null(this_arg))
          argv_strings[i] = 'NULL'            
        else          
          {
            print("Error: unsupported data type in argv");
            return( NA );
          }
      }
    
    if (!is.null(names(argv))) # take care of named args if exist
      {
        named = names(argv) != "";  
        argv_strings[named] = paste(names(argv)[named], " = ", argv_strings[named], sep = ""); 
      }

    out = sprintf('%s\n%s(%s)\n', out, func, paste(argv_strings, collapse = ",\n "));
    out
  }

#######################
# first.bigger
# 
# Given two <ordered> vectors a and b, returns vector v such that v[i] = index of first element in b that is greater than a[i]
#
#######################
first.bigger = function(a, b)
  {
    j = 0;
    out = vector(mode = "numeric", length = length(a))+NA;

    for (i in 1:length(a))
      {
        while (j<length(b) & b[j+1] <= a[i]) j=j+1;
        if (j < length(b)) out[i] = j+1;
      }
    return(out)
  }

uind = function(list)
  {
    uel = unique(list)
    sapply(uel, function(x){list(which(list==x))})
  }

rand = function(x, y=1)
{
  if (is.numeric(x))
    {
      vdim = x*y;
      array(runif(vdim), dim = c(x, y))	
    }
  else
    {
      ## is a matrix
      array(runif(prod(dim(x))), dim = dim(x))
    }
}

#############################
# readAmpDel = function(fn)
#
# reads GISTIC peak output
############################
readAmpDel = function(fn){
	lines <- readLines(fn, n = 10)
	items = noquote(strsplit(lines[grep("wide peak boundaries", lines)], '\t')[[1]])
	parsed = lapply(items[2:length(items)], function(x){ gsub('^chr(..?):([0-9]+)-([0-9]+)$', '\\1\t\\2\t\\3', x, perl = TRUE)})
	parsed2 = rapply(parsed[lapply(parsed, length)>0] , function(x){strsplit(x, '\t')}[[1]]);
	CNAS = data.frame(t(matrix(parsed2, ncol = length(parsed2)/3, nrow = 3)), stringsAsFactors = FALSE)
	names(CNAS) = c('chr', 'startbp', 'endbp')
	CNAS$startbp = as.numeric(CNAS$startbp)
	CNAS$endbp = as.numeric(CNAS$endbp)
	items = strsplit(lines[grep("cytoband", lines)], '\t')[[1]]
        CNAS$name = items[2:length(items)];
        items = strsplit(lines[grep("q value", lines)], '\t')[[1]]
        CNAS$q.val = as.numeric(items[2:length(items)]);
        items = strsplit(lines[grep("residual q value", lines)], '\t')[[1]]
        CNAS$res.q.val = as.numeric(items[2:length(itemips)]);        
	CNAS
}

#############################
# rrbind = function(df1, df2, [df3 ... etc], )
#
# like rbind, but takes the intersecting columns of the dfs
#
# if union flag is used then will take union of columns (and put NA's for columns of df1 not in df2 and vice versa)
############################
rrbind = function(..., union = T)
  {
    dfs = list(...);  # gets list of data frames
    names.list = lapply(dfs, names);
    cols = unique(unlist(names.list));
    unshared = lapply(names.list, function(x) setdiff(cols, x));
    expanded.dfs = lapply(1:length(dfs), function(x) {dfs[[x]][, unshared[[x]]] = NA; return(dfs[[x]])})
    out = do.call('rbind', expanded.dfs);

    if (!union)
      {
        shared = setdiff(cols, unique(unlist(unshared)))
        out = out[, shared];
      }
     
   return(out)
  }

#############################
# read.delim.cat(paths, ... [read.delim arguments])
#
# takes a vector of tab delimited file paths and concatenates them into a
# single data frame (takin union of identically named / numbered columns as a default)
############################
read.delim.cat = function(paths, skip = rep(0, length(paths)), cols = NULL, include.paths = T, cores = NULL, ...)	
  {
    paths[is.na(paths)] = ""; ;
    does.not.exist = !file.exists(paths);
     
    if (any(does.not.exist))
      warning(sprintf('Ignoring %s paths that do not exist on the file system.', length(which(does.not.exist))))
    paths = paths[!does.not.exist];

    if (length(skip) ==1)
      skip = rep(skip, length(paths))

    if (is.null(names(paths)))
      names(paths) = paths;
    
    # scope out files to filter out those with 0 rows and find common columns
    if (!is.null(cores))
      dfs = mclapply(1:length(paths), function(x) {tmp.df = read.delim(paths[x], skip = skip[x], ...);
    										   if (nrow(tmp.df) != 0) cbind(paths[x], tmp.df) else data.frame() }, mc.cores = cores)
    else
      dfs = lapply(1:length(paths), function(x) {tmp.df = read.delim(paths[x], skip = skip[x], ...);
    										   if (nrow(tmp.df) != 0) cbind(paths[x], tmp.df) else data.frame() })
    
    dfs = dfs[sapply(dfs, nrow)!=0];    
    
    if (length(dfs)==0)
      return(NULL);
   
    out = do.call('rrbind', dfs)

    if (!is.null(cols))
      out = cbind(out[,1], out[, cols]);

    if (include.paths)
      names(out)[1] = 'source.path'
    else
      out = out[,-1]
   
    return(out)
  }


####################
# Retrieves all snps in r2 "ldthresh" of snps in region defined by data.frame region with fields chr pos1 (and optional pos2) fields
# 
# Example:
#
# ldQuery(data.frame(chr ='chr1', pos1 = 32323232, pos2 = 32343232), 0.9);
#
#
# RMysql installation cmd
#  install.packages('RMySQL', type='source', lib = '~/R/x86_64-unknown-linux-gnu-library/2.8/', configure.args = "--with-mysql-inc='/local/util/include/mysql/' --with-mysql-lib='/local/util/lib/mysql/'")
####################
ldQuery = function(region = NULL, ldthresh)
{
  require('RMySQL')
  HOST = "vmfe7-85b"
  HOST = "guest124";
  con = dbConnect(MySQL(),
    user = "hapmap", password = "qwerasdf", dbname="hapmap", host=HOST)

  if (is.null(region$pos2))
      region$pos2 = region$pos1;

  if (!any(grep('chr', region$chr)))
    region$chr = paste('chr', region$chr, sep = "");

  q1 = sprintf('select * from ld_CEU where chr = "%s" and r2 > %f and pos1 between %s and %s', 
    region$chr, ldthresh, region$pos1, region$pos2)
  
  q2 = sprintf('select * from ld_CEU where chr = "%s" and r2 > %f and pos2 between %s and %s;', 
    region$chr, ldthresh, region$pos1, region$pos2)
  
  out = rbind(dbGetQuery(con, q1), dbGetQuery(con, q2))

  dbDisconnect(con)

  return(out[!duplicated(out),])  
}	

####################
# Retrieves all pairwise CEU ld values for snps in list snp
#
# ie ldSNP(snps)
#
####################
ldSNP = function(snps)
{
  require('RMySQL')
  con = dbConnect(MySQL(),
    user = "hapmap", password = "qwerasdf", dbname="hapmap", host="vmfe7-85b")

  q = sprintf('select * from ld_CEU where snp1 in (\"%s\") and snp2 in (\"%s\")', 
     paste(snps, collapse = '","'), paste(snps, collapse = '","'))

print(q)
  
  out = dbGetQuery(con, q);

  dbDisconnect(con)

  return(out)
}

####################
# ldSNP = function(snps, pop)
## Retrieves mafs for snp vector "snps" in HapMap pop "pop"
#
####################
af_SNP = function(snps, pop = 'CEU')
{
  require('RMySQL')
  con = dbConnect(MySQL(),
    user = "hapmap", password = "qwerasdf", dbname="hapmap", host="vmfe7-85b")

  q = sprintf('select * from af where snp in (\"%s\") and pop = \"%s\"', 
     paste(snps, collapse = '","'), pop)
  
  out = dbGetQuery(con, q);
  rownames(out) = out$snp;
  out = out[snps,];
  out$snp = snps;
  rownames(out) = NULL;
  
  dbDisconnect(con)

  return(out)
}

###############################
# fisher.plot
#
# Plots fisher contingency table 
#
###############################
fisher.plot = function(O)
  {
    require(plotrix)
    fish = fisher.test(O)
    plot.new();
    par(usr = c(0, 1, 0, 1));
    addtable2plot(0,0.5, O, display.colnames = T, display.rownames = T);
    text(0.5, 0.44, sprintf(paste('P = %0.', floor(-log10(fish$p.value))+2, 'f\nOR = %0.2f [%0.2f-%0.2f]', sep  = ""), fish$p.value, fish$estimate, fish$conf.int[[1]], fish$conf.int[[2]]))
  }

###############################
# fisher.combined
#
# Computes fisher combined p value for a matrix of p values where the columns correspond to individual tests
# rows correspond to hypotheses.
#
###############################
fisher.combined = function(Ps)
{
  if (is.vector(Ps))
    return(Ps)

  return(pchisq(rowSums(-2*log(Ps)), 2*ncol(Ps), lower.tail = F))
}

####################
# retrieves the region from CEU hapmap with r2>ld to the query region
####################
ldHood = function(region, ld)
  {
    if (is.null(region$pos2))
       {
         region$pos2 = region$pos1;
       }

     temp = ldQuery(region, ld);
     
     if (dim(temp)[1]>0)
       {
         pos = rbind(temp$pos1, temp$pos2);         
         out = data.frame(chr=temp$chr[1], pos1=min(pos), pos2=max(pos))
       }
     else
       {
         out = data.frame(chr=region$chr, pos1=region$pos1, pos2=region$pos2)
       }

    return(out)
  }

lsf.status = function()
  {
    TEMP.FILE = 'tmp.432124kiijoij1oij2oii1ijoi202082892h4iojh12oij12';
#    system(paste("bjobs -u all | awk \'{print", paste('$', c(1:7), sep = "", collapse = "\"\t\""), "\"\t\"", paste('$', c(8:10), sep = "", collapse = "\" \""), "}\' > ", TEMP.FILE, sep = ""))
    system(paste("bjobs -u all > ", TEMP.FILE, sep = ""))
        
    header = readLines(TEMP.FILE,1);
    headers = strsplit(header[[1]], "\\s+");
    headers = sub("\\_", "\\.", headers[[1]]);
    
    w = get.field.widths(header[[1]]);
    w[length(w)] = 1000;
   
#    bout = read.delim(TEMP.FILE, stringsAsFactors = FALSE, header = FALSE, skip = 1, col.names = c('JOBID','USER','STAT','QUEUE','FROM_HOST','EXEC_HOST','JOB_NAME','SUBMIT_TIME'));
    bout = read.fwf(TEMP.FILE, w, stringsAsFactors = FALSE)
    system(paste('rm ', TEMP.FILE));

    for (j in 1:dim(bout)[2])
        bout[,j] = trim(bout[,j]);       

    bout$SUBMIT.TIME = as.POSIXct(paste(bout$SUBMIT.TIME,  format(Sys.Date(), "%Y")), format =  '%b %d %H:%M %Y');
    
    bout
  }

lsf.usersum = function(lsfstat, attr = "STAT")
  {
    tab = table(lsfstat$USER, lsfstat[, attr]);
    tab = tab[order(-rowSums(tab)),]
    tab
  }

get.fwf.widths = function(file, skip=0)
  {
    l = readLines(file,skip+1);

    w = get.field.widths(l[[length(l)]]);    
  }


############################
# dev.all.off
#
# kills all windows
#
#############################
dev.all.off = function()
  {
    sapply(dev.list(), dev.off)
  }


# given logical vector returns data frame of starts and end indices of TRUE run lengths
run.lengths = function(vec)
  {
    out = data.frame();
    
    vec = c(vec, FALSE)
    last.bit = FALSE;
    this.bit = FALSE;
    last.true = as.double(this.bit);
    
    for (i in 1:length(vec))
      {
        last.bit = this.bit;
        this.bit = vec[i];
        
        if (!this.bit & last.bit)
          {
            out = rbind(out, data.frame(start = last.true, end = i-1))
          }
        
        if (this.bit & !last.bit)
          {
            last.true = i;
          }
      }
    
    out
  }

get.field.widths = function(str)
  {
    spl = strsplit(str, "")
    non.space = is.na(match(spl[[1]], " "));

    runs = run.lengths(non.space);

    out = c();
    
    if (dim(runs)[1]==1)
      {
        out = length(spl[[1]]);
      }
    else if (dim(runs)[1]>1)
      {         
        out = c(runs$start[2:dim(runs)[1]], length(spl[[1]])+1) -
          c(1, runs$start[2:dim(runs)[1]]);
      }

    out
  }

trim = function(str)
  {
    str = gsub("^\\s+", "", str);
    str = gsub("\\s+$", "", str);
    str
  }


chiblocks = function(b1, b2)
{
	tally = matrix(0, nrow=dim(hapblocks)[1], ncol =2); 
	tally[match(b1, hapblocks$BLOCKNAME), 1] = 1; 
	tally[match(b2, hapblocks$BLOCKNAME), 2] = 1;
	tab = table(as.factor(tally[,1]), as.factor(tally[,2])); 
	chisq.test(tab)
}

linblocks = function(b1, b2)
{
	tally = data.frame(matrix(0, nrow=dim(hapblocks)[1], ncol =2))
	names(tally) = c('b1', 'b2')
	
	ind = match(b1, hapblocks$BLOCKNAME)
	tally[ind[!is.na(ind)], 'b1'] = 1 

	ind = match(b2, hapblocks$BLOCKNAME)
	tally[ind[!is.na(ind)], 'b2'] = 1 

	lm('b2 ~ b1', tally)
}


listind = function(query, thelist)
{
     out = vector("list", length(query))
     for (i in (1:length(thelist)))
       {
          m = which(!is.na(match(query, thelist[[i]])))
          if (length(m)>0) {
            for (j in (1:length(m)))
              {
                out[[m[j]]] = c(out[[m[j]]], i)
              }
          }
       }
     out
}

#################
# write.tab - writes tab delimited no quotes without row names table (passes remaining arguments to write.table)
#
# = write.table(sep = "\t", quote = F, row.names = F)
#################
write.tab = function(..., sep = "\t", quote = F, row.names = F)
  {
    write.table(..., sep = sep, quote = quote, row.names = row.names)
  }


#################
# write.htab - writes data frame to pretty formatted rtable
#
#
#################
write.htab = function(tab, file,
  title = NULL, # text to be written in bold above the table  
  footer = NULL, # text to be writen in bold below the table
  highlight = NULL,  #vector of row indices of the table to highlight
  row.names = TRUE,  # includes row labels
  col.names = TRUE, # includes col labels
  high.color = 'yellow', # highlight color to use 
  row.colors = c('lightgray', 'white'), # alternating colors to shade data rows  
  header.colors = c('#4A4A4A', 'white'), # two element vector specifying background and text colors for header row, respectively,
  data.size = 15, # font size in px for data, title, and footer
  title.size = 15, footer.size = 20, header.size = round(1.1*data.size))
  {    
    require(hwriter)
    require(gplots)
    
    if (class(tab) != 'data.frame')
      tab = as.data.frame(tab)

    if (is.null(rownames(tab)))
      row.names = F;

    tab[is.na(tab)] = '';
    tab = tab[1:nrow(tab), , drop = FALSE];  #not sure why this is necessary, but deflects occasional weird R bug
    
    dir.create(dirname(file), recursive=T, showWarnings = F)
    p = openPage(file, link.css = 'hwriter.css')
    if (!is.null(title))
      hwrite(title, p, style = sprintf('font-weight:bold; font-size:%spx; margin-top;50px', title.size), center = TRUE, div = TRUE, br = TRUE);

    row.bgcolor = as.list(as.character(col2hex(row.colors)[(1:nrow(tab))%%length(row.colors)+1]));
    names(row.bgcolor) = rownames(tab)
    if (!is.null(highlight))
      row.bgcolor[rownames(tab[highlight,, drop = FALSE])] = list(col2hex(high.color));

    row.bgcolor = c(col2hex(header.colors[1]), row.bgcolor)

#    if (row.names)
      col.bgcolor = col2hex(header.colors[1])
    
    col.style = sprintf('font-weight:bold; font-size:%spx; color:%s; text-align:center', header.size, col2hex(header.colors[2]));
    
    row.style = rep(sprintf('font-size:%spx; text-align:center', data.size), nrow(tab))
    names(row.style) = rownames(tab)
    row.style = c(list(sprintf('font-weight:bold; font-size:%spx; color:%s; text-align:center', header.size, col2hex(header.colors[2]))), row.style)
    
    hwrite(tab, p, row.style = row.style, col.style = col.style, col.bgcolor = col.bgcolor, row.names = row.names, col.names = col.names,
           row.bgcolor = row.bgcolor, table.frame = 'void', table.style = 'margin-left: 30px; margin-top: 30px', br = TRUE)
    if (!is.null(footer))
      hwrite(footer, p, style = sprintf('font-weight:bold; text-align:center; font-size:%spx; margin-top;50px', footer.size), center = TRUE, div = TRUE);
    closePage(p)
  }

###############
# writecols
#
# Takes character vector and dumps into k delimited columns
###############
writeCols = function(v, k = 3, sep = "\t", file = "")
  {
    rows = ceiling(length(v)/k)
    out = matrix(ncol = k, nrow = rows, byrow = TRUE)
    out[1:length(v)] = v;
    out[is.na(out)] = "";
    write.table(as.data.frame(out), file = file, sep = sep, quote = F, row.names = F, col.names = F)
  }


##########################
# col.scale
#
# Assigns rgb colors to numeric data values in vector "x".. maps scalar values
# in val.range (default c(0,1)) to a linear color scale of between col.min (default white)
# and col.max (default black), each which are length 3 vectors or characters.  RGB values are scaled between 0 and 1. 
#
# Values below and above val.min and val.max are mapped to col.max and col.max respectively
##########################
col.scale = function(x, val.range = c(0, 1), col.min = 'white', col.max = 'black', na.col = 'white',
  invert = F # if T flips rgb.min and rgb.max
  )
  {
    if (!is.numeric(col.min))
      if (is.character(col.min))
        col.min = col2rgb(col.min)/255
      else
        error('Color should be either length 3 vector or character')

    if (!is.numeric(col.max))
      if (is.character(col.max))
        col.max = col2rgb(col.max)/255
      else
        error('Color should be either length 3 vector or character')

    col.min = as.numeric(col.min);
    col.max = as.numeric(col.max);

    x = (pmax(val.range[1], pmin(val.range[2], x))-val.range[1])/diff(val.range);
    col.min = pmax(0, pmin(1, col.min))
    col.max = pmax(0, pmin(1, col.max))

    if (invert)
      {
        tmp = col.max
        col.max = col.min
        col.min = tmp
      }

    nna = !is.na(x);
    
    out = rep(na.col, length(x))
    out[nna] = rgb((col.max[1]-col.min[1])*x[nna] + col.min[1],
        (col.max[2]-col.min[2])*x[nna] + col.min[2],
        (col.max[3]-col.min[3])*x[nna] + col.min[3])
    
    return(out)           
  }

##########################
# plot.blank
#
# Shortcut for making blank plot with no axes
##########################
plot.blank = function(xlim = c(0, 1), ylim = c(0,1), xlab = "", ylab = "", axes = F, ...)
  {
    plot(0, type = "n", axes = axes, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
#    par(usr = c(xlim, ylim))
  }

###########################
# axis.subplot
#
# function to draw simple horizontal and vertical axes at arbitrary positions on plot with
# decoupled data and plot coordinates (useful when there are many subplots within a single layout item)
#
###########################
axis.subplot = function(side = 1, # can be 1 (horizontal) or 2 (vertical)
  data.range, # (data coordinates) min and max data coordinate value for ends of axis line
              # corresponding to ypos (or xpos) plot coordinates in case of vertical (or horizontal) plot, respectively
  n.ticks = 5, # optimal number of ticks using "pretty" function
  at = NULL, # (data coordinates) tick locations same as in original axis function, over-rides n.ticks 
  ypos = NULL, # (plot coordinates) single value for horizontal, two values (range) for vertical correponding to ypos of axis line
  xpos = NULL, # (plot coordinates) single value for vertical, two values (range) for horizontal corresponding to xpos of axis line
  lwd = 1,  
  lwd.ticks = lwd,
  col = 'black',
  col.ticks = col,
  tick.length = 0.1, # in plot coordinates
  tick.orient = 1, # if >0 then facing "positive" direction, otherwise "negative"
  srt.tick = 0, # rotation of tick.labels relative to "default" (which is 0 for horizontal and 90 for vertical axis)
  cex.tick = 1, # size of tick labels
  adj.tick = c(0.5, 0.5), # justification of tick.labels (overrides defaults)
  ...) # arguments to line)
{  
  if (side == 1 & !(length(xpos)==2 & length(ypos)==1))
    stop('ypos must be of length 1 and xpos of length 2, for side = 1')
  else if (side == 2 & !(length(xpos)==1 & length(ypos)==2))
    stop('ypos must be of length 2 and xpos of length 1, for side = 2')
  else if (!(side %in% c(1, 2)))
    stop('side must be = 1 or = 2')
  
  axis.lim = data.frame(x = xpos, y = ypos);

  if (is.null(at))
    at = pretty(data.range, n = n.ticks)

  # draw axis line
  lines(axis.lim$x, axis.lim$y, lwd = lwd, col = col, ...)

  # draw ticks
  if (tick.length <= 0)
    tick.length = -tick.length;
  
  if (side ==1)
    {
      at.plot = coord.subplot(at, data.range = data.range, plot.range = xpos)
      segments(at.plot, rep(ypos, length(at)), y1 = rep(ypos, length(at))+tick.length, lwd = lwd, col = col, ...)
      text(rep(at.plot, 2), rep(ypos, length(at))+1.5*tick.length, at, adj = adj.tick, cex = cex.tick, srt = 0+srt.tick)
    }
  else if (side ==2)
    {      
      at.plot = coord.subplot(at, data.range = data.range, plot.range = ypos)
      segments(rep(xpos, length(at)), at.plot, x1 = rep(xpos, length(at))+tick.length, lwd = lwd, col = col, ...)
      text(rep(xpos, length(at))+1.5*tick.length, rep(at.plot, 2), at, adj = adj.tick, cex = cex.tick, srt = -90+srt.tick)      
    }      
}

##############################
# coord.subplot
#
# Scales data coordinates in given "data.range" to plot coordinates in "plot.range".  Can be used for x or y coordinates.
# Data coordinates that map outside of the data.range are given NA values. 
#
# Useful for doing many subplots on a single plot axis. 
# e.g. if we want to plot genomic positions 54,100,000 to 56,200,000 on plot coordinates 0.1 to 0.4
##############################
coord.subplot = function(data, # any length data vector  
              data.range, # length 2 vector 
              plot.range # length 2 vector
              )
{
  data[data<data.range[1] | data>data.range[2]] = NA;
  return(((data-data.range[1])/diff(data.range))*diff(plot.range)+plot.range[1])
}

############################
# strwrap2
#
# Wrapper around strwrap that returns a vector of characters the same length as input with new lines
# inserted in word breaks to obey a provided string width.
#
#############################
strwrap2 = function(x, width, sep = NULL, newline = '\n', ...)
  {
    return(sapply(strwrap(x, width = width, simplify = F),
                  function(x) paste(x, collapse = newline)))
  }

############################
# Capitalize first letter of each character element of vector "string"
# 
##############################
capitalize = function(string, un = FALSE)
{
    if (!un)
      {
        capped <- grep("^[^A-Z].*$", string, perl = TRUE)
        substr(string[capped], 1, 1) <- toupper(substr(string[capped],1, 1))
      }
    else
      {
        capped <- grep("^[A-Z].*$", string, perl = TRUE)
        substr(string[capped], 1, 1) <- tolower(substr(string[capped],1, 1))
      }
    
    return(string)
}

############################
# Performs fisher test on cols of matrix / df x vs cols of matrix / df y
#
# returns list with ncol(x) by ncol(y) matrices $p and $or denoting the p value and odds ratio of the result of the
# fisher test on col i of x and col j of y
#
# If y is not provided, will correlate rows of x with themselves. 
#############################
fisher.pairwise = function(x, y = x)
  {
    p = or = matrix(NA, nrow = ncol(x), ncol = ncol(y), dimnames = list(colnames(x), colnames(y)))

    if (nrow(x) != nrow(y))
      stop('x and y must have the same number of rows')
    
    logical.x = which(sapply(1:ncol(x), function(i) is.logical(x[,i])))
    logical.y = which(sapply(1:ncol(y), function(i) is.logical(y[,i])))

    if (length(logical.x)==0 | length(logical.y)==0)
      warning('No logical columns found')

    for (i in logical.x)
          for (j in logical.y)
            {
              O = table(x[,i], y[,j])
              if (min(dim(O))>1)
                {
                  res = fisher.test(O)          
                  or[i, j] = res$estimate
                  p[i, j] = res$p.value
                }
            }

    out = list(p = p, or = or)
    return(out)
  }


#################
# lighten
#
# lightens / darkens colors by brighness factor f (in -255 .. 255) that will make lighter if > 0 and darker < 0
################
lighten = function(col, f)
  {
    M = col2rgb(col)
    return(apply(matrix(pmax(0, pmin(255, M + f*matrix(rep(1, length(M)), nrow = nrow(M)))), ncol = length(col))/255, 2, function(x) rgb(x[1], x[2], x[3])))
  }

################
# fin.lines
#
# returns par('fin') in terms of lines - helps margin computation
#
################
fin.lines = function() par('fin')*(par('mar')/par('mai'))[1];  # help us make sure we don't encroach on figure margins


################
# strsplit2
#
# Strsplit when there are two layers of separators (sep1, sep2) and one needs to extract
# a collapsed vector of subitem j for all items i.
#
# Takes in a character vector and outputs a list of "separated" items
################
strsplit2 = function(x, sep1 = ",", sep2 = " ", j = 1)
  {
    return(lapply(strsplit(x, sep1), function(y) sapply(strsplit(y, sep2), function(z) z[[j]])))
  }


################ 
# timestamp
#
################
timestamp = function()
  {
    return(gsub('[^\\d]', '', as.character(Sys.time()), perl = T))
  }


###############
# img.link
#
# Returns vector of html image links to files "file" with text "caption"
#
# if embed = T, then will make img link, and additional arguments may be supplied to image tag (eg height, width)
################
img.link = function(file, caption = NULL, embed = F, ...) {
	if (is.null(caption)) {
		caption = ''
	}
	
	if (!embed) {
		return(paste('<a href = \"', file, '\">', caption, '</a>'))
	} else {
		args = list(...);
		if (length(args)>0) {
			to.quote = is.numeric(unlist(args))+1
			more.args = paste((paste(', ', names(args), " = ", c('\"', '')[to.quote], unlist(args), c('\"', '')[to.quote], sep = "")), collapse="")
		} else {
			more.args = ''
		}
		
		parts = strsplit(basename(file), "\\.")[[1]]
		file.ext = parts[length(parts)]
		if (file.ext == "tif") { # handles tif images which won't display in chrome with a regular img src tag
			out = paste('<embed src = \"', file, '\", type = "image/tiff"', more.args, ' negative=yes>', sep = "")
		} else {
			out = paste('<img src = \"', file, '\", alt = \"', caption, '\"', more.args, '>', sep = "")
		}
		
		return(out)
	}
}



#########
# html.tag
#
# makes a open and close html tag with optional text inside and optional (named) vector of name=value pairs contained inside of opening tag
#########
html.tag = function(tag, text = NULL, collapse = '\n', ## collapse is what will separate the tag and each "line" of text, alternate can be " "
  ...)
  {
    flags = unlist(list(...))
      
    if (!is.null(flags))
      {
        if (is.null(names(flags)))
          flag.str = ""
        else
          flag.str = paste(" ", paste(paste(names(flags), paste('"', flags, '"', sep = ""), sep = "=")), collapse = "")
      }
    else
      flag.str = ""
    
    return(paste(sprintf('<%s%s>', tag, flag.str), paste(text, collapse = collapse), sprintf('</%s>', tag), sep = collapse))
  }

################
# set.comp
#
# Compares two sets and outputs data frame with "left", "middle", "right" members
################
set.comp = function(s1, s2)
  {
    universe = sort(union(s1, s2));
    out = data.frame(stringsAsFactors = F)       
    tmp.comp = list(left = sort(setdiff(s1, s2)), middle = sort(intersect(s1, s2)), right = sort(setdiff(s2, s1)))
    lapply(names(tmp.comp), function(x) out[1:length(tmp.comp[[x]]),x] <<- tmp.comp[[x]])
    out[is.na(out)] = ''
    return(out)
  }


#####################################
# heatmap.plus
#
# Additional features:
#    allows several label tracks on top, bottom, left, and right, with separate top and bottom legend frames tohouse each
#    allows use of coloredData in tracks
# 
###################################
heatmap.plus = function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL,
  show.rdend = TRUE, # show row dendogram flag
  show.cdend = TRUE,  # show column dendrogram flag
  # these next four args can be coloredData or list of coloredData with category labels
  # for each column / row data point.  If named will be indexed by corresponding row/column names in X
  # otherwise they will be indexed in order.  These args can alternatively can be vector or list of vectors (which
  # will be mapped to default colormaps)
  topColAttr = NULL, # see above
  bottomColAttr = NULL, # ...
  leftRowAttr = NULL, # ...
  rightRowAttr = NULL, # ...
  leg.args = NULL,  ## legend will be populated by color mappings from coloredData and additional args (excluding position) here
  dim.heatmap = c(4, 4), ## use this to make heatmap tall or fat (instead of margins arg)
  distfun = dist, hclustfun = hclust,
  reorderfun = function(d, w) reorder(d, w),
  add.expr,  
  symm = FALSE,
  revC = identical(Colv, "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE,
  margins = c(1, 1), 
  cexRow = 0.2 +  1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
  add.grid = F,
  col.grid = 'gray',
  lwd.grid = 1,
  size.legend.panel = 0.4,
  size.feature.panel = 0.2,
  col = topo.colors(100),
  optimal.leaf = T,
  return.clust = F,    
  labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, las.col = 2,
  verbose = getOption("verbose"), ...) 
{
  require(cba)
  hcr = hcc = NULL

  ### set up top and bottom color palettes
  if (!is.null(topColAttr) | !is.null(bottomColAttr) | !is.null(leftRowAttr) | !is.null(rightRowAttr))
    {
      require(RColorBrewer)
      
      brewer.palettes = brewer.pal.info
      brewer.palettes = brewer.palettes[order(brewer.palettes$category), ]
  
      if (!is.null(topColAttr) & !is.list(topColAttr))
        topColAttr = list(topColAttr)
      
      if (!is.null(bottomColAttr) & !is.list(bottomColAttr))
        bottomColAttr = list(bottomColAttr)
      
      if (!is.null(leftRowAttr) & !is.list(leftRowAttr))
        leftRowAttr = list(leftRowAttr)
      
      if (!is.null(rightRowAttr) & !is.list(rightRowAttr))
        rightRowAttr = list(rightRowAttr)

      last.palette = 0;

      .convert.to.cData = function(x) {        
        if (!is.null(x) & class(x) != 'coloredData')
          {
            x = as.vector(x);
            last.palette <<- last.palette + 1;
            if (is.factor(x))
              uval = levels(x)
            else
              uval = unique(x)
            cmap = brewer.pal(min(brewer.palettes[last.palette, 'maxcolors'], length(uval)), rownames(brewer.palettes)[last.palette])

            if (length(uval)>length(cmap))
              {
                warning('Number of colors exceeded for colormap: duplicate colors will be created. ')
                cmap = cmap[((1:length(uval))%%length(cmap))+1]  ## we will repeat colors if the number of unique items is larger than the colormap
              }
            names(cmap) = uval;
            return(coloredData(data = x, colormap = cmap))
          }
        else
          return(x)
      }
      topColAttr = lapply(topColAttr, .convert.to.cData) ## if any attributes are not already colored data, transform them using coloredData using pre-applied palettes
      bottomColAttr = lapply(bottomColAttr, .convert.to.cData) ## if any attributes are not already colored data, transform them using coloredData using pre-applied palettes
      leftRowAttr = lapply(leftRowAttr, .convert.to.cData) ## if any attributes are not already colored data, transform them using coloredData using pre-applied palettes
      bottomColAttr = lapply(bottomColAttr, .convert.to.cData) ## if any attributes are not already colored data, transform them using coloredData using pre-applied palettes      
    }

  if (!is.null(colnames(x)))
    colNames = colnames(x)
  else
    colNames = 1:ncol(x)
  
  if (!is.null(rownames(x)))
    rowNames = rownames(x)
  else
    rowNames = 1:nrow(x)
   
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("'x' must be a numeric matrix")
    nr <- di[1L]
    nc <- di[2L]
    if (nr <= 1 || nc <= 1) 
        stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2L) 
        stop("'margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (!doRdend && identical(Colv, "Rowv")) 
        doCdend <- FALSE
    if (is.null(Rowv)) 
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv)) 
        Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram")) 
            ddr <- Rowv
        else {
            d = distfun(x);
            hcr = hclustfun(d);
            
            if (optimal.leaf) ## MARCIN ADDED
            {
              tmp <- order.optimal(d, hcr$merge)
              hcr$order = tmp$order;
              hcr$merge = tmp$merge;
              ddr <- as.dendrogram(hcr)
            }
          else
            {
              ddr <- as.dendrogram(hcr)
              if (!is.logical(Rowv) || Rowv) 
                ddr <- reorderfun(ddr, Rowv)
            }
        }
        if (nr != length(rowInd <- order.dendrogram(ddr))) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1L:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram")) 
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc) 
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
          d = distfun(if (symm) x else t(x));
          hcc = hclustfun(d);
            
          if (optimal.leaf) ## MARCIN ADDED
            {
              tmp <- order.optimal(d, hcc$merge)
              hcc$order = tmp$order;
              hcc$merge = tmp$merge;
              ddc <- as.dendrogram(hcc)
            }
          else
            {
              ddc <- as.dendrogram(hcc)
              if (!is.logical(Colv) || Colv) 
                ddc <- reorderfun(ddc, Colv)
            }
        }
        if (nc != length(colInd <- order.dendrogram(ddc))) 
          stop("column dendrogram ordering gave index of wrong length")
      }
    else colInd <- 1L:nc
    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow)) 
        if (is.null(rownames(x))) 
            (1L:nr)[rowInd]
        else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol)) 
        if (is.null(colnames(x))) 
            (1L:nc)[colInd]
        else colnames(x)
    else labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 1L, sd, na.rm = na.rm)
        x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
    }
    else if (scale == "column") {
        x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 2L, sd, na.rm = na.rm)
        x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
    }
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, dim.heatmap[1])
    lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
        dim.heatmap[2])
  
  ## COL LABEL LAYOUTS  
  col.panel.height = size.feature.panel;
  row.panel.width = size.feature.panel;
  core.panel.ind = c(2,2);
  
  if (!is.null(topColAttr))  ## add to top
    lapply(topColAttr, function(x)
           {
             lmat <<- rbind(lmat[1, ]+1, c(NA, 1), lmat[2:nrow(lmat), ]+1)             
             lhei <<- c(lhei[1L], col.panel.height, lhei[2:length(lhei)])
             core.panel.ind[1] <<- core.panel.ind[1]+1
           })

  if (!is.null(bottomColAttr)) # add to bottom
    lapply(bottomColAttr, function(x)
           {
             lmat <<- rbind(lmat+1, c(NA, 1))
             lhei <<- c(lhei, col.panel.height)
           })

  ## ROW LABEL LAYOUTS
  if (!is.null(leftRowAttr))
    lapply(leftRowAttr, function(x)
           {             
             lmat <<- cbind(lmat[, 1] + 1, c(rep(NA, core.panel.ind[1]-1),  1, rep(NA, nrow(lmat)-core.panel.ind[1])), lmat[, 2:ncol(lmat)] + 1)
             lwid <<- c(lwid[1L], row.panel.width, lwid[2L])
             core.panel.ind[2] <<- core.panel.ind[2]+1;
           })

  if (!is.null(rightRowAttr))
    lapply(rightRowAttr, function(x)
           {
             lmat <<- cbind(lmat + 1, c(rep(NA, core.panel.ind[1]-1),  1, rep(NA, nrow(lmat)-core.panel.ind[1])))
             lwid <<- c(lwid, row.panel.width)
           })

  ## ADD ROW COL LEGEND PANELS

  # top legend panel
  new.row = rep(NA, ncol(lmat)); new.row[core.panel.ind[2]] = max(lmat, na.rm = T)+1;
  core.panel.ind[1] = core.panel.ind[1]+1
  lhei = c(size.legend.panel, lhei)
  lmat = rbind(new.row, lmat)
  
  # bottom legend panel
  new.row = rep(NA, ncol(lmat)); new.row[core.panel.ind[2]] = max(lmat, na.rm = T)+1;
  lmat = rbind(lmat, new.row)  
  lhei = c(lhei, size.legend.panel)
  
  # left legend panel
  new.col = rep(NA, nrow(lmat)); new.col[core.panel.ind[1]] = max(lmat, na.rm = T)+1;
  core.panel.ind[2] = core.panel.ind[2]+1
  lmat = cbind(new.col, lmat)
  lwid = c(size.legend.panel, lwid)
  
  # right legend panel
  new.col = rep(NA, nrow(lmat)); new.col[core.panel.ind[1]] = max(lmat, na.rm = T)+1;
  lmat = cbind(lmat, new.col)
  lwid = c(lwid, size.legend.panel)
 
  lmat[is.na(lmat)] <- 0
  if (verbose) {
    cat("layout: widths = ", lwid, ", heights = ", lhei, 
        "; lmat=\n")
    print(lmat)
    }
  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  ## LAYOUT 
  layout(lmat, widths = lwid, heights = lhei, respect = TRUE) 

  if (verbose)
    print(lmat)
  
  pad.feature.panel = 0.1
  
  ## DRAW ROW LABEL PANELS
  if (!is.null(rowAttr <- c(leftRowAttr, rightRowAttr)))
    if (length(rowAttr)>0)
      lapply(rev(rowAttr), function(x)
             {
               if (is.null(x)) return
               par(mar = c(margins[1L]/2, pad.feature.panel, margins[1L]/2, pad.feature.panel))
               if (is.null(names(getData(x)))) ## we assume attributes are ordered
                 image(1, 1:nr, z = rbind(1L:nr), col = getColors(x)[rowInd], axes = FALSE, xlab = "", ylab = "")
               else   ## otherwise will use attribute names to properly order (assuming that data rownames are specified)r
                 image(1, 1:nr, z = rbind(1L:nr), col = getColors(x)[rowNames[rowInd]], axes = FALSE, xlab = "", ylab = "")
             })
  
  ## DRAW COLUMN LABEL PANELS 
  if (!is.null(colAttr <- c(topColAttr, bottomColAttr)))
    if (length(colAttr)>0)
      lapply(rev(colAttr), function(x)
             {
               if (is.null(x)) return
               par(mar = c(pad.feature.panel, margins[2L]/2, pad.feature.panel, margins[2L]/2))               
               if (is.null(names(getData(x)))) ## we assume attributes are ordered
                 image(1:nc, 1, z = cbind(1L:nc), col = getColors(x)[colInd], axes = FALSE, xlab = "", ylab = "")
               else ## otherwise will use attribute names to properly order (assuming that data colnames are specified)
                 image(1:nc, 1, z = cbind(1L:nc), col = getColors(x)[colNames[colInd]], axes = FALSE, xlab = "", ylab = "")
             })  
  
  par(mar = c(margins[1L]/2, margins[2L]/2, margins[1L]/2, margins[2L]/2))
  if (!symm || scale != "none") 
    x <- t(x)
  if (revC) {
    iy <- nr:1
    if (doRdend) 
      ddr <- rev(ddr)
    x <- x[, iy]
  }
  else iy <- 1L:nr
  
 ## CENTRAL PANEL (heatmap)
  xlim = 0.5 + c(0, nc);
  ylim = 0.5 + c(0, nr);
    image(1L:nc, 1L:nr, x, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "", ylab = "", col = col, ...)
  
    if (las.col == 2)
      axis(1, 1L:nc, labels = labCol, las = las.col, line = -0.5, tick = 0, cex.axis = cexCol)
    else
      axis(1, 1L:nc, labels = labCol, las = las.col, padj = 1, line = -0.5, tick = 0, cex.axis = cexCol)
 
     if (add.grid)
     {
       segments(-.5 + 1:(nc+1), -.5, -.5 + 1:(nc+1), .5 + nr, col= col.grid, lwd=lwd.grid)
       segments(-.5, -.5 + 1:(nr+1), .5 + nc, -.5 + 1:(nr+1), col= col.grid, lwd= lwd.grid)          
     }
  
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1L] - 1.25)

  
   axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2L] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))

  ## LEFT DENDROGRAM
  par(mar = c(margins[1L]/2, 0, margins[1L]/2, 0))
    if (doRdend & show.rdend)
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else
      if (!is.null(rownames(x)))
        {
          plot.blank(ylim = ylim, xlim = c(0, 1))
          text(rep(0.96, length(labCol)), (1:nr), labRow, srt = 0, adj = c(1, 0.5), cex = cexRow)
        }
      else
        frame()
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))

  ## TOP DENDROGRAM
  par(mar = c(0, margins[2L]/2, 0, margins[2L]/2))
    if (doCdend & show.cdend) 
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else
      if (!is.null(colnames(x)))
        {
          plot.blank(xlim = xlim, ylim = c(0, 1))
          text((1:nc), rep(0.04, length(labCol)), labCol, srt = 90, adj = c(0, 0.5), cex = cexCol)
        }
      else
        frame()
  
    if (!is.null(main)) {
        par(xpd = NA)
        title(main, cex.main = 1.5 * op[["cex.main"]])
    }

  ### LEGENDS

  # set defaults
  if (is.null(leg.args$y)) leg.args$y = 0.5;
  if (is.null(leg.args$x)) leg.args$x = 0.5;
  if (is.null(leg.args$adj)) leg.args$adj = c(0, 0.5);
  if (is.null(leg.args$xjust)) leg.args$xjust = 0.5;
  if (is.null(leg.args$yjust)) leg.args$yjust = 0.5;
  
  ## TOP
  par(mar = c(1,0,0,0))
  par(xpd = NA);
  plot(0:1,0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
  if (length(topColAttr)>0)
    {
      these.leg.args = lapply(1:length(topColAttr),  function(x) {
        out = leg.args;
        out$x = x/(length(topColAttr)+1)
        out$legend = names(getColormap(topColAttr[[x]]))
        out$fill = getColormap(topColAttr[[x]])
        return(out)
      })
      sapply(these.leg.args, function(x) do.call('legend', x)) ## call all the legends
    }

  ## BOTTOM
  par(mar = c(0,0,1,0))
  plot(0:1,0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
  if (length(bottomColAttr)>0)
    {
      these.leg.args = lapply(1:length(bottomColAttr),  function(x) {
        out = leg.args;
        out$x = x/(length(bottomColAttr)+1)
        out$legend = names(getColormap(bottomColAttr[[x]]))
        out$fill = getColormap(bottomColAttr[[x]])
        return(out)
      })
      sapply(these.leg.args, function(x) do.call('legend', x)) ## call all the legends
    }                       

  ## LEFT
  par(mar = c(0,0,0,1))
  plot(0:1,0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
  if (length(leftRowAttr)>0)
    {
      these.leg.args = lapply(1:length(leftRowAttr),  function(x) {
        out = leg.args;
        out$y = x/(length(leftRowAttr)+1)
        out$legend = names(getColormap(leftRowAttr[[x]]))
        out$fill = getColormap(leftRowAttr[[x]])
        return(out)
      })
      sapply(these.leg.args, function(x) do.call('legend', x)) ## call all the legends
    }

  ## RIGHT
  par(mar = c(0,1,0,0))
  plot(0:1,0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
  if (length(rightRowAttr)>0)
    {
      these.leg.args = lapply(1:length(rightRowAttr),  function(x) {
        out = leg.args;
        out$y = x/(length(rightRowAttr)+1)
        out$legend = names(getColormap(rightRowAttr[[x]]))
        out$fill = getColormap(rightRowAttr[[x]])
        return(out)
      })
      sapply(these.leg.args, function(x) do.call('legend', x)) ## call all the legends
    }                       
  
  invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
                                                     doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
  
  if (return.clust)
    return(list(row = hcr, col = hcc))

}

##########################
# col.scale
#
# Assigns rgb colors to numeric data values in vector "x".. maps scalar values
# in val.range (default c(0,1)) to a linear color scale of between col.min (default white)
# and col.max (default black), each which are length 3 vectors or characters.  RGB values are scaled between 0 and 1. 
#
# Values below and above val.min and val.max are mapped to col.max and col.max respectively
##########################
col.scale = function(x, val.range = c(0, 1), col.min = 'white', col.max = 'black', na.col = 'white',
  invert = F # if T flips rgb.min and rgb.max
  )
  {
    if (!is.numeric(col.min))
      if (is.character(col.min))
        col.min = col2rgb(col.min)/255
      else
        error('Color should be either length 3 vector or character')

    if (!is.numeric(col.max))
      if (is.character(col.max))
        col.max = col2rgb(col.max)/255
      else
        error('Color should be either length 3 vector or character')

    col.min = as.numeric(col.min);
    col.max = as.numeric(col.max);

    x = (pmax(val.range[1], pmin(val.range[2], x))-val.range[1])/diff(val.range);
    col.min = pmax(0, pmin(1, col.min))
    col.max = pmax(0, pmin(1, col.max))

    if (invert)
      {
        tmp = col.max
        col.max = col.min
        col.min = tmp
      }

    nna = !is.na(x);
    
    out = rep(na.col, length(x))
    out[nna] = rgb((col.max[1]-col.min[1])*x[nna] + col.min[1],
        (col.max[2]-col.min[2])*x[nna] + col.min[2],
        (col.max[3]-col.min[3])*x[nna] + col.min[3])
    
    return(out)           
  }


########################
# install.packages.bioc
#
# shortcut to install bioconductor packages
########################
install.packages.bioc = function(pkg)
  {
    source('http://bioconductor.org/biocLite.R')
    sapply(pkg, biocLite)
  }

##########################
# install.packages.github
#
# shortcut to install github packages
##########################
install.packages.github = function(pkg, username, branch)
  {
    library('devtools')
    install_github(repo = pkg, username = username, branch = branch)    
  }

########################
# coloredData
#
# simple object with data (e.g. vector or matrix of categorical, real numbers) + a colormap
# (playing with S4 classes)
#
# colormap is a (named) vector mapping factor levels / unique values in data into colors,
# or otherwise assigning a color range to numeric data.
########################
setClass('coloredData', representation(data = 'array', colormap = 'vector', type = 'character', data.names = 'character'),
         prototype(data = matrix(NA), colormap = NA, type = '', data.names = NULL)
         )
setMethod('initialize', 'coloredData', function(.Object, data, colormap, upright = T)
         {
           .Object = callNextMethod()
           .Object@type = class(data)

           if (!is.vector(colormap) | !any(!is.na(colormap)))
             stop('colormap must be a vector with non NA entries')
           
           .Object@colormap = colormap
           if (is.vector(data) & !is.list(data))
             {
               if (!is.null(names(data)))
                 .Object@data.names = names(data)
               if (upright)               
                 data = array(data, dim = c(length(data), 1), dimnames = list(names(data), NULL))
               else
                 data = array(data, dim = c(1, length(data)), dimnames = list(NULL, names(data)))
             }
           else if (is.matrix(data))
               data = array(as.vector(data), dim = dim(data), dimnames = dimnames(data))
           else if (!is.array(data))
             stop('Only vectors, matrices, and arrays supported')
           .Object@data = data
           if (!is.numeric(.Object@data)) {
             if (is.factor(.Object@data))
               uval = levels(.Object@data)
             else if (is.character(.Object@data))
               uval = unique(.Object@data)
             else
               uval = NULL;
             
             if (!is.null(uval))
               {             
                 if (is.null(names(.Object@colormap)))
                   {
                     ix = 1:min(length(uval), length(.Object@colormap))              
                     names(.Object@colormap)[ix] = uval[ix]
                   }
                 if (length(leftover <- setdiff(uval, c(NA, names(.Object@colormap))))>0)
                   warning(sprintf('The following factors are unmapped in the colormap: %s.', paste(leftover, collapse = ",")))
               }
           } 
           .Object
         })
setGeneric('getColormap', function(.Object) standardGeneric('getColormap'))
setGeneric('getData', function(.Object) standardGeneric('getData'))
setGeneric('getColors', function(.Object) standardGeneric('getColors'))
setMethod('getColormap', signature('coloredData'), function(.Object) .Object@colormap)
setMethod('getColors', signature('coloredData'), function(.Object)
          {
            dat = getData(.Object);
            cmap = getColormap(.Object);
            dat[1:length(dat)] = cmap[dat];
            return(dat)
          })
setMethod('getData', signature(.Object = 'coloredData'), function(.Object)
          {
            out = as(.Object@data, .Object@type)
            if (!is.null(.Object@data.names))
              names(out) = .Object@data.names
            return(out)
          })                            
coloredData = function(...) new('coloredData', ...)
  

##########################
# write.delim
#
# shortcut to write out a table as a tab delimited file
##########################
write.delim <- function(df, filename){
	write.table(df, filename, sep="\t", row.names= F, quote=F)
}
