
#####################################################################
#IF RSTRINGMOL LIBRARY IS AVAILABLE...
#suppressPackageStartupMessages(library(Rstringmol,quietly = TRUE));

#####################################################################
#ELSE, WE NED TO DEFINE FUNCTIONS...

suppressPackageStartupMessages(library(Matrix, quietly=TRUE))

message('loading popdy-class...')
#source('popdy-class.R')
############
sink("tmp.txt")
setClass("popdy", representation(tracks="dgCMatrix", times="numeric"))
setClass("popdat", contains="popdy")

setMethod("show", signature(object="popdy"), function(object){
	cat(sprintf("stringmol popdy object containing %d species across %d times (%d-%d)\n", ncol(object@tracks), nrow(object@tracks), min(object@times), max(object@times)))
})

setMethod("show", signature(object="popdat"), function(object){
	cat(sprintf("stringmol popdat object containing data for %d species across %d times (%d-%d)\n", ncol(object@tracks), nrow(object@tracks), min(object@times), max(object@times)))
})


setMethod("[", signature(x="popdy", i="ANY", j="missing"), function (x, i){
	return(x@tracks[,i])
})

setMethod("names", signature(x="popdy"), function(x){
	return(colnames(x@tracks))
})

setGeneric("plot")
setGeneric("lines")

setMethod("lines", signature(x="popdy"), function(x, col="black", plot.zero=FALSE, ...){
	col <- rep(col, ncol(x@tracks))
	for(i in 1:ncol(x@tracks)){
		# TODO: THIS IS A CHEAP HACK CURRENTLY!
		# CHECK THIS!
		y <- x@tracks[,i]
		if(identical(plot.zero, FALSE)) y[y==0] <- NA
		lines(x=times(x), y=y, col=col[i], ...)
	}
})

setMethod("plot", signature(x="popdy", y="missing"), function(x, col.fun=rainbow, add.lines=TRUE, add.points=FALSE, epochs=NA, epoch.border=NA, epoch.col.fun=cm.colors, epoch.alpha=0.25, axis.exponent=0, draw.box=TRUE, plot.zero=FALSE){
	x.times <- times(x)
	x.species <- speciesnames(x)
	x.range <- range(x.times)
	y.range <- range(tracks(x), na.rm=TRUE) # THIS IS MODIFIED TO ALLOW EA FAKE POPDY OBJECTS HERE.
	species.colours <- col.fun(length(x.species))
	if(!missing(epochs)) epoch.bg <- rgb(t(col2rgb(epoch.col.fun(nrow(epochs)))), alpha=(epoch.alpha * 255), maxColorValue=255)
	plot(NA, xlim=x.range, ylim=y.range, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
	if(!missing(epochs)) rect(xleft=epochs$start, ybottom=y.range[1], xright=epochs$end, ytop=y.range[2], col=epoch.bg, border=epoch.border)
	if(identical(add.points, TRUE)) lines(x, col=species.colours, plot.zero=plot.zero, type='p', pch=16, cex=0.5)
	if(identical(add.lines, TRUE)) lines(x, col=species.colours, plot.zero=plot.zero)
	axis(2, las=2)
	title(ylab='population')
	x.time.positions <- axTicks(1)
	if(missing(axis.exponent) || identical(axis.exponent, 0)){
		axis(1, at=x.time.positions, labels=x.time.positions)
		title(xlab=expression(time))
	} else {
		axis.labels <- x.time.positions / (10 ^ axis.exponent)
		axis(1, at=x.time.positions, labels=axis.labels)
		title(xlab=as.expression(substitute(time%*%10^e, list(e=axis.exponent))))
	}
	if(identical(draw.box, TRUE)) box()
})
sink(NULL)


############




read.popdy <- function(file){
	file.data <- read.table(file, sep=",")
	species.ids <- unique(as.character(file.data[,2]))
	times <- sort(unique(as.numeric((file.data[,1]))))
	time.index <- match(file.data[,1], times)
	species.index <- match(file.data[,2], species.ids)
	tracks <- as(new("dgTMatrix", i=as.integer(time.index - 1), j=as.integer(species.index - 1), x=as.double(file.data[,3]), Dim=c(length(times), length(species.ids))), "dgCMatrix")
	colnames(tracks) <- species.ids
	return(new("popdy", tracks=tracks, times=times))
}


times <- function(x){
	return(x@times)
}


nspecies <- function(x){
	return(ncol(x@tracks))
}


qnn.summarize <- function(x){
	return(sum(apply((x@tracks), 2, sum, na.rm=TRUE), na.rm=TRUE))
}


qnn.activity <- function(x, scale=FALSE){
	population.per.time <- apply(x@tracks, 1, sum)
	species.proportions <- sweep(x@tracks, 1, population.per.time, FUN="/")
	expected.species.proportions.per.time <- rBind(rep(NA, nspecies(x)), species.proportions[1:(nrow(species.proportions) - 1),])
	diff.expected <- (species.proportions - expected.species.proportions.per.time)
	diff.expected[diff.expected < 0] <- 0
	diff.expected <- diff.expected ^ 2
	if(identical(scale, TRUE)) diff.expected <- diff.expected * population.per.time
	output <- as(x, "popdat")
	output@tracks <- as(diff.expected, "dgCMatrix")
	return(output)
}

#FI
#####################################################################




message("Calculating qnn...\n");
x <- read.popdy("popdy.dat");

#xqnn<-qnn.activity(x,TRUE);
xqnn<-qnn.activity(x);
qnn<-qnn.summarize(xqnn);

#This was just a hack to see if we could evolve novel-making runs after maximising length..
#qnn<-1000000+qnn

#the following is a test to see if we can influence *anything* about the run using the GA
#qnn<-max(times(x));

message(sprintf("qnn is %e, runtime is %d, %f, %e",qnn, max(times(x)),qnn,qnn));

sink("qnn.txt")
cat(    sprintf(  "%e\n",qnn))
sink()

q("no");x
