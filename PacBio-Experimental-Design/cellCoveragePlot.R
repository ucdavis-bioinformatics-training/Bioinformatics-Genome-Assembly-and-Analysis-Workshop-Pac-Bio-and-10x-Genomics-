#
# This is an R script that produces a figure 
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#
## Requires packages installed
#install.packages(c("rjson",devtools","plotly"))

library(rjson)
library(plotly)

infile = "input.fofn.stats"
genome_size=1.5e9
target_raw_coverage = 70
est_per_smrt_cell = 1e9

## SMRT cell still needed before achieving estimate.
smrt_cells=ceiling(genome_size*target_raw_coverage/est_per_smrt_cell)

print(infile)
## read in an input file
json <- fromJSON(file=infile)

clens <- tapply(sapply(json$file_data,'[[','read_lengths'),sapply(json$file_data,'[[','cell_barcode'), function(x) do.call("c",x))

## Histogram of reads lengths
hist(do.call('c',clens),breaks=200,xlab="Read Length",main="Read lengths of all cells")

## read mean length
mlen <- sapply(clens,function(x){
  mean(x)
})
barplot(mlen,names.arg=1:length(mlen),xlab="Cell Number",ylab="Mean Length",main="Mean length by cell")

## read N50 length
n50 <- sapply(clens,function(x){
  sort(x)[which(cumsum(sort(x))/sum(sort(x)) >= 0.5)[1]]
})  
barplot(n50,names.arg=1:length(n50),xlab="Cell Number",ylab="N50 Length",main="N50 length by cell")

## Number of subreads
barplot(sapply(clens,length),names.arg=1:length(clens),xlab="Cell Number",ylab="sub-read count",main= "Subread count by Cell")

## Number of 'effective' ZMWs
zmws <- tapply(sapply(json$file_data,'[[','zmw'),sapply(json$file_data,'[[','cell_barcode'), function(x) do.call("c",x))
barplot(sapply(zmws,function(x) length(unique(x))),names.arg=1:length(zmws),xlab="Cell Number",ylab="zmw count",main= "ZMW count by Cell")

##################################################
## Coverage by cell
##################################################
k0 <- sapply(clens,function(x){
  sum(x)
})  

k8 <- sapply(clens,function(x){
  sum(x[x>8e3])
})  

k10 <- sapply(clens,function(x){
  sum(x[x>10e3])
})  

k12 <- sapply(clens,function(x){
  sum(x[x>12e3])
})  

plot(x=1:length(k0),y=cumsum(k0)/genome_size,type='l',col="black",ylim=c(0,target_raw_coverage), xlim=c(0,smrt_cells), ylab="Coverage",xlab="Cell Number",main="Coverage increase by cell")
lines(x=1:length(k0),y=cumsum(k8)/genome_size,type='l',col="red")
lines(x=1:length(k0),y=cumsum(k10)/genome_size,type='l',col="blue")
lines(x=1:length(k0),y=cumsum(k12)/genome_size,type='l',col="green")

esmrt <- smrt_cells - length(k0)
esmrt_sample <- sample(1:length(k0),size=esmrt,replace=T)

ek0 <- sapply(clens[esmrt_sample],function(x){
  sum(x)
})  

ek8 <- sapply(clens[esmrt_sample],function(x){
  sum(x[x>8e3])
})  

ek10 <- sapply(clens[esmrt_sample],function(x){
  sum(x[x>10e3])
})  

ek12 <- sapply(clens[esmrt_sample],function(x){
  sum(x[x>12e3])
})  

lines(x=(length(k0)+1):smrt_cells,y=(cumsum(ek0)+tail(cumsum(k0),1))/genome_size,type='l',lty=2,col="black")
lines(x=(length(k8)+1):smrt_cells,y=(cumsum(ek8)+tail(cumsum(k8),1))/genome_size,type='l',lty=2,col="red")
lines(x=(length(k10)+1):smrt_cells,y=(cumsum(ek10)+tail(cumsum(k10),1))/genome_size,type='l',lty=2,col="blue")
lines(x=(length(k12)+1):smrt_cells,y=(cumsum(ek12)+tail(cumsum(k12),1))/genome_size,type='l',lty=2,col="green")

####################################################################


###########
## Start Plotly version
covDat = data.frame(cellNumber = rep(1:length(k0), 4), 
                    coverage = c(cumsum(k0)/genome_size, cumsum(k8)/genome_size, cumsum(k10)/genome_size, cumsum(k12)/genome_size),
                    cutoff =factor(c(rep("0 kb", length(k0)), rep("8 kb", length(k0)), rep("10 kb", length(k0)), rep("12 kb", length(k0)))), stringsAsFactors = F)
covDat$cutoff = factor(covDat$cutoff, levels = levels(covDat$cutoff)[order(as.numeric(sapply(levels(covDat$cutoff), function(x) unlist(strsplit(x, " "))[1])))] )

covDatPred = data.frame(cellNumber = c((length(k0)+1):smrt_cells, (length(k8)+1):smrt_cells, (length(k10)+1):smrt_cells, (length(k12)+1):smrt_cells), 
                    coverage = c((cumsum(ek0)+tail(cumsum(k0),1))/genome_size, (cumsum(ek8)+tail(cumsum(k8),1))/genome_size, (cumsum(ek10)+tail(cumsum(k10),1))/genome_size, (cumsum(ek12)+tail(cumsum(k12),1))/genome_size),
                    cutoff =factor(c(rep("0 kb", smrt_cells-length(k0)), rep("8 kb", smrt_cells-length(k0)), rep("10 kb", smrt_cells-length(k0)), rep("12 kb", smrt_cells-length(k0)))), stringsAsFactors = F)
covDatPred$cutoff = factor(covDatPred$cutoff, levels = levels(covDatPred$cutoff)[order(as.numeric(sapply(levels(covDatPred$cutoff), function(x) unlist(strsplit(x, " "))[1])))] )

p = plot_ly(covDat, x=~cellNumber, y=~coverage, color=~cutoff, mode="lines", type="scatter") %>% layout(xaxis = list(range = c(0,smrt_cells)), yaxis = list(range = c(0,target_raw_coverage))) %>%
  add_lines(x=covDatPred$cellNumber, y=covDatPred$coverage, color=covDatPred$cutoff, line = list(dash="dash"), showlegend = FALSE)
p
htmlwidgets::saveWidget(as.widget(p), "coverageByCell.html")

# /////////
