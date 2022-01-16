# script to functionalize micromap development
# R 3.6.1

suppressPackageStartupMessages(library(tidyverse))#1.3.0
suppressPackageStartupMessages(library(sp)) #1.3-2
suppressPackageStartupMessages(library(rgdal)) #1.4-8
suppressPackageStartupMessages(library(micromap))#1.9.3

basinssmooth_sp <- readOGR('originalData/GIS','VAbasins_smoothNoChesPeeDee')# basin shapefile for micromaps (sp)


# VLOOKUP (Excel function hack) by Julin Maloof
vlookup <- function(ref, #the value or values that you want to look for
                    table, #the table where you want to look for it; will look in first column
                    column, #the column that you want the return data to come from,
                    range=FALSE, #if there is not an exact match, return the closest?
                    larger=FALSE) #if doing a range lookup, should the smaller or larger key be used?)
{
  # 2020 addition, make tibbles dataframes
  table <- as.data.frame(table)
  
  if(!is.numeric(column) & !column %in% colnames(table)) {
    stop(paste("can't find column",column,"in table"))
  }
  if(range) {
    if(!is.numeric(table[,1])) {
      stop(paste("The first column of table must be numeric when using range lookup"))
    }
    table <- table[order(table[,1]),] 
    index <- findInterval(ref,table[,1])
    if(larger) {
      index <- ifelse(ref %in% table[,1],index,index+1)
    }
    output <- table[index,column]
    output[!index <= dim(table)[1]] <- NA
    
  } else {
    output <- table[match(ref,table[,1]),column]
    output[!ref %in% table[,1]] <- NA #not needed?
  }
  dim(output) <- dim(ref)
  output
}


# CDF data
dat <- read.csv('processedData/allCDF.csv')# CDF results


## plot with 25th-75th percentiles

chowan <-  filter(dat,Subpopulation=='Chowan'&Indicator=='VSCIVCPMI')%>%
  select(Estimate.P,Value, StdError.P, NResp,everything())
rappahannock <-  filter(dat,Subpopulation=='Rappahannock'&Indicator=='VSCIVCPMI')%>%
  select(Estimate.P,Value, StdError.P, NResp,everything())
york <-  filter(dat,Subpopulation=='York'&Indicator=='VSCIVCPMI')%>%   
  select(Estimate.P,Value, StdError.P, NResp,everything())
potomac <-  filter(dat,Subpopulation=='Potomac'&Indicator=='VSCIVCPMI')%>%   
  select(Estimate.P,Value, StdError.P, NResp,everything())
shenandoah <-  filter(dat,Subpopulation=='Shenandoah'&Indicator=='VSCIVCPMI')%>%  
  select(Estimate.P,Value, StdError.P, NResp,everything())
roanoke <-  filter(dat,Subpopulation=='Roanoke Basin'&Indicator=='VSCIVCPMI')%>%   
  select(Estimate.P,Value, StdError.P, NResp,everything())
james <-  filter(dat,Subpopulation=='James Basin'&Indicator=='VSCIVCPMI')%>%   
  select(Estimate.P,Value, StdError.P, NResp,everything())
new <-  filter(dat,Subpopulation=='New'&Indicator=='VSCIVCPMI')%>%  
  select(Estimate.P,Value, StdError.P, NResp,everything())
bigsandy <-  filter(dat,Subpopulation=='Big Sandy'&Indicator=='VSCIVCPMI')%>% 
  select(Estimate.P,Value, StdError.P, NResp,everything())
cp <-  filter(dat,Subpopulation=='Clinch-Powell'&Indicator=='VSCIVCPMI')%>% 
  select(Estimate.P,Value, StdError.P, NResp,everything())
holston <-  filter(dat,Subpopulation=='Holston'&Indicator=='VSCIVCPMI')%>% 
  select(Estimate.P,Value, StdError.P, NResp,everything())




map.table <- create_map_table(basinssmooth_sp,'BASIN')
# pull median VSCI scores and upper/lower bounds for micromap
stats <- data.frame(Basin=c('Chowan','Rappahannock','York','Potomac','Shenandoah','Roanoke','James','New',
                            'Big Sandy','Clinch-Powell', 'Holston'),Indicator='VSCIVCPMI')
# Switch lookup column to search for median
chowan <- select(chowan,Estimate.P,everything())
rappahannock <- select(rappahannock,Estimate.P,everything())
york <- select(york,Estimate.P,everything())
potomac <- select(potomac,Estimate.P,everything())
shenandoah <- select(shenandoah,Estimate.P,everything())
roanoke <- select(roanoke,Estimate.P,everything())
james <- select(james,Estimate.P,everything())
new <- select(new,Estimate.P,everything())
bigsandy <- select(bigsandy,Estimate.P,everything())
cp <- select(cp,Estimate.P,everything())
holston <- select(holston,Estimate.P,everything())

statsbasin <- rbind(
  data.frame(x25=vlookup(25,chowan,2,TRUE),x50=vlookup(50,chowan,2,TRUE),
             x75=vlookup(75,chowan,2,TRUE),n=max(chowan$NResp)),
  data.frame(x25=vlookup(25,rappahannock,2,TRUE),x50=vlookup(50,rappahannock,2,TRUE),
             x75=vlookup(75,rappahannock,2,TRUE),n=max(rappahannock$NResp)),
  data.frame(x25=vlookup(25,york,2,TRUE),x50=vlookup(50,york,2,TRUE),
             x75=vlookup(75,york,2,TRUE),n=max(york$NResp)),
  data.frame(x25=vlookup(25,potomac,2,TRUE),x50=vlookup(50,potomac,2,TRUE),
             x75=vlookup(75,potomac,2,TRUE),n=max(potomac$NResp)),
  data.frame(x25=vlookup(25,shenandoah,2,TRUE),x50=vlookup(50,shenandoah,2,TRUE),
             x75=vlookup(75,shenandoah,2,TRUE),n=max(shenandoah$NResp)),
  data.frame(x25=vlookup(25,roanoke,2,TRUE),x50=vlookup(50,roanoke,2,TRUE),
             x75=vlookup(75,roanoke,2,TRUE),n=max(roanoke$NResp)),
  data.frame(x25=vlookup(25,james,2,TRUE),x50=vlookup(50,james,2,TRUE),
             x75=vlookup(75,james,2,TRUE),n=max(james$NResp)),
  data.frame(x25=vlookup(25,new,2,TRUE),x50=vlookup(50,new,2,TRUE),
             x75=vlookup(75,new,2,TRUE),n=max(new$NResp)),
  data.frame(x25=vlookup(25,bigsandy,2,TRUE),x50=vlookup(50,bigsandy,2,TRUE),
             x75=vlookup(75,bigsandy,2,TRUE),n=max(bigsandy$NResp)),
  data.frame(x25=vlookup(25,cp,2,TRUE),x50=vlookup(50,cp,2,TRUE),
             x75=vlookup(75,cp,2,TRUE),n=max(cp$NResp)),
  data.frame(x25=vlookup(25,holston,2,TRUE),x50=vlookup(50,holston,2,TRUE),
             x75=vlookup(75,holston,2,TRUE),n=max(holston$NResp)))

stats <- cbind(stats,statsbasin)



suppressWarnings(
  mmplot(stat.data=stats,
         map.data=map.table,
         map.link=c("Basin", "ID"),
         panel.types=c('dot_legend', 'labels','labels', 'dot_cl', 'map'),
         panel.data=list(NA,'Basin','n',list('x50', 'x25', 'x75'),NA),
         ord.by='x50',
         grouping=3,
         median.row=F,
         plot.height=7,
         plot.width=8,
         colors=brewer.pal(3, "Spectral"),
         rev.ord=T,
         panel.att=list(list(1, point.type=20, point.border=TRUE, point.size=2),
                        list(2, header='Basin', panel.width=.5, 
                             align='left', text.size=.9),
                        list(3,header='n',panel.width=.2,align='left',text.size=.9),
                        list(4, header='Estimated Median VSCI/VCPMI Score and \nAssociated Interquartile Range',
                             graph.bgcolor='lightgray', point.size=1.5,
                             xaxis.ticks=list(40,50,60,70,80), xaxis.labels=list(40,50,60,70,80)
                             ,add.line=60,add.line.col='black',add.line.typ='dashed',
                             xaxis.title='VSCI Score'),
                        list(5, header='Light Gray Means\nPreviously Displayed',
                             map.all=TRUE, fill.regions='aggregate',
                             active.border.color='black', active.border.size=1.0,
                             inactive.border.color=gray(.7), inactive.border.size=1, 
                             panel.width=1.0))) )

