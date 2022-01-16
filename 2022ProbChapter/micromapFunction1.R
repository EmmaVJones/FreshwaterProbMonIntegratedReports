indicatorRanges <- list(DO = list(breaks = list(7,8,9,10,11),
                                  units = 'mg/L',
                                  prettyName = 'Dissolved Oxygen'),
                        pH = list(breaks = list(6, 7, 8, 9), 
                                  units = 'unitless',
                                  prettyName = 'pH'),
                        SpCond = list(0, 200, 400, 600, 800),
                        TN = list(0, 0.5, 1.0, 1.5),
                        TP = list(0, 0.02, 0.04,0.06, 0.08),
                        TDS = list(0, 50, 100, 150, 200, 250),
                        NH4 =list(0, 0.01, 0.02, 0.03, 0.04, 0.05),              
                        NO3 = list(0, 0.25, 0.5, 0.75, 1.0),
                        TKN = list(0, 0.2, 0.4, 0.6),              
                        Ortho_P =list(0, 0.015, 0.03, 0.045),            
                        Turb = list(0, 5, 10, 15, 20),
                        TSS = list(0, 5, 10, 15, 20,25),               
                        Na =list(0, 5, 10, 15, 20),                
                        K = list(0, 1, 2, 3, 4),                 
                        Cl = list(0, 10, 20, 30, 40),  
                        Sf = list(0, 25, 50, 75, 100, 125,150),               
                        X70331VFine = list(0, 25, 50, 75, 100),      
                        SSCCOARSE =list(0, 1, 2, 3, 4),         
                        SSCFINE = list(0, 2, 4, 6, 8, 10),     
                        SSCTOTAL = list(0, 5, 10, 15, 20),    
                        LRBS = list(-2.5, -1.5, -0.5, 0, 0.5),             
                        Slope = list(0, 1, 2, 3),            
                        FN_PCT = list(0, 10, 20, 30),
                        SA_PCT = list(0, 10, 20, 30, 40, 50),
                        SA_FN_PCT  = list(0, 25, 50, 75, 100),   
                        LSUB_DMM = list(-1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0),      
                        BL_CB_GR_Embed_PCT = list(0, 25, 50, 75, 100),
                        Embed_PCT = list(0, 25, 50, 75, 100),      
                        TotHab = list(100, 125, 150, 175, 200),   
                        TotTaxa = list(0, 5, 10, 15, 20),       
                        EPTTax = list(0, 3, 6, 9, 12), 
                        VEphem = list(0, 10, 20, 30, 40),   
                        VPTHydropsychidae  = list(0, 10, 20, 30, 40),
                        VScrap = list(0, 10, 20, 30, 40),           
                        VChiro  = list(0, 10, 20, 30, 40), 
                        V2Dom  = list(20, 40, 60, 80),  
                        HBI  = list(2, 3, 4, 5, 6),    
                        EPTInd   = list(2, 3, 4, 5, 6),      
                        VSCIVCPMI  = list(0, 25, 50, 75, 100),
                        MetalCCU = list(0, 0.5, 1.0, 1.5),           
                        CALCIUM = list(0, 5, 10, 15),         
                        MAGNESIUM=list(0, 5, 10, 15, 20),      
                        ARSENIC = list(0, 0.2, 0.4, 0.6),
                        BARIUM  = list(0, 15, 30, 45),
                        BERYLLIUM = list(0, 0.02, 0.04,0.06, 0.08, 0.1),
                        CADMIUM  = list(0, 0.02, 0.04,0.06, 0.08, 0.1),  
                        CHROMIUM = list(0, .10, .20, .30, .40, .5, .6),
                        COPPER = list(0, 0.02, 0.04,0.06, 0.08, 0.1),   
                        IRON = list(0, 200, 400,600, 800),        
                        LEAD = list(0, 0.5, 0.10, 0.15, 0.20),     
                        MANGANESE = list(0, 20, 40, 60, 80),
                        THALLIUM = list(0, 0.02, 0.04,0.06, 0.08, 0.1), 
                        NICKEL = list(0, 0.02, 0.04,0.06, 0.08, 0.1), 
                        SILVER = list(0, 0.02, 0.04,0.06),
                        ZINC = list(0, 0.10, 0.20, 0.3),
                        ANTIMONY  = list(0, 0.10, 0.20, 0.30),  
                        ALUMINUM   = list(0, 10, 20, 30, 40, 50, 60),
                        SELENIUM = list(0, 0.10, 0.20, 0.30, 0.40, 0.5, 0.60),  
                        HARDNESS  = list(0, 50, 100, 150, 200),   
                        MERCURY  = list(0, 0.5, 1.00, 1.50, 2.00, 2.5) , 
                        wshdImpPCT  = list(0, 0.5, 1.00, 1.50, 2.00, 2.5, 3.0) ) 

# function to build dataset for interquartile range micromap
dataOrgInterquartileRange <- function(dat, # CDF data with all indicators and subpopulations
                                      subpopulations, # vector of subpopulations to pull 
                                      indicator # indicator to pull
){
  statsbasin <- tibble(Subpopulation = as.character(NA),
                       Indicator = as.character(NA),
                           x25= as.numeric(NA), x50= as.numeric(NA),
                           x75= as.numeric(NA), n= as.numeric(NA))
  for(i in 1:length(subpopulations)){
    subpopData <- filter(dat, Subpopulation == !! subpopulations[i] & Indicator == !! indicator)%>%
      select(Estimate.P,Value, StdError.P, NResp,everything())
    
    statsbasin[i,] <- tibble(Subpopulation = subpopulations[i],
                             Indicator = indicator, 
                             x25=vlookup(25,subpopData,2,TRUE),
                             x50=vlookup(50,subpopData,2,TRUE),
                             x75=vlookup(75,subpopData,2,TRUE),
                             n=max(subpopData$NResp))
    }
  return(statsbasin)
}
# organizedStats <- dataOrgInterquartileRange(dat, 
#                                     c('Chowan', 'Rappahannock', 'York', 'Potomac',
#                                            'Shenandoah', 'Roanoke Basin', 'James Basin',
#                                            'New', 'Big Sandy', 'Clinch-Powell', 'Holston', 
#                                       "Virginia"),
#                                     'DO'#'VSCIVCPMI'
#                                     ) %>% 
#   mutate(Subpopulation = case_when(Subpopulation == 'Roanoke Basin' ~ "Roanoke",
#                                    Subpopulation == 'James Basin' ~ "James",
#                                    TRUE ~ as.character(Subpopulation)),
#     Subpopulation = as.factor(Subpopulation),
#          Indicator = as.factor(Indicator))


# function to plot interquartile range by subpopulation and indicator
micromapInterquartileRange <- function(organizedStats, indicatorRanges, dropVA){
  # maxX <- round(max(organizedStats$x75)*1.2, digits = 2)
  # if(unique(organizedStats$Indicator) != 'LRBS'){minX <- 0
  # }else{minX <- round(min(organizedStats$x25)-(min(organizedStats$x25)*1.2), digits = 2)}
  # 
  # Xscale <- round(seq(from = minX, to = maxX, length.out = 5), digits = 2)
  
  plotRanges <- as.list(indicatorRanges[[as.character(unique(organizedStats$Indicator))]]$breaks)
  
  statewideMedian <- filter(organizedStats, Subpopulation == 'Virginia')$x50
  
  # drop VA from micromap plot if desired
  if(dropVA == T){ organizedStats <- filter(organizedStats, Subpopulation != 'Virginia')} 
  
  suppressWarnings(
    mmplot(stat.data=as.data.frame(organizedStats),
           map.data=map.table,
           map.link=c("Subpopulation", "ID"),
           panel.types=c('dot_legend', 'labels','labels', 'dot_cl', 'map'),
           panel.data=list(NA,'Subpopulation','n',list('x50', 'x25', 'x75'),NA),
           ord.by='x50',
           grouping=3,
           median.row=F,
           plot.height=7,
           plot.width=8,
           colors=brewer.pal(3, "Spectral"),
           rev.ord=T,
           panel.att=list(list(1, point.type=20, point.border=TRUE, point.size=2),
                          list(2, header='Subpopulation', panel.width=.5, 
                               align='left', text.size=.9),
                          list(3,header='n',panel.width=.2,align='left',text.size=.9),
                          list(4, header=paste0('Estimated Median ',
                                                indicatorRanges[[as.character(unique(organizedStats$Indicator))]]$prettyName,
                                                #unique(organizedStats$Indicator),
                                                ' and \nAssociated Interquartile Range'),
                               graph.bgcolor='lightgray', point.size=1.5,
                               xaxis.ticks=plotRanges, xaxis.labels=plotRanges,
                               #xaxis.ticks=list(Xscale), xaxis.labels=list(Xscale)
                               add.line=statewideMedian,
                               add.line.col='black',add.line.typ='dashed',
                               xaxis.title=indicatorRanges[[as.character(unique(organizedStats$Indicator))]]$units),
                                 #as.character(unique(organizedStats$Indicator))),
                          list(5, header='Light Gray Means\nPreviously Displayed',
                               map.all=TRUE, fill.regions='aggregate',
                               active.border.color='black', active.border.size=1.0,
                               inactive.border.color=gray(.7), inactive.border.size=1, 
                               panel.width=1.0))) )
}

micromapInterquartileRange(dataOrgInterquartileRange(dat, 
                                                     c('Chowan', 'Rappahannock', 'York', 'Potomac',
                                                       'Shenandoah', 'Roanoke Basin', 'James Basin',
                                                       'New', 'Big Sandy', 'Clinch-Powell', 'Holston', 
                                                       "Virginia"),
                                                     'DO'#'VSCIVCPMI'#'TotHab'#'TP'#'VSCIVCPMI'
                                                     ) %>% 
                             mutate(Subpopulation = case_when(Subpopulation == 'Roanoke Basin' ~ "Roanoke",
                                                              Subpopulation == 'James Basin' ~ "James",
                                                              TRUE ~ as.character(Subpopulation)),
                                    Subpopulation = as.factor(Subpopulation),
                                    Indicator = as.factor(Indicator)),
                           indicatorRanges = indicatorRanges,
                           dropVA = F)

