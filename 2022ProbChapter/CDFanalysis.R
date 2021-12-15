## Subpopulation estimates by parameter
subpopEstimate <- function(finalData,parameter, subpopulationCategory, subpopulation, altName, specialWeight){
  # Build in catch in case whole population needed
  if(is.na(subpopulationCategory) & subpopulation == "Virginia"){
    subpopData <- finalData
    sites.ext <- select(subpopData, siteID) %>% mutate(Use = TRUE)
    subpop.ext <- select(subpopData,siteID) %>% mutate(Region = 'Virginia')
  }else{
    subpopData <- finalData[ finalData[[subpopulationCategory]] == subpopulation, ] 
    if(nrow(subpopData) == nrow(finalData)){ # special catch for IR windows not filtering correctly
      subpopData <- finalData[ !is.na(finalData[[subpopulationCategory]] ), ] 
    }}
  
  # If no data in subpopulation, keep moving
  if(nrow(subpopData) !=0){
    sites.ext <- select(subpopData, siteID) %>% mutate(Use = TRUE)
    # Special Cases to match existing terminology for each Subpopulation
    if(is.na(altName)){
      subpop.ext <- select(subpopData,siteID) %>% mutate(Region = subpopulation)
    }else{
      subpop.ext <- select(subpopData,siteID) %>% mutate(Region = altName)
    }
    
    
    
    # Choose correct final weight and filter to stratum = 1
    if(specialWeight==FALSE){
      finalweight <- subpopData$finalweight_all
    }else{
      finalweight <- as.numeric(as.matrix(subpopData[,specialWeight]))
    }
    
    
    design.ext <- mutate(subpopData,siteID = siteID, stratum = "1", wgt = finalweight, xcoord = xmarinus, ycoord = ymarinus) %>%
      select(siteID, stratum, wgt, xcoord, ycoord)
    
    data.cont.ext <- select(finalData, siteID, parameter)
    
    subpopStatus <- cont.analysis(sites = data.frame(sites.ext), subpop = data.frame(subpop.ext), design = data.frame(design.ext), 
                                  data.cont = data.frame(data.cont.ext),
                                  pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
                                  conf=95, vartype="Local")
    return(subpopStatus)}else{return(list())}
}

# Now build it all into a single function that returns a list of lists (one list object per estimate, 
# each estimate is a list of of 4 dataframes where the important ones ($CDF and $Pct) can easily be 
# reached with list calls noted at bottom of page)
listOfResults <- function(popstatus.est, finalData, parameterName){
  list(
    popstatus.est = popstatus.est,
    dataAnalyzed = finalData,
    # All Virginia
    estimateVirginia = subpopEstimate(finalData,parameterName,NA, 'Virginia',NA, specialWeight=FALSE),
    # By Basin
    estimateRoanoke = subpopEstimate(finalData,parameterName, 'Basin', 'Roanoke','Roanoke Basin', specialWeight=FALSE),
    estimateJames = subpopEstimate(finalData,parameterName, 'Basin', 'James','James Basin', specialWeight=FALSE) ,
    estimatePotomacShenandoah = subpopEstimate(finalData,parameterName, 'Basin', 'Potomac-Shenandoah',NA, specialWeight=FALSE) ,
    estimateRappahannockYork = subpopEstimate(finalData,parameterName, 'Basin', 'Rappahannock-York',NA, specialWeight=FALSE),
    estimateNew = subpopEstimate(finalData,parameterName, 'Basin', 'New',NA, specialWeight=FALSE) ,
    estimateChowan = subpopEstimate(finalData,parameterName, 'Basin', 'Chowan',NA, specialWeight=FALSE), 
    estimateTennessee = subpopEstimate(finalData,parameterName, 'Basin', 'Tennessee',NA, specialWeight=FALSE), 
    estimateHolston = subpopEstimate(finalData,parameterName, 'SubBasin', 'Holston',NA, specialWeight=FALSE) ,
    estimateBigSandy = subpopEstimate(finalData,parameterName, 'SubBasin', 'Big Sandy',NA, specialWeight=FALSE) ,
    estimateClinchPowell = subpopEstimate(finalData,parameterName, 'SubBasin', 'Clinch-Powell',NA, specialWeight=FALSE) ,
    estimatePotomac = subpopEstimate(finalData,parameterName, 'SubBasin', 'Potomac',NA, specialWeight=FALSE) ,
    estimateShenandoah = subpopEstimate(finalData,parameterName, 'SubBasin', 'Shenandoah',NA, specialWeight=FALSE) ,
    estimateRappahannock = subpopEstimate(finalData,parameterName, 'SubBasin', 'Rappahannock',NA, specialWeight=FALSE) ,
    estimateYork = subpopEstimate(finalData,parameterName, 'SubBasin', 'York',NA, specialWeight=FALSE) ,
    
    # change here each new cycle
    # By VAHUSB, using abbreviations to make the category names cleaner here
    estimateRU = subpopEstimate(finalData,parameterName, 'VAHUSB', 'RU','Roanoke River, Upper', specialWeight=FALSE),
    estimateJM = subpopEstimate(finalData,parameterName, 'VAHUSB', 'JM','James River, Middle (Piedmont)', specialWeight=FALSE),
    estimateNE = subpopEstimate(finalData,parameterName, 'VAHUSB', 'NE','New River', specialWeight=FALSE),
    estimateRD = subpopEstimate(finalData,parameterName, 'VAHUSB', 'RD','Roanoke River- Dan River', specialWeight=FALSE),
    estimateYO = subpopEstimate(finalData,parameterName, 'VAHUSB', 'YO','York River', specialWeight=FALSE),
    estimatePL = subpopEstimate(finalData,parameterName, 'VAHUSB', 'PL','Potomac River, Lower', specialWeight=FALSE),
    estimateCU = subpopEstimate(finalData,parameterName, 'VAHUSB', 'CU','Chowan River, Upper', specialWeight=FALSE),
    estimateRA = subpopEstimate(finalData,parameterName, 'VAHUSB', 'RA','Rappahannock River', specialWeight=FALSE),
    estimateBS = subpopEstimate(finalData,parameterName, 'VAHUSB', 'BS','Big Sandy River', specialWeight=FALSE),
    estimateTH = subpopEstimate(finalData,parameterName, 'VAHUSB', 'TH','Tennessee-Holston River', specialWeight=FALSE),
    estimatePS = subpopEstimate(finalData,parameterName, 'VAHUSB', 'PS','Potomac River-Shenandoah River', specialWeight=FALSE),
    estimateJA = subpopEstimate(finalData,parameterName, 'VAHUSB', 'JA','James River- Appomattox River', specialWeight=FALSE),
    estimateCM = subpopEstimate(finalData,parameterName, 'VAHUSB', 'CM','Chowan River-Meherrin River', specialWeight=FALSE),
    estimateTC = subpopEstimate(finalData,parameterName, 'VAHUSB', 'TC','Tennessee-Clinch River', specialWeight=FALSE),
    
    
    # By Ecoregion
    estimatePiedmont = subpopEstimate(finalData,parameterName, 'EcoRegion', 'Piedmont',NA, specialWeight=FALSE) ,
    estimateNorthernPiedmont = subpopEstimate(finalData,parameterName, 'EcoRegion', 'Northern Piedmont',NA, specialWeight=FALSE) ,
    estimateCARV = subpopEstimate(finalData,parameterName, 'EcoRegion', 'Central Appalachian Ridges and Valleys',NA, specialWeight=FALSE) ,
    estimateSEplains = subpopEstimate(finalData,parameterName, 'EcoRegion', 'Southeastern Plains',NA, specialWeight=FALSE) ,
    estimateBRM = subpopEstimate(finalData,parameterName, 'EcoRegion', 'Blue Ridge Mountains',NA, specialWeight=FALSE) ,
    #estimateMACP = subpopEstimate(finalData,parameterName, 'EcoRegion', 'Middle Atlantic Coastal Plain',NA, specialWeight=FALSE) , # need 25-30 samples to run this analysis
    estimateCentralApps = subpopEstimate(finalData,parameterName, 'EcoRegion', 'Central Appalachians',NA, specialWeight=FALSE) ,
    # By Bioregion
    estimateMountainBioregion = subpopEstimate(finalData,parameterName, 'BioRegion', 'Mountain','Mountain Bioregion', specialWeight=FALSE),
    estimatePiedmontBioregion = subpopEstimate(finalData,parameterName, 'BioRegion', 'Piedmont','Piedmont Bioregion', specialWeight=FALSE) ,
    estimateCoastBioregion = subpopEstimate(finalData,parameterName, 'BioRegion', 'Coast', 'Coast Bioregion', specialWeight=FALSE),
    # By Order
    estimateFirstOrder = subpopEstimate(finalData,parameterName, 'Order', '1','First Order', specialWeight=FALSE) ,
    estimateSecondOrder = subpopEstimate(finalData,parameterName, 'Order', '2','Second Order', specialWeight=FALSE) ,
    estimateThirdOrder = subpopEstimate(finalData,parameterName, 'Order', '3','Third Order', specialWeight=FALSE) ,
    estimateFourthOrder = subpopEstimate(finalData,parameterName, 'Order', '4','Fourth Order', specialWeight=FALSE) ,
    estimateFifthOrder = subpopEstimate(finalData,parameterName, 'Order', '5','Fifth Order', specialWeight=FALSE) ,
    # By Basin Size
    estimateBasin1 = subpopEstimate(finalData,parameterName, 'BasinSize', '1','<1 square mile', specialWeight=FALSE) ,
    estimateBasin2 = subpopEstimate(finalData,parameterName, 'BasinSize', '2','1 to 10 square mile', specialWeight=FALSE) ,
    estimateBasin3 = subpopEstimate(finalData,parameterName, 'BasinSize', '3','10 to 50 square mile', specialWeight=FALSE) ,
    estimateBasin4 = subpopEstimate(finalData,parameterName, 'BasinSize', '4','>50 square mile', specialWeight=FALSE) ,
    # By Year
    estimate2001 = subpopEstimate(finalData,parameterName, 'Year', '2001','Year 2001', specialWeight='finalweight_Year') ,
    estimate2002 = subpopEstimate(finalData,parameterName, 'Year', '2002','Year 2002', specialWeight='finalweight_Year'), 
    estimate2003 = subpopEstimate(finalData,parameterName, 'Year', '2003','Year 2003', specialWeight='finalweight_Year') ,
    estimate2004 = subpopEstimate(finalData,parameterName, 'Year', '2004','Year 2004', specialWeight='finalweight_Year') ,
    estimate2005 = subpopEstimate(finalData,parameterName, 'Year', '2005','Year 2005', specialWeight='finalweight_Year') ,
    estimate2006 = subpopEstimate(finalData,parameterName, 'Year', '2006','Year 2006', specialWeight='finalweight_Year') ,
    estimate2007 = subpopEstimate(finalData,parameterName, 'Year', '2007','Year 2007', specialWeight='finalweight_Year') ,
    estimate2008 = subpopEstimate(finalData,parameterName, 'Year', '2008','Year 2008', specialWeight='finalweight_Year') ,
    estimate2009 = subpopEstimate(finalData,parameterName, 'Year', '2009','Year 2009', specialWeight='finalweight_Year') ,
    estimate2010 = subpopEstimate(finalData,parameterName, 'Year', '2010','Year 2010', specialWeight='finalweight_Year') ,
    estimate2011 = subpopEstimate(finalData,parameterName, 'Year', '2011','Year 2011', specialWeight='finalweight_Year') ,
    estimate2012 = subpopEstimate(finalData,parameterName, 'Year', '2012','Year 2012', specialWeight='finalweight_Year') ,
    estimate2013 = subpopEstimate(finalData,parameterName, 'Year', '2013','Year 2013', specialWeight='finalweight_Year') ,
    estimate2014 = subpopEstimate(finalData,parameterName, 'Year', '2014','Year 2014', specialWeight='finalweight_Year') ,
    estimate2015 = subpopEstimate(finalData,parameterName, 'Year', '2015','Year 2015', specialWeight='finalweight_Year') ,
    estimate2016 = subpopEstimate(finalData,parameterName, 'Year', '2016','Year 2016', specialWeight='finalweight_Year') ,
    estimate2017 = subpopEstimate(finalData,parameterName, 'Year', '2017','Year 2017', specialWeight='finalweight_Year') ,
    estimate2018 = subpopEstimate(finalData,parameterName, 'Year', '2018','Year 2018', specialWeight='finalweight_Year') ,
    # change here each new cycle
    estimate2019 = subpopEstimate(finalData,parameterName, 'Year', '2019','Year 2019', specialWeight='finalweight_Year') ,
    estimate2020 = subpopEstimate(finalData,parameterName, 'Year', '2020','Year 2020', specialWeight='finalweight_Year') ,
    
    # change here each new cycle
    
    # By Phase
    estimatePhase1 = subpopEstimate(finalData,parameterName, 'Panel', 'Phase1','Phase One 2001-2010', specialWeight='finalweight_Panel1') ,
    estimatePhase2 = subpopEstimate(finalData,parameterName, 'Panel', 'Phase2','Phase Two 2011-2020', specialWeight='finalweight_Panel2') ,
    
    
    # change here each new cycle
    # Bay/ NonBay
    estimateBay = subpopEstimate(finalData,parameterName, 'BayShed', 'Bay','Bay Watersheds 2001-2020', specialWeight=FALSE),
    estimateNonBay = subpopEstimate(finalData,parameterName, 'BayShed', 'NonBay','Non-Bay Watersheds 2001-2020', specialWeight=FALSE),
    estimateBayPhase1 = subpopEstimate(finalData,parameterName, 'BayPanel', 'BayPhase1','Bay Watersheds 2001-2010', specialWeight='finalweight_Panel1'),
    estimateBayPhase2 = subpopEstimate(finalData,parameterName, 'BayPanel', 'BayPhase2','Bay Watersheds 2011-2020', specialWeight='finalweight_Panel2'),
    estimateNonBayPhase1 = subpopEstimate(finalData,parameterName, 'BayPanel', 'NonBayPhase1','Non-Bay Watersheds 2001-2010', specialWeight='finalweight_Panel1'),
    estimateNonBayPhase2 = subpopEstimate(finalData,parameterName, 'BayPanel', 'NonBayPhase2','Non-Bay Watersheds 2011-2020', specialWeight='finalweight_Panel2'),
    
    # change here each new cycle
    
    # Sampling Phase
    estimateBioPanelPhase1 = subpopEstimate(finalData,parameterName, 'BioPanel', 'Phase1','VSCI Scores 2001-2005', specialWeight='finalweight_BioPanel1'),
    estimateBioPanelPhase2 = subpopEstimate(finalData,parameterName, 'BioPanel', 'Phase2','VSCI Scores 2006-2010', specialWeight='finalweight_BioPanel2'),
    estimateBioPanelPhase3 = subpopEstimate(finalData,parameterName, 'BioPanel', 'Phase3','VSCI Scores 2011-2015', specialWeight='finalweight_BioPanel3'),
    estimateBioPanelPhase4 = subpopEstimate(finalData,parameterName, 'BioPanel', 'Phase4','VSCI Scores 2016-2020', specialWeight='finalweight_BioPanel4'),
    # IR window
    estimateIR2008 = subpopEstimate(finalData,parameterName, 'IR2008', '2008','IR2008', specialWeight='finalweight_IR2008'),
    estimateIR2010 = subpopEstimate(finalData,parameterName, 'IR2010', '2010','IR2010', specialWeight='finalweight_IR2010'),
    estimateIR2012= subpopEstimate(finalData,parameterName, 'IR2012', '2012','IR2012', specialWeight='finalweight_IR2012'),
    estimateIR2014 = subpopEstimate(finalData,parameterName, 'IR2014', '2014','IR2014', specialWeight='finalweight_IR2014'),
    estimateIR2016 = subpopEstimate(finalData,parameterName, 'IR2016', '2016','IR2016', specialWeight='finalweight_IR2016'),
    estimateIR2018 = subpopEstimate(finalData,parameterName, 'IR2018', '2018','IR2018', specialWeight='finalweight_IR2018'),
    estimateIR2020 = subpopEstimate(finalData,parameterName, 'IR2020', '2020','IR2020', specialWeight='finalweight_IR2020'),
    # change here each new cycle
    estimateIR2022 = subpopEstimate(finalData,parameterName, 'IR2022', '2022','IR2022', specialWeight='finalweight_IR2022'),
    
    
    # Stream size
    estimateSmall = subpopEstimate(finalData,parameterName, 'StreamSizeCat','Small',NA, specialWeight=FALSE),
    estimateMedium = subpopEstimate(finalData,parameterName, 'StreamSizeCat','Medium',NA, specialWeight=FALSE),
    estimateLarge = subpopEstimate(finalData,parameterName, 'StreamSizeCat','Large',NA, specialWeight=FALSE),
    # Stream size and sample phase
    estimatePhase1Small = subpopEstimate(finalData,parameterName, 'StreamSizeCatPhase','Phase1Small',NA, specialWeight='finalweight_Panel1'),
    estimatePhase2Small = subpopEstimate(finalData,parameterName, 'StreamSizeCatPhase','Phase2Small',NA, specialWeight='finalweight_Panel2'),
    estimatePhase1Medium = subpopEstimate(finalData,parameterName, 'StreamSizeCatPhase','Phase1Medium',NA, specialWeight='finalweight_Panel1'),
    estimatePhase2Medium = subpopEstimate(finalData,parameterName, 'StreamSizeCatPhase','Phase2Medium',NA, specialWeight='finalweight_Panel2'),
    estimatePhase1Large = subpopEstimate(finalData,parameterName, 'StreamSizeCatPhase','Phase1Large',NA, specialWeight='finalweight_Panel1'),
    estimatePhase2Large = subpopEstimate(finalData,parameterName, 'StreamSizeCatPhase','Phase2Large',NA, specialWeight='finalweight_Panel2')
  )
}

allCDFresults <- function(designStatus, surveyData, parameter){
  
  # Organize new parameter data
  parameterData <- data.frame(select(surveyData,siteID,parameter))
  parameterData[,2] <- as.numeric(as.character(parameterData[,2])) # change to numeric in case it was saved as anything else
  
  
  # Change status if particular parameter wasnt sampled
  initialWeightData <- left_join(designStatus,parameterData,by='siteID') %>%
    mutate(parameter_status=ifelse(status == "TS" & is.na(!!as.name(parameter)),"OT", # if TS but not sampled, change to OT
                                   ifelse(status == "OT" & is.na(!!as.name(parameter)),"OT", # if OT but not sampled, keep OT
                                          ifelse(status == "PD", "PD", # if PD, keep it
                                                 ifelse(status == 'NT','NT', # if NT, keep it
                                                        ifelse(status == "TS",'TS', NA))))), # if TS AND sampled, keep as TS
           designweightoriginal = as.factor(`strahler order`), # easier to change as factor
           designweightoriginal = dplyr::recode(designweightoriginal,`1`="3790.5165999999999",
                                                `2`="947.62919999999997",
                                                `3`="541.50239999999997",
                                                `4`="315.87639999999999", 
                                                `5`="140.3895", 
                                                `6`="140.3895")) # overwrite all design weights back to original
  
  # Adjust initial weights for trend stations across various data windows
  # First calculate number of times sampled in each window
  
  trendWeightAdjustments <- mutate(initialWeightData,
                                   designweightoriginal = as.numeric(as.character(designweightoriginal)), #factor to numeric
                                   siteIDoriginal=gsub("_.*$", "", siteID)) %>% # get rid of concatenated year for trends to make calculations easier
    # Full window
    group_by(siteIDoriginal) %>% 
    mutate(nYearsSampled = ifelse(is.na(siteIDoriginal),NA,n())) %>%
    ungroup()%>%
    
    # change here each new cycle
    # 2022 IR
    group_by(siteIDoriginal, IR2022) %>%
    mutate(nYearsSampled_IR2022 = ifelse(is.na(IR2022),NA,n())) %>%
    ungroup() %>%
    
    
    # 2020 IR
    group_by(siteIDoriginal, IR2020) %>%
    mutate(nYearsSampled_IR2020 = ifelse(is.na(IR2020),NA,n())) %>%
    ungroup() %>%
    # 2018 IR
    group_by(siteIDoriginal, IR2018) %>%
    mutate(nYearsSampled_IR2018 = ifelse(is.na(IR2018),NA,n())) %>%
    ungroup() %>%
    # 2016 IR
    group_by(siteIDoriginal, IR2016) %>%
    mutate(nYearsSampled_IR2016 = ifelse(is.na(IR2016),NA,n())) %>%
    ungroup()%>%
    # 2014 IR
    group_by(siteIDoriginal, IR2014) %>%
    mutate(nYearsSampled_IR2014 = ifelse(is.na(IR2014),NA,n()))%>%
    ungroup()%>%
    # 2012 IR
    group_by(siteIDoriginal, IR2012) %>%
    mutate(nYearsSampled_IR2012 = ifelse(is.na(IR2012),NA,n())) %>%
    ungroup()%>%
    # 2010 IR
    group_by(siteIDoriginal, IR2010) %>%
    mutate(nYearsSampled_IR2010 = ifelse(is.na(IR2010),NA,n())) %>%
    ungroup()%>%
    # 2008 IR
    group_by(siteIDoriginal, IR2008) %>%
    mutate(nYearsSampled_IR2008 = ifelse(is.na(IR2008),NA,n())) %>%
    # Panels
    group_by(siteIDoriginal, Panel1) %>%
    mutate(nYearsSampled_Panel1 = ifelse(is.na(Panel1),NA,n())) %>%
    ungroup()%>%
    group_by(siteIDoriginal, Panel2) %>%
    mutate(nYearsSampled_Panel2 = ifelse(is.na(Panel2),NA,n())) %>%
    ungroup()%>%
    # Biopanels
    group_by(siteIDoriginal, BioPanel1) %>%
    mutate(nYearsSampled_BioPanel1 = ifelse(is.na(BioPanel1),NA,n())) %>%
    ungroup() %>% 
    group_by(siteIDoriginal, BioPanel2) %>%
    mutate(nYearsSampled_BioPanel2 = ifelse(is.na(BioPanel2),NA,n())) %>%
    ungroup() %>% 
    group_by(siteIDoriginal, BioPanel3) %>%
    mutate(nYearsSampled_BioPanel3 = ifelse(is.na(BioPanel3),NA,n())) %>%
    ungroup() %>% 
    group_by(siteIDoriginal, BioPanel4) %>%
    mutate(nYearsSampled_BioPanel4 = ifelse(is.na(BioPanel4),NA,n())) %>%
    ungroup() %>% 
    # Divide out weight by years sampled in each window
    # if the site wasnt sampled in a window, give it the original weight but it will be adjusted according to an OT status later
    mutate(designweight_all= designweightoriginal/nYearsSampled,
           # change here each new cycle
           designweight_IR2022= ifelse(is.na(nYearsSampled_IR2022),designweightoriginal,designweightoriginal/nYearsSampled_IR2022),
           
           designweight_IR2020= ifelse(is.na(nYearsSampled_IR2020),designweightoriginal,designweightoriginal/nYearsSampled_IR2020),
           designweight_IR2018= ifelse(is.na(nYearsSampled_IR2018),designweightoriginal,designweightoriginal/nYearsSampled_IR2018),
           designweight_IR2016= ifelse(is.na(nYearsSampled_IR2016),designweightoriginal,designweightoriginal/nYearsSampled_IR2016),
           designweight_IR2014= ifelse(is.na(nYearsSampled_IR2014),designweightoriginal,designweightoriginal/nYearsSampled_IR2014),
           designweight_IR2012= ifelse(is.na(nYearsSampled_IR2012),designweightoriginal,designweightoriginal/nYearsSampled_IR2012),
           designweight_IR2010= ifelse(is.na(nYearsSampled_IR2010),designweightoriginal,designweightoriginal/nYearsSampled_IR2010),
           designweight_IR2008= ifelse(is.na(nYearsSampled_IR2008),designweightoriginal,designweightoriginal/nYearsSampled_IR2008),
           designweight_Panel1= ifelse(is.na(nYearsSampled_Panel1),designweightoriginal,designweightoriginal/nYearsSampled_Panel1),
           designweight_Panel2= ifelse(is.na(nYearsSampled_Panel2),designweightoriginal,designweightoriginal/nYearsSampled_Panel2),
           designweight_BioPanel1= ifelse(is.na(nYearsSampled_BioPanel1),designweightoriginal,designweightoriginal/nYearsSampled_BioPanel1),
           designweight_BioPanel2= ifelse(is.na(nYearsSampled_BioPanel2),designweightoriginal,designweightoriginal/nYearsSampled_BioPanel2),
           designweight_BioPanel3= ifelse(is.na(nYearsSampled_BioPanel3),designweightoriginal,designweightoriginal/nYearsSampled_BioPanel3),
           designweight_BioPanel4= ifelse(is.na(nYearsSampled_BioPanel4),designweightoriginal,designweightoriginal/nYearsSampled_BioPanel4))
  
  
  ## Adjust design weights to get final weights
  
  # Initial sample frame inputs
  # List stream order by kilometer it represents
  sframe <- c('1st'=51210, '2nd'=13680, '3rd'=7781.08, '4th'=4448.257, 
              '5th'=1731.302, '6th'=163.901, '7th'=14.7099 )
  
  # recode to factor to make sframe match up to stream order
  trendWeightAdjustments$`strahler order` <- as.factor(trendWeightAdjustments$`strahler order`)
  levels(trendWeightAdjustments$`strahler order`) <- c('1st','2nd','3rd','4th','5th','6th','7th')
  
  finalWeights <- select(trendWeightAdjustments,
                         siteID:IR2022, # change here each new cycle
                         siteIDoriginal,designweightoriginal,
                         designweight_all:designweight_BioPanel4,
                         parameter_status,!!as.name(parameter))
  finalWeights$finalweight_Year <- adjwgt(rep(TRUE,length(finalWeights$designweightoriginal)),finalWeights$designweightoriginal,
                                          finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_all <- adjwgt(rep(TRUE,length(finalWeights$designweight_all)),finalWeights$designweight_all,
                                         finalWeights$`strahler order`, sframe)
  
  # change here each new cycle
  finalWeights$finalweight_IR2022 <- adjwgt(rep(TRUE,length(finalWeights$designweight_IR2022)),finalWeights$designweight_IR2022,
                                            finalWeights$`strahler order`, sframe)
  
  finalWeights$finalweight_IR2020 <- adjwgt(rep(TRUE,length(finalWeights$designweight_IR2020)),finalWeights$designweight_IR2020,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_IR2018 <- adjwgt(rep(TRUE,length(finalWeights$designweight_IR2018)),finalWeights$designweight_IR2018,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_IR2016 <- adjwgt(rep(TRUE,length(finalWeights$designweight_IR2016)),finalWeights$designweight_IR2016,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_IR2014 <- adjwgt(rep(TRUE,length(finalWeights$designweight_IR2014)),finalWeights$designweight_IR2014,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_IR2012 <- adjwgt(rep(TRUE,length(finalWeights$designweight_IR2012)),finalWeights$designweight_IR2012,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_IR2010 <- adjwgt(rep(TRUE,length(finalWeights$designweight_IR2010)),finalWeights$designweight_IR2010,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_IR2008 <- adjwgt(rep(TRUE,length(finalWeights$designweight_IR2008)),finalWeights$designweight_IR2008,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_Panel1 <- adjwgt(rep(TRUE,length(finalWeights$designweight_Panel1)),finalWeights$designweight_Panel1,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_Panel2 <- adjwgt(rep(TRUE,length(finalWeights$designweight_Panel2)),finalWeights$designweight_Panel2,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_BioPanel1 <- adjwgt(rep(TRUE,length(finalWeights$designweight_BioPanel1)),finalWeights$designweight_BioPanel1,
                                               finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_BioPanel2 <- adjwgt(rep(TRUE,length(finalWeights$designweight_BioPanel2)),finalWeights$designweight_BioPanel2,
                                               finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_BioPanel3 <- adjwgt(rep(TRUE,length(finalWeights$designweight_BioPanel3)),finalWeights$designweight_BioPanel3,
                                               finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_BioPanel4 <- adjwgt(rep(TRUE,length(finalWeights$designweight_BioPanel4)),finalWeights$designweight_BioPanel4,
                                               finalWeights$`strahler order`, sframe)
  
  
  
  # Change decimal degree coordinates to marinus for equal area projection (needed for spsurvey::cat.analysis())
  marinus <- spsurvey::marinus(finalWeights$`Latitude-DD`,finalWeights$`Longitude-DD`)
  finalWeights <- cbind(finalWeights,data.frame(xmarinus=marinus$x,ymarinus=marinus$y))
  rm(marinus)
  
  ####Stream Extent and Status Estimate
  
  siteExtent <- data.frame(siteID=finalWeights$siteID, Use=rep(TRUE,nrow(finalWeights)) )
  subpopExtent <- data.frame(siteID=finalWeights$siteID,Region=rep('Virginia',nrow(finalWeights)) )
  designExtent <- data.frame(siteID=finalWeights$siteID,stratum=rep(1,nrow(finalWeights)),
                             wgt=finalWeights$finalweight_all, xcoord=finalWeights$xmarinus,ycoord=finalWeights$ymarinus)
  
  StatusTNT <- finalWeights$status
  levels(StatusTNT) <- list(T=c('TS','PD','OT'), NT=c('NT') )
  
  data.cat.ext <- data.frame(finalWeights[,c('siteID','status')],StatusTNT)
  
  popstatus.est <- cat.analysis(sites = siteExtent, subpop = subpopExtent, design = designExtent,
                                data.cat = data.cat.ext, conf=95, vartype="Local")
  
  dataAnalyzed <- filter(finalWeights, parameter_status=='TS') %>%
    select(siteID,!!as.name(parameter),parameter_status,everything())
  
  output <- listOfResults(popstatus.est, dataAnalyzed, parameter)
  
  return(output)
}