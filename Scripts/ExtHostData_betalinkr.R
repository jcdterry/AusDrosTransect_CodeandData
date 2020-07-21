### Updated version of betalinkr to better assign interactions to shared or not. 

## Please note I have stripped out most of the defensive programming and options from the function. 
## and 'hard-coded' the options used in the paper. 

# i.e:

# proportions = FALSE
# index = "jaccard"
# binary = TRUE
# partitioning = "commondenom"
# distofempty = 'zero'
# function.dist = "vegdist"

# Example inputs: : 
###   webarray<- webs2array(QualWebs[c(1,2)])

# 
# FullData %>% 
#   count(Site, Host)%>%
#   filter(!is.na(Host))%>%
#   select(-n)%>%
#   mutate(TRUE)-> AllHostRecords

## AllHostRecordsmust be dataframe of the host occurances that need to be added. col 1 = Sites. col2 = Host. Like this:

#   Site  Host            
#  <chr> <chr>           
# 1 K070  bipectinata     
# 2 K070  bunnanda        
# 3 K070  Drosophila sp. 1
# 4 K070  pallidifrons     
# 5 K070  pandora         
# 6 K070  pseudoananassae 
# 7 K070  rubida          
# 8 K070  sulfurigaster   


# betalinkr_multi_Exthostinfo() # is a wrapper that repeatedly calls the function to populate a data.frame

betalinkr_Exthostinfo <-function (webarray,
                                  AllHostRecords) {
  
  WebNames<- dimnames(webarray)[[3]]
  
  
  Web1Pres<- dimnames(webarray)[[1]] %in%   AllHostRecords$Host[AllHostRecords$Site==WebNames[1]]
  Web2Pres<- dimnames(webarray)[[1]] %in%   AllHostRecords$Host[AllHostRecords$Site==WebNames[2]]
  
  webarray <- webarray[Web1Pres|Web2Pres, # This keeps any hosts that appear in at least one of the networks
                       apply(webarray, 2, sum) > 0, ,
                       drop = FALSE]
  
  Web1Pres<- dimnames(webarray)[[1]] %in%   AllHostRecords$Host[AllHostRecords$Site==WebNames[1]] # Reassess if present or not
  Web2Pres<- dimnames(webarray)[[1]] %in%   AllHostRecords$Host[AllHostRecords$Site==WebNames[2]] # Makes logic the same length
  
  
  
  ## 
  array.sharedsp <- webarray
  linkmx <- array2linkmx(webarray) # Base link matrix
  
  
  ###  Set to zero any that are not shared, then make into a link matrix of shared species (i.e.those not found in both)
  #array.sharedsp[rowSums(apply(webarray, MARGIN = c(1, 3),  sum) > 0) != 2, , ] <- 0 
  array.sharedsp[   !(Web1Pres&Web2Pres), , ] <- 0 
  array.sharedsp[, rowSums(apply(webarray, MARGIN = c(2, 3),sum) > 0) != 2, ] <- 0   # Paras as per original function
  linkmx.sharedsp <- array2linkmx(array.sharedsp)
  
  
  specmx.all <- (cbind( rbind(Web1Pres,Web2Pres)[,Web1Pres|Web2Pres], # host presence from lists not interaction web
                       apply(webarray, c(3, 2), sum))>0)*1
  # print(specmx.all)
  b_s      <- vegdist(specmx.all,    method = "jaccard", binary = TRUE)   # dissimilarity in species lists # not used further
  b_s_host <- vegdist(rbind(Web1Pres,Web2Pres)[,Web1Pres|Web2Pres],    method = "jaccard", binary = TRUE)   # dissimilarity in host species lists

  linkmx.sharedli <- linkmx
  linkmx.sharedli[, colSums(linkmx.sharedli > 0) == 1] <- 0
  linkmx.rewiring <- linkmx.sharedsp - linkmx.sharedli   # Shared species but not shared link
  linkmx.RewSha   <- linkmx.rewiring + linkmx.sharedli   # Rewiring and shared
  linkmx.uniquesp <- linkmx - linkmx.sharedsp            # Not shared species
  linkmx.UniSha   <- linkmx.uniquesp + linkmx.sharedli   # 
  
  b_wn     <- vegdist(linkmx,        method = "jaccard", binary = TRUE)   # whole network 
  b_os.raw <- vegdist(linkmx.RewSha, method = "jaccard", binary = TRUE)   # interaction change
  b_st.raw <- vegdist(linkmx.UniSha, method = "jaccard", binary = TRUE)   # species turnover
  
  b_os <- b_os.raw * b_wn/(b_os.raw + b_st.raw)
  b_st <- b_st.raw * b_wn/(b_os.raw + b_st.raw)
  return(c(S = b_s, OS = b_os, WN = b_wn, ST = b_st, S_host = b_s_host))
}

betalinkr_multi_Exthostinfo <-  function(webarray, AllHostRecords){
  data.out <- data.frame(i = integer(0), j = integer(0), S = numeric(0), 
                         OS = numeric(0), WN = numeric(0), ST = numeric(0), S_host = numeric(0))
  webnames <- dimnames(webarray)[[3]]
  if (is.null(webnames)) 
    webnames <- 1:dim(webarray)[3]
  for (i in 1:(dim(webarray)[3] - 1)) {
    for (j in (i + 1):dim(webarray)[3]) {
      data.out[nrow(data.out) + 1, ] <- c(webnames[i], 
                                          webnames[j],
                                          betalinkr_Exthostinfo(webarray = webarray[,, c(i, j)], 
                                                    AllHostRecords))
    }
  }
  return(data.out)
}




