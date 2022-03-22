# flashfm-ivis version_2022_03_10


####----Libraries----
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
library(dplyr)
library(ggplot2)
library(shiny)
library(shinydashboard) 
library(plotly)
library(DT)
library(gaston)
#library(RColorBrewer)
#display.brewer.all(type = 'qual')
library(randomcoloR)
library(igraph)
library(networkD3)
library(stringr)
library(shinybusy)
# remotes::install_github("rstudio/webshot2")
# library(webshot2)



####---- Users' input files and guideline----
#1. Finemapwithflashfm:Need
#---> .z files: GWAS
#use trait as file name, one z file for each trait to get info for Regional association plots
#---> .ld file: LD matrix
#orginal fld="finemap.ld", it will be r^2 in the internal wrapper function
#or check whether there is negative value in the matrix, if not, then already r^2 is done
#---> .RData files: mpp.pp and snpGroups
#output from PPsummarise and makeSNPgroup2 – these are the mpp.pp and snpGroups objects 
#suggest that the user saves the mpp.pp and snpGroups objects to a single RData file, 
#but should give an example of the command to do this: save(mpp.pp,snpGroups,file="anyfilename.RData")
#Processing of the files is in a wrapper function at the end of the finemap-example.R file.
#2. FLASHFMwithJAM:Need
#---> flashfmwithjam.RData (this is just the output from FLASHFMwithJAM)
#---> snpinfo.RData (this will have to be given as an example and the user has to construct it
#from their data. We should fix the first 4 names and allow the names for the pvalue columns to 
#be flexible, but restrict that pvalue columns have to be given in same order as the betas 
#submitted to flashfm. There could be 2-6 columns for p-values.


####----pre-loaded data as examples----
#users will overwrite it if they upload new files
#load("data/R_data.RData")


####----wrapper functions----
#1. Finemapwithflashfm:Need
# below fz is a vector of z files for each trait; fld is ld file;
# e.g. fz <- c("fname1.z","fname2.z","fname3.z","fname4.z")
# mpp.pp,snpGroups are from the RData files I sent awhile ago - in current use for tool
# fz <- numeric(4); for(i in 1:4) fz[i] <- paste0("/net/bliss-04-nfs/scratch/jennifer/flashfm/scripts/Example/Uganda-APOE-LDL/finemap",i,".z") 
format.flashfm.finemap <- function(fz, fld, mpp.pp, snpGroups) {
  M <- length(fz)
  GWAS_0 <- vector("list",M)
  GWAS_1 <- vector("list",M)
  GWAS <- vector("list",M)
  for(i in 1:M) {
    GWAS_0[[i]] <- read.table(fz[i], header = TRUE, as.is = TRUE, sep = " ")
    GWAS_1[[i]]$rs <- GWAS_0[[i]]$rsid
    GWAS_1[[i]]$chr <- GWAS_0[[i]]$chromosome
    GWAS_1[[i]]$ps <- GWAS_0[[i]]$position
    GWAS_1[[i]]$allele1 <- GWAS_0[[i]]$allele1
    GWAS_1[[i]]$allele0 <- GWAS_0[[i]]$allele2
    GWAS_1[[i]]$beta <- GWAS_0[[i]]$beta
    GWAS_1[[i]]$se <- GWAS_0[[i]]$se
    #GWAS[[i]]$p <-  2*(1-pnorm(abs(GWAS[[i]]$beta/GWAS[[i]]$se)))
    GWAS_1[[i]]$pval <-  2*(pnorm(abs(GWAS_1[[i]]$beta/GWAS_1[[i]]$se),lower.tail = F)) #check!!!
    GWAS_1[[i]]$af <- 1-GWAS_0[[i]]$maf
    GWAS_1[[i]]$no <- GWAS_0[[i]]$se*0
    GWAS[[i]] = as.data.frame(GWAS_1[[i]])
  }
  r <- read.table(fld, header = FALSE, as.is = TRUE, sep = " ")   
  r2 <- r^2
  colnames(r2) <- GWAS_0[[1]]$rsid
  rownames(r2) <- GWAS_0[[1]]$rsid
  
  return(list(GWAS=GWAS,
              LD=r2,
              MPP=mpp.pp$MPP,
              MPPg=mpp.pp$MPPg,
              PP=mpp.pp$PP,
              PPg=mpp.pp$PPg,
              snpGroups=snpGroups))
}
#1. Finemapwithflashfm:Need
#inputs_used = format.flashfm.finemap(fz=fz, fld=fld, mpp.pp=mpp.pp, snpGroups = snpGroups)
#2. FLASHFMwithJAM:Need 
#inputs_used = format.flashfm.finemap(fz=fz, fld=fld, mpp.pp=fm$mpp.pp, snpGroups = fm$snpGroups)
#3. Pre-loaded data (need match the variable names)
#inputs_used = format.flashfm.finemap(...)
######


## snp group network: ppg_network_links and nodes function----
ppg_network_links_nodes <- function(mpp.pp) {
  M = length(mpp.pp$PPg)
  ppg_network_links_fm_0 = c()
  ppg_network_links_flashfm_0 = c()
  ppg_network_nodes_0 = c()
  for (m in 1:M){
    ppg_rownames = rownames(mpp.pp$PPg[[m]])
    ppg_rownames = str_split(ppg_rownames, "%", simplify = TRUE)
    ppg_rownames[ppg_rownames==""] <- 'NA'
    #fm----
    ppg_rownames_link_fm = c()
    for (i in 1:dim(ppg_rownames)[1]){
      if (length(ppg_rownames[i,][ppg_rownames[i,] != 'NA'])==1){
        ppg_rownames[i,2] = ppg_rownames[i,1]
      }
      tmp = combn(ppg_rownames[i,][ppg_rownames[i,] != 'NA'], 2)
      tmp = as.data.frame(t(tmp))
      tmp$level = i
      tmp$value = mpp.pp$PPg[[m]][i,1]
      tmp$trait_linecolor = names(mpp.pp$PPg)[m]
      # if(m==1){tmp$color = "chocolate"}
      # if(m==2){tmp$color = "pink"}
      # if(m==3){tmp$color = "lightgreen"}
      # if(m==4){tmp$color = "dodgerblue"}
      # if(m==5){tmp$color = "darkred"}
      # if(m==6){tmp$color = "gold"}
      if(m==1){tmp$color = "darkorange"}
      if(m==2){tmp$color = "forestgreen"}
      if(m==3){tmp$color = "lightgreen"}
      if(m==4){tmp$color = "pink"}
      if(m==5){tmp$color = "darkred"} ##CHECK COLOR WHEN RUNNING 6 TRAITS
      if(m==6){tmp$color = "gold"}    ##CHECK COLOR WHEN RUNNING 6 TRAITS
      ppg_rownames_link_fm = rbind(ppg_rownames_link_fm, tmp)
    }
    colnames(ppg_rownames_link_fm)[1:2] = c('source', 'target')
    ppg_network_links_fm_0 = rbind(ppg_network_links_fm_0, ppg_rownames_link_fm)
    #flashfm----
    ppg_rownames_link_flashfm = c()
    for (i in 1:dim(ppg_rownames)[1]){
      if (length(ppg_rownames[i,][ppg_rownames[i,] != 'NA'])==1){
        ppg_rownames[i,2] = ppg_rownames[i,1]
      }
      tmp = combn(ppg_rownames[i,][ppg_rownames[i,] != 'NA'], 2)
      tmp = as.data.frame(t(tmp))
      tmp$level = i
      tmp$value = mpp.pp$PPg[[m]][i,2]
      tmp$trait_linecolor = names(mpp.pp$PPg)[m]
      # if(m==1){tmp$color = "chocolate"}
      # if(m==2){tmp$color = "pink"}
      # if(m==3){tmp$color = "lightgreen"}
      # if(m==4){tmp$color = "dodgerblue"}
      # if(m==5){tmp$color = "darkred"}
      # if(m==6){tmp$color = "gold"}
      if(m==1){tmp$color = "darkorange"}
      if(m==2){tmp$color = "forestgreen"}
      if(m==3){tmp$color = "lightgreen"}
      if(m==4){tmp$color = "pink"}
      if(m==5){tmp$color = "darkred"} ##CHECK COLOR WHEN RUNNING 6 TRAITS
      if(m==6){tmp$color = "gold"}    ##CHECK COLOR WHEN RUNNING 6 TRAITS
      ppg_rownames_link_flashfm = rbind(ppg_rownames_link_flashfm, tmp)
    }
    colnames(ppg_rownames_link_flashfm)[1:2] = c('source', 'target')
    ppg_network_links_flashfm_0 = rbind(ppg_network_links_flashfm_0, ppg_rownames_link_flashfm)
    #nodes
    assign(paste0("ppg_rownames_",m),ppg_rownames)
    ppg_network_nodes_0 = c(ppg_network_nodes_0, ppg_rownames)
  }
  ppg_network_nodes_1 = unique(ppg_network_nodes_0)
  ppg_network_nodes_1 = ppg_network_nodes_1[!ppg_network_nodes_1=='NA']
  ppg_network_nodes_1 = as.data.frame(ppg_network_nodes_1)
  colnames(ppg_network_nodes_1) <- 'name'
  ppg_network_nodes_1$id = c(0:(length(ppg_network_nodes_1[,1])-1))
  ppg_network_nodes_3 = c()  
  for (m in 1:M){  
    ppg_network_nodes_2 = ppg_network_nodes_1
    ppg_network_nodes_2[, names(mpp.pp$PPg)[m]] = 0
    for (j in 1:length(ppg_network_nodes_2[,1])){
      ppg_network_nodes_2[, names(mpp.pp$PPg)[m]][j] = sum(ppg_network_nodes_2$name[j] == get(paste0("ppg_rownames_",m)))
    }
    ppg_network_nodes_3 = cbind(ppg_network_nodes_3, ppg_network_nodes_2[, names(mpp.pp$PPg)[m]])
  }
  ppg_network_nodes = cbind(ppg_network_nodes_1, ppg_network_nodes_3)
  colnames(ppg_network_nodes) <- c('name', 'id', names(mpp.pp$PPg))
  ppg_network_nodes$role_size_total = rowSums(ppg_network_nodes[,3:dim(ppg_network_nodes)[2]])
  ppg_network_nodes$group_color = 'trait:--'
  for (i in 1:length(ppg_network_nodes[,1])){
    for (j in 3:(2+M)){
      if (ppg_network_nodes[i,j]>0){
        tmp = colnames(ppg_network_nodes)[j]
        ppg_network_nodes$group_color[i] = paste(c(ppg_network_nodes$group_color[i],tmp), collapse="+")
      }
    }
  }
  #consistent with the original tables
  ppg_network_links_flashfm = ppg_network_links_flashfm_0[,c(3,1,2,4,5,6)]
  colnames(ppg_network_links_flashfm) = c('level', 'source', 'target', 
                                          "importance_linesize", "trait_linecolor", "color")
  ppg_network_links_fm = ppg_network_links_fm_0[,c(3,1,2,4,5,6)]
  colnames(ppg_network_links_fm) = c('level', 'source', 'target', 
                                     "importance_linesize", "trait_linecolor", "color")
  ppg_network_nodes$role_size = ppg_network_nodes$role_size_total
  ppg_network_nodes2 = ppg_network_nodes
  # create data: group fm
  plt_links_group_fm <- data.frame(
    source=ppg_network_links_fm$source,
    target=ppg_network_links_fm$target,
    l_color =ppg_network_links_fm$color,
    l_size  =ppg_network_links_fm$importance_linesize*10
  )
  plt_nodes_group_fm <- data.frame(
    name =ppg_network_nodes2$name,
    n_color=ppg_network_nodes2$group_color,
    n_size =ppg_network_nodes2$role_size^3   #updated20220208
  )
  for (i in 1:length(plt_links_group_fm$source)){
    plt_links_group_fm$source_id[i] = ppg_network_nodes2$id[which(plt_links_group_fm$source[i]==ppg_network_nodes2$name)]
    plt_links_group_fm$target_id[i] = ppg_network_nodes2$id[which(plt_links_group_fm$target[i]==ppg_network_nodes2$name)]
  }
  # create data: group flashfm
  plt_links_group_flashfm <- data.frame(
    source=ppg_network_links_flashfm$source,
    target=ppg_network_links_flashfm$target,
    l_color =ppg_network_links_flashfm$color,
    l_size  =ppg_network_links_flashfm$importance_linesize*10
  )
  plt_nodes_group_flashfm <- data.frame(
    name =ppg_network_nodes2$name,
    n_color=ppg_network_nodes2$group_color,
    n_size =ppg_network_nodes2$role_size^3    #updated20220208
  )
  for (i in 1:length(plt_links_group_flashfm$source)){
    plt_links_group_flashfm$source_id[i] = ppg_network_nodes2$id[which(plt_links_group_flashfm$source[i]==ppg_network_nodes2$name)]
    plt_links_group_flashfm$target_id[i] = ppg_network_nodes2$id[which(plt_links_group_flashfm$target[i]==ppg_network_nodes2$name)]
  }
  
  return(list(plt_links_group_fm = plt_links_group_fm,
              plt_nodes_group_fm = plt_nodes_group_fm,
              plt_links_group_flashfm = plt_links_group_flashfm,
              plt_nodes_group_flashfm = plt_nodes_group_flashfm))
}
# ##CHECK--IMPORTANT
# ppg_network_links_nodes = ppg_network_links_nodes(mpp.pp)
# plt_links_group_fm = ppg_network_links_nodes$plt_links_group_fm
# plt_nodes_group_fm = ppg_network_links_nodes$plt_nodes_group_fm
# plt_links_group_flashfm = ppg_network_links_nodes$plt_links_group_flashfm
# plt_nodes_group_flashfm = ppg_network_links_nodes$plt_nodes_group_flashfm
ColourScale_group <- 'd3.scaleOrdinal().range(["lightblue", "lightcoral", "lightgoldenrod", "lightgreen",
                                               "lightpink", "lightsalmon", "lightyellow", "lightslateblue"]);'
### ALL DONE: snp group network


## snp individual network: pp_network_links and nodes function----
#Step 1: define credible sets based on a specific significant level (i.e. 0.99)
credset <- function(modPP,cred=0.99) {
  tmp <- modPP[order(modPP, decreasing = TRUE)]
  cpp <- cumsum(tmp)
  wh <- which(cpp <= cred)
  if (!length(wh)) wh <- 1
  wh <- c(wh, max(wh) + 1)
  keepmodPP <- tmp[wh]
  mods <- names(keepmodPP)
  usnps <- unique(unlist(strsplit(mods,"%")))
  return(usnps)
}

# M <- length(mpp.pp$PP) # number of traits 
# cs1_selected <- csM_selected <- vector("list",M)
# for(i in 1:M){
#   tmp <- mpp.pp$PP[[i]][,1] # trait i, single-trait fine-mapping
#   cs1_selected[[i]] <- credset(tmp,0.99)
#   tmp <- mpp.pp$PP[[i]][,2] # trait i, multi-trait fine-mapping
#   csM_selected[[i]] <- credset(tmp,0.99)
# }

pp_network_links_nodes <- function(mpp.pp, cs1, csM) {
  fun_check_cs <- function(x,cs){
    if (x %in% cs){
      return(x)
    }else{
      return("NA")
    }
  }
  #pp_network_links_nodes <- function(mpp.pp, cs1, csM) {
  M = length(mpp.pp$PP)
  pp_rownames_cs1_link_all = c()
  pp_cs1_node_0 = c()
  pp_rownames_csM_link_all = c()
  pp_csM_node_0 = c()
  for (m in 1:M){
    pp_rownames = rownames(mpp.pp$PP[[m]])
    pp_rownames = str_split(pp_rownames, "%", simplify = TRUE)
    pp_rownames[pp_rownames==""] <- 'NA'
    pp_rownames_cs1 = apply(pp_rownames, 1:2, fun_check_cs, cs = cs1[[m]])
    #pp_rownames_cs1 = pp_rownames_cs1[rowSums(pp_rownames_cs1 == "NA")<dim(pp_rownames_cs1)[2]-1,]
    pp_rownames_csM = apply(pp_rownames, 1:2, fun_check_cs, cs = csM[[m]])
    #pp_rownames_csM = pp_rownames_csM[rowSums(pp_rownames_csM == "NA")<dim(pp_rownames_csM)[2]-1,]
    assign(paste0("pp_rownames_cs1_",m),pp_rownames_cs1)
    assign(paste0("pp_rownames_csM_",m),pp_rownames_csM)
    ##cs1
    pp_rownames_cs1_link = c()
    for (i in 1:dim(pp_rownames_cs1)[1]){
      # if (length(pp_rownames_cs1[i,][pp_rownames_cs1[i,] != 'NA'])==1){
      #   pp_rownames_cs1[i,2] = pp_rownames_cs1[i,1]
      # }
      if(sum(pp_rownames_cs1[i,] != 'NA')>1){
        tmp = combn(pp_rownames_cs1[i,][pp_rownames_cs1[i,] != 'NA'], 2)
        tmp = as.data.frame(t(tmp))
        tmp$level = i
        tmp$value = mpp.pp$PP[[m]][i,1] #check: flashfm will be mpp.pp$PP[[m]][i,2]
        tmp$trait_linecolor = names(mpp.pp$PP)[m]
        if(m==1){tmp$color = "lightgreen"}
        if(m==2){tmp$color = "pink"}
        if(m==3){tmp$color = "chocolate"}
        if(m==4){tmp$color = "orchid"}
        if(m==5){tmp$color = "darkred"}
        if(m==6){tmp$color = "gold"}   #check: our maximum is 6 traits
        pp_rownames_cs1_link = rbind(pp_rownames_cs1_link, tmp)
      }
    }
    colnames(pp_rownames_cs1_link)[1:2] = c('source', 'target')
    pp_rownames_cs1_link_all = rbind(pp_rownames_cs1_link_all,pp_rownames_cs1_link)
    #nodes
    pp_cs1_node_0 = c(pp_cs1_node_0, cs1[[m]])
    ##csM
    pp_rownames_csM_link = c()
    for (i in 1:dim(pp_rownames_csM)[1]){
      # if (length(pp_rownames_csM[i,][pp_rownames_csM[i,] != 'NA'])==1){
      #   pp_rownames_csM[i,2] = pp_rownames_csM[i,1]
      # }
      if (sum(pp_rownames_csM[i,] != 'NA')>1){
        tmp = combn(pp_rownames_csM[i,][pp_rownames_csM[i,] != 'NA'], 2)
        tmp = as.data.frame(t(tmp))
        tmp$level = i
        tmp$value = mpp.pp$PP[[m]][i,2] #check: flashfm will be mpp.pp$PP[[m]][i,2]
        tmp$trait_linecolor = names(mpp.pp$PP)[m]
        if(m==1){tmp$color = "lightgreen"}
        if(m==2){tmp$color = "pink"}
        if(m==3){tmp$color = "chocolate"}
        if(m==4){tmp$color = "orchid"}
        if(m==5){tmp$color = "darkred"}
        if(m==6){tmp$color = "gold"}   #check: our maximum is 6 traits
        pp_rownames_csM_link = rbind(pp_rownames_csM_link, tmp)
      }
    }
    colnames(pp_rownames_csM_link)[1:2] = c('source', 'target')
    pp_rownames_csM_link_all = rbind(pp_rownames_csM_link_all,pp_rownames_csM_link)
    #nodes
    pp_csM_node_0 = c(pp_csM_node_0, csM[[m]])
  } #checked.
  ##cs1
  pp_cs1_link = pp_rownames_cs1_link_all
  pp_cs1_node = pp_cs1_node_0 
  pp_cs1_node = as.data.frame(unique(pp_cs1_node))
  colnames(pp_cs1_node) <- 'name'
  pp_cs1_node$id = c(0:(length(pp_cs1_node[,1])-1))
  
  pp_cs1_node_sum = c()
  for (m in 1:M){
    pp_cs1_node[paste0('role_size_',m)] = 0
    for (i in 1:length(pp_cs1_node[,1])){
      pp_cs1_node[ ,paste0('role_size_',m)][i] = sum(pp_cs1_node$name[i] == get(paste0("pp_rownames_cs1_",m)))
    }
    pp_cs1_node_sum = cbind(pp_cs1_node_sum, pp_cs1_node[ ,paste0('role_size_',m)])
  }
  pp_cs1_node$role_size_all = rowSums(pp_cs1_node_sum)
  colnames(pp_cs1_node) <- c('name', 'id', names(mpp.pp$PP), 'role_size_total')
  pp_cs1_node$group_color = 'trait:--'
  for (i in 1:length(pp_cs1_node[,1])){
    for (j in 3:(2+length(names(mpp.pp$PP)))){
      if (pp_cs1_node[i,j]>0){
        tmp = colnames(pp_cs1_node)[j]
        pp_cs1_node$group_color[i] = paste(c(pp_cs1_node$group_color[i],tmp), collapse="+")
      }
    }
  }
  ### ALL DONE
  # create data:
  plt_links_snp <- data.frame(
    source=pp_cs1_link$source,
    target=pp_cs1_link$target,
    l_color =pp_cs1_link$color,
    l_size  =pp_cs1_link$value*10
  )
  plt_nodes_snp <- data.frame(
    name =pp_cs1_node$name,
    n_color=pp_cs1_node$group_color,
    n_size =pp_cs1_node$role_size_total
  )
  for (i in 1:length(plt_links_snp$source)){
    plt_links_snp$source_id[i] = pp_cs1_node$id[which(plt_links_snp$source[i]==pp_cs1_node$name)]
    plt_links_snp$target_id[i] = pp_cs1_node$id[which(plt_links_snp$target[i]==pp_cs1_node$name)]
  }
  ##csM
  pp_csM_link = pp_rownames_csM_link_all
  pp_csM_node = pp_csM_node_0 
  pp_csM_node = as.data.frame(unique(pp_csM_node))
  colnames(pp_csM_node) <- 'name'
  pp_csM_node$id = c(0:(length(pp_csM_node[,1])-1))
  pp_csM_node_sum = c()
  for (m in 1:M){
    pp_csM_node[paste0('role_size_',m)] = 0
    for (i in 1:length(pp_csM_node[,1])){
      pp_csM_node[ ,paste0('role_size_',m)][i] = sum(pp_csM_node$name[i] == get(paste0("pp_rownames_csM_",m)))
    }
    pp_csM_node_sum = cbind(pp_csM_node_sum, pp_csM_node[ ,paste0('role_size_',m)])
  }
  pp_csM_node$role_size_all = rowSums(pp_csM_node_sum)
  colnames(pp_csM_node) <- c('name', 'id', names(mpp.pp$PP), 'role_size_total')
  pp_csM_node$group_color = 'trait:--'
  for (i in 1:length(pp_csM_node[,1])){
    for (j in 3:(2+length(names(mpp.pp$PP)))){
      if (pp_csM_node[i,j]>0){
        tmp = colnames(pp_csM_node)[j]
        pp_csM_node$group_color[i] = paste(c(pp_csM_node$group_color[i],tmp), collapse="+")
      }
    }
  }
  ### ALL DONE 
  # create data:
  plt_links_snp_flashfm <- data.frame(
    source=pp_csM_link$source,
    target=pp_csM_link$target,
    l_color =pp_csM_link$color,
    l_size  =pp_csM_link$value*10
  )
  plt_nodes_snp_flashfm <- data.frame(
    name =pp_csM_node$name,
    n_color=pp_csM_node$group_color,
    n_size =pp_csM_node$role_size_total
  )
  
  for (i in 1:length(plt_links_snp_flashfm$source)){
    plt_links_snp_flashfm$source_id[i] = pp_csM_node$id[which(plt_links_snp_flashfm$source[i]==pp_csM_node$name)]
    plt_links_snp_flashfm$target_id[i] = pp_csM_node$id[which(plt_links_snp_flashfm$target[i]==pp_csM_node$name)]
  }
  ##IMPORTANT
  return(list(plt_links_snp = plt_links_snp,
              plt_nodes_snp = plt_nodes_snp,
              plt_links_snp_flashfm = plt_links_snp_flashfm,
              plt_nodes_snp_flashfm = plt_nodes_snp_flashfm))
}
## testing the whole function
# pp_network_links_nodes = pp_network_links_nodes(mpp.pp, cs1_selected, csM_selected)
# plt_links_snp = pp_network_links_nodes$plt_links_snp
# plt_nodes_snp = pp_network_links_nodes$plt_nodes_snp
# plt_links_snp_flashfm = pp_network_links_nodes$plt_links_snp_flashfm
# plt_nodes_snp_flashfm = pp_network_links_nodes$plt_nodes_snp_flashfm
#color
ColourScale_snp <- 'd3.scaleOrdinal().range(["lightblue", "lightcoral", "lightgoldenrod", "lightgreen",
                                             "lightpink", "lightsalmon", "lightyellow", "lightslateblue"]);'
### ALL DONE: snp individual network



#generate a vector of color based on snpGroups----
# set.seed(1)
# n <- length(colnames(snpGroups$group.sizes))
# col_vector = distinctColorPalette(n)
# pal <- c("gray", col_vector)
# pal <- setNames(pal, c("0", colnames(snpGroups$group.sizes)))


#generating the full table GWAS_all_final----
GWAS_all_final_table <- function(GWAS, mpp.pp, cs1, csM, snpGroups){
  M <- length(names(mpp.pp$MPP))
  for (m in 1:M){
    GWAS_0 = GWAS[[m]]
    #IMPORTANT: the order must be [1]"rs" [2]"chr" [3]"ps" [4]"allele1" [5]"allele0" [6] "maf/af" [7] "no." 
    GWAS_0 = GWAS_0[ , c(1,2,3,4,5,9,10,6:8)]
    #then: [6] "beta" [7]"se" [8]"pval"
    colnames(GWAS_0)[8:10] = c(paste0('beta_',m), paste0('se_',m), paste0('pval_',m))
    colnames(GWAS_0)[6:7] = c('af','no')
    GWAS_0[,7] = 0 #updated20220208
    if (m == 1){
      GWAS_all_2 = GWAS_0
    }else {
      GWAS_all_2 = merge(GWAS_all_2, GWAS_0, by = colnames(GWAS_all_2[,1:7]))
    }
    # Credible set
    cs1_0 = as.data.frame(cs1[[m]])
    csM_0 = as.data.frame(csM[[m]])
    GWAS_all_2[paste0("cs1_",m)] = 0
    GWAS_all_2[paste0("csM_",m)] = 0
    for (i in 1:length(GWAS_all_2[,1])){
      if (sum(GWAS_all_2$rs[i] == cs1_0)>0){
        GWAS_all_2[ ,paste0("cs1_",m)][i] = 1
      }
      if (sum(GWAS_all_2$rs[i] == csM_0)>0){
        GWAS_all_2[ ,paste0("csM_",m)][i] = 1
      }  
    }
    # MPP
    temp = as.data.frame(mpp.pp$MPP[[m]])
    temp$rs = rownames(temp)
    colnames(temp)[1:2] = c(paste0("MPP_",names(mpp.pp$MPP)[m],"_Single"),
                            paste0("MPP_",names(mpp.pp$MPP)[m],"_Multi"))
    GWAS_all_2 = merge(GWAS_all_2, temp, by = c("rs"), all.x=TRUE)
  }
  # snpGroups
  GWAS_all_final_2 = GWAS_all_2
  GWAS_all_final_2$snpGroups_fm = 0
  GWAS_all_final_2$snpGroups_flashfm = 0
  GWAS_all_final_2$snpGroups_fm_size = 0
  GWAS_all_final_2$snpGroups_flashfm_size = 0
  for (j in 1:length(snpGroups$groups.fm)){
    for (h in 1:length(GWAS_all_2[,1])){
      if (sum(GWAS_all_final_2$rs[h] == snpGroups$groups.fm[[j]])>0){
        GWAS_all_final_2$snpGroups_fm[h] = names(snpGroups$groups.fm)[j]
        GWAS_all_final_2$snpGroups_fm_size[h] = length(snpGroups$groups.fm[[j]]) ##Check!!!
      }
    }
  }
  for (j in 1:length(snpGroups$groups.flashfm)){
    for (h in 1:length(GWAS_all_2[,1])){
      if (sum(GWAS_all_final_2$rs[h] == snpGroups$groups.flashfm[[j]])>0){
        GWAS_all_final_2$snpGroups_flashfm[h] = names(snpGroups$groups.flashfm)[j]
        GWAS_all_final_2$snpGroups_flashfm_size[h] = length(snpGroups$groups.flashfm[[j]]) ##Check!!!
      }
    }
  }
  
  for (m in 1:M){  #check!!!_updated20220218
    # MPPg_fm
    temp_mppg_fm = as.data.frame(mpp.pp$MPPg[[m]][,1])
    temp_mppg_fm$snpGroups_fm = rownames(temp_mppg_fm)
    colnames(temp_mppg_fm)[1] = paste0("MPPg_",names(mpp.pp$MPPg)[m],"_Single")
    GWAS_all_final_2 = merge(GWAS_all_final_2, temp_mppg_fm, all.x=TRUE)
    
    # MPPg_flashfm
    temp_mppg_flashfm = as.data.frame(mpp.pp$MPPg[[m]][,2])
    temp_mppg_flashfm$snpGroups_flashfm = rownames(temp_mppg_flashfm)
    colnames(temp_mppg_flashfm)[1] = paste0("MPPg_",names(mpp.pp$MPPg)[m],"_Multi")
    GWAS_all_final_2 = merge(GWAS_all_final_2, temp_mppg_flashfm, all.x=TRUE)
  }
  
  # GWAS_all_final_2 %>% count(snpGroups_fm)      #check groups of fm
  # GWAS_all_final_2 %>% count(snpGroups_flashfm) #check groups of flashfm
  # GWAS_all_final_2$rs[which(GWAS_all_final$snpGroups_fm=='B')] #check individual snps in the group
  # GWAS_all_final_2$rs[which(GWAS_all_final$snpGroups_flashfm=='E')]
  #update p-values
  for (m in 1:M){
    GWAS_all_final_2[paste0("pval_log_",m)] = -log10(GWAS_all_final_2[,paste0('pval_',m)])
    GWAS_all_final_2[paste0("cs_trait_",m)] = GWAS_all_final_2[,paste0("cs1_",m)] + GWAS_all_final_2[,paste0("csM_",m)]
  }
  GWAS_all_final_2[is.na(GWAS_all_final_2)] <- 0
  GWAS_all_final_2 = GWAS_all_final_2[ , -which(names(GWAS_all_final_2) %in% c("no"))]
  return(GWAS_all_final_2)
  #End of the full table GWAS_all_final
}

#GWAS_all_final <- GWAS_all_final_table(GWAS, mpp.pp, cs1, csM)


####----Max size of input file----
options(shiny.maxRequestSize = 100*1024^2) 

#


####----Shiny user interface----
ui <- dashboardPage(
  
  #Dashboard header----
  dashboardHeader(
    title = "flashfm-ivis R package",
    # Drop-down menu for messages
    dropdownMenu(type = "messages", badgeStatus = "success",
                 headerText = 'Contact info:',
                 messageItem("Contact 1: Feng",
                             "feng.zhou[AT]mrc-bsu.cam.ac.uk"
                 ),
                 messageItem("Contact 2: Jenn",
                             "jennifer.asimit[AT]mrc-bsu.cam.ac.uk"
                 )
    ),
    
    # Drop-down menu for notifications
    dropdownMenu(type = "notifications", #badgeStatus = "warning",
                 headerText = 'Update and version info:',
                 notificationItem(icon = icon("users"), status = "info",
                                  "Update: this is Version 0.3"
                 ),
                 notificationItem(icon = icon("users"), status = "info",
                                  "Update: modified on 08-Mar-2022"
                 )
    ),
    
    # Drop-down menu for tasks, with progress bar
    dropdownMenu(type = "tasks", badgeStatus = "danger",
                 headerText = 'Source code:',
                 taskItem(value = 100, color = "aqua",
                          text = "Source code #1: flashfm", 
                          href = 'https://github.com/jennasimit/flashfm'
                 ),
                 taskItem(value = 100, color = "green",
                          text = "Source code #2: flashfm-ivis", 
                          href = 'https://github.com/fz-cambridge/flashfm-ivis'
                 ),
                 taskItem(value = 100, color = "yellow",
                          text = "Source code #3: finemap-ivis", 
                          href = 'http://shiny.mrc-bsu.cam.ac.uk/apps/finemap-ivis/'
                 )
    )
    ),
  
  #Dashboard sidebar----
  dashboardSidebar(
    #Data inputs--
    radioButtons(inputId = "control_widgets_0", 
                 label = h4("Select data source:"),
                 choices = list("Pre-loaded example" = 0, 
                                "FINEMAP with flashfm" = 1,
                                "FLASHFM with JAM" = 2,
                                "Single RData file" = 3,
                                "Single csv file" = 4), 
                 selected = 0),
    
    fileInput(inputId = "inFiles",
              label = h4("Upload external data"),
              multiple = TRUE, accept = NULL, buttonLabel = "Browse..."),
    
    #Control widgets--
    radioButtons(inputId = "control_widgets_1", 
                 label = h4("Credible sets:"),
                 choices = list("Show all SNPs (default)" = 3, 
                                "Show cs1 and csM" = 2,
                                "Show cs1 or csM" = 1, 
                                "Show cs1 but not csM" = 1.1, 
                                "Show csM but not cs1" = 1.2, 
                                "Show no cs" = 0 ), 
                 selected = 3),
    
    sliderInput(inputId = "control_widgets_2", 
                label = h4("Define MPP value"), 
                min = 0, max = 1, value = c(0, 1)),   
    
    #Dashboard Menu--
    sidebarMenu(
        menuItem("Dashboard", icon = icon("dashboard"), startExpanded = TRUE,
              menuSubItem("Checking input values", tabName = "input_checking", 
                          icon = icon("angle-right")), 
              menuSubItem("Ch_0 Input values", tabName = "ch_0", 
                          icon = icon("angle-right")),   
              menuSubItem("Ch_1 Single trait", tabName = "ch_1", 
                          icon = icon("angle-right")),
              menuSubItem("Ch_2 Multi trait", tabName = "ch_2", 
                          icon = icon("angle-right")),
              menuSubItem("Ch_3 Comparison (Linked)", tabName = "ch_3", 
                          icon = icon("angle-right")),
              menuSubItem("Ch_4 Comparison (Methods)", tabName = "ch_4", 
                          icon = icon("angle-right")),
              menuSubItem("Ch_5 Comparison (Traits)", tabName = "ch_5", 
                          icon = icon("angle-right")),              
              menuSubItem("Data and table details", tabName = "data_detail", 
                          icon = icon("angle-right")),
              menuSubItem("README", tabName = "readme", 
                          icon = icon("angle-right")),
              menuSubItem(text = "Go to finemap-ivis", 
                          href = "http://shiny.mrc-bsu.cam.ac.uk/apps/finemap-ivis/", 
                          icon = icon("angle-right"))
        )
    )
  ),
  
  #Dashboard body----
  dashboardBody(  
    tabItems(
      #Input_checking----
      tabItem(tabName = "input_checking",
              #check all input files
              h3("Checking input files"),
              #h3("Check files in the current temporary server"),
              fluidRow(column(12, verbatimTextOutput("value"))),
              br(),
              h4("--- Key files' location in the server ---"),
              fluidRow(column(12, DT::dataTableOutput("file_locations"))),
              br(),
              h4("--- Information about uploading users' external data ---"),
              br(),
              h5("Case 1. Preloaded examples: see README"),
              br(),
              h5("Case 2. FINEMAP with flashfm: see README"),
              br(),
              h5("Case 3. FLASHFM with JAM: see README"),
              br(),
              h5("Case 4. Single RData file: see README"),
              br(),
              h5("Case 5. Single csv file: see README"),
      ),
      #Dashboard_ch_0 GWAS Summary Statistics and Input Values----
      tabItem(tabName = "ch_0",
              h2("Dashboard Charts 0: GWAS Summary Statistics and Input Values"),
                         fluidPage(
                           fluidRow(
                             column(# width should be between 1 and 12
                               width=12,
                               tabBox(
                                 title = "Single-trait fine-mapping (FINEMAP)",
                                 # The id lets us use input$tabset1 on the server to find the current tab
                                 id = "GWAS_tabset_1", height = "550px",
                                 tabPanel("Trait 1", 
                                          plotlyOutput("manhattan_1a"),
                                          radioButtons(inputId = "control_widgets_ch_0_a1", 
                                                       label = h4("Show a dashed line at the maximum value:"),
                                                       inline=T,
                                                       choices = list("Yes" = 0, 
                                                                      "No" = 1), 
                                                       selected = 0)),
                                 tabPanel("Trait 2", plotlyOutput("manhattan_2a"),
                                          radioButtons(inputId = "control_widgets_ch_0_a2", 
                                                       label = h4("Show a dashed line at the maximum value:"),
                                                       inline=T,
                                                       choices = list("Yes" = 0, 
                                                                      "No" = 1), 
                                                       selected = 0)),
                                 #tabPanel("Trait 3", plotlyOutput("manhattan_3a")),
                                 tabPanel("Trait 3", 
                                          selectInput("shouldShow_3a", "show plot (if you have 3 traits), otherwise showing error", c("yes", "no"), "no"),
                                          conditionalPanel(
                                            condition = "input.shouldShow_3a == 'yes'",
                                            plotlyOutput("manhattan_3a")
                                          )),
                                 #tabPanel("Trait 4", plotlyOutput("manhattan_4a")),
                                 tabPanel("Trait 4", 
                                          selectInput("shouldShow_4a", "show plot (if you have 4 traits), otherwise showing error", c("yes", "no"), "no"),
                                          conditionalPanel(
                                            condition = "input.shouldShow_4a == 'yes'",
                                            plotlyOutput("manhattan_4a")
                                          )),
                                 tabPanel("Trait 5", 
                                          selectInput("shouldShow_5a", "show plot (if you have 5 traits), otherwise showing error", c("yes", "no"), "no"),
                                          conditionalPanel(
                                            condition = "input.shouldShow_5a == 'yes'",
                                            plotlyOutput("manhattan_5a")
                                          )),
                                 tabPanel("Trait 6", 
                                          selectInput("shouldShow_6a", "show plot (if you have 6 traits), otherwise showing error", c("yes", "no"), "no"),
                                          conditionalPanel(
                                            condition = "input.shouldShow_6a == 'yes'",
                                            plotlyOutput("manhattan_6a")
                                          )),         
                                 width=NULL
                               ),
                               tabBox(
                                 title = "Multi-trait fine-mapping (flashfm)",
                                 # The id lets us use input$tabset1 on the server to find the current tab
                                 id = "GWAS_tabset_2", height = "550px",
                                 tabPanel("Trait 1", 
                                          plotlyOutput("manhattan_1b"),
                                          radioButtons(inputId = "control_widgets_ch_0_b1", 
                                                       label = h4("Show a dashed line at the maximum value:"),
                                                       inline=T,
                                                       choices = list("Yes" = 0, 
                                                                      "No" = 1), 
                                                       selected = 0)),
                                 tabPanel("Trait 2", 
                                          plotlyOutput("manhattan_2b"),
                                          radioButtons(inputId = "control_widgets_ch_0_b2", 
                                                       label = h4("Show a dashed line at the maximum value:"),
                                                       inline=T,
                                                       choices = list("Yes" = 0, 
                                                                      "No" = 1), 
                                                       selected = 0)),
                                 #tabPanel("Trait 3", plotlyOutput("manhattan_3b")),
                                 tabPanel("Trait 3", 
                                          selectInput("shouldShow_3b", "show plot (if you have 3 traits), otherwise showing error", c("yes", "no"), "no"),
                                          conditionalPanel(
                                            condition = "input.shouldShow_3b == 'yes'",
                                            plotlyOutput("manhattan_3b")
                                          )),
                                 #tabPanel("Trait 4", plotlyOutput("manhattan_4b")),
                                 tabPanel("Trait 4", 
                                          selectInput("shouldShow_4b", "show plot (if you have 4 traits), otherwise showing error", c("yes", "no"), "no"),
                                          conditionalPanel(
                                            condition = "input.shouldShow_4b == 'yes'",
                                            plotlyOutput("manhattan_4b")
                                          )),
                                 tabPanel("Trait 5", 
                                          selectInput("shouldShow_5b", "show plot (if you have 5 traits), otherwise showing error", c("yes", "no"), "no"),
                                          conditionalPanel(
                                            condition = "input.shouldShow_5b == 'yes'",
                                            plotlyOutput("manhattan_5b")
                                          )),
                                 tabPanel("Trait 6", 
                                          selectInput("shouldShow_6b", "show plot (if you have 6 traits), otherwise showing error", c("yes", "no"), "no"),
                                          conditionalPanel(
                                            condition = "input.shouldShow_6b == 'yes'",
                                            plotlyOutput("manhattan_6b")
                                          )),  
                                 width=NULL
                               ),
                             ),  
                             column(
                               width=12,
                               box(div(style = 'overflow-y: scroll; overflow-x: scroll; ', DT::dataTableOutput('ch0_table0')),
                                   title = "snpGroup$group.sizes",
                                   solidHeader = TRUE,
                                   collapsible = TRUE, width=NULL)
                             ),   
                             column(
                               width=12,
                               box(plotlyOutput("Plot_LD1"),
                                   title = "LD Matrix", 
                                   solidHeader = TRUE, collapsible = TRUE,
                                   width=NULL),
                               box(plotOutput("Plot_LD2"),
                                   downloadLink("downloadPlot", "Download Plot"),
                                   title = "LD Plot", 
                                   solidHeader = TRUE, collapsible = TRUE,
                                   width=NULL)
                             )
                           )
                         )
      ),
      #Dashboard_ch_1 forceNetwork (Single)----
      tabItem(tabName = "ch_1",
              h2("Dashboard Chart 1: Single Trait Fine Mapping"),
              fluidRow(
                column(# width should be between 1 and 12
                  width=12,
                  tabBox(
                    id = "tabset1", height = "1000px",
                    tabPanel("Network: Group", 
                             sliderInput(inputId = "control_widgets_3a_fm", 
                                         label = h4("Define thresholds for group model PP (PPg)"), 
                                         #min = 0, max = max(input$tt()$plt_links_group_fm$l_size)/10, value = c(0, max(tt()$plt_links_group_fm$l_size)/10)), 
                                         min = 0, max = 1, value = c(0, 1)), 
                             forceNetworkOutput(outputId = "ch1_network_group_fm", height = "800px"),
                             downloadButton('d_ch_1_s1', 'Download network as html')
                             ),
                    tabPanel("Network: SNPs",     
                             sliderInput(inputId = "control_widgets_3b_fm", 
                                         label = h4("Define thresholds for SNP model PP (PP)"), 
                                         #min = 0, max = max(tt()$plt_links_snp$l_size)/10, value = c(0.1, max(tt()$plt_links_snp$l_size)/10)),   
                                         min = 0, max = 1, value = c(0.1, 1)),   
                             forceNetworkOutput(outputId = "ch1_network_snp", height = "800px"),
                             downloadButton('d_ch_1_s2', 'Download network as html')
                             ),
                    width=NULL
                  )
                )
              )
      ),
      #Dashboard_ch_2 forceNetwork (Multi)----
      tabItem(tabName = "ch_2",
              h2("Dashboard Chart 2: Multi Trait Fine Mapping"),
              fluidRow(
                column(# width should be between 1 and 12
                  width=12,
                  tabBox(
                    id = "tabset1", height = "1000px",
                    tabPanel("Network: Group",
                             sliderInput(inputId = "control_widgets_3a_flashfm",
                                         label = h4("Define thresholds for group model PP (PPg)"),
                                         #min = 0, max = max(plt_links_group_flashfm$l_size)/10, value = c(0, max(plt_links_group_flashfm$l_size)/10)),
                                         min = 0, max = 1, value = c(0, 1)),
                             forceNetworkOutput(outputId = "ch1_network_group_flashfm", height = "800px"),
                             downloadButton('d_ch_2_m1', 'Download network as html')
                             ),
                    tabPanel("Network: SNPs",
                             sliderInput(inputId = "control_widgets_3b_flashfm",
                                         label = h4("Define thresholds for SNP model PP (PP)"),
                                         #min = 0, max = max(plt_links_snp_flashfm$l_size)/10, value = c(0.1, max(plt_links_snp_flashfm$l_size)/10)),
                                         min = 0, max = 1, value = c(0.1, 1)),
                             forceNetworkOutput(outputId = "ch1_network_snp_flashfm", height = "800px"),
                             downloadButton('d_ch_2_m2', 'Download network as html')
                             ),
                    width=NULL
                  )
                )
              )
       ),
      #Dashboard_ch_3 Linked Regional association plots----
      tabItem(tabName = "ch_3",
              h2("Dashboard Chart 3: Comparison Regional Association Plots of Different Traits and Methods"),
              fluidPage(
                fluidRow(
                  column(
                    width=12,
                    box(div(style = 'overflow-y: scroll; overflow-x: scroll; ', DT::dataTableOutput('ch0_table0_ch3')),
                        title = "snpGroup$group.sizes",
                        solidHeader = TRUE,
                        collapsible = TRUE, width=NULL)
                  ),
                  column(
                    width=12,
                    radioButtons(inputId = "control_widgets_ch_3", 
                                 inline = T,
                                 label = h4("How many traits: (NOTE - left panel is fm, right panel is flashfm, rows show individual traits)"),
                                 choices = list("2 traits" = 2, 
                                                "3 traits" = 3,
                                                "4 traits" = 4, 
                                                "5 traits" = 5, 
                                                "6 traits" = 6), 
                                 selected = 2),
                    radioButtons(inputId = "control_widgets_ch_3_2", 
                                 label = h4("Show a dashed line at the maximum value:"),
                                 inline=T,
                                 choices = list("Yes" = 0, 
                                                "No" = 1), 
                                 selected = 0)
                  ),
                  column(
                    width=12,
                    box(plotlyOutput("manhattan_ch3_linked"), 
                        height = 1400, width=NULL),
                  )
                ) 
              )
      ),
      #Dashboard_ch_4 Radar chart, Venn diagrams----
      tabItem(tabName = "ch_4",
              h2("Dashboard Chart 4: Comparison Results of Different Methods"),
              fluidRow(
                column(
                  width=12,
                  radioButtons(inputId = "control_widgets_ch_4", 
                               inline = T,
                               label = h4("How many traits:"),
                               choices = list("3 traits" = 3,
                                              "4 traits" = 4, 
                                              "5 traits" = 5, 
                                              "6 traits" = 6), 
                               selected = 3),
                ),
                column(# width should be between 1 and 12
                  width=12,
                  tabBox(
                    #title = "Result_1 charts",
                    # The id lets us use input$tabset1 on the server to find the current tab
                    id = "tabset1", height = "800px",
                    tabPanel("Radar chart: Credible Sets", plotlyOutput("fig_ch4_tab1_radar")),
                    tabPanel("Venn diagram 1: cs1", 
                             plotlyOutput("fig_ch4_tab2_venn_cs1"),
                             #tableOutput("table_ch4_tab2_venn_cs1"),
                             box(div(style = 'overflow-y: scroll; overflow-x: scroll; ', 
                                     DT::dataTableOutput('table_ch4_tab2_venn_cs1_table')),
                                     title = "Table of credible set cs_1",
                                     #status = "warning", 
                                     solidHeader = TRUE,
                                     collapsible = TRUE, width=NULL)
                             ),
                    tabPanel("Venn diagram 2: csM", 
                             plotlyOutput("fig_ch4_tab2_venn_csM"),
                             #tableOutput("table_ch4_tab2_venn_cs1"),
                             box(div(style = 'overflow-y: scroll; overflow-x: scroll; ', 
                                     DT::dataTableOutput('table_ch4_tab2_venn_csM_table')),
                                 title = "Table of credible set cs_M",
                                 #status = "warning", 
                                 solidHeader = TRUE,
                                 collapsible = TRUE, width=NULL)
                             ),
                    width=NULL
                  )
                )
              )
      ),
      #Dashboard_ch_5 Sankey Diagrams of Different Traits----
      tabItem(tabName = "ch_5",
              h2("Dashboard Chart 5: Sankey Diagrams of Different Traits"),
              fluidRow(
                column(# width should be between 1 and 12
                  width=12,
                  tabBox(
                    #title = "Result_1 charts",
                    # The id lets us use input$tabset1 on the server to find the current tab
                    id = "tabset1", height = "2200px",
                    tabPanel("All traits", plotlyOutput("ch5_sankey_all", height = "2000px")),
                    tabPanel("Trait_1", plotlyOutput("ch5_sankey_t1", height = "1400px")),
                    tabPanel("Trait_2", plotlyOutput("ch5_sankey_t2", height = "1400px")),
                    #tabPanel("Trait_3", plotlyOutput("ch5_sankey_t3", height = "800px")),
                    tabPanel("Trait_3", 
                             selectInput("shouldShow_ch5a", "show plot (if you have 3 traits), otherwise showing error", c("yes", "no"), "no"),
                             conditionalPanel(
                               condition = "input.shouldShow_ch5a == 'yes'",
                               plotlyOutput("ch5_sankey_t3", height = "1500px")
                             )),
                    #tabPanel("Trait_4", plotlyOutput("ch5_sankey_t4", height = "800px")),  
                    tabPanel("Trait_4", 
                             selectInput("shouldShow_ch5b", "show plot (if you have 4 traits), otherwise showing error", c("yes", "no"), "no"),
                             conditionalPanel(
                               condition = "input.shouldShow_ch5b == 'yes'",
                               plotlyOutput("ch5_sankey_t4", height = "1500px")
                             )),
                    #tabPanel("Trait_5", plotlyOutput("ch5_sankey_t5", height = "800px")),   
                    tabPanel("Trait_5", 
                             selectInput("shouldShow_ch5c", "show plot (if you have 5 traits), otherwise showing error", c("yes", "no"), "no"),
                             conditionalPanel(
                               condition = "input.shouldShow_ch5c == 'yes'",
                               plotlyOutput("ch5_sankey_t5", height = "1500px")
                             )),
                    #tabPanel("Trait_6", plotlyOutput("ch5_sankey_t6", height = "800px")),   
                    tabPanel("Trait_6", 
                             selectInput("shouldShow_ch5d", "show plot (if you have 6 traits), otherwise showing error", c("yes", "no"), "no"),
                             conditionalPanel(
                               condition = "input.shouldShow_ch5d == 'yes'",
                               plotlyOutput("ch5_sankey_t6", height = "1500px")
                             )),
                    width=NULL
                  )
                )
              )
      ),
      #Dashboard_ch_6----
      tabItem(tabName = "data_detail",
              h2("Data details"),
              # h4("---file name and location---"),
              # tableOutput("web_inputs_1"),              
              h4("---cs1---"),
              verbatimTextOutput("web_inputs_cs1"),
              h4("---csM---"),
              verbatimTextOutput("web_inputs_csM"),
              h4("---GWAS---"),
              verbatimTextOutput("web_inputs_GWAS"),
              h4("---LD---"),
              verbatimTextOutput("web_inputs_LD"),
              h4("---mpp.pp---"),
              verbatimTextOutput("web_inputs_mpp.pp"),
              h4("---snpGroups---"),
              verbatimTextOutput("web_inputs_snpGroups"),
              h4("---GWAS_all_final---"),
              div(style = 'overflow-x: scroll', DT::dataTableOutput('GWAS_all_final_table_done')),
              downloadButton("downloadData_3", "Download"),
              h4("---LD---"),
              div(style = 'height:500px; overflow-y: scroll; overflow-x: scroll; ', DT::dataTableOutput('LD_data'))
      ),
      #Dashboard_ch_7----
      tabItem(tabName = "readme",
              h2("README"),
              br(),
              h3("About this package:"),
              br(),
              h4("GitHub: "),
              br(),
              h4("flashfm R package: https://github.com/jennasimit/flashfm"),
              h4("flashfm-ivis R package: https://github.com/fz-cambridge/flashfm-ivis"),
              br(),
              h4("YouTube: "),
              #br(),
              #HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/videoseries?list=PLcy5X5WM9r0AqRrEb5vUTfRkcBRUoOc8T" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')
              fluidRow(
                fluidRow(
                column(6, 
                       box(
                         width = NULL, 
                         title = "flashfm-ivis and finemap-ivis overview", 
                         HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/hDbV9qNzsZo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')
                       )),
                column(6, 
                       box(
                         width = NULL, 
                         title = "finemap-ivis user data input", 
                         HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/dTRdfsMumyY" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')
                       )),
                column(6, 
                       box(
                         width = NULL, 
                         title = "flashfm-ivis user data input: FINEMAP with Flashfm", 
                         HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/dxcWc2fuIHI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')
                       )),
                column(6, 
                       box(
                         width = NULL, 
                         title = "flashfm-ivis user data input: FLASHFM with JAM", 
                         HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/IYmJJCqfz0A" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')
                       )),
                column(6, 
                       box(
                         width = NULL, 
                         title = "flashfm-ivis user data input: Single RData file", 
                         HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/tQaEhFfcbKE" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')
                       )),
                )
              )
              
      )
    )
  )
) #End of user interface
####----End of user interface----





####----Shiny server----
server <- function(input, output) {
  
  #dashboard 0 input file information -done-
  output$value <- renderPrint({
    if (input$control_widgets_0==0){
      print("You are using pre-loaded examples")
    }else{
      str(input$inFiles)}
      })
  #dashboard 0 input file location -done-
  output$file_locations <- DT::renderDataTable({
    if (input$control_widgets_0==0){
      f_config <- f_ld <- f_snp <- f_z <- f_rdata <- f_csv <- "You are using pre-loaded examples, therefore no external file here"
    }else if(input$control_widgets_0==1){
      file <- input$inFiles
      req(file)
      f_config <- f_ld <- f_snp <- f_z <- f_rdata <- f_csv <- 'You selected data source "Finemap with flashfm", 
                                                               it requires all files with the following extensions: .rdata; .z; .ld.
                                                               (e.g. mpp.pp and snpGroups in .RData, LD matrix in .ld, finemap.z files)'
      for (i in 1:length(file$datapath)){
        ext_check <- tools::file_ext(file$datapath[i])
        if (ext_check == "config"){f_config = file$datapath[i]} 
        if (ext_check == "ld"){f_ld = file$datapath[i]} 
        if (ext_check == "snp"){f_snp = file$datapath[i]} 
        if (ext_check == "z"){f_z = file$datapath[i]} 
        if (ext_check == "RData"){f_rdata = file$datapath[i]} 
        if (ext_check == "csv"){f_csv = file$datapath[i]} 
      }
    }else if(input$control_widgets_0==2){
      file <- input$inFiles
      req(file)
      f_config <- f_ld <- f_snp <- f_z <- f_rdata <- f_csv <- 'You selected data source "FLASHFMwithJAM", 
                                                               it requires all files with the following extensions: .rdata files only. 
                                                               (e.g. flashfmwithjam.RData and snpinfo.RData)'
      for (i in 1:length(file$datapath)){
        ext_check <- tools::file_ext(file$datapath[i])
        if (ext_check == "config"){f_config = file$datapath[i]} 
        if (ext_check == "ld"){f_ld = file$datapath[i]} 
        if (ext_check == "snp"){f_snp = file$datapath[i]} 
        if (ext_check == "z"){f_z = file$datapath[i]} 
        if (ext_check == "RData"){f_rdata = file$datapath[i]} 
        if (ext_check == "csv"){f_csv = file$datapath[i]} 
      }
    }else if(input$control_widgets_0==3){
      file <- input$inFiles
      req(file)
      f_config <- f_ld <- f_snp <- f_z <- f_rdata <- f_csv <- 'You selected data source "Single RData file", 
                                                               it requires ONE .rdata file only that contains all data.'
      for (i in 1:length(file$datapath)){
        ext_check <- tools::file_ext(file$datapath[i])
        if (ext_check == "config"){f_config = file$datapath[i]} 
        if (ext_check == "ld"){f_ld = file$datapath[i]} 
        if (ext_check == "snp"){f_snp = file$datapath[i]} 
        if (ext_check == "z"){f_z = file$datapath[i]} 
        if (ext_check == "RData"){f_rdata = file$datapath[i]} 
        if (ext_check == "csv"){f_csv = file$datapath[i]} 
      }
    }
    location = rbind(f_rdata, f_csv, f_config, f_ld, f_snp, f_z)
    colnames(location) = c('location')
    DT::datatable(location, options = list(dom = 't'))
  })
  
  
  tt <- reactive({
    withProgress(message = 'Calculation in progress',
                 detail = '...Step 1/10 done!', value = 0, {
                   for (i in 1:10) {
                     incProgress(1/10)
                     Sys.sleep(0.025)
                   }
                 })
    if (input$control_widgets_0 == 0){
      load("data/R_data_preloaded.RData")
      withProgress(message = 'Calculation in progress',
                   detail = '...Step 2/10 done!', value = 0, {
                     for (i in 1:10) {
                       incProgress(1/10)
                       Sys.sleep(0.025)
                     }
                   })
      # #snp group network: ppg_network_links and nodes function---
      # ppg_network_links_nodes = ppg_network_links_nodes(mpp.pp)
      # plt_links_group_fm = ppg_network_links_nodes$plt_links_group_fm
      # plt_nodes_group_fm = ppg_network_links_nodes$plt_nodes_group_fm
      # plt_links_group_flashfm = ppg_network_links_nodes$plt_links_group_flashfm
      # plt_nodes_group_flashfm = ppg_network_links_nodes$plt_nodes_group_flashfm
      # # ColourScale_group <- 'd3.scaleOrdinal().range(["lightblue", "lightcoral", "lightgoldenrod", "lightgreen",
      # #                                                "lightpink", "lightsalmon", "lightyellow", "lightslateblue"]);'
      # #snp individual network: pp_network_links and nodes function---
      # M <- length(mpp.pp$PP) # number of traits 
      # cs1_selected <- csM_selected <- vector("list",M)
      # for(i in 1:M){
      #   tmp <- mpp.pp$PP[[i]][,1] # trait i, single-trait fine-mapping
      #   cs1_selected[[i]] <- credset(tmp,0.99)
      #   tmp <- mpp.pp$PP[[i]][,2] # trait i, multi-trait fine-mapping
      #   csM_selected[[i]] <- credset(tmp,0.99)
      # }
      # pp_network_links_nodes = pp_network_links_nodes(mpp.pp, cs1_selected, csM_selected)
      # plt_links_snp = pp_network_links_nodes$plt_links_snp
      # plt_nodes_snp = pp_network_links_nodes$plt_nodes_snp
      # plt_links_snp_flashfm = pp_network_links_nodes$plt_links_snp_flashfm
      # plt_nodes_snp_flashfm = pp_network_links_nodes$plt_nodes_snp_flashfm
      # #generate a vector of color based on snpGroups---
      # set.seed(1)
      # n <- length(colnames(snpGroups$group.sizes))
      # col_vector = distinctColorPalette(n)
      # pal <- c("gray", col_vector)
      # pal <- setNames(pal, c("0", colnames(snpGroups$group.sizes)))
      # #generating the full table GWAS_all_final---
      # GWAS_all_final <- GWAS_all_final_table(GWAS, mpp.pp, cs1, csM, snpGroups)
      # #outputs
      # tt = list(plt_links_group_fm = plt_links_group_fm,
      #           plt_nodes_group_fm = plt_nodes_group_fm,
      #           plt_links_group_flashfm = plt_links_group_flashfm,
      #           plt_nodes_group_flashfm = plt_nodes_group_flashfm,
      #           plt_links_snp = plt_links_snp,
      #           plt_nodes_snp = plt_nodes_snp,
      #           plt_links_snp_flashfm = plt_links_snp_flashfm,
      #           plt_nodes_snp_flashfm = plt_nodes_snp_flashfm,
      #           col_vector = col_vector,
      #           pal = pal,
      #           LD = LD,
      #           cs1 = cs1,
      #           csM = csM,
      #           GWAS = GWAS,
      #           mpp.pp = mpp.pp,
      #           snpGroups = snpGroups,
      #           cs1_selected = cs1_selected,
      #           csM_selected = csM_selected,
      #           GWAS_all_final = GWAS_all_final)
    }else if (input$control_widgets_0 == 1){                       ##CHECK!!!--DONE
      #1. Finemapwithflashfm: four .z + one .ld files and mpp.pp+snpGroups from the .Rdata file
      file <- input$inFiles
      req(file)
      f_z <- numeric(1);
      j <- 1;
      for (i in 1:length(file$datapath)){
        ext_check <- tools::file_ext(file$datapath[i])
        if (ext_check == "ld"){f_ld_1 = file$datapath[i]}
        if (ext_check == "z"){
          f_z[j] = file$datapath[i]
          j=j+1
          }
        if (ext_check == "RData" || ext_check == "rdata" ){f_rdata_1 = file$datapath[i]}
      }
      load(f_rdata_1) #loading the RData file (i.e. mpp.pp and snpGroups)
      wrapped_1 = format.flashfm.finemap(fz=f_z, fld=f_ld_1, mpp.pp=mpp.pp, snpGroups = snpGroups)
      GWAS = wrapped_1$GWAS
      LD = as.matrix(wrapped_1$LD)
      withProgress(message = 'Calculation in progress',
                   detail = '...Step 2/10 done!', value = 0, {
                     for (i in 1:10) {
                       incProgress(1/10)
                       Sys.sleep(0.025)
                     }
                   })
    }else if (input$control_widgets_0 == 2){                       ##CHECK!!!ONLY NEED UPDATE HERE
      #2. FLASHFMwithJAM: four .z + one .ld files and mpp.pp+snpGroups from the two different .Rdata files
      file <- input$inFiles
      req(file)
      f_z <- numeric(1);
      j <- 1;
      for (i in 1:length(file$datapath)){
        ext_check <- tools::file_ext(file$datapath[i])
        if (ext_check == "ld"){f_ld_1 = file$datapath[i]}
        if (ext_check == "z"){
          f_z[j] = file$datapath[i]
          j=j+1
        }
        if (ext_check == "RData" || ext_check == "rdata" ){f_rdata_1 = file$datapath[i]}
      }
      load(f_rdata_1) #loading the RData file (i.e. mpp.pp and snpGroups)
      mpp.pp = fm$mpp.pp
      snpGroups = fm$snpGroups
      wrapped_1 = format.flashfm.finemap(fz=f_z, fld=f_ld_1, mpp.pp=mpp.pp, snpGroups = snpGroups)
      GWAS = wrapped_1$GWAS
      LD = as.matrix(wrapped_1$LD)
      withProgress(message = 'Calculation in progress',
                   detail = '...Step 2/10 done!', value = 0, {
                     for (i in 1:10) {
                       incProgress(1/10)
                       Sys.sleep(0.025)
                     }
                   })
    }else if (input$control_widgets_0 == 3){      ##CHECK!!!
      file <- input$inFiles
      req(file)
      load(file$datapath) 
      withProgress(message = 'Calculation in progress',
                   detail = '...Step 2/10 done!', value = 0, {
                     for (i in 1:10) {
                       incProgress(1/10)
                       Sys.sleep(0.025)
                     }
                   })
    }
    #snp group network: ppg_network_links and nodes function---
    ppg_network_links_nodes = ppg_network_links_nodes(mpp.pp)
    withProgress(message = 'Calculation in progress',
                 detail = '...Step 3/10 done!', value = 0, {
                   for (i in 1:10) {
                     incProgress(1/10)
                     Sys.sleep(0.025)
                   }
                 })
    plt_links_group_fm = ppg_network_links_nodes$plt_links_group_fm
    plt_nodes_group_fm = ppg_network_links_nodes$plt_nodes_group_fm
    plt_links_group_flashfm = ppg_network_links_nodes$plt_links_group_flashfm
    plt_nodes_group_flashfm = ppg_network_links_nodes$plt_nodes_group_flashfm
    # ColourScale_group <- 'd3.scaleOrdinal().range(["lightblue", "lightcoral", "lightgoldenrod", "lightgreen",
    #                                                "lightpink", "lightsalmon", "lightyellow", "lightslateblue"]);'
    #snp individual network: pp_network_links and nodes function---
    
    withProgress(message = 'Calculation in progress', value = 0, {
    M <- length(mpp.pp$PP) # number of traits 
    cs1_selected <- csM_selected <- vector("list",M)
    for(i in 1:M){
      tmp <- mpp.pp$PP[[i]][,1] # trait i, single-trait fine-mapping
      cs1_selected[[i]] <- credset(tmp,0.99)
      tmp <- mpp.pp$PP[[i]][,2] # trait i, multi-trait fine-mapping
      csM_selected[[i]] <- credset(tmp,0.99)
      # Increment the progress bar, and update the detail text.
      incProgress(1/M, detail = paste("...Step 4/10 done!", i))
      # Pause for 0.1 seconds to simulate a long computation.
      # Sys.sleep(1)
    }
    })
    cs1 = cs1_selected
    csM = csM_selected
    withProgress(message = 'Calculation in progress',
                 detail = '...Step 5/10 done!', value = 0, {
                   for (i in 1:10) {
                     incProgress(1/10)
                     Sys.sleep(0.1)
                   }
                 })
    pp_network_links_nodes = pp_network_links_nodes(mpp.pp, cs1_selected, csM_selected)
    withProgress(message = 'Calculation in progress',
                 detail = '...Step 6/10 done!', value = 0, {
                   for (i in 1:10) {
                     incProgress(1/10)
                     Sys.sleep(0.1)
                   }
                 })
    plt_links_snp = pp_network_links_nodes$plt_links_snp
    plt_nodes_snp = pp_network_links_nodes$plt_nodes_snp
    plt_links_snp_flashfm = pp_network_links_nodes$plt_links_snp_flashfm
    plt_nodes_snp_flashfm = pp_network_links_nodes$plt_nodes_snp_flashfm
    #generate a vector of color based on snpGroups---
    set.seed(1)
    n <- length(colnames(snpGroups$group.sizes))
    col_vector = distinctColorPalette(n)
    pal <- c("gray", col_vector)
    pal <- setNames(pal, c("0", colnames(snpGroups$group.sizes)))
    withProgress(message = 'Calculation in progress',
                 detail = '...Step 7/10 done!', value = 0, {
                   for (i in 1:10) {
                     incProgress(1/10)
                     Sys.sleep(0.025)
                   }
                 })
    #generating the full table GWAS_all_final---
    GWAS_all_final <- GWAS_all_final_table(GWAS, mpp.pp, cs1, csM, snpGroups)
    #withProgress_1
    withProgress(message = 'Calculation in progress',
                 detail = '...Step 8/10 done!', value = 0, {
                   for (i in 1:10) {
                     incProgress(1/10)
                     Sys.sleep(0.025)
                   }
                 })
    
    if (min(LD)<0){    #updated20220208
      LD = LD^2
    }
    #outputs
    tt = list(plt_links_group_fm = plt_links_group_fm,
              plt_nodes_group_fm = plt_nodes_group_fm,
              plt_links_group_flashfm = plt_links_group_flashfm,
              plt_nodes_group_flashfm = plt_nodes_group_flashfm,
              plt_links_snp = plt_links_snp,
              plt_nodes_snp = plt_nodes_snp,
              plt_links_snp_flashfm = plt_links_snp_flashfm,
              plt_nodes_snp_flashfm = plt_nodes_snp_flashfm,
              col_vector = col_vector,
              pal = pal,
              LD = LD,
              cs1 = cs1,
              csM = csM,
              GWAS = GWAS,
              mpp.pp = mpp.pp,
              snpGroups = snpGroups,
              cs1_selected = cs1_selected,
              csM_selected = csM_selected,
              GWAS_all_final = GWAS_all_final)
  })
  
  
  
  #####GWAS table trait 1a single----  
  GWAS_select_t1_a_cs <- reactive({    
    withProgress(message = 'Calculation in progress',
                 detail = '...Step 9/10 done!', value = 0, {
                   for (i in 1:10) {
                     incProgress(1/10)
                     Sys.sleep(0.025)
                   }
                 })
    GWAS_all_final = tt()$GWAS_all_final
    if (input$control_widgets_1 == 0){
      GWAS_all_final[GWAS_all_final$cs_trait_1==0, ]
    }else if (input$control_widgets_1 == 1){
      GWAS_all_final[GWAS_all_final$cs1_1==1, ]
    }else if (input$control_widgets_1 == 2){
      GWAS_all_final[GWAS_all_final$cs_trait_1==2, ]
    }else if (input$control_widgets_1 == 1.1){
      GWAS_all_final[GWAS_all_final$cs_trait_1==1 & GWAS_all_final$cs1_1==1, ]
    }else if (input$control_widgets_1 == 1.2){
      GWAS_all_final[GWAS_all_final$cs_trait_1==1 & GWAS_all_final$csM_1==1, ]  
    }else{
      GWAS_all_final
    }
  })
  GWAS_select_t1_a <- reactive({    
    mpp.pp = tt()$mpp.pp
    GWAS_select_t1_a_cs()[GWAS_select_t1_a_cs()[, paste0("MPP_",names(mpp.pp$MPP)[1],"_Single")]>=input$control_widgets_2[1] 
                        & GWAS_select_t1_a_cs()[, paste0("MPP_",names(mpp.pp$MPP)[1],"_Single")]<=input$control_widgets_2[2]+0.001, ]
  })
  
  #####GWAS table trait 1b multi----
  GWAS_select_t1_b_cs <- reactive({
    withProgress(message = 'Calculation in progress',
                 detail = '...Step 10/10 done!', value = 0, {
                   for (i in 1:10) {
                     incProgress(1/10)
                     Sys.sleep(0.025)
                   }
                 })
    GWAS_all_final = tt()$GWAS_all_final
    if (input$control_widgets_1 == 0){
      GWAS_all_final[GWAS_all_final$cs_trait_1==0, ]
    }else if (input$control_widgets_1 == 1){
      GWAS_all_final[GWAS_all_final$csM_1==1, ]
    }else if (input$control_widgets_1 == 2){
      GWAS_all_final[GWAS_all_final$cs_trait_1==2, ]
    }else if (input$control_widgets_1 == 1.1){
      GWAS_all_final[GWAS_all_final$cs_trait_1==1 & GWAS_all_final$cs1_1==1, ]
    }else if (input$control_widgets_1 == 1.2){
      GWAS_all_final[GWAS_all_final$cs_trait_1==1 & GWAS_all_final$csM_1==1, ]
    }else{
      GWAS_all_final
    }
  })
  GWAS_select_t1_b <- reactive({
    mpp.pp = tt()$mpp.pp
    GWAS_select_t1_b_cs()[GWAS_select_t1_b_cs()[, paste0("MPP_",names(mpp.pp$MPP)[1],"_Multi")]>=input$control_widgets_2[1] 
                        & GWAS_select_t1_b_cs()[, paste0("MPP_",names(mpp.pp$MPP)[1],"_Multi")]<=input$control_widgets_2[2]+0.001, ]
  })
  
  #####GWAS table trait 2a single----  
  GWAS_select_t2_a_cs <- reactive({ 
    GWAS_all_final = tt()$GWAS_all_final
    if (input$control_widgets_1 == 0){
      GWAS_all_final[GWAS_all_final$cs_trait_2==0, ]
    }else if (input$control_widgets_1 == 1){
      GWAS_all_final[GWAS_all_final$cs1_2==1, ]
    }else if (input$control_widgets_1 == 2){
      GWAS_all_final[GWAS_all_final$cs_trait_2==2, ]
    }else if (input$control_widgets_1 == 1.1){
      GWAS_all_final[GWAS_all_final$cs_trait_2==1 & GWAS_all_final$cs1_2==1, ]
    }else if (input$control_widgets_1 == 1.2){
      GWAS_all_final[GWAS_all_final$cs_trait_2==1 & GWAS_all_final$csM_2==1, ]  
    }else{
      GWAS_all_final
    }
  })
  GWAS_select_t2_a <- reactive({    
    mpp.pp = tt()$mpp.pp
    GWAS_select_t2_a_cs()[GWAS_select_t2_a_cs()[, paste0("MPP_",names(mpp.pp$MPP)[2],"_Single")]>=input$control_widgets_2[1] 
                        & GWAS_select_t2_a_cs()[, paste0("MPP_",names(mpp.pp$MPP)[2],"_Single")]<=input$control_widgets_2[2]+0.001, ]
  })
  #####GWAS table trait 2b multi----
  GWAS_select_t2_b_cs <- reactive({   
    GWAS_all_final = tt()$GWAS_all_final
    if (input$control_widgets_1 == 0){
      GWAS_all_final[GWAS_all_final$cs_trait_2==0, ]
    }else if (input$control_widgets_1 == 1){
      GWAS_all_final[GWAS_all_final$csM_2==1, ]
    }else if (input$control_widgets_1 == 2){
      GWAS_all_final[GWAS_all_final$cs_trait_2==2, ]
    }else if (input$control_widgets_1 == 1.1){
      GWAS_all_final[GWAS_all_final$cs_trait_2==1 & GWAS_all_final$cs1_2==1, ]
    }else if (input$control_widgets_1 == 1.2){
      GWAS_all_final[GWAS_all_final$cs_trait_2==1 & GWAS_all_final$csM_2==1, ]
    }else{
      GWAS_all_final
    }
  })
  GWAS_select_t2_b <- reactive({  
    mpp.pp = tt()$mpp.pp
    GWAS_select_t2_b_cs()[GWAS_select_t2_b_cs()[, paste0("MPP_",names(mpp.pp$MPP)[2],"_Multi")]>=input$control_widgets_2[1] 
                        & GWAS_select_t2_b_cs()[, paste0("MPP_",names(mpp.pp$MPP)[2],"_Multi")]<=input$control_widgets_2[2]+0.001, ]
  })    
  
  #####GWAS table trait 3a single----  
  GWAS_select_t3_a_cs <- reactive({  
    GWAS_all_final = tt()$GWAS_all_final
    if (input$control_widgets_1 == 0){
      GWAS_all_final[GWAS_all_final$cs_trait_3==0, ]
    }else if (input$control_widgets_1 == 1){
      GWAS_all_final[GWAS_all_final$cs1_3==1, ]
    }else if (input$control_widgets_1 == 2){
      GWAS_all_final[GWAS_all_final$cs_trait_3==2, ]
    }else if (input$control_widgets_1 == 1.1){
      GWAS_all_final[GWAS_all_final$cs_trait_3==1 & GWAS_all_final$cs1_3==1, ]
    }else if (input$control_widgets_1 == 1.2){
      GWAS_all_final[GWAS_all_final$cs_trait_3==1 & GWAS_all_final$csM_3==1, ]  
    }else{
      GWAS_all_final
    }
  })
  GWAS_select_t3_a <- reactive({    
    mpp.pp = tt()$mpp.pp
    GWAS_select_t3_a_cs()[GWAS_select_t3_a_cs()[, paste0("MPP_",names(mpp.pp$MPP)[3],"_Single")]>=input$control_widgets_2[1] 
                        & GWAS_select_t3_a_cs()[, paste0("MPP_",names(mpp.pp$MPP)[3],"_Single")]<=input$control_widgets_2[2]+0.001, ]
  })
  #####GWAS table trait 3b multi----
  GWAS_select_t3_b_cs <- reactive({   
    GWAS_all_final = tt()$GWAS_all_final
    if (input$control_widgets_1 == 0){
      GWAS_all_final[GWAS_all_final$cs_trait_3==0, ]
    }else if (input$control_widgets_1 == 1){
      GWAS_all_final[GWAS_all_final$csM_3==1, ]
    }else if (input$control_widgets_1 == 2){
      GWAS_all_final[GWAS_all_final$cs_trait_3==2, ]
    }else if (input$control_widgets_1 == 1.1){
      GWAS_all_final[GWAS_all_final$cs_trait_3==1 & GWAS_all_final$cs1_3==1, ]
    }else if (input$control_widgets_1 == 1.2){
      GWAS_all_final[GWAS_all_final$cs_trait_3==1 & GWAS_all_final$csM_3==1, ]
    }else{
      GWAS_all_final
    }
  })
  GWAS_select_t3_b <- reactive({  
    mpp.pp = tt()$mpp.pp
    GWAS_select_t3_b_cs()[GWAS_select_t3_b_cs()[, paste0("MPP_",names(mpp.pp$MPP)[3],"_Multi")]>=input$control_widgets_2[1] 
                        & GWAS_select_t3_b_cs()[, paste0("MPP_",names(mpp.pp$MPP)[3],"_Multi")]<=input$control_widgets_2[2]+0.001, ]
  })    
  
  #####GWAS table trait 4a single----  
  GWAS_select_t4_a_cs <- reactive({   
    GWAS_all_final = tt()$GWAS_all_final
    if (input$control_widgets_1 == 0){
      GWAS_all_final[GWAS_all_final$cs_trait_4==0, ]
    }else if (input$control_widgets_1 == 1){
      GWAS_all_final[GWAS_all_final$cs1_4==1, ]
    }else if (input$control_widgets_1 == 2){
      GWAS_all_final[GWAS_all_final$cs_trait_4==2, ]
    }else if (input$control_widgets_1 == 1.1){
      GWAS_all_final[GWAS_all_final$cs_trait_4==1 & GWAS_all_final$cs1_4==1, ]
    }else if (input$control_widgets_1 == 1.2){
      GWAS_all_final[GWAS_all_final$cs_trait_4==1 & GWAS_all_final$csM_4==1, ]  
    }else{
      GWAS_all_final
    }
  })
  GWAS_select_t4_a <- reactive({    
    mpp.pp = tt()$mpp.pp
    GWAS_select_t4_a_cs()[GWAS_select_t4_a_cs()[, paste0("MPP_",names(mpp.pp$MPP)[4],"_Single")]>=input$control_widgets_2[1] 
                        & GWAS_select_t4_a_cs()[, paste0("MPP_",names(mpp.pp$MPP)[4],"_Single")]<=input$control_widgets_2[2]+0.001, ]
  })
  #####GWAS table trait 4b multi----
  GWAS_select_t4_b_cs <- reactive({   
    GWAS_all_final = tt()$GWAS_all_final
    if (input$control_widgets_1 == 0){
      GWAS_all_final[GWAS_all_final$cs_trait_4==0, ]
    }else if (input$control_widgets_1 == 1){
      GWAS_all_final[GWAS_all_final$csM_4==1, ]
    }else if (input$control_widgets_1 == 2){
      GWAS_all_final[GWAS_all_final$cs_trait_4==2, ]
    }else if (input$control_widgets_1 == 1.1){
      GWAS_all_final[GWAS_all_final$cs_trait_4==1 & GWAS_all_final$cs1_4==1, ]
    }else if (input$control_widgets_1 == 1.2){
      GWAS_all_final[GWAS_all_final$cs_trait_4==1 & GWAS_all_final$csM_4==1, ]
    }else{
      GWAS_all_final
    }
  })
  GWAS_select_t4_b <- reactive({  
    mpp.pp = tt()$mpp.pp
    GWAS_select_t4_b_cs()[GWAS_select_t4_b_cs()[, paste0("MPP_",names(mpp.pp$MPP)[4],"_Multi")]>=input$control_widgets_2[1] 
                        & GWAS_select_t4_b_cs()[, paste0("MPP_",names(mpp.pp$MPP)[4],"_Multi")]<=input$control_widgets_2[2]+0.001, ]
  })   
  
  #####GWAS table trait 5a single----  
  GWAS_select_t5_a_cs <- reactive({    
    GWAS_all_final = tt()$GWAS_all_final
    if (input$control_widgets_1 == 0){
      GWAS_all_final[GWAS_all_final$cs_trait_5==0, ]
    }else if (input$control_widgets_1 == 1){
      GWAS_all_final[GWAS_all_final$cs1_5==1, ]
    }else if (input$control_widgets_1 == 2){
      GWAS_all_final[GWAS_all_final$cs_trait_5==2, ]
    }else if (input$control_widgets_1 == 1.1){
      GWAS_all_final[GWAS_all_final$cs_trait_5==1 & GWAS_all_final$cs1_5==1, ]
    }else if (input$control_widgets_1 == 1.2){
      GWAS_all_final[GWAS_all_final$cs_trait_5==1 & GWAS_all_final$csM_5==1, ]  
    }else{
      GWAS_all_final
    }
  })
  GWAS_select_t5_a <- reactive({    
    mpp.pp = tt()$mpp.pp
    GWAS_select_t5_a_cs()[GWAS_select_t5_a_cs()[, paste0("MPP_",names(mpp.pp$MPP)[5],"_Single")]>=input$control_widgets_2[1] 
                        & GWAS_select_t5_a_cs()[, paste0("MPP_",names(mpp.pp$MPP)[5],"_Single")]<=input$control_widgets_2[2]+0.001, ]
  })
  #####GWAS table trait 5b multi----
  GWAS_select_t5_b_cs <- reactive({   
    GWAS_all_final = tt()$GWAS_all_final
    if (input$control_widgets_1 == 0){
      GWAS_all_final[GWAS_all_final$cs_trait_5==0, ]
    }else if (input$control_widgets_1 == 1){
      GWAS_all_final[GWAS_all_final$csM_5==1, ]
    }else if (input$control_widgets_1 == 2){
      GWAS_all_final[GWAS_all_final$cs_trait_5==2, ]
    }else if (input$control_widgets_1 == 1.1){
      GWAS_all_final[GWAS_all_final$cs_trait_5==1 & GWAS_all_final$cs1_5==1, ]
    }else if (input$control_widgets_1 == 1.2){
      GWAS_all_final[GWAS_all_final$cs_trait_5==1 & GWAS_all_final$csM_5==1, ]
    }else{
      GWAS_all_final
    }
  })
  GWAS_select_t5_b <- reactive({  
    mpp.pp = tt()$mpp.pp
    GWAS_select_t5_b_cs()[GWAS_select_t5_b_cs()[, paste0("MPP_",names(mpp.pp$MPP)[5],"_Multi")]>=input$control_widgets_2[1] 
                        & GWAS_select_t5_b_cs()[, paste0("MPP_",names(mpp.pp$MPP)[5],"_Multi")]<=input$control_widgets_2[2]+0.001, ]
  })   
  
  #####GWAS table trait 6a single----  
  GWAS_select_t6_a_cs <- reactive({    
    GWAS_all_final = tt()$GWAS_all_final
    if (input$control_widgets_1 == 0){
      GWAS_all_final[GWAS_all_final$cs_trait_6==0, ]
    }else if (input$control_widgets_1 == 1){
      GWAS_all_final[GWAS_all_final$cs1_6==1, ]
    }else if (input$control_widgets_1 == 2){
      GWAS_all_final[GWAS_all_final$cs_trait_6==2, ]
    }else if (input$control_widgets_1 == 1.1){
      GWAS_all_final[GWAS_all_final$cs_trait_6==1 & GWAS_all_final$cs1_6==1, ]
    }else if (input$control_widgets_1 == 1.2){
      GWAS_all_final[GWAS_all_final$cs_trait_6==1 & GWAS_all_final$csM_6==1, ]  
    }else{
      GWAS_all_final
    }
  })
  GWAS_select_t6_a <- reactive({    
    mpp.pp = tt()$mpp.pp
    GWAS_select_t6_a_cs()[GWAS_select_t6_a_cs()[, paste0("MPP_",names(mpp.pp$MPP)[6],"_Single")]>=input$control_widgets_2[1] 
                        & GWAS_select_t6_a_cs()[, paste0("MPP_",names(mpp.pp$MPP)[6],"_Single")]<=input$control_widgets_2[2]+0.001, ]
  })
  #####GWAS table trait 6b multi----
  GWAS_select_t6_b_cs <- reactive({   
    GWAS_all_final = tt()$GWAS_all_final
    if (input$control_widgets_1 == 0){
      GWAS_all_final[GWAS_all_final$cs_trait_6==0, ]
    }else if (input$control_widgets_1 == 1){
      GWAS_all_final[GWAS_all_final$csM_6==1, ]
    }else if (input$control_widgets_1 == 2){
      GWAS_all_final[GWAS_all_final$cs_trait_6==2, ]
    }else if (input$control_widgets_1 == 1.1){
      GWAS_all_final[GWAS_all_final$cs_trait_6==1 & GWAS_all_final$cs1_6==1, ]
    }else if (input$control_widgets_1 == 1.2){
      GWAS_all_final[GWAS_all_final$cs_trait_6==1 & GWAS_all_final$csM_6==1, ]
    }else{
      GWAS_all_final
    }
    
  })
  GWAS_select_t6_b <- reactive({  
    mpp.pp = tt()$mpp.pp
    GWAS_select_t6_b_cs()[GWAS_select_t6_b_cs()[, paste0("MPP_",names(mpp.pp$MPP)[6],"_Multi")]>=input$control_widgets_2[1] 
                        & GWAS_select_t6_b_cs()[, paste0("MPP_",names(mpp.pp$MPP)[6],"_Multi")]<=input$control_widgets_2[2]+0.001, ]
  })   
  
  
  #ch0: Regional association_plot_1a----
  output$manhattan_1a <- renderPlotly({    
    #withProgress_2
    withProgress(message = 'Calculation in progress',
                 detail = '...done: manhattan_1a', value = 0, {
                   for (i in 1:10) {
                     incProgress(1/10)
                     Sys.sleep(0.05)
                   }
                 })
    
    mpp.pp = tt()$mpp.pp
    pal = tt()$pal
    
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
           line = list(color = color, dash = 'dash'))}
    
    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
           line = list(color = color))}
    
    fig_1a <- plot_ly(GWAS_select_t1_a(), x = ~ps, y = ~pval_log_1, 
                      type = 'scatter', mode = 'markers', 
                      #size = ~MPP_lowdlipo_Single, color = ~snpGroups_fm, 
                      size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[1],"_Single"))), 
                      color = ~snpGroups_fm,
                      colors = pal, 
                      sizes = c(5, 50),
                      fill = ~'',
                      marker = list(opacity = 1, sizemode = 'diameter'),
                      hoverinfo = 'text',
                      text = ~paste('SNP_ID: ', rs,
                                    '<br>SNP_position: ', ps,
                                    '<br>Allele_1: ', allele1,
                                    '<br>Allele_0 (or Allele_2): ', allele0,
                                    '<br>af: ', signif(af,4),
                                    #'<br>ps:', ps,
                                    '<br>pval: ', signif(pval_1,4),
                                    #'<br>-log10(pval):', signif(pval_log_1,4),
                                    '<br>MPP: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[1],"_Single"))),4),
                                    '<br>MPPg: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[1],"_Single"))),4), " (if MPPg=0: negligible)",
                                    '<br>snpGroups_fm: ', snpGroups_fm, ' (',snpGroups_fm_size, ' SNPs)',
                                    '<br>snpGroups_flashfm: ', snpGroups_flashfm, ' (',snpGroups_flashfm_size, ' SNPs)',
                                    '<br>Credible_set_cs1: ', cs1_1, 
                                    '<br>Credible_set_csM: ', csM_1,
                                    sep = ""))
    fig_1a <- fig_1a %>% layout(#title = 'Interactive Regional association plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             title = paste0('Chromosome_', 
                                                            GWAS_select_t1_a()$chr[1], 
                                                            '_position \n(NOTE: Size = MPP, Color = snpGroups)'),
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_select_t1_a()$ps), 
                                                       max(GWAS_select_t1_a()$ps)),
                                            #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0-10, max(GWAS_select_t1_a()$pval_log_1)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                showlegend = TRUE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    
    M_ch0_a1=input$control_widgets_ch_0_a1  ##CHECK
    if (M_ch0_a1 == 0){
      fig_1a <- fig_1a %>% layout(shapes = list(vline(GWAS_select_t1_a()$ps[which(GWAS_select_t1_a()$pval_log_1==max(GWAS_select_t1_a()$pval_log_1))]), hline(7.3)))
    }else{
      fig_1a <- fig_1a %>% layout(shapes = list(hline(7.3)))
    }
    #fig_1a <- fig_1a %>% layout(shapes = list(vline(GWAS_select_t1_a()$ps[which(GWAS_select_t1_a()$pval_log_1==max(GWAS_select_t1_a()$pval_log_1))]), hline(7.3)))
    
    # fig_1a <- fig_1a %>% layout(annotations = list( 
    # list(x = 0.001 , y = -0.05, align = 'left', 
    #      text = "Figure trait_1_FINEMAP: Size = MPP, Color = snpGroups", 
    #      showarrow = F, xref='paper', yref='paper')) )
    fig_1a
    })
  
  #ch0: Manhattan_plot_1b----
  output$manhattan_1b <- renderPlotly({
    #withProgress_3
    withProgress(message = 'Calculation in progress',
                 detail = '...done: manhattan_1b', value = 0, {
                   for (i in 1:10) {
                     incProgress(1/10)
                     Sys.sleep(0.05)
                   }
                 })
    
    mpp.pp = tt()$mpp.pp
    pal = tt()$pal
    
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
           line = list(color = color, dash = 'dash'))}

    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
           line = list(color = color))}

    fig_1b <- plot_ly(GWAS_select_t1_b(), x = ~ps, y = ~pval_log_1, 
                      type = 'scatter', mode = 'markers', 
                      size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[1],"_Multi"))), 
                      color = ~snpGroups_flashfm, 
                      colors = pal, 
                      sizes = c(5, 50),
                      fill = ~'',
                      marker = list(opacity = 1, sizemode = 'diameter'),
                      hoverinfo = 'text',
                      text = ~paste('SNP_ID: ', rs,
                                    '<br>SNP_position: ', ps,
                                    '<br>Allele_1: ', allele1,
                                    '<br>Allele_0 (or Allele_2): ', allele0,
                                    '<br>af: ', signif(af,4),
                                    #'<br>ps:', ps,
                                    '<br>pval: ', signif(pval_1,4),
                                    #'<br>-log10(pval):', signif(pval_log_1,4),
                                    '<br>MPP: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[1],"_Multi"))),4),
                                    '<br>MPPg: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[1],"_Multi"))),4), " (if MPPg=0: negligible)",
                                    '<br>snpGroups_fm: ', snpGroups_fm, ' (',snpGroups_fm_size, ' SNPs)',
                                    '<br>snpGroups_flashfm: ', snpGroups_flashfm, ' (',snpGroups_flashfm_size, ' SNPs)',
                                    '<br>Credible_set_cs1: ', cs1_1, 
                                    '<br>Credible_set_csM: ', csM_1,
                                    sep = ""))
    
    fig_1b <- fig_1b %>% layout(#title = 'Interactive Regional association plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                            title = paste0('Chromosome_', 
                                                            GWAS_select_t1_b()$chr[1], 
                                                           '_position \n(NOTE: Size = MPP, Color = snpGroups)'),
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_select_t1_b()$ps), 
                                                       max(GWAS_select_t1_b()$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0-10, max(GWAS_select_t1_b()$pval_log_1)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                showlegend = TRUE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    
    M_ch0_b1=input$control_widgets_ch_0_b1  ##CHECK
    if (M_ch0_b1 == 0){
      fig_1b <- fig_1b %>% layout(shapes = list(vline(GWAS_select_t1_b()$ps[which(GWAS_select_t1_b()$pval_log_1==max(GWAS_select_t1_b()$pval_log_1))]), hline(7.3)))
    }else{
      fig_1b <- fig_1b %>% layout(shapes = list(hline(7.3)))
    }
    # fig_1b <- fig_1b %>% layout(annotations = list( 
    # list(x = 0.05 , y = 0.95, align = 'left', 
    #      text = "Figure trait_1_flashfm: \nY-axis = pval_lowdlipo, \nX-axis = ps, \nSize = MPP_lowdlipo_Multi, \nColor = snpGroups_flashfm", 
    #      showarrow = F, xref='paper', yref='paper')) )
    fig_1b
  })  
  
  #ch0: Manhattan_plot_2a----
  output$manhattan_2a <- renderPlotly({    
    mpp.pp = tt()$mpp.pp
    pal = tt()$pal
    
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
           line = list(color = color, dash = 'dash'))}
    
    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
           line = list(color = color))}
    
    fig_2a <- plot_ly(GWAS_select_t2_a(), x = ~ps, y = ~pval_log_2, 
                      type = 'scatter', mode = 'markers', 
                      size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[2],"_Single"))), 
                      color = ~snpGroups_fm, 
                      colors = pal, 
                      sizes = c(5, 50),
                      fill = ~'',
                      marker = list(opacity = 1, sizemode = 'diameter'),
                      hoverinfo = 'text',
                      text = ~paste('SNP_ID: ', rs,
                                    '<br>SNP_position: ', ps,
                                    '<br>Allele_1: ', allele1,
                                    '<br>Allele_0 (or Allele_2): ', allele0,
                                    '<br>af: ', signif(af,4),
                                    #'<br>ps:', ps,
                                    '<br>pval: ', signif(pval_2,4),
                                    #'<br>-log10(pval):', signif(pval_log_2,4),
                                    '<br>MPP: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[2],"_Single"))),4),
                                    '<br>MPPg: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[2],"_Single"))),4), " (if MPPg=0: negligible)",
                                    '<br>snpGroups_fm: ', snpGroups_fm, ' (',snpGroups_fm_size, ' SNPs)',
                                    '<br>snpGroups_flashfm: ', snpGroups_flashfm, ' (',snpGroups_flashfm_size, ' SNPs)',
                                    '<br>Credible_set_cs1: ', cs1_2, 
                                    '<br>Credible_set_csM: ', csM_2,
                                    sep = ""))
    fig_2a <- fig_2a %>% layout(#title = 'Interactive Regional association plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             title = paste0('Chromosome_', 
                                                             GWAS_select_t2_a()$chr[1], 
                                                            '_position \n(NOTE: Size = MPP, Color = snpGroups)'),
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_select_t2_a()$ps), 
                                                       max(GWAS_select_t2_a()$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0-10, max(GWAS_select_t2_a()$pval_log_2)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                showlegend = TRUE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    
    M_ch0_a2=input$control_widgets_ch_0_a2  ##CHECK
    if (M_ch0_a2 == 0){
      fig_2a <- fig_2a %>% layout(shapes = list(vline(GWAS_select_t2_a()$ps[which(GWAS_select_t2_a()$pval_log_2==max(GWAS_select_t2_a()$pval_log_2))]), hline(7.3)))
    }else{
      fig_2a <- fig_2a %>% layout(shapes = list(hline(7.3)))
    }
    # fig_2a <- fig_2a %>% layout(annotations = list( 
    #   list(x = 0.05 , y = 0.95, align = 'left', 
    #        text = "Figure trait_2_FINEMAP: \nY-axis = pval_cholesterol, \nX-axis = ps, \nSize = MPP_cholesterol_Single, \nColor = snpGroups_fm", 
    #        showarrow = F, xref='paper', yref='paper')) )
    fig_2a
  })  
  
  #ch0: Manhattan_plot_2b----
  output$manhattan_2b <- renderPlotly({  
    mpp.pp = tt()$mpp.pp
    pal = tt()$pal
    
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
           line = list(color = color, dash = 'dash'))}
    
    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
           line = list(color = color))}
    
    fig_2b <- plot_ly(GWAS_select_t2_b(), x = ~ps, y = ~pval_log_2, 
                      type = 'scatter', mode = 'markers', 
                      size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[2],"_Multi"))), 
                      color = ~snpGroups_flashfm, 
                      colors = pal, 
                      sizes = c(5, 50),
                      fill = ~'',
                      marker = list(opacity = 1, sizemode = 'diameter'),
                      hoverinfo = 'text',
                      text = ~paste('SNP_ID: ', rs,
                                    '<br>SNP_position: ', ps,
                                    '<br>Allele_1: ', allele1,
                                    '<br>Allele_0 (or Allele_2): ', allele0,
                                    '<br>af: ', signif(af,4),
                                    #'<br>ps:', ps,
                                    '<br>pval: ', signif(pval_2,4),
                                    #'<br>-log10(pval):', signif(pval_log_2,4),
                                    '<br>MPP: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[2],"_Multi"))),4),
                                    '<br>MPPg: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[2],"_Multi"))),4), " (if MPPg=0: negligible)",
                                    '<br>snpGroups_fm: ', snpGroups_fm, ' (',snpGroups_fm_size, ' SNPs)',
                                    '<br>snpGroups_flashfm: ', snpGroups_flashfm, ' (',snpGroups_flashfm_size, ' SNPs)',
                                    '<br>Credible_set_cs1: ', cs1_2, 
                                    '<br>Credible_set_csM: ', csM_2,
                                    sep = ""))
    fig_2b <- fig_2b %>% layout(#title = 'Interactive Regional association plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             title = paste0('Chromosome_', 
                                                             GWAS_select_t2_b()$chr[1], 
                                                            '_position \n(NOTE: Size = MPP, Color = snpGroups)'),
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_select_t2_b()$ps), 
                                                       max(GWAS_select_t2_b()$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0-10, max(GWAS_select_t2_b()$pval_log_2)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                showlegend = TRUE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    
    M_ch0_b2=input$control_widgets_ch_0_b2  ##CHECK
    if (M_ch0_b2 == 0){
      fig_2b <- fig_2b %>% layout(shapes = list(vline(GWAS_select_t2_b()$ps[which(GWAS_select_t2_b()$pval_log_2==max(GWAS_select_t2_b()$pval_log_2))]), hline(7.3)))
    }else{
      fig_2b <- fig_2b %>% layout(shapes = list(hline(7.3)))
    }
    # fig_2b <- fig_2b %>% layout(annotations = list( 
    #   list(x = 0.05 , y = 0.95, align = 'left', 
    #        text = "Figure trait_2_flashfm: \nY-axis = pval_cholesterol, \nX-axis = ps, \nSize = MPP_cholesterol_Multi, \nColor = snpGroups_flashfm", 
    #        showarrow = F, xref='paper', yref='paper')) )
    fig_2b
  })  

  #ch0: Manhattan_plot_3a----
  output$manhattan_3a <- renderPlotly({   
    mpp.pp = tt()$mpp.pp
    pal = tt()$pal
    
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
           line = list(color = color, dash = 'dash'))}
    
    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
           line = list(color = color))}
    
    fig_3a <- plot_ly(GWAS_select_t3_a(), x = ~ps, y = ~pval_log_3, 
                      type = 'scatter', mode = 'markers', 
                      size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[3],"_Single"))), 
                      color = ~snpGroups_fm, 
                      colors = pal, 
                      sizes = c(5, 50),
                      fill = ~'',
                      marker = list(opacity = 1, sizemode = 'diameter'),
                      hoverinfo = 'text',
                      text = ~paste('SNP_ID: ', rs,
                                    '<br>SNP_position: ', ps,
                                    '<br>Allele_1: ', allele1,
                                    '<br>Allele_0 (or Allele_2): ', allele0,
                                    '<br>af: ', signif(af,4),
                                    #'<br>ps:', ps,
                                    '<br>pval: ', signif(pval_3,4),
                                    #'<br>-log10(pval):', signif(pval_log_3,4),
                                    '<br>MPP: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[3],"_Single"))),4),
                                    '<br>MPPg: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[3],"_Single"))),4), " (if MPPg=0: negligible)",
                                    '<br>snpGroups_fm: ', snpGroups_fm, ' (',snpGroups_fm_size, ' SNPs)',
                                    '<br>snpGroups_flashfm: ', snpGroups_flashfm, ' (',snpGroups_flashfm_size, ' SNPs)',
                                    '<br>Credible_set_cs1: ', cs1_3, 
                                    '<br>Credible_set_csM: ', csM_3,
                                    sep = ""))
    fig_3a <- fig_3a %>% layout(#title = 'Interactive Regional association plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             title = paste0('Chromosome_', 
                                                             GWAS_select_t3_a()$chr[1], 
                                                            '_position \n(NOTE: Size = MPP, Color = snpGroups)'),
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_select_t3_a()$ps), 
                                                       max(GWAS_select_t3_a()$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0-10, max(GWAS_select_t3_a()$pval_log_3)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                showlegend = TRUE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    
    #fig_3a <- fig_3a %>% layout(shapes = list(vline(GWAS_select_t3_a()$ps[which(GWAS_select_t3_a()$pval_log_3==max(GWAS_select_t3_a()$pval_log_3))]), hline(7.3)))
    fig_3a <- fig_3a %>% layout(shapes = list(hline(7.3)))
    # fig_3a <- fig_3a %>% layout(annotations = list( 
    #   list(x = 0.05 , y = 0.95, align = 'left', 
    #        text = "Figure trait_3_FINEMAP: \nY-axis = pval_triglycerides, \nX-axis = ps, \nSize = MPP_triglycerides_Single, \nColor = snpGroups_fm", 
    #        showarrow = F, xref='paper', yref='paper')) )
    fig_3a
  })      

  #ch0: Manhattan_plot_3b----
  output$manhattan_3b <- renderPlotly({    
    mpp.pp = tt()$mpp.pp
    pal = tt()$pal
    
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
           line = list(color = color, dash = 'dash'))}
    
    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
           line = list(color = color))}
    
    fig_3b <- plot_ly(GWAS_select_t3_b(), x = ~ps, y = ~pval_log_3, 
                      type = 'scatter', mode = 'markers', 
                      size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[3],"_Multi"))), 
                      color = ~snpGroups_flashfm, 
                      colors = pal, 
                      sizes = c(5, 50),
                      fill = ~'',
                      marker = list(opacity = 1, sizemode = 'diameter'),
                      hoverinfo = 'text',
                      text = ~paste('SNP_ID: ', rs,
                                    '<br>SNP_position: ', ps,
                                    '<br>Allele_1: ', allele1,
                                    '<br>Allele_0 (or Allele_2): ', allele0,
                                    '<br>af: ', signif(af,4),
                                    #'<br>ps:', ps,
                                    '<br>pval: ', signif(pval_3,4),
                                    #'<br>-log10(pval):', signif(pval_log_3,4),
                                    '<br>MPP: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[3],"_Multi"))),4),
                                    '<br>MPPg: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[3],"_Multi"))),4), " (if MPPg=0: negligible)",
                                    '<br>snpGroups_fm: ', snpGroups_fm, ' (',snpGroups_fm_size, ' SNPs)',
                                    '<br>snpGroups_flashfm: ', snpGroups_flashfm, ' (',snpGroups_flashfm_size, ' SNPs)',
                                    '<br>Credible_set_cs1: ', cs1_3, 
                                    '<br>Credible_set_csM: ', csM_3,
                                    sep = ""))
    fig_3b <- fig_3b %>% layout(#title = 'Interactive Regional association plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             title = paste0('Chromosome_', 
                                                             GWAS_select_t3_b()$chr[1], 
                                                            '_position \n(NOTE: Size = MPP, Color = snpGroups)'),
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_select_t3_b()$ps), 
                                                       max(GWAS_select_t3_b()$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0-10, max(GWAS_select_t3_b()$pval_log_3)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                showlegend = TRUE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    #fig_3b <- fig_3b %>% layout(shapes = list(vline(GWAS_select_t3_b()$ps[which(GWAS_select_t3_b()$pval_log_3==max(GWAS_select_t3_b()$pval_log_3))]), hline(7.3)))
    fig_3b <- fig_3b %>% layout(shapes = list(hline(7.3)))
    # fig_3b <- fig_3b %>% layout(annotations = list( 
    #   list(x = 0.05 , y = 0.95, align = 'left', 
    #        text = "Figure trait_3_flashfm: \nY-axis = pval_triglycerides, \nX-axis = ps, \nSize = MPP_triglycerides_Multi, \nColor = snpGroups_flashfm", 
    #        showarrow = F, xref='paper', yref='paper')) )
    fig_3b
  })   

  #ch0: Manhattan_plot_4a----
  output$manhattan_4a <- renderPlotly({    
    mpp.pp = tt()$mpp.pp
    pal = tt()$pal
    
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
           line = list(color = color, dash = 'dash'))}
    
    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
           line = list(color = color))}
    
    fig_4a <- plot_ly(GWAS_select_t4_a(), x = ~ps, y = ~pval_log_4, 
                      type = 'scatter', mode = 'markers', 
                      size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[4],"_Single"))), 
                      color = ~snpGroups_fm, 
                      colors = pal, 
                      sizes = c(5, 50),
                      fill = ~'',
                      marker = list(opacity = 1, sizemode = 'diameter'),
                      hoverinfo = 'text',
                      text = ~paste('SNP_ID: ', rs,
                                    '<br>SNP_position: ', ps,
                                    '<br>Allele_1: ', allele1,
                                    '<br>Allele_0 (or Allele_2): ', allele0,
                                    '<br>af: ', signif(af,4),
                                    #'<br>ps:', ps,
                                    '<br>pval: ', signif(pval_4,4),
                                    #'<br>-log10(pval):', signif(pval_log_4,4),
                                    '<br>MPP: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[4],"_Single"))),4),
                                    '<br>MPPg: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[4],"_Single"))),4), " (if MPPg=0: negligible)",
                                    '<br>snpGroups_fm: ', snpGroups_fm, ' (',snpGroups_fm_size, ' SNPs)',
                                    '<br>snpGroups_flashfm: ', snpGroups_flashfm, ' (',snpGroups_flashfm_size, ' SNPs)',
                                    '<br>Credible_set_cs1: ', cs1_4, 
                                    '<br>Credible_set_csM: ', csM_4,
                                    sep = ""))
    fig_4a <- fig_4a %>% layout(#title = 'Interactive Regional association plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             title = paste0('Chromosome_', 
                                                             GWAS_select_t4_a()$chr[1], 
                                                            '_position \n(NOTE: Size = MPP, Color = snpGroups)'),
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_select_t4_a()$ps), 
                                                       max(GWAS_select_t4_a()$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0-10, max(GWAS_select_t4_a()$pval_log_4)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                showlegend = TRUE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    #fig_4a <- fig_4a %>% layout(shapes = list(vline(GWAS_select_t4_a()$ps[which(GWAS_select_t4_a()$pval_log_4==max(GWAS_select_t4_a()$pval_log_4))]), hline(7.3)))
    fig_4a <- fig_4a %>% layout(shapes = list(hline(7.3)))
    # fig_4a <- fig_4a %>% layout(annotations = list( 
    #   list(x = 0.05 , y = 0.95, align = 'left', 
    #        text = "Figure trait_4_FINEMAP: \nY-axis = pval_highdlipo, \nX-axis = ps, \nSize = MPP_highdlipo_Single, \nColor = snpGroups_fm", 
    #        showarrow = F, xref='paper', yref='paper')) )
    fig_4a
  })   

  #ch0: Manhattan_plot_4b----
  output$manhattan_4b <- renderPlotly({    
    mpp.pp = tt()$mpp.pp
    pal = tt()$pal
    
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
           line = list(color = color, dash = 'dash'))}
    
    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
           line = list(color = color))}
    
    fig_4b <- plot_ly(GWAS_select_t4_b(), x = ~ps, y = ~pval_log_4, 
                      type = 'scatter', mode = 'markers', 
                      size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[4],"_Multi"))), 
                      color = ~snpGroups_flashfm, 
                      colors = pal, 
                      sizes = c(5, 50),
                      fill = ~'',
                      marker = list(opacity = 1, sizemode = 'diameter'),
                      hoverinfo = 'text',
                      text = ~paste('SNP_ID: ', rs,
                                    '<br>SNP_position: ', ps,
                                    '<br>Allele_1: ', allele1,
                                    '<br>Allele_0 (or Allele_2): ', allele0,
                                    '<br>af: ', signif(af,4),
                                    #'<br>ps:', ps,
                                    '<br>pval: ', signif(pval_4,4),
                                    #'<br>-log10(pval):', signif(pval_log_4,4),
                                    '<br>MPP: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[4],"_Multi"))),4),
                                    '<br>MPPg: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[4],"_Multi"))),4), " (if MPPg=0: negligible)",
                                    '<br>snpGroups_fm: ', snpGroups_fm, ' (',snpGroups_fm_size, ' SNPs)',
                                    '<br>snpGroups_flashfm: ', snpGroups_flashfm, ' (',snpGroups_flashfm_size, ' SNPs)',
                                    '<br>Credible_set_cs1: ', cs1_4, 
                                    '<br>Credible_set_csM: ', csM_4,
                                    sep = ""))
    fig_4b <- fig_4b %>% layout(#title = 'Interactive Regional association plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             title = paste0('Chromosome_', 
                                                             GWAS_select_t4_b()$chr[1], 
                                                            '_position \n(NOTE: Size = MPP, Color = snpGroups)'),
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_select_t4_b()$ps), 
                                                       max(GWAS_select_t4_b()$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0-10, max(GWAS_select_t4_b()$pval_log_4)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                showlegend = TRUE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    #fig_4b <- fig_4b %>% layout(shapes = list(vline(GWAS_select_t4_b()$ps[which(GWAS_select_t4_b()$pval_log_4==max(GWAS_select_t4_b()$pval_log_4))]), hline(7.3)))
    fig_4b <- fig_4b %>% layout(shapes = list(hline(7.3)))
    # fig_4b <- fig_4b %>% layout(annotations = list( 
    #   list(x = 0.05 , y = 0.95, align = 'left', 
    #        text = "Figure trait_4_flashfm: \nY-axis = pval_highdlipo, \nX-axis = ps, \nSize = MPP_highdlipo_Multi, \nColor = snpGroups_flashfm", 
    #        showarrow = F, xref='paper', yref='paper')) )
    fig_4b
  })    
  
  #ch0: Manhattan_plot_5a----
  output$manhattan_5a <- renderPlotly({    
    mpp.pp = tt()$mpp.pp
    pal = tt()$pal
    
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
           line = list(color = color, dash = 'dash'))}
    
    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
           line = list(color = color))}
    
    fig_5a <- plot_ly(GWAS_select_t5_a(), x = ~ps, y = ~pval_log_5, 
                      type = 'scatter', mode = 'markers', 
                      size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[5],"_Single"))), 
                      color = ~snpGroups_fm, 
                      colors = pal, 
                      sizes = c(5, 50),
                      fill = ~'',
                      marker = list(opacity = 1, sizemode = 'diameter'),
                      hoverinfo = 'text',
                      text = ~paste('SNP_ID: ', rs,
                                    '<br>SNP_position: ', ps,
                                    '<br>Allele_1: ', allele1,
                                    '<br>Allele_0 (or Allele_2): ', allele0,
                                    '<br>af: ', signif(af,4),
                                    #'<br>ps:', ps,
                                    '<br>pval: ', signif(pval_5,4),
                                    #'<br>-log10(pval):', signif(pval_log_5,4),
                                    '<br>MPP: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[5],"_Single"))),4),
                                    '<br>MPPg: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[5],"_Single"))),4), " (if MPPg=0: negligible)",
                                    '<br>snpGroups_fm: ', snpGroups_fm, ' (',snpGroups_fm_size, ' SNPs)',
                                    '<br>snpGroups_flashfm: ', snpGroups_flashfm, ' (',snpGroups_flashfm_size, ' SNPs)',
                                    '<br>Credible_set_cs1: ', cs1_5, 
                                    '<br>Credible_set_csM: ', csM_5,
                                    sep = ""))
    fig_5a <- fig_5a %>% layout(#title = 'Interactive Regional association plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             title = paste0('Chromosome_', 
                                                             GWAS_select_t5_a()$chr[1], 
                                                            '_position \n(NOTE: Size = MPP, Color = snpGroups)'),
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_select_t5_a()$ps), 
                                                       max(GWAS_select_t5_a()$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0-10, max(GWAS_select_t5_a()$pval_log_5)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                showlegend = TRUE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    #fig_5a <- fig_5a %>% layout(shapes = list(vline(GWAS_select_t5_a()$ps[which(GWAS_select_t5_a()$pval_log_5==max(GWAS_select_t5_a()$pval_log_5))]), hline(7.3)))
    fig_5a <- fig_5a %>% layout(shapes = list(hline(7.3)))
    # fig_5a <- fig_5a %>% layout(annotations = list( 
    #   list(x = 0.05 , y = 0.95, align = 'left', 
    #        text = "Figure trait_5_FINEMAP: \nY-axis = pval_highdlipo, \nX-axis = ps, \nSize = MPP_highdlipo_Single, \nColor = snpGroups_fm", 
    #        showarrow = F, xref='paper', yref='paper')) )
    fig_5a
  })   
  
  #ch0: Manhattan_plot_5b----
  output$manhattan_5b <- renderPlotly({    
    mpp.pp = tt()$mpp.pp
    pal = tt()$pal
    
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
           line = list(color = color, dash = 'dash'))}
    
    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
           line = list(color = color))}
    
    fig_5b <- plot_ly(GWAS_select_t5_b(), x = ~ps, y = ~pval_log_5, 
                      type = 'scatter', mode = 'markers', 
                      size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[5],"_Multi"))), 
                      color = ~snpGroups_flashfm, 
                      colors = pal, 
                      sizes = c(5, 50),
                      fill = ~'',
                      marker = list(opacity = 1, sizemode = 'diameter'),
                      hoverinfo = 'text',
                      text = ~paste('SNP_ID: ', rs,
                                    '<br>SNP_position: ', ps,
                                    '<br>Allele_1: ', allele1,
                                    '<br>Allele_0 (or Allele_2): ', allele0,
                                    '<br>af: ', signif(af,4),
                                    #'<br>ps:', ps,
                                    '<br>pval: ', signif(pval_5,4),
                                    #'<br>-log10(pval):', signif(pval_log_5,4),
                                    '<br>MPP: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[5],"_Multi"))),4),
                                    '<br>MPPg: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[5],"_Multi"))),4), " (if MPPg=0: negligible)",
                                    '<br>snpGroups_fm: ', snpGroups_fm, ' (',snpGroups_fm_size, ' SNPs)',
                                    '<br>snpGroups_flashfm: ', snpGroups_flashfm, ' (',snpGroups_flashfm_size, ' SNPs)',
                                    '<br>Credible_set_cs1: ', cs1_5, 
                                    '<br>Credible_set_csM: ', csM_5,
                                    sep = ""))
    fig_5b <- fig_5b %>% layout(#title = 'Interactive Regional association plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             title = paste0('Chromosome_', 
                                                             GWAS_select_t5_b()$chr[1], 
                                                            '_position \n(NOTE: Size = MPP, Color = snpGroups)'),
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_select_t5_b()$ps), 
                                                       max(GWAS_select_t5_b()$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0-10, max(GWAS_select_t5_b()$pval_log_5)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                showlegend = TRUE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    #fig_5b <- fig_5b %>% layout(shapes = list(vline(GWAS_select_t5_b()$ps[which(GWAS_select_t5_b()$pval_log_5==max(GWAS_select_t5_b()$pval_log_5))]), hline(7.3)))
    fig_5b <- fig_5b %>% layout(shapes = list(hline(7.3)))
    # fig_5b <- fig_5b %>% layout(annotations = list( 
    #   list(x = 0.05 , y = 0.95, align = 'left', 
    #        text = "Figure trait_5_flashfm: \nY-axis = pval_highdlipo, \nX-axis = ps, \nSize = MPP_highdlipo_Multi, \nColor = snpGroups_flashfm", 
    #        showarrow = F, xref='paper', yref='paper')) )
    fig_5b
  })    
  
  
  #ch0: Manhattan_plot_6a----
  output$manhattan_6a <- renderPlotly({    
    mpp.pp = tt()$mpp.pp
    pal = tt()$pal
    
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
           line = list(color = color, dash = 'dash'))}
    
    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
           line = list(color = color))}
    
    fig_6a <- plot_ly(GWAS_select_t6_a(), x = ~ps, y = ~pval_log_6, 
                      type = 'scatter', mode = 'markers', 
                      size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[6],"_Single"))), 
                      color = ~snpGroups_fm, 
                      colors = pal, 
                      sizes = c(5, 50),
                      fill = ~'',
                      marker = list(opacity = 1, sizemode = 'diameter'),
                      hoverinfo = 'text',
                      text = ~paste('SNP_ID: ', rs,
                                    '<br>SNP_position: ', ps,
                                    '<br>Allele_1: ', allele1,
                                    '<br>Allele_0 (or Allele_2): ', allele0,
                                    '<br>af: ', signif(af,4),
                                    #'<br>ps:', ps,
                                    '<br>pval: ', signif(pval_6,4),
                                    #'<br>-log10(pval):', signif(pval_log_6,4),
                                    '<br>MPP: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[6],"_Single"))),4),
                                    '<br>MPPg: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[6],"_Single"))),4), " (if MPPg=0: negligible)",
                                    '<br>snpGroups_fm: ', snpGroups_fm, ' (',snpGroups_fm_size, ' SNPs)',
                                    '<br>snpGroups_flashfm: ', snpGroups_flashfm, ' (',snpGroups_flashfm_size, ' SNPs)',
                                    '<br>Credible_set_cs1: ', cs1_6, 
                                    '<br>Credible_set_csM: ', csM_6,
                                    sep = ""))
    fig_6a <- fig_6a %>% layout(#title = 'Interactive Regional association plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             title = paste0('Chromosome_', 
                                                             GWAS_select_t6_a()$chr[1], 
                                                            '_position \n(NOTE: Size = MPP, Color = snpGroups)'),
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_select_t6_a()$ps), 
                                                       max(GWAS_select_t6_a()$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0-10, max(GWAS_select_t6_a()$pval_log_6)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                showlegend = TRUE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    #fig_6a <- fig_6a %>% layout(shapes = list(vline(GWAS_select_t6_a()$ps[which(GWAS_select_t6_a()$pval_log_6==max(GWAS_select_t6_a()$pval_log_6))]), hline(7.3)))
    fig_6a <- fig_6a %>% layout(shapes = list(hline(7.3)))
    # fig_6a <- fig_6a %>% layout(annotations = list( 
    #   list(x = 0.05 , y = 0.95, align = 'left', 
    #        text = "Figure trait_6_FINEMAP: \nY-axis = pval_highdlipo, \nX-axis = ps, \nSize = MPP_highdlipo_Single, \nColor = snpGroups_fm", 
    #        showarrow = F, xref='paper', yref='paper')) )
    fig_6a
  })   
  
  #ch0: Manhattan_plot_6b----
  output$manhattan_6b <- renderPlotly({    
    mpp.pp = tt()$mpp.pp
    pal = tt()$pal
    
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x,
           line = list(color = color, dash = 'dash'))}
    
    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y,
           line = list(color = color))}
    
    fig_6b <- plot_ly(GWAS_select_t6_b(), x = ~ps, y = ~pval_log_6, 
                      type = 'scatter', mode = 'markers', 
                      size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[1],"_Multi"))), 
                      color = ~snpGroups_flashfm, 
                      colors = pal, 
                      sizes = c(5, 50),
                      fill = ~'',
                      marker = list(opacity = 1, sizemode = 'diameter'),
                      hoverinfo = 'text',
                      text = ~paste('SNP_ID: ', rs,
                                    '<br>SNP_position: ', ps,
                                    '<br>Allele_1: ', allele1,
                                    '<br>Allele_0 (or Allele_2): ', allele0,
                                    '<br>af: ', signif(af,4),
                                    #'<br>ps:', ps,
                                    '<br>pval: ', signif(pval_6,4),
                                    #'<br>-log10(pval):', signif(pval_log_6,4),
                                    '<br>MPP: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[6],"_Multi"))),4),
                                    '<br>MPPg: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[6],"_Multi"))),4), " (if MPPg=0: negligible)",
                                    '<br>snpGroups_fm: ', snpGroups_fm, ' (',snpGroups_fm_size, ' SNPs)',
                                    '<br>snpGroups_flashfm: ', snpGroups_flashfm, ' (',snpGroups_flashfm_size, ' SNPs)',
                                    '<br>Credible_set_cs1: ', cs1_6, 
                                    '<br>Credible_set_csM: ', csM_6,
                                    sep = ""))
    fig_6b <- fig_6b %>% layout(#title = 'Interactive Regional association plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             title = paste0('Chromosome_', 
                                                             GWAS_select_t6_b()$chr[1], 
                                                            '_position \n(NOTE: Size = MPP, Color = snpGroups)'),
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_select_t6_b()$ps), 
                                                       max(GWAS_select_t6_b()$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0-10, max(GWAS_select_t6_b()$pval_log_6)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                showlegend = TRUE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    #fig_6b <- fig_6b %>% layout(shapes = list(vline(GWAS_select_t6_b()$ps[which(GWAS_select_t6_b()$pval_log_6==max(GWAS_select_t6_b()$pval_log_6))]), hline(7.3)))
    fig_6b <- fig_6b %>% layout(shapes = list(hline(7.3)))
    # fig_6b <- fig_6b %>% layout(annotations = list( 
    #   list(x = 0.05 , y = 0.95, align = 'left', 
    #        text = "Figure trait_6_flashfm: \nY-axis = pval_highdlipo, \nX-axis = ps, \nSize = MPP_highdlipo_Multi, \nColor = snpGroups_flashfm", 
    #        showarrow = F, xref='paper', yref='paper')) )
    fig_6b
  })    
  
  #ch0: LD matrix----  
  # output$Plot_LD <- renderPlot({
  #   LD.plot(LD[1:20,1:20], snp.positions = GWAS[[1]]$ps[1:20])
  # })
  output$Plot_LD1 <- renderPlotly({
    #withProgress_4
    withProgress(message = 'Calculation in progress',
                 detail = '...done: LD_heatmap', value = 0, {
                   for (i in 1:10) {
                     incProgress(1/10)
                     Sys.sleep(0.05)
                   }
                 })
    
    LD = tt()$LD
    LD_selected = LD[ ,GWAS_select_t1_a()$rs]
    LD_selected = LD_selected[GWAS_select_t1_a()$rs, ]
    fig_ld <- plot_ly(x = colnames(LD_selected), y = rownames(LD_selected), z = LD_selected, colors = "Greys", type = "heatmap")
    fig_ld
  })
  
  output$Plot_LD2 <- renderPlot({
    LD = tt()$LD
    LD_selected = LD[ ,GWAS_select_t1_a()$rs]
    LD_selected = LD_selected[GWAS_select_t1_a()$rs, ]
    if (input$control_widgets_1 == 0){
      hist(LD_selected)
    }else if (input$control_widgets_1 == 3){
      hist(LD_selected)
    }else{
      LD.plot(LD_selected, snp.positions = GWAS_select_t1_a()$ps)
    }
    #withProgress_5
    withProgress(message = 'Calculation in progress',
                 detail = '...done: LD_plot', value = 0, {
                   for (i in 1:10) {
                     incProgress(1/10)
                     Sys.sleep(0.05)
                   }
                 })
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function(){paste("LD_", Sys.Date(), ".pdf", sep="")},
    content = function(file){
      pdf(file)
      LD = tt()$LD
      LD_selected = LD[ ,GWAS_select_t1_a()$rs]
      LD_selected = LD_selected[GWAS_select_t1_a()$rs, ]
      if (input$control_widgets_1 == 0){
        hist(LD_selected)
      }else if (input$control_widgets_1 == 3){
        hist(LD_selected)
      }else{
        LD.plot(LD_selected, snp.positions = GWAS_select_t1_a()$ps)
      }
      dev.off()
    })
  
  
  
  #ch0: snpGroup$group.sizes table----
  output$ch0_table0 <- DT::renderDataTable({
    snpGroups = tt()$snpGroups
    col_vector = tt()$col_vector
    #generate a vector of color based on snpGroups
    #n <- length(colnames(snpGroups$group.sizes))
    #col_vector = distinctColorPalette(n)
    DT::datatable(data.frame(rbind(colnames(snpGroups$group.sizes), snpGroups$group.sizes), row.names = c("Color", 'FINEMAP', 'flashfm')),
                  options = list(dom = 't')) %>% formatStyle(
                    colnames(snpGroups$group.sizes),
                    backgroundColor = styleEqual(colnames(snpGroups$group.sizes), col_vector
                    )
                  )
  })
  
  output$ch0_table0_ch3 <- DT::renderDataTable({
    snpGroups = tt()$snpGroups
    col_vector = tt()$col_vector
    #generate a vector of color based on snpGroups
    #n <- length(colnames(snpGroups$group.sizes))
    #col_vector = distinctColorPalette(n)
    DT::datatable(data.frame(rbind(colnames(snpGroups$group.sizes), snpGroups$group.sizes), row.names = c("Color", 'FINEMAP', 'flashfm')),
                  options = list(dom = 't')) %>% formatStyle(
                    colnames(snpGroups$group.sizes),
                    backgroundColor = styleEqual(colnames(snpGroups$group.sizes), col_vector
                    )
                  )
  })
  
  #ch1: single trait networks (group vs. snp)----
  plt_links_group_fm_selected <- reactive({    
    plt_links_group_fm = tt()$plt_links_group_fm
    plt_links_group_fm[plt_links_group_fm$l_size>=input$control_widgets_3a_fm[1]*10 & plt_links_group_fm$l_size<=input$control_widgets_3a_fm[2]*10+0.001, ]
  })
  network_ch1_network_group_fm <- reactive({forceNetwork(Links = plt_links_group_fm_selected(), 
                                                              #Nodes = plt_nodes_group_fm, 
                                                              Nodes = tt()$plt_nodes_group_fm, 
                                                              Source = 'source_id', 
                                                              Target = 'target_id', 
                                                              Value = 'l_size', 
                                                              NodeID = 'name',
                                                              Nodesize = 'n_size', 
                                                              Group = 'n_color', 
                                                              #height = 500,
                                                              #width = 500, 
                                                              fontSize = 20,
                                                              fontFamily = "serif",
                                                              linkDistance = 100,
                                                              charge = -300, #numeric value indicating either the strength of the node repulsion (negative value) or attraction (positive value).
                                                              linkColour = plt_links_group_fm_selected()$l_color,
                                                              #radiusCalculation = "d.nodesize+10",
                                                              #radiusCalculation = " Math.sqrt(d.nodesize)+10",
                                                              opacity = 0.8, 
                                                              zoom = T,
                                                              legend = T, 
                                                              bounded = T,
                                                              #colourScale = JS(ColourScale),
                                                              opacityNoHover = 1)
  })
  
  output$ch1_network_group_fm <- renderForceNetwork({
    network_ch1_network_group_fm()
  })
  
  # output$d_ch_1_s1 <- downloadHandler(
  #   filename = function() {
  #     paste("ch1_network_group_fm_", Sys.Date(), ".html", sep="")
  #   },
  #   content = function(file) {
  #     saveNetwork(network_ch1_network_group_fm(), file)
  #   }
  # )
  
  output$d_ch_1_s1 <- downloadHandler(
    filename = function() {
      paste("ch1_network_group_fm_", Sys.Date(), ".html", sep="")
    },
    content = function(file) {
      saveNetwork(network_ch1_network_group_fm(), file)
      #webshot::install_phantomjs(force = TRUE)
      #webshot(file,"testing.png", vwidth = 1000, vheight = 900)
    }
  )
  
  plt_links_snp_selected <- reactive({   
    plt_links_snp = tt()$plt_links_snp
    plt_links_snp[plt_links_snp$l_size>=input$control_widgets_3b_fm[1]*10 & plt_links_snp$l_size<=input$control_widgets_3b_fm[2]*10+0.001, ]
  })
  network_ch1_network_snp_fm <- reactive({forceNetwork(Links = plt_links_snp_selected(), 
                                                            #Nodes = plt_nodes_snp, 
                                                            Nodes = tt()$plt_nodes_snp, 
                                                            Source = 'source_id', 
                                                            Target = 'target_id', 
                                                            Value = 'l_size', 
                                                            NodeID = 'name',
                                                            Nodesize = 'n_size', 
                                                            Group = 'n_color', 
                                                            #height = 500,
                                                            #width = 500, 
                                                            fontSize = 20,
                                                            fontFamily = "serif",
                                                            linkDistance = 100,
                                                            charge = -300, #numeric value indicating either the strength of the node repulsion (negative value) or attraction (positive value).
                                                            linkColour = plt_links_snp_selected()$l_color,
                                                            #radiusCalculation = "d.nodesize+10",
                                                            #radiusCalculation = " Math.sqrt(d.nodesize)+10",
                                                            opacity = 0.8, 
                                                            zoom = T,
                                                            legend = T, 
                                                            bounded = F,
                                                            #colourScale = JS(ColourScale),
                                                            opacityNoHover = 1)
  })
  
  output$ch1_network_snp <- renderForceNetwork({
    network_ch1_network_snp_fm()
  })
  
  output$d_ch_1_s2 <- downloadHandler(
    filename = function() {
      paste("ch1_network_snp_fm_", Sys.Date(), ".html", sep="")
    },
    content = function(file) {
      saveNetwork(network_ch1_network_snp_fm(), file)
    }
  )
  
  #ch2: multi trait networks (group vs. snp)----
  plt_links_group_flashfm_selected <- reactive({
    plt_links_group_flashfm = tt()$plt_links_group_flashfm
    plt_links_group_flashfm[plt_links_group_flashfm$l_size>=input$control_widgets_3a_flashfm[1]*10 & plt_links_group_flashfm$l_size<=input$control_widgets_3a_flashfm[2]*10+0.001, ]
  })
  network_ch1_network_group_flashfm <- reactive({forceNetwork(Links = plt_links_group_flashfm_selected(),
                                                                      #Nodes = plt_nodes_group_flashfm,
                                                                      Nodes = tt()$plt_nodes_group_flashfm,
                                                                      Source = 'source_id',
                                                                      Target = 'target_id',
                                                                      Value = 'l_size',
                                                                      NodeID = 'name',
                                                                      Nodesize = 'n_size',
                                                                      Group = 'n_color',
                                                                      #height = 500,
                                                                      #width = 500,
                                                                      fontSize = 20,
                                                                      fontFamily = "serif",
                                                                      linkDistance = 100,
                                                                      charge = -300, #numeric value indicating either the strength of the node repulsion (negative value) or attraction (positive value).
                                                                      linkColour = plt_links_group_flashfm_selected()$l_color,
                                                                      #radiusCalculation = "d.nodesize+10",
                                                                      #radiusCalculation = " Math.sqrt(d.nodesize)+10",
                                                                      opacity = 0.8,
                                                                      zoom = T,
                                                                      legend = T,
                                                                      bounded = T,
                                                                      #colourScale = JS(ColourScale),
                                                                      opacityNoHover = 1)
  })
  
  output$ch1_network_group_flashfm <- renderForceNetwork({
    network_ch1_network_group_flashfm()
  })
  
  output$d_ch_2_m1 <- downloadHandler(
    filename = function() {
      paste("ch2_network_group_flashfm_", Sys.Date(), ".html", sep="")
    },
    content = function(file) {
      saveNetwork(network_ch1_network_group_flashfm(), file)
    }
  )
  
  

  plt_links_snp_selected_flashfm <- reactive({ 
    plt_links_snp_flashfm = tt()$plt_links_snp_flashfm
    plt_links_snp_flashfm[plt_links_snp_flashfm$l_size>=input$control_widgets_3b_flashfm[1]*10 & plt_links_snp_flashfm$l_size<=input$control_widgets_3b_flashfm[2]*10+0.001, ]
  })
  network_ch1_network_snp_flashfm <- reactive({forceNetwork(Links = plt_links_snp_selected_flashfm(), 
                                                                    #Nodes = plt_nodes_snp, 
                                                                    Nodes = tt()$plt_nodes_snp,
                                                                    Source = 'source_id', 
                                                                    Target = 'target_id', 
                                                                    Value = 'l_size', 
                                                                    NodeID = 'name',
                                                                    Nodesize = 'n_size', 
                                                                    Group = 'n_color', 
                                                                    #height = 500,
                                                                    #width = 500, 
                                                                    fontSize = 20,
                                                                    fontFamily = "serif",
                                                                    linkDistance = 100,
                                                                    charge = -300, #numeric value indicating either the strength of the node repulsion (negative value) or attraction (positive value).
                                                                    linkColour = plt_links_snp_selected_flashfm()$l_color,
                                                                    #radiusCalculation = "d.nodesize+10",
                                                                    #radiusCalculation = " Math.sqrt(d.nodesize)+10",
                                                                    opacity = 0.8, 
                                                                    zoom = T,
                                                                    legend = T, 
                                                                    bounded = F,
                                                                    #colourScale = JS(ColourScale),
                                                                    opacityNoHover = 1)
  })  
  
  output$ch1_network_snp_flashfm <- renderForceNetwork({
    network_ch1_network_snp_flashfm()
  })
  
  output$d_ch_2_m2 <- downloadHandler(
    filename = function() {
      paste("ch2_network_snp_flashfm_", Sys.Date(), ".html", sep="")
    },
    content = function(file) {
      saveNetwork(network_ch1_network_snp_flashfm(), file)
    }
  )
  
  #ch3: Regional association plots link selection----
  output$manhattan_ch3_linked <- renderPlotly({  
    GWAS_all_final = tt()$GWAS_all_final
    mpp.pp = tt()$mpp.pp
    pal = tt()$pal
    
    M_ch3_2 = input$control_widgets_ch_3_2
    
    shared_data = GWAS_all_final %>% highlight_key()
    baseplot = plot_ly(shared_data, height = 1200)
    
    ### New plots here----
    # Figure 1: Regional association plots
    vline <- function(x = 0, color = "gray") {
      list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x, 
           line = list(color = color, dash = 'dash'))}
    
    hline <- function(y = 0, color = "red") {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y, 
           line = list(color = color))}
    
    # Figure fig_1a: pval_log_1, size = ~MPP_lowdlipo_Single, color = ~snpGroups_fm
    fig_1a_ch3 = baseplot %>%
      add_markers(x = ~ps, y = ~pval_log_1, type = 'scatter', 
                  size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[1],"_Single"))), 
                  color = ~snpGroups_fm, 
                  #colors = pal,
                  colors = tt()$pal,
                  legendgroup = ~snpGroups_fm,
                  sizes = c(5, 30),
                  fill = ~'',
                  marker = list(opacity = 1, sizemode = 'diameter'),
                  hoverinfo = 'text',
                  text = ~paste('SNP_ID: ', rs,
                                '<br>SNP_position: ', ps,
                                '<br>Allele_1: ', allele1, ', Allele_0: ', allele0,
                                #'<br>Allele 0:', allele0,
                                '<br>Group_fm: ', snpGroups_fm, ' (', snpGroups_fm_size, ' SNPs)', 
                                '<br>MPP_trait1_Single: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[1],"_Single"))),4),
                                '<br>MPPg_trait1_Single: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[1],"_Single"))),4), 
                                '<br>(If MPPg=0: it is negligible)',
                                '<br>af: ', signif(af,4),
                                '<br>pval: ', signif(pval_1,4),
                                '<br>Credible_set_cs1: ', cs1_1, 
                                '<br>Credible_set_csM: ', csM_1,
                                sep = ""))
    
    fig_1a_ch3 <- fig_1a_ch3 %>% layout(#title = 'Interactive Regional association plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_all_final$ps), max(GWAS_all_final$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0, max(GWAS_all_final$pval_log_1)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                #showlegend = FALSE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    if (M_ch3_2 == 0){
      fig_1a_ch3 <- fig_1a_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_1==max(GWAS_all_final$pval_log_1))]), hline(7.3)))
    }else{
      fig_1a_ch3 <- fig_1a_ch3 %>% layout(shapes = list(hline(7.3)))
    }
    #fig_1a_ch3 <- fig_1a_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_1==max(GWAS_all_final$pval_log_1))]), hline(7.3)))
    #fig_1a_ch3 <- fig_1a_ch3 %>% layout(annotations = list( 
    #  list(x = 0.05 , y = 0.05, align = 'left',  text = "Figure fig_1a: pval_lowdlipo, \n size = ~MPP_lowdlipo_Single, \n color = ~snpGroups_fm", showarrow = F, xref='paper', yref='paper')) )
    
    #fig_1a_ch3 
    
    # Figure fig_1b: pval_log_1, size = ~MPP_lowdlipo_Multi, color = ~snpGroups_flashfm
    fig_1b_ch3 = baseplot %>%
      add_markers(x = ~ps, y = ~pval_log_1, type = 'scatter', 
                  size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[1],"_Multi"))), 
                  color = ~snpGroups_flashfm, 
                  #colors = pal,
                  colors = tt()$pal,
                  legendgroup = ~snpGroups_flashfm,
                  sizes = c(5, 30),
                  fill = ~'',
                  marker = list(opacity = 1, sizemode = 'diameter'),
                  hoverinfo = 'text',
                  text = ~paste('SNP_ID: ', rs,
                                '<br>SNP_position: ', ps,
                                '<br>Allele_1: ', allele1, ', Allele_0: ', allele0,
                                #'<br>Allele 0:', allele0,
                                '<br>Group_flashfm: ', snpGroups_flashfm, ' (', snpGroups_flashfm_size,' SNPs)',
                                '<br>MPP_trait1_Multi: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[1],"_Multi"))),4),
                                '<br>MPPg_trait1_Multi: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[1],"_Multi"))),4),
                                '<br>(If MPPg=0: it is negligible)',
                                '<br>af: ', signif(af,4),
                                '<br>pval: ', signif(pval_1,4),
                                '<br>Credible_set_cs1: ', cs1_1, 
                                '<br>Credible_set_csM: ', csM_1,
                                sep = ""))
    
    fig_1b_ch3 <- fig_1b_ch3 %>% layout(title = 'Regional Association Plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_all_final$ps), max(GWAS_all_final$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0, max(GWAS_all_final$pval_log_1)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                #showlegend = FALSE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    if (M_ch3_2 == 0){
      fig_1b_ch3 <- fig_1b_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_1==max(GWAS_all_final$pval_log_1))]), hline(7.3)))
    }else{
      fig_1b_ch3 <- fig_1b_ch3 %>% layout(shapes = list(hline(7.3)))
    }
    #fig_1b_ch3 <- fig_1b_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_1==max(GWAS_all_final$pval_log_1))]), hline(7.3)))
    #fig_1b_ch3 <- fig_1b_ch3 %>% layout(annotations = list( 
    #  list(x = 0.05 , y = 0.95, align = 'left',  text = "Figure fig_1b: pval_lowdlipo, \n size = ~MPP_lowdlipo_Multi, \n color = ~snpGroups_flashfm", showarrow = F, xref='paper', yref='paper')) )
    
    #fig_1b_ch3
    
    # Figure fig_2a: pval_log_2, size = ~MPP_cholesterol_Single, color = ~snpGroups_fm
    fig_2a_ch3 = baseplot %>%
      add_markers(x = ~ps, y = ~pval_log_2, type = 'scatter', 
                  size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[2],"_Single"))), 
                  color = ~snpGroups_fm, 
                  #colors = pal,
                  colors = tt()$pal,
                  legendgroup = ~snpGroups_fm,
                  sizes = c(5, 30),
                  fill = ~'',
                  marker = list(opacity = 1, sizemode = 'diameter'),
                  hoverinfo = 'text',
                  text = ~paste('SNP_ID: ', rs,
                                '<br>SNP_position: ', ps,
                                '<br>Allele_1: ', allele1, ', Allele_0: ', allele0,
                                #'<br>Allele 0:', allele0,
                                '<br>Group_fm: ', snpGroups_fm, ' (', snpGroups_fm_size,' SNPs)',
                                '<br>MPP_trait2_Single: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[2],"_Single"))),4),
                                '<br>MPPg_trait2_Single: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[2],"_Single"))),4),
                                '<br>(If MPPg=0: it is negligible)',
                                '<br>af: ', signif(af,4),
                                '<br>pval: ', signif(pval_2,4),
                                '<br>Credible_set_cs1: ', cs1_2, 
                                '<br>Credible_set_csM: ', csM_2,
                                sep = ""))
    
    fig_2a_ch3 <- fig_2a_ch3 %>% layout(title = 'Regional Association Plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_all_final$ps), max(GWAS_all_final$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0, max(GWAS_all_final$pval_log_2)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                #showlegend = FALSE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    if (M_ch3_2 == 0){
      fig_2a_ch3 <- fig_2a_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_2==max(GWAS_all_final$pval_log_2))]), hline(7.3)))
    }else{
      fig_2a_ch3 <- fig_2a_ch3 %>% layout(shapes = list(hline(7.3)))
    }
    #fig_2a_ch3 <- fig_2a_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_2==max(GWAS_all_final$pval_log_2))]), hline(7.3)))
    #fig_2a_ch3 <- fig_2a_ch3 %>% layout(annotations = list( 
    #  list(x = 0.05 , y = 0.95, align = 'left', text = "Figure fig_2a: pval_holesterol, \n size = ~MPP_holesterol_Single, \n color = ~snpGroups_fm", showarrow = F, xref='paper', yref='paper')) )
    
    #fig_2a_ch3
    
    # Figure fig_2b: pval_log_2, size = ~MPP_cholesterol_Multi, color = ~snpGroups_flashfm
    fig_2b_ch3 = baseplot %>%
      add_markers(x = ~ps, y = ~pval_log_2, type = 'scatter', 
                  size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[2],"_Multi"))), 
                  color = ~snpGroups_flashfm, 
                  #colors = pal,
                  colors = tt()$pal,
                  legendgroup = ~snpGroups_flashfm,
                  sizes = c(5, 30),
                  fill = ~'',
                  marker = list(opacity = 1, sizemode = 'diameter'),
                  hoverinfo = 'text',
                  text = ~paste('SNP_ID: ', rs,
                                '<br>SNP_position: ', ps,
                                '<br>Allele_1: ', allele1, ', Allele_0: ', allele0,
                                #'<br>Allele 0:', allele0,
                                '<br>Group_flashfm: ', snpGroups_flashfm, ' (', snpGroups_flashfm_size,' SNPs)',
                                '<br>MPP_trait2_Multi: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[2],"_Multi"))),4),
                                '<br>MPPg_trait2_Multi: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[2],"_Multi"))),4),
                                '<br>(If MPPg=0: it is negligible)',
                                '<br>af: ', signif(af,4),
                                '<br>pval: ', signif(pval_2,4),
                                '<br>Credible_set_cs1: ', cs1_2, 
                                '<br>Credible_set_csM: ', csM_2,
                                sep = ""))
    
    fig_2b_ch3 <- fig_2b_ch3 %>% layout(title = 'Regional Association Plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_all_final$ps), max(GWAS_all_final$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0, max(GWAS_all_final$pval_log_2)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                #showlegend = FALSE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    if (M_ch3_2 == 0){
      fig_2b_ch3 <- fig_2b_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_2==max(GWAS_all_final$pval_log_2))]), hline(7.3)))
    }else{
      fig_2b_ch3 <- fig_2b_ch3 %>% layout(shapes = list(hline(7.3)))
    }
    #fig_2b_ch3 <- fig_2b_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_2==max(GWAS_all_final$pval_log_2))]), hline(7.3)))
    #fig_2b_ch3 <- fig_2b_ch3 %>% layout(annotations = list( 
    #  list(x = 0.05 , y = 0.95, align = 'left',  text = "Figure fig_2b: pval_holesterol, \n size = ~MPP_holesterol_Multi, \n color = ~snpGroups_flashfm", showarrow = F, xref='paper', yref='paper')) )
    
    #fig_2b_ch3
    
    # Figure fig_3a: pval_log_3, size = ~MPP_triglycerides_Single, color = ~snpGroups_fm
    fig_3a_ch3 = baseplot %>%
      add_markers(x = ~ps, y = ~pval_log_3, type = 'scatter', 
                  size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[3],"_Single"))), 
                  color = ~snpGroups_fm, 
                  #colors = pal,
                  colors = tt()$pal,
                  legendgroup = ~snpGroups_fm,
                  sizes = c(5, 30),
                  fill = ~'',
                  marker = list(opacity = 1, sizemode = 'diameter'),
                  hoverinfo = 'text',
                  text = ~paste('SNP_ID: ', rs,
                                '<br>SNP_position: ', ps,
                                '<br>Allele_1: ', allele1, ', Allele_0: ', allele0,
                                #'<br>Allele 0:', allele0,
                                '<br>Group_fm: ', snpGroups_fm, ' (', snpGroups_fm_size,' SNPs)',
                                '<br>MPP_trait3_Single: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[3],"_Single"))),4),
                                '<br>MPPg_trait3_Single: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[3],"_Single"))),4),
                                '<br>(If MPPg=0: it is negligible)',
                                '<br>af: ', signif(af,4),
                                '<br>pval: ', signif(pval_3,4),
                                '<br>Credible_set_cs1: ', cs1_3, 
                                '<br>Credible_set_csM: ', csM_3,
                                sep = ""))
    
    fig_3a_ch3 <- fig_3a_ch3 %>% layout(title = 'Regional Association Plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_all_final$ps), max(GWAS_all_final$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0, max(GWAS_all_final$pval_log_3)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                #showlegend = FALSE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    if (M_ch3_2 == 0){
      fig_3a_ch3 <- fig_3a_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_3==max(GWAS_all_final$pval_log_3))]), hline(7.3)))
    }else{
      fig_3a_ch3 <- fig_3a_ch3 %>% layout(shapes = list(hline(7.3)))
    }
    #fig_3a_ch3 <- fig_3a_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_3==max(GWAS_all_final$pval_log_3))]), hline(7.3)))
    #fig_3a_ch3 <- fig_3a_ch3 %>% layout(annotations = list( 
    #  list(x = 0.05 , y = 0.95, align = 'left',  text = "Figure fig_3a: pval_triglycerides, \n size = ~MPP_triglycerides_Single, \n color = ~snpGroups_fm", showarrow = F, xref='paper', yref='paper')) )
    
    #fig_3a_ch3
    
    # Figure fig_3b: pval_log_3, size = ~MPP_triglycerides_Multi, color = ~snpGroups_flashfm
    fig_3b_ch3 = baseplot %>%
      add_markers(x = ~ps, y = ~pval_log_3, type = 'scatter', 
                  size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[3],"_Multi"))), 
                  color = ~snpGroups_flashfm, 
                  #colors = pal,
                  colors = tt()$pal,
                  legendgroup = ~snpGroups_flashfm,
                  sizes = c(5, 30),
                  fill = ~'',
                  marker = list(opacity = 1, sizemode = 'diameter'),
                  hoverinfo = 'text',
                  text = ~paste('SNP_ID: ', rs,
                                '<br>SNP_position: ', ps,
                                '<br>Allele_1: ', allele1, ', Allele_0: ', allele0,
                                #'<br>Allele 0:', allele0,
                                '<br>Group_flashfm: ', snpGroups_flashfm, ' (', snpGroups_flashfm_size,' SNPs)',
                                '<br>MPP_trait3_Multi: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[3],"_Multi"))),4),
                                '<br>MPPg_trait3_Multi: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[3],"_Multi"))),4),
                                '<br>(If MPPg=0: it is negligible)',
                                '<br>af: ', signif(af,4),
                                '<br>pval: ', signif(pval_3,4),
                                '<br>Credible_set_cs1: ', cs1_3, 
                                '<br>Credible_set_csM: ', csM_3,
                                sep = ""))
    
    fig_3b_ch3 <- fig_3b_ch3 %>% layout(title = 'Regional Association Plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_all_final$ps), max(GWAS_all_final$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0, max(GWAS_all_final$pval_log_3)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                #showlegend = FALSE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    if (M_ch3_2 == 0){
      fig_3b_ch3 <- fig_3b_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_3==max(GWAS_all_final$pval_log_3))]), hline(7.3)))
    }else{
      fig_3b_ch3 <- fig_3b_ch3 %>% layout(shapes = list(hline(7.3)))
    }
    #fig_3b_ch3 <- fig_3b_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_3==max(GWAS_all_final$pval_log_3))]), hline(7.3)))
    #fig_3b_ch3 <- fig_3b_ch3 %>% layout(annotations = list( 
    #  list(x = 0.05 , y = 0.95, align = 'left',  text = "Figure fig_3b: pval_triglycerides, \n size = ~MPP_triglycerides_Multi, \n color = ~snpGroups_flashfm", showarrow = F, xref='paper', yref='paper')) )
    
    #fig_3b_ch3
    
    # Figure fig_4a: pval_log_4, size = ~MPP_highdlipo_Single, color = ~snpGroups_fm
    fig_4a_ch3 = baseplot %>%
      add_markers(x = ~ps, y = ~pval_log_4, type = 'scatter', 
                  size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[4],"_Single"))), 
                  color = ~snpGroups_fm, 
                  #colors = pal,
                  colors = tt()$pal,
                  legendgroup = ~snpGroups_fm,
                  sizes = c(5, 30),
                  fill = ~'',
                  marker = list(opacity = 1, sizemode = 'diameter'),
                  hoverinfo = 'text',
                  text = ~paste('SNP_ID: ', rs,
                                '<br>SNP_position: ', ps,
                                '<br>Allele_1: ', allele1, ', Allele_0: ', allele0,
                                #'<br>Allele 0:', allele0,
                                '<br>Group_fm: ', snpGroups_fm, ' (', snpGroups_fm_size,' SNPs)',
                                '<br>MPP_trait4_Single: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[4],"_Single"))),4),
                                '<br>MPPg_trait4_Single: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[4],"_Single"))),4),
                                '<br>(If MPPg=0: it is negligible)',
                                '<br>af: ', signif(af,4),
                                '<br>pval: ', signif(pval_4,4),
                                '<br>Credible_set_cs1: ', cs1_4, 
                                '<br>Credible_set_csM: ', csM_4,
                                sep = ""))
    
    fig_4a_ch3 <- fig_4a_ch3 %>% layout(title = 'Regional Association Plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_all_final$ps), max(GWAS_all_final$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0, max(GWAS_all_final$pval_log_4)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                #showlegend = FALSE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    if (M_ch3_2 == 0){
      fig_4a_ch3 <- fig_4a_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_4==max(GWAS_all_final$pval_log_4))]), hline(7.3)))
    }else{
      fig_4a_ch3 <- fig_4a_ch3 %>% layout(shapes = list(hline(7.3)))
    }
    #fig_4a_ch3 <- fig_4a_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_4==max(GWAS_all_final$pval_log_4))]), hline(7.3)))
    #fig_4a_ch3 <- fig_4a_ch3 %>% layout(annotations = list( 
    #  list(x = 0.05 , y = 0.95, align = 'left',  text = "Figure fig_4a: pval_highdlipo, \n size = ~MPP_highdlipo_Single, \n color = ~snpGroups_fm", showarrow = F, xref='paper', yref='paper')) )
    
    #fig_4a_ch3
    
    # Figure fig_4b: pval_log_4, size = ~MPP_highdlipo_Multi, color = ~snpGroups_flashfm
    fig_4b_ch3 = baseplot %>%
      add_markers(x = ~ps, y = ~pval_log_4, type = 'scatter', 
                  size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[4],"_Multi"))), 
                  color = ~snpGroups_flashfm, 
                  #colors = pal,
                  colors = tt()$pal,
                  legendgroup = ~snpGroups_flashfm,
                  sizes = c(5, 30),
                  fill = ~'',
                  marker = list(opacity = 1, sizemode = 'diameter'),
                  hoverinfo = 'text',
                  text = ~paste('SNP_ID: ', rs,
                                '<br>SNP_position: ', ps,
                                '<br>Allele_1: ', allele1, ', Allele_0: ', allele0,
                                #'<br>Allele 0:', allele0,
                                '<br>Group_flashfm: ', snpGroups_flashfm, ' (', snpGroups_flashfm_size,' SNPs)',
                                '<br>MPP_trait4_Multi: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[4],"_Multi"))),4),
                                '<br>MPPg_trait4_Multi: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[4],"_Multi"))),4),
                                '<br>(If MPPg=0: it is negligible)',
                                '<br>af: ', signif(af,4),
                                '<br>pval: ', signif(pval_4,4),
                                '<br>Credible_set_cs1: ', cs1_4, 
                                '<br>Credible_set_csM: ', csM_4,
                                sep = ""))
    
    fig_4b_ch3 <- fig_4b_ch3 %>% layout(title = 'Regional Association Plots',
                                xaxis = list(#title = 'Chromosome 19 position',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(min(GWAS_all_final$ps), max(GWAS_all_final$ps)),
                                             #type = 'log',
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2,
                                             showgrid = FALSE),
                                yaxis = list(title = '-log10(p)',
                                             gridcolor = 'rgb(255, 255, 255)',
                                             range = c(0, max(GWAS_all_final$pval_log_4)+20),
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2,
                                             showgrid = FALSE),
                                #showlegend = FALSE,
                                paper_bgcolor = 'rgb(243, 243, 243)',
                                plot_bgcolor = 'rgb(243, 243, 243)')
    if (M_ch3_2 == 0){
      fig_4b_ch3 <- fig_4b_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_4==max(GWAS_all_final$pval_log_4))]), hline(7.3)))
    }else{
      fig_4b_ch3 <- fig_4b_ch3 %>% layout(shapes = list(hline(7.3)))
    }
    #fig_4b_ch3 <- fig_4b_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_4==max(GWAS_all_final$pval_log_4))]), hline(7.3)))
    #fig_4b_ch3 <- fig_4b_ch3 %>% layout(annotations = list( 
    #  list(x = 0.05 , y = 0.95, align = 'left',  text = "Figure fig_4b: pval_highdlipo, \n size = ~MPP_highdlipo_Multi, \n color = ~snpGroups_flashfm", showarrow = F, xref='paper', yref='paper')) )
    
    #fig_4b_ch3
    
    # Figure fig_5a: pval_log_4, size = ~MPP_highdlipo_Single, color = ~snpGroups_fm
    fig_5a_ch3 = baseplot %>%
      add_markers(x = ~ps, y = ~pval_log_5, type = 'scatter', 
                  size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[5],"_Single"))), 
                  color = ~snpGroups_fm, 
                  #colors = pal,
                  colors = tt()$pal,
                  legendgroup = ~snpGroups_fm,
                  sizes = c(5, 30),
                  fill = ~'',
                  marker = list(opacity = 1, sizemode = 'diameter'),
                  hoverinfo = 'text',
                  text = ~paste('SNP_ID: ', rs,
                                '<br>SNP_position: ', ps,
                                '<br>Allele_1: ', allele1, ', Allele_0: ', allele0,
                                #'<br>Allele 0:', allele0,
                                '<br>Group_fm: ', snpGroups_fm, ' (', snpGroups_fm_size, ' SNPs)',
                                '<br>MPP_trait5_Single: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[5],"_Single"))),4),
                                '<br>MPPg_trait5_Single: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[5],"_Single"))),4),
                                '<br>(If MPPg=0: it is negligible)',
                                '<br>af: ', signif(af,4),
                                '<br>pval: ', signif(pval_5,4),
                                '<br>Credible_set_cs1: ', cs1_5, 
                                '<br>Credible_set_csM: ', csM_5,
                                sep = ""))
    
    fig_5a_ch3 <- fig_5a_ch3 %>% layout(title = 'Regional Association Plots',
                                        xaxis = list(#title = 'Chromosome 19 position',
                                                     gridcolor = 'rgb(255, 255, 255)',
                                                     range = c(min(GWAS_all_final$ps), max(GWAS_all_final$ps)),
                                                     #type = 'log',
                                                     zerolinewidth = 1,
                                                     ticklen = 5,
                                                     gridwidth = 2,
                                                     showgrid = FALSE),
                                        yaxis = list(title = '-log10(p)',
                                                     gridcolor = 'rgb(255, 255, 255)',
                                                     range = c(0, max(GWAS_all_final$pval_log_5)+20),
                                                     zerolinewidth = 1,
                                                     ticklen = 5,
                                                     gridwith = 2,
                                                     showgrid = FALSE),
                                        #showlegend = FALSE,
                                        paper_bgcolor = 'rgb(243, 243, 243)',
                                        plot_bgcolor = 'rgb(243, 243, 243)')
    if (M_ch3_2 == 0){
      fig_5a_ch3 <- fig_5a_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_5==max(GWAS_all_final$pval_log_5))]), hline(7.3)))
    }else{
      fig_5a_ch3 <- fig_5a_ch3 %>% layout(shapes = list(hline(7.3)))
    }
    #fig_5a_ch3 <- fig_5a_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_5==max(GWAS_all_final$pval_log_5))]), hline(7.3)))
    #fig_5a_ch3 <- fig_5a_ch3 %>% layout(annotations = list( 
    #  list(x = 0.05 , y = 0.95, align = 'left',  text = "Figure fig_5a: pval_highdlipo, \n size = ~MPP_highdlipo_Single, \n color = ~snpGroups_fm", showarrow = F, xref='paper', yref='paper')) )
    
    #fig_5a_ch3
    
    # Figure fig_5b: pval_log_4, size = ~MPP_highdlipo_Multi, color = ~snpGroups_flashfm
    fig_5b_ch3 = baseplot %>%
      add_markers(x = ~ps, y = ~pval_log_5, type = 'scatter', 
                  size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[5],"_Multi"))), 
                  color = ~snpGroups_flashfm, 
                  #colors = pal,
                  colors = tt()$pal,
                  legendgroup = ~snpGroups_flashfm,
                  sizes = c(5, 30),
                  fill = ~'',
                  marker = list(opacity = 1, sizemode = 'diameter'),
                  hoverinfo = 'text',
                  text = ~paste('SNP_ID: ', rs,
                                '<br>SNP_position: ', ps,
                                '<br>Allele_1: ', allele1, ', Allele_0: ', allele0,
                                #'<br>Allele 0:', allele0,
                                '<br>Group_flashfm: ', snpGroups_flashfm, ' (', snpGroups_flashfm_size,' SNPs)',
                                '<br>MPP_trait5_Multi: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[5],"_Multi"))),4),
                                '<br>MPPg_trait5_Multi: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[5],"_Multi"))),4),
                                '<br>(If MPPg=0: it is negligible)',
                                '<br>af: ', signif(af,4),
                                '<br>pval: ', signif(pval_5,4),
                                '<br>Credible_set_cs1: ', cs1_5, 
                                '<br>Credible_set_csM: ', csM_5,
                                sep = ""))
    
    fig_5b_ch3 <- fig_5b_ch3 %>% layout(title = 'Regional Association Plots',
                                        xaxis = list(#title = 'Chromosome 19 position',
                                                     gridcolor = 'rgb(255, 255, 255)',
                                                     range = c(min(GWAS_all_final$ps), max(GWAS_all_final$ps)),
                                                     #type = 'log',
                                                     zerolinewidth = 1,
                                                     ticklen = 5,
                                                     gridwidth = 2,
                                                     showgrid = FALSE),
                                        yaxis = list(title = '-log10(p)',
                                                     gridcolor = 'rgb(255, 255, 255)',
                                                     range = c(0, max(GWAS_all_final$pval_log_5)+20),
                                                     zerolinewidth = 1,
                                                     ticklen = 5,
                                                     gridwith = 2,
                                                     showgrid = FALSE),
                                        #showlegend = FALSE,
                                        paper_bgcolor = 'rgb(243, 243, 243)',
                                        plot_bgcolor = 'rgb(243, 243, 243)')
    if (M_ch3_2 == 0){
      fig_5b_ch3 <- fig_5b_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_5==max(GWAS_all_final$pval_log_5))]), hline(7.3)))
    }else{
      fig_5b_ch3 <- fig_5b_ch3 %>% layout(shapes = list(hline(7.3)))
    }
    #fig_5b_ch3 <- fig_5b_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_5==max(GWAS_all_final$pval_log_5))]), hline(7.3)))
    #fig_5b_ch3 <- fig_5b_ch3 %>% layout(annotations = list( 
    #  list(x = 0.05 , y = 0.95, align = 'left',  text = "Figure fig_5b: pval_highdlipo, \n size = ~MPP_highdlipo_Multi, \n color = ~snpGroups_flashfm", showarrow = F, xref='paper', yref='paper')) )
    
    #fig_5b_ch3
    
    # Figure fig_6a: pval_log_4, size = ~MPP_highdlipo_Single, color = ~snpGroups_fm
    fig_6a_ch3 = baseplot %>%
      add_markers(x = ~ps, y = ~pval_log_6, type = 'scatter', 
                  size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[6],"_Single"))), 
                  color = ~snpGroups_fm, 
                  #colors = pal,
                  colors = tt()$pal,
                  legendgroup = ~snpGroups_fm,
                  sizes = c(5, 30),
                  fill = ~'',
                  marker = list(opacity = 1, sizemode = 'diameter'),
                  hoverinfo = 'text',
                  text = ~paste('SNP_ID: ', rs,
                                '<br>SNP_position: ', ps,
                                '<br>Allele_1: ', allele1, ', Allele_0: ', allele0,
                                #'<br>Allele 0:', allele0,
                                '<br>Group_fm: ', snpGroups_fm, ' (', snpGroups_fm_size, ' SNPs)',
                                '<br>MPP_trait6_Single: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[6],"_Single"))),4),
                                '<br>MPPg_trait6_Single: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[6],"_Single"))),4),
                                '<br>(If MPPg=0: it is negligible)',
                                '<br>af: ', signif(af,4),
                                '<br>pval: ', signif(pval_6,4),
                                '<br>Credible_set_cs1: ', cs1_6, 
                                '<br>Credible_set_csM: ', csM_6,
                                sep = ""))
    
    fig_6a_ch3 <- fig_6a_ch3 %>% layout(title = 'Regional Association Plots',
                                        xaxis = list(#title = 'Chromosome 19 position',
                                                     gridcolor = 'rgb(255, 255, 255)',
                                                     range = c(min(GWAS_all_final$ps), max(GWAS_all_final$ps)),
                                                     #type = 'log',
                                                     zerolinewidth = 1,
                                                     ticklen = 5,
                                                     gridwidth = 2,
                                                     showgrid = FALSE),
                                        yaxis = list(title = '-log10(p)',
                                                     gridcolor = 'rgb(255, 255, 255)',
                                                     range = c(0, max(GWAS_all_final$pval_log_6)+20),
                                                     zerolinewidth = 1,
                                                     ticklen = 5,
                                                     gridwith = 2,
                                                     showgrid = FALSE),
                                        #showlegend = FALSE,
                                        paper_bgcolor = 'rgb(243, 243, 243)',
                                        plot_bgcolor = 'rgb(243, 243, 243)')
    if (M_ch3_2 == 0){
      fig_6a_ch3 <- fig_6a_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_6==max(GWAS_all_final$pval_log_6))]), hline(7.3)))
    }else{
      fig_6a_ch3 <- fig_6a_ch3 %>% layout(shapes = list(hline(7.3)))
    }
    #fig_6a_ch3 <- fig_6a_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_6==max(GWAS_all_final$pval_log_6))]), hline(7.3)))
    #fig_6a_ch3 <- fig_6a_ch3 %>% layout(annotations = list( 
    #  list(x = 0.05 , y = 0.95, align = 'left',  text = "Figure fig_6a: pval_highdlipo, \n size = ~MPP_highdlipo_Single, \n color = ~snpGroups_fm", showarrow = F, xref='paper', yref='paper')) )
    
    #fig_6a_ch3
    
    # Figure fig_6b: pval_log_4, size = ~MPP_highdlipo_Multi, color = ~snpGroups_flashfm
    fig_6b_ch3 = baseplot %>%
      add_markers(x = ~ps, y = ~pval_log_6, type = 'scatter', 
                  size = ~eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[6],"_Multi"))), 
                  color = ~snpGroups_flashfm, 
                  #colors = pal,
                  colors = tt()$pal,
                  legendgroup = ~snpGroups_flashfm,
                  sizes = c(5, 30),
                  fill = ~'',
                  marker = list(opacity = 1, sizemode = 'diameter'),
                  hoverinfo = 'text',
                  text = ~paste('SNP_ID: ', rs,
                                '<br>SNP_position: ', ps,
                                '<br>Allele_1: ', allele1, ', Allele_0: ', allele0,
                                #'<br>Allele 0:', allele0,
                                '<br>Group_flashfm: ', snpGroups_flashfm, ' (', snpGroups_flashfm_size, ' SNPs)',
                                '<br>MPP_trait6_Multi: ', signif(eval(parse(text = paste0("MPP_",names(mpp.pp$MPP)[6],"_Multi"))),4),
                                '<br>MPPg_trait6_Multi: ', signif(eval(parse(text = paste0("MPPg_",names(mpp.pp$MPPg)[6],"_Multi"))),4),
                                '<br>(If MPPg=0: it is negligible)',
                                '<br>af: ', signif(af,4),
                                '<br>pval: ', signif(pval_6,4),
                                '<br>Credible_set_cs1: ', cs1_6, 
                                '<br>Credible_set_csM: ', csM_6,
                                sep = ""))
    
    fig_6b_ch3 <- fig_6b_ch3 %>% layout(title = 'Regional Association Plots',
                                        xaxis = list(#title = 'Chromosome 19 position',
                                                     gridcolor = 'rgb(255, 255, 255)',
                                                     range = c(min(GWAS_all_final$ps), max(GWAS_all_final$ps)),
                                                     #type = 'log',
                                                     zerolinewidth = 1,
                                                     ticklen = 5,
                                                     gridwidth = 2,
                                                     showgrid = FALSE),
                                        yaxis = list(title = '-log10(p)',
                                                     gridcolor = 'rgb(255, 255, 255)',
                                                     range = c(0, max(GWAS_all_final$pval_log_6)+20),
                                                     zerolinewidth = 1,
                                                     ticklen = 5,
                                                     gridwith = 2,
                                                     showgrid = FALSE),
                                        #showlegend = FALSE,
                                        paper_bgcolor = 'rgb(243, 243, 243)',
                                        plot_bgcolor = 'rgb(243, 243, 243)')
    if (M_ch3_2 == 0){
      fig_6b_ch3 <- fig_6b_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_6==max(GWAS_all_final$pval_log_6))]), hline(7.3)))
    }else{
      fig_6b_ch3 <- fig_6b_ch3 %>% layout(shapes = list(hline(7.3)))
    }
    #fig_6b_ch3 <- fig_6b_ch3 %>% layout(shapes = list(vline(GWAS_all_final$ps[which(GWAS_all_final$pval_log_6==max(GWAS_all_final$pval_log_6))]), hline(7.3)))
    #fig_6b_ch3 <- fig_6b_ch3 %>% layout(annotations = list( 
    #  list(x = 0.05 , y = 0.95, align = 'left',  text = "Figure fig_6b: pval_highdlipo, \n size = ~MPP_highdlipo_Multi, \n color = ~snpGroups_flashfm", showarrow = F, xref='paper', yref='paper')) )
    
    #fig_6b_ch3
    
    #selection
    M_ch3=input$control_widgets_ch_3  ##CHECK
    
    if (M_ch3 == 2){
      fig_all_ch3_total <- subplot(fig_1a_ch3,
                                   style(fig_1b_ch3, showlegend = F), 
                                   style(fig_2a_ch3, showlegend = F),
                                   style(fig_2b_ch3, showlegend = F),
                                   nrows = 2) %>% 
        layout(dragmode = "select",
               plot_bgcolor='#e5ecf6', 
               xaxis = list(zerolinecolor = '#ffff', 
                            zerolinewidth = 2, 
                            gridcolor = 'ffff'), 
               yaxis = list(zerolinecolor = '#ffff', 
                            zerolinewidth = 2, 
                            gridcolor = 'ffff')) %>% 
        layout(autosize = T 
               #width = 1000, 
               #height = 1000
               ) %>% 
        # highlight(on = "plotly_click")
        highlight(on = "plotly_selected",
                  off = "plotly_deselect")
    }else if (M_ch3 == 3){
      fig_all_ch3_total <- subplot(fig_1a_ch3, 
                                   style(fig_1b_ch3, showlegend = F), 
                                   style(fig_2a_ch3, showlegend = F),
                                   style(fig_2b_ch3, showlegend = F), 
                                   style(fig_3a_ch3, showlegend = F), 
                                   style(fig_3b_ch3, showlegend = F), 
                                   nrows = 3) %>% 
        layout(dragmode = "select",
               plot_bgcolor='#e5ecf6', 
               xaxis = list(zerolinecolor = '#ffff', 
                            zerolinewidth = 2, 
                            gridcolor = 'ffff'), 
               yaxis = list(zerolinecolor = '#ffff', 
                            zerolinewidth = 2, 
                            gridcolor = 'ffff')) %>% 
        layout(autosize = T
               #width = 1000, 
               #height = 1000
               )%>% 
        # highlight(on = "plotly_click")
        highlight(on = "plotly_selected",
                  off = "plotly_deselect")
    } else if (M_ch3 == 4){
      fig_all_ch3_total <- subplot(fig_1a_ch3, 
                                   style(fig_1b_ch3, showlegend = F), 
                                   style(fig_2a_ch3, showlegend = F),
                                   style(fig_2b_ch3, showlegend = F), 
                                   style(fig_3a_ch3, showlegend = F), 
                                   style(fig_3b_ch3, showlegend = F),  
                                   style(fig_4a_ch3, showlegend = F), 
                                   style(fig_4b_ch3, showlegend = F),  
                                   nrows = 4) %>% 
        layout(dragmode = "select",
               plot_bgcolor='#e5ecf6', 
               xaxis = list(zerolinecolor = '#ffff', 
                            zerolinewidth = 2, 
                            gridcolor = 'ffff'), 
               yaxis = list(zerolinecolor = '#ffff', 
                            zerolinewidth = 2, 
                            gridcolor = 'ffff')) %>% 
        layout(autosize = T 
               #width = 1000, 
               #height = 1000
               ) %>% 
        # highlight(on = "plotly_click")
        highlight(on = "plotly_selected",
                  off = "plotly_deselect")
    }else if (M_ch3 == 5){
      fig_all_ch3_total <- subplot(fig_1a_ch3, 
                                   style(fig_1b_ch3, showlegend = F), 
                                   style(fig_2a_ch3, showlegend = F),
                                   style(fig_2b_ch3, showlegend = F), 
                                   style(fig_3a_ch3, showlegend = F), 
                                   style(fig_3b_ch3, showlegend = F),  
                                   style(fig_4a_ch3, showlegend = F), 
                                   style(fig_4b_ch3, showlegend = F), 
                                   style(fig_5a_ch3, showlegend = F), 
                                   style(fig_5b_ch3, showlegend = F),  
                                   nrows = 5) %>% 
        layout(dragmode = "select",
               plot_bgcolor='#e5ecf6', 
               xaxis = list(zerolinecolor = '#ffff', 
                            zerolinewidth = 2, 
                            gridcolor = 'ffff'), 
               yaxis = list(zerolinecolor = '#ffff', 
                            zerolinewidth = 2, 
                            gridcolor = 'ffff')) %>% 
        layout(autosize = T 
               #width = 1000, 
               #height = 1000
               ) %>% 
        # highlight(on = "plotly_click")
        highlight(on = "plotly_selected",
                  off = "plotly_deselect")
    }
    else if (M_ch3 == 6){
      fig_all_ch3_total <- subplot(fig_1a_ch3, 
                                   style(fig_1b_ch3, showlegend = F), 
                                   style(fig_2a_ch3, showlegend = F),
                                   style(fig_2b_ch3, showlegend = F), 
                                   style(fig_3a_ch3, showlegend = F), 
                                   style(fig_3b_ch3, showlegend = F),  
                                   style(fig_4a_ch3, showlegend = F), 
                                   style(fig_4b_ch3, showlegend = F), 
                                   style(fig_5a_ch3, showlegend = F), 
                                   style(fig_5b_ch3, showlegend = F), 
                                   style(fig_6a_ch3, showlegend = F), 
                                   style(fig_6b_ch3, showlegend = F), 
                                   nrows = 6) %>% 
        layout(dragmode = "select",
               plot_bgcolor='#e5ecf6', 
               xaxis = list(zerolinecolor = '#ffff', 
                            zerolinewidth = 2, 
                            gridcolor = 'ffff'), 
               yaxis = list(zerolinecolor = '#ffff', 
                            zerolinewidth = 2, 
                            gridcolor = 'ffff')) %>% 
        layout(autosize = T 
               #width = 1000, 
               #height = 1000
               ) %>% 
        # highlight(on = "plotly_click")
        highlight(on = "plotly_selected",
                  off = "plotly_deselect")
    }
    
    #withProgress_5
    withProgress(message = 'Calculation in progress',
                 detail = '...done: Linked_Manhattan', value = 0, {
                   for (i in 1:10) {
                     incProgress(1/10)
                     Sys.sleep(0.05)
                   }
                 })
    
    fig_all_ch3_total
  })
  
  #ch4: ----
  #Radar chart----
  output$fig_ch4_tab1_radar <- renderPlotly({ 
    
    GWAS_all_final = tt()$GWAS_all_final
    mpp.pp = tt()$mpp.pp
    
    fig_ch4_radar <- plot_ly(type = 'scatterpolar',
                             fill = 'toself',
                             mode = 'markers') 
    
    #selection
    M_ch4=input$control_widgets_ch_4  ##CHECK
    
    if (M_ch4 == 3){
      fig_ch4_radar <- fig_ch4_radar %>%
        add_trace(r = c(sum(GWAS_all_final$cs1_1), 
                        sum(GWAS_all_final$cs1_2),
                        sum(GWAS_all_final$cs1_3),
                        sum(GWAS_all_final$cs1_1)),
                  # theta = c(paste0('Trait_1_',names(mpp.pp$MPP)[1]),
                  #           paste0('Trait_2_',names(mpp.pp$MPP)[2]),
                  #           paste0('Trait_3_',names(mpp.pp$MPP)[3]), 
                  #           paste0('Trait_1_',names(mpp.pp$MPP)[1])),
                  theta = c(names(mpp.pp$MPP)[1],
                            names(mpp.pp$MPP)[2],
                            names(mpp.pp$MPP)[3], 
                            names(mpp.pp$MPP)[1]),
                  name = 'cs_1') 
      fig_ch4_radar <- fig_ch4_radar %>%
        add_trace(r = c(sum(GWAS_all_final$csM_1), 
                        sum(GWAS_all_final$csM_2),
                        sum(GWAS_all_final$csM_3),
                        sum(GWAS_all_final$csM_1)),
                  # theta = c(paste0('Trait_1_',names(mpp.pp$MPP)[1]),
                  #           paste0('Trait_2_',names(mpp.pp$MPP)[2]),
                  #           paste0('Trait_3_',names(mpp.pp$MPP)[3]), 
                  #           paste0('Trait_1_',names(mpp.pp$MPP)[1])),
                  theta = c(names(mpp.pp$MPP)[1],
                            names(mpp.pp$MPP)[2],
                            names(mpp.pp$MPP)[3], 
                            names(mpp.pp$MPP)[1]),
                  name = 'cs_M')
    }else if (M_ch4 == 4){
      fig_ch4_radar <- fig_ch4_radar %>%
        add_trace(r = c(sum(GWAS_all_final$cs1_1), 
                        sum(GWAS_all_final$cs1_2),
                        sum(GWAS_all_final$cs1_3),
                        sum(GWAS_all_final$cs1_4),
                        sum(GWAS_all_final$cs1_1)),
                  # theta = c(paste0('Trait_1_',names(mpp.pp$MPP)[1]),
                  #           paste0('Trait_2_',names(mpp.pp$MPP)[2]),
                  #           paste0('Trait_3_',names(mpp.pp$MPP)[3]), 
                  #           paste0('Trait_4_',names(mpp.pp$MPP)[4]), 
                  #           paste0('Trait_1_',names(mpp.pp$MPP)[1])),
                  theta = c(names(mpp.pp$MPP)[1],
                            names(mpp.pp$MPP)[2],
                            names(mpp.pp$MPP)[3], 
                            names(mpp.pp$MPP)[4], 
                            names(mpp.pp$MPP)[1]),
                  name = 'cs_1') 
      fig_ch4_radar <- fig_ch4_radar %>%
        add_trace(r = c(sum(GWAS_all_final$csM_1), 
                        sum(GWAS_all_final$csM_2),
                        sum(GWAS_all_final$csM_3),
                        sum(GWAS_all_final$csM_4),
                        sum(GWAS_all_final$csM_1)),
                  # theta = c(paste0('Trait_1_',names(mpp.pp$MPP)[1]),
                  #           paste0('Trait_2_',names(mpp.pp$MPP)[2]),
                  #           paste0('Trait_3_',names(mpp.pp$MPP)[3]), 
                  #           paste0('Trait_4_',names(mpp.pp$MPP)[4]), 
                  #           paste0('Trait_1_',names(mpp.pp$MPP)[1])),
                  theta = c(names(mpp.pp$MPP)[1],
                            names(mpp.pp$MPP)[2],
                            names(mpp.pp$MPP)[3], 
                            names(mpp.pp$MPP)[4], 
                            names(mpp.pp$MPP)[1]),
                  name = 'cs_M')
    }else if (M_ch4 == 5){
      fig_ch4_radar <- fig_ch4_radar %>%
        add_trace(r = c(sum(GWAS_all_final$cs1_1), 
                        sum(GWAS_all_final$cs1_2),
                        sum(GWAS_all_final$cs1_3),
                        sum(GWAS_all_final$cs1_4),
                        sum(GWAS_all_final$cs1_5),
                        sum(GWAS_all_final$cs1_1)),
                  # theta = c(paste0('Trait_1_',names(mpp.pp$MPP)[1]),
                  #           paste0('Trait_2_',names(mpp.pp$MPP)[2]),
                  #           paste0('Trait_3_',names(mpp.pp$MPP)[3]), 
                  #           paste0('Trait_4_',names(mpp.pp$MPP)[4]),
                  #           paste0('Trait_5_',names(mpp.pp$MPP)[5]), 
                  #           paste0('Trait_1_',names(mpp.pp$MPP)[1])),
                  theta = c(names(mpp.pp$MPP)[1],
                            names(mpp.pp$MPP)[2],
                            names(mpp.pp$MPP)[3], 
                            names(mpp.pp$MPP)[4],
                            names(mpp.pp$MPP)[5], 
                            names(mpp.pp$MPP)[1]),
                  name = 'cs_1') 
      fig_ch4_radar <- fig_ch4_radar %>%
        add_trace(r = c(sum(GWAS_all_final$csM_1), 
                        sum(GWAS_all_final$csM_2),
                        sum(GWAS_all_final$csM_3),
                        sum(GWAS_all_final$csM_4),
                        sum(GWAS_all_final$csM_5),
                        sum(GWAS_all_final$csM_1)),
                  # theta = c(paste0('Trait_1_',names(mpp.pp$MPP)[1]),
                  #           paste0('Trait_2_',names(mpp.pp$MPP)[2]),
                  #           paste0('Trait_3_',names(mpp.pp$MPP)[3]), 
                  #           paste0('Trait_4_',names(mpp.pp$MPP)[4]), 
                  #           paste0('Trait_5_',names(mpp.pp$MPP)[5]), 
                  #           paste0('Trait_1_',names(mpp.pp$MPP)[1])),
                  theta = c(names(mpp.pp$MPP)[1],
                            names(mpp.pp$MPP)[2],
                            names(mpp.pp$MPP)[3], 
                            names(mpp.pp$MPP)[4], 
                            names(mpp.pp$MPP)[5], 
                            names(mpp.pp$MPP)[1]),
                  name = 'cs_M')
    }else if (M_ch4 == 6){
    fig_ch4_radar <- fig_ch4_radar %>%
                     add_trace(r = c(sum(GWAS_all_final$cs1_1), 
                                     sum(GWAS_all_final$cs1_2),
                                     sum(GWAS_all_final$cs1_3),
                                     sum(GWAS_all_final$cs1_4),
                                     sum(GWAS_all_final$cs1_5),
                                     sum(GWAS_all_final$cs1_6),
                                     sum(GWAS_all_final$cs1_1)),
                                     # theta = c(paste0('Trait_1_',names(mpp.pp$MPP)[1]),
                                     #           paste0('Trait_2_',names(mpp.pp$MPP)[2]),
                                     #           paste0('Trait_3_',names(mpp.pp$MPP)[3]), 
                                     #           paste0('Trait_4_',names(mpp.pp$MPP)[4]), 
                                     #           paste0('Trait_5_',names(mpp.pp$MPP)[5]), 
                                     #           paste0('Trait_6_',names(mpp.pp$MPP)[6]), 
                                     #           paste0('Trait_1_',names(mpp.pp$MPP)[1])),
                               theta = c(names(mpp.pp$MPP)[1],
                                         names(mpp.pp$MPP)[2],
                                         names(mpp.pp$MPP)[3], 
                                         names(mpp.pp$MPP)[4], 
                                         names(mpp.pp$MPP)[5], 
                                         names(mpp.pp$MPP)[6], 
                                         names(mpp.pp$MPP)[1]),
                                     name = 'cs_1') 
    fig_ch4_radar <- fig_ch4_radar %>%
                     add_trace(r = c(sum(GWAS_all_final$csM_1), 
                                     sum(GWAS_all_final$csM_2),
                                     sum(GWAS_all_final$csM_3),
                                     sum(GWAS_all_final$csM_4),
                                     sum(GWAS_all_final$csM_5),
                                     sum(GWAS_all_final$csM_6),
                                     sum(GWAS_all_final$csM_1)),
                                     # theta = c(paste0('Trait_1_',names(mpp.pp$MPP)[1]),
                                     #           paste0('Trait_2_',names(mpp.pp$MPP)[2]),
                                     #           paste0('Trait_3_',names(mpp.pp$MPP)[3]), 
                                     #           paste0('Trait_4_',names(mpp.pp$MPP)[4]),
                                     #           paste0('Trait_5_',names(mpp.pp$MPP)[5]), 
                                     #           paste0('Trait_6_',names(mpp.pp$MPP)[6]),
                                     #           paste0('Trait_1_',names(mpp.pp$MPP)[1])),
                               theta = c(names(mpp.pp$MPP)[1],
                                         names(mpp.pp$MPP)[2],
                                         names(mpp.pp$MPP)[3], 
                                         names(mpp.pp$MPP)[4],
                                         names(mpp.pp$MPP)[5], 
                                         names(mpp.pp$MPP)[6],
                                         names(mpp.pp$MPP)[1]),
                                     name = 'cs_M') 
    }
    fig_ch4_radar <- fig_ch4_radar %>%
                     layout(polar = list(radialaxis = list(visible = T, range = c(0,max(c(sum(GWAS_all_final$csM_1), 
                                                                                          sum(GWAS_all_final$csM_2),
                                                                                          sum(GWAS_all_final$csM_3),
                                                                                          sum(GWAS_all_final$csM_4),
                                                                                          sum(GWAS_all_final$csM_5),
                                                                                          sum(GWAS_all_final$csM_6),
                                                                                          sum(GWAS_all_final$cs1_1), 
                                                                                          sum(GWAS_all_final$cs1_2),
                                                                                          sum(GWAS_all_final$cs1_3),
                                                                                          sum(GWAS_all_final$cs1_4),
                                                                                          sum(GWAS_all_final$cs1_5),
                                                                                          sum(GWAS_all_final$cs1_6)))))))
    fig_ch4_radar
  })
  #Venn diagram cs_1----
  fig_ch4_tab2_venn_cs1_data <- reactive({
    GWAS_all_final = tt()$GWAS_all_final
    mpp.pp = tt()$mpp.pp
    
    M_ch4=input$control_widgets_ch_4  ##CHECK
    x_cs1_1 <- x_cs1_2 <- x_cs1_3 <- x_cs1_4 <- x_cs1_5 <- x_cs1_6 <- c();
    if (M_ch4 == 3){
      for (i in 1:length(GWAS_all_final$rs)){
        if (GWAS_all_final$cs1_1[i]==1){ x_cs1_1 = c(x_cs1_1, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$cs1_2[i]==1){ x_cs1_2 = c(x_cs1_2, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$cs1_3[i]==1){ x_cs1_3 = c(x_cs1_3, GWAS_all_final$rs[i]) }
      }
      x_cs1_all <- list(cs1_1_a = x_cs1_1, 
                        cs1_2_a = x_cs1_2, 
                        cs1_3_a = x_cs1_3)
      # names(x_cs1_all) <- c(paste0('cs1_trait1_',names(mpp.pp$MPP)[1]),
      #                       paste0('cs1_trait2_',names(mpp.pp$MPP)[2]),
      #                       paste0('cs1_trait3_',names(mpp.pp$MPP)[3]))
      names(x_cs1_all) <- c(names(mpp.pp$MPP)[1],
                            names(mpp.pp$MPP)[2],
                            names(mpp.pp$MPP)[3])
      x_cs1_all
    }else if (M_ch4 == 4){
      for (i in 1:length(GWAS_all_final$rs)){
        if (GWAS_all_final$cs1_1[i]==1){ x_cs1_1 = c(x_cs1_1, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$cs1_2[i]==1){ x_cs1_2 = c(x_cs1_2, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$cs1_3[i]==1){ x_cs1_3 = c(x_cs1_3, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$cs1_4[i]==1){ x_cs1_4 = c(x_cs1_4, GWAS_all_final$rs[i]) }
      }
      x_cs1_all <- list(cs1_1 = x_cs1_1,  
                        cs1_2 = x_cs1_2, 
                        cs1_3 = x_cs1_3, 
                        cs1_4 = x_cs1_4)
      # names(x_cs1_all) <- c(paste0('cs1_trait1_',names(mpp.pp$MPP)[1]),
      #                       paste0('cs1_trait2_',names(mpp.pp$MPP)[2]),
      #                       paste0('cs1_trait3_',names(mpp.pp$MPP)[3]),
      #                       paste0('cs1_trait4_',names(mpp.pp$MPP)[4]))
      names(x_cs1_all) <- c(names(mpp.pp$MPP)[1],
                            names(mpp.pp$MPP)[2],
                            names(mpp.pp$MPP)[3],
                            names(mpp.pp$MPP)[4])
      x_cs1_all
    }else if (M_ch4 == 5){
      for (i in 1:length(GWAS_all_final$rs)){
        if (GWAS_all_final$cs1_1[i]==1){ x_cs1_1 = c(x_cs1_1, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$cs1_2[i]==1){ x_cs1_2 = c(x_cs1_2, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$cs1_3[i]==1){ x_cs1_3 = c(x_cs1_3, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$cs1_4[i]==1){ x_cs1_4 = c(x_cs1_4, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$cs1_5[i]==1){ x_cs1_5 = c(x_cs1_5, GWAS_all_final$rs[i]) }
      }
      x_cs1_all <- list(cs1_1 = x_cs1_1,  
                        cs1_2 = x_cs1_2, 
                        cs1_3 = x_cs1_3, 
                        cs1_4 = x_cs1_4, 
                        cs1_5 = x_cs1_5)
      # names(x_cs1_all) <- c(paste0('cs1_trait1_',names(mpp.pp$MPP)[1]),
      #                       paste0('cs1_trait2_',names(mpp.pp$MPP)[2]),
      #                       paste0('cs1_trait3_',names(mpp.pp$MPP)[3]),
      #                       paste0('cs1_trait4_',names(mpp.pp$MPP)[4]),
      #                       paste0('cs1_trait5_',names(mpp.pp$MPP)[5]))
      names(x_cs1_all) <- c(names(mpp.pp$MPP)[1],
                            names(mpp.pp$MPP)[2],
                            names(mpp.pp$MPP)[3],
                            names(mpp.pp$MPP)[4],
                            names(mpp.pp$MPP)[5])
      x_cs1_all
    }else if (M_ch4 == 6){
      for (i in 1:length(GWAS_all_final$rs)){
        if (GWAS_all_final$cs1_1[i]==1){ x_cs1_1 = c(x_cs1_1, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$cs1_2[i]==1){ x_cs1_2 = c(x_cs1_2, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$cs1_3[i]==1){ x_cs1_3 = c(x_cs1_3, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$cs1_4[i]==1){ x_cs1_4 = c(x_cs1_4, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$cs1_5[i]==1){ x_cs1_5 = c(x_cs1_5, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$cs1_6[i]==1){ x_cs1_6 = c(x_cs1_6, GWAS_all_final$rs[i]) }
      }
      x_cs1_all <- list(cs1_1 = x_cs1_1,
                        cs1_2 = x_cs1_2,
                        cs1_3 = x_cs1_3, 
                        cs1_4 = x_cs1_4, 
                        cs1_5 = x_cs1_5, 
                        cs1_6 = x_cs1_6)
      # names(x_cs1_all) <- c(paste0('cs1_trait1_',names(mpp.pp$MPP)[1]),
      #                       paste0('cs1_trait2_',names(mpp.pp$MPP)[2]),
      #                       paste0('cs1_trait3_',names(mpp.pp$MPP)[3]),
      #                       paste0('cs1_trait4_',names(mpp.pp$MPP)[4]),
      #                       paste0('cs1_trait5_',names(mpp.pp$MPP)[5]),
      #                       paste0('cs1_trait6_',names(mpp.pp$MPP)[6]))
      names(x_cs1_all) <- c(names(mpp.pp$MPP)[1],
                            names(mpp.pp$MPP)[2],
                            names(mpp.pp$MPP)[3],
                            names(mpp.pp$MPP)[4],
                            names(mpp.pp$MPP)[5],
                            names(mpp.pp$MPP)[6])
      x_cs1_all
    }
  })
  
  output$fig_ch4_tab2_venn_cs1 <- renderPlotly({     
    ggVennDiagram(fig_ch4_tab2_venn_cs1_data(),
                  show_intersect = TRUE)
    # ggplotly(ggVennDiagram(fig_ch4_tab2_venn_cs1_data()) + 
    #            scale_color_brewer(palette = "Paired"))
  })
  
  output$table_ch4_tab2_venn_cs1_table <- DT::renderDataTable({
    check_x_cs1_table = process_region_data(Venn(fig_ch4_tab2_venn_cs1_data()))[,c(1,2,5,4,3)]
    DT::datatable(data.frame(check_x_cs1_table), 
                  extensions = 'Buttons',
                  options = list(
                    paging = FALSE,
                    searching = TRUE,
                    fixedColumns = TRUE,
                    autoWidth = TRUE,
                    ordering = TRUE,
                    dom = 'tB',
                    buttons = list( 
                      list(extend = 'csv',   filename =  "cs1_table"),
                      list(extend = 'excel', filename =  "cs1_table"))
                  ),
                  class = "display"
                  )
  })
  
  #Venn diagram cs_M----
  fig_ch4_tab2_venn_csM_data <- reactive({
    GWAS_all_final = tt()$GWAS_all_final
    mpp.pp = tt()$mpp.pp
    
    M_ch4=input$control_widgets_ch_4  ##CHECK
    x_csM_1 <- x_csM_2 <- x_csM_3 <- x_csM_4 <- x_csM_5 <- x_csM_6 <- c();
    if (M_ch4 == 3){
      for (i in 1:length(GWAS_all_final$rs)){
        if (GWAS_all_final$csM_1[i]==1){ x_csM_1 = c(x_csM_1, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$csM_2[i]==1){ x_csM_2 = c(x_csM_2, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$csM_3[i]==1){ x_csM_3 = c(x_csM_3, GWAS_all_final$rs[i]) }
      }
      x_csM_all <- list(csM_1_a = x_csM_1, 
                        csM_2_a = x_csM_2, 
                        csM_3_a = x_csM_3)
      # names(x_csM_all) <- c(paste0('csM_trait1_',names(mpp.pp$MPP)[1]),
      #                       paste0('csM_trait2_',names(mpp.pp$MPP)[2]),
      #                       paste0('csM_trait3_',names(mpp.pp$MPP)[3]))
      names(x_csM_all) <- c(names(mpp.pp$MPP)[1],
                            names(mpp.pp$MPP)[2],
                            names(mpp.pp$MPP)[3])
      x_csM_all
    }else if (M_ch4 == 4){
      for (i in 1:length(GWAS_all_final$rs)){
        if (GWAS_all_final$csM_1[i]==1){ x_csM_1 = c(x_csM_1, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$csM_2[i]==1){ x_csM_2 = c(x_csM_2, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$csM_3[i]==1){ x_csM_3 = c(x_csM_3, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$csM_4[i]==1){ x_csM_4 = c(x_csM_4, GWAS_all_final$rs[i]) }
      }
      x_csM_all <- list(csM_1 = x_csM_1,  
                        csM_2 = x_csM_2, 
                        csM_3 = x_csM_3, 
                        csM_4 = x_csM_4)
      # names(x_csM_all) <- c(paste0('csM_trait1_',names(mpp.pp$MPP)[1]),
      #                       paste0('csM_trait2_',names(mpp.pp$MPP)[2]),
      #                       paste0('csM_trait3_',names(mpp.pp$MPP)[3]),
      #                       paste0('csM_trait4_',names(mpp.pp$MPP)[4]))
      names(x_csM_all) <- c(names(mpp.pp$MPP)[1],
                            names(mpp.pp$MPP)[2],
                            names(mpp.pp$MPP)[3],
                            names(mpp.pp$MPP)[4])
      x_csM_all
    }else if (M_ch4 == 5){
      for (i in 1:length(GWAS_all_final$rs)){
        if (GWAS_all_final$csM_1[i]==1){ x_csM_1 = c(x_csM_1, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$csM_2[i]==1){ x_csM_2 = c(x_csM_2, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$csM_3[i]==1){ x_csM_3 = c(x_csM_3, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$csM_4[i]==1){ x_csM_4 = c(x_csM_4, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$csM_5[i]==1){ x_csM_5 = c(x_csM_5, GWAS_all_final$rs[i]) }
      }
      x_csM_all <- list(csM_1 = x_csM_1,  
                        csM_2 = x_csM_2, 
                        csM_3 = x_csM_3, 
                        csM_4 = x_csM_4, 
                        csM_5 = x_csM_5)
      # names(x_csM_all) <- c(paste0('csM_trait1_',names(mpp.pp$MPP)[1]),
      #                       paste0('csM_trait2_',names(mpp.pp$MPP)[2]),
      #                       paste0('csM_trait3_',names(mpp.pp$MPP)[3]),
      #                       paste0('csM_trait4_',names(mpp.pp$MPP)[4]),
      #                       paste0('csM_trait5_',names(mpp.pp$MPP)[5]))
      names(x_csM_all) <- c(names(mpp.pp$MPP)[1],
                            names(mpp.pp$MPP)[2],
                            names(mpp.pp$MPP)[3],
                            names(mpp.pp$MPP)[4],
                            names(mpp.pp$MPP)[5])
      x_csM_all
    }else if (M_ch4 == 6){
      for (i in 1:length(GWAS_all_final$rs)){
        if (GWAS_all_final$csM_1[i]==1){ x_csM_1 = c(x_csM_1, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$csM_2[i]==1){ x_csM_2 = c(x_csM_2, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$csM_3[i]==1){ x_csM_3 = c(x_csM_3, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$csM_4[i]==1){ x_csM_4 = c(x_csM_4, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$csM_5[i]==1){ x_csM_5 = c(x_csM_5, GWAS_all_final$rs[i]) }
        if (GWAS_all_final$csM_6[i]==1){ x_csM_6 = c(x_csM_6, GWAS_all_final$rs[i]) }
      }
      x_csM_all <- list(csM_1 = x_csM_1,
                        csM_2 = x_csM_2,
                        csM_3 = x_csM_3, 
                        csM_4 = x_csM_4, 
                        csM_5 = x_csM_5, 
                        csM_6 = x_csM_6)
      # names(x_csM_all) <- c(paste0('csM_trait1_',names(mpp.pp$MPP)[1]),
      #                       paste0('csM_trait2_',names(mpp.pp$MPP)[2]),
      #                       paste0('csM_trait3_',names(mpp.pp$MPP)[3]),
      #                       paste0('csM_trait4_',names(mpp.pp$MPP)[4]),
      #                       paste0('csM_trait5_',names(mpp.pp$MPP)[5]),
      #                       paste0('csM_trait6_',names(mpp.pp$MPP)[6]))
      names(x_csM_all) <- c(names(mpp.pp$MPP)[1],
                            names(mpp.pp$MPP)[2],
                            names(mpp.pp$MPP)[3],
                            names(mpp.pp$MPP)[4],
                            names(mpp.pp$MPP)[5],
                            names(mpp.pp$MPP)[6])
      x_csM_all
    }
  })
  
  output$fig_ch4_tab2_venn_csM <- renderPlotly({     
    ggVennDiagram(fig_ch4_tab2_venn_csM_data(), show_intersect = TRUE)
  })
  
  output$table_ch4_tab2_venn_csM_table <- DT::renderDataTable({
    check_x_csM_table = process_region_data(Venn(fig_ch4_tab2_venn_csM_data()))[,c(1,2,5,4,3)]
    DT::datatable(data.frame(check_x_csM_table), 
                  extensions = 'Buttons',
                  options = list(
                    paging = FALSE,
                    searching = TRUE,
                    fixedColumns = TRUE,
                    autoWidth = TRUE,
                    scrollX=TRUE,
                    ordering = TRUE,
                    dom = 'tB',
                    buttons = list( 
                      list(extend = 'csv',   filename =  "csM_table"),
                      list(extend = 'excel', filename =  "csM_table"))
                  ),
                  class = "display"
                  )
  })
  
  
  
  
  
  #ch5: ----
  #Sankey Diagram in R and plotly (all 4 traits)----
  output$ch5_sankey_all <- renderPlotly({ 
    GWAS_all_final = tt()$GWAS_all_final
    mpp.pp = tt()$mpp.pp
    
    node_labels_snp             = GWAS_all_final$rs
    node_labels_groupfm         = unique(GWAS_all_final$snpGroups_fm)
    node_labels_groupfm_2       = paste("fm", node_labels_groupfm, sep="_")
    node_labels_groupflashfm    = unique(GWAS_all_final$snpGroups_flashfm)
    node_labels_groupflashfm_2  = paste("flashfm", node_labels_groupflashfm, sep="_")
    
    node_labels   = c(node_labels_snp, node_labels_groupfm, node_labels_groupflashfm)
    node_labels_2 = c(node_labels_snp, node_labels_groupfm_2, node_labels_groupflashfm_2)
    
    node_x      = c(rep(0.01,times=length(node_labels_snp)), 
                    rep(0.4,times=length(node_labels_groupfm)),
                    rep(0.8,times=length(node_labels_groupflashfm)))
    node_y      = c(seq(from = 1, to = length(node_labels_snp), by = 1)/length(node_labels_snp), 
                    seq(from = 1, to = length(node_labels_groupfm), by = 1)/length(node_labels_groupfm), 
                    seq(from = 1, to = length(node_labels_groupflashfm), by = 1)/length(node_labels_groupflashfm))
    
    node_table  = as.data.frame(cbind(node_labels, node_x, node_y));
    
    link_source_fm = c();
    link_target_fm = c();
    link_value_fm  = c();
    link_source_flashfm = c();
    link_target_flashfm = c();
    link_value_flashfm  = c();
    #group_fm
    for (i in 1:length(node_labels_snp)){
      if (GWAS_all_final$snpGroups_fm[i] != '0'){
        link_source_fm = c(link_source_fm, i-1)
        link_value_fm  = c(link_value_fm, 
                           mean(c(GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[1],"_Single")][i],
                                  GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[2],"_Single")][i],
                                  GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[3],"_Single")][i],
                                  GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[4],"_Single")][i]),
                                na.rm=F) #updated20220208
        )
        for (ii in 1:length(node_labels_groupfm)){
          if (node_table$node_labels[length(node_labels_snp)+ii]==GWAS_all_final$snpGroups_fm[i]){
            link_target_fm = c(link_target_fm, length(node_labels_snp)+ii-1)
          }
        }
      }
    }
    #group_flashfm
    for (i in 1:length(node_labels_snp)){
      if (GWAS_all_final$snpGroups_flashfm[i] != '0'){
        link_source_flashfm = c(link_source_flashfm, i-1)
        link_value_flashfm  = c(link_value_flashfm, 
                                mean(c(GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[1],"_Multi")][i],
                                       GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[2],"_Multi")][i],
                                       GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[3],"_Multi")][i],
                                       GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[4],"_Multi")][i]),
                                     na.rm=F)  #updated20220208
        )
        for (ii in 1:length(node_labels_groupflashfm)){
          if (node_table$node_labels[length(c(node_labels_snp,node_labels_groupfm))+ii]==GWAS_all_final$snpGroups_flashfm[i]){
            link_target_flashfm = c(link_target_flashfm, length(c(node_labels_snp,node_labels_groupfm))+ii-1)
          }
        }
      }
    }
    #check direction carefully: fm->flashfm
    # link_source = c(link_source_fm, link_source_flashfm)
    # link_target = c(link_target_fm, link_target_flashfm)
    # link_value  = c(link_value_fm, link_value_flashfm)
    #check direction carefully: flashfm->fm
    link_source = c(link_source_flashfm, link_source_fm)
    link_target = c(link_target_flashfm, link_target_fm)
    link_value  = c(link_value_flashfm, link_value_fm)
    
    link_source_updated = link_source
    for (i in 2:length(link_source)){
      if (sum(link_source[i]==link_source[1:(i-1)])!=0){
        link_source_updated[i] = link_target[which(link_source[i]==link_source[1:(i-1)])]
      }
    }
    
    fig_ch5_sankey_all <- plot_ly(type = "sankey", arrangement = "perpendicular",  #arrangement = "snap", 
                                  domain = list(x =  c(0,1),y =  c(0,1)),
                                  orientation = "h",
                                  valueformat = "0.0000",
                                  #valuesuffix = "TWh",
                          node = list(label = node_labels_2,
                                      pad = 15,
                                      thickness = 15,
                                      line = list(color = "black", width = 0.5)),
                          link = list(source = link_source_updated,      #c(0, 0, 1, 2, 5, 4, 3, 5),
                                      target = link_target,      #c(5, 3, 4, 3, 0, 2, 2, 3),
                                      value =  link_value))       #c(1, 2, 1, 1, 1, 1, 1, 2)))
    fig_ch5_sankey_all <- fig_ch5_sankey_all %>% layout(title = "Sankey Diagram of SNP groups (all traits)")
    
    fig_ch5_sankey_all
    
  })
  
  #Sankey Diagram in R and plotly (trait_1)----
  output$ch5_sankey_t1 <- renderPlotly({ 
    GWAS_all_final = tt()$GWAS_all_final
    mpp.pp = tt()$mpp.pp
    
    node_labels_snp             = GWAS_all_final$rs
    node_labels_groupfm         = unique(GWAS_all_final$snpGroups_fm)
    node_labels_groupfm_2       = paste("fm", node_labels_groupfm, sep="_")
    node_labels_groupflashfm    = unique(GWAS_all_final$snpGroups_flashfm)
    node_labels_groupflashfm_2  = paste("flashfm", node_labels_groupflashfm, sep="_")
    
    node_labels   = c(node_labels_snp, node_labels_groupfm, node_labels_groupflashfm)
    node_labels_2 = c(node_labels_snp, node_labels_groupfm_2, node_labels_groupflashfm_2)
    
    node_x      = c(rep(0.01,times=length(node_labels_snp)), 
                    rep(0.4,times=length(node_labels_groupfm)),
                    rep(0.8,times=length(node_labels_groupflashfm)))
    node_y      = c(seq(from = 1, to = length(node_labels_snp), by = 1)/length(node_labels_snp), 
                    seq(from = 1, to = length(node_labels_groupfm), by = 1)/length(node_labels_groupfm), 
                    seq(from = 1, to = length(node_labels_groupflashfm), by = 1)/length(node_labels_groupflashfm))
    
    node_table  = as.data.frame(cbind(node_labels, node_x, node_y));
    
    link_source_fm = c();
    link_target_fm = c();
    link_value_fm  = c();
    link_source_flashfm = c();
    link_target_flashfm = c();
    link_value_flashfm  = c();
    #group_fm
    for (i in 1:length(node_labels_snp)){
      if (GWAS_all_final$snpGroups_fm[i] != '0'){
        link_source_fm = c(link_source_fm, i-1)
        link_value_fm  = c(link_value_fm, 
                           mean(c(GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[1],"_Single")][i]), na.rm=F) #updated20220208
        )
        for (ii in 1:length(node_labels_groupfm)){
          if (node_table$node_labels[length(node_labels_snp)+ii]==GWAS_all_final$snpGroups_fm[i]){
            link_target_fm = c(link_target_fm, length(node_labels_snp)+ii-1)
          }
        }
      }
    }
    #group_flashfm
    for (i in 1:length(node_labels_snp)){
      if (GWAS_all_final$snpGroups_flashfm[i] != '0'){
        link_source_flashfm = c(link_source_flashfm, i-1)
        link_value_flashfm  = c(link_value_flashfm, 
                                mean(c(GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[1],"_Multi")][i]), na.rm=F) #updated20220208
        )
        for (ii in 1:length(node_labels_groupflashfm)){
          if (node_table$node_labels[length(c(node_labels_snp,node_labels_groupfm))+ii]==GWAS_all_final$snpGroups_flashfm[i]){
            link_target_flashfm = c(link_target_flashfm, length(c(node_labels_snp,node_labels_groupfm))+ii-1)
          }
        }
      }
    }
    #check direction carefully: fm->flashfm
    # link_source = c(link_source_fm, link_source_flashfm)
    # link_target = c(link_target_fm, link_target_flashfm)
    # link_value  = c(link_value_fm, link_value_flashfm)
    #check direction carefully: flashfm->fm
    link_source = c(link_source_flashfm, link_source_fm)
    link_target = c(link_target_flashfm, link_target_fm)
    link_value  = c(link_value_flashfm, link_value_fm)
    
    link_source_updated = link_source
    for (i in 2:length(link_source)){
      if (sum(link_source[i]==link_source[1:(i-1)])!=0){
        link_source_updated[i] = link_target[which(link_source[i]==link_source[1:(i-1)])]
      }
    }
    
    fig_ch5_sankey_t1 <- plot_ly(type = "sankey", arrangement = "snap", 
                                  domain = list(x =  c(0,1),y =  c(0,1)),
                                  orientation = "h",
                                 valueformat = "0.0000",
                                 #valuesuffix = "TWh",
                                  node = list(label = node_labels_2,
                                              pad = 15,
                                              thickness = 15,
                                              line = list(color = "black", width = 0.5)),
                                  link = list(source = link_source_updated,      #c(0, 0, 1, 2, 5, 4, 3, 5),
                                              target = link_target,      #c(5, 3, 4, 3, 0, 2, 2, 3),
                                              value =  link_value))       #c(1, 2, 1, 1, 1, 1, 1, 2)))
    fig_ch5_sankey_t1 <- fig_ch5_sankey_t1 %>% layout(title = paste0("Sankey Diagram of SNP groups (trait 1): ", names(mpp.pp$MPP)[1]))
    
    fig_ch5_sankey_t1
    
  })

  #Sankey Diagram in R and plotly (trait_2)----
  output$ch5_sankey_t2 <- renderPlotly({ 
    GWAS_all_final = tt()$GWAS_all_final
    mpp.pp = tt()$mpp.pp
    
    node_labels_snp             = GWAS_all_final$rs
    node_labels_groupfm         = unique(GWAS_all_final$snpGroups_fm)
    node_labels_groupfm_2       = paste("fm", node_labels_groupfm, sep="_")
    node_labels_groupflashfm    = unique(GWAS_all_final$snpGroups_flashfm)
    node_labels_groupflashfm_2  = paste("flashfm", node_labels_groupflashfm, sep="_")
    
    node_labels   = c(node_labels_snp, node_labels_groupfm, node_labels_groupflashfm)
    node_labels_2 = c(node_labels_snp, node_labels_groupfm_2, node_labels_groupflashfm_2)
    
    node_x      = c(rep(0.01,times=length(node_labels_snp)), 
                    rep(0.4,times=length(node_labels_groupfm)),
                    rep(0.8,times=length(node_labels_groupflashfm)))
    node_y      = c(seq(from = 1, to = length(node_labels_snp), by = 1)/length(node_labels_snp), 
                    seq(from = 1, to = length(node_labels_groupfm), by = 1)/length(node_labels_groupfm), 
                    seq(from = 1, to = length(node_labels_groupflashfm), by = 1)/length(node_labels_groupflashfm))
    
    node_table  = as.data.frame(cbind(node_labels, node_x, node_y));
    
    link_source_fm = c();
    link_target_fm = c();
    link_value_fm  = c();
    link_source_flashfm = c();
    link_target_flashfm = c();
    link_value_flashfm  = c();
    #group_fm
    for (i in 1:length(node_labels_snp)){
      if (GWAS_all_final$snpGroups_fm[i] != '0'){
        link_source_fm = c(link_source_fm, i-1)
        link_value_fm  = c(link_value_fm, 
                           mean(c(GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[2],"_Single")][i]),
                           na.rm=F) #updated20220208
        )
        for (ii in 1:length(node_labels_groupfm)){
          if (node_table$node_labels[length(node_labels_snp)+ii]==GWAS_all_final$snpGroups_fm[i]){
            link_target_fm = c(link_target_fm, length(node_labels_snp)+ii-1)
          }
        }
      }
    }
    #group_flashfm
    for (i in 1:length(node_labels_snp)){
      if (GWAS_all_final$snpGroups_flashfm[i] != '0'){
        link_source_flashfm = c(link_source_flashfm, i-1)
        link_value_flashfm  = c(link_value_flashfm, 
                                mean(c(GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[2],"_Multi")][i]),
                                na.rm=F) #updated20220208
        )
        for (ii in 1:length(node_labels_groupflashfm)){
          if (node_table$node_labels[length(c(node_labels_snp,node_labels_groupfm))+ii]==GWAS_all_final$snpGroups_flashfm[i]){
            link_target_flashfm = c(link_target_flashfm, length(c(node_labels_snp,node_labels_groupfm))+ii-1)
          }
        }
      }
    }
    #check direction carefully: fm->flashfm
    # link_source = c(link_source_fm, link_source_flashfm)
    # link_target = c(link_target_fm, link_target_flashfm)
    # link_value  = c(link_value_fm, link_value_flashfm)
    #check direction carefully: flashfm->fm
    link_source = c(link_source_flashfm, link_source_fm)
    link_target = c(link_target_flashfm, link_target_fm)
    link_value  = c(link_value_flashfm, link_value_fm)
    
    link_source_updated = link_source
    for (i in 2:length(link_source)){
      if (sum(link_source[i]==link_source[1:(i-1)])!=0){
        link_source_updated[i] = link_target[which(link_source[i]==link_source[1:(i-1)])]
      }
    }
    
    fig_ch5_sankey_t2 <- plot_ly(type = "sankey", arrangement = "snap", 
                                 domain = list(x =  c(0,1),y =  c(0,1)),
                                 orientation = "h",
                                 valueformat = "0.0000",
                                 #valuesuffix = "TWh",
                                 node = list(label = node_labels_2,
                                             pad = 15,
                                             thickness = 15,
                                             line = list(color = "black", width = 0.5)),
                                 link = list(source = link_source_updated,      #c(0, 0, 1, 2, 5, 4, 3, 5),
                                             target = link_target,      #c(5, 3, 4, 3, 0, 2, 2, 3),
                                             value =  link_value))       #c(1, 2, 1, 1, 1, 1, 1, 2)))
    fig_ch5_sankey_t2 <- fig_ch5_sankey_t2 %>% layout(title = paste0("Sankey Diagram of SNP groups (trait 2): ", names(mpp.pp$MPP)[2]))
    
    fig_ch5_sankey_t2
    
  })  
  
  #Sankey Diagram in R and plotly (trait_3)----
  output$ch5_sankey_t3 <- renderPlotly({ 
    GWAS_all_final = tt()$GWAS_all_final
    mpp.pp = tt()$mpp.pp
    
    node_labels_snp             = GWAS_all_final$rs
    node_labels_groupfm         = unique(GWAS_all_final$snpGroups_fm)
    node_labels_groupfm_2       = paste("fm", node_labels_groupfm, sep="_")
    node_labels_groupflashfm    = unique(GWAS_all_final$snpGroups_flashfm)
    node_labels_groupflashfm_2  = paste("flashfm", node_labels_groupflashfm, sep="_")
    
    node_labels   = c(node_labels_snp, node_labels_groupfm, node_labels_groupflashfm)
    node_labels_2 = c(node_labels_snp, node_labels_groupfm_2, node_labels_groupflashfm_2)
    
    node_x      = c(rep(0.01,times=length(node_labels_snp)), 
                    rep(0.4,times=length(node_labels_groupfm)),
                    rep(0.8,times=length(node_labels_groupflashfm)))
    node_y      = c(seq(from = 1, to = length(node_labels_snp), by = 1)/length(node_labels_snp), 
                    seq(from = 1, to = length(node_labels_groupfm), by = 1)/length(node_labels_groupfm), 
                    seq(from = 1, to = length(node_labels_groupflashfm), by = 1)/length(node_labels_groupflashfm))
    
    node_table  = as.data.frame(cbind(node_labels, node_x, node_y));
    
    link_source_fm = c();
    link_target_fm = c();
    link_value_fm  = c();
    link_source_flashfm = c();
    link_target_flashfm = c();
    link_value_flashfm  = c();
    #group_fm
    for (i in 1:length(node_labels_snp)){
      if (GWAS_all_final$snpGroups_fm[i] != '0'){
        link_source_fm = c(link_source_fm, i-1)
        link_value_fm  = c(link_value_fm, 
                           mean(c(GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[3],"_Single")][i]),
                           na.rm=F) #updated20220208
        )
        for (ii in 1:length(node_labels_groupfm)){
          if (node_table$node_labels[length(node_labels_snp)+ii]==GWAS_all_final$snpGroups_fm[i]){
            link_target_fm = c(link_target_fm, length(node_labels_snp)+ii-1)
          }
        }
      }
    }
    #group_flashfm
    for (i in 1:length(node_labels_snp)){
      if (GWAS_all_final$snpGroups_flashfm[i] != '0'){
        link_source_flashfm = c(link_source_flashfm, i-1)
        link_value_flashfm  = c(link_value_flashfm, 
                                mean(c(GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[3],"_Multi")][i]),
                                na.rm=F) #updated20220208
        )
        for (ii in 1:length(node_labels_groupflashfm)){
          if (node_table$node_labels[length(c(node_labels_snp,node_labels_groupfm))+ii]==GWAS_all_final$snpGroups_flashfm[i]){
            link_target_flashfm = c(link_target_flashfm, length(c(node_labels_snp,node_labels_groupfm))+ii-1)
          }
        }
      }
    }
    #check direction carefully: fm->flashfm
    # link_source = c(link_source_fm, link_source_flashfm)
    # link_target = c(link_target_fm, link_target_flashfm)
    # link_value  = c(link_value_fm, link_value_flashfm)
    #check direction carefully: flashfm->fm
    link_source = c(link_source_flashfm, link_source_fm)
    link_target = c(link_target_flashfm, link_target_fm)
    link_value  = c(link_value_flashfm, link_value_fm)
    
    link_source_updated = link_source
    for (i in 2:length(link_source)){
      if (sum(link_source[i]==link_source[1:(i-1)])!=0){
        link_source_updated[i] = link_target[which(link_source[i]==link_source[1:(i-1)])]
      }
    }
    
    fig_ch5_sankey_t3 <- plot_ly(type = "sankey", arrangement = "snap", 
                                 domain = list(x =  c(0,1),y =  c(0,1)),
                                 orientation = "h",
                                 valueformat = "0.0000",
                                 #valuesuffix = "TWh",
                                 node = list(label = node_labels_2,
                                             pad = 15,
                                             thickness = 15,
                                             line = list(color = "black", width = 0.5)),
                                 link = list(source = link_source_updated,      #c(0, 0, 1, 2, 5, 4, 3, 5),
                                             target = link_target,      #c(5, 3, 4, 3, 0, 2, 2, 3),
                                             value =  link_value))       #c(1, 2, 1, 1, 1, 1, 1, 2)))
    fig_ch5_sankey_t3 <- fig_ch5_sankey_t3 %>% layout(title = paste0("Sankey Diagram of SNP groups (trait 3): ", names(mpp.pp$MPP)[3]))
    
    fig_ch5_sankey_t3
    
  })
  
  #Sankey Diagram in R and plotly (trait_4)----
  output$ch5_sankey_t4 <- renderPlotly({ 
    GWAS_all_final = tt()$GWAS_all_final
    mpp.pp = tt()$mpp.pp
    
    node_labels_snp             = GWAS_all_final$rs
    node_labels_groupfm         = unique(GWAS_all_final$snpGroups_fm)
    node_labels_groupfm_2       = paste("fm", node_labels_groupfm, sep="_")
    node_labels_groupflashfm    = unique(GWAS_all_final$snpGroups_flashfm)
    node_labels_groupflashfm_2  = paste("flashfm", node_labels_groupflashfm, sep="_")
    
    node_labels   = c(node_labels_snp, node_labels_groupfm, node_labels_groupflashfm)
    node_labels_2 = c(node_labels_snp, node_labels_groupfm_2, node_labels_groupflashfm_2)
    
    node_x      = c(rep(0.01,times=length(node_labels_snp)), 
                    rep(0.4,times=length(node_labels_groupfm)),
                    rep(0.8,times=length(node_labels_groupflashfm)))
    node_y      = c(seq(from = 1, to = length(node_labels_snp), by = 1)/length(node_labels_snp), 
                    seq(from = 1, to = length(node_labels_groupfm), by = 1)/length(node_labels_groupfm), 
                    seq(from = 1, to = length(node_labels_groupflashfm), by = 1)/length(node_labels_groupflashfm))
    
    node_table  = as.data.frame(cbind(node_labels, node_x, node_y));
    
    link_source_fm = c();
    link_target_fm = c();
    link_value_fm  = c();
    link_source_flashfm = c();
    link_target_flashfm = c();
    link_value_flashfm  = c();
    #group_fm
    for (i in 1:length(node_labels_snp)){
      if (GWAS_all_final$snpGroups_fm[i] != '0'){
        link_source_fm = c(link_source_fm, i-1)
        link_value_fm  = c(link_value_fm, 
                           mean(c(GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[4],"_Single")][i]),
                           na.rm=F) #updated20220208
        )
        for (ii in 1:length(node_labels_groupfm)){
          if (node_table$node_labels[length(node_labels_snp)+ii]==GWAS_all_final$snpGroups_fm[i]){
            link_target_fm = c(link_target_fm, length(node_labels_snp)+ii-1)
          }
        }
      }
    }
    #group_flashfm
    for (i in 1:length(node_labels_snp)){
      if (GWAS_all_final$snpGroups_flashfm[i] != '0'){
        link_source_flashfm = c(link_source_flashfm, i-1)
        link_value_flashfm  = c(link_value_flashfm, 
                                mean(c(GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[4],"_Multi")][i]),
                                na.rm=F) #updated20220208
        )
        for (ii in 1:length(node_labels_groupflashfm)){
          if (node_table$node_labels[length(c(node_labels_snp,node_labels_groupfm))+ii]==GWAS_all_final$snpGroups_flashfm[i]){
            link_target_flashfm = c(link_target_flashfm, length(c(node_labels_snp,node_labels_groupfm))+ii-1)
          }
        }
      }
    }
    #check direction carefully: fm->flashfm
    # link_source = c(link_source_fm, link_source_flashfm)
    # link_target = c(link_target_fm, link_target_flashfm)
    # link_value  = c(link_value_fm, link_value_flashfm)
    #check direction carefully: flashfm->fm
    link_source = c(link_source_flashfm, link_source_fm)
    link_target = c(link_target_flashfm, link_target_fm)
    link_value  = c(link_value_flashfm, link_value_fm)
    
    link_source_updated = link_source
    for (i in 2:length(link_source)){
      if (sum(link_source[i]==link_source[1:(i-1)])!=0){
        link_source_updated[i] = link_target[which(link_source[i]==link_source[1:(i-1)])]
      }
    }
    
    fig_ch5_sankey_t4 <- plot_ly(type = "sankey", arrangement = "snap", 
                                 domain = list(x =  c(0,1),y =  c(0,1)),
                                 orientation = "h",
                                 valueformat = "0.0000",
                                 #valuesuffix = "TWh",
                                 node = list(label = node_labels_2,
                                             pad = 15,
                                             thickness = 15,
                                             line = list(color = "black", width = 0.5)),
                                 link = list(source = link_source_updated,      #c(0, 0, 1, 2, 5, 4, 3, 5),
                                             target = link_target,      #c(5, 3, 4, 3, 0, 2, 2, 3),
                                             value =  link_value))       #c(1, 2, 1, 1, 1, 1, 1, 2)))
    fig_ch5_sankey_t4 <- fig_ch5_sankey_t4 %>% layout(title = paste0("Sankey Diagram of SNP groups (trait 4): ", names(mpp.pp$MPP)[4]))
    
    fig_ch5_sankey_t4
    
  })
  
  
  #Sankey Diagram in R and plotly (trait_5)----
  output$ch5_sankey_t5 <- renderPlotly({ 
    GWAS_all_final = tt()$GWAS_all_final
    mpp.pp = tt()$mpp.pp
    
    node_labels_snp             = GWAS_all_final$rs
    node_labels_groupfm         = unique(GWAS_all_final$snpGroups_fm)
    node_labels_groupfm_2       = paste("fm", node_labels_groupfm, sep="_")
    node_labels_groupflashfm    = unique(GWAS_all_final$snpGroups_flashfm)
    node_labels_groupflashfm_2  = paste("flashfm", node_labels_groupflashfm, sep="_")
    
    node_labels   = c(node_labels_snp, node_labels_groupfm, node_labels_groupflashfm)
    node_labels_2 = c(node_labels_snp, node_labels_groupfm_2, node_labels_groupflashfm_2)
    
    node_x      = c(rep(0.01,times=length(node_labels_snp)), 
                    rep(0.4,times=length(node_labels_groupfm)),
                    rep(0.8,times=length(node_labels_groupflashfm)))
    node_y      = c(seq(from = 1, to = length(node_labels_snp), by = 1)/length(node_labels_snp), 
                    seq(from = 1, to = length(node_labels_groupfm), by = 1)/length(node_labels_groupfm), 
                    seq(from = 1, to = length(node_labels_groupflashfm), by = 1)/length(node_labels_groupflashfm))
    
    node_table  = as.data.frame(cbind(node_labels, node_x, node_y));
    
    link_source_fm = c();
    link_target_fm = c();
    link_value_fm  = c();
    link_source_flashfm = c();
    link_target_flashfm = c();
    link_value_flashfm  = c();
    #group_fm
    for (i in 1:length(node_labels_snp)){
      if (GWAS_all_final$snpGroups_fm[i] != '0'){
        link_source_fm = c(link_source_fm, i-1)
        link_value_fm  = c(link_value_fm, 
                           mean(c(GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[5],"_Single")][i]),
                           na.rm=F) #updated20220208
        )
        for (ii in 1:length(node_labels_groupfm)){
          if (node_table$node_labels[length(node_labels_snp)+ii]==GWAS_all_final$snpGroups_fm[i]){
            link_target_fm = c(link_target_fm, length(node_labels_snp)+ii-1)
          }
        }
      }
    }
    #group_flashfm
    for (i in 1:length(node_labels_snp)){
      if (GWAS_all_final$snpGroups_flashfm[i] != '0'){
        link_source_flashfm = c(link_source_flashfm, i-1)
        link_value_flashfm  = c(link_value_flashfm, 
                                mean(c(GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[5],"_Multi")][i]),
                                na.rm=F) #updated20220208
        )
        for (ii in 1:length(node_labels_groupflashfm)){
          if (node_table$node_labels[length(c(node_labels_snp,node_labels_groupfm))+ii]==GWAS_all_final$snpGroups_flashfm[i]){
            link_target_flashfm = c(link_target_flashfm, length(c(node_labels_snp,node_labels_groupfm))+ii-1)
          }
        }
      }
    }
    #check direction carefully: fm->flashfm
    # link_source = c(link_source_fm, link_source_flashfm)
    # link_target = c(link_target_fm, link_target_flashfm)
    # link_value  = c(link_value_fm, link_value_flashfm)
    #check direction carefully: flashfm->fm
    link_source = c(link_source_flashfm, link_source_fm)
    link_target = c(link_target_flashfm, link_target_fm)
    link_value  = c(link_value_flashfm, link_value_fm)
    
    link_source_updated = link_source
    for (i in 2:length(link_source)){
      if (sum(link_source[i]==link_source[1:(i-1)])!=0){
        link_source_updated[i] = link_target[which(link_source[i]==link_source[1:(i-1)])]
      }
    }
    
    fig_ch5_sankey_t5 <- plot_ly(type = "sankey", arrangement = "snap", 
                                 domain = list(x =  c(0,1),y =  c(0,1)),
                                 orientation = "h",
                                 valueformat = "0.0000",
                                 #valuesuffix = "TWh",
                                 node = list(label = node_labels_2,pad = 15,
                                             thickness = 15,
                                             line = list(color = "black", width = 0.5)),
                                 link = list(source = link_source_updated,      #c(0, 0, 1, 2, 5, 4, 3, 5),
                                             target = link_target,      #c(5, 3, 4, 3, 0, 2, 2, 3),
                                             value =  link_value))       #c(1, 2, 1, 1, 1, 1, 1, 2)))
    fig_ch5_sankey_t5 <- fig_ch5_sankey_t5 %>% layout(title = paste0("Sankey Diagram of SNP groups (trait 5): ", names(mpp.pp$MPP)[5]))
    
    fig_ch5_sankey_t5
    
  })
  
  
  #Sankey Diagram in R and plotly (trait_6)----
  output$ch5_sankey_t6 <- renderPlotly({ 
    GWAS_all_final = tt()$GWAS_all_final
    mpp.pp = tt()$mpp.pp
    
    node_labels_snp             = GWAS_all_final$rs
    node_labels_groupfm         = unique(GWAS_all_final$snpGroups_fm)
    node_labels_groupfm_2       = paste("fm", node_labels_groupfm, sep="_")
    node_labels_groupflashfm    = unique(GWAS_all_final$snpGroups_flashfm)
    node_labels_groupflashfm_2  = paste("flashfm", node_labels_groupflashfm, sep="_")
    
    node_labels   = c(node_labels_snp, node_labels_groupfm, node_labels_groupflashfm)
    node_labels_2 = c(node_labels_snp, node_labels_groupfm_2, node_labels_groupflashfm_2)
    
    node_x      = c(rep(0.01,times=length(node_labels_snp)), 
                    rep(0.4,times=length(node_labels_groupfm)),
                    rep(0.8,times=length(node_labels_groupflashfm)))
    node_y      = c(seq(from = 1, to = length(node_labels_snp), by = 1)/length(node_labels_snp), 
                    seq(from = 1, to = length(node_labels_groupfm), by = 1)/length(node_labels_groupfm), 
                    seq(from = 1, to = length(node_labels_groupflashfm), by = 1)/length(node_labels_groupflashfm))
    
    node_table  = as.data.frame(cbind(node_labels, node_x, node_y));
    
    link_source_fm = c();
    link_target_fm = c();
    link_value_fm  = c();
    link_source_flashfm = c();
    link_target_flashfm = c();
    link_value_flashfm  = c();
    #group_fm
    for (i in 1:length(node_labels_snp)){
      if (GWAS_all_final$snpGroups_fm[i] != '0'){
        link_source_fm = c(link_source_fm, i-1)
        link_value_fm  = c(link_value_fm, 
                           mean(c(GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[6],"_Single")][i]),
                           na.rm=F)  #updated20220208
        )
        for (ii in 1:length(node_labels_groupfm)){
          if (node_table$node_labels[length(node_labels_snp)+ii]==GWAS_all_final$snpGroups_fm[i]){
            link_target_fm = c(link_target_fm, length(node_labels_snp)+ii-1)
          }
        }
      }
    }
    #group_flashfm
    for (i in 1:length(node_labels_snp)){
      if (GWAS_all_final$snpGroups_flashfm[i] != '0'){
        link_source_flashfm = c(link_source_flashfm, i-1)
        link_value_flashfm  = c(link_value_flashfm, 
                                mean(c(GWAS_all_final[ , paste0("MPP_",names(mpp.pp$MPP)[6],"_Multi")][i]),
                                na.rm=F)  #updated20220208
        )
        for (ii in 1:length(node_labels_groupflashfm)){
          if (node_table$node_labels[length(c(node_labels_snp,node_labels_groupfm))+ii]==GWAS_all_final$snpGroups_flashfm[i]){
            link_target_flashfm = c(link_target_flashfm, length(c(node_labels_snp,node_labels_groupfm))+ii-1)
          }
        }
      }
    }
    #check direction carefully: fm->flashfm
    # link_source = c(link_source_fm, link_source_flashfm)
    # link_target = c(link_target_fm, link_target_flashfm)
    # link_value  = c(link_value_fm, link_value_flashfm)
    #check direction carefully: flashfm->fm
    link_source = c(link_source_flashfm, link_source_fm)
    link_target = c(link_target_flashfm, link_target_fm)
    link_value  = c(link_value_flashfm, link_value_fm)
    
    link_source_updated = link_source
    for (i in 2:length(link_source)){
      if (sum(link_source[i]==link_source[1:(i-1)])!=0){
        link_source_updated[i] = link_target[which(link_source[i]==link_source[1:(i-1)])]
      }
    }
    
    fig_ch5_sankey_t6 <- plot_ly(type = "sankey", arrangement = "snap", 
                                 domain = list(x =  c(0,1),y =  c(0,1)),
                                 orientation = "h",
                                 valueformat = "0.0000",
                                 #valuesuffix = "TWh",
                                 node = list(label = node_labels_2,
                                             pad = 15,
                                             thickness = 15,
                                             line = list(color = "black", width = 0.5)),
                                 link = list(source = link_source_updated,      #c(0, 0, 1, 2, 5, 4, 3, 5),
                                             target = link_target,      #c(5, 3, 4, 3, 0, 2, 2, 3),
                                             value =  link_value))       #c(1, 2, 1, 1, 1, 1, 1, 2)))
    fig_ch5_sankey_t6 <- fig_ch5_sankey_t6 %>% layout(title = paste0("Sankey Diagram of SNP groups (trait 6): ", names(mpp.pp$MPP)[6]))
    
    fig_ch5_sankey_t6
    
  })
  
  #ch6: tables----
  output$web_inputs_cs1 <- renderPrint({
    cs1 = tt()$cs1
    str(cs1)
  })
  output$web_inputs_csM <- renderPrint({
    csM = tt()$csM
    str(csM)
  })
  output$web_inputs_GWAS <- renderPrint({
    GWAS = tt()$GWAS
    str(GWAS)
  })
  output$web_inputs_LD <- renderPrint({
    LD = tt()$LD
    str(LD)
  })
  output$web_inputs_mpp.pp <- renderPrint({
    mpp.pp = tt()$mpp.pp
    str(mpp.pp)
  })
  output$web_inputs_snpGroups <- renderPrint({
    snpGroups = tt()$snpGroups
    str(snpGroups)
  })
  output$GWAS_all_final_table_done <- DT::renderDataTable({
    GWAS_all_final = tt()$GWAS_all_final
    DT::datatable(head(data.frame(GWAS_all_final),10), 
                  options = list(dom = 't'))
  })
  output$LD_data<- DT::renderDataTable({
    LD = tt()$LD
    DT::datatable(data.frame(LD[1:5,1:5]),
                  options = list(dom = 't'))
  })

 
   
  # Downloadable csv of selected dataset ----
  output$downloadData_3 <- downloadHandler(
    filename = function() {
      paste("GWAS_all_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(tt()$GWAS_all_final, file, row.names = FALSE)
    }
  )
  
  # output$downloadData <- downloadHandler(
  #   filename = function() {
  #     paste("data-", Sys.Date(), ".csv", sep="")
  #   },
  #   content = function(file) {
  #     write.csv(data, file)
  #   }
  # )
 
  
} #End of Shiny server ##


####----Shiny app----
shinyApp(ui, server)