#'
#' all R code used to analysis data and produce figures for the manuscript:
#' "State-transition Modeling of Blood Transcriptome Predicts Disease Evolution and Treatment Response in Chronic Myeloid Leukemia"
#' https://www.biorxiv.org/content/10.1101/2023.10.11.561908v2
#' 
#' 
#' @name analysis+figures_reproducibility
#' @source https://github.com/cohmathonc/CML_mRNA_state-transition/R/
#' @author David Frankhouser
#'



###
### SETUP ### 
###
{
  #
  # load libraries, define helper functions, load data, define plot colors
  #
  
  ### load required libraries ###
  {
    library("ggplot2")
    library("SummarizedExperiment")
    library("DESeq2")
    library("pheatmap")
    library("tximport")
    library("msigdbr")
    library("fgsea")
    library("stringr")
    library("dplyr")
    library("enrichR")
    library("msigdb")
    library("reshape")
    library("GGally")
    library("VennDiagram")
    library("RColorBrewer")
    library("tidyverse")
    library(viridis)
    library("readxl")
    library(EnhancedVolcano)
    library(STRINGdb)
    library(nVennR)
    library("biomaRt")
    library(ComplexHeatmap)
  }
  
  ### define helper functions ###
  {
    ###get gene locations from data ###
    get_genes_ind <- function(in.se, in.genes, gene.key = "gene_name", all.matches=F) {
      out.sel <- c()
      for (g in in.genes) {
        m <- which(rowData(in.se)[[gene.key]] == g)
        if (length(m)==1) {
          out.sel <- c(out.sel, m)
        } else if (length(m)==0) {
          out.sel <- c(out.sel, NA)
          print(paste(":::WARNING::: No gene matched with name: ",g,sep=""))
        } else if (length(m) > 1) {
          if (all.matches) {
            out.sel <- c(out.sel, m)
          } else {
            # ::note:: this could be modified to get the longest transcript or prefer "protein coding" transcripts
            out.sel <- c(out.sel, m[0])
            print(paste(":::WARNING::: Keeping first of multiple matches for gene with name: ",g,sep=""))
          }
        }
      }  
      return(out.sel)
    }
    
    
    ### get matching genes from expression
    #   return index of genes (e.g. expression table) that were matched from a gene set (e.g. pathway)
    #       "toMatch": list of genes for which an index will be returned for all matched genes in "geneSet"
    #       "geneSet": character list of genes to be searched for in "toMatcn"
    ### returns only first match
    match_genes <- function(matchTable, geneSet, return.names=F) {
      toMatch <- matchTable
      if (length(toMatch) < length(geneSet)) { #print warning if gene sets might be reversed...
        print(paste(":::WARNING::: input gene set is larger than gene list to get matches from. Be sure the order of arguements are correct!!!",sep=""))
      }
      
      csel <- c()
      gnames <- c()
      for (g in geneSet) {
        m <- which(toMatch==g)
        if (length(m) == 1) {
          csel <- c(csel, m)
          gnames <- c(gnames, g)
        } else if (length(m) > 1) {
          print(":::WARNING::: matched multiple genes; keeping first match only")
          csel <- c(csel, m[1])
          gnames <- c(gnames, g)
        }
      }
      if (return.names) {
        return(gnames) 
      } else {
        return(csel)
      }
    } #:::end match genes function
    #returns all matches
    match_all_genes <- function(matchTable, geneSet, return.names=F) {
      toMatch <- matchTable
      if (length(toMatch) < length(geneSet)) { #print warning if gene sets might be reversed...
        print(paste(":::WARNING::: input gene set is larger than gene list to get matches from. Be sure the order of arguements are correct!!!",sep=""))
      }
      
      csel <- c()
      gnames <- c()
      for (g in geneSet) {
        m <- which(toMatch==g)
        if (length(m) == 1) {
          csel <- c(csel, m)
          gnames <- c(gnames, g)
        } else if (length(m) > 1) {
          csel <- c(csel, m)
          gnames <- c(gnames, rep(g,length(m)) )
        }
      }
      if (return.names) {
        return(csel, gnames) 
      } else {
        return(csel)
      }
    } #:::end match genes function
    
    
    ###
    ### fGSEA FUNCTIONs
    ###
    run_fgsea <- function(inFCs, in_dir) {
      cur_path_out <- paste(in_dir,sep="")
      dir.create(cur_path_out, showWarnings = F)
      cur.lfc <- inFCs
      cur.lfc <- sort(cur.lfc)
      #set pathway loop objects
      cats <- c("H", "C2", "C2", "C3", "C3", "C5", "C5", "C5")
      subcats <- c(NA,  "CP:KEGG", "CP:WIKIPATHWAYS", "TFT:GTRD", "MIR:MIRDB", "GO:BP", "GO:MF", "GO:CC")
      #cats <- c("H", "C2", "C5")
      #subcats <- c(NA, "CP:WIKIPATHWAYS", "GO:BP")
      sort_lfc <- cur.lfc
      for (i in 1:length(cats)) {
        curc <- cats[i]
        curs <- subcats[i]
        if ( curc == "H" ) {
          curpath = msigdbr(species = "mouse", category = curc)
          pathname <- "Hallmark"
        } else {
          curpath = msigdbr(species = "mouse", category = curc, subcategory=curs)
          pathname <- gsub(":","-",curs)
        }
        print(paste("Processing: ",pathname,sep=""))
        path_list = split(x = curpath$gene_symbol, f = curpath$gs_name)
        cur_fgsea <- fgsea(pathways=path_list, stats=sort_lfc )
        data.table::fwrite(cur_fgsea[order(cur_fgsea$padj),], file=paste(cur_path_out,"/fgsea_",pathname,"_table.tsv",sep=""), sep="\t" )
        sig.sel <- which(cur_fgsea$padj<=0.1)
        sig_path <- list()
        if (length(sig.sel)==0) {next}
        le_out <- paste(cur_path_out,"/leadingEdge",sep="")
        dir.create(le_out, showWarnings = F)
        for (s in 1:length(sig.sel)) {
          curp <- cur_fgsea[s,]$pathway
          sig_path[[curp]] <- path_list[[curp]]
          png(paste(le_out,"/fgsea_",pathname,"_Obj-",curp,"_leadingEdge.png",sep=""), res=plot_res, units="in", height=6, width=6)
          p <- plotEnrichment(path_list[[curp]],sort_lfc) + labs(title=curp)
          print(p)
          graphics.off()
        }
        plottab <- cur_fgsea[sig.sel,]
        if (dim(plottab)[1] > 20) { 
          plottab <- plottab[sort.int(plottab$padj, index.return = T)$ix,][seq(1,20),] 
          plot.gl <- list()
          for (p in plottab$pathway) {
            plot.gl[[p]] <- path_list[[p]]
          }
        }
        #%>% arrange(factor(pathway, levels = plottab[["pathway"]][order(plottab[["padj"]], decreasing=F)]))[seq(1,20),] }
        png(paste(cur_path_out,"/fgsea_",pathname,"_sigTable.png",sep=""), res=plot_res, units="in", height=8, width=6)
        p <- plotGseaTable(plot.gl, sort_lfc, plottab)
        print(p)
        graphics.off()
        # dot plot #
        sigsel <- which(cur_fgsea[["padj"]]<0.05)
        sigdf <- data.frame("Pathway"= unlist(lapply(cur_fgsea[["pathway"]][sigsel], function(x) gsub("_", " ", x) )),
                            "qvalue"=cur_fgsea[["padj"]][sigsel], "Count"=cur_fgsea[["size"]][sigsel], "GeneRatio"=cur_fgsea[["ES"]][sigsel])
        sigdf$type <- "Upregulated"
        sigdf$type[sigdf$GeneRatio < 0 ] <- "Downregulated"
        # format pathways names
        fpath <- unlist(lapply(sigdf$Pathway, function(x) str_to_title( gsub("HALLMARK ","", x) ) ))
        sigdf$Pathway <- fpath
        sigdf <- sigdf %>% arrange(factor(Pathway, levels = sigdf[["Pathway"]][order(sigdf[["qvalue"]], decreasing=F)]))
        #ggplot(data=sigdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=16)
        #orddf <- sigdf %>% mutate(Pathway = fct_reorder(Pathway, qvalue, .desc=T)) #idk what this does...
        plotdf <- sigdf
        if (dim(plotdf)[1]>20) { plotdf <- plotdf[seq(1,20),]}
        plotdf$Pathway <- factor(plotdf$Pathway, levels=rev(plotdf$Pathway))
        png(paste(cur_path_out,"/fgsea_",pathname,"_dotPlot.png",sep=""), res=plot_res, units="in", height=8, width=12)
        p <- ggplot(plotdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=18) + scale_color_gradient(limits=c(0,.05), low="red", high="blue")
        print(p)
        dev.off()
      } #end fgsea for loop
    } #end fgsea function
    run_quick_fgsea <- function(inFCs, in_dir) {
      
      cur_path_out <- paste(in_dir,sep="")
      dir.create(cur_path_out, showWarnings = F)
      cur.lfc <- inFCs
      cur.lfc <- sort(cur.lfc)
      #set pathway loop objects
      #cats <- c("H", "C2", "C2", "C3", "C3", "C5", "C5", "C5")
      #subcats <- c(NA,  "CP:KEGG", "CP:WIKIPATHWAYS", "TFT:GTRD", "MIR:MIRDB", "GO:BP", "GO:MF", "GO:CC")
      #cats <- c("H", "C2", "C5")
      #subcats <- c(NA, "CP:WIKIPATHWAYS", "GO:BP")
      cats <- c("H", "C2")
      subcats <- c(NA, "CP:WIKIPATHWAYS")
      sort_lfc <- cur.lfc
      for (i in 1:length(cats)) {
        curc <- cats[i]
        curs <- subcats[i]
        if ( curc == "H" ) {
          curpath = msigdbr(species = "mouse", category = curc)
          pathname <- "Hallmark"
        } else {
          curpath = msigdbr(species = "mouse", category = curc, subcategory=curs)
          pathname <- gsub(":","-",curs)
        }
        print(paste("Processing: ",pathname,sep=""))
        path_list = split(x = curpath$gene_symbol, f = curpath$gs_name)
        cur_fgsea <- fgsea(pathways=path_list, stats=sort_lfc )
        data.table::fwrite(cur_fgsea[order(cur_fgsea$padj),], file=paste(cur_path_out,"/fgsea_",pathname,"_table.tsv",sep=""), sep="\t" )
        sig.sel <- which(cur_fgsea$padj<=0.1)
        sig_path <- list()
        if (length(sig.sel)==0) {next}
        le_out <- paste(cur_path_out,"/leadingEdge",sep="")
        dir.create(le_out, showWarnings = F)
        #if (curs != "GO:BP") {}
        for (s in 1:length(sig.sel)) {
          curp <- cur_fgsea[s,]$pathway
          sig_path[[curp]] <- path_list[[curp]]
          png(paste(le_out,"/fgsea_",pathname,"_Obj-",curp,"_leadingEdge.png",sep=""), res=plot_res, units="in", height=6, width=6)
          p <- plotEnrichment(path_list[[curp]],sort_lfc) + labs(title=curp)
          print(p)
          graphics.off()
        }
        
        plottab <- cur_fgsea[sig.sel,]
        plot.gl <- list()
        if (dim(plottab)[1] > 20) { 
          plottab <- plottab[sort.int(plottab$padj, index.return = T)$ix,][seq(1,20),] 
          for (p in plottab$pathway) {
            plot.gl[[p]] <- path_list[[p]]
          }
        }
        #exit if nothing to print; not sure why this would be empty, but errors...
        if (length(plot.gl)>0) {
          #%>% arrange(factor(pathway, levels = plottab[["pathway"]][order(plottab[["padj"]], decreasing=F)]))[seq(1,20),] }
          png(paste(cur_path_out,"/fgsea_",pathname,"_sigTable.png",sep=""), res=plot_res, units="in", height=8, width=6)
          p <- plotGseaTable(plot.gl, sort_lfc, plottab)
          print(p)
          graphics.off()
        } else  {
          if (dim(plottab)[1] > 0 ) {
            png(paste(cur_path_out,"/fgsea_",pathname,"_sigTable.png",sep=""), res=plot_res, units="in", height=8, width=6)
            p <- plotGseaTable(plottab, sort_lfc, plottab)
            print(p)
            graphics.off()
          }
        }
        # dot plot #
        sigsel <- which(cur_fgsea[["padj"]]<0.05)
        if (length(sigsel)==0) { #if none significant, take top 5
          sigsel <- order(cur_fgsea[["padj"]])[seq(1,5)]
        }
        sigdf <- data.frame("Pathway"= unlist(lapply(cur_fgsea[["pathway"]][sigsel], function(x) gsub("_", " ", x) )),
                            "qvalue"=cur_fgsea[["padj"]][sigsel], "Count"=cur_fgsea[["size"]][sigsel], "GeneRatio"=cur_fgsea[["ES"]][sigsel], "NES"=cur_fgsea[["NES"]][sigsel])
        sigdf$type <- "Upregulated"
        sigdf$type[sigdf$GeneRatio < 0 ] <- "Downregulated"
        # format pathways names
        fpath <- unlist(lapply(sigdf$Pathway, function(x) str_to_title( gsub("HALLMARK ","", x) ) ))
        sigdf$Pathway <- fpath
        sigdf <- sigdf %>% arrange(factor(Pathway, levels = sigdf[["Pathway"]][order(sigdf[["qvalue"]], decreasing=F)]))
        #ggplot(data=sigdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=16)
        #orddf <- sigdf %>% mutate(Pathway = fct_reorder(Pathway, qvalue, .desc=T)) #idk what this does...
        plotdf <- sigdf
        if (dim(plotdf)[1]>20) { plotdf <- plotdf[seq(1,20),]}
        plotdf$Pathway <- factor(plotdf$Pathway, levels=rev(plotdf$Pathway))
        png(paste(cur_path_out,"/fgsea_",pathname,"_dotPlot.png",sep=""), res=plot_res, units="in", height=8, width=12)
        p <- ggplot(plotdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=18) + scale_color_gradient(limits=c(0,.05), low="red", high="blue")
        print(p)
        dev.off()
        png(paste(cur_path_out,"/fgsea_",pathname,"_NES_dotPlot.png",sep=""), res=plot_res, units="in", height=8, width=12)
        p <- ggplot(plotdf, aes(x=NES, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=18) + scale_color_gradient(limits=c(0,.05), low="red", high="blue")
        print(p)
        dev.off()
        png(paste(cur_path_out,"/fgsea_",pathname,"_NES_ge-abs-2_dotPlot.png",sep=""), res=plot_res, units="in", height=8, width=12)
        p <- ggplot(plotdf[which(abs(plotdf$NES)>=2),], aes(x=NES, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=18) + scale_color_gradient(limits=c(0,.05), low="red", high="blue")
        print(p)
        dev.off()
      } #end fgsea for loop
    } #end fgsea function
    single_sample_fgsea <- function(inFCs, in_dir, prefix) {
      cur_path_out <- paste(in_dir,sep="")
      dir.create(cur_path_out, showWarnings = F)
      cur.lfc <- inFCs
      cur.lfc <- sort(cur.lfc)
      #set pathway loop objects
      #cats <- c("H", "C2", "C2", "C3", "C3", "C5", "C5", "C5")
      #subcats <- c(NA,  "CP:KEGG", "CP:WIKIPATHWAYS", "TFT:GTRD", "MIR:MIRDB", "GO:BP", "GO:MF", "GO:CC")
      cats <- c("H", "C2", "C5")
      subcats <- c(NA, "CP:WIKIPATHWAYS", "GO:BP")
      sort_lfc <- cur.lfc
      #for (i in 1:length(cats)) {
      for (i in 1:1) {
        curc <- cats[i]
        curs <- subcats[i]
        if ( curc == "H" ) {
          curpath = msigdbr(species = "mouse", category = curc)
          pathname <- "Hallmark"
        } else {
          curpath = msigdbr(species = "mouse", category = curc, subcategory=curs)
          pathname <- gsub(":","-",curs)
        }
        print(paste("Processing: ",pathname,sep=""))
        path_list = split(x = curpath$gene_symbol, f = curpath$gs_name)
        cur_fgsea <- fgsea(pathways=path_list, stats=sort_lfc, nPermSimple = 10000 )
        data.table::fwrite(cur_fgsea[order(cur_fgsea$padj),], file=paste(cur_path_out,"/fgsea_",prefix,"_db-",pathname,"_table.tsv",sep=""), sep="\t" )
        sig.sel <- which(cur_fgsea$padj<=0.1)
        sig_path <- list()
        if (length(sig.sel)==0) {next}
        #le_out <- paste(cur_path_out,"/leadingEdge",sep="")
        #dir.create(le_out, showWarnings = F)
        #for (s in 1:length(sig.sel)) {
        #  curp <- cur_fgsea[s,]$pathway
        #  sig_path[[curp]] <- path_list[[curp]]
        #  png(paste(le_out,"/fgsea_",pathname,"_Obj-",curp,"_leadingEdge.png",sep=""), res=plot_res, units="in", height=6, width=6)
        #  p <- plotEnrichment(path_list[[curp]],sort_lfc) + labs(title=curp)
        #  print(p)
        #  graphics.off()
        #}
        plottab <- cur_fgsea[sig.sel,]
        if (dim(plottab)[1] > 20) { 
          plottab <- plottab[sort.int(plottab$padj, index.return = T)$ix,][seq(1,20),] 
          plot.gl <- list()
          for (p in plottab$pathway) {
            plot.gl[[p]] <- path_list[[p]]
          }
        }
        #%>% arrange(factor(pathway, levels = plottab[["pathway"]][order(plottab[["padj"]], decreasing=F)]))[seq(1,20),] }
        #png(paste(cur_path_out,"/fgsea_",prefix,"_db-",pathname,"_sigTable.png",sep=""), res=plot_res, units="in", height=8, width=6)
        #p <- plotGseaTable(plot.gl, sort_lfc, plottab)
        #print(p)
        #graphics.off()
        # dot plot #
        sigsel <- which(cur_fgsea[["padj"]]<0.05)
        sigdf <- data.frame("Pathway"= unlist(lapply(cur_fgsea[["pathway"]][sigsel], function(x) gsub("_", " ", x) )),
                            "qvalue"=cur_fgsea[["padj"]][sigsel], "Count"=cur_fgsea[["size"]][sigsel], "GeneRatio"=cur_fgsea[["ES"]][sigsel])
        sigdf$type <- "Upregulated"
        sigdf$type[sigdf$GeneRatio < 0 ] <- "Downregulated"
        # format pathways names
        fpath <- unlist(lapply(sigdf$Pathway, function(x) str_to_title( gsub("HALLMARK ","", x) ) ))
        sigdf$Pathway <- fpath
        sigdf <- sigdf %>% arrange(factor(Pathway, levels = sigdf[["Pathway"]][order(sigdf[["qvalue"]], decreasing=F)]))
        #ggplot(data=sigdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=16)
        #orddf <- sigdf %>% mutate(Pathway = fct_reorder(Pathway, qvalue, .desc=T)) #idk what this does...
        plotdf <- sigdf
        if (dim(plotdf)[1]>20) { plotdf <- plotdf[seq(1,20),]}
        plotdf$Pathway <- factor(plotdf$Pathway, levels=rev(plotdf$Pathway))
        png(paste(cur_path_out,"/fgsea_",prefix,"_db-",pathname,"_dotPlot.png",sep=""), res=plot_res, units="in", height=8, width=12)
        p <- ggplot(plotdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=18) + scale_color_gradient(limits=c(0,.05), low="red", high="blue")
        print(p)
        dev.off()
      } #end fgsea for loop
    } #end fgsea function
    single_sample_fgsea_orig <- function(inFCs, in_dir, prefix) {
      cur_path_out <- paste(in_dir,sep="")
      dir.create(cur_path_out, showWarnings = F)
      cur.lfc <- inFCs
      cur.lfc <- sort(cur.lfc)
      #set pathway loop objects
      #cats <- c("H", "C2", "C2", "C3", "C3", "C5", "C5", "C5")
      #subcats <- c(NA,  "CP:KEGG", "CP:WIKIPATHWAYS", "TFT:GTRD", "MIR:MIRDB", "GO:BP", "GO:MF", "GO:CC")
      cats <- c("H", "C2", "C5")
      subcats <- c(NA, "CP:WIKIPATHWAYS", "GO:BP")
      sort_lfc <- cur.lfc
      #for (i in 1:length(cats)) {
      for (i in 1:1) {
        curc <- cats[i]
        curs <- subcats[i]
        if ( curc == "H" ) {
          curpath = msigdbr(species = "mouse", category = curc)
          pathname <- "Hallmark"
        } else {
          curpath = msigdbr(species = "mouse", category = curc, subcategory=curs)
          pathname <- gsub(":","-",curs)
        }
        print(paste("Processing: ",pathname,sep=""))
        path_list = split(x = curpath$gene_symbol, f = curpath$gs_name)
        cur_fgsea <- fgsea(pathways=path_list, stats=sort_lfc )
        data.table::fwrite(cur_fgsea[order(cur_fgsea$padj),], file=paste(cur_path_out,"/fgsea_",prefix,"_db-",pathname,"_table.tsv",sep=""), sep="\t" )
        sig.sel <- which(cur_fgsea$padj<=0.1)
        sig_path <- list()
        if (length(sig.sel)==0) {next}
        #le_out <- paste(cur_path_out,"/leadingEdge",sep="")
        #dir.create(le_out, showWarnings = F)
        #for (s in 1:length(sig.sel)) {
        #  curp <- cur_fgsea[s,]$pathway
        #  sig_path[[curp]] <- path_list[[curp]]
        #  png(paste(le_out,"/fgsea_",pathname,"_Obj-",curp,"_leadingEdge.png",sep=""), res=plot_res, units="in", height=6, width=6)
        #  p <- plotEnrichment(path_list[[curp]],sort_lfc) + labs(title=curp)
        #  print(p)
        #  graphics.off()
        #}
        plottab <- cur_fgsea[sig.sel,]
        if (dim(plottab)[1] > 20) { 
          plottab <- plottab[sort.int(plottab$padj, index.return = T)$ix,][seq(1,20),] 
          plot.gl <- list()
          for (p in plottab$pathway) {
            plot.gl[[p]] <- path_list[[p]]
          }
        }
        #%>% arrange(factor(pathway, levels = plottab[["pathway"]][order(plottab[["padj"]], decreasing=F)]))[seq(1,20),] }
        #png(paste(cur_path_out,"/fgsea_",prefix,"_db-",pathname,"_sigTable.png",sep=""), res=plot_res, units="in", height=8, width=6)
        #p <- plotGseaTable(plot.gl, sort_lfc, plottab)
        #print(p)
        #graphics.off()
        # dot plot #
        sigsel <- which(cur_fgsea[["padj"]]<0.05)
        sigdf <- data.frame("Pathway"= unlist(lapply(cur_fgsea[["pathway"]][sigsel], function(x) gsub("_", " ", x) )),
                            "qvalue"=cur_fgsea[["padj"]][sigsel], "Count"=cur_fgsea[["size"]][sigsel], "GeneRatio"=cur_fgsea[["ES"]][sigsel])
        sigdf$type <- "Upregulated"
        sigdf$type[sigdf$GeneRatio < 0 ] <- "Downregulated"
        # format pathways names
        fpath <- unlist(lapply(sigdf$Pathway, function(x) str_to_title( gsub("HALLMARK ","", x) ) ))
        sigdf$Pathway <- fpath
        sigdf <- sigdf %>% arrange(factor(Pathway, levels = sigdf[["Pathway"]][order(sigdf[["qvalue"]], decreasing=F)]))
        #ggplot(data=sigdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=16)
        #orddf <- sigdf %>% mutate(Pathway = fct_reorder(Pathway, qvalue, .desc=T)) #idk what this does...
        plotdf <- sigdf
        if (dim(plotdf)[1]>20) { plotdf <- plotdf[seq(1,20),]}
        plotdf$Pathway <- factor(plotdf$Pathway, levels=rev(plotdf$Pathway))
        png(paste(cur_path_out,"/fgsea_",prefix,"_db-",pathname,"_dotPlot.png",sep=""), res=plot_res, units="in", height=8, width=12)
        p <- ggplot(plotdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=18) + scale_color_gradient(limits=c(0,.05), low="red", high="blue")
        print(p)
        dev.off()
      } #end fgsea for loop
    } #end fgsea function
    single_sample_fgsea_simple <- function(inFCs, in_dir, prefix) { #explicitly set multilevel
      cur_path_out <- paste(in_dir,sep="")
      dir.create(cur_path_out, showWarnings = F)
      cur.lfc <- inFCs
      cur.lfc <- sort(cur.lfc)
      #set pathway loop objects
      #cats <- c("H", "C2", "C2", "C3", "C3", "C5", "C5", "C5")
      #subcats <- c(NA,  "CP:KEGG", "CP:WIKIPATHWAYS", "TFT:GTRD", "MIR:MIRDB", "GO:BP", "GO:MF", "GO:CC")
      cats <- c("H", "C2", "C5")
      subcats <- c(NA, "CP:WIKIPATHWAYS", "GO:BP")
      sort_lfc <- cur.lfc
      #for (i in 1:length(cats)) {
      for (i in 1:1) {
        curc <- cats[i]
        curs <- subcats[i]
        if ( curc == "H" ) {
          curpath = msigdbr(species = "mouse", category = curc)
          pathname <- "Hallmark"
        } else {
          curpath = msigdbr(species = "mouse", category = curc, subcategory=curs)
          pathname <- gsub(":","-",curs)
        }
        print(paste("Processing: ",pathname,sep=""))
        path_list = split(x = curpath$gene_symbol, f = curpath$gs_name)
        cur_fgsea <- fgseaSimple(pathways=path_list, stats=sort_lfc, nperm = 10000, nproc=4 )
        data.table::fwrite(cur_fgsea[order(cur_fgsea$padj),], file=paste(cur_path_out,"/fgsea_",prefix,"_db-",pathname,"_table.tsv",sep=""), sep="\t" )
        sig.sel <- which(cur_fgsea$padj<=0.1)
        sig_path <- list()
        if (length(sig.sel)==0) {next}
        #le_out <- paste(cur_path_out,"/leadingEdge",sep="")
        #dir.create(le_out, showWarnings = F)
        #for (s in 1:length(sig.sel)) {
        #  curp <- cur_fgsea[s,]$pathway
        #  sig_path[[curp]] <- path_list[[curp]]
        #  png(paste(le_out,"/fgsea_",pathname,"_Obj-",curp,"_leadingEdge.png",sep=""), res=plot_res, units="in", height=6, width=6)
        #  p <- plotEnrichment(path_list[[curp]],sort_lfc) + labs(title=curp)
        #  print(p)
        #  graphics.off()
        #}
        plottab <- cur_fgsea[sig.sel,]
        if (dim(plottab)[1] > 20) { 
          plottab <- plottab[sort.int(plottab$padj, index.return = T)$ix,][seq(1,20),] 
          plot.gl <- list()
          for (p in plottab$pathway) {
            plot.gl[[p]] <- path_list[[p]]
          }
        }
        #%>% arrange(factor(pathway, levels = plottab[["pathway"]][order(plottab[["padj"]], decreasing=F)]))[seq(1,20),] }
        #png(paste(cur_path_out,"/fgsea_",prefix,"_db-",pathname,"_sigTable.png",sep=""), res=plot_res, units="in", height=8, width=6)
        #p <- plotGseaTable(plot.gl, sort_lfc, plottab)
        #print(p)
        #graphics.off()
        # dot plot #
        sigsel <- which(cur_fgsea[["padj"]]<0.05)
        sigdf <- data.frame("Pathway"= unlist(lapply(cur_fgsea[["pathway"]][sigsel], function(x) gsub("_", " ", x) )),
                            "qvalue"=cur_fgsea[["padj"]][sigsel], "Count"=cur_fgsea[["size"]][sigsel], "GeneRatio"=cur_fgsea[["ES"]][sigsel])
        sigdf$type <- "Upregulated"
        sigdf$type[sigdf$GeneRatio < 0 ] <- "Downregulated"
        # format pathways names
        fpath <- unlist(lapply(sigdf$Pathway, function(x) str_to_title( gsub("HALLMARK ","", x) ) ))
        sigdf$Pathway <- fpath
        sigdf <- sigdf %>% arrange(factor(Pathway, levels = sigdf[["Pathway"]][order(sigdf[["qvalue"]], decreasing=F)]))
        #ggplot(data=sigdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=16)
        #orddf <- sigdf %>% mutate(Pathway = fct_reorder(Pathway, qvalue, .desc=T)) #idk what this does...
        plotdf <- sigdf
        if (dim(plotdf)[1]>20) { plotdf <- plotdf[seq(1,20),]}
        plotdf$Pathway <- factor(plotdf$Pathway, levels=rev(plotdf$Pathway))
        png(paste(cur_path_out,"/fgsea_",prefix,"_db-",pathname,"_dotPlot.png",sep=""), res=plot_res, units="in", height=8, width=12)
        p <- ggplot(plotdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=18) + scale_color_gradient(limits=c(0,.05), low="red", high="blue")
        print(p)
        dev.off()
      } #end fgsea for loop
    } #end fgsea function
    
    
    # enrichR setup #
    setEnrichrSite("Enrichr")
    websiteLive <- TRUE
    #dbs <- listEnrichrDbs() #get all dbs
    #dbs[grep("Hallmark", dbs$libraryName),]
    test_dbs <- c("GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021", "KEGG_2019_Mouse", 
                  "WikiPathways_2019_Mouse", "MSigDB_Hallmark_2020", "MSigDB_Oncogenic_Signatures", 
                  "TF_Perturbations_Followed_by_Expression", "TRRUST_Transcription_Factors_2019", "ChEA_2022")
    output_enrichr_results <- function(inres, name, path_out) {
      for ( db in names(inres)) {
        if (dim(inres[[db]])[1]==0) {next} #skip output
        write.table(inres[[db]], paste(path_out,"/",db,"_",name,"_table.tsv", sep="" ), sep="\t", row.names=F )
        png(paste(path_out,"/",db,"_",name,"_plot.png", sep="" ), res=plot_res, units="in", height=12, width=8 )
        p <- plotEnrich(inres[[db]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
        print(p)
        graphics.off()
      }
    }
    
    
    ### makes list of 3-way-venn sector infomation
    # values:
    #   "sets" contains all objects represented by each sector of the 3 way venn
    #   "sizes" contains the size of each sector of the 3 way venn
    venn3_from_named_list <- function(inlist) {
      ### three way venn ###
      n1 <- names(inlist)[1]
      n2 <- names(inlist)[2]
      n3 <- names(inlist)[3]
      #get lists for each category
      l1 <- inlist[[n1]]
      l2 <- inlist[[n2]]
      l3 <- inlist[[n3]]
      #get lenghts of lists
      sa1 <- length(l1)
      sa2 <- length(l2)
      sa3 <- length(l3)
      
      #get objects in each venn sector
      a1 <- setdiff(l1, c(l2,l3))
      a2 <- setdiff(l2, c(l1,l3))
      a3 <- setdiff(l3, c(l2,l1))
      a12 <- intersect(l1, l2)
      a23 <- intersect(l2, l3)
      a13 <- intersect(l1, l3)
      a123 <- intersect( intersect(l1,l2), l3 )
      #get size of sectors
      n1 <- length(a1)
      n2 <- length(a2)
      n3 <- length(a3)
      n12 <- length(a12)
      n23 <- length(a23)
      n13 <- length(a13)
      n123 <- length(a123)
      #make output
      out <- list()
      #append objects
      setnames <- c("a1", "a2", "a3", "a12", "a13", "a23", "a123")
      sets <- list(a1, a2, a3, a12, a13, a23, a123)
      for (si in 1:length(setnames) ) {
        s <- setnames[si]
        out[["sets"]][[s]] <- sets[[si]]
      }
      sizenames <- c("n1", "n2", "n3", "n12", "n13", "n23", "n123")
      sizes <- list(n1, n2, n3, n12, n13, n23, n123)
      for (si in 1:length(sizenames) ) {
        s <- sizenames[si]
        out[["sizes"]][[s]] <- sizes[[si]]
      }
      
      return(out)
    }
    
    
    ### stringDB HELPER FUNCTIONS
    {
      #function to keep only stringdb gene ids in an input list
      rmAllUnmatchedGenes <- function(in.int, toMatch) {
        out.int <- c()
        for (i in 1:dim(in.int)[1]) {
          curline <- in.int[i,]
          #match base gene "from"
          mf <- which(toMatch == curline$from)
          #match target gene "to"
          mt <- which(toMatch == curline$to)
          #remove unless both matched
          if ( length(mf) > 0 & length(mt) > 0 ) {
            out.int <- rbind(out.int, curline)
          }
        }
        print(paste("retained ",dim(out.int)[1]," of ",dim(in.int)[1]," gene interactions"))
        return(out.int)
      }
      #remove only targets not in list
      rmUnmatchedTargetGenes <- function(in.int, toMatch) {
        out.int <- c()
        for (i in 1:dim(in.int)[1]) {
          curline <- in.int[i,]
          #match base gene "from"
          mf <- which(toMatch == curline$from)
          #match target gene "to"
          mt <- which(toMatch == curline$to)
          #remove unless both matched
          if ( length(mt) > 0 ) {
            out.int <- rbind(out.int, curline)
          }
        }
        print(paste("retained ",dim(out.int)[1]," of ",dim(in.int)[1]," gene interactions"))
        return(out.int)
      }
      #function to add logFC column to table used for stringDB
      add_l2fc <- function(in.tab, deg.tab) {
        out.tab <- c()
        for (i in 1:dim(in.tab)[1]) {
          cgene <- in.tab$gene[i]
          m <- which(deg.tab$gene_name==cgene)
          #get l2fc value
          if (length(m)==1) {
            l2fc <- deg.tab$log2FoldChange[m]
          } else if (length(m) > 1 ) {
            fcs <- deg.tab$log2FoldChange[m]
            if ( length(which(fcs>0))==length(fcs) | length(which(fcs<0))==length(fcs) ) { #handle all same dir
              l2fc <- mean(fcs)
            } else {
              print(paste(":::WARNING::: multi-direction gene: ",cgene,sep=""))
              l2fc <- fcs[which(abs(fcs)==max(abs(fcs)))]
            }
          } else {
            print(paste(":::WARNING::: could not match: ",cgene,sep=""))
            l2fc <- 0
          }
          #add l2FC to output
          out.tab <- rbind(out.tab, c(unlist(in.tab[i,]), l2fc) )
        }
        #add colnames
        colnames(out.tab) <- c(colnames(in.tab), "l2FC")
        
        return(data.frame(out.tab))
      }
      add_l2fc_wColor <- function(in.tab, deg.tab, minColor="blue", maxColor="red", lfcMax=8, lfcMin=-8 ) {
        # in.tab <- iec.c2.orig
        # deg.tab <- all.degs
        ###make color function
        # Define the color range
        color_min <- lfcMin
        color_max <- lfcMax
        # Define the color scale
        color_map <- colorRamp(c(minColor, maxColor))
        color_map <- colorRamp(rev(brewer.pal(11,"RdBu")))
        #use function to get color; sets a reasonable color scale and handles expreme log2FCs
        get_col <- function(lfc, c_map=color_map, minVal=color_min, maxVal=color_max, maxC =maxColor, minC=minColor) {
          norm_lfc <- ( lfc - minVal) / (maxVal - minVal)
          if (norm_lfc > 1) {
            col <- maxC
          } else if (norm_lfc < 0) {
            col <- minC
          } else {
            col <- rgb(c_map(norm_lfc), maxColorValue = 255)  
          }
          return (col)
        }
        
        out.tab <- c()
        for (i in 1:dim(in.tab)[1]) {
          cgene <- in.tab$gene[i]
          m <- which(deg.tab$gene_name==cgene)
          #get l2fc value
          if (length(m)==1) {
            l2fc <- as.numeric(deg.tab$log2FoldChange[m])
            col <- get_col(l2fc)
          } else if (length(m) > 1 ) {
            fcs <- as.numeric(deg.tab$log2FoldChange[m])
            if ( length(which(fcs>0))==length(fcs) | length(which(fcs<0))==length(fcs) ) { #handle all same dir
              l2fc <- mean(fcs)
              col <- get_col(l2fc)
            } else {
              print(paste(":::WARNING::: multi-direction gene: ",cgene,sep=""))
              l2fc <- fcs[which(abs(fcs)==max(abs(fcs)))]
              col <- "black"
            }
          } else {
            print(paste(":::WARNING::: could not match: ",cgene,sep=""))
            l2fc <- 0
            col <- "grey60"
          }
          #add l2FC  + color to output
          out.tab <- rbind(out.tab, c(unlist(in.tab[i,]), l2fc, col) )
        }
        #add colnames
        dim(in.tab)
        colnames(out.tab) <- c(colnames(in.tab), "l2FC", "color")
        out.df <- data.frame(out.tab)
        
        return(out.df)
      }
      
      
    }
    
    
    
    # Create a matrix with the desired dimensions
    matrix_out <- function(my_list, num_cols) {
      num_items <- length(my_list)
      num_rows <- ceiling(num_items / num_cols)
      my_matrix <- matrix("", nrow = num_rows, ncol = num_cols)
      my_matrix[1:num_items] <- my_list
      return(my_matrix) 
    }
    
    #
    # function to get DEGs using gene_name 
    #
    match_deg_genes <- function(deg.names, deg.obj) {
      match.ind <- c()
      for (g in deg.names) {
        m <- which(deg.obj$gene_name==g & deg.obj$padj <= 0.05)
        if ( length(m) < 1 ) {
          print(paste("::Warning:: No DEG matched for: ",g,sep="") )
        } else if (length(m) >1 ) {
          print(paste("::Warning:: MULTIPLE DEGs matched and kept for: ",g,sep="") )
          match.ind <- c(match.ind, m)
        } else {
          match.ind <- c(match.ind, m)
        }
      }
      return(match.ind)
    }
    
    #
    # Boxplot function for DEGs
    #
    plot_deg_dir <- function(degs.up, degs.down, outname) {
      cdf <- data.frame("Count"=c(length(degs.up),length(degs.down) ), "Direction"=c("Up", "Down"))
      png(outname, res=plot_res, units="in", height=3, width=5 )
      p <- ggplot(cdf, aes(x=Direction, y=Count, fill=Direction)) + geom_bar(stat="identity") +
        theme_bw(base_size=16) + theme(legend.position="none") + scale_fill_manual(values=list("Up"="red", "Down"="dodgerblue"))
      print(p)
      graphics.off()
    }
    
    
    ### significant pathway dot plot from fgsea table ###
    plot_sig_fGSEA <- function(cur_fgsea, bounds = c(-2,2), outname="unnamed_dotplot.png", ht=8, sigpval = 0.05) {
      sigsel <- which(cur_fgsea[["padj"]]<sigpval)
      sigdf <- data.frame("Pathway"= unlist(lapply(cur_fgsea[["pathway"]][sigsel], function(x) gsub("_", " ", x) )),
                          "qvalue"=cur_fgsea[["padj"]][sigsel], "Count"=cur_fgsea[["size"]][sigsel], "NES"=cur_fgsea[["NES"]][sigsel])
      sigdf$type <- "Upregulated"
      sigdf$type[sigdf$NES < 0 ] <- "Downregulated"
      # format pathways names
      fpath <- unlist(lapply(sigdf$Pathway, function(x) str_to_title( gsub("HALLMARK ","", x) ) ))
      sigdf$Pathway <- fpath
      sigdf <- sigdf %>% arrange(factor(Pathway, levels = sigdf[["Pathway"]][order(sigdf[["qvalue"]], decreasing=F)]))
      #ggplot(data=sigdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=16)
      #orddf <- sigdf %>% mutate(Pathway = fct_reorder(Pathway, qvalue, .desc=T)) #idk what this does...
      plotdf <- sigdf
      #if (dim(plotdf)[1]>20) { plotdf <- plotdf[seq(1,20),]} #limit output if too large
      # #resize if too many pathways
      # ht <- 8
      # if (dim(plotdf)[1]>20 & dim(plotdf)[1]>30) { 
      #   ht <- 12
      # } else {
      #   ht <- 18
      # }
      plotdf$Pathway <- factor(plotdf$Pathway, levels=rev(plotdf$Pathway))
      print(paste("Outputting to: ",outname," with dims: ", paste(dim(plotdf),collapse=" x "),sep=""))
      png(outname, res=plot_res, units="in", height=ht, width=12)
      p <- ggplot(plotdf, aes(x=NES, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=18) + 
        scale_color_gradient(limits=c(0,sigpval), low="red", high="blue") 
        coord_cartesian(xlim=bounds)
      print(p)
      dev.off()
    }
    
    
    ###
    # function to plot a data.frame of CML_contribution type; used to make CML contirbution plots for subsets of genes
    ###
    plot_CML_contribution <- function(plotName, in.df, plot.all=F) {
      t2.df <- in.df[which(in.df$diff.lab=="up" | in.df$diff.lab=="down" ),]
      tup <- sort.int(t2.df$CML_contribution, index.return = T, decreasing = T)$ix[seq(1,10)]
      tdown <- sort.int(t2.df$CML_contribution, index.return = T, decreasing = F)$ix[seq(1,5)]
      sel.pts <- c(tup, tdown)
      if (plot.all) { sel.pts <- seq(1,dim(t2.df)[1])}
      cl.up <- length(which(in.df$CML_contribution>0 & in.df$current_lfc>0))
      cl.down <- length(which(in.df$CML_contribution>0 & in.df$current_lfc<0))
      al.up <- length(which(in.df$CML_contribution<0 & in.df$current_lfc>0))
      al.down <- length(which(in.df$CML_contribution<0 & in.df$current_lfc<0))
      png(plotName, res=plot_res, units="in", height=8, width=10)
      p <- ggplot(in.df, aes(x=CML_contribution, y= PC1, color=diff.lab )) +
        geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", linewidth=1.25) +
        geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
        geom_text(data=t2.df[sel.pts,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene)) +
        ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
      #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" ))
      print(p)
      graphics.off()
    }
    
    
    
    ###
    ### LOADING VALUE CONTRIBUTION FUNCTION
    ###
    # :::REQUIREMENTS:::
    # "out_name" - label for comparison used for output filenames
    # "deg_list" - input gene set for which CML contribution will be calculated
    # "deg_comp" - label of DEG comparison in "gene.df" that will be used for log2FC; this should be a pro-CML comparison (i.e. first "test" term is the more CML-like group [ex. c5 vs controls])
    # :::inputs with default values:::
    # "gene.df" - (default input) object that holds loading values AND log2FC from any DEG comparison being used; alternative object may be specified
    # "enrichOut" - (default: True) should enrichR be performed on pro- and anti-CML genes separately? Set to "False" to save time.
    # "all.gene.labs" - (default: False) should all gene labels be output in plot? Only set to "True" if the gene set is small in size
    gene_set_CML_contribution <- function( out_name, deg_list, deg_comp, 
                                           gene.df = gene.df, enrichOut = T, all.gene.labs=F ) { #all
      
      curdeg <- out_name #deg_name[[i]] #NEEDED:change curpath to curdeg
      cgenes.id <- deg_list
      curcomp <- deg_comp
      
      sel <- match_genes( gene.df$gene_id, cgenes.id)
      c.df <- gene.df[sel,]
      
      m.1 <- c.df[[paste(curcomp,"_mean.",1,sep="")]]
      m.2 <- c.df[[paste(curcomp,"_mean.",2,sep="")]]
      c.diff <- m.1 - m.2
      c.ld <- c.df$eigengenes
      up.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]
      dn.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]
      c.all <- c.diff * c.ld #get contribution for all genes in pathway
      c.up <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
      c.down <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
      ctab <- matrix( c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", sum(c.ld), sum(up.ld), sum(dn.ld), sum(c.all,na.rm=T), sum(c.up), sum(c.down) ), nrow=3, ncol=5) 
      #add lfc contribution; need to figure out why "Inf" is produced
      c.lfc <- log2((m.1+amin)/(m.2+amin) )
      cl.all <- c.ld * c.lfc #get contribution for all genes in pathway based on LFC
      cl.up <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
      cl.down <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
      ctab <- matrix(c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", 
                        length(which(!is.na(c.all))), length(which(!is.na(c.up))), length(which(!is.na(c.down))), # count of genes contributing (non NA genes)
                        sum(c.ld), sum(up.ld), sum(dn.ld), #eigenvalue
                        sum(c.all,na.rm=T), sum(c.up, na.rm=T), sum(c.down, na.rm=T), #contribution w/ exp diff
                        sum(cl.all,na.rm=T), sum(cl.up, na.rm=T), sum(cl.down, na.rm=T) ), # contrib. w/ log2FC
                     nrow=3, ncol=7)
      if (length(which(is.na(ctab))) != 0) {
        print(paste(":::ERROR::: NAs produced for comp: ",curcomp," pathway: ",curdeg,sep=""))
        print(ctab)
        stop("NA encountered")
      }
      cont.tab <- rbind(cont.tab, ctab)
      #plot space
      c.df[["exp.diff"]] <- c.diff
      c.df[["cont.all"]] <- c.all
      dif <- c.df$exp.diff
      # dif[which(c.df$exp.diff>0)] <- "up"
      # dif[which(c.df$exp.diff<0)] <- "down"
      dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]>0)] <- "up"
      dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]<0)] <- "down"
      c.df[["diff.lab"]] <- dif
      #build loading contribution
      "eigen.up" = mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      "eigen.down" <- mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
      "PC1.up" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      "PC1.down" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      
      png(paste(plot_out,"/eigenDEGs_",curdeg,"_comparison-",curcomp,"_eigengene.vs.PC1_pts.png",sep=""), res=plot_res, units="in", height=8, width=10)
      p <- ggplot(c.df, aes(x=eigengenes, y= PC1, color=diff.lab )) + geom_point() + theme_bw() +
        geom_segment(aes(x = 0, y = 0, xend = eigen.down, yend = PC1.down), color="#06398c", #old="#1ebd80"
                     arrow = arrow(length = unit(.65, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = eigen.up, yend = PC1.up), color= "#ab1844", #old= "#c97f10"
                     arrow = arrow(length = unit(.65, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        scale_x_continuous(trans = "reverse") + scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" ))
      print(p)
      graphics.off()
      
      ### make "eigen-contribution"
      # positive value = positive contribution (neg eig + pos exp OR pos eig + neg exp)
      # neg contribution (neg eig + neg exp OR pos eig + pos exp)
      cml.c <- c.df$eigengenes
      # cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff > 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff > 0)] * -1 #pro-CML up (neg eig. -> adjust)
      # cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff < 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff < 0)] * 1 #pro-CML down ( pos eig. -> NO adjust)
      # cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff < 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff < 0)] * 1 # anti-CML down (neg -> NO adjust)
      # cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff > 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff > 0)] * -1 # anti-CML up (pos -> adjust)
      cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 #pro-CML up (neg eig. -> adjust)
      cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 #pro-CML down ( pos eig. -> NO adjust)
      cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 # anti-CML down (neg -> NO adjust)
      cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 # anti-CML up (pos -> adjust)
      
      #sums
      # neg.up <- sum(cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff > 0)], na.rm = F)
      # neg.up.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff > 0 & !is.na(c.df$eigengenes) )] )
      # pos.down <- sum(cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff < 0)], na.rm = F)
      # pos.down.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff < 0 & !is.na(c.df$eigengenes) )] )
      # neg.down <- sum(cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff < 0)], na.rm = F)
      # neg.down.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff < 0 & !is.na(c.df$eigengenes) )] )
      # pos.up <- sum(cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff > 0)], na.rm = F)
      # pos.up.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff > 0 & !is.na(c.df$eigengenes) )] )
      neg.up <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
      neg.up.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
      pos.down <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
      pos.down.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
      neg.down <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
      neg.down.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
      pos.up <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
      pos.up.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
      cmat <- matrix( c( rep(curdeg, 4), rep(curcomp, 4), "neg|up", "pos|down", "neg|down", "pos|up", 
                         neg.up, pos.down, neg.down, pos.up,
                         neg.up.l, pos.down.l, neg.down.l, pos.up.l),
                      nrow=4, ncol=5) 
      cml.tab <- rbind(cml.tab, cmat )
      #cml label
      cml.lab <- rep(NA, length(cml.c))
      cml.lab[which(cml.c > 0)] <- "pro-CML"
      cml.lab[which(cml.c < 0)] <- "anti-CML"
      #full matrix
      # :::note::: the mean contribution of each class of genes must be added!
      #(neg.up  / neg.up.l) + (pos.down / pos.down.l), (neg.down / neg.down.l) + (pos.up / pos.up.l)
      n.u.m <- (neg.up  / neg.up.l) 
      p.d.m <- (pos.down / pos.down.l)
      n.d.m <- (neg.down / neg.down.l) 
      p.u.m <- (pos.up / pos.up.l)
      fmat <- matrix( c( rep(curdeg, 2), rep(curcomp, 2), "pro-CML", "anti-CML",
                         sum(c(n.u.m, p.d.m), na.rm=T), sum(c(n.d.m, p.u.m), na.rm=T),
                         neg.up.l + pos.down.l, neg.down.l + pos.up.l),
                      nrow=2, ncol=5) 
      cml.full <- rbind(cml.full, fmat )
      ### eigengene value vs PC1 ###
      #counts
      # cl.up <- length(which(c.df$eigengenes < 0 & c.df$exp.diff > 0) )
      # cl.down <- length(which(c.df$eigengenes > 0 & c.df$exp.diff < 0) )
      # al.up <- length(which(c.df$eigengenes < 0 & c.df$exp.diff < 0) )
      # al.down <- length(which(c.df$eigengenes > 0 & c.df$exp.diff < 0) )
      cl.up <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
      cl.down <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
      al.up <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
      al.down <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
      c.df[["CML_contribution"]] <- cml.c
      cml.up = mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      cml.down <- mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
      PC1.up <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      PC1.down <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      cml.tot <- mean(cml.c, na.rm=T) 
      PC1.tot <- mean(c.df$PC1, na.rm=T)
      #make df for labeling genes without NA genes
      l.df <- c.df[which(c.df$diff.lab=="up" | c.df$diff.lab=="down" ),]
      tup <- sort.int(l.df$CML_contribution, index.return = T, decreasing = T)$ix[seq(1,10)]
      tdown <- sort.int(l.df$CML_contribution, index.return = T, decreasing = F)$ix[seq(1,5)]
      png(paste(plot_out,"/eigenDEGs_",curdeg,"_comparison-",curcomp,"_pts+label.png",sep=""), res=plot_res, units="in", height=8, width=10)
      p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab )) + 
        geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", size=1.25) +
        geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.5, 
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
        geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene)) +
        geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene)) +
        ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
      #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" )) 
      print(p)
      graphics.off()
      png(paste(plot_out,"/eigenDEGs_",curdeg,"_comparison-",curcomp,"_pts+label+tot.png",sep=""), res=plot_res, units="in", height=8, width=10)
      p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab )) + 
        geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", linewidth=1.25) +
        geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.5, 
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
        geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene)) +
        geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene)) +
        ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
      #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" )) 
      print(p)
      graphics.off()
      #all names
      if (all.gene.labs) {
        png(paste(plot_out,"/eigenDEGs_",curdeg,"_comparison-",curcomp,"_pts+label+tot_allGenes.png",sep=""), res=plot_res, units="in", height=8, width=10)
        p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab )) +
          geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", linewidth=1.25) +
          geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.5,
                       arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
          geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.5,
                       arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
          geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.5,
                       arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
          scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
          geom_text(data=l.df, aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene)) +
          geom_text(data=l.df, aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene)) +
          ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
        #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" ))
        print(p)
        graphics.off()
      }
      
      ###add current test columns
      cols <- grep(paste(curdeg,"_",sep=""), colnames(c.df), value=T) #get all columns labeled with the comparison being processed
      for (c_lab in cols ){ #add identical columns with "current" used as label
        c.df[[gsub(curdeg, "current", c_lab)]] <- c.df[[c_lab]]
      }
      write.table(c.df, paste(plot_out,"/eigenDEGs_",curdeg,"_comparison-",curcomp,"_dataTable.tsv",sep=""), sep="\t", row.names=F) 
      
      
      ### vs log2FC ###
      
      lfc.up <- mean(c.df[[paste(curcomp,"_lfc",sep="")]][which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      lfc.down <- mean(c.df[[paste(curcomp,"_lfc",sep="")]][which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      lfc.tot <- mean(c.df[[paste(curcomp,"_lfc",sep="")]], na.rm=T)
      png(paste(plot_out,"/eigenDEGs_",curdeg,"_vs.log2FC_comparison-",curcomp,"_pts+label+tot.png",sep=""), res=plot_res, units="in", height=8, width=10)
      p <- ggplot(c.df, aes(x=CML_contribution, y= .data[[paste(curcomp,"_lfc",sep="")]], color=diff.lab )) + 
        geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", linewidth=1.25) +
        geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = lfc.down), color="#06398c", size=1.5, 
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = lfc.up), color= "#ab1844", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = lfc.tot), color= "black", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
        geom_text(data=l.df[tup,], aes(x=CML_contribution, y=.data[[paste(curcomp,"_lfc",sep="")]]+0.0007, color = diff.lab, label=gene)) +
        geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=.data[[paste(curcomp,"_lfc",sep="")]]+.0007, color = diff.lab, label=gene)) +
        ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
      #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" )) 
      print(p)
      graphics.off()
      
      if (enrichOut) {
        ### ENRICHR on pro- & anti-CML genes ###
        cur_path_out <- paste("enrichr_CML-cont_",curcomp,sep="")
        dir.create(cur_path_out, showWarnings = F)
        degs <- c.df$gene
        degs.pro <- c.df$gene[which(cml.c > 0)]
        print(paste("Pro-CML genes ",length(degs.pro)))
        degs.anti <- c.df$gene[which(cml.c < 0)]
        print(paste("Pro-CML genes ",length(degs.anti)))
        enriched.up <- enrichr(degs.up, test_dbs)
        enriched.down <- enrichr(degs.down, test_dbs)
        output_enrichr_results(enriched.up, "pro-CML", cur_path_out)
        output_enrichr_results(enriched.down, "anti-CML", cur_path_out)
      }
      
    } #end gene set contribution function
    #set prefix of output; :NOTE: can include directories but they must exist!
    gene_set_CML_contribution_outpre <- function( out.prefix, out_name, deg_list, deg_comp, 
                                                  gene.df = gene.df, enrichOut = T, all.gene.labs=F ) { #all
      curdeg <- out_name #deg_name[[i]] #NEEDED:change curpath to curdeg
      cgenes.id <- deg_list
      curcomp <- deg_comp
      
      sel <- match_genes( gene.df$gene_id, cgenes.id)
      c.df <- gene.df[sel,]
      
      m.1 <- c.df[[paste(curcomp,"_mean.",1,sep="")]]
      m.2 <- c.df[[paste(curcomp,"_mean.",2,sep="")]]
      c.diff <- m.1 - m.2
      c.ld <- c.df$eigengenes
      up.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]
      dn.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]
      c.all <- c.diff * c.ld #get contribution for all genes in pathway
      c.up <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
      c.down <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
      ctab <- matrix( c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", sum(c.ld), sum(up.ld), sum(dn.ld), sum(c.all,na.rm=T), sum(c.up), sum(c.down) ), nrow=3, ncol=5) 
      #add lfc contribution; need to figure out why "Inf" is produced
      c.lfc <- log2((m.1+amin)/(m.2+amin) )
      cl.all <- c.ld * c.lfc #get contribution for all genes in pathway based on LFC
      cl.up <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
      cl.down <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
      ctab <- matrix(c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", 
                        length(which(!is.na(c.all))), length(which(!is.na(c.up))), length(which(!is.na(c.down))), # count of genes contributing (non NA genes)
                        sum(c.ld), sum(up.ld), sum(dn.ld), #eigenvalue
                        sum(c.all,na.rm=T), sum(c.up, na.rm=T), sum(c.down, na.rm=T), #contribution w/ exp diff
                        sum(cl.all,na.rm=T), sum(cl.up, na.rm=T), sum(cl.down, na.rm=T) ), # contrib. w/ log2FC
                     nrow=3, ncol=7)
      if (length(which(is.na(ctab))) != 0) {
        print(paste(":::ERROR::: NAs produced for comp: ",curcomp," pathway: ",curdeg,sep=""))
        print(ctab)
        stop("NA encountered")
      }
      cont.tab <- rbind(cont.tab, ctab)
      #plot space
      c.df[["exp.diff"]] <- c.diff
      c.df[["cont.all"]] <- c.all
      dif <- c.df$exp.diff
      # dif[which(c.df$exp.diff>0)] <- "up"
      # dif[which(c.df$exp.diff<0)] <- "down"
      dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]>0)] <- "up"
      dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]<0)] <- "down"
      c.df[["diff.lab"]] <- dif
      #build loading contribution
      "eigen.up" = mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      "eigen.down" <- mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
      "PC1.up" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      "PC1.down" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      
      png(paste(out.prefix,"_eigenDEGs","_comparison-",curcomp,"_eigengene.vs.PC1_pts.png",sep=""), res=plot_res, units="in", height=8, width=10)
      p <- ggplot(c.df, aes(x=eigengenes, y= PC1, color=diff.lab )) + geom_point() + theme_bw() +
        geom_segment(aes(x = 0, y = 0, xend = eigen.down, yend = PC1.down), color="#06398c", #old="#1ebd80"
                     arrow = arrow(length = unit(.65, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = eigen.up, yend = PC1.up), color= "#ab1844", #old= "#c97f10"
                     arrow = arrow(length = unit(.65, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        scale_x_continuous(trans = "reverse") + scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" ))
      print(p)
      graphics.off()
      
      ### make "eigen-contribution"
      # positive value = positive contribution (neg eig + pos exp OR pos eig + neg exp)
      # neg contribution (neg eig + neg exp OR pos eig + pos exp)
      cml.c <- c.df$eigengenes
      # cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff > 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff > 0)] * -1 #pro-CML up (neg eig. -> adjust)
      # cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff < 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff < 0)] * 1 #pro-CML down ( pos eig. -> NO adjust)
      # cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff < 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff < 0)] * 1 # anti-CML down (neg -> NO adjust)
      # cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff > 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff > 0)] * -1 # anti-CML up (pos -> adjust)
      cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 #pro-CML up (neg eig. -> adjust)
      cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 #pro-CML down ( pos eig. -> NO adjust)
      cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 # anti-CML down (neg -> NO adjust)
      cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 # anti-CML up (pos -> adjust)
      
      #sums
      # neg.up <- sum(cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff > 0)], na.rm = F)
      # neg.up.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff > 0 & !is.na(c.df$eigengenes) )] )
      # pos.down <- sum(cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff < 0)], na.rm = F)
      # pos.down.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff < 0 & !is.na(c.df$eigengenes) )] )
      # neg.down <- sum(cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff < 0)], na.rm = F)
      # neg.down.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df$exp.diff < 0 & !is.na(c.df$eigengenes) )] )
      # pos.up <- sum(cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff > 0)], na.rm = F)
      # pos.up.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df$exp.diff > 0 & !is.na(c.df$eigengenes) )] )
      neg.up <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
      neg.up.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
      pos.down <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
      pos.down.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
      neg.down <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
      neg.down.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
      pos.up <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
      pos.up.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
      cmat <- matrix( c( rep(curdeg, 4), rep(curcomp, 4), "neg|up", "pos|down", "neg|down", "pos|up", 
                         neg.up, pos.down, neg.down, pos.up,
                         neg.up.l, pos.down.l, neg.down.l, pos.up.l),
                      nrow=4, ncol=5) 
      cml.tab <- rbind(cml.tab, cmat )
      #cml label
      cml.lab <- rep(NA, length(cml.c))
      cml.lab[which(cml.c > 0)] <- "pro-CML"
      cml.lab[which(cml.c < 0)] <- "anti-CML"
      #full matrix
      # :::note::: the mean contribution of each class of genes must be added!
      #(neg.up  / neg.up.l) + (pos.down / pos.down.l), (neg.down / neg.down.l) + (pos.up / pos.up.l)
      n.u.m <- (neg.up  / neg.up.l) 
      p.d.m <- (pos.down / pos.down.l)
      n.d.m <- (neg.down / neg.down.l) 
      p.u.m <- (pos.up / pos.up.l)
      fmat <- matrix( c( rep(curdeg, 2), rep(curcomp, 2), "pro-CML", "anti-CML",
                         sum(c(n.u.m, p.d.m), na.rm=T), sum(c(n.d.m, p.u.m), na.rm=T),
                         neg.up.l + pos.down.l, neg.down.l + pos.up.l),
                      nrow=2, ncol=5) 
      cml.full <- rbind(cml.full, fmat )
      ### eigengene value vs PC1 ###
      #counts
      # cl.up <- length(which(c.df$eigengenes < 0 & c.df$exp.diff > 0) )
      # cl.down <- length(which(c.df$eigengenes > 0 & c.df$exp.diff < 0) )
      # al.up <- length(which(c.df$eigengenes < 0 & c.df$exp.diff < 0) )
      # al.down <- length(which(c.df$eigengenes > 0 & c.df$exp.diff < 0) )
      cl.up <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
      cl.down <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
      al.up <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
      al.down <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
      c.df[["CML_contribution"]] <- cml.c
      cml.up = mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      cml.down <- mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
      PC1.up <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      PC1.down <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      cml.tot <- mean(cml.c, na.rm=T) 
      PC1.tot <- mean(c.df$PC1, na.rm=T)
      #make df for labeling genes without NA genes
      l.df <- c.df[which(c.df$diff.lab=="up" | c.df$diff.lab=="down" ),]
      tup <- sort.int(l.df$CML_contribution, index.return = T, decreasing = T)$ix[seq(1,10)]
      tdown <- sort.int(l.df$CML_contribution, index.return = T, decreasing = F)$ix[seq(1,5)]
      png(paste(out.prefix,"_eigenDEGs","_comparison-",curcomp,"_pts+label.png",sep=""), res=plot_res, units="in", height=8, width=10)
      p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab )) + 
        geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", size=1.25) +
        geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.5, 
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
        geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene)) +
        geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene)) +
        ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
      #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" )) 
      print(p)
      graphics.off()
      png(paste(out.prefix,"_eigenDEGs","_comparison-",curcomp,"_pts+label+tot.png",sep=""), res=plot_res, units="in", height=8, width=10)
      p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab )) + 
        geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", linewidth=1.25) +
        geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.5, 
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
        geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene)) +
        geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene)) +
        ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
      #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" )) 
      print(p)
      graphics.off()
      #all names
      if (all.gene.labs) {
        png(paste(out.prefix,"_eigenDEGs","_comparison-",curcomp,"_pts+label+tot_allGenes.png",sep=""), res=plot_res, units="in", height=8, width=10)
        p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab )) +
          geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", linewidth=1.25) +
          geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.5,
                       arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
          geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.5,
                       arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
          geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.5,
                       arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
          scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
          geom_text(data=l.df, aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene)) +
          geom_text(data=l.df, aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene)) +
          ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
        #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" ))
        print(p)
        graphics.off()
      }
      
      ###add current test columns
      cols <- grep(paste(curdeg,"_",sep=""), colnames(c.df), value=T) #get all columns labeled with the comparison being processed
      for (c_lab in cols ){ #add identical columns with "current" used as label
        c.df[[gsub(curdeg, "current", c_lab)]] <- c.df[[c_lab]]
      }
      write.table(c.df, paste(out.prefix,"_eigenDEGs","_comparison-",curcomp,"_dataTable.tsv",sep=""), sep="\t", row.names=F) 
      
      
      ### vs log2FC ###
      
      lfc.up <- mean(c.df[[paste(curcomp,"_lfc",sep="")]][which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      lfc.down <- mean(c.df[[paste(curcomp,"_lfc",sep="")]][which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      lfc.tot <- mean(c.df[[paste(curcomp,"_lfc",sep="")]], na.rm=T)
      png(paste(out.prefix,"_eigenDEGs","_vs.log2FC_comparison-",curcomp,"_pts+label+tot.png",sep=""), res=plot_res, units="in", height=8, width=10)
      p <- ggplot(c.df, aes(x=CML_contribution, y= .data[[paste(curcomp,"_lfc",sep="")]], color=diff.lab )) + 
        geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", linewidth=1.25) +
        geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = lfc.down), color="#06398c", size=1.5, 
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = lfc.up), color= "#ab1844", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = lfc.tot), color= "black", size=1.5,
                     arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
        geom_text(data=l.df[tup,], aes(x=CML_contribution, y=.data[[paste(curcomp,"_lfc",sep="")]]+0.0007, color = diff.lab, label=gene)) +
        geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=.data[[paste(curcomp,"_lfc",sep="")]]+.0007, color = diff.lab, label=gene)) +
        ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
      #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" )) 
      print(p)
      graphics.off()
      
      if (enrichOut) {
        ### ENRICHR on pro- & anti-CML genes ###
        cur_path_out <- paste("enrichr_CML-cont_",curcomp,sep="")
        dir.create(cur_path_out, showWarnings = F)
        degs <- c.df$gene
        degs.pro <- c.df$gene[which(cml.c > 0)]
        print(paste("Pro-CML genes ",length(degs.pro)))
        degs.anti <- c.df$gene[which(cml.c < 0)]
        print(paste("Pro-CML genes ",length(degs.anti)))
        enriched.up <- enrichr(degs.up, test_dbs)
        enriched.down <- enrichr(degs.down, test_dbs)
        output_enrichr_results(enriched.up, "pro-CML", cur_path_out)
        output_enrichr_results(enriched.down, "anti-CML", cur_path_out)
      }
      
    } #end gene set contribution function
    
  
  } #end helper function section
  
  ### Set output directory for plots ###
  plot_out <- "plots"
  dir.create(plot_out, showWarnings = F)
  
  manFig_out <- "manuscript_figures"
  dir.create(manFig_out, showWarnings = F)
  
  ### load data ###
  dat.se <- readRDS("Robj_dat.se_final.rds")
  #:notes on data:
  # CML state-space coordinates are stored in: colData(dat.se)[["CML_space"]] 
  # HSA gene expression stored in "colData" are taken from "abundance"
  
  ### hardcoded critical points values (in case state-space construction is not run)
  # :note: these will be over writen if state-space construction section is run
  ss.cps <- c(-208.51554, -145.15660,  -43.88014, 111.85604, 138.55368)
  
  
  ### plotting parameters ###
  {
    ### set plot resolution
    plot_res = 150
    
    ### set variable ploting scheme (color, shapes, etc) ###
    {

      state_rx.cols <- c("Ctrl"="#000000","c1"="#1b2afa","c2"="#fab71b","c3"="#fc497c","A.post"="#45e9f5","D.Rx"="#6507b3","D.post"="#d955c9")
      trt.cols <- list("C" = "black", "A.post"= "#45e9f5", "D.Rx" = "#6507b3", "D.post"="#d955c9")
      group.cols <- c("A"="#4bd2de", "B"="#eb176c", "C"="#000000","D"="#6507b3")
      group_name.cols <- c("Tet-Off-On"="#10cc81", "CML"="#eb176c", "Control"="#000000","TKI"="#1745eb")
      state_chr.cols <- c("c1"="#1b2afa","c2"="#fab71b","c3" = "#fa1b2e")
      deg.cols <- list("no" = "grey", "yes" =  "#54d676")
      
      ###shapes###
      state_chr.shapes <- c("Ctrl"=21,"c1"=21,"c2"=22,"c3"=23,"A.post"=24,"D.Rx"=22,"D.post"=25)
      
    }
  }

}




###
### CML State-space construction on expression data ###
###
# :NOTE: included for reproducibility
#   CML state-space coordinates are included in: colData(dat.se)[["CML_space"]]
{
  #
  # perform svd on the space; rotate so that Ctrl fit line is has slope = 0; project treatment mice; determine critical points
  #
  
  
  ### perform SVD on CML and Control mice
  {  
    #get index of genes that are zero in all samples; these genes are not included in state-space construction
    lowExpRM <- which(rowSums(data.matrix(assay( dat.se, "abundance")))>0)
    
    minexp <- min(data.matrix(assay(dat.se, "counts"))[which(data.matrix(assay(dat.se, "counts"))>0 )]) # get counts zero offset
    minabd <- min(data.matrix(assay(dat.se, "abundance"))[which(data.matrix(assay(dat.se, "abundance"))>0 )]) # get abundance zero offset
    dat.mc <- scale( t(log2(data.matrix(assay(dat.se, "abundance")) + minabd)), scale=F )
    assays(dat.se)$almc <- t(dat.mc) #ADD abundance log mean-centered (almc) assay to SE object
    dat.alog <-  log2(data.matrix(assay(dat.se, "abundance")) + minabd)
    assays(dat.se)$alog <- dat.alog #ADD abundance log abundance (alog) assay to SE object
    #!!!NOTE !!! use CML ("Group"="B") and Control ("Group"="C") ONLY to construct svd state-space
    ldat <- log2(data.matrix(assay(dat.se, "abundance")) + minabd)
    cbsel <- which( colData(dat.se)$Group=="C" | colData(dat.se)$Group=="B" )
    cb.mc <- rowMeans(ldat[,cbsel])
    cb.almc <- sweep(ldat[lowExpRM,cbsel], 1, cb.mc[lowExpRM], FUN="-") #remove low expression genes
    bsvd <- svd(t(cb.almc))
    bV <- bsvd$v
    bU <- bsvd$u
    bD <- bsvd$d
    bDm <- diag(bD)
    png(paste(plot_out,"/CML_space-unrotated_screeplot.png",sep=""), res=plot_res, units="in", height=3, width=5)
    barplot(bD/sum(bD)*100, col="black", border=NA)
    graphics.off()
    rownames(bV) <- rowData(dat.se)$gene_name[lowExpRM]
    rownames(bU) <- colnames(dat.se)[cbsel]
    #CB only
    sdf <- data.frame("PC1"=bU[,1], "PC2"=bU[,2], "PC3"=bU[,3], "PC4"=bU[,4], "PC5"=bU[,5], "PC6"=bU[,6], "PC7"=bU[,7], 
                      "Treatment"=colData(dat.se)$treatment[cbsel], "Time"=as.numeric(colData(dat.se)$timepoint[cbsel]),
                      "state_chr"=colData(dat.se)$state_chr[cbsel], "Myeloid"=colData(dat.se)$myeloid[cbsel], "myeloid_chr"=colData(dat.se)$myeloid_chr[cbsel],
                      "mouse_id"=colData(dat.se)$mouse_id[cbsel], "Group"=colData(dat.se)$Group[cbsel])
    
    
    ### ROTATE SPACE ###
    { 
      # :notes:
      # CML state-space is rotated so that the linear fit of all control samples has a slope = 0
      # PC2 can be flipped so that CML state is associated with positive PC2 values instead of negative
      
      #if needed, arrange state-space so CML is "down" (i.e., negative PC2 values)
      if (mean(bU[,2][which(sdf$Group=="B")]) > 0) {
        bU[,2] <- bU[,2] * -1
      }

      
      ### process to find slope of controls closest to zero
      {
        test_deg <- seq(18, 22, by=0.001) 
        mlist <- c()
        for (deg in test_deg) {
          theta <- deg * pi / 180
          RotM <- matrix(c(c(cos(theta), sin(theta) ), c(-1*sin(theta), cos(theta)) ), ncol=2, nrow=2 )
          rU <- bU[,c(1,2)] %*% bDm[c(1,2), c(1,2)]  %*% RotM
          r.exp.df <- data.frame("PC1"=rU[,1], "PC2"=rU[,2], 
                                 "Treatment"=colData(dat.se)$treatment[cbsel], "Time"=as.numeric(colData(dat.se)$sample_weeks[cbsel]),
                                 "mouse_id"=colData(dat.se)$mouse_id[cbsel], "Sex"=colData(dat.se)$sex[cbsel])
          r.ctrl.fit <- lm(r.exp.df[which(r.exp.df$Treatment=="TET_ON_C"),]$PC2 ~ r.exp.df[which(r.exp.df$Treatment=="TET_ON_C"),]$PC1 )
          mlist <- c(mlist, r.ctrl.fit$coefficients[[2]]) #record slope of current rotation
        }
        #set theta to when slope of controls is closest to zero
        minM <- which(abs(mlist)==min(abs(mlist)) ) #get smallest slope
        deg <- test_deg[minM]
        theta <- deg * pi / 180
      }
      
      ### rotate state-space by angle "theta"
      RotM <- matrix(c(c(cos(theta), sin(theta) ), c(-1*sin(theta), cos(theta)) ), ncol=2, nrow=2 )
      rU <- bU[,c(1,2)] %*% bDm[c(1,2), c(1,2)]  %*% RotM #scale space by eigenvalues
      rV <- bV[,c(1,2)]  %*% RotM
      rownames(rV) <- rowData(dat.se[lowExpRM,])$gene_name
      r.exp.df <- data.frame("PC1"=rU[,1], "PC2"=rU[,2], 
                             "Treatment"=colData(dat.se)$treatment[cbsel], "Time"=as.numeric(colData(dat.se)$sample_weeks[cbsel]),
                             "mouse_id"=colData(dat.se)$mouse_id[cbsel], "Sex"=colData(dat.se)$sex[cbsel],
                             "group_name"=colData(dat.se)$group_name[cbsel])
      
      r.ctrl.fit <- lm(r.exp.df[which(r.exp.df$Treatment=="TET_ON_C"),]$PC2 ~ r.exp.df[which(r.exp.df$Treatment=="TET_ON_C"),]$PC1 )
      png(paste(plot_out,"/CML_space-FINAL_PC1.vs.PC2_control_fit_line.png",sep=""), res=plot_res, units="in", height=4, width=4)
      plot(r.exp.df$PC1, r.exp.df$PC2, pch=16)
      abline(a=r.ctrl.fit$coefficients[[1]], b=r.ctrl.fit$coefficients[[2]])
      graphics.off()
    }
    
    #plot final rotated state-space
    {
      png(paste(plot_out,"/CML_space-FINAL_PC1.vs.PC2_BC-only.png",sep=""), res=plot_res, units="in", height=6, width=8)
      ggplot(r.exp.df, aes(x=PC1, y=PC2, color=group_name)) + geom_point(size=2) + theme_bw(base_size=16) +
        scale_color_manual(values=group_name.cols)
      graphics.off()
      png(paste(plot_out,"/CML_space-FINAL_time.vs.PC2_BC-only.png",sep=""), res=plot_res, units="in", height=6, width=8)
      ggplot(r.exp.df, aes(x=Time, y=PC2, color=as.character(group_name), group=mouse_id)) + geom_point() + geom_line() + theme_bw(base_size=16)
      graphics.off()
    }
  }
  
  
  ### project treatment mice into space ###
  {
    #project TKI ("Group"="D") and TOTO ("Group"="A") mice
    adsel <- which( colData(dat.se)$Group=="A" | colData(dat.se)$Group=="D" )
    ad.almc <- sweep(ldat[, adsel], 1, cb.mc, FUN="-") #mean-center using means of samples used to construct the state-space (i.e., CML + Control mice)
    rU.treat <- t(ad.almc[lowExpRM, ])  %*% rV[,c(1,2)] # get state-space corrodinates for treatment mice using expression of genes included in state-space

    ### order samples to match column order of "dat.se"
    cbcol <- colData(dat.se)[cbsel,] #colData that matches "rU"
    adcol <- colData(dat.se)[adsel,] #colData that matches "rU.treat"
    cml_space <- c()
    samp_test <-c()
    for (i in 1:dim(colData(dat.se))[1]) {
      #get sample
      cursamp <- colData(dat.se)$sample[i]
      samp_test <- c(samp_test, cursamp)
      #match to correct group; CB or AD
      if (colData(dat.se)$Group[i]=="C" | colData(dat.se)$Group[i]=="B" ) {
        m <- which(cbcol$sample==cursamp)
        if (length(m)!=1) {
          print(paste(":ERROR: matching sample ",cursamp,sep=""))
        } else {
          cml_space <- c(cml_space, rU[m,2])
        }
      } else if (colData(dat.se)$Group[i]=="A" | colData(dat.se)$Group[i]=="D" ) {
        m <- which(adcol$sample==cursamp)
        if (length(m)!=1) {
          print(paste(":ERROR: matching sample ",cursamp,sep=""))
        } else {
          cml_space <- c(cml_space, rU.treat[m,2])
        }
      }else {
        print(":::ERROR:::")
      }
    } #end cml_state construction
    
    # state-space coordinates calculated above stored as: "CML_space.calc"
    # state-space coordinates used in manuscript strored as: "CML_space.manuscript"
    tdf <- data.frame("CML_space.calc"=cml_space, "Time"=as.numeric(colData(dat.se)$sample_weeks),
                      "group_name"=colData(dat.se)$group_name, "mouse_id"=colData(dat.se)$mouse_id,
                      "CML_space.manuscript"=colData(dat.se)[["CML_space"]] ) 
    #plot calculated state-space
    png(paste(plot_out,"/CML_space-FINAL_calculated_PC1.vs.PC2_all_mouse_trajectories.png",sep=""), res=plot_res, units="in", height=3, width=5)
    ggplot(tdf, aes(x=Time, y=CML_space.calc, color=group_name, group=mouse_id))  + geom_path(alpha=.5) + geom_point(size=2) + theme_bw(base_size=16) +
      scale_color_manual(values=group_name.cols)
    graphics.off()
    #plot manuscript state-space
    png(paste(plot_out,"/CML_space-FINAL_manuscript_PC1.vs.PC2_all_mouse_trajectories.png",sep=""), res=plot_res, units="in", height=3, width=5)
    ggplot(tdf, aes(x=Time, y=CML_space.manuscript, color=group_name, group=mouse_id))  + geom_path(alpha=.5) + geom_point(size=2) + theme_bw(base_size=16) +
      scale_color_manual(values=group_name.cols)
    graphics.off()
    
  }
  
  
  ### check that the above process results in the state-space coordinates ###
  # CML state-space coordinates are stored in dat.se object as: colData(dat.se)[["CML_space"]]
  # :note: difference in matrix multiplication results are not identical
  #   check that difference are insignificant
  any( abs(tdf[["CML_space.calc"]] - tdf[["CML_space.manuscript"]]) > 10^-8 )
  
  
  
  #state-space plots
  {
    col.df <- data.frame(colData(dat.se))
    
    png(paste(plot_out,"/CML_space-FINAL-rot_PC1.vs.PC2_BC-only.png",sep=""), res=plot_res, units="in", height=6, width=8)
    ggplot(r.exp.df, aes(x=PC1, y=PC2, color=group_name)) + geom_point(size=2) + theme_bw(base_size=16) +
      scale_color_manual(values=group_name.cols)
    graphics.off()
    png(paste(plot_out,"/CML_space-FINAL-rot_time.vs.PC2_BC-only.png",sep=""), res=plot_res, units="in", height=6, width=8)
    ggplot(r.exp.df, aes(x=Time, y=PC2, color=as.character(group_name), group=mouse_id)) + geom_point() + geom_line() + theme_bw(base_size=16)
    graphics.off()
    
    png(paste(plot_out,"/CML_space-FINAL-rot_Group-A_time.vs.CML_space_line+point.png",sep=""), res=plot_res, units="in", height=6, width=10)
    p <- ggplot(col.df, aes(x=sample_weeks, y=CML_space )) + geom_point(color="grey", alpha=.5) +
      geom_line(data=col.df[which(col.df$Group=="A"),], aes(x=sample_weeks, y=CML_space, color=group_name, group=mouse_id)) +
      geom_point(data=col.df[which(col.df$Group=="A"),], aes(x=sample_weeks, y=CML_space, color=group_name, group=mouse_id)) +
      theme_bw(base_size=18)
    print(p)
    graphics.off()
    
    png(paste(plot_out,"/CML_space-FINAL-rot_Group-D_time.vs.CML_space_line+point.png",sep=""), res=plot_res, units="in", height=6, width=10)
    p <- ggplot(col.df, aes(x=sample_weeks, y=CML_space )) + geom_point(color="grey", alpha=.5) +
      geom_line(data=col.df[which(col.df$Group=="D"),], aes(x=sample_weeks, y=CML_space, color=group_name, group=mouse_id)) +
      geom_point(data=col.df[which(col.df$Group=="D"),], aes(x=sample_weeks, y=CML_space, color=group_name, group=mouse_id)) +
      theme_bw(base_size=18)
    print(p)
    graphics.off()
    
    
  }
  
  
  
  ### DEFINE CRITICAL POINTS ###
  ### kernel density to get critical points zeros ###
  {
    png(paste(plot_out,"/CML_space_CML-mouse_density.png",sep=""), res=plot_res, units="in", height=3, width=5)
    ggplot(data = data.frame(colData(dat.se)[which(colData(dat.se)$Group=="B"),]), aes(x = CML_space)) +
      geom_histogram(aes(x=CML_space, y=after_stat(density)), bins=20,  fill="grey", color="grey") + 
      geom_density(adjust=.75) +theme_bw()
    graphics.off()

    #get density function
    dens <- density(data.frame(colData(dat.se)[which(colData(dat.se)$Group=="B"),])$CML_space, adjust=.75, 
                    bw="nrd0", kernel="gaussian")
    # dy/dx first derivative
    first<-diff(dens$y)/diff(dens$x)
    # Second derivative
    second<-diff(first)/diff(dens$x[1:511])
    # Condition for inflection point
    flections<-c()
    cps <- c()
    for(i in 2:length(second)){
      if(sign(first[i])!=sign(first[i-1])){
        flections<-c(flections,i)
      }
      if (first[i]==0) {
        cps <- c(cps, i)
      }
    }
    png(paste(plot_out,"/CML_space_CML-mouse_density+zeros.png",sep=""), res=plot_res, units="in", height=3, width=5)
    plot(dens, main="CML mouse density")
    abline(v=dens$x[flections], lty=2, col="red")
    graphics.off()
    
    
    ### define CML critical points ###
    ss.cps <- dens$x[flections]
    # :NOTE: critical points found by fitting a polynomial to the critical points; 
    #   polynomial could only be fit by removing one of the critical points (c4)
    #   for details see manuscript and Matlab code:
    #   https://github.com/cohmathonc/CML_mRNA_state-transition/Matlab/XXX
    ss.cps[2] <- -145.1566 #manually replace c4 with fitted value from Matlab
    
    png(paste(plot_out,"/CML_space_Group-B_CML_space_dens+hist.png",sep=""), res=plot_res, units="in", height=3, width=5)
    ggplot(data = data.frame(colData(dat.se)[which(colData(dat.se)$Group=="B"),]), aes(x = CML_space)) +
      geom_histogram(aes(x=CML_space, y=..density..), bins=20, fill="grey", color="grey") + 
      geom_density(adjust=.75, color="black") + theme_bw() +
      geom_vline(xintercept = dens$x[flections], color="pink", linetype="dashed") +
      geom_vline(xintercept = ss.cps, color="red", linetype="dashed") 
    graphics.off()
    
    
    
    ### control density ###
    png(paste(plot_out,"/CML_space_Control-mouse_density.png",sep=""), res=plot_res, units="in", height=3, width=5)
    ggplot(data = data.frame(colData(dat.se)[which(colData(dat.se)$Group=="C"),]), aes(x=CML_space, fill=Group)) +
      geom_histogram(aes(x=CML_space, y=..density..),  fill="grey", color="grey") +
      geom_density(adjust=2, color="red", alpha=.5)
    graphics.off()
    c.dens <- density(data.frame(colData(dat.se)[which(colData(dat.se)$Group=="C"),])$CML_space, adjust=2, 
                      bw="nrd0", kernel="gaussian")
    # dy/dx first derivative
    first<-diff(c.dens$y)/diff(c.dens$x)
    # Second derivative
    second<-diff(first)/diff(c.dens$x[1:511])
    # Condition for inflection point
    ctl.flections<-c()
    ctl.cps <- c()
    for(i in 2:length(second)){
      if(sign(first[i])!=sign(first[i-1])){
        ctl.flections<-c(ctl.flections,i)
      }
      if (first[i]==0) {
        ctl.cps <- c(ctl.cps, i)
      }
    }
    png(paste(plot_out,"/CML_space_Control-mouse_density+zeros.png",sep=""), res=plot_res, units="in", height=3, width=5)
    plot(c.dens, xlab="CML_space", main="Controls mouse density")
    abline(v=c.dens$x[ctl.flections], lty=2, col="red")
    graphics.off()
    c.dens$x[ctl.flections]
    c.dens$y[ctl.flections]
    
  }
    
    
  
  
  ### MSD - mean square displacement analysis ###
  {
    
    ### rescale state-space from [c1,c5] to [0,1]
    traj.list <- list()
    time_limit <- 19
    for (g in unique(colData(dat.se)$Group)) {
      cur.dat <- c()
      for (m in unique(colData(dat.se)$mouse_id)) {
        ctime <- as.numeric(colData(dat.se)$sample_weeks[which(colData(dat.se)$Group==g & colData(dat.se)$mouse_id==m)])
        cspace <- colData(dat.se)$CML_space[which(colData(dat.se)$Group==g & colData(dat.se)$mouse_id==m)]
        clab <- rep(m, length(ctime))
        cgroup <- rep(g, length(ctime))
        #scale data
        a <- ss.cps[1]
        b <- ss.cps[5]
        c <- 0
        d <- 1
        s.space <- (cspace - a) * ((d - c) / (b - a)) + c
      
        cdf <- data.frame("x"=as.numeric(ctime), "y"=s.space, "mouse_id"=clab)
        
        #add data
        cur.dat <- rbind(cur.dat, cdf)
      }    
      #process data 
      cur.dat[["mouse"]] <- paste("mid_",cur.dat[["mouse_id"]],sep="")
      #cur.dat[["label"]] <- NULL
      cur.dat[["time"]] <- cur.dat$x + 1 #make 1-indexed
      traj.list[[g]] <- cur.dat
    }
    png(paste(plot_out,"/CML_space_Group-B_time.vs.scaled_CML_space.png",sep=""), res=plot_res, units="in", height=3, width=5)
    ggplot(traj.list[["B"]], aes(x=x, y=y,color=mouse, group=mouse)) + geom_line()
    graphics.off()
    
    
    ### MSD analysis using time step = 1 week
    #:note: TOTO mice (Group==A) was not samples weekly after treatment
    msd.dat <- c()
    for (grp in unique(colData(dat.se)$Group)) {
      df <- traj.list[[grp]]  
      for (m in unique(df$mouse_id)) {
        pdat <- df[which(df$mouse_id==m),]
        max_tp <- max(pdat$time)
        if (grp=="A") { max_tp <- 6}
        #:note: if time is ordered for each sample, then "selts" is the same as the time AND indexes the y values
        if ( !identical(sort(pdat$time),pdat$time) ) {
          print(paste(":::ERROR::: skipping mouse ",m, " where time points are not sorted!"))
          next
        }
        #handle missing time points
        npdat <- c()
        for (tps in 1:max_tp) {
          mch <- which(pdat$time==tps)
          if (length(mch)>0) {
            npdat <- rbind(npdat, pdat[mch,])
          } else {
            npdat <- rbind(npdat, c(tps-1, NA, m, pdat$mouse[1], tps))
          }
        }
        colnames(npdat) <- colnames(pdat)
        pdat <- data.frame(npdat)
        pdat <- pdat[order(as.numeric(pdat$time)),]
        #calculate MSD for each time step (dt)
        for ( dt in 1:(max_tp-1)) { # remove last time point because there can't be an interval of that
          ts <- dt + 1 #account for 1-indexed time
          #calculate square displacements for each possible path
          sq_disp <- c()
          selt_i <- which(pdat$time==1)
          selt_f <- which(pdat$time==ts)

          dx <- rep(dt, length(selt_i) ) #x is the time which is defined as the timestep (dt); divid by time_limit to scale time from 0-1
          #:note: if time is ordered for each sample, then "selts" is the same as the time AND indexes the y values
          y_i <- as.numeric(pdat$y[selt_i])
          y_f <- as.numeric(pdat$y[selt_f])
          dy <- y_f - y_i
          #sq_disp <- c(sq_disp, dx^2 + dy^2 ) #add displacement
          sq_disp <- c(sq_disp, dy^2 )
          msd.dat <- rbind(msd.dat, c(grp, m, dt, mean(sq_disp, na.rm=T)))
        } #end time step
      } #end mouse
    }# end group
    #create df
    colnames(msd.dat) <- c("group", "mouse_id", "time_interval", "MSD")
    msd.df <- data.frame(msd.dat)  
    msd.df[["mouse_id"]] <- as.character(msd.df[["mouse_id"]])
    msd.df[["MSD"]] <- as.numeric(msd.df[["MSD"]])
    msd.df[["time_interval"]] <- as.numeric(msd.df[["time_interval"]])
    
    ### fit the MSD for each mouse
    fit.dat <- c()
    for (m in unique(msd.df$mouse_id)) {
      cdat <- msd.df[which(msd.df$mouse_id==m),]
      fit <- lm(cdat$MSD ~ cdat$time_interval)
      sd <- stats::sigma(fit)
      cm <- summary(fit)$coefficients[[2]]
      cb <- summary(fit)$coefficients[[1]]
      fit.dat <- rbind(fit.dat, c(cdat$group[1], m, cm, cb, sd))
    }
    colnames(fit.dat) <- c("group", "mouse", "slope", "intercept", "sd")
    write.table(fit.dat, paste(plot_out,"/CML_space_MSD_linearFitTable.tsv",sep=""), sep="\t", row.names=F )
    fit.df <- data.frame(fit.dat)
    mean.fits <- c()
    lm.mean <- c()
    for (grp in unique(msd.df$group)) {
      #get mean fit
      grp.sel <- which(fit.df$group==grp)
      mean.fits <- rbind(mean.fits, c(grp, mean(as.numeric(fit.df$slope)[grp.sel]), mean(as.numeric(fit.df$intercept)[grp.sel])) )
      
      #fit group
      cdat <- msd.df[which(msd.df$group==grp),]
      fit <- lm(cdat$MSD ~ cdat$time_interval)
      sd <- stats::sigma(fit)
      cm <- summary(fit)$coefficients[[2]]
      cb <- summary(fit)$coefficients[[1]]
      lm.mean <- rbind(lm.mean, c(grp, cm, cb, sd))
      
      #plots
      png(paste(plot_out,"/CML_space_MSD_grp-",grp,"_linearFit.png",sep=""), res=plot_res, units="in", height=3, width=5)
      p <- ggplot(msd.df[which(msd.df$group==grp),], aes(x = time_interval, y = MSD, color = mouse_id)) + geom_point(alpha=.5)+
         geom_smooth(method="lm", linetype="solid", linewidth=2, se=F) +
        labs(title = "Mean Squared Displacement",
             x = "Time Lag", y = "MSD") +
        theme_minimal()
      print(p)
      graphics.off()
      
      png(paste(plot_out,"/CML_space_MSD_grp-",grp,"_MSD_point+line.png",sep=""), res=plot_res, units="in", height=3, width=5)
      p <- ggplot(msd.df[which(msd.df$group==grp),], aes(x = time_interval, y = MSD, color = mouse_id)) + geom_point()+
        geom_line() +
        labs(title = "Mean Squared Displacement",
             x = "Time Lag", y = "MSD") + 
        theme_minimal()
      print(p)
      graphics.off()
    }
    
    write.table(fit.dat, paste(plot_out,"/CML_space_MSD_linearFitTable.tsv",sep=""), sep="\t", row.names=F )
    
    #all fits
    colnames(mean.fits) <- c("group", "slope", "intercept")
    mean.df <- data.frame(mean.fits)
    mean.df$slope <- as.numeric(mean.df$slope)
    png(paste(plot_out,"/CML_space_MSD_all-Grps_linearFit.png",sep=""), res=plot_res, units="in", height=6, width=6)
    ggplot(msd.df, aes(x = time_interval, y = MSD, color = group)) + geom_point(alpha=.5)+
      geom_smooth(method="lm", linetype="solid", linewidth=2, se=F) +
      labs(title = paste("CML = ",round(mean.df$slope[which(mean.df$group=="B")],digits=4),"\n",
                         "Ctl = ",round(mean.df$slope[which(mean.df$group=="C")],digits=4),"\n",
                         "TOTO  ",round(mean.df$slope[which(mean.df$group=="A")],digits=4),"\n",
                         "TKI = ",round(mean.df$slope[which(mean.df$group=="D")],digits=4),"\n",sep=""),
           x = "Time", y = "MSD") +
      theme_minimal()
    graphics.off()
    
    
    ### ZZZ potentially not needed 
    {  
      ### all CML + treatment pre-week 6 ###
      {
        msdt6.dat <- c()
        for (grp in c("B","D","A")) {
          df <- traj.list[[grp]]  
          for (m in unique(df$mouse_id)) {
            pdat <- df[which(df$mouse_id==m),]
            max_tp <- max(pdat$time)
            if (grp=="A" | grp == "D") { max_tp <- 6}
            #:note: if time is ordered for each sample, then "selts" is the same as the time AND indexes the y values
            if ( !identical(sort(pdat$time),pdat$time) ) {
              print(paste(":::ERROR::: skipping mouse ",m, " where time points are not sorted!"))
              next
            }
            #handle missing time points
            npdat <- c()
            for (tps in 1:max_tp) {
              mch <- which(pdat$time==tps)
              if (length(mch)>0) {
                npdat <- rbind(npdat, pdat[mch,])
              } else {
                npdat <- rbind(npdat, c(tps-1, NA, m, pdat$mouse[1], tps))
              }
            }
            colnames(npdat) <- colnames(pdat)
            pdat <- data.frame(npdat)
            pdat <- pdat[order(as.numeric(pdat$time)),]
            #calculate msdt6 for each time step (dt)
            for ( dt in 1:(max_tp-1)) { # remove last time point because there can't be an interval of that
              ts <- dt + 1 #account for 1-indexed time
              #calculate square displacements for each possible path
              sq_disp <- c()
              selt_i <- which(pdat$time==1)
              selt_f <- which(pdat$time==ts)
              
              dx <- rep(dt, length(selt_i) ) #x is the time which is defined as the timestep (dt); divid by time_limit to scale time from 0-1
              #:note: if time is ordered for each sample, then "selts" is the same as the time AND indexes the y values
              y_i <- as.numeric(pdat$y[selt_i])
              y_f <- as.numeric(pdat$y[selt_f])
              dy <- y_f - y_i
              #sq_disp <- c(sq_disp, dx^2 + dy^2 ) #add displacement
              sq_disp <- c(sq_disp, dy^2 )
              msdt6.dat <- rbind(msdt6.dat, c(grp, m, dt, mean(sq_disp, na.rm=T)))
            } #end time step
          } #end mouse
        }# end group
        #create df
        colnames(msdt6.dat) <- c("group", "mouse_id", "time_interval", "msdt6")
        msdt6.df <- data.frame(msdt6.dat)  
        msdt6.df[["mouse_id"]] <- as.character(msdt6.df[["mouse_id"]])
        msdt6.df[["msdt6"]] <- as.numeric(msdt6.df[["msdt6"]])
        msdt6.df[["time_interval"]] <- as.numeric(msdt6.df[["time_interval"]])
        
        #fits
        fit.dat <- c()
        for (m in unique(msdt6.df$mouse_id)) {
          cdat <- msdt6.df[which(msdt6.df$mouse_id==m),]
          fit <- lm(cdat$msdt6 ~ cdat$time_interval)
          cm <- summary(fit)$coefficients[[2]]
          cb <- summary(fit)$coefficients[[1]]
          fit.dat <- rbind(fit.dat, c(cdat$group[1], m, cm, cb))
        }
        colnames(fit.dat) <- c("group", "mouse", "slope", "intercept")
        #write.table(fit.dat, paste(plot_out,"/CML_space_msdt6_linearFitTable.tsv",sep=""), sep="\t", row.names=F )
        fit.df <- data.frame(fit.dat)
        mean.fits <- c()
        for (grp in unique(msdt6.df$group)) {
          #get mean fit
          grp.sel <- which(fit.df$group==grp)
          mean.fits <- rbind(mean.fits, c(grp, mean(as.numeric(fit.df$slope)[grp.sel]), mean(as.numeric(fit.df$intercept)[grp.sel])) )
        }
      } #end CML + week 6
      
      ### TKI Rx ###
      {
        msd.tki.dat <- c()
        df <- traj.list[["D"]]  
        for (grp in c("D-Rx", "D-postRx")) {
          
          for (m in unique(df$mouse_id)) {
            pdat <- df[which(df$mouse_id==m),]
            
            if (grp=="D-Rx") { 
              pdat <- pdat[which(pdat$time>=6 & pdat$time<10 ),]
              max_tp <- max(pdat$time)  
              min_tp <- min(pdat$time)
            } else {
              pdat <- pdat[which( pdat$time>=10 ),]
              max_tp <- max(pdat$time)
              min_tp <- min(pdat$time)
            }
            #:note: if time is ordered for each sample, then "selts" is the same as the time AND indexes the y values
            if ( !identical(sort(pdat$time),pdat$time) ) {
              print(paste(":::ERROR::: skipping mouse ",m, " where time points are not sorted!"))
              next
            }
            #handle missing time points
            npdat <- c()
            for (tps in min_tp:max_tp) {
              mch <- which(pdat$time==tps)
              if (length(mch)>0) {
                npdat <- rbind(npdat, pdat[mch,])
              } else {
                npdat <- rbind(npdat, c(tps-1, NA, m, pdat$mouse[1], tps))
              }
            }
            colnames(npdat) <- colnames(pdat)
            pdat <- data.frame(npdat)
            pdat <- pdat[order(as.numeric(pdat$time)),]
            #calculate msd.tki for each time step (dt)
            for ( dt in min_tp:(max_tp-1)) { # remove last time point because there can't be an interval of that
              ts <- dt + 1 #account for 1-indexed time
              #calculate square displacements for each possible path
              sq_disp <- c()
              selt_i <- which(pdat$time==min_tp)
              selt_f <- which(pdat$time==ts)
              
              dx <- rep(dt, length(selt_i) ) #x is the time which is defined as the timestep (dt); divid by time_limit to scale time from 0-1
              #:note: if time is ordered for each sample, then "selts" is the same as the time AND indexes the y values
              y_i <- as.numeric(pdat$y[selt_i])
              y_f <- as.numeric(pdat$y[selt_f])
              dy <- y_f - y_i
              #sq_disp <- c(sq_disp, dx^2 + dy^2 ) #add displacement
              sq_disp <- c(sq_disp, dy^2 )
              msd.tki.dat <- rbind(msd.tki.dat, c(grp, m, dt, mean(sq_disp, na.rm=T)))
            } #end time step
          } #end mouse
        }# end group
        #create df
        colnames(msd.tki.dat) <- c("group", "mouse_id", "time_interval", "msd.tki")
        msd.tki.df <- data.frame(msd.tki.dat)  
        msd.tki.df[["mouse_id"]] <- as.character(msd.tki.df[["mouse_id"]])
        msd.tki.df[["msd.tki"]] <- as.numeric(msd.tki.df[["msd.tki"]])
        msd.tki.df[["time_interval"]] <- as.numeric(msd.tki.df[["time_interval"]])
        
        #fits
        fit.dat <- c()
        for (grp in c("D-Rx", "D-postRx")) {
          for (m in unique(msd.tki.df$mouse_id)) {
            cdat <- msd.tki.df[which(msd.tki.df$mouse_id==m & msd.tki.df$group==grp),]
            fit <- lm(cdat$msd.tki ~ cdat$time_interval)
            cm <- summary(fit)$coefficients[[2]]
            cb <- summary(fit)$coefficients[[1]]
            fit.dat <- rbind(fit.dat, c(grp, m, cm, cb))
          }
        }
        colnames(fit.dat) <- c("group", "mouse", "slope", "intercept")
        write.table(fit.dat, paste(plot_out,"/CML_space_msd.tki_linearFitTable.tsv",sep=""), sep="\t", row.names=F )
        fit.df <- data.frame(fit.dat)
        mean.fits <- c()
        for (grp in c("D-Rx", "D-postRx")) {
          #get mean fit
          grp.sel <- which(fit.df$group==grp)
          mean.fits <- rbind(mean.fits, c(grp, mean(as.numeric(fit.df$slope)[grp.sel], na.rm=T), mean(as.numeric(fit.df$intercept)[grp.sel])) )
          
          # #plots
          png(paste(plot_out,"/CML_space_msd.TKI-TRT_grp-",grp,"_linearFit.png",sep=""), res=plot_res, units="in", height=4, width=6)
          p <- ggplot(msd.tki.df[which(msd.tki.df$group==grp),], aes(x = time_interval, y = msd.tki, color = mouse_id)) + geom_point(alpha=.5)+
            geom_smooth(method="lm", linetype="solid", linewidth=2, se=F) +
            labs(title = "Mean Squared Displacement",
                 x = "Time Lag", y = "msd.tki") +
            theme_minimal()
          print(p)
          graphics.off()
  
          png(paste(plot_out,"/CML_space_msd.TKI-TRT_grp-",grp,"_linearFit.png",sep=""), res=plot_res, units="in", height=4, width=6)
          p <- ggplot(msd.tki.df[which(msd.tki.df$group==grp),], aes(x = time_interval, y = msd.tki, color = mouse_id)) + geom_point()+
            geom_line() +
            labs(title = "Mean Squared Displacement",
                 x = "Time Lag", y = "msd.tki") +
            theme_minimal()
          print(p)
          graphics.off()
        }
        write.table(fit.dat, paste(plot_out,"/CML_space_msd.tki_linearFitTable.tsv",sep=""), sep="\t", row.names=F )
        
        #all fits
        colnames(mean.fits) <- c("group", "slope", "intercept")
        mean.df <- data.frame(mean.fits)
        mean.df$slope <- as.numeric(mean.df$slope)
        
      }
    }
    
  } #end MSD


  
    
} #end state-space construction


###
### SURVIVAL ANALYSIS & STATE-TRANSITION MODELING
### :note: See "Matlab" directory all survival analysis performed in Matlab
### https://github.com/cohmathonc/CML_mRNA_state-transition/Matlab/
###


###
### Figure 1-2 plots (Fig. S1-2)
###
{
  #setup for plotting
  hg <- "HSA_BCR_ABL1_gene"
  col.df <- data.frame(colData(dat.se))
  b.sel <- which(colData(dat.se)$Group=="B")
  c.sel <- which(colData(dat.se)$Group=="C")
  
  
  ### 
  # Figure 1A
  # created in biorender
  ###

  
  ### 
  # Figure 1B
  ###
  {
  # BCR-ABL vs CML - Fig 1B
  png(paste(manFig_out,"/Fig-1A_BCR-ABL.vs.CML_B-ONLY_space_criticalPoints_line+point-MANUSCRIPT.png",sep=""), res=plot_res, units="in", height=4, width=3)
  ggplot( col.df[b.sel,], aes(x=HSA_BCR_ABL1_gene, y=CML_space, color=Group)) + geom_point() +
     geom_line(data=col.df[b.sel,], aes(x=HSA_BCR_ABL1_gene, y=CML_space, group=mouse_id)) + 
    geom_smooth( method="lm", se=F, col="dodgerblue", linewidth=2) + scale_color_manual(values=unlist(group.cols)) +
    #geom_smooth( method="lm", se=F, col="grey20", linetype="dashed", linewidth=2) 
    scale_y_continuous(limits = c(-320, 200)) +
    theme_bw(base_size=18) + theme(legend.position="none") 
  graphics.off()
  # CML _space
  png(paste(manFig_out,"/Fig-1A_time.vs.CML_space_BC-ONLY_space_criticalPoints_line+point-MANUSCRIPT.png",sep=""), res=plot_res, units="in", height=4, width=3)
  ggplot( col.df[c(b.sel, c.sel),], aes(x=sample_weeks, y=CML_space, color=Group, group=mouse_id)) + geom_point() +
    geom_hline(yintercept=ss.cps[c(2,4)], color="dimgrey", linetype="dashed") +
    geom_line() + theme_bw(base_size=18) + theme(legend.position="none") + scale_color_manual(values=unlist(group.cols))
  graphics.off()
  

  }
  
  
  ###
  #  Table S1: BCR-ABL correlation with PCs
  ###
  {
    outtab <- c()
    #rotated
    for (pc in 1:dim(rU)[2]) {
      curpc <- rU[,pc]
      ba.fit <- lm(curpc[b.sel] ~ col.df$HSA_BCR_ABL1_gene[b.sel], na.action = na.omit)
      rsq <- summary(ba.fit)$adj.r.squared
      pv <- summary(ba.fit)$coefficients[,4][2]
      #log
      lba.fit <- lm(curpc[b.sel] ~ log(col.df$HSA_BCR_ABL1_gene[b.sel]+minabd), na.action = na.omit)
      lrsq <- summary(lba.fit)$adj.r.squared
      lpv <- summary(lba.fit)$coefficients[,4][2]
      #myeloid
      mba.fit <- lm(curpc[b.sel] ~ col.df$myeloid[b.sel], na.action = na.omit)
      mrsq <- summary(mba.fit)$adj.r.squared
      mpv <- summary(mba.fit)$coefficients[,4][2]
      #output
      outtab <- rbind(outtab, c(paste("rotated-",pc,sep=""), rsq, pv, lrsq, lpv, mrsq, mpv))
    }
    #original
    for (pc in 1:dim(bU)[2]) {
      curpc <- bU[,pc]
      ba.fit <- lm(curpc[b.sel] ~ col.df$HSA_BCR_ABL1_gene[b.sel], na.action = na.omit)
      rsq <- summary(ba.fit)$adj.r.squared
      pv <- summary(ba.fit)$coefficients[,4][2]
      #log
      lba.fit <- lm(curpc[b.sel] ~ log(col.df$HSA_BCR_ABL1_gene[b.sel]+minabd), na.action = na.omit)
      lrsq <- summary(lba.fit)$adj.r.squared
      lpv <- summary(lba.fit)$coefficients[,4][2]
      #myeloid
      mba.fit <- lm(curpc[b.sel] ~ col.df$myeloid[b.sel], na.action = na.omit)
      mrsq <- summary(mba.fit)$adj.r.squared
      mpv <- summary(mba.fit)$coefficients[,4][2]
      #output
      outtab <- rbind(outtab, c(paste("orig-",pc,sep=""), rsq, pv, lrsq, lpv, mrsq, mpv))
    }
    colnames(outtab) <- c("PC", "Rsquared" ,"pvalue", "log-Rsquared", "log-pvalue", "myeloid-Rsquared", "myeloid-pvalue")
    write.table(outtab, paste(manFig_out,"/Tab-S1_BCR-ABL.vs.CML_space_table.tsv",sep=""), sep="\t", row.names=F)
  }
  
  ###
  #  Fig. 1C
  ###
  {
    png(paste(manFig_out,"/Fig-1C_HSA-gene-",hg,"_group-B-ONLY_timepoint-weeks-0-1_boxplot_MANUSCRIPT.png",sep=""), res=plot_res, units="in", height=2, width=3)
    ggplot(data.frame(colData(dat.se)[which(colData(dat.se)$Group=="B" & colData(dat.se)$sample_weeks<=1),]), aes(x=as.character(sample_weeks), y=.data[[hg]], group=sample_weeks)) + 
      geom_boxplot(fill="dimgrey") + theme_bw(base_size=16) + theme(legend.position="none") + scale_y_continuous(limits=c(0,1))
    graphics.off() 
    
    png(paste(manFig_out,"/Fig-1C_CML-space_group-B-ONLY_timepoint-weeks-0-1_boxplot_MANUSCRIPT.png",sep=""), res=plot_res, units="in", height=4, width=3)
    ggplot(data.frame(colData(dat.se)[which(colData(dat.se)$Group=="B" & colData(dat.se)$sample_weeks<=1),]), aes(x=as.character(sample_weeks), y=CML_space, group=sample_weeks, fill=Group)) + 
      geom_boxplot() + theme_bw(base_size=16) +scale_fill_manual(values=group.cols) + theme(legend.position="none") +
      scale_y_continuous(limits=c(-300,200))
    graphics.off() 
    
    #p-value
    t.test(colData(dat.se)$CML_space[which(colData(dat.se)$Group=="B"&colData(dat.se)$sample_weeks==0)], colData(dat.se)$CML_space[which(colData(dat.se)$Group=="B"&colData(dat.se)$sample_weeks==1)], paired=T )

    }
  
  
  ###
  #  Figure 2A: 
  ###
    # mechanistic model diagram constructed in powerpoint
  
  
  ###
  #  Figure 2B: 
  ###
    # see Matlab code: XXX
  
  
  ###
  # Figure 2C: potentials and sample density 
  ###
  {
    
    ### "potential wells" using sample density ###
    {

      ### CTL ###
      # Use the spline function to interpolate the curve
      c.interp_dens <- spline(c.dens$x, c.dens$y, n = length(c.dens$x) * 10)
      # Define a function that interpolates the curve
      c.dens_function <- function(x) {
        approx(c.interp_dens$x, c.interp_dens$y, xout = x)$y
      }
      png(paste(manFig_out,"/Fig-2C_potential_C-ONLY_reverse-manusript.png",sep=""), res=plot_res, units="in", height=4, width=4) 
      plot(-1*c.dens$x, -1*c.dens$y, main="Ctl samples", pch=19, col="black", type="l", ylab="", yaxt="n", xlab="",xaxt="n", bty="n")
      points(-1*colData(dat.se)$CML_space[c.sel], -1*c.dens_function(colData(dat.se)$CML_space[c.sel])+0.00005, col="black", pch=19, cex=1.5)
      graphics.off()
      
      ### CML ###
      # Use the spline function to interpolate the curve
      interp_dens <- spline(dens$x, dens$y, n = length(dens$x) * 10)
      # Define a function that interpolates the curve
      dens_function <- function(x) {
        approx(interp_dens$x, interp_dens$y, xout = x)$y
      }
      #use critical points
      png(paste(manFig_out,"/Fig-2C_potential_B-ONLY_reverse-manusript.png",sep=""), res=plot_res, units="in", height=3, width=6) 
      plot(-1*dens$x, -1*dens$y, main="CML samples", pch=19, col="red", type="l", ylab="", yaxt="n", xlab="",xaxt="n", bty="n")
      points(-1*colData(dat.se)$CML_space[b.sel], -1*dens_function(colData(dat.se)$CML_space[b.sel])+0.00005, col="red", pch=19, cex=1.5)
      graphics.off()
    }
    
    ### Histograms ###
    png(paste(manFig_out,"/Fig-2C_CML_space_B-ONLY_histogram-manuscript.png",sep=""), res=plot_res, units="in", height=3, width=5)
    ggplot(data = data.frame(colData(dat.se)[which(colData(dat.se)$Group=="B"),]), aes(x = CML_space)) +
      geom_histogram(aes(x=CML_space, y=..density.. ), bins=10,  fill="grey", color="grey") + geom_density(aes(color=Group), adjust=.75, linewidth=2) +
      geom_vline(xintercept=ss.cps, color="black", linetype="dashed") + scale_x_reverse() + 
      scale_color_manual(values=unlist(group.cols)) +theme_bw(base_size=18) + theme(legend.position="none")
    graphics.off()
    png(paste(manFig_out,"/Fig-2C_CML_space_C-ONLY_histogram-manuscript.png",sep=""), res=plot_res, units="in", height=3, width=5)
    ggplot(data = data.frame(colData(dat.se)[which(colData(dat.se)$Group=="C"),]), aes(x = CML_space)) +
      geom_histogram(aes(x=CML_space, y=..density.. ), bins=10,  fill="grey", color="grey") + geom_density(aes(color=Group), adjust=2, linewidth=2)  +
      geom_vline(xintercept=ss.cps, color="black", linetype="dashed") + scale_x_reverse() + 
      scale_color_manual(values=unlist(group.cols)) +theme_bw(base_size=18) + theme(legend.position="none")
    graphics.off()

  } #end Fig 2C
  
  ### SUPPLEMENTAL ###
  
  ###
  # Figure S1: myeloid vs CML
  ###
  {
    # myeloid vs CML - Fig 1C
    png(paste(manFig_out,"/Fig-S1_myeloid.vs.CML_B-ONLY_space_criticalPoints_line+point.png",sep=""), res=plot_res, units="in", height=4, width=3)
    ggplot( col.df[b.sel,], aes(x=myeloid, y=CML_space, color=Group)) + geom_point() +
      geom_path(data=col.df[b.sel,], aes(x=myeloid, y=CML_space, group=mouse_id)) + 
      geom_smooth( method="lm", se=F, col="grey40", linewidth=1.5) + scale_color_manual(values=unlist(group.cols)) +
      geom_smooth( method="lm", se=F, col="black", linetype="dashed", linewidth=1.5) + scale_y_continuous(limits = c(-320, 200)) +
      theme_bw() + theme(legend.position="none") 
    graphics.off()
  }
  
  ###
  # Figure S1: log-BCR-ABL and CML vs all weeks: boxplot and trajectories
  ###
  {
    
    png(paste(manFig_out,"/Fig-S1_log-BCR-ABL.vs.CML_B-ONLY_space_criticalPoints_line+point-MANUSCRIPT.png",sep=""), res=plot_res, units="in", height=4, width=3)
    ggplot( col.df[b.sel,], aes(x=log(HSA_BCR_ABL1_gene), y=CML_space, color=Group)) + geom_point() +
      geom_line(data=col.df[b.sel,], aes(x=log(HSA_BCR_ABL1_gene), y=CML_space, group=mouse_id)) + 
      geom_smooth( method="lm", se=F, col="grey40", linewidth=1.5) + scale_color_manual(values=unlist(group.cols)) +
      geom_smooth( method="lm", se=F, col="black", linetype="dashed", linewidth=1.5) + scale_y_continuous(limits = c(-320, 200)) +
      theme_bw() + theme(legend.position="none") 
    graphics.off()
    
    png(paste(manFig_out,"/Fig-S1_HSA-gene-",hg,"_group-B-ONLY_timepoint-weeks_boxplot_MANUSCRIPT.png",sep=""), res=plot_res, units="in", height=4, width=5)
    ggplot(data.frame(colData(dat.se)[which(colData(dat.se)$Group=="B" ),]), aes(x=sample_weeks, y=.data[[hg]], group=sample_weeks)) + 
      geom_boxplot(fill="dimgrey") + theme_bw(base_size=16) + theme(legend.position="none") + scale_y_continuous(limits=c(0,50))
    graphics.off() 
    
    png(paste(manFig_out,"/Fig-S1_HSA-gene-",hg,"_group-B-ONLY_timepoint-weeks_boxplot-log.png",sep=""), res=plot_res, units="in", height=6, width=10)
    ggplot(data.frame(colData(dat.se)[which(colData(dat.se)$Group=="B"),]), aes(x=sample_weeks, y=log(.data[[hg]]), group=sample_weeks)) + 
      geom_boxplot(fill=group.cols[["B"]]) + theme_bw(base_size=16) 
    graphics.off() 

  }
  
  ###
  # Fig S2
  ### 
  {
    png(paste(manFig_out,"/Fig-S2A_MSD_all-Grps_samplePath.png",sep=""), res=plot_res, units="in", height=4, width=6)
    ggplot(msd.df, aes(x = time_interval, y = MSD, color = group, group=mouse_id)) + geom_point(alpha=.5)+
      geom_path() +
      labs(x = "Time", y = "MSD") +
      theme_minimal()  + scale_color_manual(values=group.cols)
    graphics.off()

    ### all Rx pre-Rx only ###
    prerx.msd.df <- msd.df
    prerx.msd.df[which(!(prerx.msd.df$group=="D" & prerx.msd.df$time_interval>5)),]
    png(paste(manFig_out,"/Fig-S2B_MSD_all-Grps-pre-Rx-Only_linearFit.png",sep=""), res=plot_res, units="in", height=4, width=6)
    ggplot(prerx.msd.df[which(!(prerx.msd.df$group=="D" & prerx.msd.df$time_interval>5)),], aes(x = time_interval, y = MSD, color = group)) + geom_point(alpha=.5)+
      geom_smooth(method="lm", linetype="solid", linewidth=2, se=F) +
      labs(x = "Time", y = "MSD") + scale_color_manual(values=group.cols)+
      theme_minimal()
    graphics.off()
    
  }
  
}



###
### Growth rate analysis
###
{
  d_myeloid <- c()
  pom_myeloid <- c()
  dcml_myeloid <- c()
  s.d_myeloid <- c()
  s.d_cml <- c()
  samps <- c()
  ### build derivative ###
  for (s in unique(colData(dat.se)$mouse_id)) {
    ssel <- which(colData(dat.se)$mouse_id==s)
    ctime <- colData(dat.se)$sample_weeks[ssel]
    cmye <- c()
    csamp <- c()
    ccml <- c()
    for (t in ctime) { #build myeloid percentage
      cmye <- c(cmye, colData(dat.se)$myeloid[intersect(ssel, which(colData(dat.se)$sample_weeks==t) )] )
      ccml <- c(ccml, colData(dat.se)$CML_space[intersect(ssel, which(colData(dat.se)$sample_weeks==t) )] )
      csamp <- c(csamp, colData(dat.se)$sample[intersect(ssel, which(colData(dat.se)$sample_weeks==t) )] )
    }
    
    #spline fit
    s.time <- smooth.spline(ctime, cmye)
    s.cml <- smooth.spline(ccml, cmye)
    
    # Calculate the differences between consecutive x and y values
    dx <- diff(as.numeric(ctime))
    dy <- diff(as.numeric(cmye))
    dcml <- diff(as.numeric(ccml))
    
    # Calculate the derivative as the ratio of the differences
    derivative <- dy / dx
    d_cml <- dy / dcml
    # Add NA at the beginning to align with the original x values
    derivative <- c(0, derivative)
    d_cml <- c(0, d_cml)
    s.dm <- predict(s.time, ctime, deriv = 1)$y
    s.dcml <- predict(s.cml, ccml, deriv = 1)$y

    #add values to lists
    d_myeloid <- c(d_myeloid, derivative)
    pom_myeloid <- c(pom_myeloid, derivative / max(derivative))
    samps <- c(samps, csamp)
    dcml_myeloid <- c(dcml_myeloid, d_cml)
    s.d_myeloid <- c(s.d_myeloid, s.dm)
    s.d_cml <- c(s.d_cml, s.dcml)
  }
  
  tmp.df <- data.frame("d_myeloid" = d_myeloid, "sample"=samps, "pom_myeloid"=pom_myeloid,"dcml_myeloid"=dcml_myeloid,
                       "s.d_myeloid"=s.d_myeloid, "s.d_cml"=s.d_cml)
  meta.df <- merge(colData(dat.se), tmp.df, by="sample")
  meta.df <- data.frame(meta.df)
  # 5 STATE
  s5 <- c()
  for (i in 1:1-length(ss.cps)) {
    s5 <- c(s5, ss.cps[i] + (ss.cps[i+1]-ss.cps[i]) )
  }
  state.5 <- rep(NA, dim(dat.se)[2])
  state.5[which(colData(dat.se)$CML_space < s5[1])] <- "s5"
  state.5[which(colData(dat.se)$CML_space >= s5[1] & colData(dat.se)$CML_space < s5[2])] <- "s4"
  state.5[which(colData(dat.se)$CML_space >= s5[2] & colData(dat.se)$CML_space < s5[3])] <- "s3"
  state.5[which(colData(dat.se)$CML_space >= s5[3] & colData(dat.se)$CML_space < s5[4])] <- "s2"
  state.5[which(colData(dat.se)$CML_space >= s5[4] )] <- "s1"
  meta.df[["state.5"]] <- state.5
  

  ### spline ttest ###
  #c3 vs ctl
  t.test(meta.df$s.d_myeloid[which(meta.df$state_rx=="c3"&meta.df$Group=="B")], meta.df$s.d_myeloid[which(meta.df$Group=="C")])
  #c3 vs c1
  t.test(meta.df$s.d_myeloid[which(meta.df$state_rx=="c3"&meta.df$Group=="B")], meta.df$s.d_myeloid[which(meta.df$state_rx=="c1"&meta.df$Group=="B")])
  #c2 vs c1
  t.test(meta.df$s.d_myeloid[which(meta.df$state_rx=="c2"&meta.df$Group=="B")], meta.df$s.d_myeloid[which(meta.df$state_rx=="c1"&meta.df$Group=="B")])
  #c3 vs c2
  t.test(meta.df$s.d_myeloid[which(meta.df$state_rx=="c3"&meta.df$Group=="B")], meta.df$s.d_myeloid[which(meta.df$state_rx=="c2"&meta.df$Group=="B")])
  

  
  ### ADD DATA TO "dat.se"
  # !!!using splined data as THE derivative!!! i.e. s.d_myeloid is input as d_myeloid
  # :note: already included in "Robj_dat.se_final.rds"
  # identical(meta.df$sample, rownames(colData(dat.se)))
  # cur.meta <- colData(dat.se)
  # colnames(meta.df)
  # mye.meta <- data.frame(meta.df$s.d_myeloid)
  # colnames(mye.meta) <- "d_myeloid"
  # mye.meta[["sample"]] <- meta.df$sample
  # new.meta <- merge(cur.meta, mye.meta, by="sample")
  # identical(new.meta$sample, rownames(colData(dat.se)))
  # new.ord <- data.frame(new.meta) %>% arrange(factor(sample, levels=colData(dat.se)$sample))
  # identical(new.ord$sample, rownames(colData(dat.se)))
  # dim(colData(dat.se))
  # dim(new.ord) #should have one more column
  # colData(dat.se)[["d_myeloid"]] <- new.ord$d_myeloid #!!!WARNING!!! make sure this is correct before executing!!!
}



###
### DESeq2 differential expression analysis
###
{
  deg_out <- "DEG_output"
  dir.create(deg_out, showWarnings = F)
  
  ### CML vs Control DEGs ###
  {
    comps <- c("B-state-c3.vs.C", "B-state-c2.vs.C", "B-state-c1.vs.C", "B-state-c3.vs.B-state-c1", "B-state-c3.vs.B-state-c2", "B-state-c2.vs.B-state-c1",
               "C_early.vs.late", "C_late.vs.early")
    #objects to use for looping
    trt1sel <- c("TET_OFF_B", "TET_OFF_B", "TET_OFF_B", "TET_OFF_B", "TET_OFF_B", "TET_OFF_B",
                 "TET_ON_C", "TET_ON_C"  )
    trt2sel <- c("TET_ON_C", "TET_ON_C", "TET_ON_C", "TET_OFF_B", "TET_OFF_B", "TET_OFF_B", 
                 "TET_ON_C", "TET_ON_C")
    class1name <- c("state_chr", "state_chr", "state_chr", "state_chr", "state_chr", "state_chr",
                    "sample_weeks", "sample_weeks")
    class2name <- c(NA, NA, NA, "state_chr", "state_chr", "state_chr",
                    "sample_weeks", "sample_weeks")
    class1sel <- c("c3", "c2", "c1", "c3", "c3", "c2",  "early", "late")
    class2sel <- c(NA, NA, NA, "c1", "c2", "c1",  "late", "early")
    
    
    ### loop through each comparison ###
    for (ind in 1:length(comps) ) {
    #for (ind in 1:1 ) {  
      compname <- comps[ind]
      trt1 <- trt1sel[ind]
      trt2 <- trt2sel[ind]
      class1 <- class1sel[ind]
      class2 <- class2sel[ind]
      cname1 <- class1name[ind]
      cname2 <- class2name[ind]
      if (trt1 == trt2) {
        lab1 <- paste(trt1, class1, sep="_")
        lab2 <- paste(trt2, class2, sep="_")
      } else {
        lab1 <- trt1
        lab2 <- trt2
      }
      print(paste("Processing: ",compname,sep=""))
      # class 1
      if (class1 == "early" | class1=="late") {
        if (class1=="early") {
          sel1 <- which(colData(dat.se)$treatment==trt1 & colData(dat.se)[[cname1]]<=4 )
        } else if (class1=="late") {
          sel1 <- which(colData(dat.se)$treatment==trt1 & colData(dat.se)[[cname1]]>=14 )
        }
      } else {
        sel1 <- which(colData(dat.se)$treatment==trt1 & colData(dat.se)[[cname1]]==class1 )
      }
      # class 2
      if (is.na(class2)) {
        sel2 <- which(colData(dat.se)$treatment==trt2 )
      } else if (class2 == "early" | class2=="late") {
        if (class2=="early") {
          sel2 <- which(colData(dat.se)$treatment==trt2 & colData(dat.se)[[cname2]]<=4 )
        } else if (class2=="late") {
          sel2 <- which(colData(dat.se)$treatment==trt2 & colData(dat.se)[[cname2]]>=14 )
        } 
      } else {
        sel2 <- which(colData(dat.se)$treatment==trt2 & colData(dat.se)[[cname2]]==class2 )  
      }
      comp.se <- dat.se[, c(sel1, sel2) ]
      #add label to comp.se
      colData(comp.se)[["label"]] <- c( rep(lab1, length(sel1)), rep(lab2, length(sel2)))
      col.df <- data.frame(colData(comp.se))
      png(paste(plot_out,"/DESeq-sexFactor_",compname,"_line+point.png",sep=""), res=plot_res, units="in", height=6, width=10)
      p <- ggplot(data.frame(colData(dat.se)), aes(x=sample_weeks, y=CML_space, group=mouse_id )) + geom_line(color="grey", alpha=.5) + geom_point(color="grey", alpha=.5, size=.75) + 
        geom_point(data=col.df, aes(x=sample_weeks, y=CML_space, color=label, group=mouse_id), size=2) +
        theme_bw(base_size=18) + geom_hline(yintercept=ss.cps, color="black", linetype="dashed", alpha=.75)
      print(p)
      graphics.off()
      png(paste(plot_out,"/DESeq-sexFactor_",compname,"_space-boxplot.png",sep=""), res=plot_res, units="in", height=6, width=8)
      p <- ggplot(col.df, aes(x=label, y=CML_space, fill=label )) + geom_boxplot() + theme_bw() + 
        theme(axis.title.x = element_blank(), axis.text.x=element_blank()) + coord_cartesian(ylim = c(-350, 225) ) 
      print(p)
      graphics.off()
      
      #txi object
      lens <- matrix(rep(rowData(comp.se)$basepairs, dim(comp.se)[2]), nrow=dim(comp.se)[1], ncol=dim(comp.se)[2] ) 
      colnames(lens) <- colnames(comp.se)
      rownames(lens) <- rownames(comp.se)
      txi.comp <- list( "abundance" = data.matrix(assay(comp.se, "abundance")), "counts"=data.matrix(assay(comp.se, "counts")), "length"=lens,
                        "countsFromAbundance"="no" )
      #condition
      coldata <- data.frame("treatment"=colData(comp.se)$treatment, "sex"=colData(comp.se)$sex, "timepoint"=colData(comp.se)$timepoint, 
                            "weeks"=colData(comp.se)$sample_weeks, "label" = colData(comp.se)$label)
      #colnames(coldata) <- c("treatment", "sex"," timepoint","weeks")
      coldata$condition <- factor(coldata$label, levels=c(lab2, lab1)) #set so trt2 is reference
      #DESeq results
      comp.dss <- DESeqDataSetFromTximport(txi = txi.comp,
                                           colData = coldata,
                                           design = ~ sex + condition)
      comp.dss <- DESeq(comp.dss)
      comp.res <- results(comp.dss)
      comp.dss <- estimateSizeFactors(comp.dss)
      comp.norm <- counts(comp.dss, normalized=T)
      sort(comp.res$log2FoldChange)
      #compare DEGs
      length(which(comp.res$padj<0.05))
      cmean <- rowMeans(comp.norm[, which(coldata$condition==lab2)]) 
      bmean <- rowMeans(comp.norm[, which(coldata$condition==lab1)])
      write.table(cbind(rowData(dat.se), comp.res, cmean, bmean), paste(deg_out,"/DESeq-sexFactor_",compname,"_table-means.tsv",sep=""),sep="\t", row.names=F)
      write.table(cbind(rowData(dat.se), comp.res, cmean, bmean)[which(!is.na(comp.res$padj)) ,], paste(deg_out,"/DESeq-sexFactor_",compname,"_table-means_noNAs.tsv",sep=""), row.names=F,sep="\t")
      deg.sel <- which( comp.res$padj <0.05 & abs(comp.res$log2FoldChange) >= 2 )
      write.table(cbind(rowData(dat.se), comp.res, cmean, bmean)[deg.sel,], paste(deg_out,"/DESeq-sexFactor_",compname,"_DEG-table-means.tsv",sep=""),sep="\t", row.names=F)
      compdegs <- rep("no", dim(comp.res)[1])
      compdegs[deg.sel] <- "yes"
      ddf <- data.frame("log2FC"=comp.res$log2FoldChange, "pval" = -1*log10(comp.res$padj), "DEG" = compdegs )
      
      # enchancedVolcano
      no.na <- which(!is.na(comp.res$padj))
      keyvals <- ifelse(
        comp.res$log2FoldChange[no.na] < -2, "#05c9e3",
        ifelse(comp.res$log2FoldChange[no.na] > 2, "#e30562",
               "dimgrey"))
      keyvals[is.na(keyvals)] <- 'dimgrey'
      names(keyvals)[keyvals == "#e30562"] <- 'high'
      names(keyvals)[keyvals == 'dimgrey'] <- 'mid'
      names(keyvals)[keyvals == "#05c9e3"] <- 'low'
      
      png(paste(plot_out,"/DEG-sexFactor_",compname,"_volcano.png",sep=""), res=plot_res, units="in", height=7, width=6)
      p <- EnhancedVolcano(comp.res[no.na,],
                           lab = rowData(comp.se)$gene_name[no.na],
                           x = 'log2FoldChange',
                           y = 'padj',
                           title = compname, 
                           subtitle = paste(length(which(comp.res$padj<0.05 & comp.res$log2FoldChange>=2))," Up; ",length(which(comp.res$padj<0.05 & comp.res$log2FoldChange<=-2))," Down",sep=""),
                           pCutoff = 0.05,
                           FCcutoff = 2,
                           legendPosition = 'none',
                           pointSize = 3.0,
                           labSize = 6.0,
                           shape = c(1, 1, 1, 4),
                           col=c('dimgrey', 'dimgrey', 'dimgrey', "#05c9e3"),
                           colCustom = keyvals,
                           selectLab = rowData(comp.se)$gene_name[no.na][which(names(keyvals) %in% c('high', 'low'))])
      print(p)
      graphics.off()
      
      tdat <- data.matrix(assay(comp.se, "counts"))[grep("HSA",rowData(comp.se)$gene_name),]
      colnames(tdat) <- paste(colData(comp.se)$Group,colnames(comp.se),sep="_")
      
      
      ### fGSEA ###
      cur_path_out <- paste(deg_out,"/fgsea_sexFactor_",compname,sep="")
      dir.create(cur_path_out, showWarnings = F)
      cur.lfc <- comp.res$log2FoldChange
      names(cur.lfc) <- rowData(dat.se)[,2]
      cur.lfc <- sort(cur.lfc)
      run_quick_fgsea(cur.lfc, cur_path_out)



      ### enrichR ###
      cur_path_out <- paste(deg_out,"/enrichr_sexFactor_",compname,sep="")
      dir.create(cur_path_out, showWarnings = F)
      degs <- rowData(comp.se)$gene_name[deg.sel]
      degs.up <- rowData(comp.se)$gene_name[which(comp.res$log2FoldChange > 2 & comp.res$padj < 0.05)]
      degs.down <- rowData(comp.se)$gene_name[which(comp.res$log2FoldChange < -2 & comp.res$padj < 0.05)]
      enriched.up <- enrichr(degs.up, test_dbs)
      enriched.down <- enrichr(degs.down, test_dbs)
      output_enrichr_results(enriched.up, "Up-DEGs", cur_path_out)
      output_enrichr_results(enriched.down, "Down-DEGs", cur_path_out)
      
      
    } # end DEG comparison loop
    
  } #end CML vs control DEGs
  
  
  
  ### TREATMENT (Rx) Group DEGs; include SEX AS FACTOR ###
  # :note: treatment starts AFTER week 6
  # :note: group D nilotnib ends at week 9; week 10 is first postT time point
  {
    comps <- c("A-postRx.vs.C", 
               "D-postRx.vs.C", "D-Rx.vs.C", 
               "A-postRx.vs.A-preRx-c1", "D-postRx.vs.D-preRx-c3", 
               "D-Rx.vs.D-preT-c1", "D-Rx.vs.D-preT-c2", "D-Rx.vs.D-preT-c3",
               "D-Rx.vs.A-postRx", "D-postRx.vs.A-postRx", "D-postRx.vs.D-preT-c2", "D-postRx.vs.D-preRx", "D-postRx.vs.B-c3",
               "A-postRx.vs.B-state-c3", "A-postRx.vs.B-state-c2", "A-postRx.vs.B-state-c1",
               "D-postRx.vs.B-state-c3", "D-postRx.vs.B-state-c2", "D-postRx.vs.B-state-c1",
               "D-Rx.vs.B-state-c3", "D-Rx.vs.B-state-c2", "D-Rx.vs.B-state-c1",
               "A-postRx.vs.C-late")
    #objects to use for looping
    trt1sel <- c("TET_OFF_ON_A", 
                 "TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D",
                 "TET_OFF_ON_A", "TET_OFF_NIL_ON_D",
                 "TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D",
                 "TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D", "TET_OFF_NIL_ON_D", "TET_OFF_NIL_ON_D", "TET_OFF_NIL_ON_D",
                 "TET_OFF_ON_A", "TET_OFF_ON_A", "TET_OFF_ON_A", #vs CML
                 "TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D", #D-postRx vs CML
                 "TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D", #D-Rx vs CML
                 "TET_OFF_ON_A"
    )
    trt2sel <- c("TET_ON_C", "TET_ON_C", "TET_ON_C", 
                 "TET_OFF_ON_A","TET_OFF_NIL_ON_D",
                 "TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D",
                 "TET_OFF_ON_A", "TET_OFF_ON_A", "TET_OFF_NIL_ON_D", "TET_OFF_NIL_ON_D", "TET_OFF_B", 
                 "TET_OFF_B", "TET_OFF_B", "TET_OFF_B", # A postRx vs B
                 "TET_OFF_B", "TET_OFF_B", "TET_OFF_B", # D postRx vs B
                 "TET_OFF_B", "TET_OFF_B", "TET_OFF_B", # D Rx vs B
                 "TET_ON_C"
    )
    class1name <- c("rx_class", 
                    "rx_class", "rx_class", 
                    "rx_class", "rx_class", 
                    "rx_class","rx_class", "rx_class", 
                    "rx_class", "rx_class", "rx_class", "rx_class", "rx_class",
                    "rx_class","rx_class", "rx_class", #A postRX vs B
                    "rx_class","rx_class", "rx_class", #D postRx vs B
                    "rx_class","rx_class", "rx_class", #D Rx vs B
                    "rx_class"
    )
    class2name <- c(NA,NA,NA, 
                    "state_chr", "state_chr", 
                    "state_chr", "state_chr", "state_chr", 
                    "rx_class", "rx_class", "state_chr", "rx_class", "state_chr", 
                    "state_chr", "state_chr", "state_chr", #A postRx vs B
                    "state_chr", "state_chr", "state_chr", #D postRx vs B
                    "state_chr", "state_chr", "state_chr", #D Rx vs B
                    "sample_weeks"
    )
    class1sel <- c("A.post", "D.post", "D.Rx", "A", "D.post", "D.Rx", "D.Rx", "D.Rx", "D.Rx", "D.post" , "D.post", "D.post", "D.post",
                   "A.post", "A.post", "A.post",
                   "D.post", "D.post", "D.post", 
                   "D.Rx", "D.Rx", "D.Rx", 
                   "A.post")
    class2sel <- c(NA,NA,NA, 
                   "c1","c3",
                   "c1","c2","c3", 
                   "A.post","A.post", "c2", "D.pre", "c3",
                   "c3","c2","c1", 
                   "c3","c2","c1",
                   "c3","c2","c1",
                   "late")
    #pre or post treatment samples; NA for non-Rx groups (B,C)
    #:this isnt needed with new "rx_class"
    #rx1sel <- c("post","post","rx", "post","post","rx", "rx","rx", "rx", "post", "post", "post", "post")
    #rx2sel <- c(NA,NA,NA, "pre","pre", "pre","pre","pre", "post", "post", "pre", "pre", NA)
    
    
    for (ind in 1:length(comps) ) {
      compname <- comps[ind]
      trt1 <- trt1sel[ind]
      trt2 <- trt2sel[ind]
      class1 <- class1sel[ind]
      class2 <- class2sel[ind]
      cname1 <- class1name[ind]
      cname2 <- class2name[ind]
      if (!is.na(class2)) {
        lab1 <- paste(trt1, class1,  sep="_")
        lab2 <- paste(trt2, class2,  sep="_")
      } else if (trt1 == trt2) {
        lab1 <- paste(trt1, class1, sep="_")
        lab2 <- paste(trt2, class2, sep="_")
      } else {
        lab1 <- trt1
        lab2 <- trt2
      }
      
      print(paste("Processing: ",compname,sep=""))
      
      ### condition 1; no exception for class == NA
      if (class1 == "early" | class1=="late") {
        if (class1=="early") {
          sel1 <- which(colData(dat.se)$treatment==trt1 & colData(dat.se)[[cname1]]<=4 )
        } else if (class1=="late") {
          sel1 <- which(colData(dat.se)$treatment==trt1 & colData(dat.se)[[cname1]]>=14 )
        }
      } else {
        sel1 <- which(colData(dat.se)$treatment==trt1 & colData(dat.se)[[cname1]]==class1 )
      }
      ### condition 2
      if (is.na(class2)  ) { #niether class2  have values
        sel2 <- which(colData(dat.se)$treatment==trt2 )
      }  else if (class2 == "early" | class2=="late") {
        if (class2=="early") {
          sel2 <- which(colData(dat.se)$treatment==trt2 & colData(dat.se)[[cname2]]<=4 )
        } else if (class2=="late") {
          sel2 <- which(colData(dat.se)$treatment==trt2 & colData(dat.se)[[cname2]]>=14 )
        } 
      } else {
        sel2 <- which(colData(dat.se)$treatment==trt2 & colData(dat.se)[[cname2]]==class2)
      }
      
      ###check conditions have matched some samples
      if (length(sel1)==0 | length(sel2)==0 ) {
        print(paste(":::ERROR::: unable to match samples to both conditions for comparison! SKIPPING!",sep=""))
        next
      }
      
      comp.se <- dat.se[, c(sel1, sel2) ]
      #add label to comp.se
      colData(comp.se)[["label"]] <- c( rep(lab1, length(sel1)), rep(lab2, length(sel2)))
      comp.df <- data.frame(colData(comp.se))
      head(comp.df$time)
      png(paste(plot_out,"/DESeq-sexFactor_",compname,"_line+point.png",sep=""), res=plot_res, units="in", height=6, width=10)
      p <- ggplot(data.frame(colData(dat.se)), aes(x=sample_weeks, y=CML_space, group=mouse_id )) + geom_line(color="grey", alpha=.5) + geom_point(color="grey", alpha=.5, size=.75) + 
        geom_point(data=comp.df, aes(x=sample_weeks, y=CML_space, color=label, group=mouse_id), size=2) +
        theme_bw(base_size=18) + geom_hline(yintercept=ss.cps, color="black", linetype="dashed", alpha=.75)
      print(p)
      graphics.off()
      png(paste(plot_out,"/DESeq-sexFactor_",compname,"_space-boxplot.png",sep=""), res=plot_res, units="in", height=6, width=8)
      p <- ggplot(comp.df, aes(x=label, y=CML_space, fill=label )) + geom_boxplot() + theme_bw() + 
        theme(axis.title.x = element_blank(), axis.text.x=element_blank()) + coord_cartesian(ylim = c(-350, 225) )
      print(p)
      graphics.off()
      
      
      #txi object
      lens <- matrix(rep(rowData(comp.se)$basepairs, dim(comp.se)[2]), nrow=dim(comp.se)[1], ncol=dim(comp.se)[2] ) 
      colnames(lens) <- colnames(comp.se)
      rownames(lens) <- rownames(comp.se)
      txi.comp <- list( "abundance" = data.matrix(assay(comp.se, "abundance")), "counts"=data.matrix(assay(comp.se, "counts")), "length"=lens,
                        "countsFromAbundance"="no" )
      #condition
      coldata <- data.frame("treatment"=colData(comp.se)$treatment, "sex"=colData(comp.se)$sex, "timepoint"=colData(comp.se)$timepoint, 
                            "weeks"=colData(comp.se)$sample_weeks, "label" = colData(comp.se)$label)
      #colnames(coldata) <- c("treatment", "sex"," timepoint","weeks")
      coldata$condition <- factor(coldata$label, levels=c(lab2, lab1)) #set so trt2 is reference
      write.table(cbind(comp.df, as.character(coldata$condition)), paste(deg_out,"/DESeq-sexFactor_",compname,"_sampleInfoTable.tsv",sep=""),sep="\t" )
      #DESeq results
      comp.dss <- DESeqDataSetFromTximport(txi = txi.comp,
                                           colData = coldata,
                                           design = ~ sex + condition)
      comp.dss <- DESeq(comp.dss)
      comp.res <- results(comp.dss)
      comp.dss <- estimateSizeFactors(comp.dss)
      comp.norm <- counts(comp.dss, normalized=T)
      sort(comp.res$log2FoldChange)
      #compare DEGs
      length(which(comp.res$padj<0.05))
      cmean <- rowMeans(comp.norm[, which(coldata$condition==lab2)]) 
      bmean <- rowMeans(comp.norm[, which(coldata$condition==lab1)])
      write.table(cbind(rowData(dat.se), comp.res, cmean, bmean), paste(deg_out,"/DESeq-sexFactor_",compname,"_table-means.tsv",sep=""),sep="\t", row.names=F)
      write.table(cbind(rowData(dat.se), comp.res, cmean, bmean)[which(!is.na(comp.res$padj)) ,], paste(deg_out,"/DESeq-sexFactor_",compname,"_table-means_noNAs.tsv",sep=""),sep="\t", row.names=F)
      deg.sel <- which( comp.res$padj <0.05 & abs(comp.res$log2FoldChange) >= 2 )
      write.table(cbind(rowData(dat.se), comp.res, cmean, bmean)[deg.sel,], paste(deg_out,"/DESeq-sexFactor_",compname,"_DEG-table-means.tsv",sep=""),sep="\t", row.names=F)
      compdegs <- rep("no", dim(comp.res)[1])
      compdegs[deg.sel] <- "yes"
      ddf <- data.frame("log2FC"=comp.res$log2FoldChange, "pval" = -1*log10(comp.res$padj), "DEG" = compdegs )
      # png(paste(plot_out,"/DEG-sexFactor_",compname,"_volcano.png",sep=""), res=plot_res, units="in", height=6, width=6)
      # p <- ggplot(ddf, aes(x=log2FC, y=pval, color=DEG, shape=DEG)) + geom_point(alpha=.75) + theme_bw(base_size=16) +
      #   scale_color_manual(values=c("dimgrey", "#05c9e3" )) + scale_shape_manual(values=c(1, 19))
      # print(p)
      # graphics.off()
      # NEW enchancedVolcano
      no.na <- which(!is.na(comp.res$padj))
      keyvals <- ifelse(
        comp.res$log2FoldChange[no.na] < -2, "#05c9e3",
        ifelse(comp.res$log2FoldChange[no.na] > 2, "#e30562",
               "dimgrey"))
      keyvals[is.na(keyvals)] <- 'dimgrey'
      names(keyvals)[keyvals == "#e30562"] <- 'high'
      names(keyvals)[keyvals == 'dimgrey'] <- 'mid'
      names(keyvals)[keyvals == "#05c9e3"] <- 'low'
      
      png(paste(plot_out,"/DEG-sexFactor_",compname,"_volcano.png",sep=""), res=plot_res, units="in", height=7, width=6)
      p <- EnhancedVolcano(comp.res[no.na,],
                           lab = rowData(comp.se)$gene_name[no.na],
                           x = 'log2FoldChange',
                           y = 'padj',
                           title = compname,
                           subtitle = paste(length(which(comp.res$padj<0.05 & comp.res$log2FoldChange>=2))," Up; ",length(which(comp.res$padj<0.05 & comp.res$log2FoldChange<=-2))," Down",sep=""),
                           pCutoff = 0.05,
                           FCcutoff = 2,
                           legendPosition = 'none',
                           pointSize = 3.0,
                           labSize = 6.0,
                           shape = c(1, 1, 1, 4),
                           col=c('dimgrey', 'dimgrey', 'dimgrey', "#05c9e3"),
                           colCustom = keyvals,
                           selectLab = rowData(comp.se)$gene_name[no.na][which(names(keyvals) %in% c('high', 'low'))])
      print(p)
      graphics.off()
      
      
      tdat <- data.matrix(assay(comp.se, "counts"))[grep("HSA",rowData(comp.se)$gene_name),]
      colnames(tdat) <- paste(colData(comp.se)$Group,colnames(comp.se),sep="_")
      
      
      ### fGSEA ###
      #unique(msigdbr(species = "mouse", category = "C5")$gs_subcat) #check included pathways
      #txi object
      cur_path_out <- paste(deg_out,"/fgsea_sexFactor_",compname,sep="")
      dir.create(cur_path_out, showWarnings = F)
      cur.lfc <- comp.res$log2FoldChange
      names(cur.lfc) <- rowData(dat.se)[,2]
      cur.lfc <- sort(cur.lfc)
      # run fgsea
      #run_fgsea(cur.lfc, cur_path_out)
      run_quick_fgsea(cur.lfc, cur_path_out)
      
      
      
      
      ### enrichR ###
      cur_path_out <- paste(deg_out,"/enrichr_sexFactor_",compname,sep="")
      dir.create(cur_path_out, showWarnings = F)
      degs <- rowData(comp.se)$gene_name[deg.sel]
      degs.up <- rowData(comp.se)$gene_name[which(comp.res$log2FoldChange > 2 & comp.res$padj < 0.05)]
      degs.down <- rowData(comp.se)$gene_name[which(comp.res$log2FoldChange < -2 & comp.res$padj < 0.05)]
      enriched.up <- enrichr(degs.up, test_dbs)
      enriched.down <- enrichr(degs.down, test_dbs)
      output_enrichr_results(enriched.up, "Up-DEGs", cur_path_out)
      output_enrichr_results(enriched.down, "Down-DEGs", cur_path_out)
      
      
    } #end treatment group DEGs
  } #end treatment DEGs
  
  
  
  ###
  ### setup DEG comparisons and Venns and UPSETs
  ###
  {
    
    ### inter-CML ###    
    b1.vs.b0 <- read.table(paste(deg_out,"/DESeq-sexFactor_B-state-c2.vs.B-state-c1_table-means_noNAs.tsv",sep=""), sep="\t", header=T)
    b2.vs.b0 <- read.table(paste(deg_out,"/DESeq-sexFactor_B-state-c3.vs.B-state-c1_table-means_noNAs.tsv",sep=""), sep="\t", header=T)
    b2.vs.b1 <- read.table(paste(deg_out,"/DESeq-sexFactor_B-state-c3.vs.B-state-c2_table-means_noNAs.tsv",sep=""), sep="\t", header=T)
    b1v0.degs <- which(b1.vs.b0$padj<0.05 & abs(b1.vs.b0$log2FoldChange)>2)
    b2v0.degs <- which(b2.vs.b0$padj<0.05 & abs(b2.vs.b0$log2FoldChange)>2)
    b2v1.degs <- which(b2.vs.b1$padj<0.05 & abs(b2.vs.b1$log2FoldChange)>2)
    #venn diagram
    set1 <- list(set1=b1.vs.b0[b1v0.degs,]$gene_name)
    set2 <- list(set2=b2.vs.b1[b2v1.degs,]$gene_name)
    set3 <- list(set3=b2.vs.b0[b2v0.degs,]$gene_name)
    
    # vs CTRL + other
    b0.vs.c <- read.table(paste(deg_out,"/DESeq-sexFactor_B-state-c1.vs.C_table-means_noNAs.tsv",sep=""), sep="\t", header=T)
    b1.vs.c <- read.table(paste(deg_out,"/DESeq-sexFactor_B-state-c2.vs.C_table-means_noNAs.tsv",sep=""), sep="\t", header=T)
    b2.vs.c <- read.table(paste(deg_out,"/DESeq-sexFactor_B-state-c3.vs.C_table-means_noNAs.tsv",sep=""), sep="\t", header=T)
    apost.vs.c <- read.table(paste(deg_out,"/DESeq-sexFactor_A-postRx.vs.C_table-means_noNAs.tsv",sep=""), sep="\t", header=T)
    dpost.vs.c <- read.table(paste(deg_out,"/DESeq-sexFactor_D-postRx.vs.C_table-means_noNAs.tsv",sep=""), sep="\t", header=T)
    drx.vs.c <- read.table(paste(deg_out,"/DESeq-sexFactor_D-Rx.vs.C_table-means_noNAs.tsv",sep=""), sep="\t", header=T)
    dpost.vs.b3 <- read.table(paste(deg_out,"/DESeq-sexFactor_D-postRx.vs.B-state-c3_table-means_noNAs.tsv",sep=""), sep="\t", header=T)
    c.early.vs.late <- read.table(paste(deg_out,"/DESeq-sexFactor_C_early.vs.late_table-means_noNAs.tsv",sep=""), sep="\t", header = T)
    c.late.vs.early <- read.table(paste(deg_out,"/DESeq-sexFactor_C_late.vs.early_table-means_noNAs.tsv",sep=""), sep="\t", header = T)
    apost.vs.clate <- read.table(paste(deg_out,"/DESeq-sexFactor_A-postRx.vs.C-late_table-means_noNAs.tsv",sep=""), sep="\t",header=T)
    drx.vs.apost <- read.table(paste(deg_out,"/DESeq-sexFactor_D-Rx.vs.A-postRx_table-means_noNAs.tsv",sep=""), sep="\t", header=T)
    
    
    b0vc.degs <- which(b0.vs.c$padj<0.05 & abs(b0.vs.c$log2FoldChange)>2)
    b1vc.degs <- which(b1.vs.c$padj<0.05 & abs(b1.vs.c$log2FoldChange)>2)
    b2vc.degs <- which(b2.vs.c$padj<0.05 & abs(b2.vs.c$log2FoldChange)>2)
    apost.degs <- which(apost.vs.c$padj<0.05 & abs(apost.vs.c$log2FoldChange)>2)
    dpost.degs <- which(dpost.vs.c$padj<0.05 & abs(dpost.vs.c$log2FoldChange)>2)
    drx.degs <- which(drx.vs.c$padj<0.05 & abs(drx.vs.c$log2FoldChange)>2)
    dpost.b3.degs <- which(dpost.vs.b3$padj<0.05 & abs(dpost.vs.b3$log2FoldChange)>2)
    c.evl.degs <- which(c.early.vs.late$padj<0.05 & abs(c.early.vs.late$log2FoldChange)>2)
    c.lve.degs <- which(c.late.vs.early$padj<0.05 & abs(c.late.vs.early$log2FoldChange)>2)
    apostvcL.degs <- which(apost.vs.clate$padj<0.05 & abs(apost.vs.clate$log2FoldChange)>2)
    rx.da.degs <- which(drx.vs.apost$padj<0.05 & abs(drx.vs.apost$log2FoldChange)>2)
    
    #DEG ids
    b0vc.genes.id <- b0.vs.c[b0vc.degs,]$gene_id
    b1vc.genes.id <- b1.vs.c[b1vc.degs,]$gene_id
    b2vc.genes.id <- b2.vs.c[b2vc.degs,]$gene_id
    rx.da.genes.id <- drx.vs.apost$gene_id[rx.da.degs]
    
    
    
    ### gene lists ###
    degs.cml_ctl <- list("c5.vs.ctl"=b2.vs.c[b2vc.degs, ]$gene_id, "c3.vs.ctl"=b1.vs.c[b1vc.degs, ]$gene_id, 
                         "c1.vs.ctl"=b0.vs.c[b0vc.degs, ]$gene_id)
    degs.cml6 <- list("c5.vs.ctl"=b2.vs.c[b2vc.degs, ]$gene_id, "c3.vs.ctl"=b1.vs.c[b1vc.degs, ]$gene_id, 
                      "c1.vs.ctl"=b0.vs.c[b0vc.degs, ]$gene_id, "c3.vs.c1"=b1.vs.b0[b1v0.degs, ]$gene_id,
                      "c5.vs.c3"=b2.vs.b1[b2v1.degs, ]$gene_id,"c5.vs.c0"=b2.vs.b0[b2v0.degs, ]$gene_id)
    ### UPREG ###
    b0vc.degs.up <- which(b0.vs.c$padj<0.05 & b0.vs.c$log2FoldChange>2)
    b1vc.degs.up <- which(b1.vs.c$padj<0.05 & b1.vs.c$log2FoldChange>2)
    b2vc.degs.up <- which(b2.vs.c$padj<0.05 & b2.vs.c$log2FoldChange>2)
    b1v0.degs.up <- which(b1.vs.b0$padj<0.05 & b1.vs.b0$log2FoldChange>2)
    b2v0.degs.up <- which(b2.vs.b0$padj<0.05 & b2.vs.b0$log2FoldChange>2)
    b2v1.degs.up <- which(b2.vs.b1$padj<0.05 & b2.vs.b1$log2FoldChange>2)
    degs.cml6.up <- list("c5.vs.ctl"=b2.vs.c[b2vc.degs.up, ]$gene_id, "c3.vs.ctl"=b1.vs.c[b1vc.degs.up, ]$gene_id, 
                         "c1.vs.ctl"=b0.vs.c[b0vc.degs.up, ]$gene_id, "c3.vs.c1"=b1.vs.b0[b1v0.degs.up, ]$gene_id,
                         "c5.vs.c3"=b2.vs.b1[b2v1.degs.up, ]$gene_id,"c5.vs.c0"=b2.vs.b0[b2v0.degs.up, ]$gene_id)
    ### DN REG ###
    b0vc.degs.dn <- which(b0.vs.c$padj<0.05 & b0.vs.c$log2FoldChange < -2)
    b1vc.degs.dn <- which(b1.vs.c$padj<0.05 & b1.vs.c$log2FoldChange < -2)
    b2vc.degs.dn <- which(b2.vs.c$padj<0.05 & b2.vs.c$log2FoldChange < -2)
    b1v0.degs.dn <- which(b1.vs.b0$padj<0.05 & b1.vs.b0$log2FoldChange < -2)
    b2v0.degs.dn <- which(b2.vs.b0$padj<0.05 & b2.vs.b0$log2FoldChange < -2)
    b2v1.degs.dn <- which(b2.vs.b1$padj<0.05 & b2.vs.b1$log2FoldChange < -2)
    degs.cml6.dn <- list("c5.vs.ctl"=b2.vs.c[b2vc.degs.dn, ]$gene_id, "c3.vs.ctl"=b1.vs.c[b1vc.degs.dn, ]$gene_id, 
                         "c1.vs.ctl"=b0.vs.c[b0vc.degs.dn, ]$gene_id, "c3.vs.c1"=b1.vs.b0[b1v0.degs.dn, ]$gene_id,
                         "c5.vs.c3"=b2.vs.b1[b2v1.degs.dn, ]$gene_id,"c5.vs.c0"=b2.vs.b0[b2v0.degs.dn, ]$gene_id)
    degs.Rx <- list("c5.vs.ctl"=b2.vs.c[b2vc.degs, ]$gene_id, "c3.vs.ctl"=b1.vs.c[b1vc.degs, ]$gene_id, 
                    "c1.vs.ctl"=b0.vs.c[b0vc.degs, ]$gene_id, "TKI_Rx"=drx.vs.c[drx.degs,]$gene_id, 
                    "TKI_postRx"=dpost.vs.c[dpost.degs,]$gene_id, "TOTO_postRx"=apost.vs.c[apost.degs,]$gene_id)
    
    
    ### UPSET OBJECTS ###
    cm.cml_ctls = ComplexHeatmap::make_comb_mat(degs.cml_ctl)
    cm.cml6 = ComplexHeatmap::make_comb_mat(degs.cml6)
    # Up- vs Down-regulated upset objects #
    cm.cml6.up = ComplexHeatmap::make_comb_mat(degs.cml6.up)
    cm.cml6.dn = ComplexHeatmap::make_comb_mat(degs.cml6.dn)
    
    
    
    ## Rx- & CML-vs-ctl
    cm.Rx = ComplexHeatmap::make_comb_mat(degs.Rx)
    ### treatment upset objects ###
    # UP #
    apost.degs.up <- which(apost.vs.c$padj<0.05 & apost.vs.c$log2FoldChange>2)
    dpost.degs.up <- which(dpost.vs.c$padj<0.05 & dpost.vs.c$log2FoldChange>2)
    drx.degs.up <- which(drx.vs.c$padj<0.05 & drx.vs.c$log2FoldChange>2)
    degs.Rx.up <- list("c5.vs.ctl"=b2.vs.c[b2vc.degs.up, ]$gene_id, "c3.vs.ctl"=b1.vs.c[b1vc.degs.up, ]$gene_id, 
                       "c1.vs.ctl"=b0.vs.c[b0vc.degs.up, ]$gene_id, "TKI_Rx"=drx.vs.c[drx.degs.up,]$gene_id, 
                       "TKI_postRx"=dpost.vs.c[dpost.degs.up,]$gene_id, "TOTO_postRx"=apost.vs.c[apost.degs.up,]$gene_id)
    cm.Rx.up = ComplexHeatmap::make_comb_mat(degs.Rx.up)
    
    # DOWN #
    apost.degs.dn <- which(apost.vs.c$padj<0.05 & apost.vs.c$log2FoldChange < -2)
    dpost.degs.dn <- which(dpost.vs.c$padj<0.05 & dpost.vs.c$log2FoldChange < -2)
    drx.degs.dn <- which(drx.vs.c$padj<0.05 & drx.vs.c$log2FoldChange < -2)
    degs.Rx.dn <- list("c5.vs.ctl"=b2.vs.c[b2vc.degs.dn, ]$gene_id, "c3.vs.ctl"=b1.vs.c[b1vc.degs.dn, ]$gene_id, 
                       "c1.vs.ctl"=b0.vs.c[b0vc.degs.dn, ]$gene_id, "TKI_Rx"=drx.vs.c[drx.degs.dn,]$gene_id, 
                       "TKI_postRx"=dpost.vs.c[dpost.degs.dn,]$gene_id, "TOTO_postRx"=apost.vs.c[apost.degs.dn,]$gene_id)
    cm.Rx.dn = ComplexHeatmap::make_comb_mat(degs.Rx.dn)
    
    
    ### RX comparison vs CML-Ctl ###
    rx.comp.cml <- list("TKI.vs.TOTO" = drx.vs.apost$gene_name[rx.da.degs], "c5.vs.ctl"=b2.vs.c[b2vc.degs, ]$gene_name, "c3.vs.ctl"=b1.vs.c[b1vc.degs, ]$gene_name, 
                        "c1.vs.ctl"=b0.vs.c[b0vc.degs, ]$gene_name)
    cm.Rx.cml = ComplexHeatmap::make_comb_mat(rx.comp.cml)
  } 
  
  ###
  ### DEG lists
  ###
  {
    #deg gene names
    dg.b3.c <- b2.vs.c[b2vc.degs, ]$gene_name
    dg.b2.c <- b1.vs.c[b1vc.degs, ]$gene_name 
    dg.b1.c <- b0.vs.c[b0vc.degs, ]$gene_name
    dg.b2.b1 <- b1.vs.b0[b1v0.degs,]$gene_name
    dg.b3.b2 <- b2.vs.b1[b2v1.degs,]$gene_name
    dg.b3.b1 <- b2.vs.b0[b2v0.degs,]$gene_name
    dg.apost.c <- apost.vs.c$gene_name[apost.degs]
    dg.dpost.c <- dpost.vs.c$gene_name[dpost.degs]
    dg.drx.c <- drx.vs.c$gene_name[drx.degs]
    dg.dpost.b3 <- dpost.vs.b3$gene_name[dpost.b3.degs]
    dg.c.evl <- c.early.vs.late$gene_name[c.evl.degs]
    dg.c.lve <- c.late.vs.early$gene_name[c.lve.degs]
    dg.apost.cL <- apost.vs.clate$gene_name[apostvcL.degs]
    dg.rx.da <- drx.vs.apost$gene_name[rx.da.degs]
    
    #gene ids
    ig.b3.c <- b2.vs.c[b2vc.degs, ]$gene_id
    ig.b2.c <- b1.vs.c[b1vc.degs, ]$gene_id 
    ig.b1.c <- b0.vs.c[b0vc.degs, ]$gene_id
    ig.b2.b1 <- b1.vs.b0[b1v0.degs,]$gene_id
    ig.b3.b2 <- b2.vs.b1[b2v1.degs,]$gene_id
    ig.b3.b1 <- b2.vs.b0[b2v0.degs,]$gene_id
    ig.apost.c <- apost.vs.c$gene_id[apost.degs]
    ig.dpost.c <- dpost.vs.c$gene_id[dpost.degs]
    ig.drx.c <- drx.vs.c$gene_id[drx.degs]
    ig.dpost.b3 <- dpost.vs.b3$gene_id[dpost.b3.degs]
    ig.c.evl <- c.early.vs.late$gene_id[c.evl.degs]
    ig.c.lve <- c.late.vs.early$gene_id[c.lve.degs]
    ig.apost.cL <- apost.vs.clate$gene_id[apostvcL.degs]
    ig.rx.da <- drx.vs.apost$gene_id[rx.da.degs]
  }
  
  ###
  ### fgsea results
  ###
  {
    ### vs CTRL
    f.b1.c <- read.table(paste(deg_out,"/fgsea_sexFactor_B-state-c1.vs.C/fgsea_Hallmark_table.tsv",sep=""), header=T)
    f.b2.c <- read.table(paste(deg_out,"/fgsea_sexFactor_B-state-c2.vs.C/fgsea_Hallmark_table.tsv",sep=""), header=T)
    f.b3.c <- read.table(paste(deg_out,"/fgsea_sexFactor_B-state-c3.vs.C/fgsea_Hallmark_table.tsv",sep=""), header=T)
    f.apost.c <- read.table(paste(deg_out,"/fgsea_sexFactor_A-postRx.vs.C/fgsea_Hallmark_table.tsv",sep=""), header=T)
    f.dpost.c <- read.table(paste(deg_out,"/fgsea_sexFactor_D-postRx.vs.C/fgsea_Hallmark_table.tsv",sep=""), header=T)
    f.drx.c <- read.table(paste(deg_out,"/fgsea_sexFactor_D-Rx.vs.C/fgsea_Hallmark_table.tsv",sep=""), header=T)
    f.dpost.b3 <- read.csv(paste(deg_out,"/fgsea_sexFactor_D-postRx.vs.B-state-c3/fgsea_Hallmark_table.tsv",sep=""), header=T, sep="\t")
    f.c.evl <- read.csv(paste(deg_out,"/fgsea_sexFactor_C_early.vs.late/fgsea_Hallmark_table.tsv",sep=""), header=T, sep="\t")
    f.c.lve <- read.csv(paste(deg_out,"/fgsea_sexFactor_C_late.vs.early/fgsea_Hallmark_table.tsv",sep=""), header=T, sep="\t")
    f.apost.clate <- read.table(paste(deg_out,"/fgsea_sexFactor_A-postRx.vs.C-late/fgsea_Hallmark_table.tsv",sep=""), header=T)
    f.rx.da <- read.table(paste(deg_out,"/fgsea_sexFactor_D-Rx.vs.A-postRx/fgsea_Hallmark_table.tsv",sep=""), header=T)
    
    ### vs State
    f.b2.b0 <- read.table(paste(deg_out,"/fgsea_sexFactor_B-state-c3.vs.B-state-c1/fgsea_Hallmark_table.tsv",sep=""), header=T)
    f.b2.b1 <- read.table(paste(deg_out,"/fgsea_sexFactor_B-state-c3.vs.B-state-c2/fgsea_Hallmark_table.tsv",sep=""), header=T)
    f.b1.b0 <- read.table(paste(deg_out,"/fgsea_sexFactor_B-state-c2.vs.B-state-c1/fgsea_Hallmark_table.tsv",sep=""), header=T)
    
  }
  
  ###
  ### early, transition, late DEG definitions
  ###
  {
    #early
    early.deg.genes <- setdiff(b0.vs.c[b0vc.degs, ]$gene_name, union( b1.vs.b0[b1v0.degs, ]$gene_name, b2.vs.b1[b2v1.degs, ]$gene_name ))
    early.deg.genes.id <- setdiff(b0.vs.c[b0vc.degs, ]$gene_id, union( b1.vs.b0[b1v0.degs, ]$gene_id, b2.vs.b1[b2v1.degs, ]$gene_id ))
    #sel.early <- match(early.deg.genes, b1.vs.c$gene_name )
    sel.early <- match_deg_genes( early.deg.genes, b0.vs.c)
    
    #transition
    trans.deg.genes <- setdiff(b1.vs.c[b1vc.degs, ]$gene_name, union(b0.vs.c[b0vc.degs, ]$gene_name, b2.vs.b1[b2v1.degs, ]$gene_name ))
    trans.deg.genes.id <- setdiff(b1.vs.c[b1vc.degs, ]$gene_id, union(b0.vs.c[b0vc.degs, ]$gene_id, b2.vs.b1[b2v1.degs, ]$gene_id ))
    sel.trans <- match_deg_genes( trans.deg.genes, b1.vs.c)
    
    #late
    late.deg.genes <- setdiff(b2.vs.c[b2vc.degs, ]$gene_name, union(b1.vs.c[b1vc.degs, ]$gene_name, b0.vs.c[b0vc.degs, ]$gene_name))
    late.deg.genes.id <- setdiff(b2.vs.c[b2vc.degs, ]$gene_id, union(b1.vs.c[b1vc.degs, ]$gene_id, b0.vs.c[b0vc.degs, ]$gene_id))
    sel.late <- match_deg_genes( late.deg.genes, b2.vs.c)
  }
  
  
}# end DESeq section



###
### Loading value analysis
###
{
  #
  # :notes: 
  #   - this section use the eigengenes to calculate the CML contribution for each gene
  #   - analysis/prepartion for Figure 3
  #
  
  
  # make new tables for some reason....
  tab.b3.c <- b2.vs.c
  tab.b2.c <- b1.vs.c
  tab.b1.c <- b0.vs.c
  tab.b2.b1 <- b1.vs.b0
  tab.b3.b2 <- b2.vs.b1
  tab.b3.b1 <- b2.vs.b0
  
  #build data
  # for each sample group comparison, a list of genes, the DEG table, a name, and a GSEA result object are stored in lists
  #:note: each DEG comparison only contains data for the genes that were able to be tested; untested genes have pvalue, l2fc, and means -> NA
  ### !!!NOTE!!! CML contribution ASSUMES that all comparisons use the more control-like set as reference
  ###     Expression direction is assumed to be oriented so that increses in expression move away from controls (!ONLY TRUE IF CONTROL-LIKE SAMPLE IS THE REFERENCE!)
  rownames(rV) <- rowData(dat.se[lowExpRM,])$gene_name #add row names
  amin <- minexp #set offset value
  gene.df <- data.frame("gene"=rownames(rV), "gene_id"=rowData(dat.se[lowExpRM,])$gene_id, "eigengenes"=rV[,2], "PC1"=rV[,1])
  deglist <- list(dg.b3.c, dg.b2.c, dg.b1.c, dg.b2.b1, dg.b3.b2, dg.b3.b1, dg.apost.c, dg.dpost.c, dg.drx.c, dg.dpost.b3, dg.c.evl, dg.c.lve, dg.apost.cL, dg.rx.da  )
  degids <- list(ig.b3.c, ig.b2.c, ig.b1.c, ig.b2.b1, ig.b3.b2, ig.b3.b1, ig.apost.c, ig.dpost.c, ig.drx.c, ig.dpost.b3, ig.c.evl, ig.c.lve, ig.apost.cL, ig.rx.da )
  degtab <- list(tab.b3.c, tab.b2.c, tab.b1.c, tab.b2.b1, tab.b3.b2, tab.b3.b1, apost.vs.c, dpost.vs.c, drx.vs.c, dpost.vs.b3, c.early.vs.late, c.late.vs.early, apost.vs.clate, drx.vs.apost )
  degcomp <- list("B3.vs.C", "B2.vs.C", "B1.vs.C", "B2.vs.B1", "B3.vs.B2" ,"B3.vs.B1", "Apost.vs.C", "Dpost.vs.C", "Drx.vs.C", "Dpost.vs.B3", "C_early.vs.late", "C_late.vs.early", "Apost.vs.Clate", "Drx.vs.Apost" )
  gsea.objs <- list(f.b3.c, f.b2.c, f.b1.c, f.b1.b0, f.b2.b1, f.b2.b0, f.apost.c, f.dpost.c, f.drx.c, f.dpost.b3, f.c.evl, f.c.lve, f.apost.clate, f.rx.da   )
  combcomp <- c(rep("NA", dim(gene.df)[1]))
  for (i in 1:length(degcomp)) {
    curgenes <- deglist[[i]]
    curids <- degids[[i]] 
    curcomp <- degcomp[[i]]
    # sel <- match_genes(gene.df$gene, curgenes)
    sel <- match_genes(gene.df$gene_id, curids)
    curlab <- rep("no", dim(gene.df)[1])
    curlab[sel] <- "yes"
    gene.df[[curcomp]] <- curlab
    #pvalue
    curtab <- degtab[[i]]
    psel <- match_genes(gene.df$gene, curtab$gene_name)
    curpv <- curtab$padj
    plab <- rep(1, dim(gene.df)[1])
    plab[psel] <- curpv
    gene.df[[paste(curcomp,"_padj",sep="")]] <- as.numeric(plab)
    #log2FC
    curfc <- curtab$log2FoldChange
    fclab <- rep(NA, dim(gene.df)[1])
    fclab[psel] <- curfc
    gene.df[[paste(curcomp,"_lfc",sep="")]] <- as.numeric(fclab)
    #group means
    m1lab <- rep(NA, dim(gene.df)[1] ) 
    m1lab[psel] <- curtab$bmean
    gene.df[[paste(curcomp,"_mean.",1,sep="")]] <- as.numeric(m1lab)
    
    m2lab <- rep(NA, dim(gene.df)[1] )
    m2lab[psel] <- curtab$cmean
    gene.df[[paste(curcomp,"_mean.",2,sep="")]] <- as.numeric(m2lab)
    
    #append to combined label
    combcomp[sel] <- paste(combcomp[sel],curcomp, sep="|")
  }
  comblab <- unlist(lapply(combcomp, function(x) gsub("NA\\|", "", x)))
  gene.df[["degs"]] <- comblab

  
  
  
  ###
  ### make pathway contribution plots (loading val * exp change) ###
  # :note: only genes included in DEG testing will contribute for each comparison
  #   NAs are produced when gene were not tested (low expression, etc.)
  ###
  {
    msigdb.mm = getMsigdb(org = 'mm', id = 'SYM', version = '7.4')
    mm.h <- subsetCollection(msigdb.mm, "h")
    cont.tab <- c() #holds eigengene, contribution (diff.exp and l2fc) using keys all, up, down
    cml.tab <- c() # holds pro-CML and anti CML using pos|up, pos|down, neg|up, neg|down
    cml.full <- c()
    for (curpath in names(mm.h) ){
      cgenes <- geneIds(mm.h[[curpath]])
      sel <- match_all_genes( gene.df$gene, cgenes)
      c.df <- gene.df[sel,]
      for (i in 1:length(degcomp)) {
        curcomp <- degcomp[[i]]
        m.1 <- c.df[[paste(curcomp,"_mean.",1,sep="")]]
        m.2 <- c.df[[paste(curcomp,"_mean.",2,sep="")]]
        c.diff <- m.1 - m.2
        c.ld <- c.df$eigengenes
        up.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]
        dn.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]
        c.all <- c.diff * c.ld #get contribution for all genes in pathway
        c.up <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
        c.down <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
        ctab <- matrix( c( rep(curpath, 3), rep(curcomp, 3), "all", "up", "down", sum(c.ld), sum(up.ld), sum(dn.ld), sum(c.all,na.rm=T), sum(c.up), sum(c.down) ), nrow=3, ncol=5) 
        #add lfc contribution; need to figure out why "Inf" is produced
        c.lfc <- log2((m.1+amin)/(m.2+amin) )
        cl.all <- c.ld * c.lfc #get contribution for all genes in pathway based on LFC
        cl.up <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
        cl.down <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
        ctab <- matrix(c( rep(curpath, 3), rep(curcomp, 3), "all", "up", "down", 
                          length(which(!is.na(c.all))), length(which(!is.na(c.up))), length(which(!is.na(c.down))), # count of genes contributing (non NA genes)
                          sum(c.ld), sum(up.ld), sum(dn.ld), #eigenvalue
                          sum(c.all,na.rm=T), sum(c.up, na.rm=T), sum(c.down, na.rm=T), #contribution w/ exp diff
                          sum(cl.all,na.rm=T), sum(cl.up, na.rm=T), sum(cl.down, na.rm=T) ), # contrib. w/ log2FC
                       nrow=3, ncol=7)
        if (length(which(is.na(ctab))) != 0) {
          print(paste(":::ERROR::: NAs produced for comp: ",curcomp," pathway: ",curpath,sep=""))
          print(ctab)
          stop("NA encountered")
        }
        cont.tab <- rbind(cont.tab, ctab)
        #plot space
        c.df[["exp.diff"]] <- c.diff
        c.df[["cont.all"]] <- c.all
        dif <- c.df$exp.diff
        # dif[which(c.df$exp.diff>0)] <- "up"
        # dif[which(c.df$exp.diff<0)] <- "down"
        dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]>0)] <- "up"
        dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]<0)] <- "down"
        c.df[["diff.lab"]] <- dif
        #build loading contribution
        "eigen.up" = mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
        "eigen.down" <- mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
        "PC1.up" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
        "PC1.down" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
        # uncomment to get all plots
        # png(paste(plot_out,"/eigenpath_",curpath,"_color-",curcomp,"_eigengene.vs.PC1_pts.png",sep=""), res=plot_res, units="in", height=8, width=10)
        # p <- ggplot(c.df, aes(x=eigengenes, y= PC1, color=diff.lab )) + geom_point() + theme_bw() +
        #   geom_segment(aes(x = 0, y = 0, xend = eigen.down, yend = PC1.down), color="#06398c", #old="#1ebd80"
        #                arrow = arrow(length = unit(.65, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        #   geom_segment(aes(x = 0, y = 0, xend = eigen.up, yend = PC1.up), color= "#ab1844", #old= "#c97f10"
        #                arrow = arrow(length = unit(.65, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        #   scale_x_continuous(trans = "reverse") + scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" ))
        # print(p)
        # graphics.off()
        
        ### make "eigen-contribution"
        cml.c <- c.df$eigengenes
        cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 #pro-CML up (neg eig. -> adjust)
        cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 #pro-CML down ( pos eig. -> NO adjust)
        cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 # anti-CML down (neg -> NO adjust)
        cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 # anti-CML up (pos -> adjust)
        
        #sums
        neg.up <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
        neg.up.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
        pos.down <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
        pos.down.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
        neg.down <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
        neg.down.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
        pos.up <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
        pos.up.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
        
        cmat <- matrix( c( rep(curpath, 4), rep(curcomp, 4), "neg|up", "pos|down", "neg|down", "pos|up",
                           neg.up, pos.down, neg.down, pos.up,
                           neg.up.l, pos.down.l, neg.down.l, pos.up.l),
                        nrow=4, ncol=5) 
        cml.tab <- rbind(cml.tab, cmat )
        #full matrix
        # :::note::: the mean contribution of each class of genes must be added!
        n.u.m <- (neg.up  / neg.up.l) 
        p.d.m <- (pos.down / pos.down.l)
        n.d.m <- (neg.down / neg.down.l) 
        p.u.m <- (pos.up / pos.up.l)
        fmat <- matrix( c( rep(curpath, 2), rep(curcomp, 2), "pro-CML", "anti-CML",
                           sum(c(n.u.m, p.d.m), na.rm=T), sum(c(n.d.m, p.u.m), na.rm=T),
                           neg.up.l + pos.down.l, neg.down.l + pos.up.l),
                        nrow=2, ncol=5) 
        cml.full <- rbind(cml.full, fmat )
        
        ### uncomment to get a bunch of plots...
        {
        #counts
        # cl.up <- length(which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
        # cl.down <- length(which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
        # al.up <- length(which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
        # al.down <- length(which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
        # c.df[["CML_contribution"]] <- cml.c
        # cml.up = mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
        # cml.down <- mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
        # PC1.up <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
        # PC1.down <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
        # cml.tot <- mean(cml.c, na.rm=T) 
        # PC1.tot <- mean(c.df$PC1, na.rm=T)
        # #make df for labeling genes without NA genes
        # l.df <- c.df[which(c.df$diff.lab=="up" | c.df$diff.lab=="down" ),]
        # tup <- sort.int(l.df$CML_contribution, index.return = T, decreasing = T)$ix[seq(1,10)]
        # tdown <- sort.int(l.df$CML_contribution, index.return = T, decreasing = F)$ix[seq(1,5)]
        # png(paste(plot_out,"/eigenpath_",curpath,"_color-",curcomp,"_pts+label.png",sep=""), res=plot_res, units="in", height=8, width=10)
        # p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab )) + 
        #   geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", size=1.25) +
        #   geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.5, 
        #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        #   geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.5,
        #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        #   geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.5,
        #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        #   scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
        #   geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene)) +
        #   geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene)) +
        #   ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
        # #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" )) 
        # print(p)
        # graphics.off()
        # png(paste(plot_out,"/eigenpath_",curpath,"_color-",curcomp,"_pts+label+tot.png",sep=""), res=plot_res, units="in", height=8, width=10)
        # p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab )) + 
        #   geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", size=1.25) +
        #   geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.5, 
        #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        #   geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.5,
        #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        #   geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.5,
        #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        #   scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
        #   geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene)) +
        #   geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene)) +
        #   ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
        # #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" )) 
        # print(p)
        # graphics.off()
        # write.table(c.df, paste(plot_out,"/eigenpath_",curpath,"_dataTable.tsv",sep=""), sep="\t", row.names=F) 
        # 
        # ### vs log2FC ###
        # 
        # lfc.up <- mean(c.df[[paste(curcomp,"_lfc",sep="")]][which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
        # lfc.down <- mean(c.df[[paste(curcomp,"_lfc",sep="")]][which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
        # lfc.tot <- mean(c.df[[paste(curcomp,"_lfc",sep="")]], na.rm=T)
        # png(paste(plot_out,"/eigenpath_",curpath,"_vs.log2FC_color-",curcomp,"_pts+label+tot.png",sep=""), res=plot_res, units="in", height=8, width=10)
        # p <- ggplot(c.df, aes(x=CML_contribution, y= .data[[paste(curcomp,"_lfc",sep="")]], color=diff.lab )) + 
        #   geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", linewidth=1.25) +
        #   geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = lfc.down), color="#06398c", size=1.5, 
        #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        #   geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = lfc.up), color= "#ab1844", size=1.5,
        #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        #   geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = lfc.tot), color= "black", size=1.5,
        #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
        #   scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
        #   geom_text(data=l.df[tup,], aes(x=CML_contribution, y=.data[[paste(curcomp,"_lfc",sep="")]]+0.0007, color = diff.lab, label=gene)) +
        #   geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=.data[[paste(curcomp,"_lfc",sep="")]]+.0007, color = diff.lab, label=gene)) +
        #   ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
        # #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" )) 
        # print(p)
        # graphics.off()
        }
        
        
      }
    }#end pathway contribution
  }
  colnames(cont.tab) <- c("pathway","comp","DEGs", "gene_counts", "eigenpath", "contribution", "cont.l2fc")
  cont.df <- data.frame(cont.tab)
  colnames(cml.tab) <- c("pathway", "comp", "dir", "CML.sum", "gene.count")
  cml.df <- data.frame(cml.tab)
  colnames(cml.full) <- c("pathway", "comp", "dir", "CML.mean", "gene.count")
  full.df <- data.frame(cml.full)
  for (curpath in unique(full.df$pathway)){
    for (curcomp in unique(full.df$comp)) {
      curdat <- full.df[which(full.df$pathway==curpath & full.df$comp==curcomp),]
      cdat <- c(curpath, curcomp, "total", sum(as.numeric(curdat$CML.mean)), sum(as.numeric(curdat$gene.count)))
      cml.full <- rbind(cml.full, cdat)
    }
  }
  full.df <- data.frame(cml.full)
  #add mean contribution
  cml.df[["CML.mean"]] <- as.numeric(cml.df$CML.sum) / as.numeric(cml.df$gene.count)
  cml.df$CML.mean[which(is.na(cml.df$CML.mean))] <- 0
  # write.table(cont.df, paste(plot_out,"/eigenpath_all-Pathway_dataTable.tsv",sep=""), sep="\t", row.names=F)
  # write.table(cml.df, paste(plot_out,"/eigenpath_all-Pathway_CML-eigenpath_dataTable.tsv",sep=""), sep="\t", row.names=F)
  # write.table(full.df, paste(plot_out,"/eigenpath_all-Pathway_CML-eigenpath-totalContribution_dataTable.tsv",sep=""), sep="\t", row.names=F)
  
  
  ###
  ### build CML contribution tables: DEGs and GSEA
  ###
  # :::note::: outputs used to construct Table S2
  {
    cctab_out <- "CML_contribution_tables"
    dir.create(cctab_out, showWarnings = F)
    ### DEGS ###
    #read each DEG file, add CML contribution to file
    {
      deg.files <- list.files(path=deg_out, pattern=glob2rx("DESeq-sexFactor_*_DEG-table-means.tsv"), full.names=T )
      for (fpath in deg.files) {
        f <- strsplit(fpath, "/")[[1]][2]
        comp <- strsplit(f,"_")[[1]][2]
        grp1 <- strsplit(comp,".vs.")[[1]][1]
        grp2 <- strsplit(comp,".vs.")[[1]][2]
        glab1 <- gsub("c2", "c3", gsub("c3", "c5", gsub("D-", "TKI-", gsub("A-", "TOTO-", gsub("B-state-", "", gsub( "C", "Ctl", grp1)) ))))
        glab2 <- gsub("c2", "c3", gsub("c3", "c5", gsub("D-", "TKI-", gsub("A-", "TOTO-", gsub("B-state-", "", gsub( "C", "Ctl", grp2)) ))))
        compout <- gsub("c2", "c3", gsub("c3", "c5", comp))
        cur.deg <- read.table(fpath,header=T)
        if (dim(cur.deg)[1]==0) { 
          print(paste("::SKIPPING: ",comp,":::",sep=""))
          next }
        colnames(cur.deg)[which(colnames(cur.deg)=="cmean")] <- paste(glab2,"_mean",sep="") 
        colnames(cur.deg)[which(colnames(cur.deg)=="bmean")] <- paste(glab1,"_mean",sep="")
        #get eigenvalue and CML contribution
        eig <- c()
        cc <- c()
        cclab <- c()
        rev.warn <- F
        for (gi in 1:dim(cur.deg)[1]) {
          g <- cur.deg$gene_name[gi]
          gid <- cur.deg$gene_id[gi]
          lfc <- cur.deg$log2FoldChange[gi]
          m <- which(gene.df$gene_id==gid)
          if ( length(m)>0) {
            if (length(m)>1) {
              print(paste(":::WARNING::: multiple matches for: ",gid,sep=""))
              m <- m[1]
            }
            e <- gene.df$eigengenes[m]
            dim(gene.df)
            length(unique(gene.df$gene_id))
            
            #ensure log2FC has correct sign; CML contribution requires that grp1 is the more CML-like (lower in space) group 
            #flip l2FC sign for the following comparisons
            rev.comps <- c("A-postRx.vs.B-state-c2", "A-postRx.vs.B-state-c3", 
                           "D-Rx.vs.B-state-c3","D-Rx.vs.D-preT-c3")
            if ( comp %in% rev.comps ) {
              if (!rev.warn) {
                print(paste(":::WARNING::: reversing log2FC for ",comp,sep=""))
                rev.warn <- T
              }
              lfc <- -1*lfc
            }
            
            ###pro - adjust to make positive value
            ###anti - adjust to make negative value
            if ( (e < 0) & (lfc > 0)) { #Pro: neg ld; pos lfc
              cml.c <- e*-1 #adj. e pos 
              lab <- "Pro"
            } else if ((e > 0) & (lfc < 0)) { #Pro: pos ld; neg lfc
              cml.c <- e #no adj. e already pos
              lab <- "Pro"
            } else if ((e < 0) & (lfc < 0)) { #Anti: neg ld; neg lfc
              cml.c <- e #no adj; e already negative
              lab <- "Anti"
            } else if ((e > 0) & (lfc > 0)) { #Anti: pos ld; pos lfc
              cml.c <- e*-1 #adj. e pos 
              lab <- "Anti"
            }
            
            
            eig <- c(eig,e)
            cc <- c(cc, cml.c)
            cclab <- c(cclab, lab)
          } else {
            print(paste(":::ERROR::: no match for: ",gid,". Outputting NA",sep=""))
            eig <- c(eig,e)
            cc <- c(cc, cml.c)
            cclab <- c(cclab, lab)
          }
        }
        length(eig)
        new.deg <- cbind(cur.deg, eig, cc, cclab)
        
        colnames(new.deg)[seq(dim(new.deg)[2]-2,dim(new.deg)[2])] <- c("eigenvalue", "eigengene_value", "CML_contribution")
        #write.table(new.deg, paste(cctab_out,"/",gsub(".tsv","_CML-cont.tsv", f),sep=""),sep="\t", row.names=F)
        write.table(new.deg, paste(cctab_out,"/CML.contribution_",compout,"_DEG_table.tsv",sep=""),sep="\t", row.names=F)
      } #end file loop
      
      
      ### output for defined DEG sets - Early, Trans, Late ###
      deg.sets <- list(b0.vs.c, b1.vs.c, b2.vs.c)
      deg.comps <- c("B-state-c1.vs.C","B-state-c2.vs.C","B-state-c3.vs.C")
      sel.genes <- list(early.deg.genes.id, trans.deg.genes.id, late.deg.genes.id)
      set.names <- c("Early", "Transition","Late")
      for (fi in 1:length(deg.sets)) {
        degdat <- deg.sets[[fi]]
        #comp <- strsplit(f,"_")[[1]][2]
        comp <- deg.comps[fi]
        curgenes <- sel.genes[[fi]]
        cursel <- match_genes(degdat$gene_id, curgenes) #get genes that match current gene set
        cur.deg <- degdat[cursel,]
        
        grp1 <- strsplit(comp,".vs.")[[1]][1]
        grp2 <- strsplit(comp,".vs.")[[1]][2]
        glab1 <- gsub("c2", "c3", gsub("c3", "c5", gsub("D-", "TKI-", gsub("A-", "TOTO-", gsub("B-state-", "", gsub( "C", "Ctl", grp1)) ))))
        glab2 <- gsub("c2", "c3", gsub("c3", "c5", gsub("D-", "TKI-", gsub("A-", "TOTO-", gsub("B-state-", "", gsub( "C", "Ctl", grp2)) ))))
        compout <- gsub("c2", "c3", gsub("c3", "c5", comp))
        setname <- set.names[fi]
        
        if (dim(cur.deg)[1]==0) { 
          print(paste("::SKIPPING: ",comp,":::",sep=""))
          next }
        colnames(cur.deg)[which(colnames(cur.deg)=="cmean")] <- paste(glab2,"_mean",sep="") 
        colnames(cur.deg)[which(colnames(cur.deg)=="bmean")] <- paste(glab1,"_mean",sep="")
        #get eigenvalue and CML contribution
        eig <- c()
        cc <- c()
        cclab <- c()
        for (gi in 1:dim(cur.deg)[1]) {
          g <- cur.deg$gene_name[gi]
          gid <- cur.deg$gene_id[gi]
          lfc <- cur.deg$log2FoldChange[gi]
          m <- which(gene.df$gene_id==gid)
          if ( length(m)>0) {
            if (length(m)>1) {
              print(paste(":::WARNING::: multiple matches for: ",gid,sep=""))
              m <- m[1]
            }
            e <- gene.df$eigengenes[m]
            dim(gene.df)
            length(unique(gene.df$gene_id))
            
            
            ###pro - adjust to make positive value
            ###anti - adjust to make negative value
            if ( (e < 0) & (lfc > 0)) { #Pro: neg ld; pos lfc
              cml.c <- e*-1 #adj. e pos 
              lab <- "Pro"
            } else if ((e > 0) & (lfc < 0)) { #Pro: pos ld; neg lfc
              cml.c <- e #no adj. e already pos
              lab <- "Pro"
            } else if ((e < 0) & (lfc < 0)) { #Anti: neg ld; neg lfc
              cml.c <- e #no adj; e already negative
              lab <- "Anti"
            } else if ((e > 0) & (lfc > 0)) { #Anti: pos ld; pos lfc
              cml.c <- e*-1 #adj. e pos 
              lab <- "Anti"
            }
            
            
            eig <- c(eig,e)
            cc <- c(cc, cml.c)
            cclab <- c(cclab, lab)
          } else {
            print(paste(":::ERROR::: no match for: ",gid,". Outputting NA",sep=""))
            eig <- c(eig,e)
            cc <- c(cc, cml.c)
            cclab <- c(cclab, lab)
          }
        }
        length(eig)
        new.deg <- cbind(cur.deg, eig, cc, cclab)
        
        colnames(new.deg)[seq(dim(new.deg)[2]-2,dim(new.deg)[2])] <- c("eigenvalue", "eigengene_value", "CML_contribution")
        ### Table S2 ###
        write.table(new.deg, paste(cctab_out,"/CML.contribution_",setname,"_DEG_table.tsv",sep=""),sep="\t", row.names=F)
      } #end file loop
      
      
    }
    
    
    ### GSEA ###
    #read GSEA file; output contribution column to each GSEA table
    {
      gsea.dirs <- list.files(path=deg_out, pattern=glob2rx("fgsea_sexFactor_*"), include.dirs=T, all.files=T, full.names=T)
      
      # gsea.dirs <- gsea.dirs[which(gsea.dirs==fdir):length(gsea.dirs)] #restart at last file if processing gets interrupted
      for (fdir in gsea.dirs) {
        f <- paste(fdir,"/fgsea_Hallmark_table.tsv",sep="")
        
        comp <- strsplit(strsplit(f,"/")[[1]][2], "_")[[1]][3]
        compout <- gsub("c2", "c3", gsub("c3", "c5", comp))
        curcomp <- gsub("-Rx", "rx", gsub("-postRx", "post", gsub("-state-c","",comp) ) ) #make compatible with "gene.df" names
        grp1 <- strsplit(comp,".vs.")[[1]][1]
        grp2 <- strsplit(comp,".vs.")[[1]][2]
        glab1 <- gsub("c2", "c3", gsub("c3", "c5", grp1))
        glab2 <- gsub("c2", "c3", gsub("c3", "c5", grp2))
        cur.tab <- read.table(f, header=T, sep="\t")
        if (dim(cur.tab)[1]==0) {next} #skip empty
        out.tab <- c()
        for (pi in 1:dim(cur.tab)[1]) {#get CML contribution for each pathway
          curpath <- cur.tab$pathway[pi]
          #leading edge
          le.genes <- strsplit(cur.tab$leadingEdge[pi],"\\|")[[1]]
          sel <- match_all_genes( gene.df$gene, le.genes)
          le.df <- gene.df[sel,]
          cml.le <- le.df$eigengenes
          cml.le[which(le.df$eigengenes < 0 & le.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.le[which(le.df$eigengenes < 0 & le.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 #pro-CML up (neg eig. -> adjust)
          cml.le[which(le.df$eigengenes > 0 & le.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.le[which(le.df$eigengenes > 0 & le.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 #pro-CML down ( pos eig. -> NO adjust)
          cml.le[which(le.df$eigengenes < 0 & le.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.le[which(le.df$eigengenes < 0 & le.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 # anti-CML down (neg -> NO adjust)
          cml.le[which(le.df$eigengenes > 0 & le.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.le[which(le.df$eigengenes > 0 & le.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 # anti-CML up (pos -> adjust)
          le.pos.g <- paste(le.df$gene[which(cml.le>0)],collapse="|")
          le.neg.g <- paste(le.df$gene[which(cml.le<0)],collapse="|")
          le.pos.cnt <- length(cml.le[which(cml.le>0)])
          le.neg.cnt <- length(cml.le[which(cml.le<0)])
          le.pos.mean <- mean(cml.le[which(cml.le>0)])
          le.neg.mean <- mean(cml.le[which(cml.le<0)])
          le.tot <- mean(cml.le)
          
          #all pathway
          cgenes <- geneIds(mm.h[[curpath]])
          sel <- match_all_genes( gene.df$gene, cgenes)
          c.df <- gene.df[sel,]
          cml.c <- c.df$eigengenes
          cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 #pro-CML up (neg eig. -> adjust)
          cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 #pro-CML down ( pos eig. -> NO adjust)
          cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 # anti-CML down (neg -> NO adjust)
          cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 # anti-CML up (pos -> adjust)
          c.pos.g <- paste(c.df$gene[which(cml.c>0)],collapse="|")
          c.neg.g <- paste(c.df$gene[which(cml.c<0)],collapse="|")
          c.pos.cnt <- length(cml.c[which(cml.c>0)])
          c.neg.cnt <- length(cml.c[which(cml.c<0)])
          c.pos.mean <- mean(cml.c[which(cml.c>0)])
          c.neg.mean <- mean(cml.c[which(cml.c<0)])
          c.tot <- mean(cml.c)
          
          #output pathway
          outline <- c(cur.tab[pi,], le.tot, le.pos.mean, le.neg.mean, le.pos.cnt, le.neg.cnt, le.pos.g, le.neg.g,
                       c.tot, c.pos.mean, c.neg.mean, c.pos.cnt, c.neg.cnt, c.pos.g, c.neg.g)
          out.tab <- rbind(out.tab, outline)             
        } #end pathway for loop
        
        ### output gsea table ###
        header <- c(colnames(cur.tab), "leadingEdge.CML_contribution", "leadingEdge.pro-CML_cont", "leadingEdge.anti-CML_cont", 
                    "leadingEdge.pro-CML_gene_count", "leadingEdge.anti-CML_gene_count", "leadingEdge.pro-CML_genes", "leadingEdge.anti-CML_genes",
                    "fullHallmark.CML_contribution", "fullHallmark.pro-CML_cont", "fullHallmark.anti-CML_cont", 
                    "fullHallmark.pro-CML_gene_count", "fullHallmark.anti-CML_gene_count", "fullHallmark.pro-CML_genes", "fullHallmark.anti-CML_genes" ) 
        colnames(out.tab) <- header
        ### Table S2 ###
        write.table(out.tab, paste(cctab_out,"/CML.contribution_",compout,"_GSEA_table.tsv",sep=""),sep="\t", row.names=F)
      }
    }# end GSEA contribution tables

    
  } #end CML contribution tables
  
  
  
}# end loading value analysis



###
### FIGURE 3 PLOTS (Fig. S3-4)
###
{
  ###
  # Figure 3A: critical point classification of trajectories
  ###
  {
    bc.sel <- which(colData(dat.se)$Group=="B" | colData(dat.se)$Group=="C")
    col.df <- data.frame(colData(dat.se))
    #state-space colors
    png(paste(manFig_out,"/Fig-3A_time.vs.CML_space_BC-ONLY_color-state_chr-B-ONLY_line+point.png",sep=""), res=plot_res, units="in", height=3, width=5)
    p <- ggplot(col.df[bc.sel,], aes(x=as.numeric(sample_weeks), y=CML_space, color=Group, group=mouse_id )) + geom_line() +
      scale_color_manual(values=unlist(group.cols) ) +
      geom_point(data=col.df[which(col.df$Group=="B"),], aes(x=as.numeric(sample_weeks), y=CML_space,  shape=state_chr, fill=state_chr, group=mouse_id, stroke=NA), size=3) +
      theme_bw(base_size=18) + geom_hline(yintercept=ss.cps[c(2,4)], color="dimgrey", linetype="dashed", alpha=.75) +
      scale_shape_manual(values=c(21,22,24)) + scale_fill_manual(values=state_rx.cols) + theme(legend.position="none")
    print(p)
    graphics.off()
  }
  
  ###
  # Figure 3B: DEG upset
  ###
  {
    ## CML vs ctl & CML ##
    degs.cml6 <- list("c5.vs.ctl"=b2.vs.c[b2vc.degs, ]$gene_id, "c3.vs.ctl"=b1.vs.c[b1vc.degs, ]$gene_id, 
                      "c1.vs.ctl"=b0.vs.c[b0vc.degs, ]$gene_id, "c3.vs.c1"=b1.vs.b0[b1v0.degs, ]$gene_id,
                      "c5.vs.c3"=b2.vs.b1[b2v1.degs, ]$gene_id,"c5.vs.c0"=b2.vs.b0[b2v0.degs, ]$gene_id)
    cm.cml6 = ComplexHeatmap::make_comb_mat(degs.cml6)
    cm_annot <- function(cm.cml6) {
      comb.annot <- c()
      for (ni in 1:length(comb_name(cm.cml6)) ){ #add unique
        n <- comb_name(cm.cml6)[ni]
        n.list <- as.numeric(strsplit(n, "")[[1]])
        #name order: "c5", "c3", "c1", "early", "late", "full"
        if (sum(n.list[1:3])>1 & sum(n.list[4:6])==0) { #vs ctl
          comb.annot <- c(comb.annot, "vs Ctl")
        } else if (sum(n.list[1:3])==0 & sum(n.list[4:6])>1) { # inter CML
          comb.annot <- c(comb.annot, "inter-CML")
        } else if (sum(n.list)==1) { # unique
          comb.annot <- c(comb.annot, "unique")
        } else {
          comb.annot <- c(comb.annot, "mixed")
        }
      }
      return(comb.annot)
    }
    cm.cml6.sel <- cm.cml6[comb_size(cm.cml6)>20]
    comb.annot <- cm_annot(cm.cml6.sel)
    comb.cols <- c( "vs Ctl" = "red", "inter-CML" = "#5000a1", "mixed"="grey80", "unique"="black")
    comb.sort <- sort.int(comb_size(cm.cml6.sel), index.return = T)$ix
    png(paste(manFig_out,"/Fig-2B_vsControl+interCML_upset-uniq.png",sep=""), res=plot_res, units="in", height=3, width=5)
    UpSet(cm.cml6.sel, set_order = c("c1.vs.ctl", "c3.vs.ctl", "c5.vs.ctl", "c3.vs.c1",  "c5.vs.c3",  "c5.vs.c0"),
          comb_order = rev(comb.sort),
          right_annotation = upset_right_annotation(cm.cml6.sel, add_numbers = TRUE,
                                                    gp = gpar(fill = c("#fa1b2e", "#fab71b", "#1b2afa", "#a9f5b2", "#facba5",  "#faa5c0" )) ),
          top_annotation = HeatmapAnnotation(
            Treatment = comb.annot,
            col = list(Treatment=comb.cols),
            "Intersection\nsize" = anno_barplot(comb_size(cm.cml6.sel), 
                                                border = FALSE, 
                                                gp = gpar(fill = "black"), 
                                                height = unit(2, "cm")
            ), 
            annotation_name_side = "left", 
            annotation_name_rot = 0) )
    
    dev.off()
  }
  
  ###
  # Figure 3C: GSEA heatmap
  ###
  {
    plimit <- 0.0001
    c1sig <- which(f.b1.c$padj <= plimit)
    c2sig <- which(f.b2.c$padj <= plimit)
    c3sig <- which(f.b3.c$padj <= plimit)
    f.b1.c[c1sig,]
    f.b2.c[c2sig,]
    f.b3.c[c3sig,]
    
    min.padj <- min(c(f.b1.c[c1sig,]$padj,
                      f.b2.c[c2sig,]$padj,
                      f.b3.c[c3sig,]$padj))
    ### ALT: dotplots ###
    {
      # plot_sig_fGSEA_man <- function(cur_fgsea, bounds = c(-2,2), outname="unnamed_dotplot.png", ht=8, sigpval = 0.05, minNES = 0) {
      #   sigsel <- which(cur_fgsea[["padj"]]<sigpval & abs(cur_fgsea[["NES"]]) > minNES)
      #   sigdf <- data.frame("Pathway"= unlist(lapply(cur_fgsea[["pathway"]][sigsel], function(x) gsub("_", " ", x) )),
      #                       "qvalue"=cur_fgsea[["padj"]][sigsel], "Count"=cur_fgsea[["size"]][sigsel], "NES"=cur_fgsea[["NES"]][sigsel])
      #   sigdf$type <- "Upregulated"
      #   sigdf$type[sigdf$NES < 0 ] <- "Downregulated"
      #   # format pathways names
      #   fpath <- unlist(lapply(sigdf$Pathway, function(x) str_to_title( gsub("HALLMARK ","", x) ) ))
      #   sigdf$Pathway <- fpath
      #   sigdf <- sigdf %>% arrange(factor(Pathway, levels = sigdf[["Pathway"]][order(sigdf[["qvalue"]], decreasing=F)]))
      #   #ggplot(data=sigdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=16)
      #   #orddf <- sigdf %>% mutate(Pathway = fct_reorder(Pathway, qvalue, .desc=T)) #idk what this does...
      #   plotdf <- sigdf
      #   #if (dim(plotdf)[1]>20) { plotdf <- plotdf[seq(1,20),]} #limit output if too large
      #   # #resize if too many pathways
      #   # ht <- 8
      #   # if (dim(plotdf)[1]>20 & dim(plotdf)[1]>30) { 
      #   #   ht <- 12
      #   # } else {
      #   #   ht <- 18
      #   # }
      #   plotdf$Pathway <- factor(plotdf$Pathway, levels=rev(plotdf$Pathway))
      #   print(paste("Outputting to: ",outname," with dims: ", paste(dim(plotdf),collapse=" x "),sep=""))
      #   png(outname, res=plot_res, units="in", height=ht, width=6)
      #   p <- ggplot(plotdf, aes(x=NES, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=16) + 
      #     scale_color_gradient(limits=c(0,sigpval), low="red", high="blue") +
      #     coord_cartesian(xlim=bounds) + theme(legend.position="none")
      #   print(p)
      #   dev.off()
      # }
      # plot_sig_fGSEA_man(f.b1.c, bounds = c(-2.8,2.8), 
      #                outname=paste(plot_out,"/fgsea-summary_NES_B1.vs.C_allSig_p-0.0001_dotplot.png",sep=""), 
      #                sigpval = 0.0001, minNES = 2,
      #                ht=4)
      # plot_sig_fGSEA_man(f.b2.c, bounds = c(-2.8,2.8), 
      #                outname=paste(plot_out,"/fgsea-summary_NES_B2.vs.C_allSig_p-0.0001_dotplot.png",sep=""), 
      #                sigpval = 0.0001, minNES = 2,
      #                ht=2)
      # plot_sig_fGSEA_man(f.b3.c, bounds = c(-2.8,2.8), 
      #                outname=paste(plot_out,"/fgsea-summary_NES_B3.vs.C_allSig_p-0.0001_dotplot.png",sep=""), 
      #                sigpval = 0.0001, minNES = 2,
      #                ht=5)
    }
    ## heatmap - all sig ##
    colnames(f.b1.c)
    
    length(c1sig)
    tmp1 <- merge(f.b1.c[c1sig,c("pathway","NES")], f.b2.c[c2sig,c("pathway","NES")],  by="pathway", suffix=c("c1","c3") )
    df.list <- list(f.b1.c[c1sig,c("pathway","NES")], f.b2.c[,c("pathway","NES")], f.b3.c[,c("pathway","NES")] )
    
    
    inc <- 50
    #mval <- max(abs(unlist(sig.nes.nz)))
    mval <- max(abs(c(f.b1.c$NES[c1sig], f.b2.c$NES[c2sig], f.b3.c$NES[c3sig])))
    blist <- seq(-1*mval, mval, length.out=inc+1)
    corcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(inc)
    all.sig <-  sort(unique(c(f.b1.c[c1sig,c("pathway")], f.b2.c[c2sig,c("pathway")], f.b3.c[c3sig,c("pathway")] )))
    df.list <- list(f.b1.c[c1sig,], f.b2.c[c2sig,], f.b3.c[c3sig,] )
    df.name <- c("c1.vs.ctl", "c3.vs.ctl", "c5.vs.ctl")
    sig.nes <- c()
    for (df in df.list) {
      out <- c()
      for (p in all.sig) {
        m <- which(df$pathway==p)
        if (length(m)==0) {
          out <- c(out,NA)
        } else {
          out <- c(out, df[["NES"]][m])
        }
      }
      sig.nes <- cbind(sig.nes, out)
    }
    colnames(sig.nes) <- df.name
    rownames(sig.nes) <- all.sig
    sig.nes.nz <- sig.nes
    sig.nes.nz[which(is.na(sig.nes))] <- 0
    #?hclust(1-cor(t(sig.nes.nz)))
    n.ph <- pheatmap::pheatmap(t(sig.nes.nz))
    col.ord <- n.ph$tree_col$order
    dim(sig.nes.nz)
    col.ord <- order(sig.nes[,3], decreasing=T)
    all.sig.name <- unlist(lapply(all.sig, function(x) gsub("_", " ", gsub("HALLMARK_", "",x ))  ))
    nes.range <- seq(-1*max(abs(sig.nes),na.rm=T), max(abs(sig.nes),na.rm=T), length.out=inc+1)
    png(paste(manFig_out,"/Fig-2C_fgsea-summary_NES_all-Significant_p-",plimit,"_pheatmap.png",sep=""), res=plot_res, units="in", height=8, width=6)
    p <- pheatmap::pheatmap(sig.nes[col.ord,], scale="none", cluster_rows = F, cluster_cols = F, breaks=nes.range, color=corcol , 
             labels_row=str_to_title(all.sig.name[col.ord]), angle_col="315", na_col="dimgrey", fontsize=18)
    print(p)
    graphics.off()  

  
  }
  
  ###
  # Figure 3D: CML contribution
  ###
  {
    amin <- minabd
    gene_set_CML_contribution_manFig <- function( out_name, deg_list, deg_comp, gene.df = gene.df ) { #all
      
      curdeg <- out_name #deg_name[[i]] #NEEDED:change curpath to curdeg
      cgenes.id <- deg_list
      curcomp <- deg_comp
      
      sel <- match_genes( gene.df$gene_id, cgenes.id)
      c.df <- gene.df[sel,]
      
      m.1 <- c.df[[paste(curcomp,"_mean.",1,sep="")]]
      m.2 <- c.df[[paste(curcomp,"_mean.",2,sep="")]]
      c.diff <- m.1 - m.2
      c.ld <- c.df$eigengenes
      up.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]
      dn.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]
      c.all <- c.diff * c.ld #get contribution for all genes in pathway
      c.up <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
      c.down <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
      ctab <- matrix( c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", sum(c.ld), sum(up.ld), sum(dn.ld), sum(c.all,na.rm=T), sum(c.up), sum(c.down) ), nrow=3, ncol=5) 
      #add lfc contribution; need to figure out why "Inf" is produced
      c.lfc <- log2((m.1+amin)/(m.2+amin) )
      cl.all <- c.ld * c.lfc #get contribution for all genes in pathway based on LFC
      cl.up <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
      cl.down <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
      ctab <- matrix(c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", 
                        length(which(!is.na(c.all))), length(which(!is.na(c.up))), length(which(!is.na(c.down))), # count of genes contributing (non NA genes)
                        sum(c.ld), sum(up.ld), sum(dn.ld), #eigenvalue
                        sum(c.all,na.rm=T), sum(c.up, na.rm=T), sum(c.down, na.rm=T), #contribution w/ exp diff
                        sum(cl.all,na.rm=T), sum(cl.up, na.rm=T), sum(cl.down, na.rm=T) ), # contrib. w/ log2FC
                     nrow=3, ncol=7)
      if (length(which(is.na(ctab))) != 0) {
        print(paste(":::ERROR::: NAs produced for comp: ",curcomp," pathway: ",curdeg,sep=""))
        print(ctab)
        stop("NA encountered")
      }
      cont.tab <- rbind(cont.tab, ctab)
      #plot space
      c.df[["exp.diff"]] <- c.diff
      c.df[["cont.all"]] <- c.all
      dif <- c.df$exp.diff
      # dif[which(c.df$exp.diff>0)] <- "up"
      # dif[which(c.df$exp.diff<0)] <- "down"
      dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]>0)] <- "up"
      dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]<0)] <- "down"
      c.df[["diff.lab"]] <- dif
      #build loading contribution
      "eigen.up" = mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      "eigen.down" <- mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
      "PC1.up" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      "PC1.down" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      
      ### make "eigen-contribution"
      # positive value = positive contribution (neg eig + pos exp OR pos eig + neg exp)
      # neg contribution (neg eig + neg exp OR pos eig + pos exp)
      cml.c <- c.df$eigengenes
      cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 #pro-CML up (neg eig. -> adjust)
      cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 #pro-CML down ( pos eig. -> NO adjust)
      cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 # anti-CML down (neg -> NO adjust)
      cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 # anti-CML up (pos -> adjust)
      
      #sums
      neg.up <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
      neg.up.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
      pos.down <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
      pos.down.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
      neg.down <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
      neg.down.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
      pos.up <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
      pos.up.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
      cmat <- matrix( c( rep(curdeg, 4), rep(curcomp, 4), "neg|up", "pos|down", "neg|down", "pos|up", 
                         neg.up, pos.down, neg.down, pos.up,
                         neg.up.l, pos.down.l, neg.down.l, pos.up.l),
                      nrow=4, ncol=5) 
      cml.tab <- rbind(cml.tab, cmat )
      #cml label
      cml.lab <- rep(NA, length(cml.c))
      cml.lab[which(cml.c > 0)] <- "pro-CML"
      cml.lab[which(cml.c < 0)] <- "anti-CML"
      #full matrix
      # :::note::: the mean contribution of each class of genes must be added!
      #(neg.up  / neg.up.l) + (pos.down / pos.down.l), (neg.down / neg.down.l) + (pos.up / pos.up.l)
      n.u.m <- (neg.up  / neg.up.l) 
      p.d.m <- (pos.down / pos.down.l)
      n.d.m <- (neg.down / neg.down.l) 
      p.u.m <- (pos.up / pos.up.l)
      fmat <- matrix( c( rep(curdeg, 2), rep(curcomp, 2), "pro-CML", "anti-CML",
                         sum(c(n.u.m, p.d.m), na.rm=T), sum(c(n.d.m, p.u.m), na.rm=T),
                         neg.up.l + pos.down.l, neg.down.l + pos.up.l),
                      nrow=2, ncol=5) 
      cml.full <- rbind(cml.full, fmat )
      ### eigengene value vs PC1 ###
      cl.up <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
      cl.down <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
      al.up <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
      al.down <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
      c.df[["CML_contribution"]] <- cml.c
      cml.up = mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      cml.down <- mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
      PC1.up <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      PC1.down <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      cml.tot <- mean(cml.c, na.rm=T) 
      PC1.tot <- mean(c.df$PC1, na.rm=T)
      #make df for labeling genes without NA genes
      l.df <- c.df[which(c.df$diff.lab=="up" | c.df$diff.lab=="down" ),]
      tup <- sort.int(l.df$CML_contribution, index.return = T, decreasing = T)$ix[seq(1,10)]
      tdown <- sort.int(l.df$CML_contribution, index.return = T, decreasing = F)$ix[seq(1,5)]
      png(paste(manFig_out,"/Fig-2D_eigenDEGs_",curdeg,"_comparison-",curcomp,"_pts+label.png",sep=""), res=plot_res, units="in", height=4, width=4)
      p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab, fill=diff.lab )) + 
        geom_point(aes(shape=diff.lab), alpha=.75) + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", linewidth=1) +
        # geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", linewidth=1.2, 
        #              arrow = arrow(length = unit(.25, "cm"), type="open"), lineend="butt", linejoin="mitre") +
        # geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", linewidth=1.2,
        #              arrow = arrow(length = unit(.25, "cm"), type="open"), lineend="butt", linejoin="mitre") +
        geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", linewidth=1.25,
                     arrow = arrow(length = unit(.30, "cm"), type="closed"), lineend="butt", linejoin="mitre") +
        scale_color_manual(values= list("up"="#f76a95", "down" = "#4287f5" )) + scale_fill_manual(values= list("up"="#f76a95", "down" = "#4287f5" )) +
        # geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene), linewidth=14) +
        # geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene), linewidth=14) +
        theme(legend.position = "none") + scale_shape_manual(values=c("up"=24, "down"=25)) +
        scale_x_continuous(limits=c(-.03, 0.058)) + scale_y_continuous(limits=c(-.025, 0.018)) 
      #ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
      #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" )) 
      print(p)
      graphics.off()
      
      
      
      
    } #end gene set contribution function
    
    deg_list <- list( early.deg.genes.id, trans.deg.genes.id, late.deg.genes.id, 
                      ig.b1.c, ig.b2.c, ig.b3.c,
                      ig.b3.b1, ig.b3.b2, ig.b2.b1,
                      ig.apost.c, ig.dpost.c, ig.drx.c, ig.c.evl, ig.c.lve, ig.apost.cL, ig.rx.da)
    deg_name <- list("Early",  "Transition", "Late", 
                     "B1.vs.C", "B2.vs.C", "B3.vs.C",
                     "B3.vs.B1", "B3.vs.B2", "B2.vs.B1",
                     "Apost.vs.C", "Dpost.vs.C", "Drx.vs.C", "C_early.vs.late", "C_late.vs.early", "Apost.vs.Clate", "Drx.vs.Apost")
    deg_comp <- list( "B1.vs.C", "B2.vs.C", "B3.vs.C", 
                      "B1.vs.C", "B2.vs.C", "B3.vs.C",
                      "B3.vs.B1", "B3.vs.B2", "B2.vs.B1",
                      "Apost.vs.C", "Dpost.vs.C", "Drx.vs.C","C_early.vs.late", "C_late.vs.early", "Apost.vs.Clate", "Drx.vs.Apost")
    for (i in 4:6){
      i<-4
      gene_set_CML_contribution_manFig( deg_name[i], deg_list[[i]], deg_comp[i], gene.df=gene.df)
    }
  }
  
  
  ###
  # Figure 3E: Gene set CML contribution
  ###
  {
    gsea.b1.c <- read.table("CML_contribution_tables/CML.contribution_B-state-c1.vs.C_GSEA_table.tsv", sep="\t", header=T)
    gsea.b2.c <- read.table("CML_contribution_tables/CML.contribution_B-state-c3.vs.C_GSEA_table.tsv", sep="\t", header=T)
    gsea.b3.c <- read.table("CML_contribution_tables/CML.contribution_B-state-c5.vs.C_GSEA_table.tsv", sep="\t", header=T)
    plimit <- 0.0001
    sel.b1.c <- which(gsea.b1.c$padj < plimit)
    sel.b2.c <- which(gsea.b2.c$padj < plimit)
    sel.b3.c <- which(gsea.b3.c$padj < plimit)
    
    # color by state-space #
    {
      all.sig <-  sort(unique(c(gsea.b1.c[sel.b1.c,c("pathway")], gsea.b2.c[sel.b2.c,c("pathway")], gsea.b3.c[sel.b3.c,c("pathway")] )))
      df.list <- list(gsea.b1.c[sel.b1.c,], gsea.b2.c[sel.b2.c,], gsea.b3.c[sel.b3.c,] )
      df.name <- c("c1.vs.ctl", "c3.vs.ctl", "c5.vs.ctl")
      sig.nes <- c()
      for (df in df.list) {
        out <- c()
        for (p in all.sig) {
          m <- which(df$pathway==p)
          if (length(m)==0) {
            out <- c(out,NA)
          } else {
            out <- c(out, df[["leadingEdge.CML_contribution"]][m])
          }
        }
        sig.nes <- cbind(sig.nes, out)
      }
      colnames(sig.nes) <- df.name
      rownames(sig.nes) <- all.sig
      sig.df <- data.frame(sig.nes)
      dim(sig.nes)
      
      state.df <- data.frame("pathway" = rep(rownames(sig.df),3), "contribution" = c( sig.df$c1.vs.ctl, sig.df$c3.vs.ctl, sig.df$c5.vs.ctl ), 
                             "state" = c(rep("c1", dim(sig.df)[1]),rep("c2", dim(sig.df)[1]), rep("c3", dim(sig.df)[1])  ) )
      state.df$pathway <-  unlist(lapply(state.df$pathway, function(x) str_to_title(gsub("_"," ", gsub( "HALLMARK_", "", x)) )) )
      state.df$state <- factor(state.df$state, levels=c("c1", "c2", "c3"))
      # line_positions holds divider line between each variable
      state.df2 <- state.df %>% mutate(line_positions = as.numeric(factor(pathway, levels = unique(pathway))), 
             line_positions = line_positions + .5,  
             line_positions = ifelse(line_positions == max(line_positions), NA, line_positions)) 
      png(paste(manFig_out,"/Fig-2E_fgsea-CML-contribution_p-",plimit,"_bar.png",sep=""), res=plot_res, units="in", height=3, width=12)
      p <- ggplot(state.df2, aes(x=pathway, y=contribution, fill=state)) + 
        geom_bar(stat = "identity", position=position_dodge(width = 0.5), width=.5) + 
        theme_bw(base_size=12) + scale_fill_manual(values=state_rx.cols) +
        theme(axis.text.x = element_text(angle = 60,  hjust=1), legend.position="none", 
              panel.grid.major.x = element_blank(), panel.grid.minor.x = element_line(colour = "black") ) +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
        geom_vline(aes(xintercept = line_positions))
      print(p)
      graphics.off()
    }
  }
  
  ### SUPPLEMENTAL ###
  
  ###
  # Figure S3 - upset plots
  ###
  {
    ## CML vs ctl & CML ##
    cm_annot <- function(cm.cml6) {
      comb.annot <- c()
      for (ni in 1:length(comb_name(cm.cml6)) ){ #add unique
        n <- comb_name(cm.cml6)[ni]
        n.list <- as.numeric(strsplit(n, "")[[1]])
        #name order: "c5", "c3", "c1", "early", "late", "full"
        if (sum(n.list[1:3])>1 & sum(n.list[4:6])==0) { #vs ctl
          comb.annot <- c(comb.annot, "vs Ctl")
        } else if (sum(n.list[1:3])==0 & sum(n.list[4:6])>1) { # inter CML
          comb.annot <- c(comb.annot, "inter-CML")
        } else if (sum(n.list)==1) { # unique
          comb.annot <- c(comb.annot, "unique")
        } else {
          comb.annot <- c(comb.annot, "mixed")
        }
      }
      return(comb.annot)
    }
    comb.annot <- cm_annot(cm.cml6)
    comb.cols <- c( "vs Ctl" = "red", "inter-CML" = "#5000a1", "mixed"="grey80", "unique"="black")
    comb.sort <- sort.int(comb_size(cm.cml6), index.return = T)$ix
    png(paste(manFig_out,"/Fig-S2_vsControl+interCML_upset-uniq.png",sep=""), res=plot_res, units="in", height=3, width=6)
    UpSet(cm.cml6, set_order = c("c1.vs.ctl", "c3.vs.ctl", "c5.vs.ctl", "c3.vs.c1",  "c5.vs.c3",  "c5.vs.c0"),
          comb_order = rev(comb.sort),
          right_annotation = upset_right_annotation(cm.cml6, add_numbers = TRUE,
                                                    gp = gpar(fill = c("#fa1b2e", "#fab71b", "#1b2afa", "#a9f5b2", "#facba5",  "#faa5c0" )) ),
          top_annotation = upset_top_annotation(cm.cml6, add_numbers = TRUE) )
    
    dev.off()
    # UP plot #
    cm.cml6.up = ComplexHeatmap::make_comb_mat(degs.cml6.up)
    comb.cols <- c( "vs Ctl" = "red", "inter-CML" = "#5000a1", "mixed"="grey80", "unique"="black")
    comb.sort <- sort.int(comb_size(cm.cml6.up), index.return = T)$ix
    comb.annot <- cm_annot(cm.cml6.up)
    png(paste(manFig_out,"/Fig-S2_vsControl+interCML_UPREGULATED_upset-uniq.png",sep=""), res=plot_res, units="in", height=3, width=6)
    UpSet(cm.cml6.up, set_order = c("c1.vs.ctl", "c3.vs.ctl", "c5.vs.ctl", "c3.vs.c1",  "c5.vs.c3",  "c5.vs.c0"),
          comb_order = rev(comb.sort),
          right_annotation = upset_right_annotation(cm.cml6.up, add_numbers = TRUE,
                                                    gp = gpar(fill = c("#fa1b2e", "#fab71b", "#1b2afa", "#a9f5b2", "#facba5",  "#faa5c0" )) ),
          top_annotation = HeatmapAnnotation(
            Comparison = comb.annot,
            col = list(Comparison=comb.cols),
            "Intersection\nsize" = anno_barplot(comb_size(cm.cml6.up), 
                                               border = FALSE, 
                                               gp = gpar(fill = "#8c060a"), 
                                               height = unit(2, "cm")
            ), 
            annotation_name_side = "left", 
            annotation_name_rot = 0) )
    
    dev.off()
    # DOWN plot #
    cm.cml6.dn = ComplexHeatmap::make_comb_mat(degs.cml6.dn)
    comb.cols <- c( "vs Ctl" = "red", "inter-CML" = "#5000a1", "mixed"="grey80", "unique"="black")
    comb.sort <- sort.int(comb_size(cm.cml6.dn), index.return = T)$ix
    comb.annot <- cm_annot(cm.cml6.dn)
    png(paste(manFig_out,"/Fig-S2_vsControl+interCML_DOWNREGULATED_upset-uniq.png",sep=""), res=plot_res, units="in", height=3, width=6)
    UpSet(cm.cml6.dn, set_order = c("c1.vs.ctl", "c3.vs.ctl", "c5.vs.ctl", "c3.vs.c1",  "c5.vs.c3",  "c5.vs.c0"),
          comb_order = rev(comb.sort),
          right_annotation = upset_right_annotation(cm.cml6.dn, add_numbers = TRUE,
                                                    gp = gpar(fill = c("#fa1b2e", "#fab71b", "#1b2afa", "#a9f5b2", "#facba5",  "#faa5c0" )) ),
          top_annotation = HeatmapAnnotation(
            Comparison = comb.annot,
            col = list(Comparison=comb.cols),
            "Intersection\nsize" = anno_barplot(comb_size(cm.cml6.dn), 
                                                 border = FALSE, 
                                                 gp = gpar(fill = "#6d99c2"), 
                                                 height = unit(2, "cm")
            ), 
            annotation_name_side = "left", 
            annotation_name_rot = 0) )
    
    dev.off()
  }
  
  ###
  # Figure S3B - inter-CML dotplots
  ###
  {
    plot_sig_fGSEA_S3B <- function(cur_fgsea, bounds = c(-2,2), outname="unnamed_dotplot.png", ht=8, sigpval = 0.05) {
      sigsel <- which(cur_fgsea[["padj"]]<sigpval)
      sigdf <- data.frame("Pathway"= unlist(lapply(cur_fgsea[["pathway"]][sigsel], function(x) gsub("_", " ", x) )),
                          "qvalue"=cur_fgsea[["padj"]][sigsel], "Count"=cur_fgsea[["size"]][sigsel], "NES"=cur_fgsea[["NES"]][sigsel])
      sigdf$type <- "Upregulated"
      sigdf$type[sigdf$NES < 0 ] <- "Downregulated"
      # format pathways names
      fpath <- unlist(lapply(sigdf$Pathway, function(x) str_to_title( gsub("HALLMARK ","", x) ) ))
      sigdf$Pathway <- fpath
      sigdf <- sigdf %>% arrange(factor(Pathway, levels = sigdf[["Pathway"]][order(sigdf[["qvalue"]], decreasing=F)]))
      #ggplot(data=sigdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=16)
      #orddf <- sigdf %>% mutate(Pathway = fct_reorder(Pathway, qvalue, .desc=T)) #idk what this does...
      plotdf <- sigdf
      #if (dim(plotdf)[1]>20) { plotdf <- plotdf[seq(1,20),]} #limit output if too large
      # #resize if too many pathways
      # ht <- 8
      # if (dim(plotdf)[1]>20 & dim(plotdf)[1]>30) { 
      #   ht <- 12
      # } else {
      #   ht <- 18
      # }
      plotdf$Pathway <- factor(plotdf$Pathway, levels=rev(plotdf$Pathway))
      print(paste("Outputting to: ",outname," with dims: ", paste(dim(plotdf),collapse=" x "),sep=""))
      png(outname, res=plot_res, units="in", height=ht, width=10)
      p <- ggplot(plotdf, aes(x=NES, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=24) + 
        scale_color_gradient(limits=c(0,sigpval), low="red", high="blue") 
      coord_cartesian(xlim=bounds)
      print(p)
      dev.off()
    }
    
    
    c2.c1.sig <- which(f.b1.b0$padj <= 0.001)
    c3.c2.sig <- which(f.b2.b1$padj <= 0.001)
    c3.c1.sig <- which(f.b2.b0$padj <= 0.001)
    sig.b1.b0 <- f.b1.b0[c2.c1.sig,]
    sig.b2.b1 <- f.b2.b1[c3.c2.sig,]
    sig.b2.b0 <- f.b2.b0[c3.c1.sig,]
    
    min.padj <- min(c(f.b1.b0[c2.c1.sig,]$padj,
                      f.b2.b1[c3.c2.sig,]$padj,
                      f.b2.b0[c3.c1.sig,]$padj))
    length(c3sig)
    
    dim(f.b2.b1)
    f.b2.b1[seq(1,25), c(1,2,3,4,5,6)]
    sig.b2.b1[which(sig.b2.b1$NES<0),c(1,2,3,4,5)]
    sig.b2.b0[which(sig.b2.b0$NES<0),c(1,2,3,4,5)]
    
    plot_sig_fGSEA_S3B(sig.b1.b0, bounds = c(-2.8,2.8), 
                   outname=paste(manFig_out,"/Fig-S3B_fgsea-summary_NES_B2.vs.B1_allSig_p-0.001_dotplot.png",sep=""), 
                   ht=length(c2.c1.sig)/2)
    plot_sig_fGSEA_S3B(sig.b2.b1, bounds = c(-2.8,2.8), 
                   outname=paste(manFig_out,"/Fig-S3B_fgsea-summary_NES_B3.vs.B2_allSig_p-0.001_dotplot.png",sep=""), 
                   ht=length(c3.c2.sig)/2)
    plot_sig_fGSEA_S3B(sig.b2.b0, bounds = c(-2.8,2.8), 
                   outname=paste(manFig_out,"/Fig-S3B_fgsea-summary_NES_B3.vs.B1_allSig_p-0.001_dotplot.png",sep=""), 
                   ht=length(c3.c1.sig)/2)
  }
  
  ###
  # Figure S3 - Early, Transition, Late unique DEGs
  ###
  {
    ### venns
    plotVenn( list(b0.vs.c[b0vc.degs, ]$gene_name, b1.vs.b0[b1v0.degs, ]$gene_name, b2.vs.b1[b2v1.degs, ]$gene_name), 
              sNames=c("Early", "c3.vs.c1","c5.vs.c3"), showPlot = T, 
               labelRegions=F, fontScale = 2.5, outFile=paste(manFig_out,"/Fig-S2_vsControl_plot_early-DEGs-ryb.svg",sep=""),
              setColors=c("dodgerblue", "grey40", "grey"))
    rsvg::rsvg_png(
      paste(manFig_out,"/Fig-S2_vsControl_plot_early-DEGs-ryb.svg",sep=""), paste(manFig_out,"/Fig-S2_vsControl_plot_early-DEGs-ryb.png",sep=""),
      width = 900, height = 700)
    
    plotVenn( list(b1.vs.c[b1vc.degs, ]$gene_name, b0.vs.c[b0vc.degs, ]$gene_name, b2.vs.b1[b2v1.degs, ]$gene_name), 
              sNames=c("Transition", "c1.vs.Ctrl","c5.vs.c3"), showPlot = T, fontScale = 2.5,
              outFile=paste(manFig_out,"/Fig-S2_vsControl_plot_transition-DEGs-ryb.svg",sep=""), labelRegions=F,
              setColors=c("goldenrod", "grey40", "grey"))
    rsvg::rsvg_png(
      paste(manFig_out,"/Fig-S2_vsControl_plot_transition-DEGs-ryb.svg",sep=""), paste(manFig_out,"/Fig-S2_vsControl_plot_transition-DEGs-ryb.png",sep=""),
      width = 900, height = 700)
    
    plotVenn( list(b2.vs.c[b2vc.degs, ]$gene_name, b1.vs.c[b1vc.degs, ]$gene_name, b0.vs.c[b0vc.degs, ]$gene_name), 
              sNames=c("c5.vs.Ctrl", "c3.vs.Ctrl","c1.vs.Ctrl"), showPlot = T, fontScale = 2.5,
              outFile=paste(manFig_out,"/Fig-S2_vsControl_plot_Late-DEGs-ryb.svg",sep=""), labelRegions=F,
              setColors=c("red", "grey40", "grey"))
    rsvg::rsvg_png(
      paste(manFig_out,"/Fig-S2_vsControl_plot_late-DEGs-ryb.svg",sep=""), paste(manFig_out,"/Fig-S2_vsControl_plot_late-DEGs-ryb.png",sep=""),
      width = 900, height = 700)
    
    ### direction boxplots
    plot_deg_dir_man <- function(degs.up, degs.down, outname) {
      cdf <- data.frame("Count"=c(length(degs.up),length(degs.down) ), "Direction"=c("Up", "Down"))
      png(outname, res=plot_res, units="in", height=3, width=3 )
      p <- ggplot(cdf, aes(x=Direction, y=Count, fill=Direction)) + geom_bar(stat="identity", width=.5) +
        theme_bw(base_size=18) + theme(legend.position="none") + scale_x_discrete(name="") + 
        scale_fill_manual(values=list("Up"="red", "Down"="dodgerblue")) + scale_y_continuous(name="DEGs")
      print(p)
      graphics.off()
    }
    
    sel.early <- match_deg_genes( early.deg.genes, b0.vs.c)
    early.degs <- b0.vs.c[sel.early,]
    degs.up <- early.degs$gene_name[which(early.degs$log2FoldChange > 2 & early.degs$padj < 0.05)]
    degs.down <- early.degs$gene_name[which(early.degs$log2FoldChange < -2 & early.degs$padj < 0.05)]
    plot_deg_dir_man(degs.up, degs.down, paste(manFig_out,"/Fig-S2_DEGs-def_Early_barplot.png", sep="") )
    
    sel.trans <- match_deg_genes( trans.deg.genes, b1.vs.c)
    trans.degs <- b1.vs.c[sel.trans,]
    degs.up <- trans.degs$gene_name[which(trans.degs$log2FoldChange > 2 & trans.degs$padj < 0.05)]
    degs.down <- trans.degs$gene_name[which(trans.degs$log2FoldChange < -2 & trans.degs$padj < 0.05)]
    plot_deg_dir_man(degs.up, degs.down, paste(manFig_out,"/Fig-S2_DEGs-def_Transition_barplot.png", sep="") )
    
    sel.late <- match_deg_genes( late.deg.genes, b2.vs.c)
    late.degs <- b2.vs.c[sel.late,]
    degs.up <- late.degs$gene_name[which(late.degs$log2FoldChange > 2 & late.degs$padj < 0.05)]
    degs.down <- late.degs$gene_name[which(late.degs$log2FoldChange < -2 & late.degs$padj < 0.05)]
    plot_deg_dir_man(degs.up, degs.down, paste(manFig_out,"/Fig-S2_DEGs-def_Late_barplot.png", sep="") )
  }
  
  
  ###
  # Figure S4: Myeloid
  ###
  {
    mean(colData(dat.se)$myeloid[which(colData(dat.se)$state_rx=="c3"&colData(dat.se)$Group=="B")])
    mean(colData(dat.se)$myeloid[which(colData(dat.se)$state_rx=="c2"&colData(dat.se)$Group=="B")])
    
    png(paste(manFig_out,"/Fig-S3_myeloid.vs.CML_space_boxplot.png",sep=""), res=plot_res, units="in", height=3, width=5)
    p <- ggplot( col.df[which(col.df$Group=="C" | col.df$Group=="B"),], aes(x=state_rx, y=myeloid, fill=state_rx)) + geom_boxplot() +
      theme_bw(base_size=18) + scale_fill_manual(values=state_rx.cols)
    print(p)
    graphics.off()
    
    ### growth >:( ###
    png(paste(manFig_out,"/Fig-S3_spline.dmdt_group-",g,"-ONLY_color-state_rx_boxplot.png",sep=""),  res=plot_res, units="in", height=3, width=5)
    p <- ggplot(meta.df[which(meta.df$Group=="C" | meta.df$Group=="B"),], aes(x=state_rx, y=s.d_myeloid, fill=state_rx)) + 
      geom_boxplot() + theme_bw(base_size=16) +
      scale_fill_manual(values=state_rx.cols) 
    print(p)
    graphics.off()
    
    png(paste(manFig_out,"/Fig-S3_spline.dmdt_group-",g,"-ONLY_color-state_rx_boxplot.png",sep=""),  res=plot_res, units="in", height=3, width=5)
    p <- ggplot(meta.df[which(meta.df$Group=="C" | meta.df$Group=="B"),], aes(x=CML_space, y=s.d_myeloid, fill=state_rx, color=state_rx, group=state_rx)) + 
      geom_point() + geom_smooth(method="lm") + theme_bw(base_size=16) +
      scale_fill_manual(values=state_rx.cols) + scale_color_manual(values=state_rx.cols) +
      scale_x_continuous(limits=c(200,-300),trans = "reverse")
    print(p)
    graphics.off()
  }
  
  ###
  # Figure S4: CML contribution
  ###
  {
    gene_set_CML_contribution_manFigS4 <- function( out_name, deg_list, deg_comp, gene.df = gene.df ) { #all
      
      curdeg <- out_name #deg_name[[i]] #NEEDED:change curpath to curdeg
      cgenes.id <- deg_list
      curcomp <- deg_comp
      
      sel <- match_genes( gene.df$gene_id, cgenes.id)
      c.df <- gene.df[sel,]
      
      m.1 <- c.df[[paste(curcomp,"_mean.",1,sep="")]]
      m.2 <- c.df[[paste(curcomp,"_mean.",2,sep="")]]
      c.diff <- m.1 - m.2
      c.ld <- c.df$eigengenes
      up.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]
      dn.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]
      c.all <- c.diff * c.ld #get contribution for all genes in pathway
      c.up <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
      c.down <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
      ctab <- matrix( c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", sum(c.ld), sum(up.ld), sum(dn.ld), sum(c.all,na.rm=T), sum(c.up), sum(c.down) ), nrow=3, ncol=5) 
      #add lfc contribution; need to figure out why "Inf" is produced
      c.lfc <- log2((m.1+amin)/(m.2+amin) )
      cl.all <- c.ld * c.lfc #get contribution for all genes in pathway based on LFC
      cl.up <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
      cl.down <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
      ctab <- matrix(c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", 
                        length(which(!is.na(c.all))), length(which(!is.na(c.up))), length(which(!is.na(c.down))), # count of genes contributing (non NA genes)
                        sum(c.ld), sum(up.ld), sum(dn.ld), #eigenvalue
                        sum(c.all,na.rm=T), sum(c.up, na.rm=T), sum(c.down, na.rm=T), #contribution w/ exp diff
                        sum(cl.all,na.rm=T), sum(cl.up, na.rm=T), sum(cl.down, na.rm=T) ), # contrib. w/ log2FC
                     nrow=3, ncol=7)
      if (length(which(is.na(ctab))) != 0) {
        print(paste(":::ERROR::: NAs produced for comp: ",curcomp," pathway: ",curdeg,sep=""))
        print(ctab)
        stop("NA encountered")
      }
      cont.tab <- rbind(cont.tab, ctab)
      #plot space
      c.df[["exp.diff"]] <- c.diff
      c.df[["cont.all"]] <- c.all
      dif <- c.df$exp.diff
      # dif[which(c.df$exp.diff>0)] <- "up"
      # dif[which(c.df$exp.diff<0)] <- "down"
      dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]>0)] <- "up"
      dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]<0)] <- "down"
      c.df[["diff.lab"]] <- dif
      #build loading contribution
      "eigen.up" = mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      "eigen.down" <- mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
      "PC1.up" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      "PC1.down" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      
      ### make "eigen-contribution"
      # positive value = positive contribution (neg eig + pos exp OR pos eig + neg exp)
      # neg contribution (neg eig + neg exp OR pos eig + pos exp)
      cml.c <- c.df$eigengenes
      cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 #pro-CML up (neg eig. -> adjust)
      cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 #pro-CML down ( pos eig. -> NO adjust)
      cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 # anti-CML down (neg -> NO adjust)
      cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 # anti-CML up (pos -> adjust)
      
      #sums
      neg.up <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
      neg.up.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
      pos.down <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
      pos.down.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
      neg.down <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
      neg.down.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
      pos.up <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
      pos.up.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
      cmat <- matrix( c( rep(curdeg, 4), rep(curcomp, 4), "neg|up", "pos|down", "neg|down", "pos|up", 
                         neg.up, pos.down, neg.down, pos.up,
                         neg.up.l, pos.down.l, neg.down.l, pos.up.l),
                      nrow=4, ncol=5) 
      cml.tab <- rbind(cml.tab, cmat )
      #cml label
      cml.lab <- rep(NA, length(cml.c))
      cml.lab[which(cml.c > 0)] <- "pro-CML"
      cml.lab[which(cml.c < 0)] <- "anti-CML"
      #full matrix
      # :::note::: the mean contribution of each class of genes must be added!
      #(neg.up  / neg.up.l) + (pos.down / pos.down.l), (neg.down / neg.down.l) + (pos.up / pos.up.l)
      n.u.m <- (neg.up  / neg.up.l) 
      p.d.m <- (pos.down / pos.down.l)
      n.d.m <- (neg.down / neg.down.l) 
      p.u.m <- (pos.up / pos.up.l)
      fmat <- matrix( c( rep(curdeg, 2), rep(curcomp, 2), "pro-CML", "anti-CML",
                         sum(c(n.u.m, p.d.m), na.rm=T), sum(c(n.d.m, p.u.m), na.rm=T),
                         neg.up.l + pos.down.l, neg.down.l + pos.up.l),
                      nrow=2, ncol=5) 
      cml.full <- rbind(cml.full, fmat )
      ### eigengene value vs PC1 ###
      cl.up <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
      cl.down <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
      al.up <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
      al.down <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
      c.df[["CML_contribution"]] <- cml.c
      cml.up = mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      cml.down <- mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
      PC1.up <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      PC1.down <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      cml.tot <- mean(cml.c, na.rm=T) 
      PC1.tot <- mean(c.df$PC1, na.rm=T)
      #make df for labeling genes without NA genes
      l.df <- c.df[which(c.df$diff.lab=="up" | c.df$diff.lab=="down" ),]
      tup <- sort.int(l.df$CML_contribution, index.return = T, decreasing = T)$ix[seq(1,10)]
      tdown <- sort.int(l.df$CML_contribution, index.return = T, decreasing = F)$ix[seq(1,5)]
      png(paste(manFig_out,"/Fig-S4D_eigenDEGs_",curdeg,"_comparison-",curcomp,"_pts+label.png",sep=""), res=plot_res, units="in", height=4, width=4)
      p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab, fill=diff.lab )) + 
        geom_point(aes(shape=diff.lab), alpha=.75) + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", size=1) +
        # geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.2, 
        #              arrow = arrow(length = unit(.25, "cm"), type="open"), lineend="butt", linejoin="mitre") +
        # geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.2,
        #              arrow = arrow(length = unit(.25, "cm"), type="open"), lineend="butt", linejoin="mitre") +
        geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.25,
                     arrow = arrow(length = unit(.30, "cm"), type="closed"), lineend="butt", linejoin="mitre") +
        scale_color_manual(values= list("up"="#f76a95", "down" = "#4287f5" )) + scale_fill_manual(values= list("up"="#f76a95", "down" = "#4287f5" )) +
        # geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene), size=14) +
        # geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene), size=14) +
        theme(legend.position = "none") + scale_shape_manual(values=c("up"=24, "down"=25)) +
        scale_x_continuous(limits=c(-.03, 0.058)) + scale_y_continuous(limits=c(-.025, 0.018)) 
      #ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
      #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" )) 
      print(p)
      graphics.off()
      # png(paste(plot_out,"/eigenDEGs_",curdeg,"_comparison-",curcomp,"_pts+label+tot.png",sep=""), res=plot_res, units="in", height=8, width=10)
      # p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab )) + 
      #   geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", linewidth=1.25) +
      #   geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.5, 
      #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
      #   geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.5,
      #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
      #   geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.5,
      #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
      #   scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
      #   geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene)) +
      #   geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene)) +
      #   ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
      # #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" )) 
      # print(p)
      # graphics.off()
      
      
      
      
    } #end gene set contribution function
    
    deg_list <- list( early.deg.genes.id, trans.deg.genes.id, late.deg.genes.id, 
                      ig.b1.c, ig.b2.c, ig.b3.c,
                      ig.b3.b1, ig.b3.b2, ig.b2.b1,
                      ig.apost.c, ig.dpost.c, ig.drx.c, ig.c.evl, ig.c.lve, ig.apost.cL, ig.rx.da)
    deg_name <- list("Early",  "Transition", "Late", 
                     "B1.vs.C", "B2.vs.C", "B3.vs.C",
                     "B3.vs.B1", "B3.vs.B2", "B2.vs.B1",
                     "Apost.vs.C", "Dpost.vs.C", "Drx.vs.C", "C_early.vs.late", "C_late.vs.early", "Apost.vs.Clate", "Drx.vs.Apost")
    deg_comp <- list( "B1.vs.C", "B2.vs.C", "B3.vs.C", 
                      "B1.vs.C", "B2.vs.C", "B3.vs.C",
                      "B3.vs.B1", "B3.vs.B2", "B2.vs.B1",
                      "Apost.vs.C", "Dpost.vs.C", "Drx.vs.C","C_early.vs.late", "C_late.vs.early", "Apost.vs.Clate", "Drx.vs.Apost")
    for (i in c(1:3,7:9)){
      gene_set_CML_contribution_manFigS4( deg_name[i], deg_list[[i]], deg_comp[i], gene.df=gene.df)
    }
  }
  
  ###
  # Figure S4: Gene set CML contribution
  ###
  {
    gsea.b2.b1 <- read.table("CML_contribution_tables/CML.contribution_B-state-c3.vs.B-state-c1_GSEA_table.tsv", sep="\t", header=T)
    gsea.b3.b2 <- read.table("CML_contribution_tables/CML.contribution_B-state-c5.vs.B-state-c3_GSEA_table.tsv", sep="\t", header=T)
    gsea.b3.b1 <- read.table("CML_contribution_tables/CML.contribution_B-state-c5.vs.B-state-c1_GSEA_table.tsv", sep="\t", header=T)
    plimit <- 0.0001
    sel.b2.b1 <- which(gsea.b2.b1$padj < plimit)
    sel.b3.b2 <- which(gsea.b3.b2$padj < plimit)
    sel.b3.b1 <- which(gsea.b3.b1$padj < plimit)
    
    # color by state-space #
    {
      all.sig <-  sort(unique(c(gsea.b2.b1[sel.b2.b1,c("pathway")], gsea.b3.b2[sel.b3.b2,c("pathway")], gsea.b3.b1[sel.b3.b1,c("pathway")] )))
      df.list <- list(gsea.b2.b1[sel.b2.b1,], gsea.b3.b2[sel.b3.b2,], gsea.b3.b1[sel.b3.b1,] )
      df.name <- c("c2.vs.c1", "c3.vs.c2", "c3.vs.c1")
      sig.nes <- c()
      for (df in df.list) {
        out <- c()
        for (p in all.sig) {
          m <- which(df$pathway==p)
          if (length(m)==0) {
            out <- c(out,NA)
          } else {
            out <- c(out, df[["leadingEdge.CML_contribution"]][m])
          }
        }
        sig.nes <- cbind(sig.nes, out)
      }
      colnames(sig.nes) <- df.name
      rownames(sig.nes) <- all.sig
      sig.df <- data.frame(sig.nes)
      dim(sig.nes)
      inter.cml.cols <- c("c2.vs.c1"="#a9f5b2", "c3.vs.c2"="#facba5", "c3.vs.c1" = "#faa5c0")
      state.df <- data.frame("pathway" = rep(rownames(sig.df),3), "contribution" = c( sig.df$c2.vs.c1, sig.df$c3.vs.c2, sig.df$c3.vs.c1 ), 
                             "state" = c(rep("c2.vs.c1", dim(sig.df)[1]),rep("c3.vs.c2", dim(sig.df)[1]), rep("c3.vs.c1", dim(sig.df)[1])  ) )
      state.df$pathway <-  unlist(lapply(state.df$pathway, function(x) str_to_title(gsub("_"," ", gsub( "HALLMARK_", "", x)) )) )
      state.df$state <- factor(state.df$state, levels=c("c2.vs.c1", "c3.vs.c2", "c3.vs.c1"))
      # line_positions holds divider line between each variable
      state.df2 <- state.df %>% mutate(line_positions = as.numeric(factor(pathway, levels = unique(pathway))), 
                                       line_positions = line_positions + .5,  
                                       line_positions = ifelse(line_positions == max(line_positions), NA, line_positions)) 
      head(state.df2)
      png(paste(manFig_out,"/Fig-S4_fgsea-CML-contribution_p-",plimit,"_bar.png",sep=""), res=plot_res, units="in", height=4, width=12)
      p <- ggplot(state.df2, aes(x=pathway, y=contribution, fill=state)) + 
        geom_bar(stat = "identity", position=position_dodge(width = 0.5), width=.5) + 
        theme_bw(base_size=12) + scale_fill_manual(values=inter.cml.cols) +
        theme(axis.text.x = element_text(angle = 60,  hjust=1), legend.position="none", 
              panel.grid.major.x = element_blank(), panel.grid.minor.x = element_line(colour = "black") ) +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
        geom_vline(aes(xintercept = line_positions))
      print(p)
      graphics.off()
      png(paste(manFig_out,"/Fig-S4_fgsea-CML-contribution_p-",plimit,"_bar-rev.png",sep=""), res=plot_res, units="in", height=10, width=4)
      p <- ggplot(state.df2, aes(y=pathway, x=contribution, fill=state)) + 
        geom_bar(stat = "identity", position=position_dodge(width = 0.5), width=.5) + 
        theme_bw(base_size=12) + scale_fill_manual(values=inter.cml.cols) +
        theme( legend.position="none", 
              panel.grid.major.y = element_blank(), panel.grid.minor.y = element_line(colour = "black") ) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
        geom_hline(aes(yintercept = line_positions))
      print(p)
      graphics.off()
    }
  }
  

}




###
### DEG-based expression dynamics
###
{
  #
  # :note: analysis/preparation for Figure 4A
  #
  #CML vs Control DEGs
  ctl.degs <- list("c1" = b0.vs.c$gene_name[b0vc.degs], "c2" = b1.vs.c$gene_name[b1vc.degs], "c3" = b2.vs.c$gene_name[b2vc.degs] )
  ctl.degs.sel <- list("c1" = b0vc.degs, "c2" = b1vc.degs, "c3" = b2vc.degs )
  ctl.names <- c("c1", "c2", "c3" )
  #CML vs CML
  cml.degs <- list("early"=b1.vs.b0$gene_name[b1v0.degs], "late"=b2.vs.b1$gene_name[b2v1.degs], "full"=b2.vs.b0$gene_name[b2v0.degs] )
  cml.degs.sel <- list("early"=b1v0.degs, "late"=b2v1.degs, "full"=b2v0.degs )
  cml.names <- c("early", "late", "full")
  #make venns
  ctl.venn <- venn3_from_named_list(ctl.degs)
  cml.venn <- venn3_from_named_list(cml.degs)
  ctl.sel.venn <- venn3_from_named_list(ctl.degs.sel)
  cml.sel.venn <- venn3_from_named_list(cml.degs.sel)
  
  # c2 vs ctrl unique DEGs
  c2.uniq.genes <- ctl.venn[["sets"]][["a2"]]
  c2.uniq.sel <- match_genes(rowData(dat.se)$gene_name, c2.uniq.genes)  
  length(c2.uniq.genes)
  
  ###
  ### CORRELATION 
  ###
  # :::note::: this should be made into a function that passes: selected genes, output name
  # NO CTRL SAMPLES
  b.sel <- which(colData(dat.se)$Group=="B")
  gene.sets <- list( early.deg.genes.id, trans.deg.genes.id, late.deg.genes.id, 
                     ig.b1.c, ig.b2.c, ig.b3.c,
                     ig.b3.b2, ig.apost.c, ig.dpost.c, ig.drx.c)
  set.names <- c("Early-DEGs", "Transition-DEGs", "Late-DEGs", 
                 'c1-DEGs', "c2-DEGs", "c3-DEGs" ,
                 "c3.vs.c2-DEGs", "A-postRx-DEGs", "D-Rx-DEGs", "D-postRx-DEGs")
  deg_comp <- list( "B1.vs.C", "B2.vs.C", "B3.vs.C", 
                    "B1.vs.C", "B2.vs.C", "B3.vs.C",
                    "B3.vs.B2", "Apost.vs.C", "Dpost.vs.C", "Drx.vs.C")
  knums <- c(4, 2, 3, 
             4, 4, 4, 
             4, 4, 4, 4)
  #for (si in 1:length(gene.sets) ) {
  for (si in 2:3 ) {
    curgenes <- gene.sets[[si]] #holds gene (unique) IDs 
    set.name <- set.names[si]
    set.sel <- match_genes(rownames(dat.se), curgenes)
    knum <- knums[si]
    #samp.dat <- scale( t(log(data.matrix(assay(dat.se, "abundance"))+amin)), scale=F )[b.sel, set.sel]
    samp.dat <- t(data.matrix(assay(dat.se, "almc")))[b.sel, set.sel]
    cgenes <- rowData(dat.se)$gene_name[set.sel] #holds gene names
    g.info <- rowData(dat.se[set.sel])
    basemean <- rowMeans(log(data.matrix(assay(dat.se[ set.sel, b.sel], "abundance"))+amin))
    g.info[["mean.log.abd"]] <- basemean
    colnames(samp.dat) <- cgenes
    scor <- cor( samp.dat )
    inc <- 20
    blist <- seq(-1, 1, length.out=inc+1)
    corcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(inc)
    phem <- pheatmap::pheatmap(scor, scale="none", color=corcol, breaks=blist, show_rownames=F, show_colnames = F)
    tcutem <- cutree(phem$tree_row, k=knum) #or h=5
    tcutem <- as.character(tcutem)
    names(tcutem) <- cgenes
    tanndf <- data.frame( "Groups" = as.character(tcutem) )
    colnames(tanndf) <- "Group_Number"
    colnames(scor) <- cgenes
    trbcol <- grDevices::rainbow(length(unique(tcutem))) 
    tanncol <- trbcol
    names(tanncol) <- as.numeric(unique(tcutem))
    tanncol <- list(Group_Number = tanncol)
    rownames(scor) <- rownames(tanndf) #match rownames; use numeric to avoid duplicates
    graphics.off()
    png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_pheatmap.png",sep=""), res=plot_res, units="in", height=8, width=10) 
    print(pheatmap(scor, show_rownames=F, show_colnames=F, breaks=blist, color=corcol, 
                   annotation_row=tanndf, annotation_names_row = F, annotation_color=tanncol) )#, cellheight=wd, cellwidth=wd)
    graphics.off()
    for (gr in unique(tcutem) ) {
      csel <- which(tcutem==gr)
      cdat <- samp.dat[,csel]
      gr.info <- g.info[csel,]
      c.cor <- scor[csel, csel]
      usel <- upper.tri(c.cor, diag=F)
      cm.cor <- mean(c.cor[usel])
      med.cor <- median(c.cor[usel])
      png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_Grp-",gr,"_phatmap.png",sep=""), res=plot_res, units="in", height=5, width=6) 
      print(pheatmap(c.cor, show_rownames=F, show_colnames=F, breaks=blist, color=corcol, 
                     main=paste("Average correlation: ",round(cm.cor, digits=4),"\n",
                                "Median correlation: ",round(med.cor, digits=4),sep="")) ) 
      graphics.off()
      write.table(gr.info, paste(plot_out,"/DEG-dynamics_set-",set.name,"_Grp-",gr,"_geneInfo.tsv",sep=""),sep="\t", "row.names"=F )
      #mean expression
      cmean <- rowMeans(cdat)
      m.df <- data.frame(cmean)
      colnames(m.df) <- "meanExp"
      m.df[["sample"]] <- rownames(cdat)
      c.info <- data.frame(colData(dat.se)[b.sel,])
      cm.df <- merge(m.df, c.info, by="sample")
      ### NEEDED ###
      #  add in controls: only controls plot AND plot with both CML and controlsl; except they'll all be smoothed together on the x axis
      for (grp in c("A","C","D")) {
        cur.sel <- which(colData(dat.se)$Group==grp)
        alt.genes <- gr.info$gene_id
        alt.sel <- match_genes(rowData(dat.se)$gene_id, alt.genes) #get genes using "gene_id"
        #alt.dat <- scale( t(log(data.matrix(assay(dat.se, "abundance"))+amin)), scale=F )[cur.sel, alt.sel ]
        alt.dat <- t(data.matrix(assay(dat.se, "almc")))[cur.sel, alt.sel ]
        altmean <- rowMeans(alt.dat)
        a.df <- data.frame(altmean)
        colnames(a.df) <- "meanExp"
        a.df[["sample"]] <- rownames(alt.dat)
        alt.info <- data.frame(colData(dat.se)[cur.sel,])
        alt.df <- merge(a.df, alt.info, by="sample")
        #C data
        ctlsel <- which(colData(dat.se)$Group=="C")
        ctl.dat <- t(data.matrix(assay(dat.se, "almc")))[cur.sel, ctlsel ]
        ctlmean <- rowMeans(ctl.dat)
        ctl.df <- data.frame(ctlmean)
        colnames(ctl.df) <- "meanExp"
        ctl.df[["sample"]] <- rownames(ctl.dat)
        ctl.info <- data.frame(colData(dat.se)[ctlsel,])
        ctl.df <- merge(ctl.df, ctl.info, by="sample")
        
        png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_SampleGroup-",grp,"_DynGrp-",gr,"_CML.space_point+smooth_x-reverse.png",sep=""), res=plot_res, units="in", height=5, width=6) 
        p <- ggplot(alt.df, aes(x=CML_space, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
          theme_bw(base_size=16)  + scale_x_continuous(limits=c(200,-300),trans = "reverse") + geom_vline(xintercept=ss.cps[c(2,4)], linetype="dashed")
        print(p)
        graphics.off()
        colnames(alt.df)
        colnames(cm.df)
        merge.df <- rbind(cm.df, alt.df)
        merge2.df <- rbind(merge.df, ctl.df)
        if (grp=="A") {
          png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_SampleGroup-",grp,"+B_postRx_DynGrp-",gr,"_CML.space_point+smooth_x-reverse.png",sep=""), res=plot_res, units="in", height=5, width=6) 
          p <- ggplot(merge.df[which(merge.df$rx_class=="A.post"|merge.df$Group=="B"),], aes(x=CML_space, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
            theme_bw(base_size=16)  + scale_x_continuous(limits=c(200,-300),trans = "reverse") + geom_vline(xintercept=ss.cps[c(2,4)], linetype="dashed")+
            scale_color_manual(values=unlist(group.cols) )
          print(p)
          graphics.off()
          #+C
          png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_SampleGroup-",grp,"+B+C_postRx_DynGrp-",gr,"_CML.space_point+smooth_x-reverse.png",sep=""), res=plot_res, units="in", height=5, width=6) 
          p <- ggplot(merge2.df[which(merge.df$rx_class=="A.post"|merge.df$Group=="B"|merge.df$Group=="C"),], aes(x=CML_space, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
            theme_bw(base_size=16)  + scale_x_continuous(limits=c(200,-300),trans = "reverse") + geom_vline(xintercept=ss.cps[c(2,4)], linetype="dashed")+
            scale_color_manual(values=unlist(group.cols) )
          print(p)
          graphics.off()
        }
        if (grp=="D") {
          png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_SampleGroup-",grp,"+B_postRx_DynGrp-",gr,"_CML.space_point+smooth_x-reverse.png",sep=""), res=plot_res, units="in", height=5, width=6) 
          p <- ggplot(merge.df[which(merge.df$rx_class=="D.post"|merge.df$Group=="B"),], aes(x=CML_space, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
            theme_bw(base_size=16)  + scale_x_continuous(limits=c(200,-300),trans = "reverse") + geom_vline(xintercept=ss.cps[c(2,4)], linetype="dashed")+
            scale_color_manual(values=unlist(group.cols) )
          print(p)
          graphics.off()
          #+C
          png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_SampleGroup-",grp,"+B+C_postRx_DynGrp-",gr,"_CML.space_point+smooth_x-reverse.png",sep=""), res=plot_res, units="in", height=5, width=6) 
          p <- ggplot(merge2.df[which(merge.df$rx_class=="D.post"|merge.df$Group=="B"|merge.df$Group=="C"),], aes(x=CML_space, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
            theme_bw(base_size=16)  + scale_x_continuous(limits=c(200,-300),trans = "reverse") + geom_vline(xintercept=ss.cps[c(2,4)], linetype="dashed")+
            scale_color_manual(values=unlist(group.cols) )
          print(p)
          graphics.off()
          png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_SampleGroup-",grp,"+B_Rx_DynGrp-",gr,"_CML.space_point+smooth_x-reverse.png",sep=""), res=plot_res, units="in", height=5, width=6) 
          p <- ggplot(merge.df[which(merge.df$rx_class=="D.Rx"|merge.df$Group=="B"),], aes(x=CML_space, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
            theme_bw(base_size=16)  + scale_x_continuous(limits=c(200,-300),trans = "reverse") + geom_vline(xintercept=ss.cps[c(2,4)], linetype="dashed") +
            scale_color_manual(values=unlist(group.cols) )
          print(p)
          graphics.off()
          #+C
          png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_SampleGroup-",grp,"+B+C_Rx_DynGrp-",gr,"_CML.space_point+smooth_x-reverse.png",sep=""), res=plot_res, units="in", height=5, width=6) 
          p <- ggplot(merge2.df[which(merge.df$rx_class=="D.Rx"|merge.df$Group=="B"|merge.df$Group=="C"),], aes(x=CML_space, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
            theme_bw(base_size=16)  + scale_x_continuous(limits=c(200,-300),trans = "reverse") + geom_vline(xintercept=ss.cps[c(2,4)], linetype="dashed")+
            scale_color_manual(values=unlist(group.cols) )
          print(p)
          graphics.off()
        }
        if (grp=="C") {
          png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_SampleGroup-",grp,"_DynGrp-",gr,"_time_point+smooth.png",sep=""), res=plot_res, units="in", height=5, width=6) 
          p <- ggplot(alt.df, aes(x=sample_weeks, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
            theme_bw(base_size=16)  +  
            scale_color_manual(values=unlist(group.cols) )
          print(p)
          graphics.off()
          png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_SampleGroup-",grp,"+B_DynGrp-",gr,"_time_point+smooth.png",sep=""), res=plot_res, units="in", height=5, width=6) 
          p <- ggplot(merge.df[which(merge.df$Group=="C"|merge.df$Group=="B"),], aes(x=sample_weeks, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
            theme_bw(base_size=16)  + 
            scale_color_manual(values=unlist(group.cols) )
          print(p)
          graphics.off()
          png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_SampleGroup-",grp,"+B_DynGrp-",gr,"_CML.space_point+smooth_x-reverse.png",sep=""), res=plot_res, units="in", height=5, width=6) 
          p <- ggplot(merge.df[which(merge.df$Group=="C"|merge.df$Group=="B"),], aes(x=CML_space, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
            theme_bw(base_size=16)  + scale_x_continuous(limits=c(200,-300),trans = "reverse") + geom_vline(xintercept=ss.cps[c(2,4)], linetype="dashed")+
            scale_color_manual(values=unlist(group.cols) )
          print(p)
          graphics.off()
        }
      }
      png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_Grp-",gr,"_CML.space_point+smooth_x-reverse.png",sep=""), res=plot_res, units="in", height=5, width=6) 
      p <- ggplot(cm.df, aes(x=CML_space, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
        theme_bw(base_size=16) + scale_x_continuous(trans = "reverse") + geom_vline(xintercept=ss.cps[c(2,4)], linetype="dashed")
      print(p)
      graphics.off()
      
      png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_Grp-",gr,"_CML.space_point+smooth.png",sep=""), res=plot_res, units="in", height=5, width=6) 
      p <- ggplot(cm.df, aes(x=CML_space, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
        theme_bw(base_size=16) + ggtitle(paste("Mean correlation: ",round(cm.cor, digits=4),sep=""))
      print(p)
      graphics.off()
      
      if (cm.cor > .5) {
        png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_CORR-SEL_Grp-",gr,"_CML.space_point+smooth_x-reverse.png",sep=""), res=plot_res, units="in", height=5, width=6) 
        p <- ggplot(cm.df, aes(x=CML_space, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
          theme_bw(base_size=16) + scale_x_continuous(trans = "reverse") + geom_vline(xintercept=ss.cps[c(2,4)], linetype="dashed")
        print(p)
        graphics.off()
      }
      

      png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_Grp-",gr,"_myeloid_point+smooth.png",sep=""), res=plot_res, units="in", height=5, width=6) 
      p <- ggplot(cm.df, aes(x=myeloid, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
        theme_bw(base_size=16)
      print(p)
      graphics.off()
      #sample color
      png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_Grp-",gr,"_CML.space_point+smooth.png",sep=""), res=plot_res, units="in", height=5, width=6) 
      p <- ggplot(cm.df, aes(x=CML_space, y=meanExp)) + geom_smooth() +
        geom_point(data=cm.df, aes(x=CML_space, y=meanExp, color=mouse_id), alpha=.75) + theme_bw(base_size=16)
      print(p)
      graphics.off()
      png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_Grp-",gr,"_CML.space_point+path-ss.png",sep=""), res=plot_res, units="in", height=5, width=6) 
      p <- ggplot(cm.df, aes(x=CML_space, y=meanExp, group=mouse_id, color=mouse_id)) + geom_path(alpha=.75) +
        geom_point(data=cm.df, aes(x=CML_space, y=meanExp, color=mouse_id), alpha=.75) + theme_bw(base_size=16)
      print(p)
      graphics.off()
      png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_Grp-",gr,"_CML.space_point+path-time.png",sep=""), res=plot_res, units="in", height=5, width=6) 
      p <- ggplot(cm.df, aes(x=CML_space, y=meanExp, group=mouse_id, color=mouse_id)) + geom_path(alpha=.75) +
        geom_point(data=cm.df, aes(x=CML_space, y=meanExp, color=mouse_id), alpha=.75) + theme_bw(base_size=16)
      print(p)
      graphics.off()

      png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_Grp-",gr,"_myeloid_point+smooth.png",sep=""), res=plot_res, units="in", height=5, width=6) 
      p <- ggplot(cm.df, aes(x=myeloid, y=meanExp)) + geom_smooth() +
        geom_point(data=cm.df, aes(x=myeloid, y=meanExp, color=mouse_id), alpha=.75) + theme_bw(base_size=16)
      print(p)
      graphics.off()

    } #end loop over gene modules
  }
  
}



###
###  Correlate expression with inflection points of the sample density curve
###
{
  #
  # :note: analysis/preparation for Figure 4B-C
  #
  
  
  ### load all CML vs Ctrl/CML DEGs into a table; these will be correlated with the potential
  {
    deg.files <- list.files(path=deg_out, pattern=glob2rx("DESeq-sexFactor_B-state-*DEG-table-means.tsv"), full.names=T )
    deg.genes <- c()
    all.degs.tab <- c() #make table of all DEGs
    for (f in deg.files) {
      curtab <- read.table(f, header=T)
      curgenes <- curtab$gene_name
      deg.genes <- c(deg.genes, curgenes)
      #make table output
      fname <- strsplit(f, "/")[[1]][2]
      comp <- gsub("_DEG-table-means.tsv","",gsub("DESeq-sexFactor_","",fname))
      outtab <- cbind(rep(comp,dim(curtab)[1]), curtab)
      all.degs.tab <- rbind(all.degs.tab, outtab)
    }
    colnames(all.degs.tab) <- c("comp", colnames(curtab))
    all.degs <- data.frame(all.degs.tab)
    deg.genes <- sort(unique(deg.genes))
    
    
    ### get expression table for deg.genes
    deg.sel <- match_genes(rowData(dat.se)$gene_name, deg.genes )
    deg.se <- dat.se[deg.sel,] #DEG only SE object
    
    ### correlate each gene expression with state-space polynomial
    # Use the spline function to interpolate the curve
    interp_dens <- spline(dens$x, dens$y, n = length(dens$x) * 10)
    
    # Define a function that interpolates the curve
    dens_function <- function(x) {
      approx(interp_dens$x, interp_dens$y, xout = x)$y
    }
  
    # get location of each sample on polynomial
    samp.dens <- as.numeric(dens_function(colData(dat.se)$CML_space))
  }
  
  ### c2 and c4 INFLECTION CORRELATION ###
  {
    for (pt.name in c("c4", "c2")) {
      if (pt.name=="c2") {
        pt.range=c(ss.cps[3], ss.cps[5])
      } else if (pt.name=="c4") {
        pt.range=c(ss.cps[1], ss.cps[3])
      }
      
      #only use CML ("B") samples that fall within the current inflection point
      pt.sel <- which(colData(dat.se)$Group=="B" & (colData(dat.se)$CML_space > pt.range[1] & colData(dat.se)$CML_space <= pt.range[2]) )
      
      ### perform correlation between gene expression and the local shape of potential ("pt.name"; either c2 or c4)
      deg.cor <- c()
      deg.lm <- c()
      for (i in 1:length(deg.sel)) {
        sel <- deg.sel[i]
        gname <- rowData(deg.se)$gene_name[i]
        l2fc <- rowData(deg.se)$log[i]
        cdat <- as.numeric(data.matrix(assay(deg.se[i,],"almc" )))
        cc <- cor(samp.dens[pt.sel], cdat[pt.sel]) #correlate only CML sample in current inflection point
        fit <- lm(cdat[pt.sel] ~ samp.dens[pt.sel])
        rs <- summary(fit)$r.squared
        ar <- summary(fit)$adj.r.squared
        pv <- summary(fit)$coefficients[2,4]
        
        deg.cor <- c(deg.cor, cc)
        deg.lm <- rbind(deg.lm, c(rs, ar, pv))
      }
      #make outdir
      cordegs <- paste("plots/DEG_potential-expression_correlation_",pt.name,sep="")
      dir.create(cordegs, showWarnings = F)
      
      cor.tab <- c()
      cor.sum <- c()
      #output gene plots to separate directory
      cordegs_genes <- paste("plots/DEG_potential-expression_correlation_",pt.name,"/gene_plots",sep="")
      dir.create(cordegs_genes, showWarnings = F)
      for (i in 1:length(deg.cor)) {
        cc <- deg.cor[i]
        pv <- deg.lm[i,3]
        if (is.na(cc)) {next } #need to figure out why this happens...
        if (abs(cc)<.5 | pv > 0.05) { next } 
        if (cc>0) {cdir <- "positive"}
        if (cc <0) {cdir <- "negative"}
        cgene <- rowData(deg.se)$gene_name[i]
        bsel <- which(colData(deg.se)$Group=="B")
        allexp <- as.numeric(data.matrix(assay(deg.se[i,],"almc" )))
        cexp <- allexp[bsel] #using only CML samples
        cdat <- cbind(rep(cgene, length(bsel)), cexp, 
                      colData(deg.se)$CML_space[bsel], rep(cdir, length(bsel)), rep(cc, length(bsel)) )
        cor.tab <- rbind(cor.tab, cdat)
        cmean <- mean(as.numeric(data.matrix(assay(deg.se[i,bsel],"abundance" ))))
        cor.sum <- rbind(cor.sum, c(cgene, cmean, "global", cdir, cc, deg.lm[i,] ))
        #plot genes; TURNED OFF TO SAVE TIME!!!
        model.pts <- rep("no", dim(deg.se)[2] )
        model.pts[pt.sel] <- "yes"
        gdf <- data.frame("expression"=cexp, "model" = model.pts[bsel], density = samp.dens[bsel],
                          "CML_space"=colData(deg.se)$CML_space[bsel], "mouse_id" = colData(deg.se)$mouse_id[bsel], state_rx = colData(deg.se)$state_rx[bsel] )
      }
      # output correlation results
      chead <- c("gene", "expression", "CML_space", "direction", "corrCoef", "lm-Rsquared", "lm-adjRsquared","lm-pvalue" )
      write.table(rbind(chead,cor.sum),paste(cordegs,"/corrTable_top-genes-",pt.name,".tsv",sep=""), sep="\t", col.names=F, row.names=F )
    }
  }
  
  ### protein-protein interaction (PPI) analysis for c2 & c4 inflection pts ###
  {
    
    #load genes
    iec.c2.tab <- read.table(paste("plots/DEG_potential-expression_correlation_c2/corrTable_top-genes-c2.tsv",sep=""), sep="\t", header=T)
    iec.c4.tab <- read.table(paste("plots/DEG_potential-expression_correlation_c4/corrTable_top-genes-c4.tsv",sep=""), sep="\t", header=T)
    
    #human: 9606; mouse: 10090
    string_db <- STRINGdb$new( species = 10090, score_threshold = 0.4, version="11.0b")
    
    #setup biomaRt
    mm.ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl") #"hsapiens_gene_ensembl" For Homo sapiens
    attributes <- c("ensembl_gene_id", "external_gene_name")
    
    ### process each inflection pt ###
    #setup variables    
    share.tab <- data.frame("gene"=intersect(iec.c2.tab$gene,iec.c4.tab$gene))
    pt.tabs <- list("c2" = iec.c2.tab, "c4" = iec.c4.tab) # "c2+c4"=share.tab
    trans.pt <- list() #list to hold genes
    for (pt.name in c( "c2",  "c4") ) {
      print(paste("processing ",pt.name, sep=""))
      cordegs <- paste("plots/DEG_potential-expression_correlation_",pt.name,sep="")
      dir.create(cordegs, showWarnings = F)
      iec.orig <- string_db$map( pt.tabs[[pt.name]], "gene", removeUnmappedRows = TRUE )
      png(paste(cordegs,"/allDEGs_stringDB_",pt.name,"_network-allGenes.png",sep=""), res=plot_res, units="in", height=10, width=10)
      string_db$plot_network(iec.orig$STRING_id)
      graphics.off()
      png(paste(cordegs,"/allDEGs_stringDB_",pt.name,"_network-min900.png",sep=""), res=plot_res, units="in", height=10, width=10)
      string_db$plot_network(iec.orig$STRING_id, required_score=900)
      graphics.off()
      
      
      ### add log2FC to table ### 
      iec <- add_l2fc(iec.orig, all.degs)
      #extended
      int <- string_db$get_interactions( iec$STRING_id )
      png(paste(cordegs,"/allDEGs_stringDB_",pt.name,"_EXTENDED_network.png",sep=""), res=plot_res, units="in", height=10, width=10)
      string_db$plot_network(int[which(int$combined_score>900),])
      graphics.off()
      
      
      #add gene expression colors
      # log2FC coloring on network
      iec.lfc <- add_l2fc_wColor(iec.orig, all.degs)
      int.lfc <- string_db$get_interactions( iec.lfc$STRING_id ) #should be same as above, but running interactions here to be safe
      lfc_id <- string_db$post_payload( iec.lfc$STRING_id, colors=iec.lfc$color )
      png(paste(cordegs,"/allDEGs_stringDB_",pt.name,"_halo-l2FC_network-allGenes.png",sep=""), res=plot_res, units="in", height=10, width=10)
      string_db$plot_network( iec.orig$STRING_id, payload_id=lfc_id )
      graphics.off()
      png(paste(cordegs,"/allDEGs_stringDB_",pt.name,"_halo-l2FC_network-min-900.png",sep=""), res=plot_res, units="in", height=10, width=10)
      string_db$plot_network( iec.orig$STRING_id, payload_id=lfc_id, required_score=900 )
      graphics.off()
      #extended
      png(paste(cordegs,"/allDEGs_stringDB_",pt.name,"_halo-l2FC_EXTENDED_network.png",sep=""), res=plot_res, units="in", height=10, width=10)
      string_db$plot_network( int.lfc[which(int.lfc$combined_score>900),], payload_id=lfc_id )
      graphics.off()
      
      
      
      ##enrichment 
      enrich <- string_db$get_enrichment( iec.orig )
      write.table(enrich, paste(cordegs,"/allDEGs_stringDB-network_enrichment_table.tsv",sep=""),sep="\t", row.names = F )
      enrich.kegg <- string_db$get_enrichment( iec.orig, category="KEGG" )
      write.table(enrich.kegg, paste(cordegs,"/allDEGs_stringDB-network_KEGG_table.tsv",sep=""),sep="\t", row.names = F )
      #extended
      enrich <- string_db$get_enrichment( int.lfc[which(int.lfc$combined_score>900),] )
      write.table(enrich, paste(cordegs,"/allDEGs_stringDB-network_EXTENDED_enrichment_table.tsv",sep=""),sep="\t", row.names = F )
      enrich.kegg <- string_db$get_enrichment( int.lfc[which(int.lfc$combined_score>900),], category="KEGG" )
      write.table(enrich.kegg, paste(cordegs,"/allDEGs_stringDB-network_EXTENDED_KEGG_table.tsv",sep=""),sep="\t", row.names = F )
      
      ##mapping of StringDB IDs to gene names
      stringdb_ids <- unique(unlist(iec.orig$STRING_id)) 
      stringdb_ids <- unlist(lapply(stringdb_ids, function(x) gsub("10090.","",x)))
      gene_mapping <- getBM(attributes = attributes, filters = "ensembl_peptide_id", values = stringdb_ids, mart = mm.ensembl)
      genes <- gene_mapping[,2]
      trans.pt[[pt.name]] <- list("full" = genes)
      #extended
      stringdb_ids <- unique(unlist(int.lfc[which(int.lfc$combined_score>900),c(1,2)])) 
      stringdb_ids <- unlist(lapply(stringdb_ids, function(x) gsub("10090.","",x)))
      gene_mapping <- getBM(attributes = attributes, filters = "ensembl_peptide_id", values = stringdb_ids, mart = mm.ensembl)
      genes_ext <- gene_mapping[,2]
      trans.pt[[pt.name]] <- list("full" = genes, "extended" = genes_ext)
      
      ### enrichR ###
      cur_path_out <- paste(cordegs,"/enrichr_",pt.name,sep="")
      dir.create(cur_path_out, showWarnings = F)
      enriched.up <- enrichr(genes, test_dbs)
      output_enrichr_results(enriched.up, paste(pt.name,"-all",sep=""), cur_path_out)
      #extended
      cur_path_out <- paste(cordegs,"/enrichr_EXTENDED_",pt.name,sep="")
      dir.create(cur_path_out, showWarnings = F)
      enriched.up <- enrichr(genes_ext, test_dbs)
      output_enrichr_results(enriched.up, paste(pt.name,"-all",sep=""), cur_path_out)
      
      
      ### CML contribution ###
      #get gene IDs
      gsel <- match_all_genes( rowData(dat.se)$gene_name, genes)
      gene_ids <- rowData(dat.se)$gene_id[gsel]
      if (pt.name == "c2") {
        gene_set_CML_contribution(paste(pt.name,".trans.pt",sep=""), gene_ids, "B2.vs.B1", gene.df = gene.df, enrichOut = F, all.gene.labs=F)
      }
      if (pt.name == "c4") {
        gene_set_CML_contribution(paste(pt.name,".trans.pt",sep=""), gene_ids, "B3.vs.B2", gene.df = gene.df, enrichOut = F, all.gene.labs=F)
      }
      gene_set_CML_contribution(paste(pt.name,".trans.pt",sep=""), gene_ids, "B3.vs.C", gene.df = gene.df, enrichOut = F, all.gene.labs=F )
      #extended
      gsel <- match_all_genes( rowData(dat.se)$gene_name, genes_ext)
      gene_ids <- rowData(dat.se)$gene_id[gsel]
      if (pt.name == "c2") {
        gene_set_CML_contribution(paste("EXTENDED_",pt.name,".trans.pt",sep=""), gene_ids, "B2.vs.B1", gene.df = gene.df, enrichOut = F, all.gene.labs=F)
      }
      if (pt.name == "c4") {
        gene_set_CML_contribution(paste("EXTENDED_",pt.name,".trans.pt",sep=""), gene_ids, "B3.vs.B2", gene.df = gene.df, enrichOut = F, all.gene.labs=F)
      }
      gene_set_CML_contribution(paste("EXTENDED_",pt.name,".trans.pt",sep=""), gene_ids, "B3.vs.C", gene.df = gene.df, enrichOut = F, all.gene.labs=F )
      
      ### clustering ###
      {  
        clustersList <- string_db$get_clusters(iec.orig)
        cl.len <- min( c(length(clustersList), 4) )
        ## pane of clusters ##
        png(paste(cordegs,"/allDEGs_stringDB_",pt.name,"_clusters_halo-l2FC_network-lg-allGenes.png",sep=""), res=plot_res, units="in", height=15, width=15)
        par(mfrow=c(2,2))
        for(i in seq(1:cl.len)){
          string_db$plot_network(clustersList[[i]], payload_id=lfc_id)
        }
        graphics.off()
        
        png(paste(cordegs,"/allDEGs_stringDB_",pt.name,"_clusters_halo-l2FC_network-lg-min-900.png",sep=""), res=plot_res, units="in", height=15, width=15)
        par(mfrow=c(2,2))
        for(i in seq(1:cl.len)){
          string_db$plot_network(clustersList[[i]], payload_id=lfc_id, required_score = 900)
        }
        graphics.off()
        
        png(paste(cordegs,"/allDEGs_stringDB_",pt.name,"_clusters_halo-l2FC_network-lg-min-700.png",sep=""), res=plot_res, units="in", height=15, width=15)
        par(mfrow=c(2,2))
        for(i in seq(1:cl.len)){
          string_db$plot_network(clustersList[[i]], payload_id=lfc_id, required_score = 700)
        }
        graphics.off()
        
        ## enrichment + CML contribution##
        for(i in seq(1:cl.len)) {
          c.enr <- string_db$get_enrichment( clustersList[[i]] )
          write.table(c.enr, paste(cordegs,"/allDEGs_stringDB-network_cluster-",i,"_enrichment_table.tsv",sep=""),sep="\t", row.names = F )
          
          # enrichR #
          stringdb_ids <- unique(unlist( clustersList[[i]] )) 
          stringdb_ids <- unlist(lapply(stringdb_ids, function(x) gsub("10090.","",x)))
          gene_mapping <- getBM(attributes = attributes, filters = "ensembl_peptide_id", values = stringdb_ids, mart = mm.ensembl)
          
          genes <- gene_mapping[,2]
          trans.pt[[pt.name]][[i]] = genes
          
          cur_path_out <- paste(cordegs,"/enrichr_",pt.name,"_cluster-",i,sep="")
          dir.create(cur_path_out, showWarnings = F)
          enriched.up <- enrichr(genes, test_dbs)
          output_enrichr_results(enriched.up, paste(pt.name,"_c-",i,sep=""), cur_path_out)
           
          ### CML cont ###
          #get gene IDs
          gsel <- match_all_genes(rowData(dat.se)$gene_name, genes)
          gene_ids <- rowData(dat.se)$gene_id[gsel]
          if (pt.name == "c2") {
            gene_set_CML_contribution(paste(pt.name,".trans.pt_cluster-",i,sep=""), gene_ids, "B2.vs.B1", gene.df = gene.df, enrichOut = F, all.gene.labs=F)
          }
          if (pt.name == "c4") {
            gene_set_CML_contribution(paste(pt.name,".trans.pt_cluster-",i,sep=""), gene_ids, "B3.vs.B2", gene.df = gene.df, enrichOut = F, all.gene.labs=F)
          }
          gene_set_CML_contribution(paste(pt.name,".trans.pt_cluster-",i,sep=""), gene_ids, "B3.vs.C", gene.df = gene.df, enrichOut = F, all.gene.labs=F)
          
        }
        
        #individual cluster plots
        for(i in seq(1:cl.len)){
          png(paste(cordegs,"/allDEGs_stringDB_",pt.name,"_cluster-",i,"-ONLY_halo-l2FC_network-lg-allGenes.png",sep=""), res=plot_res, units="in", height=10, width=10)
          string_db$plot_network(clustersList[[i]], payload_id=lfc_id)
          graphics.off()
          png(paste(cordegs,"/allDEGs_stringDB_",pt.name,"_cluster-",i,"-ONLY_halo-l2FC_network-lg-min-900.png",sep=""), res=plot_res, units="in", height=10, width=10)
          string_db$plot_network(clustersList[[i]], payload_id=lfc_id, required_score = 900)
          graphics.off()
        }
      } #end clustering 
      

    }#end inflection point for loop
    
    
    ###upset and venns for critical point DEGs vs inflection point dynamic genes ###
    dvt <- list("c4_Trans" = pt.tabs[["c4"]]$gene, "c2_Trans"=pt.tabs[["c2"]]$gene, "Early_DEGs" = early.deg.genes, 
                "Trans_DEGs" = trans.deg.genes, "Late_DEGs" = late.deg.genes )
    cm.dvt = ComplexHeatmap::make_comb_mat(dvt)
    comb.sort <- sort.int(comb_size(cm.dvt), index.return = T)$ix
    png(paste(plot_out,"/venn5_transitionPt.vs.uniq-DEGs_upset.png",sep=""), res=plot_res, units="in", height=6, width=12)
    UpSet(cm.dvt, set_order = names(dvt), comb_order = rev(comb.sort),
          right_annotation = upset_right_annotation(cm.dvt, 
                                                    gp = gpar(fill = c("#2bcf64", "#eb8a52", "dodgerblue",  "goldenrod", "red" )) ),
           )
    graphics.off()

    dvt2 <- list("c4_Trans" = pt.tabs[["c4"]]$gene, "c2_Trans"=pt.tabs[["c2"]]$gene, "c1_DEGs" = dg.b1.c, 
                "c3_DEGs" = dg.b2.c, "c5_DEGs" = dg.b3.c, "c5.c1_DEGs" = dg.b3.b1, 
                "c5.c3_DEGs" = dg.b3.b2, "c3.c1_DEGs" = dg.b3.b1 )
    cm.dvt2 = ComplexHeatmap::make_comb_mat(dvt2)
    comb.sort <- sort.int(comb_size(cm.dvt2), index.return = T)$ix
    png(paste(plot_out,"/venn8_transitionPt.vs.ctl+inter-DEGs_upset.png",sep=""), res=plot_res, units="in", height=6, width=12)
    UpSet(cm.dvt2, set_order = names(dvt2), comb_order = rev(comb.sort),
          right_annotation = upset_right_annotation(cm.dvt2, 
                                                    gp = gpar(fill = c("#2bcf64", "#eb8a52", "dodgerblue",  "goldenrod", "red", "#faa5c0", "#facba5", "#a9f5b2" )) ),
    )
    graphics.off()
    cm.dvt2.sel <- cm.dvt2[comb_size(cm.dvt2)>20]
    comb.sort <- sort.int(comb_size(cm.dvt2.sel), index.return = T)$ix
    png(paste(plot_out,"/venn8_transitionPt.vs.ctl+inter-DEGs_minSet-20_upset.png",sep=""), res=plot_res, units="in", height=6, width=8)
    UpSet(cm.dvt2.sel, set_order = names(dvt2), comb_order = rev(comb.sort),
          right_annotation = upset_right_annotation(cm.dvt2.sel, 
                                                    gp = gpar(fill = c("#2bcf64", "#eb8a52", "dodgerblue",  "goldenrod", "red", "#faa5c0", "#facba5", "#a9f5b2" )) ),
    )
    graphics.off()
    
  }
  
  ### CML contribution tables: c2, c4 Inflection Pt Dynamics ###
  # :::note::: output used to construct Table S3
  {
    #get files
    enrich.files <- c("plots/DEG_potential-expression_correlation_c2/enrichr_EXTENDED_c2/MSigDB_Hallmark_2020_c2-all_table.tsv", 
                      "plots/DEG_potential-expression_correlation_c4/enrichr_EXTENDED_c4/MSigDB_Hallmark_2020_c4-all_table.tsv")
    enrich.comps <- c("B2.vs.C", "B3.vs.C")
    enrich.names <- c("c2.transPt", "c4.transPt")
    # process each file
    for (fi in 1:length(enrich.files) ) {
      f <- enrich.files[fi]
      comp <- enrich.comps[fi]
      compout <- enrich.names[fi]
      curcomp <- gsub("-Rx", "rx", gsub("-postRx", "post", gsub("-state-c","",comp) ) ) #make compatible with "gene.df" names
      grp1 <- strsplit(comp,".vs.")[[1]][1]
      grp2 <- strsplit(comp,".vs.")[[1]][2]
      glab1 <- gsub("c2", "c3", gsub("c3", "c5", grp1))
      glab2 <- gsub("c2", "c3", gsub("c3", "c5", grp2))
      cur.tab <- read.table(f, header=T)
      if (dim(cur.tab)[1]==0) {next} #skip empty
      out.tab <- c()
      for (pi in 1:dim(cur.tab)[1]) {#get CML contribution for each pathway
        curpath <- cur.tab$Term[pi]
        if (cur.tab$Adjusted.P.value[pi] > .2) {next} #skip high pvalues
        #leading edge
        le.genes <- strsplit(cur.tab$Genes[pi],";")[[1]]
        #modify genes to all be upper case to fix matching issue between human and mouse
        df.genes <- toupper(gene.df$gene)
        sel <- match_all_genes( df.genes, le.genes)
        le.df <- gene.df[sel,]
        cml.le <- le.df$eigengenes
        cml.le[which(le.df$eigengenes < 0 & le.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.le[which(le.df$eigengenes < 0 & le.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 #pro-CML up (neg eig. -> adjust)
        cml.le[which(le.df$eigengenes > 0 & le.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.le[which(le.df$eigengenes > 0 & le.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 #pro-CML down ( pos eig. -> NO adjust)
        cml.le[which(le.df$eigengenes < 0 & le.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.le[which(le.df$eigengenes < 0 & le.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 # anti-CML down (neg -> NO adjust)
        cml.le[which(le.df$eigengenes > 0 & le.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.le[which(le.df$eigengenes > 0 & le.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 # anti-CML up (pos -> adjust)
        le.pos.g <- paste(le.df$gene[which(cml.le>0)],collapse="|")
        le.neg.g <- paste(le.df$gene[which(cml.le<0)],collapse="|")
        le.pos.cnt <- length(cml.le[which(cml.le>0)])
        le.neg.cnt <- length(cml.le[which(cml.le<0)])
        le.pos.mean <- mean(cml.le[which(cml.le>0)])
        le.neg.mean <- mean(cml.le[which(cml.le<0)])
        le.tot <- mean(cml.le)
        
        
        #append current pathway data
        outline <- c(cur.tab[pi,], le.tot, le.pos.mean, le.neg.mean, le.pos.cnt, le.neg.cnt, le.pos.g, le.neg.g )
        out.tab <- rbind(out.tab, outline)             
      } #end pathway for loop
      
      ### output gsea table ###
      if (is.null(dim(out.tab)) ) {
        print(paste("Skipping: ",compout," with no significant pathways",sep=""))
        next
      }
      header <- c(colnames(cur.tab), "CML_contribution", "pro-CML_cont", "anti-CML_cont", 
                  "pro-CML_gene_count", "anti-CML_gene_count", "pro-CML_genes", "anti-CML_genes" ) 
      colnames(out.tab) <- header
      ### Table S3 ###
      write.table(out.tab, paste(cctab_out,"/transPt.CML.contribution_",compout,"_HallmarkEnrichR_table.tsv",sep=""),sep="\t", row.names=F)
    }# end file loop
    
  }
  
} #end inflection point expression dynamics section



###
### Figure 4: expression dynamics (Fig. S5-6)
###
{
  ###
  # Figure 4A
  ###
  {
    #CML vs Control DEGs
    ctl.degs <- list("c1" = b0.vs.c$gene_name[b0vc.degs], "c2" = b1.vs.c$gene_name[b1vc.degs], "c3" = b2.vs.c$gene_name[b2vc.degs] )
    ctl.degs.sel <- list("c1" = b0vc.degs, "c2" = b1vc.degs, "c3" = b2vc.degs )
    ctl.names <- c("c1", "c2", "c3" )
    #CML vs CML
    cml.degs <- list("early"=b1.vs.b0$gene_name[b1v0.degs], "late"=b2.vs.b1$gene_name[b2v1.degs], "full"=b2.vs.b0$gene_name[b2v0.degs] )
    cml.degs.sel <- list("early"=b1v0.degs, "late"=b2v1.degs, "full"=b2v0.degs )
    cml.names <- c("early", "late", "full")
    #make venns
    ctl.venn <- venn3_from_named_list(ctl.degs)
    cml.venn <- venn3_from_named_list(cml.degs)
    ctl.sel.venn <- venn3_from_named_list(ctl.degs.sel)
    cml.sel.venn <- venn3_from_named_list(cml.degs.sel)
    
    # c2 vs ctrl unique DEGs
    c2.uniq.genes <- ctl.venn[["sets"]][["a2"]]
    c2.uniq.sel <- match_genes(rowData(dat.se)$gene_name, c2.uniq.genes)  
    length(c2.uniq.genes)
    
    ###
    ### EXPRESSION CORRELATION ANALYSIS
    ###
    {
      # Control samples are not includes
      b.sel <- which(colData(dat.se)$Group=="B")
      gene.sets <- list( early.deg.genes.id, trans.deg.genes.id, late.deg.genes.id )
      set.names <- c("Early-DEGs", "Transition-DEGs", "Late-DEGs")
      # deg_comp <- list( "B1.vs.C", "B2.vs.C", "B3.vs.C") #DEG comparisons used to get 
      deg_comp <- NULL
      knums <- c(4, 2, 3 ) #number of groups cut using "cutree()"; decided after exploring different values
      for (si in 1:length(gene.sets) ) {
        curgenes <- gene.sets[[si]] #holds gene (unique) IDs 
        set.name <- set.names[si]
        set.sel <- match_genes(rownames(dat.se), curgenes)
        knum <- knums[si]
        #samp.dat <- scale( t(log(data.matrix(assay(dat.se, "abundance"))+amin)), scale=F )[b.sel, set.sel]
        samp.dat <- t(data.matrix(assay(dat.se, "almc")))[b.sel, set.sel]
        cgenes <- rowData(dat.se)$gene_name[set.sel] #holds gene names
        g.info <- rowData(dat.se[set.sel])
        basemean <- rowMeans(log(data.matrix(assay(dat.se[ set.sel, b.sel], "abundance"))+amin))
        g.info[["mean.log.abd"]] <- basemean
        colnames(samp.dat) <- cgenes
        scor <- cor( samp.dat )
        inc <- 20
        blist <- seq(-1, 1, length.out=inc+1)
        corcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(inc)
        phem <- pheatmap::pheatmap(scor, scale="none", color=corcol, breaks=blist, show_rownames=F, show_colnames = F)
        tcutem <- cutree(phem$tree_row, k=knum) #or h=5
        tcutem <- as.character(tcutem)
        names(tcutem) <- cgenes
        tanndf <- data.frame( "Groups" = as.character(tcutem) )
        colnames(tanndf) <- "Group_Number"
        colnames(scor) <- cgenes
        trbcol <- grDevices::rainbow(length(unique(tcutem))) 
        tanncol <- trbcol
        names(tanncol) <- as.numeric(unique(tcutem))
        tanncol <- list(Group_Number = tanncol)
        rownames(scor) <- rownames(tanndf) #match rownames; use numeric to avoid duplicates
        graphics.off()
        png(paste(manFig_out,"/Fig-S3_DEG-dynamics_set-",set.name,"_pheatmap.png",sep=""), res=plot_res, units="in", height=6, width=8) 
        print(pheatmap(scor, show_rownames=F, show_colnames=F, breaks=blist, color=corcol, 
                       annotation_row=tanndf, annotation_names_row = F, annotation_color=tanncol) )#, cellheight=wd, cellwidth=wd)
        graphics.off()
        for (gr in unique(tcutem) ) {
          csel <- which(tcutem==gr)
          cdat <- samp.dat[,csel]
          gr.info <- g.info[csel,]
          c.cor <- scor[csel, csel]
          usel <- upper.tri(c.cor, diag=F)
          cm.cor <- mean(c.cor[usel])
          #mean expression
          cmean <- rowMeans(cdat)
          m.df <- data.frame(cmean)
          colnames(m.df) <- "meanExp"
          m.df[["sample"]] <- rownames(cdat)
          c.info <- data.frame(colData(dat.se)[b.sel,])
          cm.df <- merge(m.df, c.info, by="sample")
          if (si==1 ) { #print blanks
            cm.df$meanExp <- rep(NA, dim(cm.df)[1])
            gr <- 0
          }
          png(paste(manFig_out,"/Fig-3B_DEG-dynamics_set-",set.name,"_Grp-",gr,"_CML.space_point+smooth_x-reverse.png",sep=""), res=plot_res, units="in", height=4, width=6) 
          p <- ggplot(cm.df, aes(x=CML_space, y=meanExp, color=Group)) + geom_point(alpha=.75) + geom_smooth() +
            theme_bw(base_size=18) + scale_x_continuous(trans = "reverse") + geom_vline(xintercept=ss.cps[c(2,4)], linetype="dashed") +
            geom_vline(xintercept=ss.cps[c(1,3,5)], linetype="dashed", color="grey") + theme(legend.position="none")
          print(p)
          graphics.off()
          
        } #end loop over gene modules
      }
    }
  
  }
  
  ###
  # Figure 4B-C: potential well
  ###
  {
    ### CML ###
    # Use the spline function to interpolate the curve
    interp_dens <- spline(dens$x, dens$y, n = length(dens$x) * 10)
    # Define a function that interpolates the curve
    dens_function <- function(x) {
      approx(interp_dens$x, interp_dens$y, xout = x)$y
    }
    #use critical points
    png(paste(manFig_out,"/Fig-3C-D_potential_B-ONLY_reverse-manusript.png",sep=""), res=plot_res, units="in", height=4, width=5) 
    plot(-1*dens$x, -1*dens$y, pch=19, col="black", type="l", ylab="", yaxt="n", xaxt="n", xlab="", xlim=c(-220, 325), bty="n" )
    points(-1*colData(dat.se)$CML_space[b.sel], -1*dens_function(colData(dat.se)$CML_space[b.sel])+0.00001, col="red", pch=19, cex=1)
    graphics.off()
  }
  
  ###
  # Figure 4C-D: stringDB network
  ###
  {
    
    #load genes found to be locally correlated with c2 and c4
    iec.c2.tab <- read.table(paste("plots/DEG_potential-expression_correlation_c2/corrTable_top-genes-c2.tsv",sep=""), sep="\t", header=T)
    iec.c4.tab <- read.table(paste("plots/DEG_potential-expression_correlation_c4/corrTable_top-genes-c4.tsv",sep=""), sep="\t", header=T)
  
    #human: 9606; mouse: 10090
    string_db <- STRINGdb$new( species = 10090, score_threshold = 0.4, version="11.0b")
    
    ### process both inflection points ###
    #setup biomaRt
    mm.ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl") #"hsapiens_gene_ensembl" For Homo sapiens
    attributes <- c("ensembl_gene_id", "external_gene_name")
    #setup variables    
    share.tab <- data.frame("gene"=intersect(iec.c2.tab$gene,iec.c4.tab$gene))
    pt.tabs <- list("c2" = iec.c2.tab, "c4" = iec.c4.tab, "c2+c4"=share.tab)
    
    ### PPI analysis for each inflection point ###
    for (pt.name in c( "c2",  "c4") ) {
      #for (pt.name in c("c4") ) {
      print(paste("processing ",pt.name, sep=""))
      #pt.name <- "c2" #testing
      cordegs <- paste("plots/DEG_potential-expression_correlation_",pt.name,sep="")
      dir.create(cordegs, showWarnings = F)
      iec.orig <- string_db$map( pt.tabs[[pt.name]], "gene", removeUnmappedRows = TRUE )
      
      
      
      ### add log2FC to table ### 
      iec <- add_l2fc(iec.orig, all.degs)
      #extended
      int <- string_db$get_interactions( iec$STRING_id )
      
      #add gene expression colors
      # log2FC coloring on network
      iec.lfc <- add_l2fc_wColor(iec.orig, all.degs)
      int.lfc <- string_db$get_interactions( iec.lfc$STRING_id ) #should be same as above, but running interactions here to be safe
      lfc_id <- string_db$post_payload( iec.lfc$STRING_id, colors=iec.lfc$color )
      print(length(which(int.lfc$combined_score>900)))
      #extended
      png(paste("plots/manFig-3C-D_stringDB_",pt.name,"_halo-l2FC_EXTENDED_network.png",sep=""), res=plot_res, units="in", height=5, width=5)
      string_db$plot_network( int.lfc[which(int.lfc$combined_score>900),], payload_id=lfc_id )
      graphics.off()
      
    
    }#end inflection point for loop
  }
  
  ###
  # Figure 4C-D: enrichment barplot
  ###
  {
    transPt.cols <- c("transPt.c2" = "#74cf75", "transPt.c4"="#f5a45d")
    ### c2 transPt ###
    enr.pt.c2 <- read.table("CML_contribution_tables/transPt.CML.contribution_c2.transPt_HallmarkEnrichR_table.tsv", sep="\t", header=T)
    plimit <- 0.01
    sel.pt.c2 <- which(enr.pt.c2$Adjusted.P.value < plimit  )
    sig.nes <- data.frame("transPt.c2" = enr.pt.c2[["CML_contribution"]][sel.pt.c2] )
    rownames(sig.nes) <- enr.pt.c2[["Term"]][sel.pt.c2]
    sig.df <- data.frame(sig.nes)
    state.df <- data.frame("pathway" = rownames(sig.df), "contribution" = sig.df$transPt.c2, 
                           "state" = rep("transPt.c2", dim(sig.df)[1]) )
    state.df$pathway <-  unlist(lapply(state.df$pathway, function(x) str_to_title(gsub("_"," ", gsub( "HALLMARK_", "", x)) )) )
    # line_positions holds divider line between each variable
    state.df2 <- state.df %>% mutate(line_positions = as.numeric(factor(pathway, levels = unique(pathway))), 
                                     line_positions = line_positions + .5,  
                                     line_positions = ifelse(line_positions == max(line_positions), NA, line_positions)) 
    png(paste(manFig_out,"/Fig-4C_fgsea-CML-contribution_p-",plimit,"_bar.png",sep=""), res=plot_res, units="in", height=4, width=4)
    p <- ggplot(state.df2, aes(y=reorder(pathway, contribution), x=contribution, fill=state)) + 
      geom_bar(stat = "identity", position=position_dodge(width = 0.5), width=.5) + 
      theme_bw(base_size=12) + scale_fill_manual(values=transPt.cols) +
      theme( legend.position="none", axis.text.y = element_text(color = "black", size = 14)) +
      scale_x_continuous(n.breaks = 3)
    # geom_vline(aes(xintercept = line_positions))
    print(p)
    graphics.off()
    
    ### c4 transPt ###
    enr.pt.c4 <- read.table("CML_contribution_tables/transPt.CML.contribution_c4.transPt_HallmarkEnrichR_table.tsv", sep="\t", header=T)
    plimit <- 0.01
    sel.pt.c4 <- which(enr.pt.c4$Adjusted.P.value < plimit  )
    sig.nes <- data.frame("transPt.c4" = enr.pt.c4[["CML_contribution"]][sel.pt.c4] )
    rownames(sig.nes) <- enr.pt.c4[["Term"]][sel.pt.c4]
    sig.df <- data.frame(sig.nes)
    state.df <- data.frame("pathway" = rownames(sig.df), "contribution" = sig.df$transPt.c4, 
                           "state" = rep("transPt.c4", dim(sig.df)[1]) )
    state.df$pathway <-  unlist(lapply(state.df$pathway, function(x) str_to_title(gsub("_"," ", gsub( "HALLMARK_", "", x)) )) )
    # line_positions holds divider line between each variable
    state.df2 <- state.df %>% mutate(line_positions = as.numeric(factor(pathway, levels = unique(pathway))), 
                                     line_positions = line_positions + .5,  
                                     line_positions = ifelse(line_positions == max(line_positions), NA, line_positions)) 
    png(paste(manFig_out,"/Fig-4D_fgsea-CML-contribution_p-",plimit,"_bar.png",sep=""), res=plot_res, units="in", height=5, width=5)
    p <- ggplot(state.df2, aes(y=reorder(pathway, contribution), x=contribution, fill=state)) + 
      geom_bar(stat = "identity", position=position_dodge(width = 0.5), width=.5) + 
      theme_bw(base_size=12) + scale_fill_manual(values=transPt.cols) +
      theme( legend.position="none", axis.text.y = element_text(color = "black", size = 14), axis.text.x = element_text(color = "black", size = 10, angle=60, hjust=.5, vjust=.5)) +
      scale_x_continuous(n.breaks = 3)
      #scale_x_continuous(breaks = seq(0, 0.015, length.out = 3))
    # geom_vline(aes(xintercept = line_positions))
    print(p)
    graphics.off()
    
  }
  
  ### SUPPLEMENTAL ###
  
  ###
  # Fig S5: CML contribution of gene modules
  ###
  {
    b.sel <- which(colData(dat.se)$Group=="B")
    gene.sets <- list( early.deg.genes.id, trans.deg.genes.id, late.deg.genes.id, 
                       ig.b1.c, ig.b2.c, ig.b3.c,
                       ig.b3.b2, ig.apost.c, ig.dpost.c, ig.drx.c)
    set.names <- c("Early-DEGs", "Transition-DEGs", "Late-DEGs", 
                   'c1-DEGs', "c2-DEGs", "c3-DEGs" ,
                   "c3.vs.c2-DEGs", "A-postRx-DEGs", "D-Rx-DEGs", "D-postRx-DEGs")
    deg_comp <- list( "B1.vs.C", "B2.vs.C", "B3.vs.C", 
                      "B1.vs.C", "B2.vs.C", "B3.vs.C",
                      "B3.vs.B2", "Apost.vs.C", "Dpost.vs.C", "Drx.vs.C")
    knums <- c(4, 2, 3, 
               4, 4, 4, 
               4, 4, 4, 4)
    #for (si in 1:length(gene.sets) ) {
    for (si in 1:3 ) {
      curgenes <- gene.sets[[si]] #holds gene (unique) IDs 
      set.name <- set.names[si]
      set.sel <- match_genes(rownames(dat.se), curgenes)
      knum <- knums[si]
      #samp.dat <- scale( t(log(data.matrix(assay(dat.se, "abundance"))+amin)), scale=F )[b.sel, set.sel]
      samp.dat <- t(data.matrix(assay(dat.se, "almc")))[b.sel, set.sel]
      cgenes <- rowData(dat.se)$gene_name[set.sel] #holds gene names
      g.info <- rowData(dat.se[set.sel])
      basemean <- rowMeans(log(data.matrix(assay(dat.se[ set.sel, b.sel], "abundance"))+amin))
      g.info[["mean.log.abd"]] <- basemean
      colnames(samp.dat) <- cgenes
      scor <- cor( samp.dat )
      inc <- 20
      blist <- seq(-1, 1, length.out=inc+1)
      corcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(inc)
      phem <- pheatmap::pheatmap(scor, scale="none", color=corcol, breaks=blist, show_rownames=F, show_colnames = F)
      tcutem <- cutree(phem$tree_row, k=knum) #or h=5
      tcutem <- as.character(tcutem)
      names(tcutem) <- cgenes
      tanndf <- data.frame( "Groups" = as.character(tcutem) )
      colnames(tanndf) <- "Group_Number"
      colnames(scor) <- cgenes
      trbcol <- grDevices::rainbow(length(unique(tcutem))) 
      tanncol <- trbcol
      names(tanncol) <- as.numeric(unique(tcutem))
      tanncol <- list(Group_Number = tanncol)
      rownames(scor) <- rownames(tanndf) #match rownames; use numeric to avoid duplicates
      graphics.off()
      png(paste(manFig_out,"/Fig-S4_DEG-dynamics_set-",set.name,"_pheatmap.png",sep=""), res=plot_res, units="in", height=8, width=10) 
      print(pheatmap(scor, show_rownames=F, show_colnames=F, breaks=blist, color=corcol, 
                     annotation_row=tanndf, annotation_names_row = F, annotation_color=tanncol) )#, cellheight=wd, cellwidth=wd)
      graphics.off()
      for (gr in unique(tcutem) ) {
        csel <- which(tcutem==gr)
        cdat <- samp.dat[,csel]
        gr.info <- g.info[csel,]
        c.cor <- scor[csel, csel]
        usel <- upper.tri(c.cor, diag=F)
        cm.cor <- mean(c.cor[usel])
        med.cor <- median(c.cor[usel])
        png(paste(plot_out,"/DEG-dynamics_set-",set.name,"_Grp-",gr,"_phatmap.png",sep=""), res=plot_res, units="in", height=5, width=6) 
        print(pheatmap(c.cor, show_rownames=F, show_colnames=F, breaks=blist, color=corcol, 
                       main=paste("Average correlation: ",round(cm.cor, digits=4),"\n",
                                  "Median correlation: ",round(med.cor, digits=4),sep="")) ) 
        graphics.off()
        write.table(gr.info, paste(plot_out,"/DEG-dynamics_set-",set.name,"_Grp-",gr,"_geneInfo.tsv",sep=""),sep="\t", "row.names"=F )
        #mean expression
        cmean <- rowMeans(cdat)
        m.df <- data.frame(cmean)
        colnames(m.df) <- "meanExp"
        m.df[["sample"]] <- rownames(cdat)
        c.info <- data.frame(colData(dat.se)[b.sel,])
        cm.df <- merge(m.df, c.info, by="sample")
        ### NEEDED ###
        #  add in controls: only controls plot AND plot with both CML and controlsl; except they'll all be smoothed together on the x axis
        for (grp in c("A","C","D")) {
          cur.sel <- which(colData(dat.se)$Group==grp)
          alt.genes <- gr.info$gene_id
          alt.sel <- match_genes(rowData(dat.se)$gene_id, alt.genes) #get genes using "gene_id"
          #alt.dat <- scale( t(log(data.matrix(assay(dat.se, "abundance"))+amin)), scale=F )[cur.sel, alt.sel ]
          alt.dat <- t(data.matrix(assay(dat.se, "almc")))[cur.sel, alt.sel ]
          altmean <- rowMeans(alt.dat)
          a.df <- data.frame(altmean)
          colnames(a.df) <- "meanExp"
          a.df[["sample"]] <- rownames(alt.dat)
          alt.info <- data.frame(colData(dat.se)[cur.sel,])
          alt.df <- merge(a.df, alt.info, by="sample")
          #C data
          ctlsel <- which(colData(dat.se)$Group=="C")
          ctl.dat <- t(data.matrix(assay(dat.se, "almc")))[cur.sel, ctlsel ]
          ctlmean <- rowMeans(ctl.dat)
          ctl.df <- data.frame(ctlmean)
          colnames(ctl.df) <- "meanExp"
          ctl.df[["sample"]] <- rownames(ctl.dat)
          ctl.info <- data.frame(colData(dat.se)[ctlsel,])
          ctl.df <- merge(ctl.df, ctl.info, by="sample")

          colnames(alt.df)
          colnames(cm.df)
          merge.df <- rbind(cm.df, alt.df)
          merge2.df <- rbind(merge.df, ctl.df)

        ### loading value analysis ###
        # !!!!!!!
        # !!! commented out to save time while regenerating plots!!!
        # !!!!!!!
        deg.sets <- c(deg_comp[[si]], "Apost.vs.C", "Dpost.vs.C", "Drx.vs.C" )
        for (curdeg in deg.sets ) {
          #curdeg <- deg_comp[[si]] #NEEDED:change curpath to curdeg
          cgenes.id <- gr.info$gene_id
          curcomp <- curdeg
          sel <- match_genes( gene.df$gene_id, cgenes.id)
          c.df <- gene.df[sel,]

          m.1 <- c.df[[paste(curcomp,"_mean.",1,sep="")]]
          m.2 <- c.df[[paste(curcomp,"_mean.",2,sep="")]]
          c.diff <- m.1 - m.2
          c.ld <- c.df$eigengenes
          up.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]
          dn.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]
          c.all <- c.diff * c.ld #get contribution for all genes in pathway
          c.up <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
          c.down <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
          ctab <- matrix( c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", sum(c.ld), sum(up.ld), sum(dn.ld), sum(c.all,na.rm=T), sum(c.up), sum(c.down) ), nrow=3, ncol=5)
          #add lfc contribution; need to figure out why "Inf" is produced
          c.lfc <- log2((m.1+amin)/(m.2+amin) )
          cl.all <- c.ld * c.lfc #get contribution for all genes in pathway based on LFC
          cl.up <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
          cl.down <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
          ctab <- matrix(c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down",
                            length(which(!is.na(c.all))), length(which(!is.na(c.up))), length(which(!is.na(c.down))), # count of genes contributing (non NA genes)
                            sum(c.ld), sum(up.ld), sum(dn.ld), #eigenvalue
                            sum(c.all,na.rm=T), sum(c.up, na.rm=T), sum(c.down, na.rm=T), #contribution w/ exp diff
                            sum(cl.all,na.rm=T), sum(cl.up, na.rm=T), sum(cl.down, na.rm=T) ), # contrib. w/ log2FC
                         nrow=3, ncol=7)
          if (length(which(is.na(ctab))) != 0) {
            print(paste(":::ERROR::: NAs produced for comp: ",curcomp," pathway: ",curdeg,sep=""))
            print(ctab)
            stop("NA encountered")
          }
          cont.tab <- rbind(cont.tab, ctab)
          #plot space
          c.df[["exp.diff"]] <- c.diff
          c.df[["cont.all"]] <- c.all
          dif <- c.df$exp.diff
          # dif[which(c.df$exp.diff>0)] <- "up"
          # dif[which(c.df$exp.diff<0)] <- "down"
          dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]>0)] <- "up"
          dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]<0)] <- "down"
          c.df[["diff.lab"]] <- dif
          #build loading contribution
          "eigen.up" = mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)])
          "eigen.down" <- mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
          "PC1.up" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)])
          "PC1.down" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])

          ### make "eigen-contribution"
          # positive value = positive contribution (neg eig + pos exp OR pos eig + neg exp)
          # neg contribution (neg eig + neg exp OR pos eig + pos exp)
          cml.c <- c.df$eigengenes

          cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 #pro-CML up (neg eig. -> adjust)
          cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 #pro-CML down ( pos eig. -> NO adjust)
          cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 # anti-CML down (neg -> NO adjust)
          cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 # anti-CML up (pos -> adjust)

          #sums
          neg.up <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
          neg.up.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
          pos.down <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
          pos.down.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
          neg.down <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
          neg.down.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
          pos.up <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
          pos.up.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
          cmat <- matrix( c( rep(curdeg, 4), rep(curcomp, 4), "neg|up", "pos|down", "neg|down", "pos|up",
                             neg.up, pos.down, neg.down, pos.up,
                             neg.up.l, pos.down.l, neg.down.l, pos.up.l),
                          nrow=4, ncol=5)
          cml.tab <- rbind(cml.tab, cmat )
          #cml label
          cml.lab <- rep(NA, length(cml.c))
          cml.lab[which(cml.c > 0)] <- "pro-CML"
          cml.lab[which(cml.c < 0)] <- "anti-CML"
          #full matrix
          # :::note::: the mean contribution of each class of genes must be added!
          #(neg.up  / neg.up.l) + (pos.down / pos.down.l), (neg.down / neg.down.l) + (pos.up / pos.up.l)
          n.u.m <- (neg.up  / neg.up.l)
          p.d.m <- (pos.down / pos.down.l)
          n.d.m <- (neg.down / neg.down.l)
          p.u.m <- (pos.up / pos.up.l)
          fmat <- matrix( c( rep(curdeg, 2), rep(curcomp, 2), "pro-CML", "anti-CML",
                             sum(c(n.u.m, p.d.m), na.rm=T), sum(c(n.d.m, p.u.m), na.rm=T),
                             neg.up.l + pos.down.l, neg.down.l + pos.up.l),
                          nrow=2, ncol=5)
          cml.full <- rbind(cml.full, fmat )
          ### eigengene value vs PC1 ###
          #counts
          cl.up <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
          cl.down <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
          al.up <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
          al.down <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
          c.df[["CML_contribution"]] <- cml.c
          cml.up = mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)])
          cml.down <- mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
          PC1.up <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)])
          PC1.down <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
          cml.tot <- mean(cml.c, na.rm=T)
          PC1.tot <- mean(c.df$PC1, na.rm=T)
          #make df for labeling genes without NA genes
          l.df <- c.df[which(c.df$diff.lab=="up" | c.df$diff.lab=="down" ),]
          tup <- sort.int(l.df$CML_contribution, index.return = T, decreasing = T)$ix[seq(1,10)]
          tdown <- sort.int(l.df$CML_contribution, index.return = T, decreasing = F)$ix[seq(1,5)]

          png(paste(manFig_out,"/Fig-S5_DEG-dynamics_set-",set.name,"_Grp-",gr,"_color-",curcomp,"_pts+label+tot.png",sep=""), res=plot_res, units="in", height=4, width=4)
          p <- ggplot(c.df[which(!is.na(c.df$diff.lab)),], aes(x=CML_contribution, y= PC1, fill=diff.lab, color=diff.lab )) +
            geom_point(aes(shape=diff.lab), alpha=.75) + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", size=1) +
            # geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.2, 
            #              arrow = arrow(length = unit(.25, "cm"), type="open"), lineend="butt", linejoin="mitre") +
            # geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.2,
            #              arrow = arrow(length = unit(.25, "cm"), type="open"), lineend="butt", linejoin="mitre") +
            geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.25,
                         arrow = arrow(length = unit(.30, "cm"), type="closed"), lineend="butt", linejoin="mitre") +
            scale_color_manual(values= list("up"="#f76a95", "down" = "#4287f5" )) + scale_fill_manual(values= list("up"="#f76a95", "down" = "#4287f5" )) +
            # geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene), size=14) +
            # geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene), size=14) +
            theme(legend.position = "none") + scale_shape_manual(values=c("up"=24, "down"=25)) +
            scale_x_continuous(limits=c(-.03, 0.058)) + scale_y_continuous(limits=c(-.025, 0.018))
            # geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", linewidth=1.25) +
            # geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.5,
            #              arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
            # scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
            # # geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene)) +
            # # geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene)) +
            # ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") ) +
            # theme(legend.position="none")
          #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" ))
          print(p)
          graphics.off()
          #all names

 

        } #end loading value analysis
      } #end loop over gene modules
    }
  }
  
}

  ###
  # Fig S6: inflection points CML contribution
  ###
  {
    
    gene_set_CML_contribution_S6 <- function( out_name, deg_list, deg_comp, 
                                           gene.df = gene.df, enrichOut = T, all.gene.labs=F ) { #all
      
      curdeg <- out_name #deg_name[[i]] #NEEDED:change curpath to curdeg
      cgenes.id <- deg_list
      curcomp <- deg_comp
      
      sel <- match_genes( gene.df$gene_id, cgenes.id)
      c.df <- gene.df[sel,]
      
      m.1 <- c.df[[paste(curcomp,"_mean.",1,sep="")]]
      m.2 <- c.df[[paste(curcomp,"_mean.",2,sep="")]]
      c.diff <- m.1 - m.2
      c.ld <- c.df$eigengenes
      up.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]
      dn.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]
      c.all <- c.diff * c.ld #get contribution for all genes in pathway
      c.up <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
      c.down <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
      ctab <- matrix( c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", sum(c.ld), sum(up.ld), sum(dn.ld), sum(c.all,na.rm=T), sum(c.up), sum(c.down) ), nrow=3, ncol=5) 
      #add lfc contribution; need to figure out why "Inf" is produced
      c.lfc <- log2((m.1+amin)/(m.2+amin) )
      cl.all <- c.ld * c.lfc #get contribution for all genes in pathway based on LFC
      cl.up <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
      cl.down <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
      ctab <- matrix(c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", 
                        length(which(!is.na(c.all))), length(which(!is.na(c.up))), length(which(!is.na(c.down))), # count of genes contributing (non NA genes)
                        sum(c.ld), sum(up.ld), sum(dn.ld), #eigenvalue
                        sum(c.all,na.rm=T), sum(c.up, na.rm=T), sum(c.down, na.rm=T), #contribution w/ exp diff
                        sum(cl.all,na.rm=T), sum(cl.up, na.rm=T), sum(cl.down, na.rm=T) ), # contrib. w/ log2FC
                     nrow=3, ncol=7)
      if (length(which(is.na(ctab))) != 0) {
        print(paste(":::ERROR::: NAs produced for comp: ",curcomp," pathway: ",curdeg,sep=""))
        print(ctab)
        stop("NA encountered")
      }
      cont.tab <- rbind(cont.tab, ctab)
      #plot space
      c.df[["exp.diff"]] <- c.diff
      c.df[["cont.all"]] <- c.all
      dif <- c.df$exp.diff
      # dif[which(c.df$exp.diff>0)] <- "up"
      # dif[which(c.df$exp.diff<0)] <- "down"
      dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]>0)] <- "up"
      dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]<0)] <- "down"
      c.df[["diff.lab"]] <- dif
      #build loading contribution
      "eigen.up" = mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      "eigen.down" <- mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
      "PC1.up" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      "PC1.down" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      
      
      ### make "eigen-contribution"
      # positive value = positive contribution (neg eig + pos exp OR pos eig + neg exp)
      # neg contribution (neg eig + neg exp OR pos eig + pos exp)
      cml.c <- c.df$eigengenes
      
      cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 #pro-CML up (neg eig. -> adjust)
      cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 #pro-CML down ( pos eig. -> NO adjust)
      cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 # anti-CML down (neg -> NO adjust)
      cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 # anti-CML up (pos -> adjust)
      
      #sums
      
      neg.up <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
      neg.up.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
      pos.down <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
      pos.down.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
      neg.down <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
      neg.down.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
      pos.up <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
      pos.up.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
      cmat <- matrix( c( rep(curdeg, 4), rep(curcomp, 4), "neg|up", "pos|down", "neg|down", "pos|up", 
                         neg.up, pos.down, neg.down, pos.up,
                         neg.up.l, pos.down.l, neg.down.l, pos.up.l),
                      nrow=4, ncol=5) 
      cml.tab <- rbind(cml.tab, cmat )
      #cml label
      cml.lab <- rep(NA, length(cml.c))
      cml.lab[which(cml.c > 0)] <- "pro-CML"
      cml.lab[which(cml.c < 0)] <- "anti-CML"
      #full matrix
      # :::note::: the mean contribution of each class of genes must be added!
      #(neg.up  / neg.up.l) + (pos.down / pos.down.l), (neg.down / neg.down.l) + (pos.up / pos.up.l)
      n.u.m <- (neg.up  / neg.up.l) 
      p.d.m <- (pos.down / pos.down.l)
      n.d.m <- (neg.down / neg.down.l) 
      p.u.m <- (pos.up / pos.up.l)
      fmat <- matrix( c( rep(curdeg, 2), rep(curcomp, 2), "pro-CML", "anti-CML",
                         sum(c(n.u.m, p.d.m), na.rm=T), sum(c(n.d.m, p.u.m), na.rm=T),
                         neg.up.l + pos.down.l, neg.down.l + pos.up.l),
                      nrow=2, ncol=5) 
      cml.full <- rbind(cml.full, fmat )
      ### eigengene value vs PC1 ###
      #counts
      cl.up <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
      cl.down <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
      al.up <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
      al.down <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
      c.df[["CML_contribution"]] <- cml.c
      cml.up = mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      cml.down <- mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
      PC1.up <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      PC1.down <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      cml.tot <- mean(cml.c, na.rm=T) 
      PC1.tot <- mean(c.df$PC1, na.rm=T)
      #make df for labeling genes without NA genes
      l.df <- c.df[which(c.df$diff.lab=="up" | c.df$diff.lab=="down" ),]
      tup <- sort.int(l.df$CML_contribution, index.return = T, decreasing = T)$ix[seq(1,10)]
      tdown <- sort.int(l.df$CML_contribution, index.return = T, decreasing = F)$ix[seq(1,5)]

      png(paste(manFig_out,"/Fig-S6_eigenDEGs_",curdeg,"_comparison-",curcomp,"_pts+label+tot.png",sep=""), res=plot_res, units="in", height=4, width=4)
      ### OLD ###
      # p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab )) + 
      #   geom_point() + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", linewidth=1.25) +
      #   geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.5, 
      #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
      #   geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.5,
      #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
      #   geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.5,
      #                arrow = arrow(length = unit(.35, "cm"), type="closed"), lineend="butt", linejoin="bevel") +
      #   scale_color_manual(values= list("up"="#f54278", "down" = "#4287f5" )) +
      #   geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene)) +
      #   geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene)) +
      #   ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
      #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" )) 
      p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab, fill=diff.lab )) + 
        geom_point(aes(shape=diff.lab), alpha=.75) + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", size=1) +
        # geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.2, 
        #              arrow = arrow(length = unit(.25, "cm"), type="open"), lineend="butt", linejoin="mitre") +
        # geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.2,
        #              arrow = arrow(length = unit(.25, "cm"), type="open"), lineend="butt", linejoin="mitre") +
        geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.25,
                     arrow = arrow(length = unit(.30, "cm"), type="closed"), lineend="butt", linejoin="mitre") +
        scale_color_manual(values= list("up"="#f76a95", "down" = "#4287f5" )) + scale_fill_manual(values= list("up"="#f76a95", "down" = "#4287f5" )) +
        # geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene), size=14) +
        # geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene), size=14) +
        theme(legend.position = "none") + scale_shape_manual(values=c("up"=24, "down"=25)) +
        scale_x_continuous(limits=c(-.03, 0.058)) + scale_y_continuous(limits=c(-.025, 0.018)) 
      print(p)
      graphics.off()
      #all names

      
    } #end gene set contribution function
    output_enrichr_manFigS6 <- function(inres, name, path_out) {
      for ( db in names(inres)) {
        if (dim(inres[[db]])[1]==0) {next} #skip output
        write.table(inres[[db]], paste(path_out,"/",db,"_",name,"_table.tsv", sep="" ), sep="\t", row.names=F )
        png(paste(path_out,"/Fig-S6_",db,"_",name,"_plot.png", sep="" ), res=plot_res, units="in", height=12, width=8 )
        p <- plotEnrich(inres[[db]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
        print(p)
        graphics.off()
      }
    }
    
    #load genes
    iec.c2.tab <- read.table(paste("plots/DEG_potential-expression_correlation_c2/corrTable_top-genes-c2.tsv",sep=""), sep="\t", header=T)
    iec.c4.tab <- read.table(paste("plots/DEG_potential-expression_correlation_c4/corrTable_top-genes-c4.tsv",sep=""), sep="\t", header=T)
    
    ### venn of overlap ###
    plotVenn( list(iec.c2.tab$gene, iec.c4.tab$gene), sNames=c("c2", "c4"), showPlot = T, 
              outFile=paste(manFig_out,"/Fig-S6_venn2_inflection-pt_topGenes_c2.vs.c4.svg",sep=""), labelRegions=F,
              setColors=c("#5ee06b", "#f29c33"))
    rsvg::rsvg_png(
      paste(manFig_out,"/Fig-S6_venn2_inflection-pt_topGenes_c2.vs.c4.svg",sep=""), paste(manFig_out,"/Fig-S6_venn2_inflection-pt_topGenes_c2.vs.c4.png",sep=""),
      width = 900, height = 700)
    write.table(intersect(iec.c2.tab$gene, iec.c4.tab$gene), paste(manFig_out,"/Fig-S6_venn2_inflection-pt_topGenes_c2.vs.c4_shared.tsv",sep=""), sep="\t", row.names=F, col.names = F) 
    
    #human: 9606; mouse: 10090
    string_db <- STRINGdb$new( species = 10090, score_threshold = 0.4, version="11.5")
    
    ### process both inflection points ###
    #setup biomaRt
    mm.ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl") #"hsapiens_gene_ensembl" For Homo sapiens
    attributes <- c("ensembl_gene_id", "external_gene_name")
    #setup variables    
    share.tab <- data.frame("gene"=intersect(iec.c2.tab$gene,iec.c4.tab$gene))
    pt.tabs <- list("c2" = iec.c2.tab, "c4" = iec.c4.tab)
    trans.pt <- list() #list to hold genes

    for (pt.name in c( "c2",  "c4") ) {
      print(paste("processing ",pt.name, sep=""))
      cordegs <- paste("plots/DEG_potential-expression_correlation_",pt.name,sep="")
      dir.create(cordegs, showWarnings = F)
      iec.orig <- string_db$map( pt.tabs[[pt.name]], "gene", removeUnmappedRows = TRUE )
      
      ### add log2FC to table ### 
      iec <- add_l2fc(iec.orig, all.degs)
      #extended
      int <- string_db$get_interactions( iec$STRING_id )
      
      #add gene expression colors
      # log2FC coloring on network
      iec.lfc <- add_l2fc_wColor(iec.orig, all.degs)
      int.lfc <- string_db$get_interactions( iec.lfc$STRING_id ) #should be same as above, but running interactions here to be safe
      lfc_id <- string_db$post_payload( iec.lfc$STRING_id, colors=iec.lfc$color )
      #extended
      if (pt.name == "c2") { 
        png(paste(manFig_out,"/Fig-S6_allDEGs_stringDB_",pt.name,"_halo-l2FC_EXTENDED_network.png",sep=""), res=plot_res, units="in", height=10, width=10)
        string_db$plot_network( int.lfc[which(int.lfc$combined_score>900),], payload_id=lfc_id )
        graphics.off()
      } else { #reduce number of proteins for larger c4 inflection pt
        png(paste(manFig_out,"/Fig-S6_allDEGs_stringDB_",pt.name,"_halo-l2FC_EXTENDED_network.png",sep=""), res=plot_res, units="in", height=10, width=10)
        string_db$plot_network( iec.orig$STRING_id, payload_id=lfc_id, required_score=900 )
        graphics.off()
      }
      # Get the mapping of StringDB IDs to gene names
      stringdb_ids <- unique(unlist(iec.orig$STRING_id)) 
      stringdb_ids <- unlist(lapply(stringdb_ids, function(x) gsub("10090.","",x)))
      gene_mapping <- getBM(attributes = attributes, filters = "ensembl_peptide_id", values = stringdb_ids, mart = mm.ensembl)
      genes <- gene_mapping[,2]
      trans.pt[[pt.name]] <- list("full" = genes)
      #extended
      stringdb_ids <- unique(unlist(int.lfc[which(int.lfc$combined_score>900),c(1,2)])) 
      stringdb_ids <- unlist(lapply(stringdb_ids, function(x) gsub("10090.","",x)))
      gene_mapping <- getBM(attributes = attributes, filters = "ensembl_peptide_id", values = stringdb_ids, mart = mm.ensembl)
      genes_ext <- gene_mapping[,2]
      trans.pt[[pt.name]] <- list("full" = genes, "extended" = genes_ext)
      
      ## CML contribution plots - Figure S6 ##
      #get gene IDs
      gsel <- match_all_genes( rowData(dat.se)$gene_name, genes)
      gene_ids <- rowData(dat.se)$gene_id[gsel]
      if (pt.name == "c2") {
        gene_set_CML_contribution_S6(paste(pt.name,".trans.pt",sep=""), gene_ids, "B2.vs.B1", gene.df = gene.df, enrichOut = F, all.gene.labs=F)
      }
      if (pt.name == "c4") {
        gene_set_CML_contribution_S6(paste(pt.name,".trans.pt",sep=""), gene_ids, "B3.vs.B2", gene.df = gene.df, enrichOut = F, all.gene.labs=F)
      }
      gene_set_CML_contribution_S6(paste(pt.name,".trans.pt",sep=""), gene_ids, "B3.vs.C", gene.df = gene.df, enrichOut = F, all.gene.labs=F )
      #extended
      gsel <- match_all_genes( rowData(dat.se)$gene_name, genes_ext)
      gene_ids <- rowData(dat.se)$gene_id[gsel]
      if (pt.name == "c2") {
        gene_set_CML_contribution_S6(paste("EXTENDED_",pt.name,".trans.pt",sep=""), gene_ids, "B2.vs.B1", gene.df = gene.df, enrichOut = F, all.gene.labs=F)
      }
      if (pt.name == "c4") {
        gene_set_CML_contribution_S6(paste("EXTENDED_",pt.name,".trans.pt",sep=""), gene_ids, "B3.vs.B2", gene.df = gene.df, enrichOut = F, all.gene.labs=F)
      }
      gene_set_CML_contribution_S6(paste("EXTENDED_",pt.name,".trans.pt",sep=""), gene_ids, "B3.vs.C", gene.df = gene.df, enrichOut = F, all.gene.labs=F )
      
      ### clustering ###
      {  
        clustersList <- string_db$get_clusters(iec.orig)
        cl.len <- min( c(length(clustersList), 4) )
        png(paste(manFig_out,"/Fig-S6_allDEGs_stringDB_",pt.name,"_halo-l2FC_EXTENDED_network.png",sep=""), res=plot_res, units="in", height=10, width=10)
        if (pt.name=="c4") { #only larger c4 point is clustered
          png(paste(manFig_out,"/Fig-S6_allDEGs_stringDB_",pt.name,"_clusters_halo-l2FC_network-lg-min-900.png",sep=""), res=plot_res, units="in", height=15, width=15)
          par(mfrow=c(2,2))
          for(i in seq(1:cl.len)){
            string_db$plot_network(clustersList[[i]], payload_id=lfc_id, required_score = 900)
          }
          graphics.off()
        }
        
        ## enrichment + CML contribution##
        for(i in seq(1:cl.len)) {
          c.enr <- string_db$get_enrichment( clustersList[[i]] )
          
          # enrichR #
          stringdb_ids <- unique(unlist( clustersList[[i]] )) 
          stringdb_ids <- unlist(lapply(stringdb_ids, function(x) gsub("10090.","",x)))
          gene_mapping <- getBM(attributes = attributes, filters = "ensembl_peptide_id", values = stringdb_ids, mart = mm.ensembl)
          
          genes <- gene_mapping[,2]
          trans.pt[[pt.name]][[i]] = genes
          
          cur_path_out <- manFig_out
          dir.create(cur_path_out, showWarnings = F)
          enriched.up <- enrichr(genes, test_dbs ) #up only because nearly everything is upregulated
          output_enrichr_manFigS6(enriched.up, paste(pt.name,"_c-",i,sep=""), cur_path_out)
           
          ### CML contribution - Figure S6 ###
          #get gene IDs
          gsel <- match_all_genes(rowData(dat.se)$gene_name, genes)
          gene_ids <- rowData(dat.se)$gene_id[gsel]
          if (pt.name == "c2") {
            gene_set_CML_contribution_S6(paste(pt.name,".trans.pt_cluster-",i,sep=""), gene_ids, "B2.vs.B1", gene.df = gene.df, enrichOut = F, all.gene.labs=F)
          }
          if (pt.name == "c4") {
            gene_set_CML_contribution_S6(paste(pt.name,".trans.pt_cluster-",i,sep=""), gene_ids, "B3.vs.B2", gene.df = gene.df, enrichOut = F, all.gene.labs=F)
          }
          gene_set_CML_contribution_S6(paste(pt.name,".trans.pt_cluster-",i,sep=""), gene_ids, "B3.vs.C", gene.df = gene.df, enrichOut = F, all.gene.labs=F)
          
        }
        
      } #end clustering 

    }#end inflection point for loop
    
    
    
    
  }#end Fig S6
  
  
}#end Fig 4 plots



###
### Figure 5: treatment analysis (Fig. S7)
###
{
  ###
  # Figure 5A - TREATMENT RESPONSE BOXPLOT
  ###
  {
    ### trajectories ###
    png(paste(manFig_out,"/Fig-5A_Group-A_time.vs.CML_space_line+point.png",sep=""), res=plot_res, units="in", height=3, width=4)
    p <- ggplot(col.df, aes(x=sample_weeks, y=CML_space)) + geom_point(color="grey", alpha=.5) + geom_hline(yintercept=ss.cps[c(2,4)], linetype="dashed")+
      geom_line(data=col.df[which(col.df$Group=="A"),], aes(x=sample_weeks, y=CML_space, color=Group, group=mouse_id)) +
      geom_point(data=col.df[which(col.df$Group=="A"),], aes(x=sample_weeks, y=CML_space, color=Group, group=mouse_id)) +
      theme_bw(base_size=18) + theme(legend.position="none") + coord_cartesian(ylim = c(-350, 225) ) + scale_color_manual(values=group.cols)
    print(p)
    graphics.off()
    
    png(paste(manFig_out,"/Fig-5A_Group-D_time.vs.CML_space_line+point.png",sep=""), res=plot_res, units="in", height=3, width=4)
    p <- ggplot(col.df, aes(x=sample_weeks, y=CML_space )) + geom_point(color="grey", alpha=.5) + geom_hline(yintercept=ss.cps[c(2,4)], linetype="dashed") +
      geom_line(data=col.df[which(col.df$Group=="D"),], aes(x=sample_weeks, y=CML_space, color=Group, group=mouse_id)) +
      geom_point(data=col.df[which(col.df$Group=="D"),], aes(x=sample_weeks, y=CML_space, color=Group, group=mouse_id)) +
      theme_bw(base_size=18) + theme(legend.position="none") + coord_cartesian(ylim = c(-350, 225) ) +  scale_color_manual(values=group.cols) 
    print(p)
    graphics.off()
    
    
    ### boxplot ###
    trt.df <- data.frame(colData(dat.se)[which(colData(dat.se)$Group=="C" | 
                                                 colData(dat.se)$rx_class == "A.post" | 
                                                 colData(dat.se)$rx_class == "D.Rx" | 
                                                 colData(dat.se)$rx_class == "D.post"),])
    lab <- trt.df$rx_class
    lab[which(is.na(lab))] <- "C"
    trt.df[["label"]] <- lab
    trt.df[["label"]] <- factor(trt.df[["label"]], levels=c("C", "A.post", "D.Rx", "D.post"))
    cols <- c("TOTO-postRx" = "#45e9f5", "TKI-postRx" = "#d955c9", "TKI"="#fca7cf", "TKI-Rx" = "#6507b3")
    #trt.cols <- list("C" = "black", "A.post"= "#45e9f5", "D.Rx" = "#6507b3", "D.post"="#d955c9")
    trt.cols <- list("C" = "black", "A.post"= "#45e9f5", "D.Rx" = "#6507b3", "D.post"="#d955c9")
    png(paste(manFig_out,"/Fig-4A_treatmentResponse_space-boxplot.png",sep=""), res=plot_res, units="in", height=3, width=4)
    p <- ggplot(trt.df, aes(x=label, y=CML_space, fill=label )) + geom_hline(yintercept=ss.cps[c(2,4)], linetype="dashed") + 
      geom_boxplot() + theme_bw(base_size=18) + 
      theme(axis.title.x = element_blank(), axis.text.x=element_blank()) + coord_cartesian(ylim = c(-350, 225) ) +
      scale_fill_manual(values=trt.cols) + theme(legend.position="none") 
    print(p)
    graphics.off()
  }
  
  ###
  # Figure 5B - upset
  ###
  {
    ## Rx- & CML-vs-ctl
    
    #limit combinations #
    comb.cols <- c( "TOTO-postRx" = "#45e9f5", "TKI-postRx" = "#d955c9", "TKI"="#fca7cf", "TKI-Rx" = "#6507b3", "mixed_Rx"="#014508", "CML.vs.ctl"="grey")
    selcm <- length(which(comb_size(cm.Rx) > 10))
    cm.Rx.lim <- cm.Rx[comb_size(cm.Rx) > 10]
    comb.sort <- sort.int(comb_size(cm.Rx.lim), index.return = T)$ix
    comb.annot <- c()
    set.ord <- c("c5.vs.ctl", "c3.vs.ctl", "c1.vs.ctl", "TOTO_postRx", "TKI_Rx", "TKI_postRx")
    for (ni in 1:length(comb_name(cm.Rx.lim)) ){
      n <- comb_name(cm.Rx.lim)[ni]
      n.list <- as.numeric(strsplit(n, "")[[1]])
      #name order: "c5", "c3", "c1", "TKI.Rx", "TKI-postRx", "TOTO-postRx"
      if (n.list[4]==1 & sum(n.list[4:6])==1) { #TKI RX only
        comb.annot <- c(comb.annot, "TKI-Rx")
      } else if (n.list[5]==1 & sum(n.list[4:6])==1) { #TKI postRx only
        comb.annot <- c(comb.annot, "TKI-postRx")
      } else if (n.list[6]==1 & sum(n.list[4:6])==1) { #TOTO postRx only
        comb.annot <- c(comb.annot, "TOTO-postRx")
      } else if  (n.list[5]==1 & n.list[4]==1  & sum(n.list[4:6])==2) { #any TKI
        comb.annot <- c(comb.annot, "TKI")
      } else if  (n.list[6]==1 & sum(n.list[4:6])>=2) { #mixed Rx
        comb.annot <- c(comb.annot, "mixed_Rx")
      } else {
        comb.annot <- c(comb.annot, "CML.vs.ctl")
      }
    }
    png(paste(manFig_out,"/Fig-5B_vsControl+Rx_upset-10-limit.png",sep=""), res=plot_res, units="in", height=3, width=5)
    UpSet(cm.Rx.lim, set_order = set.ord, comb_order = rev(comb.sort),
          right_annotation = upset_right_annotation(cm.Rx, add_numbers = TRUE,
                                                    gp = gpar(fill = c("red", "goldenrod", "dodgerblue", "#6507b3", "#d955c9",  "#45e9f5" )) ),
    )
    # top_annotation = HeatmapAnnotation(
    #   Treatment = comb.annot,
    #   col = list(Treatment=comb.cols),
    #   "Intersection\nsize" = anno_barplot(comb_size(cm.Rx.lim), 
    #                                       border = FALSE, 
    #                                       gp = gpar(fill = "black"), 
    #                                       height = unit(2, "cm")
    #   ), 
    # annotation_name_side = "left", 
    # annotation_name_rot = 0) )
    graphics.off()
  }
  
  ###
  # Figure 5C-D: Rx CML contribution 
  ###
  {
    gene_set_CML_contribution_manFig5 <- function( out_name, deg_list, deg_comp, gene.df = gene.df ) { #all
      
      curdeg <- out_name 
      cgenes.id <- deg_list
      curcomp <- deg_comp
      
      sel <- match_genes( gene.df$gene_id, cgenes.id)
      c.df <- gene.df[sel,]
      
      m.1 <- c.df[[paste(curcomp,"_mean.",1,sep="")]]
      m.2 <- c.df[[paste(curcomp,"_mean.",2,sep="")]]
      c.diff <- m.1 - m.2
      c.ld <- c.df$eigengenes
      up.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]
      dn.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]
      c.all <- c.diff * c.ld #get contribution for all genes in pathway
      c.up <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
      c.down <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
      ctab <- matrix( c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", sum(c.ld), sum(up.ld), sum(dn.ld), sum(c.all,na.rm=T), sum(c.up), sum(c.down) ), nrow=3, ncol=5) 
      #add lfc contribution; need to figure out why "Inf" is produced
      c.lfc <- log2((m.1+amin)/(m.2+amin) )
      cl.all <- c.ld * c.lfc #get contribution for all genes in pathway based on LFC
      cl.up <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
      cl.down <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
      ctab <- matrix(c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", 
                        length(which(!is.na(c.all))), length(which(!is.na(c.up))), length(which(!is.na(c.down))), # count of genes contributing (non NA genes)
                        sum(c.ld), sum(up.ld), sum(dn.ld), #eigenvalue
                        sum(c.all,na.rm=T), sum(c.up, na.rm=T), sum(c.down, na.rm=T), #contribution w/ exp diff
                        sum(cl.all,na.rm=T), sum(cl.up, na.rm=T), sum(cl.down, na.rm=T) ), # contrib. w/ log2FC
                     nrow=3, ncol=7)
      if (length(which(is.na(ctab))) != 0) {
        print(paste(":::ERROR::: NAs produced for comp: ",curcomp," pathway: ",curdeg,sep=""))
        print(ctab)
        stop("NA encountered")
      }
      cont.tab <- rbind(cont.tab, ctab)
      #plot space
      c.df[["exp.diff"]] <- c.diff
      c.df[["cont.all"]] <- c.all
      dif <- c.df$exp.diff
      # dif[which(c.df$exp.diff>0)] <- "up"
      # dif[which(c.df$exp.diff<0)] <- "down"
      dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]>0)] <- "up"
      dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]<0)] <- "down"
      c.df[["diff.lab"]] <- dif
      #build loading contribution
      "eigen.up" = mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      "eigen.down" <- mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
      "PC1.up" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      "PC1.down" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      
      ### make "eigen-contribution"
      # positive value = positive contribution (neg eig + pos exp OR pos eig + neg exp)
      # neg contribution (neg eig + neg exp OR pos eig + pos exp)
      cml.c <- c.df$eigengenes
      cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 #pro-CML up (neg eig. -> adjust)
      cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 #pro-CML down ( pos eig. -> NO adjust)
      cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 # anti-CML down (neg -> NO adjust)
      cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 # anti-CML up (pos -> adjust)
      
      #sums
      neg.up <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
      neg.up.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
      pos.down <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
      pos.down.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
      neg.down <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
      neg.down.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
      pos.up <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
      pos.up.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
      cmat <- matrix( c( rep(curdeg, 4), rep(curcomp, 4), "neg|up", "pos|down", "neg|down", "pos|up", 
                         neg.up, pos.down, neg.down, pos.up,
                         neg.up.l, pos.down.l, neg.down.l, pos.up.l),
                      nrow=4, ncol=5) 
      cml.tab <- rbind(cml.tab, cmat )
      #cml label
      cml.lab <- rep(NA, length(cml.c))
      cml.lab[which(cml.c > 0)] <- "pro-CML"
      cml.lab[which(cml.c < 0)] <- "anti-CML"
      #full matrix
      # :::note::: the mean contribution of each class of genes must be added!
      n.u.m <- (neg.up  / neg.up.l) 
      p.d.m <- (pos.down / pos.down.l)
      n.d.m <- (neg.down / neg.down.l) 
      p.u.m <- (pos.up / pos.up.l)
      fmat <- matrix( c( rep(curdeg, 2), rep(curcomp, 2), "pro-CML", "anti-CML",
                         sum(c(n.u.m, p.d.m), na.rm=T), sum(c(n.d.m, p.u.m), na.rm=T),
                         neg.up.l + pos.down.l, neg.down.l + pos.up.l),
                      nrow=2, ncol=5) 
      cml.full <- rbind(cml.full, fmat )
      ### eigengene value vs PC1 ###
      cl.up <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
      cl.down <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
      al.up <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
      al.down <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
      c.df[["CML_contribution"]] <- cml.c
      cml.up = mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      cml.down <- mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
      PC1.up <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      PC1.down <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      cml.tot <- mean(cml.c, na.rm=T) 
      PC1.tot <- mean(c.df$PC1, na.rm=T)
      #make df for labeling genes without NA genes
      l.df <- c.df[which(c.df$diff.lab=="up" | c.df$diff.lab=="down" ),]
      tup <- sort.int(l.df$CML_contribution, index.return = T, decreasing = T)$ix[seq(1,10)]
      tdown <- sort.int(l.df$CML_contribution, index.return = T, decreasing = F)$ix[seq(1,5)]
      png(paste(manFig_out,"/Fig-5CD_eigenDEGs_",curdeg,"_comparison-",curcomp,"_pts+label.png",sep=""), res=plot_res, units="in", height=4, width=4)
      p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab, fill=diff.lab )) + 
        geom_point(aes(shape=diff.lab), alpha=.75) + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", linewidth=1) +
        # geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.2, 
        #              arrow = arrow(length = unit(.25, "cm"), type="open"), lineend="butt", linejoin="mitre") +
        # geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.2,
        #              arrow = arrow(length = unit(.25, "cm"), type="open"), lineend="butt", linejoin="mitre") +
        geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.25,
                     arrow = arrow(length = unit(.30, "cm"), type="closed"), lineend="butt", linejoin="mitre") +
        scale_color_manual(values= list("up"="#f76a95", "down" = "#4287f5" )) + scale_fill_manual(values= list("up"="#f76a95", "down" = "#4287f5" )) +
        # geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene), size=14) +
        # geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene), size=14) +
        theme(legend.position = "none") + scale_shape_manual(values=c("up"=24, "down"=25)) +
        scale_x_continuous(limits=c(-.03, 0.058)) + scale_y_continuous(limits=c(-.025, 0.018)) 

      print(p)
      graphics.off()

      
      
      
    } #end gene set contribution function
    ### TOTO minus age-related genes ###
    ig.toto.time <- setdiff(ig.apost.c, ig.c.lve)
    
    
    ### TOTO and TKI vs Ctrl ###
    # :note: same as Fig. 2D with different for loop indexes AND file name prefix in "gene_set_CML_contribution_manFig4" function

    deg_list <- list( early.deg.genes.id, trans.deg.genes.id, late.deg.genes.id, 
                      ig.b1.c, ig.b2.c, ig.b3.c,
                      ig.b3.b1, ig.b3.b2, ig.b2.b1,
                      ig.apost.c, ig.dpost.c, ig.drx.c, ig.c.evl, ig.c.lve, ig.apost.cL, ig.rx.da, ig.toto.time )
    deg_name <- list("Early",  "Transition", "Late", 
                     "B1.vs.C", "B2.vs.C", "B3.vs.C",
                     "B3.vs.B1", "B3.vs.B2", "B2.vs.B1",
                     "Apost.vs.C", "Dpost.vs.C", "Drx.vs.C", "C_early.vs.late", "C_late.vs.early", "Apost.vs.Clate", "Drx.vs.Apost", "TOTO.minus.time")
    deg_comp <- list( "B1.vs.C", "B2.vs.C", "B3.vs.C", 
                      "B1.vs.C", "B2.vs.C", "B3.vs.C",
                      "B3.vs.B1", "B3.vs.B2", "B2.vs.B1",
                      "Apost.vs.C", "Dpost.vs.C", "Drx.vs.C","C_early.vs.late", "C_late.vs.early", "Apost.vs.Clate", "Drx.vs.Apost", "Apost.vs.C")
    for (i in c(10, 12, 15, 17)){
      gene_set_CML_contribution_manFig5( deg_name[i], deg_list[[i]], deg_comp[i], gene.df=gene.df)
    }
  }
  
  ###
  # Figure 5E: TKI gene set CML contribution
  ###
  {
    
    # CML contribution #
    {
      gsea.drx.c <- read.table("CML_contribution_tables/CML.contribution_D-Rx.vs.C_GSEA_table.tsv", sep="\t", header=T)
      gsea.dpost.c <- read.table("CML_contribution_tables/CML.contribution_D-postRx.vs.C_GSEA_table.tsv", sep="\t", header=T)
      plimit <- 0.0001
      sel.drx.c <- which(gsea.drx.c$padj < plimit & abs(gsea.drx.c$NES) >=2 )
      sel.dpost.c <- which(gsea.dpost.c$padj < plimit & abs(gsea.dpost.c$NES) >=2 )
      sig.nes <- data.frame("TKI.vs.ctl" = gsea.dpost.c[["leadingEdge.CML_contribution"]][sel.dpost.c] )
      rownames(sig.nes) <- gsea.dpost.c[["pathway"]][sel.dpost.c]
      sig.df <- data.frame(sig.nes)
      
      state_rx.cols
      state.df <- data.frame("pathway" = rownames(sig.df), "contribution" = sig.df$TKI.vs.ctl, 
                             "state" = rep("D.Rx", dim(sig.df)[1]) )
      state.df$pathway <-  unlist(lapply(state.df$pathway, function(x) str_to_title(gsub("_"," ", gsub( "HALLMARK_", "", x)) )) )
      # line_positions holds divider line between each variable
      state.df2 <- state.df %>% mutate(line_positions = as.numeric(factor(pathway, levels = unique(pathway))), 
                                       line_positions = line_positions + .5,  
                                       line_positions = ifelse(line_positions == max(line_positions), NA, line_positions)) 
      png(paste(manFig_out,"/Fig-5E_fgsea-CML-contribution_p-",plimit,"_bar.png",sep=""), res=plot_res, units="in", height=3, width=7)
      p <- ggplot(state.df2, aes(x=reorder(pathway, -contribution), y=contribution, fill=state)) + 
        geom_bar(stat = "identity", position=position_dodge(width = 0.5), width=.5) + 
        theme_bw(base_size=12) + scale_fill_manual(values=state_rx.cols) +
        theme(axis.text.x = element_text(angle = 60,  hjust=1), legend.position="none", 
              panel.grid.major.x = element_blank(), panel.grid.minor.x = element_line(colour = "black") ) +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) 
      # geom_vline(aes(xintercept = line_positions))
      print(p)
      graphics.off()
    }
    
    # horizontal dot plot #
    {
      gsea.drx.c <- read.table("CML_contribution_tables/CML.contribution_D-Rx.vs.B-state-c3_GSEA_table.tsv", sep="\t", header=T)
      plimit <- 0.0001
      sel.drx.c <- which(gsea.drx.c$padj < plimit & abs(gsea.drx.c$NES) >=2 )
      gsea.drx.c$Pathway <- unlist(lapply(gsea.drx.c$pathway, function(x) str_to_title(gsub("_"," ", gsub( "HALLMARK_", "", gsub("_PATHWAY","",x))) )) )
      gsea.drx.c[sel.drx.c,]$padj
      #function to add space to margins
      margin_spacer <- function(x) {
        # where x is the column in your dataset
        left_length <- nchar(levels(factor(x)))[1]
        if (left_length > 8) {
          return((left_length - 8) * 4)
        }
        else
          return(0)
      }
      png(paste(manFig_out,"/Fig-5E_fgsea_p-",plimit,"_dotplot.png",sep=""), res=plot_res, units="in", height=3, width=7)
      ggplot(gsea.drx.c[sel.drx.c,], aes(y=NES, x=reorder(Pathway, padj), color=padj, size=size)) + geom_point() + 
        theme_bw(base_size=12) + scale_color_gradient(limits=c(0,.00001), low="red", high="blue") + 
        scale_x_discrete(labels = function(x) str_wrap(x, width = 25)) +
        theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 60,  hjust=1), 
              plot.margin = margin(l = 0 + margin_spacer(gsea.drx.c[sel.drx.c,]$Pathway)) ) +
        expand_limits(y = c(-2.5, 3.5))
      graphics.off()
      png(paste(manFig_out,"/Fig-5E_fgsea_p-",plimit,"_dotplot-legendOnly.png",sep=""), res=plot_res, units="in", height=6, width=6)
      ggplot(gsea.drx.c[sel.drx.c,], aes(x=NES, y=reorder(Pathway, padj), color=padj, size=size)) + geom_point() + 
        theme_bw(base_size=12) + scale_color_gradient(limits=c(0,.00001), low="red", high="blue") + 
        scale_x_discrete(labels = function(x) str_wrap(x, width = 25)) +
        theme(axis.text.x = element_text(angle = 60,  hjust=1), plot.margin = margin(l = 0 + margin_spacer(gsea.drx.c[sel.drx.c,]$Pathway)) )
      graphics.off()
    }
  }
  
  ###
  # Figure 5F: TOTO
  ###
  
  
  ###
  # Figure S7 - full upset
  ###
  {    
    cm.Rx = ComplexHeatmap::make_comb_mat(degs.Rx)
    comb.annot <- c()
    #n <- comb_name(cm.Rx)[1]
    for (ni in 1:length(comb_name(cm.Rx)) ){
      n <- comb_name(cm.Rx)[ni]
      n.list <- as.numeric(strsplit(n, "")[[1]])
      #name order: "c5", "c3", "c1", "TKI.Rx", "TKI-postRx", "TOTO-postRx"
      if (n.list[4]==1 & sum(n.list[4:6])==1) { #TKI RX only
        comb.annot <- c(comb.annot, "TKI-Rx")
      } else if (n.list[5]==1 & sum(n.list[4:6])==1) { #TKI postRx only
        comb.annot <- c(comb.annot, "TKI-postRx")
      } else if (n.list[6]==1 & sum(n.list[4:6])==1) { #TOTO postRx only
        comb.annot <- c(comb.annot, "TOTO-postRx")
      } else if  (n.list[5]==1 & n.list[4]==1  & sum(n.list[4:6])==2) { #any TKI
        comb.annot <- c(comb.annot, "TKI")
      } else if  (n.list[6]==1 & sum(n.list[4:6])>=2) { #mixed Rx
        comb.annot <- c(comb.annot, "mixed_Rx")
      } else {
        comb.annot <- c(comb.annot, "CML.vs.ctl")
      }
    }
    comb.cols <- c( "TOTO-postRx" = "#45e9f5", "TKI-postRx" = "#d955c9", "TKI"="#fca7cf", "TKI-Rx" = "#6507b3", "mixed_Rx"="#014508", "CML.vs.ctl"="grey")
    comb.sort <- sort.int(comb_size(cm.Rx), index.return = T)$ix
    png(paste(manFig_out,"/Fig-S7_vsControl+Rx_upset.png",sep=""), res=plot_res, units="in", height=3, width=6)
    UpSet(cm.Rx, set_order = names(degs.Rx), comb_order = rev(comb.sort),
          right_annotation = upset_right_annotation(cm.Rx,  add_numbers = TRUE,
                                                    gp = gpar(fill = c("red", "goldenrod", "dodgerblue", "#6507b3", "#d955c9",  "#45e9f5" )) ),
          top_annotation = HeatmapAnnotation(
            Treatment = comb.annot,
            col = list(Treatment=comb.cols),
            "Intersection\nsize" = anno_barplot(comb_size(cm.Rx), 
                                                border = FALSE, 
                                                gp = gpar(fill = "black"), 
                                                height = unit(2, "cm")
            ), 
            annotation_name_side = "left", 
            annotation_name_rot = 0) )
    graphics.off()
  }
  
  ###
  # Figure S7: Myeloid
  ###
  {
    b.rx.sel <- which(col.df$treatment=="TET_OFF_B" | col.df$state_rx=="A.post" | col.df$state_rx=="D.post" | col.df$state_rx=="D.Rx")
    png(paste(manFig_out,"/Fig-S5_myeloid.vs.CML_space_boxplot.png",sep=""), res=plot_res, units="in", height=3, width=5)
    p <- ggplot( col.df[b.rx.sel,], aes(x=state_rx, y=myeloid, fill=state_rx)) + geom_boxplot() +
      theme_bw(base_size=18) + scale_fill_manual(values=state_rx.cols)
    print(p)
    graphics.off()
  }
  
  ###
  # Fig S7A
  ###
  {
    col.df <- data.frame(colData(dat.se))
    png(paste(manFig_out,"/Fig-S7A_TOTO_myeloid.vs.CML_space_line+point.png",sep=""), res=plot_res, units="in", height=3, width=4)
    p <- ggplot( col.df[which(col.df$Group=="A"),], aes(y=myeloid, x=sample_weeks, color=Group, group=mouse_id)) + geom_line() + geom_point() +
      theme_bw(base_size=18) + theme(legend.position="none") + scale_color_manual(values=unlist(group.cols))
    print(p)
    graphics.off()
    png(paste(manFig_out,"/Fig-S7A_TKI_myeloid.vs.CML_space_line+point.png",sep=""), res=plot_res, units="in", height=3, width=4)
    p <- ggplot( col.df[which(col.df$Group=="D"),], aes(y=myeloid, x=sample_weeks, color=Group, group=mouse_id)) + geom_line() + geom_point() +
      theme_bw(base_size=18) + theme(legend.position="none") + scale_color_manual(values=unlist(group.cols))
    print(p)
    graphics.off()
    #BCR-ABL
    png(paste(manFig_out,"/Fig-S7A_TOTO_BCR-ABL.vs.CML_space_line+point.png",sep=""), res=plot_res, units="in", height=3, width=4)
    p <- ggplot( col.df[which(col.df$Group=="A"),], aes(y=HSA_BCR_ABL1_gene, x=sample_weeks, color=Group, group=mouse_id)) + geom_line() + geom_point() +
      theme_bw(base_size=18) + theme(legend.position="none") + scale_color_manual(values=unlist(group.cols))
    print(p)
    graphics.off()
    png(paste(manFig_out,"/Fig-S7A_TKI_BCR-ABL.vs.CML_space_line+point.png",sep=""), res=plot_res, units="in", height=3, width=4)
    p <- ggplot( col.df[which(col.df$Group=="D"),], aes(y=HSA_BCR_ABL1_gene, x=sample_weeks, color=Group, group=mouse_id)) + geom_line() + geom_point() +
      theme_bw(base_size=18) + theme(legend.position="none") + scale_color_manual(values=unlist(group.cols))
    print(p)
    graphics.off()
  }
  
  ###
  # Fig S7B
  ###
  {
    ### UPSET plots for up- and down-reg. separately ###
    # UP #
    apost.degs.up <- which(apost.vs.c$padj<0.05 & apost.vs.c$log2FoldChange>2)
    dpost.degs.up <- which(dpost.vs.c$padj<0.05 & dpost.vs.c$log2FoldChange>2)
    drx.degs.up <- which(drx.vs.c$padj<0.05 & drx.vs.c$log2FoldChange>2)
    degs.Rx.up <- list("c5.vs.ctl"=b2.vs.c[b2vc.degs.up, ]$gene_id, "c3.vs.ctl"=b1.vs.c[b1vc.degs.up, ]$gene_id, 
                       "c1.vs.ctl"=b0.vs.c[b0vc.degs.up, ]$gene_id, "TKI_Rx"=drx.vs.c[drx.degs.up,]$gene_id, 
                       "TKI_postRx"=dpost.vs.c[dpost.degs.up,]$gene_id, "TOTO_postRx"=apost.vs.c[apost.degs.up,]$gene_id)
    cm.Rx.up = ComplexHeatmap::make_comb_mat(degs.Rx.up)
    comb.annot <- c()
    for (ni in 1:length(comb_name(cm.Rx.up)) ){
      n <- comb_name(cm.Rx)[ni]
      n.list <- as.numeric(strsplit(n, "")[[1]])
      #name order: "c5", "c3", "c1", "TKI.Rx", "TKI-postRx", "TOTO-postRx"
      if (n.list[4]==1 & sum(n.list[4:6])==1) { #TKI RX only
        comb.annot <- c(comb.annot, "TKI-Rx")
      } else if (n.list[5]==1 & sum(n.list[4:6])==1) { #TKI postRx only
        comb.annot <- c(comb.annot, "TKI-postRx")
      } else if (n.list[6]==1 & sum(n.list[4:6])==1) { #TOTO postRx only
        comb.annot <- c(comb.annot, "TOTO-postRx")
      } else if  (n.list[5]==1 & n.list[4]==1  & sum(n.list[4:6])==2) { #any TKI
        comb.annot <- c(comb.annot, "TKI")
      } else if  (n.list[6]==1 & sum(n.list[4:6])>=2) { #mixed Rx
        comb.annot <- c(comb.annot, "mixed_Rx")
      } else {
        comb.annot <- c(comb.annot, "CML.vs.ctl")
      }
    }
    comb.cols <- c( "TOTO-postRx" = "#45e9f5", "TKI-postRx" = "#d955c9", "TKI"="#fca7cf", "TKI-Rx" = "#6507b3", "mixed_Rx"="#014508", "CML.vs.ctl"="grey")
    comb.sort <- sort.int(comb_size(cm.Rx.up), index.return = T)$ix
    png(paste(manFig_out,"/Fig-S7B_vsControl+Rx_UPREGULATED_upset.png",sep=""), res=plot_res, units="in", height=6, width=12)
    UpSet(cm.Rx.up, set_order = rev(names(degs.Rx.up)), comb_order = rev(comb.sort),
          right_annotation = upset_right_annotation(cm.Rx.up, 
                                                    gp = gpar(fill = c("red", "goldenrod", "dodgerblue", "#6507b3", "#d955c9",  "#45e9f5" )) ),
          top_annotation = HeatmapAnnotation(
            Comparisons = comb.annot,
            col = list(Comparisons=comb.cols),
            "Upregulated\nDEGs" = anno_barplot(comb_size(cm.Rx.up), 
                                               border = FALSE, 
                                               gp = gpar(fill = "#8c060a"), 
                                               height = unit(2, "cm")
            ), 
            annotation_name_side = "left", 
            annotation_name_rot = 0) )
    graphics.off()
    # DOWN #
    apost.degs.dn <- which(apost.vs.c$padj<0.05 & apost.vs.c$log2FoldChange < -2)
    dpost.degs.dn <- which(dpost.vs.c$padj<0.05 & dpost.vs.c$log2FoldChange < -2)
    drx.degs.dn <- which(drx.vs.c$padj<0.05 & drx.vs.c$log2FoldChange < -2)
    degs.Rx.dn <- list("c5.vs.ctl"=b2.vs.c[b2vc.degs.dn, ]$gene_id, "c3.vs.ctl"=b1.vs.c[b1vc.degs.dn, ]$gene_id, 
                       "c1.vs.ctl"=b0.vs.c[b0vc.degs.dn, ]$gene_id, "TKI_Rx"=drx.vs.c[drx.degs.dn,]$gene_id, 
                       "TKI_postRx"=dpost.vs.c[dpost.degs.dn,]$gene_id, "TOTO_postRx"=apost.vs.c[apost.degs.dn,]$gene_id)
    cm.Rx.dn = ComplexHeatmap::make_comb_mat(degs.Rx.dn)
    comb.annot <- c()
    for (ni in 1:length(comb_name(cm.Rx.dn)) ){
      n <- comb_name(cm.Rx)[ni]
      n.list <- as.numeric(strsplit(n, "")[[1]])
      #name order: "c5", "c3", "c1", "TKI.Rx", "TKI-postRx", "TOTO-postRx"
      if (n.list[4]==1 & sum(n.list[4:6])==1) { #TKI RX only
        comb.annot <- c(comb.annot, "TKI-Rx")
      } else if (n.list[5]==1 & sum(n.list[4:6])==1) { #TKI postRx only
        comb.annot <- c(comb.annot, "TKI-postRx")
      } else if (n.list[6]==1 & sum(n.list[4:6])==1) { #TOTO postRx only
        comb.annot <- c(comb.annot, "TOTO-postRx")
      } else if  (n.list[5]==1 & n.list[4]==1  & sum(n.list[4:6])==2) { #any TKI
        comb.annot <- c(comb.annot, "TKI")
      } else if  (n.list[6]==1 & sum(n.list[4:6])>=2) { #mixed Rx
        comb.annot <- c(comb.annot, "mixed_Rx")
      } else {
        comb.annot <- c(comb.annot, "CML.vs.ctl")
      }
    }
    comb.cols <- c( "TOTO-postRx" = "#45e9f5", "TKI-postRx" = "#d955c9", "TKI"="#fca7cf", "TKI-Rx" = "#6507b3", "mixed_Rx"="#014508", "CML.vs.ctl"="grey")
    comb.sort <- sort.int(comb_size(cm.Rx.dn), index.return = T)$ix
    png(paste(manFig_out,"/Fig-S7B_vsControl+Rx_DOWNREGULATED_upset.png",sep=""), res=plot_res, units="in", height=6, width=12)
    UpSet(cm.Rx.dn, set_order = rev(names(degs.Rx.dn)), comb_order = rev(comb.sort),
          right_annotation = upset_right_annotation(cm.Rx.dn, 
                                                    gp = gpar(fill = c("red", "goldenrod", "dodgerblue", "#6507b3", "#d955c9",  "#45e9f5" )) ),
          top_annotation = HeatmapAnnotation(
            Comparisons = comb.annot,
            col = list(Comparisons=comb.cols),
            "Downregulated\nDEGs" = anno_barplot(comb_size(cm.Rx.dn), 
                                                 border = FALSE, 
                                                 gp = gpar(fill = "#6d99c2"), 
                                                 height = unit(2, "cm")
            ), 
            annotation_name_side = "left", 
            annotation_name_rot = 0) )
    graphics.off()
    
  }
  
  ###
  # Fig S7C
  ###
  {
    ###
    ### build counts of DEGs vs state-space difference between compared groups
    ###
    {
      ss.v.deg <- c()
      ### CML vs Control ###
      {
        comps <- c("B-state-c3.vs.C", "B-state-c2.vs.C", "B-state-c1.vs.C", "B-state-c3.vs.B-state-c1", "B-state-c3.vs.B-state-c2", "B-state-c2.vs.B-state-c1",
                   "B-state-3_M.vs.F", "B-state-2_M.vs.F", "B-state-1_M.vs.F")
        #objects to use for looping
        trt1sel <- c("TET_OFF_B", "TET_OFF_B", "TET_OFF_B", "TET_OFF_B", "TET_OFF_B", "TET_OFF_B",
                     "TET_OFF_B", "TET_OFF_B", "TET_OFF_B" )
        trt2sel <- c("TET_ON_C", "TET_ON_C", "TET_ON_C", "TET_OFF_B", "TET_OFF_B", "TET_OFF_B", 
                     "TET_OFF_B", "TET_OFF_B", "TET_OFF_B")
        class1name <- c("state_chr", "state_chr", "state_chr", "state_chr", "state_chr", "state_chr",
                        "sex", "sex", "sex")
        class2name <- c(NA, NA, NA, "state_chr", "state_chr", "state_chr",
                        "sex", "sex", "sex")
        class1sel <- c("c3", "c2", "c1", "c3", "c3", "c2", "M", "M", "M")
        class2sel <- c(NA, NA, NA, "c1", "c2", "c1", "F","F", "F")
        
        for (ind in 1:6 ) { #skip sex comparisons
          compname <- comps[ind]
          trt1 <- trt1sel[ind]
          trt2 <- trt2sel[ind]
          class1 <- class1sel[ind]
          class2 <- class2sel[ind]
          cname1 <- class1name[ind]
          cname2 <- class2name[ind]
          if (trt1 == trt2) {
            lab1 <- paste(trt1, class1, sep="_")
            lab2 <- paste(trt2, class2, sep="_")
          } else {
            lab1 <- trt1
            lab2 <- trt2
          }
          print(paste("Processing: ",compname,sep=""))
          
          sel1 <- which(colData(dat.se)$treatment==trt1 & colData(dat.se)[[cname1]]==class1 )
          if (is.na(class2)) {
            sel2 <- which(colData(dat.se)$treatment==trt2 )
          } else {
            sel2 <- which(colData(dat.se)$treatment==trt2 & colData(dat.se)[[cname2]]==class2 )  
          }
          comp.se <- dat.se[, c(sel1, sel2) ]
          #add label to comp.se
          colData(comp.se)[["label"]] <- c( rep(lab1, length(sel1)), rep(lab2, length(sel2)))
          col.df <- data.frame(colData(comp.se))
          
          #add DEG diff to table
          ss1 <- mean(col.df$CML_space[which(col.df$label==lab1)])
          ss2 <- mean(col.df$CML_space[which(col.df$label==lab2)])
          cmean <- ss1 - ss2 
          
          #add DEG count to table
          dtab <-  read.table(paste(deg_out,"/DESeq-sexFactor_",compname,"_DEG-table-means.tsv",sep=""),sep="\t",header=T) 
          dcnt <- dim(dtab)[1]
          c.up <- length(which(dtab$log2FoldChange>0))
          c.down <- length(which(dtab$log2FoldChange<0))
          
          #add to table
          ss.v.deg <- rbind(ss.v.deg, c(compname, cmean, ss1, ss2, dcnt, c.up, c.down) )
        }
        
        
        
      }
      
      ### Treatment DEGs ###
      {
        comps <- c("A-postRx.vs.C", 
                   "D-postRx.vs.C", "D-Rx.vs.C", 
                   "A-postRx.vs.A-preRx-c1", "D-postRx.vs.D-preRx-c3", 
                   "D-Rx.vs.D-preT-c1", "D-Rx.vs.D-preT-c2", "D-Rx.vs.D-preT-c3",
                   "D-Rx.vs.A-postRx", "D-postRx.vs.A-postRx", "D-postRx.vs.D-preT-c2", "D-postRx.vs.D-preRx", "D-postRx.vs.B-c3",
                   "A-postRx.vs.B-state-c3", "A-postRx.vs.B-state-c2", "A-postRx.vs.B-state-c1",
                   "D-postRx.vs.B-state-c3", "D-postRx.vs.B-state-c2", "D-postRx.vs.B-state-c1")
        #objects to use for looping
        trt1sel <- c("TET_OFF_ON_A", 
                     "TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D",
                     "TET_OFF_ON_A", "TET_OFF_NIL_ON_D",
                     "TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D",
                     "TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D", "TET_OFF_NIL_ON_D", "TET_OFF_NIL_ON_D", "TET_OFF_NIL_ON_D",
                     "TET_OFF_ON_A", "TET_OFF_ON_A", "TET_OFF_ON_A", #vs CML
                     "TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D" #vs CML
        )
        trt2sel <- c("TET_ON_C", "TET_ON_C", "TET_ON_C", 
                     "TET_OFF_ON_A","TET_OFF_NIL_ON_D",
                     "TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D","TET_OFF_NIL_ON_D",
                     "TET_OFF_ON_A", "TET_OFF_ON_A", "TET_OFF_NIL_ON_D", "TET_OFF_NIL_ON_D", "TET_OFF_B", 
                     "TET_OFF_B", "TET_OFF_B", "TET_OFF_B", # A postRx vs B
                     "TET_OFF_B", "TET_OFF_B", "TET_OFF_B" # D postRx vs B
        )
        class1name <- c("rx_class", 
                        "rx_class", "rx_class", 
                        "rx_class", "rx_class", 
                        "rx_class","rx_class", "rx_class", 
                        "rx_class", "rx_class", "rx_class", "rx_class", "rx_class",
                        "rx_class","rx_class", "rx_class", #A postRX vs B
                        "rx_class","rx_class", "rx_class" #D postRx vs B
        )
        class2name <- c(NA,NA,NA, 
                        "state_chr", "state_chr", 
                        "state_chr", "state_chr", "state_chr", 
                        "rx_class", "rx_class", "state_chr", "rx_class", "state_chr", 
                        "state_chr", "state_chr", "state_chr", #A postRx vs B
                        "state_chr", "state_chr", "state_chr" #D postRx vs B
        )
        class1sel <- c("A.post", "D.post", "D.Rx", "A", "D.post", "D.Rx", "D.Rx", "D.Rx", "D.Rx", "D.post" , "D.post", "D.post", "D.post",
                       "A.post", "A.post", "A.post",
                       "D.post", "D.post", "D.post")
        class2sel <- c(NA,NA,NA, 
                       "c1","c3",
                       "c1","c2","c3", 
                       "A.post","A.post", "c2", "D.pre", "c3",
                       "c3","c2","c1", 
                       "c3","c2","c1" )
  
        for (ind in seq(1,length(comps))) {
          compname <- comps[ind]
          trt1 <- trt1sel[ind]
          trt2 <- trt2sel[ind]
          class1 <- class1sel[ind]
          class2 <- class2sel[ind]
          cname1 <- class1name[ind]
          cname2 <- class2name[ind]
          if (!is.na(class2)) {
            lab1 <- paste(trt1, class1,  sep="_")
            lab2 <- paste(trt2, class2,  sep="_")
          } else if (trt1 == trt2) {
            lab1 <- paste(trt1, class1, sep="_")
            lab2 <- paste(trt2, class2, sep="_")
          } else {
            lab1 <- trt1
            lab2 <- trt2
          }
          
          print(paste("Processing: ",compname,sep=""))
          
          ### condition 1; no exception for class == NA
          sel1 <- which(colData(dat.se)$treatment==trt1 & colData(dat.se)[[cname1]]==class1 )
          
          ### condition 2
          if (is.na(class2)  ) { #niether class2  have values
            sel2 <- which(colData(dat.se)$treatment==trt2 )
          } else {
            sel2 <- which(colData(dat.se)$treatment==trt2 & colData(dat.se)[[cname2]]==class2)
          }
          
          ###check conditions have matched some samples
          if (length(sel1)==0 | length(sel2)==0 ) {
            print(paste(":::ERROR::: unable to match samples to both conditions for comparison! SKIPPING!",sep=""))
            next
          }
          
          comp.se <- dat.se[, c(sel1, sel2) ]
          #add label to comp.se
          colData(comp.se)[["label"]] <- c( rep(lab1, length(sel1)), rep(lab2, length(sel2)))
          comp.df <- data.frame(colData(comp.se))
          #add DEG diff to table
          ss1 <- mean(comp.df$CML_space[which(comp.df$label==lab1)])
          ss2 <- mean(comp.df$CML_space[which(comp.df$label==lab2)])
          cmean <- ss1 - ss2 
          
          #add DEG count to table
          dtab <-  read.table(paste(deg_out,"/DESeq-sexFactor_",compname,"_DEG-table-means.tsv",sep=""),sep="\t",header=T) 
          dcnt <- dim(dtab)[1]
          c.up <- length(which(dtab$log2FoldChange>0))
          c.down <- length(which(dtab$log2FoldChange<0))
          
          #add to table
          ss.v.deg <- rbind(ss.v.deg, c(compname, cmean, ss1, ss2, dcnt, c.up, c.down) )
          
        }
      }
      
      #plot ss.diff vs DEGs
      dim(ss.v.deg)
      ssvd.df <- data.frame(ss.v.deg)
      colnames(ssvd.df) <- c("comp","ss.diff", "ss.1", "ss.2", "DEGs", "DEGs.up", "DEGs.down")
    }
    
    png(paste(manFig_out,"/Fig-S7C_DEGs.vs.state-space_DEGs-all.png",sep=""), res=plot_res, units="in", height=4, width=6)
    ggplot(ssvd.df, aes(x=as.numeric(DEGs), y=abs(as.numeric(ss.diff)))) + geom_point(size=3) + theme_bw(base_size=18) + labs(x="DEGs", y="State-space difference") 
    graphics.off()
  }
  
  ###
  # Fig S7D-E: CML contribution of treatment DEGs
  ###
  {
    gene_set_CML_contribution_manFigS7 <- function( out_name, deg_list, deg_comp, gene.df = gene.df ) { #all
      
      curdeg <- out_name 
      cgenes.id <- deg_list
      curcomp <- deg_comp
      
      sel <- match_genes( gene.df$gene_id, cgenes.id)
      c.df <- gene.df[sel,]
      
      m.1 <- c.df[[paste(curcomp,"_mean.",1,sep="")]]
      m.2 <- c.df[[paste(curcomp,"_mean.",2,sep="")]]
      c.diff <- m.1 - m.2
      c.ld <- c.df$eigengenes
      up.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]
      dn.ld <- c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]
      c.all <- c.diff * c.ld #get contribution for all genes in pathway
      c.up <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
      c.down <- c.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
      ctab <- matrix( c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", sum(c.ld), sum(up.ld), sum(dn.ld), sum(c.all,na.rm=T), sum(c.up), sum(c.down) ), nrow=3, ncol=5) 
      #add lfc contribution; need to figure out why "Inf" is produced
      c.lfc <- log2((m.1+amin)/(m.2+amin) )
      cl.all <- c.ld * c.lfc #get contribution for all genes in pathway based on LFC
      cl.up <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] #get up-regulated genes
      cl.down <- cl.all[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] #get down-regulated genes
      ctab <- matrix(c( rep(curdeg, 3), rep(curcomp, 3), "all", "up", "down", 
                        length(which(!is.na(c.all))), length(which(!is.na(c.up))), length(which(!is.na(c.down))), # count of genes contributing (non NA genes)
                        sum(c.ld), sum(up.ld), sum(dn.ld), #eigenvalue
                        sum(c.all,na.rm=T), sum(c.up, na.rm=T), sum(c.down, na.rm=T), #contribution w/ exp diff
                        sum(cl.all,na.rm=T), sum(cl.up, na.rm=T), sum(cl.down, na.rm=T) ), # contrib. w/ log2FC
                     nrow=3, ncol=7)
      if (length(which(is.na(ctab))) != 0) {
        print(paste(":::ERROR::: NAs produced for comp: ",curcomp," pathway: ",curdeg,sep=""))
        print(ctab)
        stop("NA encountered")
      }
      cont.tab <- rbind(cont.tab, ctab)
      #plot space
      c.df[["exp.diff"]] <- c.diff
      c.df[["cont.all"]] <- c.all
      dif <- c.df$exp.diff
      # dif[which(c.df$exp.diff>0)] <- "up"
      # dif[which(c.df$exp.diff<0)] <- "down"
      dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]>0)] <- "up"
      dif[which(c.df[[paste(curcomp,"_lfc",sep="")]]<0)] <- "down"
      c.df[["diff.lab"]] <- dif
      #build loading contribution
      "eigen.up" = mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      "eigen.down" <- mean(c.ld[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
      "PC1.up" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      "PC1.down" <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      
      ### make "eigen-contribution"
      # positive value = positive contribution (neg eig + pos exp OR pos eig + neg exp)
      # neg contribution (neg eig + neg exp OR pos eig + pos exp)
      cml.c <- c.df$eigengenes
      cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 #pro-CML up (neg eig. -> adjust)
      cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 #pro-CML down ( pos eig. -> NO adjust)
      cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] <- cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)] * 1 # anti-CML down (neg -> NO adjust)
      cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] <- cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)] * -1 # anti-CML up (pos -> adjust)
      
      #sums
      neg.up <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
      neg.up.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
      pos.down <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
      pos.down.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
      neg.down <- sum(cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0)], na.rm = F)
      neg.down.l <- length( cml.c[which(c.df$eigengenes < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0 & !is.na(c.df$eigengenes) )] )
      pos.up <- sum(cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0)], na.rm = F)
      pos.up.l <- length( cml.c[which(c.df$eigengenes > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0 & !is.na(c.df$eigengenes) )] )
      cmat <- matrix( c( rep(curdeg, 4), rep(curcomp, 4), "neg|up", "pos|down", "neg|down", "pos|up", 
                         neg.up, pos.down, neg.down, pos.up,
                         neg.up.l, pos.down.l, neg.down.l, pos.up.l),
                      nrow=4, ncol=5) 
      cml.tab <- rbind(cml.tab, cmat )
      #cml label
      cml.lab <- rep(NA, length(cml.c))
      cml.lab[which(cml.c > 0)] <- "pro-CML"
      cml.lab[which(cml.c < 0)] <- "anti-CML"
      #full matrix
      # :::note::: the mean contribution of each class of genes must be added!
      #(neg.up  / neg.up.l) + (pos.down / pos.down.l), (neg.down / neg.down.l) + (pos.up / pos.up.l)
      n.u.m <- (neg.up  / neg.up.l) 
      p.d.m <- (pos.down / pos.down.l)
      n.d.m <- (neg.down / neg.down.l) 
      p.u.m <- (pos.up / pos.up.l)
      fmat <- matrix( c( rep(curdeg, 2), rep(curcomp, 2), "pro-CML", "anti-CML",
                         sum(c(n.u.m, p.d.m), na.rm=T), sum(c(n.d.m, p.u.m), na.rm=T),
                         neg.up.l + pos.down.l, neg.down.l + pos.up.l),
                      nrow=2, ncol=5) 
      cml.full <- rbind(cml.full, fmat )
      ### eigengene value vs PC1 ###
      cl.up <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
      cl.down <- length(which(cml.c > 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
      al.up <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] > 0) )
      al.down <- length(which(cml.c < 0 & c.df[[paste(curcomp,"_lfc",sep="")]] < 0) )
      c.df[["CML_contribution"]] <- cml.c
      cml.up = mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      cml.down <- mean(cml.c[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)]) 
      PC1.up <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] > 0)]) 
      PC1.down <- mean(c.df$PC1[which(c.df[[paste(curcomp,"_lfc",sep="")]] < 0)])
      cml.tot <- mean(cml.c, na.rm=T) 
      PC1.tot <- mean(c.df$PC1, na.rm=T)
      #make df for labeling genes without NA genes
      l.df <- c.df[which(c.df$diff.lab=="up" | c.df$diff.lab=="down" ),]
      tup <- sort.int(l.df$CML_contribution, index.return = T, decreasing = T)$ix[seq(1,10)]
      tdown <- sort.int(l.df$CML_contribution, index.return = T, decreasing = F)$ix[seq(1,5)]
      png(paste(manFig_out,"/Fig-S7_eigenDEGs_",curdeg,"_comparison-",curcomp,"_pts+label.png",sep=""), res=plot_res, units="in", height=4, width=4)
      p <- ggplot(c.df, aes(x=CML_contribution, y= PC1, color=diff.lab, fill=diff.lab )) + 
        geom_point(aes(shape=diff.lab), alpha=.75) + theme_bw() + geom_vline(xintercept = 0, linetype="dashed", color="black", linewidth=1) +
        # geom_segment(aes(x = 0, y = 0, xend = cml.down, yend = PC1.down), color="#06398c", size=1.2, 
        #              arrow = arrow(length = unit(.25, "cm"), type="open"), lineend="butt", linejoin="mitre") +
        # geom_segment(aes(x = 0, y = 0, xend = cml.up, yend = PC1.up), color= "#ab1844", size=1.2,
        #              arrow = arrow(length = unit(.25, "cm"), type="open"), lineend="butt", linejoin="mitre") +
        geom_segment(aes(x = 0, y = 0, xend = cml.tot, yend = PC1.tot), color= "black", size=1.25,
                     arrow = arrow(length = unit(.30, "cm"), type="closed"), lineend="butt", linejoin="mitre") +
        scale_color_manual(values= list("up"="#f76a95", "down" = "#4287f5" )) + scale_fill_manual(values= list("up"="#f76a95", "down" = "#4287f5" )) +
        # geom_text(data=l.df[tup,], aes(x=CML_contribution, y=PC1+0.0007, color = diff.lab, label=gene), size=14) +
        # geom_text(data=l.df[tdown,], aes(x=CML_contribution, y=PC1+.0007, color = diff.lab, label=gene), size=14) +
        theme(legend.position = "none") + scale_shape_manual(values=c("up"=24, "down"=25)) +
        scale_x_continuous(limits=c(-.03, 0.058)) + scale_y_continuous(limits=c(-.025, 0.018)) 
      #ggtitle(paste("CML: Up=",cl.up,"; Down=",cl.down,"; anti-CML: Up=",al.up,"; Down=",al.down,sep="") )
      #scale_color_manual(values= list("up"="#ab1844", "down" = "#06398c" )) 
      print(p)
      graphics.off()

  
    } #end gene set contribution function
    ### TOTO minus age-related genes ###
    ig.toto.time <- setdiff(ig.apost.c, ig.c.lve)
    
    
    ### TOTO and TKI vs Ctrl ###
    # :note: same as Fig. 2D with different for loop indexes AND file name prefix in "gene_set_CML_contribution_manFig4" function
    
    
    deg_list <- list( early.deg.genes.id, trans.deg.genes.id, late.deg.genes.id, 
                      ig.b1.c, ig.b2.c, ig.b3.c,
                      ig.b3.b1, ig.b3.b2, ig.b2.b1,
                      ig.apost.c, ig.dpost.c, ig.drx.c, ig.c.evl, ig.c.lve, ig.apost.cL, ig.rx.da, ig.toto.time )
    deg_name <- list("Early",  "Transition", "Late", 
                     "B1.vs.C", "B2.vs.C", "B3.vs.C",
                     "B3.vs.B1", "B3.vs.B2", "B2.vs.B1",
                     "Apost.vs.C", "Dpost.vs.C", "Drx.vs.C", "C_early.vs.late", "C_late.vs.early", "Apost.vs.Clate", "Drx.vs.Apost", "TOTO.minus.time")
    deg_comp <- list( "B1.vs.C", "B2.vs.C", "B3.vs.C", 
                      "B1.vs.C", "B2.vs.C", "B3.vs.C",
                      "B3.vs.B1", "B3.vs.B2", "B2.vs.B1",
                      "Apost.vs.C", "Dpost.vs.C", "Drx.vs.C","C_early.vs.late", "C_late.vs.early", "Apost.vs.Clate", "Drx.vs.Apost", "Apost.vs.C")

    for (i in c(10, 14)){
      gene_set_CML_contribution_manFigS7( deg_name[i], deg_list[[i]], deg_comp[i], gene.df=gene.df)
    }
  }
  
  ###
  # Fig S7F
  ###
  {
    gsea.tki.post <- read.table("CML_contribution_tables/CML.contribution_D-postRx.vs.C_GSEA_table.tsv", sep="\t", header=T)
    gsea.tki.rx <- read.table("CML_contribution_tables/CML.contribution_D-Rx.vs.C_GSEA_table.tsv", sep="\t", header=T)
    
    plimit <- 0.0001
    sel.rx <- which(gsea.tki.rx$padj < plimit)
    sel.post <- which(gsea.tki.post$padj < plimit)
    
    
    # color by state-space #
    {
      all.sig <-  sort(unique(c(gsea.tki.rx[sel.rx,c("pathway")], gsea.tki.post[sel.post,c("pathway")] )))
      df.list <- list(gsea.tki.rx[sel.rx,], gsea.tki.post[sel.post,] )
      df.name <- c("tki.rx.vs.c", "tki.post.vs.c")
      sig.nes <- c()
      for (df in df.list) {
        out <- c()
        for (p in all.sig) {
          m <- which(df$pathway==p)
          if (length(m)==0) {
            out <- c(out,NA)
          } else {
            out <- c(out, df[["leadingEdge.CML_contribution"]][m])
          }
        }
        sig.nes <- cbind(sig.nes, out)
      }
      colnames(sig.nes) <- df.name
      rownames(sig.nes) <- all.sig
      sig.df <- data.frame(sig.nes)
      dim(sig.nes)
      tki.cols <- c("tki.rx.vs.c"="#6507b3", "tki.post.vs.c"="#d955c9")
      state.df <- data.frame("pathway" = rep(rownames(sig.df),2), "contribution" = c( sig.df$tki.rx.vs.c, sig.df$tki.post.vs.c ), 
                             "state" = c(rep("tki.rx.vs.c", dim(sig.df)[1]),rep("tki.post.vs.c", dim(sig.df)[1])  ) )
      state.df$pathway <-  unlist(lapply(state.df$pathway, function(x) str_to_title(gsub("_"," ", gsub( "HALLMARK_", "", x)) )) )
      state.df$state <- factor(state.df$state, levels=c("tki.post.vs.c", "tki.rx.vs.c"))
      # line_positions holds divider line between each variable
      state.df2 <- state.df %>% mutate(line_positions = as.numeric(factor(pathway, levels = unique(pathway))), 
                                       line_positions = line_positions + .5,  
                                       line_positions = ifelse(line_positions == max(line_positions), NA, line_positions)) 
      png(paste(manFig_out,"/Fig-S7_fgsea-CML-contribution_p-",plimit,"_bar.png",sep=""), res=plot_res, units="in", height=12, width=8)
      p <- ggplot(state.df2, aes(y=reorder(pathway, contribution), x=contribution, fill=state)) + 
        geom_bar(stat = "identity", position=position_dodge(width = 0.5), width=.5) + 
        theme_bw(base_size=20) + scale_fill_manual(values=tki.cols) +
        theme( legend.position="none", 
               panel.grid.major.y = element_blank(), panel.grid.minor.y = element_line(colour = "black") ) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
        geom_hline(aes(yintercept = line_positions))
      print(p)
      graphics.off()
    }
  }
  
  
  ###
  # Fig S7G
  ###
  {
    for (i in c(10)) {
      curcomp <- degcomp[[i]] 
      cc.df <- cont.df[which(cont.df$comp==curcomp & cont.df$DEGs=="all"),]
      ###find top 10 pathways - smallest (large negative) values have most contribution to leukem
      
      ### make contributions the AVERAGE so that large pathways don't skew ###
      cc.df[["cont.diff.sum"]] <- as.numeric(cc.df$contribution)
      cc.df[["cont.l2fc.sum"]] <- as.numeric(cc.df$cont.l2fc)
      cc.df[["contribution"]] <- cc.df$cont.diff.sum / as.numeric(cc.df$gene_counts)
      cc.df[["cont.l2fc"]] <- cc.df$cont.l2fc.sum / as.numeric(cc.df$gene_counts)
      
      #get sig pathways - NOTE: required objects
      cur.gsea <- gsea.objs[[i]]
      cur.sig <- which(cur.gsea$padj <= 0.05)
      psort <- sort.int(cur.gsea[cur.sig,]$padj, index.return = T, decreasing = F)$ix
      nes.sort <- sort.int(cur.gsea[cur.sig,]$NES, index.return = T, decreasing = T)$ix
      
      curpaths <- cur.gsea[cur.sig,]$pathway[psort]
      nes.paths <- cur.gsea[cur.sig,]$pathway[nes.sort]
      colnames(cur.gsea)
      cur.gsea[cur.sig,c(1,6)][nes.sort,]
      cur.gsea[,6]

      ### CML - eigengene contribution ###
      etpath <- curpaths #set so the existing code can be used
      ep.all <- c()
      ep.full <- c() #for matching all full pathways so that plot limits can be set
      for (p in etpath){
        ep.all <- c(ep.all, which(cml.df$pathway==p) )
        ep.full <- c(ep.full, which(full.df$pathway==p) )
      }
      
      if (length(curpath)==0) {next}
      
      ###plots
      ht <- length(curpaths) / 2.5
      # pro-CML ONLY
      
      c.df <- cml.df[intersect(ep.all, which(cml.df$comp==curcomp)), ]
      if (dim(c.df)[1]==0) {next}
      f.df <- full.df[intersect(ep.full, which(full.df$comp==curcomp)), ]
      evals <- c( as.numeric(f.df$CML.mean), 0) #add zero so that min is always at least zero
      #order pathways
      c.df$pathway <- factor(c.df$pathway, levels=rev(curpaths))
      f.df$pathway <- factor(f.df$pathway, levels=rev(curpaths))
      
      
      ###total contribution only
      png(paste(manFig_out,"/Fig-S7G_eigenpath_",curcomp,"_allSig_CML-contribution_total-Only.png",sep=""), res=plot_res, units="in", height=4, width=14)
      cdir <- rep("Pro-CML", dim(f.df)[1])
      cdir[which(as.numeric(f.df$CML.mean)<0)] <- "Anti-CML"
      f.df[["cont.dir"]] <- cdir
      p <- ggplot(f.df[which(f.df$dir=="total"),], aes(y=pathway, x=as.numeric(CML.mean), fill=cont.dir )) +  # multiply by -1 to get pro-CML to be positive
        geom_bar(stat="identity")  + theme_bw(base_size=18) +
        xlab("CML contribution") + ylab("") + 
        scale_fill_manual(values= list("Pro-CML"="#f54278", "Anti-CML" = "#4287f5" )) +
        coord_cartesian(xlim=c(min(evals)*1.1, max(evals)*1.1))
      print(p)
      graphics.off()
      
    }
  }
}


