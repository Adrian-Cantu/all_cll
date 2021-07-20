library(gt23)
library(RMySQL)
library(dplyr)
library(geneRxCluster)
library(tidyr)
#library(ggplot2)
#source('lib.R')
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Homo.sapiens)
all_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp where Trial="CART19_CLL"')


if(! file.exists('all_cll_intSites.rds')){
intSites <- getDBgenomicFragments(samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq') %>%
  stdIntSiteFragments() %>%
  collapseReplicatesCalcAbunds() %>%
  annotateIntSites()
  intSites <- data.frame(intSites) %>% filter(seqnames %in% paste0('chr', c(1:22, 'X'))) %>%
    mutate(cellType=ifelse(cellType=='Whole blood','Whole Blood',cellType))
  saveRDS(intSites,file='all_cll_intSites.rds')
} else {
  intSites <- readRDS('all_cll_intSites.rds')
}

# make a table summary of the samples
kk2 <- group_by(data.frame(intSites), cellType,timePointMonths,GTSP) %>%
  summarise(cellsPerSample = sum(estAbund)) %>%
  arrange(cellType,timePointMonths)

group_by(data.frame(intSites),cellType,timePointMonths) %>%
  summarise(cellsPerSample = sum(estAbund)) %>% arrange(timePointMonths) %>%
  pivot_wider(names_from = cellType, values_from = cellsPerSample,values_fill=0)

group_by(data.frame(intSites),cellType,timePointMonths,GTSP) %>%
  summarise() %>% group_by(cellType,timePointMonths) %>%
  summarise(Sample_num = n()) %>%
  arrange(timePointMonths)


# count samples before and after removing samples with less than 100 cells
group_by(data.frame(intSites), GTSP) %>% mutate(cellsPerSample = sum(estAbund)) %>% summarise() %>% nrow(.)
group_by(data.frame(intSites), GTSP) %>% mutate(cellsPerSample = sum(estAbund)) %>% filter(cellsPerSample >= 100) %>% summarise() %>% nrow(.)

#filter our samples with less than 100 cells  
filtered_intSites <- group_by(data.frame(intSites), GTSP) %>% mutate(cellsPerSample = sum(estAbund)) %>% filter(cellsPerSample >= 100) %>% ungroup()

#number of insertion sites per group
filtered_intSites %>% filter(timePointMonths == 0 & cellType=='Whole Blood')  %>% nrow(.)     #%>% makeGRangesFromDataFrame(.)
filtered_intSites %>% filter(timePointMonths >= 1 & cellType=='Whole Blood' )  %>% nrow(.)    # %>% makeGRangesFromDataFrame(.)

g1 <-filtered_intSites %>% filter(timePointMonths == 0 & cellType=='Whole Blood')  %>% makeGRangesFromDataFrame(.,keep.extra.columns=TRUE)
g2 <-filtered_intSites %>% filter(timePointMonths >= 1 & cellType=='Whole Blood' ) %>% makeGRangesFromDataFrame(.,keep.extra.columns=TRUE)

#data.frame(g1) %>% nrow(.)
#data.frame(g2) %>% nrow(.)

list_genes <- function(chr,g_start,g_end) {
  temp <- intSites %>% filter(seqnames==chr & start>=g_start & start<=g_end)
  return(paste(unique(temp$nearestFeature),collapse = ' '))
}

list_all_genes <- function(seqnames,start,end,genes) {
  tmp_R <-GRanges(paste0(seqnames,':',start,'-',end,':*'))
  tmp <- subsetByOverlaps(genes, tmp_R)
  tmp_g <-  as.data.frame(org.Hs.egSYMBOL) %>% filter(gene_id %in% tmp$gene_id)
  #return(list(tmp$gene_id))
  #return(paste(tmp$gene_id,collapse = ' '))
  return(paste(tmp_g$symbol,collapse = ' '))
 }

results <- scanStats(g1,g2,gr1.label = "T0",gr2.label = "T1", kvals = "3L:20L") %>% GenomicRanges::as.data.frame(.) %>% arrange(clusterSource)


#names(subsetByOverlaps(all_genes, results))
#as.data.frame(org.Hs.egSYMBOL) %>% head
#makeGRangesFromDataFrame


p_results <- results %>% filter(width < 1500000)
p_result_anot <- p_results %>% rowwise() %>% mutate(genes_intsites=list_genes(as.character(seqnames),start,end),
                                   genes_entrez=list_all_genes(seqnames = seqnames,start = start,end = end,all_genes))
openxlsx::write.xlsx(p_result_anot, file = 'scan_table.xlsx')


stat3 <- as.data.frame(org.Hs.egSYMBOL) %>% filter(symbol=='STAT3')
stat3_r <- all_genes["6774"]
subsetByOverlaps(g1, stat3_r)
subsetByOverlaps(g2, stat3_r)

