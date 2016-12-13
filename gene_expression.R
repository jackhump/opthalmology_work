library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
options(echo=T)
iFolder="/SAN/vyplab/HuRNASeq/opthalmology_work"

# rpkm data for each eye RNA dataset
MH_fibroblast_data <- "/SAN/vyplab/NCMD/JessicaGardner/RNAseq/fibroblasts/MH/deseq2/rpkm_values.csv"
BJ_JR_fibroblast_data <- "/SAN/vyplab/NCMD/JessicaGardner/RNAseq/fibroblasts/BJ_JR/deseq2/rpkm_values.csv"
BJ_JR_eyecup_data <- ""


annotation_file="/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/biomart_annotations_human.tab"
gtex_samples=paste(iFolder, "GTEx_Data_V6_Annotations_SampleAttributesDS.txt",sep="/")
gtex_subjects=paste(iFolder, "GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct", sep = "/")
# GTEx data 
gsample <- fread(gtex_samples, header=T )
gtex <- as.data.frame( fread(gtex_subjects, header=T ) )

# insertion region is chrX:136,420,790-144,420,790 
target <- c("chrX", 136420790, 144420790)
target <- c("chr3", 63499538, 63669292)

#biomart hg38
annotation <- fread(annotation_file)
# skin fibroblast RPKMs

# eyecup RPKMs

# read in each csv file
BJ_JR_fibro <- fread(BJ_JR_fibroblast_data, header=T)
MH_fibro <- fread(MH_fibroblast_data, header=T)

quit()
# match all the RPKMs together to make one big table
fibro <- merge( BJ_JR_fibro, MH_fibro, by = "ensemblID" )
# match in the coordinates for each gene. In complex cases split the EnsemblIDs and match the first EnsemblID
fibro_ids <- str_split_fixed(fibro$ensemblID, "\\+", 2)[,1] 
fibro$chr <- annotation$chromosome_name[match(fibro_ids, annotation$EnsemblID )]
fibro$start <- annotation$start_position[match(fibro_ids, annotation$EnsemblID )]
fibro$end <- annotation$end_position[match(fibro_ids, annotation$EnsemblID )]
# throw out genes not within the target region
fibro_region <- subset(fibro, chr == "chrX" & start >= 136420790 & end <= 144420790 & external_gene_id.x != "Y_RNA" & external_gene_id.x != "snoU13" & external_gene_id.x != "SNORA18" )

# how to visualise?
# remove genes that aren't expressed?
#fibro_rpkms <- select(fibro_region, BJ001,BJ002,BJ003,JR001,JR002,JR003,MH_1,MH_2,MH_3)
#fibro_region$mean_rpkm <- apply(fibro_rpkms, MAR=1, FUN = mean)
#fibro_region <- subset(fibro_region, mean_rpkm > 0.001)

# gather into a long dataframe
fibro_gather <- select(fibro_region,ensemblID,external_gene_id.x, BJ001,BJ002,BJ003,JR001,JR002,JR003,MH_1,MH_2,MH_3)
fibro_gather <-  as.data.frame( gather(fibro_gather, sample, rpkm, -ensemblID, -external_gene_id.x ) )  
names(fibro_gather)[2] <- "gene"
fibro_gather$subject <- substr(fibro_gather$sample, 1, 2)

# try a ggplot
# put the genes in coordinate order
gene.order <- fibro_region$external_gene_id.x
p <- ggplot( fibro_gather, aes(x = gene, y = rpkm, colour = subject) ) + 
	geom_jitter(size=0.8, width = 0.4) +
	ylab("RPKM") +
	ylim(c(0,50)) +
	scale_x_discrete(limits=gene.order) +
	theme( axis.text.x  = element_text(angle=45, hjust = 1, vjust=1),
	axis.title.x = element_blank(),
	legend.justification=c(1,1), # put legend within plotting area
	legend.position=c(1,1) ) +
	 labs(
    title = paste0("Fibroblast gene expression " ,target[1],":", target[2],"-", target[3] )
    caption = "expression presented in fragments per kilobase of exon per million mapped reads (FPKM)"
  )
 
plot <- paste0(iFolder,"fibroblast_expression_plot_",target[1],".pdf")
ggsave(plot, p, width = 60, height=15, units = "cm" )

# all the skin samples
gtex_skin_samples <- subset(gsample, gsample$SMTS == "Skin" )$SAMPID
# get just the transformed fibroblast samples
gtex_fibro_samples <- subset(gsample, gsample$SMTSD == "Cells - Transformed fibroblasts")$SAMPID

gtex_datasets <- c("all_skin" ,"fibroblasts_only")
gtex_datasets_samples <- list(gtex_skin_samples, gtex_fibro_samples)

for(i in 1:2){
	g_skin <- gtex[, which(names(gtex) %in% c("Name", "Description",gtex_datasets_samples[[i]]) ) ]

	g_skin$geneID <- str_split_fixed(g_skin$Name, "\\.", 2)[,1]

	# find just the genes covered by the interval
	gtex_skin_region <- g_skin[which(g_skin$geneID %in% fibro_region$ensemblID) , ]
	# gather 
	gtex_skin_region_gather <- gather(gtex_skin_region, sample, rpkm, -Name, -Description, -geneID)
	# plot
	gtex_expression_plot <- paste0(iFolder,"/gtex_", gtex_datasets[i],"_expression_plot_",target[1],".pdf")
	p <- ggplot( gtex_skin_region_gather, aes(x = Description, y = rpkm) ) + 
	geom_boxplot() +
	ylab("RPKM") +
	ylim(c(0,50) ) +
        scale_x_discrete(limits=gene.order) +
        theme( axis.text.x  = element_text(angle=45, hjust = 1, vjust=1),
        axis.title.x = element_blank() ) +
        labs(title = paste("GTEx gene expression ",target[1],":", target[2],"-", target[3]," ",gsub("_", " ", gtex_datasets[i]) ) ,
    		caption = "expression presented in fragments per kilobase of exon per million mapped reads (FPKM)")

	ggsave(gtex_expression_plot, p, width = 60, height=15, units = "cm" )
}


