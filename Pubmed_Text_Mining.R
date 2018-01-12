if(!"RISmed" %in% rownames(installed.packages()))  install.packages("RISmed")
if(!"stringr" %in% rownames(installed.packages()))  install.packages("stringr")

library(RISmed)
library(stringr)
#data <- read.table("TFs.txt", header=F, stringsAsFactors = F)
data <- read.table("GeneSet_1.txt", header=F, stringsAsFactors = F)
gene_info <- read.delim("Mus_musculus.gene_info", header=FALSE, stringsAsFactors=FALSE)
gene_info$symbols[gene_info$V5!='-'] <- paste("", gene_info$V3[gene_info$V5!='-'], gene_info$V5[gene_info$V5!='-'], "", sep="|")
gene_info$symbols[gene_info$V5=='-'] <- paste("", gene_info$V3[gene_info$V5=='-'], "", sep="|")

alias <- vector()
count <- vector()
PMID <- vector()
query_trans <- vector()

for(i in 1:length(data$V1)){
    if(length(gene_info$symbols[grep(paste("^\\|",data$V1[i],"\\|",sep=""),gene_info$symbols, perl =T)]))
        tmp <- str_sub(gene_info$symbols[grep(paste("^\\|",data$V1[i],"\\|",sep=""),gene_info$symbols, perl =T)], 2, -2)
    else{
        if(length(gene_info$symbols[grep(paste("\\|",data$V1[i],"\\|",sep=""),gene_info$symbols, perl =T)]))
            tmp <- str_sub(gene_info$symbols[grep(paste("\\|",data$V1[i],"\\|",sep=""),gene_info$symbols, perl =T)], 2, -2)
        else
            tmp <- data$V1[i]
    }
    alias <- c(alias, tmp)
#     if(gene_info$V5[gene_info$V3 == data$V1[i]] == '-'){
#         query_term <- paste(data$V1[i], " AND ", "Tregs", " AND ", "(cancer or tumor)", sep = "")
#         query_TF <- data$V1[i]
#     }
#     else{
#         query_term <- paste("(", paste(data$V1[i],gene_info$V5[gene_info$V3 == data$V1[i]], sep="|"), ") AND ", "Tregs", " AND ", "(cancer or tumor)", sep = "")
#         query_TF <- paste(data$V1[i],gene_info$V5[gene_info$V3 == data$V1[i]], sep="|")
#     }
    #query_term <- paste(data$V1[i],"[Title/Abstract]"," AND ", "Tregs","[Title/Abstract]",sep = "")

    query_term <- paste("(", tmp, ") AND ", "T lymphocytes", sep = "")
    #query_term <- paste("(", tmp, ") AND ", "Tregs", " AND ", "T lymphocytes", sep = "")
    query_TF <- tmp

    res <- EUtilsSummary(query_term, type="esearch", db="pubmed")
    res_TF <- EUtilsSummary(query_TF, type="esearch", db="pubmed")

    if(QueryCount(res_TF)){
        count <- c(count, QueryCount(res))
        if(QueryCount(res))
            PMID <- c(PMID, paste(QueryId(res), collapse = ";"))
        else
            PMID <- c(PMID, "NA")
        query_trans <- c(query_trans, QueryTranslation(res))
    }
    else{
        count <- c(count, 0)
        PMID <- c(PMID, "NA")
        query_trans <- c(query_trans, paste(data$V1[i], "not found in PubMed"))
    }
}

results <- data.frame("Gene_name"=data$V1, "Aliases"=alias, "Query"=query_trans, "Count"=count, "PMID"=PMID)

write.table(results, "Pubmed_query_hits_alias_GeneSet_1_T_lymphocytes.csv", row.names=F, col.names=T, sep=",", quote=F)
