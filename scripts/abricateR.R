#' A function for processing blast output.
#
#' Read the documents on the github page https://github.com/maxlcummins/abricateR for more info
#' @param file Filename of your abricate input
#' @param output Prefix for object/filename of your output files
#' @param identity Nucleotide identity cut-off. Default = 90
#' @param lenth Length/coverage cut-off. Default = 90
#' @param writecsv TRUE/FALSE: Write/Dont write to csv. Default = FALSE
#
#
#
#


abricateR <- function(file, output, identity = 90, length = 90, writecsv = FALSE) {
        
        require(readr)
        require(magrittr)
        require(dplyr)

        #provides a name for downstream files based on user input
        name <- output
        #reads input file (BLAST output in tab separated format.
        #                   Note: MUST BE GENERATED USING CUSTOM
        #                   PROVIDED BLAST SUBMISSION SCRIPT)
        #Colname reassignment
        
        
        #Read in full abricate genotype data sheet
        file <- read_delim(
                file,
                "\t",
                escape_double = FALSE,
                trim_ws = TRUE,
                col_names = TRUE
        )
        
        #Remove cases where there are multiple headers from concatenation of abricate reports
        file <- file %>% filter(SEQUENCE != "SEQUENCE")
        
        #Colname reassignment
        colnames(file)[c(1, 10:11)] <-
                c("name", "perc_coverage", "perc_identity")
        df_colnames <- colnames(file)
        
        #Convert percent coverage and identity to numeric type to allow filtering
        file$perc_coverage <- as.numeric(file$perc_coverage)
        file$perc_identity <- as.numeric(file$perc_identity)
        
        # appends a 1 to the final column of each row (used later by dcast)
        rep(x = 1, times = nrow(file)) -> file$gene_present
        
        # Removes .fasta* in column sample name
        gsub(pattern = "\\.fasta.*", replacement = "", x = file$name) ->
                file$name
        
        #adds the name of the database to the start of the gene name (for later sorting)
        file$GENE <- paste0(file$DATABASE, "_", file$GENE)
        
        # Assigns all hits meeting criteria of
        # X% Nucleotide match and Y% Length Match to variable filename.NX.LY.PASS
        # (tweak arguments 3 and 4 of blastlord to change the nucleotide and length criteria)
        file1 <- subset.data.frame(file, file$perc_identity >= identity & file$perc_coverage >= length)
        assign(paste(name, paste("N", identity, sep = ""), paste("L", length, sep =""), "PASS", sep = "."), file1, envir=globalenv())
        
        # Assigns all hits NOT meeting criteria of
        # X% Nucleotide match and Y% Length Match to variable filename.NX.LY.FAIL
        # (tweak arguments 3 and 4 of blastlord to change the nucleotide and length criteria)
        file2 <-  subset.data.frame(file, file$perc_identity < identity & file$perc_coverage < length)
        assign(paste(name, paste("N", identity, sep = ""), paste("L", length, sep =""), "FAIL", sep = "."), file2, envir=globalenv())
        
        require("reshape2")
        
        # Generates a simple summary table of genes (not including allele variants) from
        # file1 (PASS file) that contains hits that passed the set criteria
        dcast(data = file1, name ~ GENE,
              value.var = 'gene_present', drop = FALSE) -> file3
        assign(paste(name, "_simple_summary_", "N", identity, "L", length, sep = ""), file3, envir =globalenv())
        
        # Generates a table showing co-occurence of all genes on a given SEQUENCE for each sample.
        file5 <- file1 %>%
                group_by(name, SEQUENCE) %>%
                summarise(same_scaff=paste(unique(GENE), collapse=" "))
        assign(paste(name, "_co-occurence_", "N", identity, "L", length, sep= ""), file5, envir=globalenv())
        
        if(writecsv == TRUE) {
                message("Writing objects to csv file...")
                write.csv(file1, paste(paste(name, paste("N", identity, sep = ""), paste("L", length, sep =""), "PASS", sep = "."),".csv", sep = ""))
                write.csv(file2, paste(paste(name, paste("N", identity, sep = ""), paste("L", length, sep =""), "FAIL", sep = "."),".csv", sep = ""))
                write.csv(file3, paste(paste(name, "_simple_summary_", "N", identity, "L", length, sep = ""),".csv", sep = ""))
                write.csv(file5, paste(paste(name, "_same_scaff_", "N", identity, "L", length, sep= ""),".csv", sep = ""))
                message("Writing complete; script finished.")
        }
        else {
                message("Files not written to disk; script finished.")
        }
}
#***Changelog***
# blastlord4.0  - added an option to write PASS, FAIL and summary files to disk in CSV format. Default of write = FALSE
# blastlord5.0  - changed read.delim to read_delim from readr package.
#                 This greatly improves read time of the original CSV file - the most time consuming step
#               - added a change log
#