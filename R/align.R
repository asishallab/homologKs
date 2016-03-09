#' Runs the pipeline to align a set of coding sequences: First translates them,
#' then validates them for premature stop codons, subsequently generates a
#' multiple sequence alignment (MSA) of amino acid (AA) sequences, then uses
#' this AA MSA as guide and aligns the coding sequences in the final step.
#'
#' @param cds an instance of Biostrings::DNAStringSet representing the coding
#' sequences that need to be aligned
#' @param aas instance of Biostrings::AAStringSet containing the amino acid
#' sequences for the genes held in 'cds'. You can use Biostrings::transate to
#' generate them, or the Java program MACSE.
#' @param work.dir the working directory to use and in which to save the
#' relevant files
#' @param gene.group.name a string being used to name the output files written
#' into work.dir. Could be something like 'fam1234'.
#' @param mafft.call The string passed to system to invoke the alignment
#' program MAFFT on the translated amino acid sequences. Defaut is 'mafft
#' --auto'. Set option 'paranomeKsR.mafft.call' to change this defaut.
#' @param pal2nal.call The string processed and then passed to system to invoke
#' the program 'pal2nal.pl' in order to generate a multiple coding sequence
#' alignment. Three substrings must be present and will be replaced with their
#' appropriate paths: #msa#, #cds#, #cds.msa.path#. Default is '#msa# #cds#
#' -nogap -nomismatch -output paml > #cds.msa.path#', set option
#' 'paranomeKsR.pal2nal.call' to change this defaut.
#'
#' @import Biostrings
#' @export
#' @return The path to the aligned coding sequences file
alignCodingSequencesPipeline <- function(cds, aas, work.dir, gene.group.name, mafft.call = getOption("paranomeKsR.mafft.call", 
    "mafft --auto"), pal2nal.call = getOption("paranomeKsR.pal2nal.call", paste(file.path(path.package("paranomeKsR"), 
    "pal2nal.pl"), "#msa# #cds# -nogap -nomismatch -output paml > #cds.msa.path#"))) {
    cds.path <- file.path(work.dir, paste(gene.group.name, "_CDS.fasta", sep = ""))
    writeXStringSet(cds, cds.path)
    aas.path <- sub("_CDS", "_AAS", cds.path, fixed = TRUE)
    writeXStringSet(aas, aas.path)
    cds.msa.path <- file.path(work.dir, paste(gene.group.name, "_CDS_MSA.fasta", 
        sep = ""))
    aas.msa.path <- sub("_CDS_", "_AAS_", cds.msa.path, fixed = TRUE)
    system(paste(mafft.call, aas.path, ">", aas.msa.path))
    system(sub("#cds.msa.path#", cds.msa.path, sub("#cds#", cds.path, sub("#msa#", 
        aas.msa.path, pal2nal.call, fixed = TRUE), fixed = TRUE), fixed = TRUE))
    cds.msa.path
} 
