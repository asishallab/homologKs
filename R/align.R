#' Uses the program MACSE to translate the coding sequences to amino acid
#' sequences.
#'
#' @param path.2.cds.fasta gives the valid file path to the FASTA file holding
#' the coding sequences.
#' @param macse.call The call with which to invoke MACSE through the
#' system(...) function. Use option 'macse.call' to modify the default.
#'
#' @return TRUE if and only if no error has occurred.
#' @export
translate2AASeqs <- function(path.2.cds.fasta, macse.call = getOption("macse.call", 
    paste("java -Xmx600m -jar ", file.path(path.package("paranomeKsR"), "macse_v1.01b.jar"), 
        " -prog translateNT2AA", sep = ""))) {
    cmd <- paste(macse.call, "-seq", path.2.cds.fasta)
    system(cmd)
    TRUE
}

#' Validates an AAString
#'
#' @param 'aa.seq' an instance of Biostrings::AAString
#'
#' @return TRUE if and only if the argument does not have a premature stop codon.
#' @export
validateAASeqs <- function(aa.seq) {
    if (class(aa.seq) != "AAString") 
        stop("Argument 'aa.seq' is not of class 'AAString'.")
    !grepl("\\S\\*\\S", aa.seq, perl = TRUE)
}

#' Validates an AAStringSet
#'
#' @param 'aa.seq' an instance of Biostrings::AAStringSet
#'
#' @return A logical vector indicating the indices of valid AA-Seqs in the argument set.
#' @export
validateAAStringSet <- function(aa.set) {
    as.logical(lapply(names(aa.set), function(acc) {
        tryCatch({
            validateAASeqs(aa.set[[acc]])
        }, error = function(e) {
            warning("Amino-Acid-Sequence ", acc, " caused an error: ", e)
            FALSE
        })
    }))
}

#' Aligns the set of amino acid sequences based on the chemical properties of
#' their respective amino acid residues. Does so using MAFFT. NOTE: mafft is
#' expected to be installed.
#'
#' @param 'path.2.aa.seqs' the valid file path to the amino acid sequences
#' FASTA
#' @param 'path.2.aa.msa' the valid file path to the amino acid multiple
#' sequence alignment file that is to be generated
#' @param 'mafft.call' the string used to build the system command which
#' invokes MAFFT.
#'
#' @return TRUE if and only of no error occurred.
#' @export
alignAAStringSet <- function(path.2.aa.seqs, path.2.aa.msa, mafft.call = getOption("mafft.call", 
    paste("mafft --thread", getOption("mc.cores", 1), "--auto"))) {
    cmd <- paste(mafft.call, path.2.aa.seqs, ">", path.2.aa.msa)
    system(cmd)
    TRUE
}

#' Uses an amino acid multiple sequence alignment to guide the alignment of
#' codon sequence 'cds'.
#'
#' @param 'cds' an instance of Biostrings::DNAString representing the codon
#' sequence to align.
#' @param 'aligned.aa.seq' an instance of Biostrings::AAString representing the
#' aligned amino acid sequence.
#' @param 'aa.gap' the string used to represent a gap in the aligned AA sequence
#' @param 'codon.gap' the string used to represent a gap in the aligned codon
#' sequence
#'
#' @return An instance of Biostrings::DNAString representing the aligned codon
#' sequence
#' @export
alignCDSWithAlignedAASeq <- function(cds, aligned.aa.seq, aa.gap = "-", codon.gap = "---") {
    aa.chars <- strsplit(toString(aligned.aa.seq), split = NULL)[[1]]
    cds.chars <- strsplit(toString(cds), split = NULL)[[1]]
    i <- 1
    DNAString(paste(lapply(aa.chars, function(aa.char) {
        if (aa.char == aa.gap) {
            codon.gap
        } else {
            start.ind <- (i - 1) * 3 + 1
            stop.ind <- i * 3
            i <<- i + 1
            paste(cds.chars[start.ind:stop.ind], collapse = "")
        }
    }), collapse = ""))
}

#' Uses an amino acid multiple sequence alignment to guide the alignment of
#' codon sequence 'cds'.
#'
#' @param 'unaligned.cds.set' an instance of Biostrings::DNAStringSet
#' representing the unaligned codon sequences
#' @param 'aligned.aa.set' an instance of Biostrings::AAMultipleAlignment
#' representing the aligned amino acid sequences used as guide
#'
#' @return An instance of Biostrings::DNAMultipleAlignment representing the
#' aligned codon sequences.
#' @export
alignCDSSetWithAlignedAAsAsGuide <- function(unaligned.cds.set, aligned.aa.set) {
    cds.set <- DNAStringSet(lapply(names(aligned.aa.set), function(acc) {
        alignCDSWithAlignedAASeq(unaligned.cds.set[[acc]], aligned.aa.set[[acc]])
    }))
    names(cds.set) <- names(aligned.aa.set)
    DNAMultipleAlignment(cds.set)
}

#' Many programs in Bioinformatics cannot handle gene accessions with special
#' characters. That is why we rename the genes of an Biostrings::XStringSet
#' with 'PROT1', 'PROT2'...
#'
#' @param 'xstring.set' an instance of Biostrings::XStringSet
#'
#' @return An instance of data.frame with two columns the 'original' and the
#' 'sanitized' names.
#' @export
sanitizeNames <- function(xstring.set) {
    df <- data.frame(original = names(xstring.set), stringsAsFactors = FALSE)
    df$sanitized <- paste("PROT", 1:nrow(df), sep = "")
    df
}

#' Returns the matching original names for the sanitized ones held in argument
#' 'san.names'.
#'
#' @param 'san.names' a character vector of sanitized names (gene accessions)
#' @param 'name.mappings' a data.frame generated by function
#' 'sanitizeNames(...)'. This data.frame is used to lookup the matching
#' original names for the argument sanitized names.
#'
#' @return A character vector in which the original names are stored in the
#' order of their matching sanitized ones as they appear in the argument
#' 'san.names'.
#' @export
replaceWithOriginal <- function(san.names, name.mappings) {
    as.character(lapply(san.names, function(x) {
        name.mappings[which(name.mappings$sanitized == x), "original"]
    }))
}

#' Runs the pipeline to align a set of coding sequences: First translates them,
#' then validates them for premature stop codons, subsequently generates a
#' multiple sequence alignment (MSA) of amino acid (AA) sequences, then uses
#' this AA MSA as guide and aligns the coding sequences in the final step.
#'
#' @param cds an instance of Biostrings::DNAStringSet representing the coding
#' sequences that need to be aligned
#' @param work.dir the working directory to use and in which to save the
#' relevant files
#' @param gene.group.name a string being used to name the output files written
#' into work.dir. Could be something like 'fam1234'.
#'
#' @return The ALIGNED and validated coding sequences as an instance of
#' Biostrings::DNAStringSet, or nothing if validation discards the rest of
#' 'cds'.
#' @export
alignCodingSequencesPipeline <- function(cds, work.dir, gene.group.name) {
    #' Sanitize the gene identifiers:
    name.maps <- sanitizeNames(cds)
    names(cds) <- name.maps$sanitized
    cds.path <- file.path(work.dir, paste(gene.group.name, "_CDS.fasta", sep = ""))
    writeXStringSet(cds, cds.path)
    name.maps.path <- file.path(work.dir, paste(gene.group.name, "_name_mappings_table.txt", 
        sep = ""))
    write.table(name.maps, name.maps.path, row.names = FALSE, sep = "\t", quote = FALSE)
    #' Convert to AA and align the AA-sequences:
    translate2AASeqs(cds.path)
    #' Remove invalid AA-Sequences, i.e. AA-Seqs with premature stop-codons:
    aa.path <- file.path(work.dir, paste(gene.group.name, "_CDS_macse_AA.fasta", 
        sep = ""))
    aas <- readAAStringSet(aa.path)
    aas.san <- aas[validateAAStringSet(aas)]
    #' Warn about removed AA-Seqs:
    cds.san <- if (length(aas.san) < length(aas)) {
        len.diff <- length(aas) - length(aas.san)
        warning(len.diff, " amino-acid-sequences had a premature stop codon and were removed from further analysis.")
        cds[names(aas.san)]
    } else cds
    #' If only a single sequence is left, we're done:
    cds.msa.orig <- if (length(cds.san) > 1) {
        #' Write out the sanitized amino acid seqs:
        aas.san.path <- file.path(work.dir, paste(gene.group.name, "_AA_sanitized.fasta", 
            sep = ""))
        writeXStringSet(aas.san, aas.san.path)
        #' Generate a multiple sequence alignment:
        aas.msa.path <- file.path(work.dir, paste(gene.group.name, "_AA_sanitized_MSA.fasta", 
            sep = ""))
        alignAAStringSet(aas.san.path, aas.msa.path)
        aas.san.msa <- readAAMultipleAlignment(aas.msa.path)
        #' Use the aligned AA-Seqs as quide to align the CDS Sequences:
        cds.msa <- alignCDSSetWithAlignedAAsAsGuide(cds.san, attr(aas.san.msa, "unmasked"))
        cds.msa.path <- file.path(work.dir, paste(gene.group.name, "_CDS_MSA.fasta", 
            sep = ""))
        writeXStringSet(attr(cds.msa, "unmasked"), cds.msa.path)
        cds.msa@unmasked
    } else {
        cds.san
    }
    #' Return the CDS MSA as an instance of Biostrings::DNAStringSet using
    #' the ORIGINAL gene identifiers:
    if (length(cds.msa.orig) > 0) 
        names(cds.msa.orig) <- replaceWithOriginal(names(cds.msa.orig), name.maps)
    cds.msa.orig
}

#' Filters a multiple coding sequence alignment discarding positions where any
#' sequence has a gap or a 'N' character.
#'
#' @param cds.msa an instance of Biostrings::XStringSet representing the
#' multiple sequence alignment
#' @param chars.2.discard.regex a string representing the regular expression
#' defining the character class of the positions to be discarded. Default is
#' '[nN-]'
#' @param seq.set.type one of either Biostrings::DNAStringSet or
#' Biostrings::AAStringSet. Default is Biostrings::DNAStringSet.
#'
#' @return The filtered XStringSet with the possibly reduced sequences.
#' @export
filterMSA <- function(cds.msa, chars.2.discard.regex = "[nN-]", seq.set.type = DNAStringSet) {
    disc.pos <- Reduce(union, lapply(cds.msa, function(x) {
        as.integer(gregexpr("[nN-]", toString(x))[[1]])
    }))
    seq.set.type(setNames(lapply(cds.msa, function(x) {
        x[setdiff(1:length(x), disc.pos)]
    }), names(cds.msa)))
}

#' Tests for function filterMSA
#'
#' @return TRUE if and only if the tests do not fail
#' @export
testFilterMSA <- function() {
    nms <- c("Prot_A", "Prot_B")
    t.msa <- DNAStringSet(setNames(list(DNAString("ACG-TNT"), DNAString("TN-TACG")), 
        nms))
    exp.msa <- DNAStringSet(setNames(list(DNAString("ATT"), DNAString("TAG")), nms))
    res.msa <- filterMSA(t.msa)
    checkEquals(toString(exp.msa[[1]]), toString(res.msa[[1]]))
    checkEquals(toString(exp.msa[[2]]), toString(res.msa[[2]]))
    checkEquals(names(exp.msa), names(t.msa))
} 
