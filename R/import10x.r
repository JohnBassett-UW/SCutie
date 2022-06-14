#' import10x
#' import 10x data from filtered feature matrices generated by Cell Ranger, and CITEseq count.
#' 	some public data sets are supported depending on file format and file naming convention.
#'
#' @param path Accepts either a directory containing compressed 'tsv' files or an 'h5' file.
#' @param type (optional) Force input file type to user input. Supported input types are "ff_bc_matrix" or "h5".
#' @param ESNG Defaults to FALSE. Logical that determines weather gene names should be used or ensembl gene id
#'
#' @importFrom Matrix Matrix readMM sparseMatrix
#' @importFrom hdf5r H5File
#' @importFrom utils read.table
#'
#' @return Sparse matrix of 10x single-cell count data. Class dgCMatrix.
#' @export
#'
#' @examples
#' path_to_file <- '/path/to/data'
#' data <- import10x(path_to_file)
#'
#' path_to_file <- '/path/to/data.h5'
#' data <- import10x(path = path_to_file, ESNG = T)
import10x <- function(path, type = "unspecified", ESNG = F, garbage_collection = T){
  supported_file_types <- c("ff_bc_matrix", "h5")
  ptm <- proc.time() #set time start

  ####Define sub routines#########################################################
  ##subRoutine##
  #extract_table_names#
  extract_table_names <- function(tables.list, check = F){
    tables.names <- names(tables.list)
    components <- c("barcodes", "features", "matrix")

    for(name in components){
      index <- grep(name, tables.names)
      names(tables.list)[index] <- name
    }

    matched = which(components %in% names(tables.list))

    if(length(matched) == 3){
      return(tables.list)
    }else if(length(matched) == 2){
      names(tables.list)[-matched] <- components[-matched]
      return(tables.list)
    }else{
      message("unable to match files to content")
      stop("file names must contain values 'barcodes', 'features', and 'matrix'")
    }
  }
  ##subRoutine##
  #compile_dgCMatrix#
  compile_dgCMatrix <- function(tables.list, attempt = 0) { #formatting function for dgCMatrices

    if(inherits(tables.list[["barcodes"]], "data.frame") == T){ #check if barcodes are data.frame and coerce to character vector
      tables.list[["barcodes"]] <- tables.list[["barcodes"]][,1]
    }

    data <- Matrix::Matrix(tables.list[["matrix"]] ,  # define dim names and ensure format
                           sparse = T ,
                           dimnames = list(tables.list[["features"]],
                                           tables.list[["barcodes"]]) ,
                           forceCheck = T,
                           doDiag = F)


    data <- as(data, Class = "dgCMatrix") # cast from 'T' to 'C' representation
  }

  ##subRoutine##
  #Check_for_file_exceptions#
  check_file_exceptions <- function(path, type) { #format and detect file type

    path <- file.path(path)

    if(!file.exists(path)){return(stop("file not found"))}

    if(!(type %in% supported_file_types)){  #check if type specified. if not, attempt to detect
      type <- strsplit(basename(path), split = "\\.")
      if(length(type[[1]]) == 1){ #check if has extension. if not assume directory.
        if(!dir_contains_ffbcmat(path)){stop("could not find files in directory")}
        type[[1]] <- path
        type[[2]] <- "ff_bc_matrix"
        return(type)
      }else if(type[[1]][length(type[[1]])] == 'h5'){ #check if extension is '.h5' and set type
        type[[2]] <- type[[1]][length(type[[1]])]
        type[[1]] <- path
        return(type)
      }else{ #any other file extension is unsupported
        try(type <- type[[1]][length(type[[1]])])
      }
    }else if(type %in% supported_file_types){ #check if specified type is accepted
      usr_def_type <- list()
      usr_def_type[[1]] <- path
      usr_def_type[[2]] <- type
      return(usr_def_type)
    }
    try(message((paste("specified file type '.", type, "' not recognized", sep = ''))))
    return(stop("invalid file or file type", call. = F))
  }

  ##subRoutine##
  #dir_contains_ffbcmat#
  dir_contains_ffbcmat <- function(path){      #check folder contents for subfiles of ff_bc_matrix

    cat("searching directory for associated files \n")
    dir.contents <- dir(path)
    extensions <- strsplit(dir(path), split = "\\.")
    components <- c("barcodes", "features", "matrix")
    recognized <- c(F,F,F)

    if(all(lengths(extensions) <= 1)){  #throw exception if files don't have extensions
      return(stop("path contains only directories or files without recognizeable extensions"))
    }

    for(i in 1:length(components)){
      recognized[i] <- length(grep(components[i], dir.contents)) > 0
    }

    cat(paste("location contains", + length(dir.contents), "files \n"))

    return(sum(recognized) >= 2)
  }

  ##SubRoutine#
  #generate_unique_names#
  #Example: (TBCE, TBCE, TBCE) -> (TBCE.1, TBCE.2, TBCE.3)
  generate_unique_names <- function(char_names, count = 0){ #accepts character vector of duplicate names and returns a vector of equal length containing unique elements
    count = count + 1
    dupes <- duplicated(char_names) # logical vector of duplicate indices
    dup_names <- char_names[dupes] # character vector of duplicate names

    if(length(dup_names) == 0){
      return(paste(char_names[!dupes], ".", count, sep = ""))
    }

    char_names[dupes] <- generate_unique_names(dup_names, count)
    char_names[!dupes] <- paste(char_names[!dupes], ".", count, sep = "") # set non-duplicate elements to enumerated name
    return(char_names)
  }

  ##subRoutine##
  #dimnames_unique#
  #checks if dimension names are unique. If UMIs are not unique, issues a warning. If gene symbols are not unique, replace duplicates with unique gene symbol names
  dimnames_unique <- function(data_10x){
    if(any(duplicated(dimnames(data_10x)[[2]]))){
      warning("duplicate UMI's found in data set")
    }

    dupes <- duplicated(dimnames(data_10x)[[1]]) #logical vector of duplicate indices

    if(any(dupes)){
      message("Warning in column names: gene symbols contain duplicate identities")
      dup_names <- dimnames(data_10x)[[1]][dupes] #character vector of duplicate names
      dimnames(data_10x)[[1]][dupes] <- generate_unique_names(dup_names)
      message("All column names have been made unique")
      return(data_10x)
    }else{
      cat("All column names verified unique \n")
      return(data_10x)
    }
  }


  ####allocate fun vars###########################################################
  tables.list <- list()
  type <- check_file_exceptions(path, type) #appropriately handle user input for type and path
  path <- type[[1]] # resets path to platform independent output from check_file_exceptions()


  cat(paste("format:", type[[2]], "\n")) #print detected file type
  ####main flow control switch for type specific importing########################
  data_10x <- switch(type[[2]],
                     "ff_bc_matrix" = { ############################################

                       #import files as elements of a list
                       cat("importing data from file: \n")
                       if(garbage_collection == T){gc()}
                       for(set in dir(path)){
                         cat(paste(set, "\n"))
                         table.index <- substr(set, 0, nchar(set)-7)
                         if(substr(set, nchar(set)-5, nchar(set)-3) == "tsv"){
                           tables.list[[table.index]] <- read.table(file.path(path, set))
                         }else{
                           message("Large matrices may take a minute to process...")
                           tables.list[[table.index]]<- Matrix::readMM(file = file.path(path, set))
                         }
                       }
                       if(garbage_collection == T){gc()}

                       tables.list <- extract_table_names(tables.list)

                       if(ESNG == T | ncol(tables.list[["features"]]) == 1){
                         tables.list[["features"]] <- tables.list[["features"]][,1] #10x Output w/ ESNG option OR CITEseq count hash matrix
                       }else if(ncol(tables.list[["features"]]) >=2){ #10x output
                         tables.list[["features"]] <- tables.list[["features"]][,2]
                       }
                       compile_dgCMatrix(tables.list) #compile imported files into a single dgCMatrix by index
                     },
                     "h5" = { ######################################################
                       #import data from h5
                       activeFile <- hdf5r::H5File$new(filename =  path, mode = 'r')
                       activeFile.datasets <- names(activeFile)
                       #support for v3 (no multimodal support)
                       if(length(activeFile.datasets) > 1){stop("h5 files with multiple datasets not supported")}

                       data <- activeFile[[activeFile.datasets]]
                       data.attrs <- names(data)
                       if(ESNG == T){data.attrs[data.attrs == "features"] <- "features/id"}else{
                         data.attrs[data.attrs == "features"]<- "features/name"
                       }
                       cat("importing data from h5: \n")
                       for(set in data.attrs){
                         cat(paste(set, "\n"))
                         tables.list[[set]] <- data[[set]][]
                       }

                       names(tables.list)[grep("features", names(tables.list))] <- "features"

                       tables.list[["matrix"]] <- Matrix::sparseMatrix(
                         i = tables.list[["indices"]] + 1,
                         p = tables.list[["indptr"]],
                         x = tables.list[["data"]],
                         dims = tables.list[["shape"]],
                         repr = "T")

                       activeFile$close_all()

                       compile_dgCMatrix(tables.list)
                     }
  )

  cat("[*_*] Done. \n time ")
  print((proc.time() - ptm)[3]) #print time elapsed


  return(dimnames_unique(data_10x))
}
