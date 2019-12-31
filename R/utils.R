# Annoying ERROR: hard-coded installation path
# .path_config_file <- file.path(
#   system.file("extdata",
#     package = "neopeptides"
#   ),
#   "path_config.yml"
# )
.path_config_file <- getOption(
  "neopeptides.config",
  default = file.path(
    Sys.getenv("HOME", unset = "~/"),
    ".neopeptide",
    "config.yml"
  )
)

#' Find Entry Path from Environment Variable or Config file
#'
#' Firstly environment variable PATH will be searched, if not
#' found, configure file then will be used to search the entry.
#'
#' @param entry a program or blast database.
#' @return a path.
#' @export
#' @seealso [set_blast_path] for setting blast path, [install_database]
#' for installing database.
find_path <- function(entry) {
  if (has_program(entry)) {
    entry
  } else {
    # Check config file
    prog_path <- load_from_config(entry)
    if (is.null(prog_path)) {
      stop("Path for ", paste(entry, collapse = " or "),
        " not found! Please make sure that the software is correctly installed and path is set.",
        call. = FALSE
      )
    } else {
      if (length(prog_path) > 1) {
        warning("More than 2 location found, only the first will be used.")
      }
      prog_path[1]
    }
  }
}

has_program <- function(program) {
  stopifnot(length(program) == 1)
  path <- Sys.which(program) %>% as.character()
  path != ""
}

add_path <- function(path) {
  current <- Sys.getenv("PATH")
  PATH <- unique(strsplit(current, ":")[[1]])
  PATH <- c(PATH, path[!path %in% PATH])
  Sys.setenv(PATH = paste(PATH, collapse = ":"))
}

# Save key, path pair to config file
save_to_config <- function(key, path) {
  if (!file.exists(.path_config_file)) {
    message("=> Config file ", .path_config_file, " does not exist, creating it...")
    if (!dir.exists(dirname(.path_config_file))) {
      dir.create(dirname(.path_config_file), recursive = TRUE)
    }
    file.create(.path_config_file)
  }
  new_to_add <- list(path)
  names(new_to_add) <- key
  configs <- yaml::yaml.load_file(.path_config_file)
  message(sprintf("=> Saving %s:%s to config file %s", key, path, .path_config_file))
  if (is.null(configs)) {
    # has no content
    yaml::write_yaml(new_to_add, file = .path_config_file)
  } else {
    # has content
    configs[[key]] <- path
    yaml::write_yaml(configs, file = .path_config_file)
  }
}

# load path from config file by key
load_from_config <- function(key) {
  if (!file.exists(.path_config_file)) {
    msg <- paste(
      sprintf(
        "Config file %s does not exist or permission denied, run the following command to chek!",
        .path_config_file
      ),
      "\tyaml::yaml.load_file(neopeptides:::.path_config_file)",
      sep = "\n"
    )
    stop(msg, call. = FALSE)
  }
  configs <- yaml::yaml.load_file(.path_config_file)
  if (!key %in% names(configs)) {
    stop("Key '", key, "' not found in config file ", .path_config_file)
  } else {
    configs[[key]]
  }
}

# Internal function to Smith-Waterman align two vectors of peptides.
SW_align <- function(pep1,
                     pep2,
                     gap_open = -11L,
                     gap_extend = -1L) {
  al <- Biostrings::pairwiseAlignment(pep1, pep2,
    substitutionMatrix = "BLOSUM62",
    gapOpening = gap_open,
    gapExtension = gap_extend,
    type = "local",
    scoreOnly = TRUE
  )

  if (length(al) == 0) al <- as.numeric(NA)

  return(al)
}

modeleR <- function(als, a = 26, k = 4.86936, dislike = FALSE) {
  be <- -k * (a - als)
  sumexp <- sum(exp(be))
  Zk <- 1 + sumexp
  R <- sumexp / Zk
  if (dislike) {
    return(1 - R)
  } else {
    return(R)
  }
}
