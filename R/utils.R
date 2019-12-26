.path_config_file <- file.path(
  system.file("extdata",
    package = "neopeptides"
  ),
  "path_config.yml"
)

find_program <- function(program) {
  path <- Sys.which(program)
  if (all(path == "")) {
    stop("programcutable for ", paste(program, collapse = " or "),
      " not found! Please make sure that the software is correctly installed or path is set.",
      call. = FALSE
    )
  }
  path[which(path != "")[1]]
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

set_path <- function(path, program, store = TRUE) {
  if (!dir.exists(path)) {
    stop(path, " does not exist, please check!")
  }
  add_path(path)
  if (store) {
    if (!file.exists(.path_config_file)) {
      file.create(.path_config_file)
    }
    new_to_add <- list(path)
    names(new_to_add) <- program
    configs <- yaml::yaml.load_file(.path_config_file)
    if (is.null(configs)) {
      # has no content
      yaml::write_yaml(new_to_add, file = .path_config_file)
    } else {
      # has content
      configs[[program]] <- path
    }
  }
}

init_path <- function() {
  paths <- yaml::yaml.load_file(.path_config_file)
  paths %>%
    purrr::reduce(base::c) %>%
    add_path()
}
