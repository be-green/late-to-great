library(haven)
library(data.table)
library(stringr)

dta_files <- list.files("data/raw/OHIE/OHIE_Data")

make_csv <- function(dta_file, inpath, outpath) {
  filename <- str_replace_all(dta_file, "\\.dta", "")
  fwrite(
    as.data.table(
      haven::read_dta(file.path(inpath, dta_file))
    ),
    file.path(outpath, paste0(filename,".csv"))
  )
}

lapply(dta_files, make_csv,
       inpath = "data/raw/OHIE/OHIE_Data",
       outpath = "data/interim/OHIE")

