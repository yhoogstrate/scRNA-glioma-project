# data per sample ----

data_per_sample <- read.table("data-raw/data_per_sample.tsv", sep = "\t", header = T) |> 
  dplyr::mutate(date_added = as.POSIXct(.data$date_added, format = "%d-%m-%Y"))
usethis::use_data(data_per_sample, compress = "xz", overwrite = T)

# run code, then:
# devtools::document()
# devtools::load_all()

