# fetch data from zenodo

library(httr)

url <- "https://zenodo.org/api/records/10139562/files-archive"


response <- GET(url, progress())

if (status_code(response) == 200) {
  writeBin(content(response, "raw"), "./tmp.zip")
  unzip(destination_zip, exdir = "./")
  cat("File downloaded successfully.\n")
} else {
  cat("Failed to download the file. Status code:", status_code(response), "\n")
}

