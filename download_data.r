# fetch data from zenodo
library(curl)

url <- "https://zenodo.org/api/records/10139562/files-archive"


curl_download(url, destfile = "./tmp.zip", quiet = FALSE)

unzip("./tmp.zip", exdir = "./")
cat("File downloaded successfully.\n")



