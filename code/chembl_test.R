#!/usr/bin/env RScript

# This scripts reads a list of compounds, identified by their chembl ids from a CSV file 
# given as an input. For each compound it find related targets, optionally filtered by
# organism. It saves a mepping between the compound and targets in the output CSV file.

#añadido
#librerias
library(jsonlite)
library(httr)
# install.packages('devtools')
library(devtools)
# devtools::install_github('rajarshi/chemblr/package')
library(chemblr)

URL_ROOT <- "https://www.ebi.ac.uk/chembl/api/data"

get_objects2 <- function(url){
  req <- GET(url)
  warn_for_status(req)
  json <- content(req, "text", encoding=parsed$encoding)
  validate(json)
  return(fromJSON(json))
}
direcc <- "https://www.ebi.ac.uk/chembl/api/data/molecule?limit=1000&max_phase=4&format=json"
pr <- get_objects2(direcc)
#hay que fusionar hasta coger los 3944
moleculas_pr <- pr[1]$molecules # moleculas
moleculas_1 <- moleculas_1[[1]] # datos primera molecula
moleculas_1$molecule_chembl_id

json_format <- "&format=json"
url_mechanism <- "mechanism?molecule_chembl_id=" #despues de esto va la id de la molécula
mecan_mol_1 <- file.path(URL_ROOT, paste(url_mechanism, moleculas_1$molecule_chembl_id, json_format, sep = ""))
mecan_1_json <- get_objects2(mecan_mol_1)
#la molecula 1 no tiene targets, buscamos con un ejemplo
mecan_ejemplo <- file.path(URL_ROOT, paste(url_mechanism, "CHEMBL998", json_format, sep = ""))
mecan_ejemplo_json <- get_objects2(mecan_ejemplo)
mecan_ejemplo_json$mechanisms[[1]]$target_chembl_id

#/añadido