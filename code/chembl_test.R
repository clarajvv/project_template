#este script se ha obtenido a partir de un código encontrado en https://gist.github.com/mnowotka/99a232900116df22be84ead82f234d9e.

#librerias
library(jsonlite)
library(httr)
# install.packages('devtools')
library(devtools)
# devtools::install_github('rajarshi/chemblr/package')
library(chemblr)

URL_ROOT <- "https://www.ebi.ac.uk/chembl/api/data"
DEFAULT_ENCODING <- "UTF-8"

get_objects2 <- function(url){
  req <- GET(url)
  warn_for_status(req)
  json <- content(req, "text", encoding=DEFAULT_ENCODING)
  validate(json)
  return(fromJSON(json))
}
direcc <- file.path(URL_ROOT, "molecule?limit=1000&max_phase=4&format=json")
#"https://www.ebi.ac.uk/chembl/api/data/molecule?limit=1000&max_phase=4&format=json" es la direccion de la primera pagina de los fármacos
pr <- get_objects2(direcc)
#hay que fusionar hasta coger los 3944
moleculas_pr <- pr[1]$molecules # moleculas
moleculas_1 <- moleculas_pr[[1]] # datos primera molecula
moleculas_1$molecule_chembl_id # el id del primer fármaco
moleculas_2 <- moleculas_pr[[2]] #datos segunda molécula

json_format <- "&format=json" # cadena para que la web te devuelva json
url_mechanism <- "mechanism?molecule_chembl_id=" # cadena para obtener los mecanismos, despues de esto va la id de la molécula
mecan_mol_1 <- file.path(URL_ROOT, paste(url_mechanism, moleculas_1$molecule_chembl_id, json_format, sep = "")) # direccion mecanismos del primer fármaco
mecan_1_json <- get_objects2(mecan_mol_1) # mecanismos del primer fármaco
mecan_1_json$mechanisms[[1]]$target_chembl_id #id del target del primer fármaco (este solo tiene uno, pero podría tener más)

#la molecula 1 no tiene targets, buscamos con un ejemplo
mecan_ejemplo <- file.path(URL_ROOT, paste(url_mechanism, "CHEMBL998", json_format, sep = ""))
mecan_ejemplo_json <- get_objects2(mecan_ejemplo)
mecan_ejemplo_json$mechanisms[[1]]$target_chembl_id
