# Este script se ha obtenido a partir de un c√≥digo encontrado en https://gist.github.com/mnowotka/99a232900116df22be84ead82f234d9e.



#################################################################
#################### LIBRERIAS ##################################
#################################################################

#Las librerias se tienen que instalar en la carpeta SOFTWARE
library(jsonlite) ## Tiene que ser la ultima version (1.7.2)
library(httr)
# install.packages('devtools')
library(devtools)
# devtools::install_github('rajarshi/chemblr/package')
library(chemblr)
library(sjmisc)





#################################################################
################## CARGA PROTEINAS ##############################
#################################################################

proteinas_targets <- scan(file = "data/proteinas_chembl.txt", what = character())

# Proteinas STRING con interaction score = 0.9 y max number of interactors to show  no more than 5 interactors
proteinas_targets_segundo_grado <- scan(file = "data/proteinas_chembl_segundo_grado.txt", what = character())

# Probar con interaction score o more than 5 interactors



#################################################################
################## FUNCIONES ####################################
#################################################################

# Le pasamos una URL que es una CONSULTA y nos devuelve la info, que se encuentra en formato json
get_objects2 <- function(url){
  req <- GET(url)
  warn_for_status(req)
  json <- content(req, "text", encoding=DEFAULT_ENCODING)
  validate(json)
  return(fromJSON(json))
}


# Obtiene los ids de todos los f·rmacos de una fase indicada y guarda en una lista toda la info de los f·rmacos 
get_all_drugs <- function(fase_farmaco){
  
  limite_farmacos <- 1000 #N˙mero de f·rmacos por p·gina
  fase_farmacos <- fase_farmaco #Fase en el que se encuentra un f·rmaco (en fase 4 ya est· aprobado)
  offset_farmacos <- 0 #PosiciÛn del f·rmaco que se muestra el primero en la p·gina (pasar p·gina)
  
  URL_farmacos <- file.path(URL_ROOT, 
                            paste("molecule?limit=",limite_farmacos,
                                  "&offset=", offset_farmacos,
                                  "&max_phase=",fase_farmacos, 
                                  json_format, 
                                  sep = ""))
  
  farmacos_pagina_1 <- get_objects2(URL_farmacos)
  
  lista_info_farmacos <- list()
  lista_info_farmacos[[1]] <- farmacos_pagina_1
  vector_farmacos_id <- get_drugs_from_page(farmacos_pagina_1)
  
  
  # Numero total de farmacos
  num_total_farmacos <- round(farmacos_pagina_1$page_meta$total_count)
  # Numero de paginas si mostramos con pagina 1000
  num_paginas <- ceiling(num_total_farmacos/limite_farmacos)
  
  #Si hay m·s de una p·gina, hacemos el proceso anterior con las que haya
  for (pagina in 1:(num_paginas-1)){
    offset_farmacos <- offset_farmacos + 1000
    URL_farmacos <- file.path(URL_ROOT, 
                              paste("molecule?limit=",limite_farmacos,
                                    "&offset=", offset_farmacos,
                                    "&max_phase=",fase_farmacos, 
                                    json_format, 
                                    sep = ""))
    farmacos_pagina <- get_objects2(URL_farmacos)
    
    lista_info_farmacos[[(1+pagina)]] <- farmacos_pagina
    
    vector_farmacos_id_intermedio <- get_drugs_from_page(farmacos_pagina)
    vector_farmacos_id <- c(vector_farmacos_id, vector_farmacos_id_intermedio)
  }
  
  # Devolvemos un vector con todos los IDs de los f·rmacos
  return(list("vector_farmacos_ids" = vector_farmacos_id, "info_farmacos" = lista_info_farmacos))
}


# Pasandole una pagina con f·rmacos, devuelve un vector con los ids de los farmacos
get_drugs_from_page <- function(farmacos_page){
  vector_farmacos_id <- vector()
  for (posicion in 1:length(farmacos_page$molecules)){
    vector_farmacos_id <- c(vector_farmacos_id, 
                           farmacos_page$molecules[[posicion]]$molecule_chembl_id)
    
  }
  return(vector_farmacos_id)
}


# Le pasamos el vector con los ids de las drogas y las proteinas targets que nos interesa
# Devuelve una lista donde el nombre corresponde al target y el contenido ser· un vector con los f·rmacos asociados
# Devuelve un dataframe con las drogas seleccionadas y la informacion relacionada
get_targets <- function(drug_vector, proteinas_targets, drug_info){
  #Para devolver: lista y dataframe
  lista_target_drug <- list()
  informacion_farmacos <- data.frame("id_farmaco" = NA,
                                     "id_target" = NA, 
                                     "fecha_aprobaciÛn" = NA,
                                     "canonical_smile" = NA, 
                                     "action_type" = NA,
                                     "mechanism_of_action" = NA)
  
  ###### Valores para mostrar como va la ejecucion
  
  valor25 <- round(length(drug_vector)*0.25)
  valor50 <- round(length(drug_vector)*0.5)
  valor75 <- round(length(drug_vector)*0.75)
  
  for (posDrug in 1:length(drug_vector)){
    
    # Comprobacion Valores para mostrar como va la ejecucion
    ##############################################
    if(posDrug == valor25){
      print("Se ha alcanzado un 25% de la carga")
    }
    if(posDrug == valor50){
      print("Se ha alcanzado un 50% de la carga")
    }
    if(posDrug == valor75){
      print("Se ha alcanzado un 75% de la carga")
    }
   ################################################ 
    # Hacemos peticiÛn para traernos el mecanismo y conseguir los ids de los target (si tiene alguno)
    id_farmaco <- drug_vector[posDrug]
    URL_targets <- file.path(URL_ROOT, 
                             paste("mechanism?molecule_chembl_id=", 
                                   id_farmaco, 
                                   json_format, 
                                   sep = ""))
    mechanism <- get_objects2(URL_targets)
    num_targets <- mechanism$page_meta$total_count
    
    if(num_targets > 0){
      for(t in 1:(num_targets)){
        target <- mechanism$mechanisms[[t]]$target_chembl_id
        
        # Comprobamos que no sea null, ya que a veces aparece la etiqueta pero a NULL
        if(!is.null(target)){
          # Comprobamos que el target corresponde a una de nuestras proteinas objetivo
          if(match(target, proteinas_targets, nomatch = 0) != 0){
            # Usamos get_drug_info para completar el dataframe
            informacion_farmacos <- get_drug_info(tablaInfo = informacion_farmacos, 
                                                  drugID = id_farmaco, 
                                                  drugPos = posDrug,
                                                  drugInfo = drug_info,
                                                  mechanism = mechanism, 
                                                  targetpos = t, 
                                                  targetID = target)
            
            #Comprobamos si ya tenemos una entrada con ese target para la lista
            if(!is.null(lista_target_drug[[target]])){
              lista_target_drug[[target]] <- c(lista_target_drug[[target]], id_farmaco)
            }
            else{
              lista_target_drug[[target]] <- id_farmaco
            }
          }
        }
      }
    }
  }
  na.omit(informacion_farmacos)
  return(list("lista_target_drug" = lista_target_drug, "df_informacion_farmacos" = informacion_farmacos))
}

# Cargamos los datos que queremos de las drogas
# AÒade al dataframe la info de un f·rmaco. Esta informaciÛn puede venir de la consulta de todos los f·rmacos o de la consulta de los mecanismos
get_drug_info <- function(tablaInfo, drugID, drugPos, drugInfo, mechanism, targetpos, targetID){
  
  ######## Info del listado de f·rmacos ###########
  pagina <- ceiling(drugPos/1000)
  posicion_en_pagina <- (drugPos-((pagina-1)*1000))
  farmaco_info <- drugInfo[[pagina]]$molecules[[posicion_en_pagina]]
  
  # Fecha de aprobaciÛn del f·rmaco
  fecha_aprobacion <- farmaco_info$first_approval
  # Canonical smile
  c_smile <- farmaco_info$molecule_structures[1]
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Mirar para aÒadir m·s !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ######## Info del mecanismo ###########
  # Tipo de accion
  action_type <- mechanism$mechanisms[[targetpos]]$action_type
  # Mecanismo de acciÛn
  mechanism_of_action <- mechanism$mechanisms[[targetpos]]$mechanism_of_action
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Mirar para aÒadir m·s !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  tablaInfo <- rbind(tablaInfo, c("id_farmaco" = drugID,
                                  "id_target" = targetID, 
                                  "fecha_aprobaciÛn" = fecha_aprobacion,
                                  "canonical_smile" = c_smile,
                                  "action_type" = action_type,
                                  "mechanism_of_action" = mechanism_of_action))
  return(tablaInfo)
  
}




#################################################################
################### URLS ########################################
#################################################################

## Ruta raiz

URL_ROOT <- "https://www.ebi.ac.uk/chembl/api/data"
DEFAULT_ENCODING <- "UTF-8"





#################################################################
################### EJECUCION ###################################
#################################################################

# Todos los f·rmacos en fase 4. Tanto para proteinas de 1 y 2 grado

farmacos_fase4 <- get_all_drugs(fase_farmaco = 4)
farmacos_ids_fase4 <- farmacos_fase4[["vector_farmacos_ids"]]
farmacos_info_fase4 <- farmacos_fase4[["info_farmacos"]]

targets_proteinas_1g <- get_targets(drug_vector = farmacos_ids_fase4, proteinas_targets = proteinas_targets, drug_info = farmacos_info_fase4)
targets_proteinas_2g <- get_targets(drug_vector = farmacos_ids_fase4, proteinas_targets = proteinas_targets_segundo_grado, drug_info = farmacos_info_fase4)

# Todos los f·rmacos en fase 3. Tanto para proteinas de 1 y 2 grado

farmacos_fase3 <- get_all_drugs(fase_farmaco = 3)
farmacos_ids_fase3 <- farmacos_fase3[["vector_farmacos_ids"]]
farmacos_info_fase3 <- farmacos_fase3[["info_farmacos"]]

targets_proteinas_1g_f3 <- get_targets(drug_vector = farmacos_ids_fase3, proteinas_targets = proteinas_targets, drug_info = farmacos_info_fase3)
targets_proteinas_2g_f3 <- get_targets(drug_vector = farmacos_ids_fase3, proteinas_targets = proteinas_targets_segundo_grado, drug_info = farmacos_info_fase3)




#################################################################
################### Pruebas(borrar) #############################
#################################################################

limite_farmacos <- 1000 #N˙mero de f·rmacos por p·gina
fase_farmacos <- fase_farmaco #Fase en el que se encuentra un f·rmaco (en fase 4 ya est· aprobado)
offset_farmacos <- 0 #PosiciÛn del f·rmaco que se muestra el primero en la p·gina (pasar p·gina)

URL_farmacos <- file.path(URL_ROOT, 
                          paste("molecule?limit=",limite_farmacos,
                                "&offset=", offset_farmacos,
                                "&max_phase=",fase_farmacos, 
                                json_format, 
                                sep = ""))

farmacos_pagina_1 <- get_objects2(URL_farmacos)
# Nos traemos los primeros 1000 farmacos
farmacos_pagina_2 <- get_objects2(URL_farmacos)
# Numero total de farmacos
num_farmacos_tot <- round(farmacos_pagina_1$page_meta$total_count)
# Numero de paginas si mostramos con pagina 1000
num_paginas <- ceiling(num_farmacos_tot/limite_farmacos)

# Chembl id de los farmacos
farmaco1_id <- farmacos_pagina_1$molecules[[1]]$molecule_chembl_id
farmaco1_molecula <- farmacos_pagina_1$molecules[[1]]


farmaco2_id <- farmacos_pagina_1$molecules[[2]]$molecule_chembl_id
# ...
farmaco1000_id <- farmacos_pagina_1$molecules[[1000]]$molecule_chembl_id

# Teniendo el id de un farmaco, podemos ver sus target
id_farmaco <- farmaco1_id
URL_targets <- file.path(URL_ROOT, 
                         paste("mechanism?molecule_chembl_id=", 
                               id_farmaco, 
                               json_format, 
                               sep = ""))
mecanismo_farmaco_1 <- get_objects2(URL_targets)
num_mecanismos_farmaco1 <- mecanismo_farmaco_1$page_meta$total_count

targets_farmaco1 <- mecanismo_farmaco_1$mechanisms[[1]]$target_chembl_id

id_farmaco <- farmaco2_id
URL_targets <- file.path(URL_ROOT, 
                         paste("mechanism?molecule_chembl_id=", 
                               id_farmaco, 
                               json_format, 
                               sep = ""))
mecanismo_farmaco_2 <- get_objects2(URL_targets)
num_targets_farmaco2 <- length(mecanismo_farmaco_2$mechanisms[[1]]$target_chembl_id)
targets_farmaco2 <- mecanismo_farmaco_2$mechanisms[[1]]$target_chembl_id

id_farmaco <- "CHEMBL9"
URL_targets <- file.path(URL_ROOT, 
                         paste("mechanism?molecule_chembl_id=", 
                               id_farmaco, 
                               json_format, 
                               sep = ""))
mecanismo_farmaco_X <- get_objects2(URL_targets)
num_targets_farmacoX <- mecanismo_farmaco_X$page_meta$total_count
mecanismo_farmaco_X$mechanisms[[2]]$target_chembl_id

