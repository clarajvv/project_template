# Este script se ha obtenido a partir de un cÃƒÂ³digo encontrado en https://gist.github.com/mnowotka/99a232900116df22be84ead82f234d9e.



#################################################################
#################### LIBRERIAS ##################################
#################################################################


dir <- getwd()
libr <- paste(dir,"/software/deps", sep = "")
.libPaths(c(libr, .libPaths()))



#Las librerias se tienen que instalar en la carpeta SOFTWARE
library(jsonlite) ## Tiene que ser la ultima version (1.7.2)
library(httr)
library(devtools)
library(chemblr)
library(sjmisc)
library(tidyverse)
library(networkD3)
library(magrittr)
library(ggplot2)

# Meter nueva!!!!!!!!!!!!!!!!!!
library(tidyverse)




#################################################################
################### URLS ########################################
#################################################################

## VARIABLES GLOBALES

URL_ROOT <- "https://www.ebi.ac.uk/chembl/api/data"
DEFAULT_ENCODING <- "UTF-8"
JSONFORMAT <- "&format=json"




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


# Obtiene los ids de todos los fÃ¡rmacos de una fase indicada y guarda en una lista toda la info de los fÃ¡rmacos 
get_all_drugs <- function(fase_farmaco){
  
  limite_farmacos <- 1000 #NÃºmero de fÃ¡rmacos por pÃ¡gina
  fase_farmacos <- fase_farmaco #Fase en el que se encuentra un fÃ¡rmaco (en fase 4 ya estÃ¡ aprobado)
  offset_farmacos <- 0 #PosiciÃ³n del fÃ¡rmaco que se muestra el primero en la pÃ¡gina (pasar pÃ¡gina)
  
  URL_farmacos <- file.path(URL_ROOT, 
                            paste("molecule?limit=",limite_farmacos,
                                  "&offset=", offset_farmacos,
                                  "&max_phase=",fase_farmacos, 
                                  JSONFORMAT, 
                                  sep = ""))
  
  farmacos_pagina_1 <- get_objects2(URL_farmacos)
  
  lista_info_farmacos <- list()
  lista_info_farmacos[[1]] <- farmacos_pagina_1
  vector_farmacos_id <- get_drugs_from_page(farmacos_pagina_1)
  
  
  # Numero total de farmacos
  num_total_farmacos <- round(farmacos_pagina_1$page_meta$total_count)
  # Numero de paginas si mostramos con pagina 1000
  num_paginas <- ceiling(num_total_farmacos/limite_farmacos)
  
  #Si hay mÃ¡s de una pÃ¡gina, hacemos el proceso anterior con las que haya
  for (pagina in 1:(num_paginas-1)){
    offset_farmacos <- offset_farmacos + 1000
    URL_farmacos <- file.path(URL_ROOT, 
                              paste("molecule?limit=",limite_farmacos,
                                    "&offset=", offset_farmacos,
                                    "&max_phase=",fase_farmacos, 
                                    JSONFORMAT, 
                                    sep = ""))
    farmacos_pagina <- get_objects2(URL_farmacos)
    
    lista_info_farmacos[[(1+pagina)]] <- farmacos_pagina
    
    vector_farmacos_id_intermedio <- get_drugs_from_page(farmacos_pagina)
    vector_farmacos_id <- c(vector_farmacos_id, vector_farmacos_id_intermedio)
  }
  
  # Devolvemos un vector con todos los IDs de los fÃ¡rmacos
  return(list("vector_farmacos_ids" = vector_farmacos_id, "info_farmacos" = lista_info_farmacos))
}


# Pasandole una pagina con fÃ¡rmacos, devuelve un vector con los ids de los farmacos
get_drugs_from_page <- function(farmacos_page){
  vector_farmacos_id <- vector()
  for (posicion in 1:length(farmacos_page$molecules)){
    vector_farmacos_id <- c(vector_farmacos_id, 
                           farmacos_page$molecules[[posicion]]$molecule_chembl_id)
    
  }
  return(vector_farmacos_id)
}


# Le pasamos el vector con los ids de las drogas y las proteinas targets que nos interesa
# Devuelve una lista donde el nombre corresponde al target y el contenido serÃ¡ un vector con los fÃ¡rmacos asociados
# Devuelve un dataframe con las drogas seleccionadas y la informacion relacionada
get_targets <- function(drug_vector, proteinas_targets, drug_info){
  #Para devolver: lista y dos dataframe
  lista_target_drug <- list()
  informacion_farmacos <- data.frame("id_farmaco" = NA,
                                     "id_target" = NA, 
                                     "fecha_aprobaciÃ³n" = NA,
                                     "canonical_smile" = NA, 
                                     "action_type" = NA,
                                     "mechanism_of_action" = NA)
  informacion_quimica <- data.frame("alogp" = NA, 
                                    "aromatic_rings" = NA,
                                    "cx_logd" = NA, 
                                    "cx_logp" = NA, 
                                    "cx_most_apka" = NA, 
                                    "cx_most_bpka" = NA, 
                                    "hba" = NA, 
                                    "hba_lipinski" = NA,
                                    "hbd" = NA, 
                                    "hbd_lipisnki" = NA, 
                                    "heavy_atoms" = NA, 
                                    "molecular_species" = NA, 
                                    "mw_freebase" = NA, 
                                    "mw_monoistopic" = NA, 
                                    "psa" = NA, 
                                    "qed_weight" = NA, 
                                    "ro3_pass" = NA, 
                                    "rtb" = NA)
  
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
    # Hacemos peticiÃ³n para traernos el mecanismo y conseguir los ids de los target (si tiene alguno)
    id_farmaco <- drug_vector[posDrug]
    URL_targets <- file.path(URL_ROOT, 
                             paste("mechanism?molecule_chembl_id=", 
                                   id_farmaco, 
                                   JSONFORMAT, 
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
            tablas <- get_drug_info(tablaInfo = informacion_farmacos,
                                                  tablaQuimica = informacion_quimica,
                                                  drugID = id_farmaco, 
                                                  drugPos = posDrug,
                                                  drugInfo = drug_info,
                                                  mechanism = mechanism, 
                                                  targetpos = t, 
                                                  targetID = target)
            informacion_farmacos <- tablas$tablaInfo
            informacion_quimica <- tablas$tablaQuimica
            
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
  informacion_farmacos <- na.omit(informacion_farmacos)
  informacion_quimica <- na.omit(informacion_quimica)
  return(list("lista_target_drug" = lista_target_drug, "df_informacion_farmacos" = informacion_farmacos, "df_informacion_quimica" = informacion_quimica))
}

# Cargamos los datos que queremos de las drogas
# AÃ±ade al dataframe la info de un fÃ¡rmaco. Esta informaciÃ³n puede venir de la consulta de todos los fÃ¡rmacos o de la consulta de los mecanismos
get_drug_info <- function(tablaInfo, tablaQuimica, drugID, drugPos, drugInfo, mechanism, targetpos, targetID){
  
  ######## Info del listado de fÃ¡rmacos ###########
  pagina <- ceiling(drugPos/1000)
  posicion_en_pagina <- (drugPos-((pagina-1)*1000))
  farmaco_info <- drugInfo[[pagina]]$molecules[[posicion_en_pagina]]
  
  # Fecha de aprobaciÃ³n del fÃ¡rmaco
  fecha_aprobacion <- farmaco_info$first_approval
  # Canonical smile
  c_smile <- farmaco_info$molecule_structures[1]
  
  ####### Info quimica ##############################
  # alogp
  alogp <- farmaco_info$molecule_properties$alogp
  
  # aromatic rings
  aro_rings <- farmaco_info$molecule_properties$aromatic_rings
  
  # cx_logd
  cx_logd <- farmaco_info$molecule_properties$cx_logd
  
  # cx_logp
  cx_logp <- farmaco_info$molecule_properties$cx_logp
  
  # cx_most_apka
  cx_most_apka <- farmaco_info$molecule_properties$cx_most_apka
  
  # cx_most_bpka
  cx_most_bpka <- farmaco_info$molecule_properties$cx_most_bpka
  
  # hba 
  hba <- farmaco_info$molecule_properties$hba
  
  # hba_lipinski
  hba_lipinski <- farmaco_info$molecule_properties$hba_lipinski
  
  # hbd
  hbd <- farmaco_info$molecule_properties$hbd
  
  # hbd_lipinski
  hbd_lipinski <- farmaco_info$molecule_properties$hbd_lipinski
  
  # heavy_atoms
  heavy_atoms <- farmaco_info$molecule_properties$heavy_atoms
  
  # molecular_species
  molecular_species <- farmaco_info$molecule_properties$molecular_species
  
  # mw_freebase
  mw_freebase <- farmaco_info$molecule_properties$mw_freebase
  
  # mw_monoistopic
  mw_monoistopic <- farmaco_info$molecule_properties$mw_monoistopic
  
  # psa
  psa <- farmaco_info$molecule_properties$psa
  
  # qed_weight
  qed_weight <- farmaco_info$molecule_properties$qed_weight
  
  # ro3_pass
  ro3_pass <- farmaco_info$molecule_properties$ro3_pass
  
  # rtb
  rtb <- farmaco_info$molecule_properties$rtb
  
  ######## Info del mecanismo ###########
  # Tipo de accion
  action_type <- mechanism$mechanisms[[targetpos]]$action_type
  # Mecanismo de acciÃ³n
  mechanism_of_action <- mechanism$mechanisms[[targetpos]]$mechanism_of_action
  
 
  
  tablaInfo <- rbind(tablaInfo, c("id_farmaco" = drugID,
                                  "id_target" = targetID, 
                                  "fecha_aprobaciÃ³n" = fecha_aprobacion,
                                  "canonical_smile" = c_smile,
                                  "action_type" = action_type,
                                  "mechanism_of_action" = mechanism_of_action))
  tablaQuimica <- rbind(tablaQuimica, c("id_farmaco" = drugID,
                                        "id_target" = targetID, 
                                        "alogp" = alogp, 
                                        "aromatic_rings" = aro_rings,
                                        "cx_logd" = cx_logd, 
                                        "cx_logp" = cx_logp, 
                                        "cx_most_apka" = cx_most_apka, 
                                        "cx_most_bpka" = cx_most_bpka, 
                                        "hba" = hba, 
                                        "hba_lipinski" = hba_lipinski,
                                        "hbd" = hbd, 
                                        "hbd_lipisnki" = hbd_lipinski, 
                                        "heavy_atoms" = heavy_atoms, 
                                        "molecular_species" = molecular_species, 
                                        "mw_freebase" = mw_freebase, 
                                        "mw_monoistopic" = mw_monoistopic, 
                                        "psa" = psa, 
                                        "qed_weight" = qed_weight, 
                                        "ro3_pass" = ro3_pass, 
                                        "rtb" = rtb))
  
  return(list("tablaInfo" = tablaInfo, "tablaQuimica" = tablaQuimica))
  
}





#################################################################
################### EJECUCION ###################################
#################################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Debe indicar la fase del fármaco", call.=FALSE)
} else if (length(args)==1) {
  if(args == 1) {
    faseFarmaco = "3"
    proteinas_targets <- scan(file = "data/proteinas_chembl.txt", what = character())
    grado_proteinas <- 1
  }
  else if(args == 2) {
    faseFarmaco = "3"
    proteinas_targets <- scan(file = "data/proteinas_chembl_secundarias.txt", what = character())
    grado_proteinas <- 2
  }
  else if (args == 3) {
    faseFarmaco = "4"
    proteinas_targets <- scan(file = "data/proteinas_chembl.txt", what = character())
    grado_proteinas <- 1
  } 
  else if (args == 4) {
    faseFarmaco = "4"
    proteinas_targets <- scan(file = "data/proteinas_chembl_secundarias.txt", what = character())
    grado_proteinas <- 2
  }
}
faseFarmaco <- 4
grado_proteinas <- 1
proteinas_targets <- scan(file = "data/proteinas_chembl.txt", what = character())

print(paste("Comienza la búsqueda de fármacos en la fase de ensayo ", faseFarmaco, sep = ""))
farmacos_fase <- get_all_drugs(fase_farmaco = faseFarmaco)
farmacos_ids_fase <- farmacos_fase[["vector_farmacos_ids"]]
farmacos_info_fase <- farmacos_fase[["info_farmacos"]]
print(paste("Fin de la búsqueda de fármacos en la fase de ensayo ", faseFarmaco, sep = ""))

print("Comienza el filtrado de fármacos para las proteínas seleccionadas")
targets_proteinas <- get_targets(drug_vector = farmacos_ids_fase, proteinas_targets = proteinas_targets, drug_info = farmacos_info_fase)
print("Fin del filtrado de fármacos para las proteínas seleccionadas")




#################################################################
#################### GUARDAR DATOS ##############################
#################################################################

# Creamos una carpeta nueva en resultados con el día y hora de la ejecución
HOY  <- format(Sys.time(), "%F_%H.%M.%S")
nombre_carpeta <- paste(HOY, "_fase", faseFarmaco, "_grado", grado_proteinas, sep = "" )
directorio_carpeta_resultados <- paste("../results", nombre_carpeta, sep = "/")
dir.create(directorio_carpeta_resultados)


save_info <- function(datos, titulo, directorio){
  if(dim(datos)[1] >= 1){
    file_name_dir <- paste(directorio, "/", titulo, ".csv", sep = "")
    write.csv(datos, file_name_dir, row.names = FALSE)
  }
  else{
     print(paste("No hay información para guardar en: ", titulo, sep = ""))
  }
}

save_info(targets_proteinas$df_informacion_farmacos, "informacionGeneralFarmacos", directorio_carpeta_resultados)
save_info(targets_proteinas$df_informacion_quimica, "informacionQuimicaFarmacos", directorio_carpeta_resultados)



#################################################################
####################### GRAFICAS ################################
#################################################################

if(dim(targets_proteinas$df_informacion_farmacos)[1] >=1){
  # Si hay datos, se hacen las gráficas
  
  ##################################### Gráfica actionType ############################
  
  dfActionType <- data.frame("actionType" = targets_proteinas$df_informacion_farmacos$action_type)
  dfActionType <- as.data.frame(dfActionType %>% group_by(actionType) %>% tally())
  
  file_name_actionType <- paste(directorio_carpeta_resultados, "/graficaTipoDeAccion.png", sep = "")
  png(file_name_actionType, width = 1050, height = 1200) 
  
  g <- ggplot(data = dfActionType, aes(x = actionType, y = n, fill = actionType)) + geom_bar(stat="identity", width = 0.6) +
    theme_minimal()+ geom_text(aes(label = n), vjust=1.6, color="white", size=3) +
    theme (text = element_text(size=30)) + ggtitle ("Frecuencia del tipo de acción")  + 
    theme(plot.title = element_text(size=rel(2),
                                    vjust=3)) +
    scale_fill_brewer(palette="Paired") +
    labs(x = "Tipo de acción", y = "Frecuencia de aparición") + 
    theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    theme(axis.title.x = element_text(vjust=-0.5, size=rel(1.5))) +
    theme(axis.title.y = element_text(vjust=1.5, size=rel(1.5))) 
  
  print(g)
  
  dev.off()
  
  ##################################### Gráfica  mechanism of action ############################
  
  dfMechanism <- data.frame("Mechanism_of_action" = targets_proteinas$df_informacion_farmacos$mechanism_of_action)
  dfMechanism <- as.data.frame(dfMechanism %>% group_by(Mechanism_of_action) %>% tally())
  
  file_name_mechanism <- paste(directorio_carpeta_resultados, "/graficaMecanismoDeAccion.png", sep = "")
  png(file_name_mechanism, width = 1250, height = 1250) 
  
  g2 <- ggplot(data = dfMechanism, aes(x = Mechanism_of_action, y = n, fill = Mechanism_of_action)) + geom_bar(stat="identity", width = 0.6) +
    theme_minimal() +
    theme (text = element_text(size=20)) + ggtitle ("Frecuencia del mecanismo de acción")  + 
    theme(plot.title = element_text(size=rel(2),
                                    vjust=3)) +
    labs(x = "Mecanismo de acción", y = "Frecuencia") + 
    scale_colour_gradient2() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(axis.title.x = element_text(vjust=1.5, size=rel(1.5))) +
    theme(axis.title.y = element_text(vjust=1.5, size=rel(1.5))) 
  
  print(g2)
  
  dev.off()
  
}else{
  print("El númrero de fármacos obtenido es insuficiente para crear gráficas")
}




#red de interacciÃ³n

# CAMBIAR EL FORMATO DE LA RED (PONERLO BONITO) Y QUITAR QUE SEAN LOS DOS GRADOS DE LA PROTEÍNA(O PRIMER GRADO O SEGUNDO GRADO)
view_Interaction <- function(type_1, type_2, day){
   proteins <- c()
   list_name <- c()
   targets <- c()
   type <- c(type_1, type_2)
   
   if(length(type) >= 2){
       for(i in 1:length(type)){
           target <- attributes(type[i])
           targets <- append(targets, target$names)
           }
           for(i in 1:length(targets)){
               for(j in 1:length(type[[i]])){
                   prot <- type[[i]][j]
                   proteins <- append(proteins, prot)
                   name <- targets[i]
                   list_name <- append(list_name, name)
               }
             }
           inter <- data.frame(list_name, proteins)
           bsk.network<-graph.data.frame(inter, directed=F)
           img = img = paste("Red-medicamento-proteina", day, sep = "-")
           dir <- file_name <- paste("../results",img, sep = "/")
           file_name <- paste(dir, "jpeg", sep = ".")
           print(file_name)
           jpeg(file_name)
           plot(bsk.network)
       }else if(length(type) <= 1){
           print("El medicamento no tiene proteínas interaccionando")
         }
}

#interact_graphic1<- view_Interaction(targets_proteinas_1g$lista_target_drug, targets_proteinas_2g$lista_target_drug, HOY)






