#/bin/bash 
mkdir software software/deps
chmod 755 software/deps

paquetesCRAN="'jsonlite' 'httr' 'devtools' 'sjmisc' 'tidyverse' 'networkD3' 'magrittr' "

paquetesDevtools="'rajarshi/chemblr/package'"

for paqueteCRAN in $paquetesCRAN
	do 
		Rscript -e 'install.packages('$paqueteCRAN', repos="https://cran.rstudio.com/", lib = "software/deps" )'
	done 

for paqueteDevtools in $paquetesDevtools
	do 
		Rscript -e 'devtools::install_github('$paqueteDevtools', lib = "software/deps" )'
	done