#/bin/bash 
mkdir software software/deps
chmod 755 software/deps

paquetesCRAN="'jsonlite' 'httr' 'devtools' 'sjmisc' 'tidyverse' 'networkD3' 'magrittr'"

for paqueteCRAN in $paquetesCRAN
	do 
		Rscript -e 'install.packages('$paqueteCRAN', repos="https://cran.rstudio.com/", lib = "software/deps" )'
	done 
