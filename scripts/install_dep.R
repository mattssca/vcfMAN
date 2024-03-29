#dependencies
install.packages("stringr", repos = 'http://cran.us.r-project.org')
install.packages("table1", repos = 'http://cran.us.r-project.org')
install.packages("dplyr", repos = 'http://cran.us.r-project.org')
install.packages("knitr",repos = 'http://cran.us.r-project.org')
install.packages("devtools", repos = 'http://cran.us.r-project.org')
install.packages("gridExtra", repos = 'http://cran.us.r-project.org')
install.packages("BiocManager", repos = 'http://cran.us.r-project.org')
install.packages("openxlsx", repos = 'http://cran.us.r-project.org')
install.packages("RCircos", repos = 'http://cran.us.r-project.org')
install.packages("webshot", repos = 'http://cran.us.r-project.org')
install.packages("psych", repos = 'http://cran.us.r-project.org')
install.packages("data.table", repos = 'http://cran.us.r-project.org')
BiocManager::install("karyoploteR", force = TRUE, dependencies = TRUE)
devtools::install_github('Mikata-Project/ggthemr')
webshot::install_phantomjs()
