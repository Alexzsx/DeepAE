
if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}

library(RCy3)
library(magrittr)
library(stringr)


cytoscapePing()
RCy3::commandsRun('wikipathways import-as-network id=WP295')
RCy3::commandsRun('wikipathways import-as-pathway id=WP295')


RCy3::commandsRun('wikipathways import-as-network id=WP4396')
RCy3::commandsRun('wikipathways import-as-pathway id=WP4396')


RCy3::commandsRun('wikipathways import-as-network id=WP3888')
RCy3::commandsRun('wikipathways import-as-pathway id=WP3888')


RCy3::commandsRun('wikipathways import-as-network id=WP307')
RCy3::commandsRun('wikipathways import-as-pathway id=WP307')


RCy3::commandsRun('wikipathways import-as-network id=WP4239')
RCy3::commandsRun('wikipathways import-as-pathway id=WP4239')



RCy3::commandsRun('wikipathways import-as-network id=WP43')
RCy3::commandsRun('wikipathways import-as-pathway id=WP43')


RCy3::commandsRun('wikipathways import-as-network id=WP43')
RCy3::commandsRun('wikipathways import-as-pathway id=WP43')


RCy3::commandsRun('wikipathways import-as-network id=WP2185')
RCy3::commandsRun('wikipathways import-as-pathway id=WP2185')


RCy3::commandsRun('wikipathways import-as-network id=WP210')
RCy3::commandsRun('wikipathways import-as-pathway id=WP210')


