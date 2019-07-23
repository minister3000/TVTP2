

#=============================
# Attach DATA to Package
#=============================

# read raw data in as data.frame
RawEconomicData <- read.csv(file="C:\\Users\\tgummer\\Desktop\\Uni HH\\Tobias\\R Package for DISS\\Project\\TVTP\\ExcludeFromCompile\\RawEconomicData.csv", header = TRUE, row.names = 1, sep=",")
# format and save the data suitable for package building under the name provided here. File will be saved in the \data folder of the R package
devtools::use_data(RawEconomicData, pkg=".", overwrite = TRUE)
# generate required data description file in the main folder. Edit that file if desired (you have to at least provide a title) and move it into the \man folder manually
prompt("RawEconomicData")
# --> Rebuild the R package in RStudio
# --> When package is loaded through 'library(TVTP)', the dataset automatically becomes available and does not need to get loaded or extracted from the package separately.
#     just type 'RawEconomicData'

# print all available datasets in the 'TVTP' package:
data(package = "TVTP")

# this will provide the 'help' information that we provided in the data description file above:
help(RawEconomicData)
