library(FactoMineR)
library(factoextra)
library(PCAtools)
###############did not work#############
install.packages('reticulate')
#loading required R libraries 
library(reticulate) #the superpower bridges python and R
#specifying which version of python to use
use_python('/usr/local/Cellar/python/3.7.7/bin/python3')
py_config()
#importing required Python libraries/modules
library(reticulate)
py_install("pandas")
py_install("seaborn")
py_install("matplotlib.pyplot")
sns <- import('seaborn')
plt <- import('matplotlib.pyplot')
pd <- import('pandas')
#building a seaborn pairplot using pairplot()
sns$pairplot(r_to_py(filtered_grouped_nona), hue = 'combined_animals')
#display the plot
plt$show()

######################################################


pairs(select_if(filtered_grouped_nona[,-1], is.numeric))
heatmap(x, scale = "column")
tiff("~/Box Sync/public_repertoire/results/pair_plot.tiff")


biplot(p)
pca_asd <- PCA(x)
library(dplyr)
library(GGally)
filtered_grouped_nona <- filtered_grouped_nona %>%
  mutate(E11_or_E16 = ifelse(grepl("^E11, E16$",combined_animals),"YES", "NO"))

ggpairs(filtered_grouped_nona, columns = 4:21, 
        ggplot2::aes(colour=E11_and_E16, alpha = 0.7)) 


mutate(specificity = str_detect("PreF")

library(ggplot2)
install.packages("GGally")

msa()



filtered_grouped %>% ggplot(aes(x = E11_and_E16, y = n, fill = spec)) + 
  geom_bar(stat = "identity", position = "dodge", width = .5)

