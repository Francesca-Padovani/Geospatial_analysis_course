setwd("/Users/FrancescaPadovani/Desktop/MAGISTRALE DATA SCIENCE/SECONDO ANNO/GEOSPATIAL ANALYSIS/PROGETTO_PROVE/Francesca_Padovani_Geospatial/shapefile/")
library(rgdal)
sf_use_s2(FALSE)
library(spdep)

#open the shapefile created with Jupyter Notebook
rome <- readOGR("final_dataset.shp")
summary(rome)
dim(rome)
rome@data
head(rome@data)

plot(rome)
text(coordinates(rome), labels=rome$ids, cex=0.6, col="blue")
title(main = "Neighbourhoods of Rome")

#compute the centroid 
centroids <- coordinates(rome)

plot(rome, border="blue") 
points(centroids, cex=0.8)


#By referring to the distances among centroids we can then define 
#the neighborhood relationships among the spatial units. 

# k-nearest neighbors

#with k = 1
knn1 <- knn2nb(knearneigh(centroids, k=1, longlat = T))
plot(rome, border="grey")
plot(knn1, centroids, col="red", add=TRUE)

#with k = 4,
knn4 <- knn2nb(knearneigh(centroids,k=4,longlat=T))
plot(rome, border="grey")
plot(knn4, centroids, col = 'red',add=TRUE) 


#Critical cut-off neighborhood criterion 

knn1RO <- knn2nb(knearneigh(centroids,k=1,longlat=T))
all.linkedT <- max(unlist(nbdists(knn1RO, centroids, longlat=T))) 
all.linkedT

#and we found that the minimum threshold distance is equal to 9.25 km
dnb10 <- dnearneigh(centroids, 0, 10, longlat=TRUE); dnb9
dnb13 <- dnearneigh(centroids, 0, 13, longlat=TRUE); dnb13
dnb15 <- dnearneigh(centroids, 0, 15, longlat=TRUE); dnb15 
dnb18 <- dnearneigh(centroids, 0, 18, longlat=TRUE); dnb18 


#as the cut-off distance increases, the number of links grows rapidly.
plot(rome, border="grey",xlab="",ylab="",xlim=NULL)
title(main="d nearest neighbours, d = 9-18") 
plot(dnb10, centroids, add=TRUE, col="blue")
plot(dnb13, centroids, add=TRUE, col="red")
plot(dnb15, centroids, add=TRUE, col="yellow")
plot(dnb18, centroids, add=TRUE, col="gold")


#contiguity-based neighbourhood
contnb_q <- poly2nb(rome, queen=T)  
contnb_q
plot(rome, border="grey")
plot(contnb_q, centroids, add=TRUE)
#########
#########



# Defining spatial weights

dnb1.listw <- nb2listw(dnb10, style = "W")
dnb2.listw <- nb2listw(dnb13, style = "W")
dnb3.listw <- nb2listw(dnb15, style = "W")
dnb4.listw <- nb2listw(dnb18, style = "W")
dnbneigh.listw <- nb2listw(knn4, style = 'W') 

listw2mat(dnb1.listw)

#Build free form matrix 
distM <- as.matrix(dist(centroids)) 
distM   

#three possible types of free form weight matrices
W1 <- 1/(1+(distM))
diag(W1) <- 0                       

W2 <- 1/(1+distM)^2
diag(W2) <- 0

W3 <- exp(-distM^2)
diag(W3) <- 0


#Row-standardize them 
W1s <- W1/rowSums(W1) 
W2s <- W2/rowSums(W2) 
W3s <- W3/rowSums(W3) 


#We can convert the weight matrix into a "listw" object
listW1s <- mat2listw(W1s)
listW2s <- mat2listw(W2s)
listW3s <- mat2listw(W3s)



### GLOBAL SPATIAL AUTOCORRELATION

#The visual inspection of the spatial quantile distribution of the growth rate may or may not suggest 
#the presence of some form of spatial dependence. 

brks <- round(quantile(rome$avg_price_, digits=3))
colours <- grey((length(brks):2)/length(brks))
plot(rome, col=colours[findInterval(rome$avg_price_, brks, all.inside=TRUE)])
title(main = "Average Airbnb price")

#in this case it seems that tehre's a bit of spatial dependence on the left neighborhoods.

#Moran's I test under the assumption of normality
moran.test(rome$avg_price_, dnb1.listw, randomisation=FALSE) 
moran.test(rome$avg_price_, dnb2.listw, randomisation=FALSE)
moran.test(rome$avg_price_, dnb3.listw, randomisation=FALSE)
moran.test(rome$avg_price_, dnb4.listw, randomisation=FALSE)
moran.test(rome$avg_price_, dnbneigh.listw, randomisation=FALSE)
moran.test(rome$avg_price_, listW1s, randomisation=FALSE)
moran.test(rome$avg_price_, listW2s, randomisation=FALSE)
moran.test(rome$avg_price_, listW3s, randomisation=FALSE)


#Moran's I test under the assumption of randomization
moran.test(rome$avg_price_, dnb1.listw, randomisation=TRUE)
moran.test(rome$avg_price_, dnb2.listw, randomisation=TRUE)
moran.test(rome$avg_price_, dnb3.listw, randomisation=TRUE)
moran.test(rome$avg_price_, dnb4.listw, randomisation=TRUE)
moran.test(rome$avg_price_, dnbneigh.listw, randomisation=TRUE)
moran.test(rome$avg_price_, listW1s, randomisation=TRUE)
moran.test(rome$avg_price_, listW2s, randomisation=TRUE)
moran.test(rome$avg_price_, listW3s, randomisation=TRUE)

#Moran's I test under permutation bootstrap
moran.mc(rome$avg_price_, dnb1.listw, nsim=999)
moran.mc(rome$avg_price_, dnb2.listw, nsim=999)
moran.mc(rome$avg_price_, dnb3.listw, nsim=999)
moran.mc(rome$avg_price_, dnb4.listw, nsim=999)
moran.mc(rome$avg_price_, dnbneigh.listw, nsim=999)
moran.mc(rome$avg_price_, listW1s, nsim=999)
moran.mc(rome$avg_price_, listW2s, nsim=999)
moran.mc(rome$avg_price_, listW3s, nsim=999)


##LOCAL SPATIAL AUTOCORRELATION Moran's Scatterplot

mplot1 <- moran.plot(rome$avg_price_, listw=dnbneigh.listw , main="1st Moran scatterplot \n from 4(k)nn matrix", xlab="Average Airbnb price", ylab="Spatial lags", return_df=F)
mplot2 <- moran.plot(rome$avg_price_, listw=dnb2.listw, main="2nd Moran scatterplot \n from critical cut-off W ", xlab="Average Airbnb price", ylab="Spatial lags", return_df=F)
mplot3 <- moran.plot(rome$avg_price_, listw=listW2s, main="3rd Moran scatterplot \n from free form W", xlab="Average Airbnb price", ylab="Spatial lags", return_df=F)

grid()
summary(mplot3) 
#we can see which are the variables that contributed most according to some criteria (Cook distance, covariance etc..)
hotspot <- as.numeric(row.names(as.data.frame(summary(mplot1))))
hotspot


#MORAN'S I INDEX to test significance
lmI1 <- localmoran(rome$avg_price_, dnbneigh.listw)
print(head(lmI1))
lmI2 <- localmoran(rome$avg_price_, dnb2.listw)
print(head(lmI2))
lmI3 <- localmoran(rome$avg_price_, listW2s)
print(head(lmI3))

