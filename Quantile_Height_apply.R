#Lectura de datos

setwd("~/Geo/Maestría/Datos_Erick/Tablas")
#With all point

library(foreign)
b <- read.dbf("int_BvFRP10p_MergeAllP.dbf")
head(b,3)

b <- write.csv(b, file = "BufferALL.csv")
bb <- read.csv("BufferALL.csv")
b <- bb[,-(3:27)]
bd <- data.frame(bb)


# Split buffers. This give me the clusters of points with points. 
s<- split(bd,bd$FID_Buffer)


#For -- Adding labels to each histogram
ind = bd[,'FID_Buffer']
ni <- unique(ind)
Ni = length(unique(ind))
ind1 = unique(ind)
indi <- sort(ind1,decreasing = F)  # Clusters with data


#Creating list of clusters, that will be useful locating clusters in file fil_Dat.csv
clus1 <- paste(c("X00"), indi[1:10], sep="")
clus2 <- paste(c("X0"), indi[11:100], sep="")
clus3 <- paste(c("X"), indi[101:675], sep="")
clus <- c(clus1, clus2, clus3)


#creating wave length Matrix for spectral firms 
datos <- read.csv("Ref_fil2.csv", header = TRUE, sep = ',')

# substacting two first rows (especie, file)
datos <- datos[-(2154),]
datos <- datos[-(1:2),]
datos <- datos[,-(775)]

write.csv(datos, file = "datos.csv")


#Reading Data line 47
dat <- read.csv("datos.csv", header = T, sep ="," )
dat <- data.frame(dat)

#Creating Matrix for implement PLS (****RIGHT HAND***)
# Selec only the clusters with data

library(dplyr)
data <- tbl_df(dat)
matrixpls <- select(dat, one_of(clus))

#Transpose. 
tmatrixpls <- t(matrixpls)
colnames(tmatrixpls) <- rownames(matrixpls)       #  I Assing the matrix colnames 

# I need to create unique vectors with height data from each cluster to create 
# Quantiles prob = seq(0.1,1,0.1)
hR <- lapply(s, function(s) (unique(sort(s[,'Height']))))
QuantilehR <- lapply(hR, function(hR) quantile(hR, probs = seq(0.1,1,0.1)))

QM <- data.frame(QuantilehR)
head(QM,3)


#Renaming names from QM to X00n, X0nn, Xnnn  

QM2<- sub("X","",names(QM))
cl1 <- paste(c("00"),QM2[1:10], sep="")
cl2 <- paste(c("0"), QM2[11:100], sep="")
cl3 <- QM2[101:675]
cl  <- c(cl1, cl2, cl3)
clu <- paste(c("X"),cl, sep = "")

names(QM) <- clu

#Transpose the matrix QM. dim(QM) =  675 x 10 
TQM <- t(QM)

#tmatrixpls = wave length matrix
library(pls)

#pls.options(parallel = 1) # Use mclapply with 4 CPUs
t <- proc.time() # Inicia el cronómetro
#library(parallel) # for the makeCluster call


height.pls <- plsr(TQM ~ tmatrixpls, 65, validation = "CV", dframe = T, method = "simpls")
#t <- proc.time() # Inicia el cronómetro
#height.pls <- plsr(tmatrixpls ~ TQM, 6, validation = "LOO", dframe = T)
proc.time()-t    # Detiene el cronómetro
coefplot(height.pls, ncomp = 1:65)
biplot(height.pls)
summary(height.pls)

Xscores<- scores(height.pls)
Xloadings <- loadings(height.pls)
Yscores <- Yscores(height.pls) 
Yloadings <- Yloadings(height.pls)

diag(cov(Xscores)) - apply(Xscores,2,var)

# components arranged in decreasing order
Yl <- apply(Yloadings, 1, function(x) order(x, decreasing = T))
# Yloadings arranged in decreasing order
Yls <- apply(Yloadings, 1, function(x) sort(x, decreasing =T))
# Ploting Yloadings
#pYlsa <- apply(t(Ylsa), 1, function(x) plot(x, main = paste("Var")))
for (i in 1:dim(Yloadings)[1]){plot(Yls[,i], main=paste("Var",i), xlab = "Components")}

#Creating the model Y = Xs*t(Yl) + med(Y)
xs <- matrix(NA, dim(Xscores)[1],dim(Xscores)[2])
Xs <- replicate(dim(Yloadings)[1],xs)
y0 <- matrix(NA, dim(Xscores)[1] ,dim(Yl)[2])

for (i in 1: dim(Yloadings)[1] ) 
{
  #Yloadings Transpose
  Xs[,,i] <- Xscores[,c(Yl[,i])]
  y0[,i] <- Xs[,,i]%*%as.matrix(Yls[,i])
}

#Mean Matrix from Y
mYls <- apply(TQM,2,mean)
M.Yls <- (matrix(1,dim(Xscores)[1],1))%*%mYls

#Model Y = Xs*t(Yl) + med(Y), that give us the Y estimated

Yest <- y0+M.Yls

#Variance
Xs.var <- var(Xscores[,1])

y <- sort((Yls[,2])^2, decreasing = T)
plot(Xs.var*(cumsum(y)))

for (i in 1: 10)
  {
   plot(Xs.var*cumsum( sort((Yls[,i])^2, decreasing=T) ), main = paste("Var",i), ylab = "CumSum", xlab = "Components")
  }

