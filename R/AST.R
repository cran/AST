

AST <- function(data.residual, spaceMatrix, par.time=0.5,par.age=1,
                weight.coverage=0.9, agecat, minyear,maxyear){

################### err #####

  if(   !is.numeric(c(weight.coverage,par.time,par.age,agecat,minyear,maxyear))     ) stop("  weight.coverage,par.time,par.age,agecat,minyear,maxyear must be number.")
  if(   prod(c(weight.coverage,par.time,par.age,agecat,minyear,maxyear))<0     ) stop("weight.coverage,par.time,par.age,agecat,minyear,maxyear must be positive numbers")
  if( weight.coverage<0 | weight.coverage>1) stop("weight.coverage must be a number between 0 and 1")
  if(maxyear<minyear) stop("min year must be lower than max year.")
  if(!is.matrix(spaceMatrix)) stop("spaceMatrix must be a matrix format.")
  if(!is.data.frame(data.residual)) stop("data.residual must be a data.frame format.")
  if( !all( c("age","year","location","residual")%in%colnames(data.residual) )  ) stop("data.residual must contain these variable names: age, year , location, residual  ")


  if( any(c(  !is.numeric(spaceMatrix) , min(spaceMatrix)<0 , max(spaceMatrix)>1   ))  ) stop("spaceMatrix must be between 0 and 1 ")
  if( is.null(rownames(spaceMatrix)  )  ) stop("row names of spaceMatrix is necessary!")
  if(   any(  is.na(
    as.numeric(rownames(spaceMatrix))
  )
  )
  )stop("row names of spaceMatrix must be matched with location(number)")


  if(  c("coverage")%in%colnames(data.residual)   ){
   t <- data.residual$coverage
   if( any(c(  !is.numeric(t)   ,   length(unique(t))!=2 , max(t)!=1 , min(t)!=0  ))   ) stop("coverage must be a binary variable ")
   data.residual <- data.residual[!is.na(data.residual$coverage),]
  }else{
   data.residual$coverage <- 1
   weight.coverage = 1
  }

 #######################
    T=(maxyear-minyear)+1
	LOC = as.numeric(rownames(spaceMatrix)	 )
	data.pred = expand.grid(year=minyear:maxyear,age=1:agecat,location=LOC)

	data.residual <- data.residual[!is.na(data.residual$year),]
	data.residual <- data.residual[!is.na(data.residual$age),]
	data.residual <- data.residual[!is.na(data.residual$location),]
	data.residual <- data.residual[!is.na(data.residual$residual),]

	data.residual <- data.residual[data.residual$age%in%c(1:agecat) & data.residual$location%in%LOC & data.residual$year%in%c(minyear:maxyear) , ]

########################### fun ###########

	  ############age weight
calcage = function (ag)
{
	z = matrix(NA , nrow=ag , ncol=ag)
	z[1:ag , 1:ag] = matrix(NA , nrow=ag , ncol=ag)
	A <- as.numeric(rownames(z) <- c(1:ag))
	B <- (colnames(z) <- c(1:ag))

	for (i in 1:ag)
		{
		for (j in 1:ag)
			{
			z[i,j]=1/exp(par.age*(abs(A[i]-B[j])))
			}
		}
return(z)
}


###############time_weight
calctime = function (T ,minyear,maxyear )
{
y = matrix(NA,nrow=T,ncol=T)
y[1:T,1:T] = matrix(NA,nrow=T,ncol=T)
A <- as.numeric(rownames(y) <- c(minyear:maxyear))
B <- (colnames(y) <- c(minyear:maxyear))

for (i in 1:T)
		{
		for (j in 1:T)
			{
			y[i,j]=(1-(abs((A[i])-(B[j]))/T)^par.time)^3
			}
		}

return (y)
}


###############final weight and rescaling

calcW = function(x,y,z)
{
	zy = kronecker(z, y)
	weight = apply(zy,1,sum,na.rm=T)
		for (i in 1:nrow(zy))
			{
			zy[i,]= zy[i,]/weight[i]
			}
	W = kronecker(x, zy)
return(W)
}

###################weight matrix colnames

calcmatchnames = function(W , Data)
{
	colnames(W) = rep (1:ncol(W))
	out <- c()
	for (i in LOC) for (j in (1:agecat)) for(k in minyear:maxyear) out= c(out,(paste(k,"-",j,"-",i, sep="")))
	colnames(W) = out
	Data <- Data[order(Data$location, Data$age, Data$year), ]
	Data$ID = paste0(Data$year,"-",Data$age,"-",Data$location)
	W_Data <-  as.matrix( W[, as.character(Data$ID)])
return(W_Data)
}

########################coverage weight
calccove = function(W_Data , Data)
{
	data <- Data[order(Data$location, Data$age, Data$year), ]
	data$ID = paste0(data$year,data$age,data$location)
	data$N = 1:nrow(data)
	d = split(data$N, as.factor(data$coverage) )
	d1 = as.vector(d$'1')
	d2 = as.vector(d$'0')
	colnames(W_Data) = rep (1:ncol(W_Data))
	W_Data[,d1] = apply(W_Data[,d1], 2, function(x) (x*weight.coverage))
	W_Data[,d2] = apply(W_Data[,d2], 2, function(x) (x*(1-weight.coverage)))
	w_all = W_Data
return(w_all)

}
########################final rescale
finalw = function(w_all)
{
	weight = apply(w_all,1,sum,na.rm=T)
		for (i in 1:nrow(w_all))
			{
			w_all[i,]= w_all[i,]/weight[i]
			}
	fw = w_all
return(fw)
}

#################residual and weight matrix
calcrespred = function (fw , Data)
{
	Data <- Data[order(Data$location, Data$age, Data$year), ]
	Data$ID = paste0(Data$year,"-",Data$age,"-",Data$location)
	resvec = Data[ , c("residual")]
	mat.vec <- c(fw %*% resvec)
		rownames(fw) = rep (1:nrow(fw))
		out.row <- c()
		for (i in LOC) for (j in (1:agecat)) for(k in minyear:maxyear) out.row= c(out.row,(paste(k,"-",j,"-",i, sep="")))
		rownames(fw) = out.row
	mat = cbind(mat.vec , out.row)
	mat = as.data.frame(mat)
return(mat)
}
####################
calcpred = function (mat , pData)
{
	pData$ID = paste0(pData$year,"-",pData$age,"-",pData$location)
	mat$ID = as.character(mat$ID)
	colnames(mat)[2] <- "ID"
	outdata= merge(mat, pData, by.x= "ID" , by.y="ID" )
return(outdata)
}

################# run #######

  timeMat = calctime(T , minyear , maxyear)
  ageMat = calcage(agecat)

   final_weight = calcW(spaceMatrix ,timeMat , ageMat)

# W=final_weight ; Data = data.residual
   names = calcmatchnames(final_weight ,data.residual)
    rm(final_weight)
  coverageW = calccove(names,data.residual)
    rm(names)
  fw = finalw(coverageW)
    rm(coverageW)
  wmat = calcrespred(fw , data.residual)
    rm(fw)
  colnames(wmat) <- c("residual_AST","ID")
  out = calcpred(wmat , data.pred)
  out$residual_AST = as.numeric( as.character(out$residual_AST)  )
  out = out[order(out$location,out$year,out$age),]

  OUT = list()
  OUT$adj.res <- out
  OUT$Age_weight <- ageMat
  OUT$time_weight <- timeMat

  return(OUT)

}
