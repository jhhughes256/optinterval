#R SCRIPT TO FIND COEFFICIENTS FOR CUBIC SPLINE INTERPOLATION
#INPUT is a dataframe in a NONMEM style format


#remove all current objects in the workspace
 rm(list=ls(all=TRUE))

#Load the required libraries
 library(splines)
 library(doBy)

#Load the file and set working directory (may be Windows specific!)
  file.name.ext <- "cubicdata.csv"
  file.name <- gsub(".csv", "", file.name.ext)
  file.name.out <- paste(file.name,".int.csv",sep="")


#Read *.fit file and attach, so column names are available
  orgdata <- read.csv(file=file.name.ext, stringsAsFactors=F, na.string=".", as.is=T)
  names(orgdata)

#Rename the first column, as C was used as a comment flag
  names(orgdata)[1] <- "ID"


#-------------------------------------------------------------------------------
#Example interpolation of the spline for the first subject
  orgdata1 <- orgdata[orgdata$ID==1,]
  spline1 <- interpSpline(COBS ~ TIME, orgdata1)      #piece-wise polynomial spline is default

  #Output is a spline object with the coefficients of the interpolation for each time point (aka knot)
  spline1
  spline1$knots  #The x-values
  spline1$coefficients   #The 4 coefficients of a cubic spline for each knot
  spline1all <- cbind(spline1$knots,spline1$coefficients)  #Joint them together if required
  
  #Important to plot the spline!
   plot( orgdata1$COBS ~ orgdata1$TIME, ylim=c(0,6) )
   splinedata <- predict(spline1, seq(0,80,length.out = 800))
   points(splinedata, type = "l", col="red")
   #Note that different assumptions can be made to control the behaviour of the spline at either end - the boundary conditions



#-------------------------------------------------------------------------------
#Example of application of interpSpline to all subjects in a NONMEM style dataframe
#Requires the doBy package - very handy in this line of work!
  splineall <- lapplyBy(~ID, data=orgdata,   function(d) interpSpline(d$COBS ~ d$TIME))
  #Output is a list of spline objects
  splineall

  #How to access the items of the spline list
  #splineall_1 <- splineall$'1'
  splineall_1 <- splineall[[1]]
  names(splineall_1)
  splineall_1$knots
  splineall_1$coefficients

  
#-------------------------------------------------------------------------------
#Example of application of interpSpline to all subjects in a NONMEM style dataframe
#Adding plotting of spline subject by subject
#Plots to screen - plotting to a file may be more useful
  
  #Define a wrapper function for interpSpline that adds plotting
  interpSplineplot <- function(ydata,xdata)
   {
    ymin <- min(ydata)*0.8
    ymax <- max(ydata)*1.2
    xmin <- min(xdata)
    xmax <- max(xdata)
    splineobj <- interpSpline(ydata ~ xdata)
    plot(ydata ~ xdata, ylim=c(ymin,ymax))
    splinedata <-  predict(splineobj, seq(xmin,xmax,length.out=length(xdata)*10))
    points(splinedata, type = "l", col="red")
    return(splineobj)
   }
  
  splineall <- lapplyBy(~ID, data=orgdata,   function(d) interpSplineplot(d$COBS,d$TIME))
  #Output remains a list of spline objects
  splineall


#-------------------------------------------------------------------------------
#Return a dataframe with the exported coefficients
#Turn the results back into a dataframe:

 #A utility function to merge a list of spline objects into a dataframe
 merge.spline.list <- function(x)
  {
  #Merges a list of spline objects into one big dataframe
  alldata <- NULL
  for (i in 1:length(x))
    {
    temp <- x[[i]]
    temp$knots <- format(temp$knots, digits=2)          #NONMEM doesn't like too many digits
    temp$coefficients <- format(temp$coefficients, digits=2) #NONMEM doesn't like too many digits
    temp <- cbind(temp$knots,temp$coefficients)
    alldata <- rbind(alldata, temp)
    alldata <- as.data.frame(alldata, stringAsFactors=F)
    }
  names(alldata) <-c("KNOT","COF1","COF2","COF3","COF4")
  return(alldata)
  }

#Apply the function
splinedf <- merge.spline.list(splineall)

#Make a final dataframe for export to NONMEM
   output <- cbind(orgdata,splinedf)
  
  #Change NA back to "."
   output[is.na(output)==T] <- "."

  #Change ID to CID
   names(output)[1] <- "CID"
   
  #Add an empty DV column for simulation
   output$DV <- "."

  #Write the file
  write.csv(output, file=file.name.out, row.names=F, quote=F)
