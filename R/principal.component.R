##########################################################################
## IPCAPS Library
## Author: Kridsadakorn Chaichoompu (kridsadakorn.chaichoompu@ulg.ac.be)
## Date: 01/06/2016
## Description:
##    This code is a part of Iterative Pruning to CApture Population 
##    Structure (IPCAPS) Library
##    
##Licence: GPL V3
## 
##    Copyright (C) 2016  Kridsadakorn Chaichoompu
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#--include--begin--#

###
# Parameters
# X         Data matrix which rows represent samples and columns represent features.
# PCscore  To determine whether scaled principal components will be calculated, 
#           otherwise the eigen vectors will be used instead. Default value is TRUE.
# no.pc     A number  of principal components (PCs) to be calculated. Default value is NA,
#           it means that all PCs will be calculated.
# data.type To determine the type of data. It can be "linear" (default) and "snp".
# weight    To determine whether weights will be used to solve correlation between
#           features.
#
# Examples:
#
# PC = cal.PC.linear(X)
# plot(PC$PC)
# plot(PC$evalue)
#
# PC = cal.PC.linear(X,PCscore = FALSE)
# plot(PC$PC)
# plot(PC$evalue)
#
# PC = cal.PC.linear(X,no.pc=3)
# plot(PC$PC)
# plot(PC$evalue)
#
# PC = cal.PC.linear(X, data.type="snp")
# plot(PC$PC)
# plot(PC$evalue)


cal.PC.linear <- function(X, PCscore = TRUE, no.pc = NA, data.type="linear", XXT=TRUE){
  
  replace.missing <- function(X,missing=NA,rep){
    if (is.na(missing)){
      idx = which(is.na(X))
      X[idx] = rep[idx]
    }else{
      idx = which(X == missing)
      X[idx] = rep[idx]
    }
    return(X)
  }
  

  #Resolve missing value by median
  #cat(paste0("Missing values are set as medians\n"))
  X.median = apply(X,2,median,na.rm=TRUE)
  missing.char=NA
  X = t(apply(X,1,replace.missing,missing=missing.char,rep=X.median))
  
  if (data.type == "snp"){
    XPi = (colSums(X)+1)/(2+(2*dim(X)[1]))
    XSD = (XPi * (1-XPi))^0.5
  }else{
    #Assume linear
    XSD = apply(X,2,sd)
  }
  
  XCM = colMeans(X)
  A = t((t(X)-XCM)/XSD)
  A = A[,colSums(is.na(A))<nrow(A)]
  
  if (is.na(no.pc)){
    no.pc = dim(A)[1]
  }
  
  if (XXT == TRUE){
    AA=A %*% t(A)
    UDV = svds(AA, k=no.pc)
    eigen.value = UDV$d
    PCs = UDV$u
  }else{
    UDV = svds(A, k=no.pc)
    eigen.value = UDV$d
    PCs = UDV$u
  }
  
  
  if (PCscore == TRUE){
    #Calculate the PCscore and the scaled-up PCs
    B.T=diag(1/sqrt(eigen.value)) %*% t(PCs) %*% A
    PCs=A %*% t(B.T)
  }
  
  return(list("PC"=PCs,"evalue"=eigen.value))
}

#--include--end--#
