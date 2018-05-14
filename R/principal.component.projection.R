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
cal.PC.projection <- function(data, id, status, pop, no.pc = 10){
  #library(rARPACK)


  X.con = data[which(status==1),]
  X.case = data[which(status==2),]
  new.pop = pop[which(status==1)]
  tmp.pop = pop[which(status==2)]
  new.pop = c(new.pop,tmp.pop)
  new.id = id[which(status==1)]
  tmp.id = id[which(status==2)]
  new.id = c(new.id,tmp.id)

  X = X.con

  XPi = (colSums(X)+1)/(2+(2*dim(X)[1]))
  XSD = (XPi * (1-XPi))^0.5

  XCM = colMeans(X)
  A = t((t(X)-XCM)/XSD)
  A = A[,colSums(is.na(A))<nrow(A)]

  AA=A %*% t(A)
  evv = eigs_sym(AA, k=no.pc)
  eigen.value = evv$values
  PCs = evv$vectors
  B.T=diag(1/sqrt(eigen.value)) %*% t(PCs) %*% A

  PCcon=A %*% t(B.T)

  X = X.case

  XPi = (colSums(X)+1)/(2+(2*dim(X)[1]))
  XSD = (XPi * (1-XPi))^0.5

  XCM = colMeans(X)
  A = t((t(X)-XCM)/XSD)
  A = A[,colSums(is.na(A))<nrow(A)]

  PCcase=A %*% t(B.T)

  label.status=rep('control',dim(PCcon)[1])
  tmp.label=rep('case',dim(PCcase)[1])
  label.status = c(label.status,tmp.label)

  PC_project = rbind(PCcon,PCcase)

  return(list("PC"=PC_project,"id"=new.id,"label"=new.pop,"status"=label.status))
}
#--include--end--#
