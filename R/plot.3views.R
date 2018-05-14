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

#col.pat.table, example
#[,1]   [,2] [,3]    
#[1,] "pop1" "x"  "yellow"
#[2,] "pop2" "u"  "blue" 

plot.3views <- function(X,labels,col.pat.table=NA,only.row=NA){  
  map.color = c("red",rgb(0,68,27,max=255),"blue",rgb(231,41,138,max=255),"darkorange","black")
  map.color = c(map.color,rgb(102,37,6,max=255),rgb(63,0,125,max=255),"green")
  map.color = c(map.color,"cyan",rgb(250,159,181,max=255),"yellow","darkgrey")
  map.color = c(map.color,rgb(116,196,118,max=255))
  
  map.pch = c(1,0,2:18,35:38,60:64,94,126)
  map.pch = c(map.pch,33:34,42,45,47,40,91,123,41,92,93,125)
  map.pch = c(map.pch,49:57,97:107,109:110,112:119,121:122)
  map.pch = c(map.pch,65:78,81:82,84:85,89)

  map.pattern = c()
  for (i in 1:length(map.pch))
    for (j in 1:length(map.color)){
      tmp = c(i,j)
      map.pattern = rbind(map.pattern,tmp)
    }

  if (class(labels) == "data.frame"){
    labels = labels[,1]
    u.label = sort(unique(labels))
  }else{
    u.label = sort(unique(labels))
  }
  
  par(mfrow=c(2,2))
  
  #Top-Left
  par(mar=c(4, 1, 1, 2))
  plot(c(min(X[,1]),max(X[,1])),c(min(X[,2]),max(X[,2])),type="n",xlab="",ylab="",main="",axes=FALSE)
  axis(side=1, labels=TRUE, line=2)
  axis(side=4, labels=TRUE, line=2)
  mtext("PC1", side=1, line=0.5)
  mtext("PC2", side=4, line=0.5)
  set_legend = NULL
  set_pch = NULL
  set_col = NULL
  for (k in 1:length(u.label)){
    if (anyNA(col.pat.table)){
      spch = map.pch[map.pattern[k,1]]
      scolor = map.color[map.pattern[k,2]]
    }else{
      idx = which(col.pat.table[,1] == u.label[k])
      spch = col.pat.table[idx,2]
      scolor = col.pat.table[idx,3]
    }
    
    if (anyNA(only.row)){
      points(X[labels %in% u.label[k],1],X[labels %in% u.label[k],2],col=scolor,pch=spch)
    }else{
      tmp.idx1=which(labels %in% u.label[k])
      tmp.idx2=intersect(tmp.idx1,only.row)
      if (length(tmp.idx2)>0){
        points(X[tmp.idx2,1],X[tmp.idx2,2],col=scolor,pch=spch)
      }
    }
    set_pch = c(set_pch,spch)
    set_col = c(set_col, scolor)
  }
  
  #Top-Right
  par(mar=c(4, 3.5, 1, 1))
  plot(c(min(X[,3]),max(X[,3])),c(min(X[,2]),max(X[,2])),type="n",xlab="",ylab="",main="",axes=FALSE)
  axis(side=1, labels=TRUE, line=2)
  mtext("PC3", side=1, line=0.5)
  for (k in 1:length(u.label)){
    if (anyNA(col.pat.table)){
      spch = map.pch[map.pattern[k,1]]
      scolor = map.color[map.pattern[k,2]]
    }else{
      idx = which(col.pat.table[,1] == u.label[k])
      spch = col.pat.table[idx,2]
      scolor = col.pat.table[idx,3]
    }
    if (anyNA(only.row)){
      points(X[labels %in% u.label[k],3],X[labels %in% u.label[k],2],col=scolor,pch=spch)
    }else{
      tmp.idx1=which(labels %in% u.label[k])
      tmp.idx2=intersect(tmp.idx1,only.row)
      if (length(tmp.idx2)>0){
        points(X[tmp.idx2,3],X[tmp.idx2,2],col=scolor,pch=spch)
      }
    }
  }
  
  #Bottom-Left
  par(mar=c(1, 1, 0.5, 2))
  plot(c(min(X[,1]),max(X[,1])),c(min(X[,3]),max(X[,3])),type="n",xlab="",ylab="",main="",axes=FALSE)
  axis(side=4, labels=TRUE, line=2)
  mtext("PC3", side=4, line=0.5)
  for (k in 1:length(u.label)){
    if (anyNA(col.pat.table)){
      spch = map.pch[map.pattern[k,1]]
      scolor = map.color[map.pattern[k,2]]
    }else{
      idx = which(col.pat.table[,1] == u.label[k])
      spch = col.pat.table[idx,2]
      scolor = col.pat.table[idx,3]
    }
    
    if (anyNA(only.row)){
      points(X[labels %in% u.label[k],1],X[labels %in% u.label[k],3],col=scolor,pch=spch)
    }else{
      tmp.idx1=which(labels %in% u.label[k])
      tmp.idx2=intersect(tmp.idx1,only.row)
      if (length(tmp.idx2)>0){
        points(X[tmp.idx2,1],X[tmp.idx2,3],col=scolor,pch=spch)
      }
    }
  }
  
  #Bottom-Right
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  if (length(u.label)>54){
    legend('center', inset=0, legend=u.label, pch=set_pch, col=set_col, ncol=4)
  }else if (length(u.label)>36){
    legend('center', inset=0, legend=u.label, pch=set_pch, col=set_col, ncol=3)
  }else if (length(u.label)>18){
    legend('center', inset=0, legend=u.label, pch=set_pch, col=set_col, ncol=2)
  }else{
    legend('center', inset=0, legend=u.label, pch=set_pch, col=set_col, ncol=1)
  } 
}

#--include--end--#
