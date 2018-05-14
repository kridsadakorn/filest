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

#raw.data is a SNP matrix while rows represent individuals and columns represent SNP
#idx.p1 is a vector containing the indexes of population 1 in raw.data
#idx.p2 is a vector containing the indexes of population 2 in raw.data
#return a global Fst as a single value
fst.hudson <-function(raw.data,idx.p1,idx.p2){
  prestep.fst.one.marker <- function(alleles,idx.p1,idx.p2){
    
    #Pop 1
    G = alleles[idx.p1]
    no.AA=length(which(G==0))
    no.AB=length(which(G==1))
    no.BB=length(which(G==2))
    n1=(no.AA+no.AB+no.BB)*2
    p.A=(no.AA*2 + no.AB)/n1
    #p.B=(no.BB*2 + no.AB)/n1
    #p1 = min(p.A,p.B,na.rm = T)
    p1 = p.A

    #Pop 2
    G = alleles[idx.p2]
    no.AA=length(which(G==0))
    no.AB=length(which(G==1))
    no.BB=length(which(G==2))
    n2=(no.AA+no.AB+no.BB)*2
    p.A=(no.AA*2 + no.AB)/n2
    #p.B=(no.BB*2 + no.AB)/n2
    #p2 = min(p.A,p.B,na.rm = T)
    p2 = p.A

    
    N = (p1 - p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
    D = p1*(1-p2) + p2*(1-p1)
    
    return(c(N,D))
  }
  
  set.fst = apply(raw.data,2,prestep.fst.one.marker,idx.p1=idx.p1, idx.p2=idx.p2)
  #fst = mean(set.fst[1,],na.rm=T) / mean(set.fst[2,],na.rm=T)
  fst = mean(set.fst[1,]/set.fst[2,],na.rm = T)

  return(fst)
}

#--include--end--#

#24/01/2016
#raw.data is a SNP matrix while rows represent individuals and columns represent SNP
#idx.p1 is a vector containing the indexes of population 1 in raw.data
#idx.p2 is a vector containing the indexes of population 2 in raw.data
#return a vector of Fst
fst.each.snp.hudson <-function(raw.data,idx.p1,idx.p2){
  prestep.fst.one.marker <- function(alleles,idx.p1,idx.p2){
    
    #Pop 1
    G = alleles[idx.p1]
    no.AA=length(which(G==0))
    no.AB=length(which(G==1))
    no.BB=length(which(G==2))
    n1=(no.AA+no.AB+no.BB)*2
    p.A=(no.AA*2 + no.AB)/n1
    #p.B=(no.BB*2 + no.AB)/n1
    #p1 = min(p.A,p.B,na.rm = T)
    p1 = p.A
    
    #Pop 2
    G = alleles[idx.p2]
    no.AA=length(which(G==0))
    no.AB=length(which(G==1))
    no.BB=length(which(G==2))
    n2=(no.AA+no.AB+no.BB)*2
    p.A=(no.AA*2 + no.AB)/n2
    #p.B=(no.BB*2 + no.AB)/n2
    #p2 = min(p.A,p.B,na.rm = T)
    p2 = p.A
    
    
    N = (p1 - p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
    D = p1*(1-p2) + p2*(1-p1)
    
    return(c(N,D))
  }
  
  set.fst = apply(raw.data,2,prestep.fst.one.marker,idx.p1=idx.p1, idx.p2=idx.p2)
  #fst = mean(set.fst[1,],na.rm=T) / mean(set.fst[2,],na.rm=T)
  vec.fst = set.fst[1,]/set.fst[2,]
  
  return(vec.fst)
}

