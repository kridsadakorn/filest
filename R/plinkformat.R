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

#For detail about BED file format, check http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

# Return object as
# object$snp - SNP matrix from bed file
# object$snp.info - SNP information from bim file
# object$ind.info - individual information from fam file

#Example:
#bed = "test.bed"
#bim = "test.bim"
#fam = "test.fam"
#snp = read.bed(bed,bim,fam)

read.bed <- function(bed,bim,fam,only.snp=FALSE){
  ret = NA
  
  if (!file.exists(bed)){
    cat(paste0("Error: BED file doesn't exist: ",bed,"\n"))
    return(ret)
  }
  if (!file.exists(bim)){
    cat(paste0("Error: BIM file doesn't exist: ",bim,"\n"))
    return(ret)
  }
  if (!file.exists(fam)){
    cat(paste0("Error: FAM file doesn't exist: ",fam,"\n"))
    return(ret)
  }
  #Read BIM file
  snp.info = read.table(bim,header=FALSE)
  colnames(snp.info) = c("chr","ID","GD","position","allele1","allele2")
  
  #Read FAM file
  ind.info = read.table(fam,header=FALSE)
  colnames(ind.info) = c("FamID","IndID","PatID","MatID","sex","phenotype")
  
  no.ind = dim(ind.info)[1]
  no.snp = dim(snp.info)[1]
  
  if (only.snp == TRUE){
    snp.info = NA
    ind.info = NA
  }
  
  #Read BED file
  fh = file(bed, 'rb')
  #Read the first three bytes to check file format
  buff = readBin(fh, what="raw",n=3)
  if (sum(buff[1:2] == c('6c','1b')) != 2){
    cat(paste0("Error: BED file is not in a correct format: ",bed,"\n"))
    return(ret)
  }
  no.byte.to.read = NA
  no.loop = NA
  if (buff[3] == '01'){
    no.byte.to.read = ceiling(no.ind/4.0)
    no.loop = no.snp
  }else{
    no.byte.to.read = ceiling(no.snp/4.0)
    no.loop = no.ind
  }

  snp = readBin(fh, what="raw",n=(no.byte.to.read*no.loop))
  close(fh)
  
  length.snp = length(snp)
  no.part=10
  sub.length.snp = round(length.snp / no.part)
  
  tmp.snp = c()
  for (i in 1:no.part){
    if (i != no.part){
      idx1 = (i-1)*sub.length.snp + 1
      idx2 = i*sub.length.snp
      sub.snp =snp[idx1:idx2]
    }else{
      idx1 = (i-1)*sub.length.snp + 1
      idx2 = length.snp
      sub.snp =snp[idx1:idx2]
    }

    sub.snp.bits = matrix(as.integer(rawToBits(sub.snp)),ncol=2,byrow=T)
    
    tmp.sub.snp.bits = sub.snp.bits[,1]+sub.snp.bits[,2]
    tmp.sub.snp.bits[which((sub.snp.bits[,2]-sub.snp.bits[,1])<0)] = NA
    tmp.snp = c(tmp.snp,tmp.sub.snp.bits)
  }
  
  snp = matrix(tmp.snp,ncol=no.snp,byrow=F)[1:no.ind,]
  
  if (only.snp == FALSE){
    return(list("snp"=snp,"snp.info"=snp.info,"ind.info"=ind.info))
  }else{
    return(list("snp"=snp))
  }
}

# The variable "object" needs to be in this format:
# object$snp - SNP matrix from bed file
# object$snp.info - SNP information from bim file
# object$ind.info - individual information from fam file
# The variable "file" is a prefix for output files
write.bed <- function(object,file){
  
  bed = paste0(file,'.bed')
  bim = paste0(file,'.bim')
  fam = paste0(file,'.fam')
  if (file.exists(bed)){
    cat(paste0("Overwrite the existed bed file\n"))
  }
  if (file.exists(bim)){
    cat(paste0("Overwrite the existed bim file\n"))
  }
  if (file.exists(fam)){
    cat(paste0("Overwrite the existed fam file\n"))
  }
  
  #Write BIM file
  if (is.null(object$snp.info)){
    cat(paste0("Couldn't find object$snp.info: a bim file was not created\n"))
  }else{
    write.table(object$snp.info,file=bim,quote=F,row.names=F,col.names=F)
  }
  
  #Write FAM file
  if (is.null(object$ind.info)){
    cat(paste0("Couldn't find object$ind.info: a fam file was not created\n"))
  }else{
    write.table(object$ind.info,file=fam,quote=F,row.names=F,col.names=F)
  }
  
  #Header of BED file
  header = c(0,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,1,0,0,0,0,0,0,0)
  
  no.ind = dim(object$ind.info)[1]
  no.snp = dim(object$snp.info)[1]
  
  fill.up = 4-(no.ind %% 4)
  tmp = NA
  if (fill.up != 4){
    tmp = matrix(rep(0,no.snp*fill.up),ncol=no.snp)
  }
  
  vec = object$snp
  if (!anyNA(tmp)){
    vec = rbind(vec,tmp)
  }
  
  vec = as.vector(vec)
  buff = matrix(rep(0,length(vec)*2),ncol=2)
  buff[which(vec == 1),1] = 0
  buff[which(vec == 1),2] = 1
  buff[which(vec == 2),1] = 1
  buff[which(vec == 2),2] = 1
  buff[which(is.na(vec)),1] = 1
  buff[which(is.na(vec)),2] = 0
  vec = NULL
  buff = c(header,as.vector(t(buff)))
  buff=packBits(as.raw(buff))
  
  fh = file(bed, 'wb')
  writeBin(buff,fh)
  close(fh)
}

#--include--end--#


