

#' (Internal function) Calculate the first margin for random numbers, internally used for
#' parallelization
#'
#' @param prob A number to specify probability
#' @param fst A number to specify Fst
#'
#' @return The first margin for random numbers

cal.margin1 = function(prob,fst){
  ret = prob*(1-fst)/fst
  return(ret)
}

#' (Internal function) Calculate the second margin for random numbers, internally used for
#' parallelization
#'
#' @param prob A number to specify probability
#' @param fst A number to specify Fst
#'
#' @return The second margin for random numbers

cal.margin2 = function(prob,fst){
  ret = (1-prob)*(1-fst)/fst
  return(ret)
}

#' (Internal function) Get random number based on beta distribution, internally used for
#' parallelization
#'
#' @param mar1 A number to specify the first margin for random function of beta
#' distribution.
#' @param mar2 A number to specify the second margin for random function of beta
#' distribution.
#'
#' @return A random number based on beta distribution
#'
#' @importFrom stats rbeta

get.random.beta = function(mar1,mar2){
  ret = rbeta(1, mar1, mar2, ncp = 0)
  return(ret)
}

#' (Internal function) Calculate probability for wild-type allele, internally used for
#' parallelization
#'
#' @param prob A number to specify probability
#'
#' @return A probability of wild-type allele

cal.prob.AA = function(prob){
  ret = (1-prob)^2
  return(ret)
}

#' (Internal function) Calculate probability for heterozygous allele, internally used for
#' parallelization
#'
#' @param prob A number to specify probability
#' @param riskratio A number to specify risk ratio for cases
#'
#' @return A probability of heterozygous allele

cal.prob.AB = function(prob,riskratio=1.0){
  ret = 2*prob*(1-prob)*riskratio
  return(ret)
}

#' (Internal function) Calculate probability for mutant-type allele, internally used for
#' parallelization
#'
#' @param prob A number to specify probability
#' @param riskratio A number to specify risk ratio for cases
#'
#' @return A probability of mutant-type allele

cal.prob.BB = function(prob,riskratio=1.0){
  ret = prob^2 * riskratio^2
  return(ret)
}

#' (Internal function) Sample SNP according to given probabilities, internally used for
#' parallelization
#'
#' @param p.AA A number to specify probability of wild-type allele
#' @param p.AB A number to specify probability of heterozygous allele
#' @param p.BB A number to specify probability of mutant-type allele
#' @param no.ind A number to specify the amount of individuals to be generated
#'
#' @return A simulated SNP in additive coding

do.sample.SNP = function(p.AA,p.AB,p.BB,no.ind){
  ret = sample(c(0,1,2),no.ind, prob=c(p.AA,p.AB,p.BB), replace=TRUE)
  return(ret)
}

#' (Internal function) Create categorical data, internally used for
#' parallelization
#'
#' @param marker A string to specify a marker
#' @param name A string to specify a name of marker
#'
#' @return A categorical value

do.create.cate = function(marker,name){
  ret = paste(name,marker,sep="_")
  return(ret)
}


#' (Internal function) Create outliers,internally used for
#' parallelization
#'
#' @param X An input matrix containing SNPs in additive coding
#' @param population A vector of labels for all individuals
#' @param outlier A number to specify the amount of outliers to be generated
#' @import rARPACK
#' @importFrom stats runif
#' @return The modified matrix with outliers

create.outlier = function(X,population,outlier){
  SVD=svd(X)
  umod=SVD$u
  idx = c()
  for (i in 1:length(population)){
    if (outlier[i]>0){
      if (i>1){
        start = sum(population[1:(i-1)])+1
        end = sum(population[1:i])
      }else{
        start = 1
        end = population[1]
      }
      idx = c(idx,start:end)
    }
  }

  no.outlier = length(idx)
  rn = runif(no.outlier,min = -0.1 ,max=0.1)
  umod[idx,1] = rn
  rn = runif(no.outlier,min = -0.2 ,max=0.2)
  umod[idx,2] = rn
  rn = runif(no.outlier,min = -0.4 ,max=0.4)
  umod[idx,3] = rn

  A = umod %*% diag(SVD$d) %*% t(SVD$v)
  A = round(A)
  A[which(A>2)] = 2
  A[which(A<0)] = 0

  return(A)
}

#' (Internal function) Calculate probabilities for categorical data, internally used for
#' parallelization
#'
#' @param prob A number to specify probability
#' @param max.category A number to specify the maximum number of categorical groups
#'
#' @return Probability for categorical data
#'
#' @importFrom stats runif

cal.prob.categorical.data = function(prob,max.category){

  pset = runif(max.category-1,min=0,max=1)
  pset = pset/(sum(pset))
  modifier = 1/(1-prob)
  pset=pset/modifier

  pset = c(prob,sort(pset))

  return(pset)
}

#' (Internal function) Sample categorical data,internally used for
#' parallelization
#'
#' @param prob A number to specify probability
#' @param max.category A number to specify the maximum number of categorical groups
#' @param no.ind A number to specify the amount of individuals to be generated
#'
#' @return A simulated data in categorical mode

do.sample.categorical.data = function(prob,max.category,no.ind){
  ret = sample(1:max.category,no.ind, prob=prob, replace=TRUE)
  return(ret)
}

#' (Internal function) Generate categorical data, internally used for
#' parallelization
#'
#' @param pop A number to specify probability
#' @param fst A number to specify Fst
#' @param no.marker A number to specify the amount of markers to be generated
#'
#' @return A marrix of categorical data
#'
#' @importFrom stats runif

generate.categorical.data =function(pop,fst,no.marker){
  ret = NULL
  if ((length(pop) !=  length(fst)) || (length(pop)<=0) || (length(fst)<=0) || (sum(pop<0)!=0) || (sum(fst<0)!=0) || no.marker<0){
    cat("Incorrect usage of generate.categorical.data()\n\te.g. generate.categorical.data(pop=c(50,50),fst=c(0.01,0.01),no.marker=100)\n")
    return(ret)
  }
  p = runif(no.marker, min = 0.1, max = 0.9)

  for (i in 1:length(pop)){
    mar1 = mapply(cal.margin1,prob=p,fst=fst[i]) #change here
    mar2 = mapply(cal.margin2,prob=p,fst=fst[i]) #change here
    p.b = mapply(get.random.beta,mar1=mar1,mar2=mar2)
    p.all = mapply(cal.prob.categorical.data,prob=p.b,max.category=10)
    data.mat = apply(p.all,2,do.sample.categorical.data,max.category=10,no.ind=pop[i])
    ret = rbind_bigmatrix(ret,data.mat)
  }
  return(ret)
}

#' (Internal function) Generate SNPs, internally used for
#' parallelization
#'
#' @param pop A number to specify the amount of populations
#' @param fst A number to specify Fst (representing genetic distance)
#' @param no.snp A number to specify the amount of SNPs
#' @param riskratio A number to specify risk ratio for cases
#'
#' @return Simulated SNPs
#'
#' @importFrom stats runif

generate.snp =function(pop,fst,no.snp,riskratio=1.0){
  ret = NULL
  if ((length(pop) !=  length(fst)) || (length(pop)<=0) || (length(fst)<=0) || (sum(pop<0)!=0) || (sum(fst<0)!=0) || no.snp<0){
    cat("Incorrect usage of generate.snp()\n\te.g. generate.snp(pop=c(50,50),fst=c(0.01,0.01),no.snp=1000)\n")
    return(ret)
  }
  p = runif(no.snp, min = 0.1, max = 0.9)

  for (i in 1:length(pop)){
    desired.fst = NA
    if (i<=2){
      desired.fst=fst[i]
    }else{
      desired.fst=2*fst[i]-fst[1]
    }
    mar1 = mapply(cal.margin1,prob=p,fst=desired.fst) #change here
    mar2 = mapply(cal.margin2,prob=p,fst=desired.fst) #change here
    p.b = mapply(get.random.beta,mar1=mar1,mar2=mar2)
    p.AA = mapply(cal.prob.AA,prob=p.b)
    p.AB = mapply(cal.prob.AB,prob=p.b,riskratio=riskratio)
    p.BB = mapply(cal.prob.BB,prob=p.b,riskratio=riskratio)
    data.mat = mapply(do.sample.SNP,p.AA=p.AA,p.AB=p.AB,p.BB=p.BB,no.ind=pop[i]) #change here
    ret = rbind_bigmatrix(ret,data.mat)
  }
  return(ret)
}


#' (Internal function) Generate labels, internally used for
#' parallelization
#'
#' @param pop A number to specify the amount of populations
#' @param outlier A number to specify the amount of outliers
#'
#' @return A vector of population labels

generate.label = function(pop,outlier){
  ret = NULL
  if ((length(pop)<=0) || (sum(pop<0)!=0)){
    cat("Incorrect usage of generate.label()\n\te.g. generate.label(pop=c(50,50))\n")
    return(ret)
  }

  for (i in 1:length(pop)){
    label = NULL
    if (outlier[i]>0){
      label = rep(paste0('outlier',i),times=pop[i]) #change here
    }else{
      label = rep(paste0('pop',i),times=pop[i]) #change here
    }
    ret = c(ret,label)
  }

  #ret=as.factor(ret)
  return(ret)
}



#' Combind two matrices by row for big data, internally used for
#' parallelization
#'
#' @param a The first matrix
#' @param b The second matrix
#'
#' @return The combined matrix by row
#'
#' @export
#' @seealso \code{\link{cbind_bigmatrix}}
#'
#' @examples
#' X <- matrix(c(1,2,0,1,2,2,1,2,0,0,1,2,1,2,2,2),ncol=4)
#' Y <- matrix(c(2,1,1,0,1,0,0,1,1,2,2,0,0,1,1,0),ncol=4)
#' Z <- rbind_bigmatrix(X,Y)
#' print(Z)

rbind_bigmatrix <- function(a,b){
  sizeA=dim(a)
  sizeB=dim(b)

  if (is.null(sizeA) && is.null(sizeB)){
    return(NULL)
  }else if (!is.null(sizeA) && is.null(sizeB)){
    return(a)
  }else if (is.null(sizeA) && !is.null(sizeB)){
    return(b)
  }else if (sizeA[2] != sizeB[2]){ #don't have a same number of column
    return(NULL)
  }
  else{ #both a and b are not null
    size=sizeB+sizeA
    ret = matrix(nrow=size[1],ncol=sizeA[2])
    ret[1:sizeA[1],] = a
    ret[(sizeA[1]+1):size[1],] = b
    return(ret)
  }
}

#' Combind two matrices by column for big data, internally used for
#' parallelization
#'
#' @param a The first matrix
#' @param b The second matrix
#'
#' @return The combined matrix by column
#'
#' @export
#' @seealso \code{\link{rbind_bigmatrix}}
#'
#' @examples
#' X <- matrix(c(1,2,0,1,2,2,1,2,0,0,1,2,1,2,2,2),ncol=4)
#' Y <- matrix(c(2,1,1,0,1,0,0,1,1,2,2,0,0,1,1,0),ncol=4)
#' Z <- cbind_bigmatrix(X,Y)
#' print(Z)
#'
cbind_bigmatrix <- function(a,b){
  sizeA=dim(a)
  sizeB=dim(b)

  if (is.null(sizeA) && is.null(sizeB)){
    return(NULL)
  }else if (!is.null(sizeA) && is.null(sizeB)){
    return(a)
  }else if (is.null(sizeA) && !is.null(sizeB)){
    return(b)
  }else if (sizeA[1] != sizeB[1]){ #don't have a same number of column
    return(NULL)
  }
  else{ #both a and b are not null
    size=sizeB+sizeA
    ret = matrix(nrow=sizeA[1],ncol=size[2])
    ret[,1:sizeA[2]] = a
    ret[,(sizeA[2]+1):size[2]] = b
    return(ret)
  }
}

#' (Internal function) Get all parameters
#'
#' @param param A string of parameters
#'
#' @return A list of parameters

get.para = function(param){
  name=NA
  population=NA
  fst=NA
  case=NA
  marker=NA
  replicate=NA
  riskratio=NA
  no.case.snp=NA
  missing=NA
  outlier=NA
  pc=NA
  fulloutput=NA
  cate=FALSE
  for (idx in 1:dim(param)[1]){
    if (param[idx,1]=='--setting'){
      name = as.character(param[idx,2])
    }
    if (param[idx,1]=='--population'){
      population = as.numeric(strsplit(as.character(param[idx,2]),split=",")[[1]])
      if (length(population[which(population<2)])>0){
        cat(paste0("Warning: replaced all numbers of population that are less than 2 by 2\n"))
        population[which(population<2)] = 2
      }
    }
    if (param[idx,1]=='--fst'){
      fst = as.double(strsplit(as.character(param[idx,2]),split=",")[[1]])
    }
    if (param[idx,1]=='--case'){
      case = as.double(strsplit(as.character(param[idx,2]),split=",")[[1]])
    }
    if (param[idx,1]=='--marker'){
      marker = as.numeric(as.character(param[idx,2]))
    }
    if (param[idx,1]=='--replicate'){
      replicate = as.numeric(as.character(param[idx,2]))
    }
    if (param[idx,1]=='--riskratio'){
      riskratio = as.numeric(as.character(param[idx,2]))
    }
    if (param[idx,1]=='--no.case.snp'){
      no.case.snp = as.numeric(as.character(param[idx,2]))
    }
    if (param[idx,1]=='--missing'){
      missing = as.numeric(as.character(param[idx,2]))
    }
    if (param[idx,1]=='--outlier'){
      outlier = as.numeric(strsplit(as.character(param[idx,2]),split=",")[[1]])
    }
    if (param[idx,1]=='--pc'){
      pc = as.character(param[idx,2])
    }
    if (param[idx,1]=='--fulloutput'){
      fulloutput = as.character(param[idx,2])
    }
    if (param[idx,1]=='--cate'){
      cate = as.character(param[idx,2])
    }
  }
  if ((length(population)+length(fst)+length(case)+length(outlier))/4.0 != length(population)){
    cat(paste0("Parameters are not correct! The given lists of population, fst, case, and outlier need to be the same length!\n"))
    quit()
  }
  ret = list("name"=name,"population"=population,"fst"=fst,"case"=case,"marker"=marker,"replicate"=replicate,
             "riskratio"=riskratio,"no.case.snp"=no.case.snp,"missing"=missing,"outlier"=outlier,"pc"=pc,
             "fulloutput"=fulloutput,"cate"=cate)
  return(ret)
}

#' Simulate data for multiple populations
#'
#' The output files are saved to the specified directory according to \code{out}.
#'
#' @param setting An absolute path to a setting file
#' @param out An absolute path for output files
#' @param thread A number to specify a maximum thread to be run in parallel
#'
#' @return NULL
#'
#' @details This function takes the specific input file containing the settings
#' for simulations. It allows multiple settings for several simulation within
#' one file. The simulation-setting file must be a text file. The line started
#' with "--" indicates the parameters for simulation, and the line started with
#' "#" are comments. Empty lines are allowed in the setting file. The parameters
#' in the setting file are listed below:
#' \itemize{
#' \item \code{--setting}	A name of setting
#' \item \code{--population} A list that indicates the numbers of population
#' separated by comma
#' \item \code{--fst}	A list that indicates the Fst values separated by comma.
#' Each Fst value represents a genetic distance of that particular population
#' and the first population. The Fst values for the first population and the
#' second population should be the same values, otherwise they will be summed up
#' and devided by two.
#' \item \code{--case} A list that indicates the ratio values of cases separated
#' by comma
#' \item \code{--outlier}	A list that indicates the logical values (0/1) whether
#' that population are outliers, separated by comma
#' \item \code{--marker} A number of SNPs
#' \item \code{--replicate}	A number of replicates
#' \item \code{--riskratio} A number of replicates
#' \item \code{--no.case.snp}	A number of case SNPs
#' \item \code{--pc} A logical value (TRUE/FALSE) whether PCs will be calculated.
#' \item \code{--fulloutput} A logical value (TRUE/FALSE) whether all
#' information will be exported.
#' }
#'
#' @export
#'
# @import doParallel
# @import parallel
# @import foreach
#' @import KRIS
#' @importFrom utils read.table write.table write.csv
#' @importFrom grDevices dev.off pdf
#'
#' @examples
#'
#' #See the examples from
#' \donttest{
#' demo.filest()
#' }
#'
filest <- function(setting, out, thread = 1){

  #To run this script
  #$ Rscript [sript].r --setting [setting file] --outdir [output directory] --thread [number of thread]
  start.time = Sys.time()
  cat(paste0("Start [S0] at ",format(start.time),"\n"))
  args = commandArgs()

  no.thread = as.integer(thread)
  out.dir.prefix = out
  fname.input = setting

  if (length(no.thread)==0){
    no.thread = 1
  }

  #registerDoMC(no.thread)
  #parallel.clusters=makeCluster(no.thread,type="SOCK")
  #registerDoParallel(parallel.clusters)

  options(scipen=15)
  options(digits=18)
  set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )

  cat(paste0("Setting file is : ",fname.input,"\n"))

  param.all = read.table(fname.input,header=F,sep="=")
  settings = which(param.all=='--setting')

  idx.rep <- NULL
  #foreach (set.no = 1:length(settings))  %dopar% {
  for (set.no in 1:length(settings)) {
    if (set.no != length(settings)){
      idx.range = settings[set.no]:(settings[(set.no+1)]-1)
    }else{
      idx.range = settings[set.no]:length(param.all[,1])
    }
    param=get.para(param.all[idx.range,])
    out.dir =  file.path(out.dir.prefix,param$name)
    cat(paste0("The simulated files will be saved in this directory: ",out.dir,"\n"))
    dir.create(out.dir, showWarnings = TRUE, mode = "0755")
    #foreach (idx.rep = 1:param$replicate) %dopar% {
    for (idx.rep in 1:param$replicate) {

      out.filename.prefix = paste0(file.path(out.dir,"simSNP_rep"),idx.rep)
      cat(paste0("Creating data file setting #",set.no," - rep #",idx.rep,"\n"))
      no.sample = sum(param$population)
      no.missing.marker = floor(as.double(no.sample) * as.double(param$missing))

      #raw.data = generate.snp(pop=c(no.pop1,no.pop2,no.pop3,no.pop4),fst=c(FST1,FST2,FST3,FST4),no.snp=no.marker)
      #raw.data = generate.categorical.data(pop=c(no.pop1,no.pop2,no.pop3,no.pop4),fst=c(FST1,FST2,FST3,FST4),no.marker=no.marker)
      raw.data = generate.snp(pop=param$population,fst=param$fst,no.snp=param$marker)

      cat(paste0("Done - ",format(Sys.time()-start.time),"\n"))

      if (sum(param$outlier)>0){
        cat(paste0("Create outliers #",set.no," - rep #",idx.rep,"\n"))
        raw.data = create.outlier(raw.data,param$population,param$outlier)
        cat(paste0("Done - ",format(Sys.time()-start.time),"\n"))
      }

      #introduce Missing Genotype
      no.row = dim(raw.data)[1]
      no.col = dim(raw.data)[2]
      if (param$missing > 0){
        cat(paste0("Create missing data #",set.no," - rep #",idx.rep,"\n"))
        no.missing = round(param$missing * sum(param$population))
        if (no.missing>0){
          for (i in 1:no.missing){
            idx.row = sample.int(dim(raw.data)[1],size=1)
            idx.col = sample.int(dim(raw.data)[2],size=1)
            raw.data[idx.row,idx.col] = -1
          }
        }
        cat(paste0("Done - ",format(Sys.time()-start.time),"\n"))
      }

      object = list()

      if (!anyNA(param$no.case.snp)){
        if (param$no.case.snp > 0){
          object.case = list()
          object.ctrl = list()
        }
      }

      cat(paste0("Writing data files setting #",set.no," - rep #",idx.rep,"\n"))
      #Numeric data - individuals as row, variables as column
      if (param$fulloutput == "TRUE"){
        filename = paste(out.filename.prefix,"_data_numMark_rowInd_colVar.txt",sep="")
        write.table(raw.data,file=filename,sep=" ",col.names=F, row.names=F, quote=F)
      }

      #Numeric data - for ipPCA's input, individuals as column, variables as row
      filename = paste(out.filename.prefix,"_data_numMark_rowVar_colInd.txt",sep="")
      write.table(t(raw.data),file=filename,sep=" ",col.names=F, row.names=F, quote=F)
      cat(paste0("Done - ",format(Sys.time()-start.time),"\n"))

      if (param$cate == "TRUE"){
        #Categorical data - factominer input, individuals as row, variables as column
        cat(paste0("Creating Categorical files setting #",set.no," - rep #",idx.rep,"\n"))
        cate.data = c()
        marker_id = paste0(rep('m',param$marker),1:param$marker)
        cate.data = t(apply(raw.data,1,do.create.cate,t(as.matrix(marker_id))))
        #for (i in 1:length(marker_id)){
        #  tmp = paste(marker_id[i],raw.data[,i],sep="_")
        #  cate.data = cbind(cate.data, tmp)
        #}
        filename = paste(out.filename.prefix,"_data_catMark_rowInd_colVar.txt",sep="")
        write.table(cate.data,file=filename,sep=" ",col.names=F, row.names=F, quote=F)

        #Numeric data encoded for categorical - transposed format , individuals as column, variables as row
        #       cat(paste0("Creating Encoded Categorical files #",set.no," - rep #",idx.rep,"\n"))
        #       encoded.cate.data = c()
        #       marker_id = paste0(rep('m',param$marker),1:param$marker)
        #       marker.variable=c()
        #       for (i in 1:length(marker_id)){
        #         for (j in sort(unique(raw.data[,i]))){
        #           variable = paste0(marker_id[i],"_",j)
        #           marker.variable = c(marker.variable,variable)
        #           tmp = rep.int(0,no.sample)
        #           tmp[which(raw.data[,i]==j)] = 1
        #           encoded.cate.data = cbind(encoded.cate.data,tmp)
        #         }
        #       }
        #
        #       filename = paste(out.filename.prefix,"_data_catEncodeMark_rowVar_colInd.txt",sep="")
        #       write.table(t(encoded.cate.data),file=filename,sep=" ",col.names=F, row.names=F, quote=F)
        #
        #       filename = paste(out.filename.prefix,"_data_catEncodeMark_markerHeader.txt",sep="")
        #       write.table(marker.variable,file=filename,sep=" ",col.names=F, row.names=F, quote=F)

        cat(paste0("Done - ",format(Sys.time()-start.time),"\n"))
      }

      cat(paste0("Creating status file setting #",set.no," - rep #",idx.rep,"\n"))

      #create disease status, 1 = uneffected, 2=affected
      status = rep('1',(no.sample))
      status.label = rep('control',(no.sample))
      no.case = param$case * param$population

      #generate the disease maker, but first generate normal SNPs
      if (!anyNA(param$no.case.snp)){
        if (param$no.case.snp > 0){
          case.snp.set = generate.snp(pop=param$population,fst=param$fst,no.snp=param$no.case.snp,riskratio = 1.0)
          case.snp.only = generate.snp(pop=no.case,fst=param$fst,no.snp=param$no.case.snp,riskratio = param$riskratio)
        }
      }

      last.idx = 0
      last.idy = 0
      for (i in 1:length(param$population)){
        if (no.case[i]>0){
          idx = (last.idx + (param$population[i]-no.case[i]) + 1):(param$population[i] + last.idx)
          status[idx] = 2
          status.label[idx] = "case"

          if (!anyNA(param$no.case.snp)){
            if (param$no.case.snp > 0){
              idy = (last.idy +  1):(no.case[i] + last.idy)
              case.snp.set[idx,] = case.snp.only[idy,]
              last.idy = last.idy + no.case[i]
            }
          }

        }
        last.idx = last.idx + param$population[i]
      }

      #create labels

      iid = 1:no.sample
      label = generate.label(pop=param$population,outlier=param$outlier)
      plot.label = paste0(label,"_",status.label)
      pid = rep('0',no.sample)
      sex = rep('1',no.sample)
      individuals = cbind(iid,label,status)
      object$ind.info = cbind(iid,iid,pid,pid,sex,status)
      idx.ctrl = which(status == 1)
      idx.case = which(status == 2)
      object$ind.info = cbind(iid,iid,pid,pid,sex,status)

      if (!anyNA(param$no.case.snp)){
        if (param$no.case.snp > 0){
          object.ctrl$ind.info = object$ind.info[idx.ctrl,]
          object.case$ind.info = object$ind.info[idx.case,]
        }
      }

      colnames(individuals) = c("individual","population","disease")
      filename = paste0(out.filename.prefix,"_individuals.txt")
      write.table(individuals,file=filename,sep=" ",col.names=F, row.names=F, quote=F)
      if (param$fulloutput == "TRUE"){
        filename = paste0(out.filename.prefix,"_individuals_with_header.txt",sep="")
        write.table(individuals,file=filename,sep=" ",col.names=T, row.names=F, quote=F)
      }
      cat(paste0("Done - ",format(Sys.time()-start.time),"\n"))

      #estimate Fst
      cat(paste0("Estimating Fst setting #",set.no," - rep #",idx.rep,"\n"))
      unique.label = unique(label)
      if (length(unique.label)>1){
        fst.tab = matrix(NA,length(unique.label),length(unique.label))
        fst.col.row.name = matrix(NA,length(unique.label),1)
        for (i in 1:(length(unique.label)-1)) {
          for (j in (i+1):length(unique.label)) {
            idx1=which(label==unique.label[i])
            idx2=which(label==unique.label[j])
            fst.col.row.name[i] = paste0(unique.label[i],'(',length(idx1),')')
            fst.col.row.name[j] = paste0(unique.label[j],'(',length(idx2),')')
            dat1=raw.data[idx1,]
            dat2=raw.data[idx2,]
            new.data=rbind(dat1,dat2)
            p1=1:length(idx1)
            p2=(length(idx1)+1):(length(idx1)+length(idx2))
            f_st=fst.hudson(new.data,p1,p2)
            fst.tab[i,j]=f_st
            fst.tab[j,i]=f_st
          }
        }
        colnames(fst.tab) = fst.col.row.name
        rownames(fst.tab) = fst.col.row.name

        filename = paste0(out.filename.prefix,"_estimated_Fst.txt")
        write.csv(fst.tab,file=filename,quote=F,row.names=T,na='')
      }

      cat(paste0("Done - ",format(Sys.time()-start.time),"\n"))
      cat(paste0("Creating maker information setting #",set.no," - rep #",idx.rep,"\n"))

      #create marker information
      marker_id = paste0(rep('marker',param$marker),1:param$marker)
      marker_position = 1:param$marker * 10000
      marker_chr = rep('1',param$marker)
      marker_morgan = rep('0',param$marker)
      marker_allele1 = rep('A',param$marker)
      marker_allele2 = rep('T',param$marker)

      if (!anyNA(param$no.case.snp)){
        if (param$no.case.snp > 0){
          caseSNP_id = paste0(rep('caseSNP',param$no.case.snp),1:param$no.case.snp)
          caseSNP_position = (param$marker * 10000) + (1:param$no.case.snp * 10000)
          caseSNP_chr = rep('1',param$no.case.snp)
          caseSNP_morgan = rep('0',param$no.case.snp)
          caseSNP_allele1 = rep('A',param$no.case.snp)
          caseSNP_allele2 = rep('T',param$no.case.snp)

          marker_id = c(marker_id,caseSNP_id)
          marker_position = c(marker_position,caseSNP_position)
          marker_chr = c(marker_chr,caseSNP_chr)
          marker_morgan = c(marker_morgan,caseSNP_morgan)
          marker_allele1 = c(marker_allele1,caseSNP_allele1)
          marker_allele2 = c(marker_allele2,caseSNP_allele2)
        }
      }

      object$snp.info = cbind(marker_chr,marker_id,marker_morgan,marker_position,marker_allele1,marker_allele2)

      if (!anyNA(param$no.case.snp)){
        if (param$no.case.snp > 0){
          object.ctrl$snp.info = object$snp.info
          object.case$snp.info = object$snp.info
        }
      }

      #Combine Null SNPs with case SNPs

      object$snp = raw.data

      if (!anyNA(param$no.case.snp)){
        if (param$no.case.snp > 0){
          object$snp = cbind(raw.data,case.snp.set)

          object.ctrl$snp = object$snp[idx.ctrl,]
          object.case$snp = object$snp[idx.case,]

          filename = paste0(out.filename.prefix,"_controls")
          write.bed(object.ctrl,file=filename)

          filename = paste0(out.filename.prefix,"_cases")
          write.bed(object.case,file=filename)
        }
      }

      filename = paste0(out.filename.prefix)
      write.bed(object,file=filename)

      cat(paste0("Done - ",format(Sys.time()-start.time),"\n"))

      if (param$pc == "TRUE"){
        cat(paste0("Generating PC scores #",set.no," - rep #",idx.rep,"\n"))
        PC = cal.pc.linear(object$snp,PCscore=TRUE,no.pc=10,data.type="snp")$PC

        filename = paste0(out.filename.prefix,"_PC10.txt")
        write.table(cbind(iid,iid,PC),file=filename,sep=" ",col.names=F, row.names=F, quote=F)

        filename = paste0(out.filename.prefix,"_PC.pdf")
        pdf(filename)
        plot3views(PC,plot.label)
        dev.off()

        cat(paste0("Generating EigenVector  #",set.no," - rep #",idx.rep,"\n"))
        PC = cal.pc.linear(object$snp,PCscore=FALSE,no.pc=10,data.type="snp")$PC

        filename = paste0(out.filename.prefix,"_eigenvector10.txt")
        write.table(cbind(iid,iid,PC),file=filename,sep=" ",col.names=F, row.names=F, quote=F)

        PC.projected = NULL
        if (sum(no.case) > 0){
          cat(paste0("Generating projected PCs  #",set.no," - rep #",idx.rep,"\n"))
          PC.projected = cal.pc.projection(object$snp, iid, status, label, no.pc = 10)

          filename = paste0(out.filename.prefix,"_PC.projected.txt")
          tmp.pc.projected = cbind(PC.projected$id,PC.projected$id,PC.projected$label,PC.projected$status,PC.projected$PC)
          #colnames(tmp.pc.projected) <- c('fid','iid','pop','status','pc1','pc2','pc3','pc4','pc5','pc6','pc7','pc8','pc9','pc10')
          #PC.projected.sorted <- tmp.pc.projected[order(fid),]
          PC.projected.sorted <- tmp.pc.projected

          write.table(PC.projected.sorted,file=filename,sep=" ",col.names=F, row.names=F, quote=F)

          filename = paste0(out.filename.prefix,"_PC.projected.pdf")
          pdf(filename)
          plot.label = paste0(PC.projected$label,"_",PC.projected$status)
          plot3views(PC.projected$PC,plot.label)
          dev.off()
        }


        #Save R object
        filename = paste0(out.filename.prefix,".RData")
        if (!anyNA(param$no.case.snp)){
          if (param$no.case.snp > 0){
            save(file=filename,PC,PC.projected,raw.data,case.snp.set,label,object)
          }else{
            save(file=filename,PC,PC.projected,raw.data,label,object)
          }
        }else{
          save(file=filename,PC,PC.projected,raw.data,label,object)
        }
      }else{
        #Save R object
        filename = paste0(out.filename.prefix,".RData")
        if (!anyNA(param$no.case.snp)){
          if (param$no.case.snp > 0){
            save(file=filename,raw.data,case.snp.set,label,object)
          }else{
            save(file=filename,raw.data,label,object)
          }
        }else{
          save(file=filename,raw.data,label,object)
        }
      }
      cat(paste0("Done - ",format(Sys.time()-start.time),"\n"))
    }
  }

  invisible(NULL)
}

#' Demonstration the filest function
#'
#' This function generates the setting file and demonstrate how to use \code{\link{filest}}.
#'
#' @return The output directory
#' @export
#'
#' @examples
#'
#' #To run this function, simply call demo.filest()
#' \donttest{
#' demo.filest()
#' }
#'
demo.filest <- function(){
  txt = "--setting=example1\n"
  txt = paste0(txt,"--population=500,500\n")
  txt = paste0(txt,"--fst=0.005,0.005\n")
  txt = paste0(txt,"--case=0,0\n")
  txt = paste0(txt,"--outlier=0,0\n")
  txt = paste0(txt,"--marker=3000\n")
  txt = paste0(txt,"--replicate=1\n")
  txt = paste0(txt,"--riskratio=1\n")
  txt = paste0(txt,"--no.case.snp=0\n")
  txt = paste0(txt,"--pc=TRUE\n")
  txt = paste0(txt,"--missing=0\n")
  txt = paste0(txt,"--fulloutput=TRUE\n")

  outdir = file.path(tempdir())

  settingfile = file.path(outdir,"example1.txt")
  cat(paste0("Creating a setting file ... ",settingfile,"\n"))
  fo = file(settingfile,"w")
  for (i in txt){ write(i,fo)}
  close(fo)

  cat(paste0("Generating the simulated data  to  ... ",outdir,"\n"))
  filest(setting = settingfile, out = outdir, thread = 1)

 return(outdir)
}
