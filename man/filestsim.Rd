\name{filestsim}
\alias{filestsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
A function to simulate SNP data which genetic distances between populations are controlable.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function operates on a Nx3 matrix, where N is the number of samples and data are collected on 3 variables.
}
\usage{
filestsim(setting, out, thread = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{setting}{
%%     ~~Describe \code{x} here~~
A setting file. See \emph{demo.filesttim} for demonstration of setting files.
}
  \item{out}{
%%     ~~Describe \code{x} here~~
To set an absolute path for FilestSim output. If the specified output directory already exists, the old simulated files will be replaced.
}
  \item{thread}{
%%     ~~Describe \code{x} here~~
To specify a number of threads or parallel processes. Default = 1.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function \emph{rubikClust} is able to take up to 3 variables (N x 3 matrix). In case, a matrix contains more than 3 columns, only the first three columns are used;  the other columns are ignored.
}
\value{

  The returned value is a vector of numbers representing cluster memberships.


}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data(simSNP) #to load simulated data set
groups <- rubikClust(PC[,1:3])
print(groups)
table(groups)

}

