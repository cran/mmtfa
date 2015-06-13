\name{mmtfa}
\alias{mmtfa}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
mmtfa: Function for Model-Based Clustering and Classification with Mixtures of Multivariate t Factor Analyzers
}
\description{
Fits mixtures of multivariate modified t-factor analyzers via the alternating expectation-conditional maximization algorithm to the given data under a clustering (default) or classification paradigm (by giving either training index or percentage of data taken to be known) in serial or parallel.
}
\usage{
mmtfa(x, Gs=1:4, Qs=1:2, clas=0, init="kmeans", scale=TRUE, models = "all", 
        dfstart=50, dfupdate="approx", gauss=FALSE, eps=0.05, known=NULL, 
        parallel.cores=FALSE)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
A numeric matrix or data frame.
}
  \item{Gs}{
%%     ~~Describe \code{Gs} here~~
An integer or integer vector indicating the number of groups to fit. Default is 1-4.
}
\item{Qs}{
%%     ~~Describe \code{Qs} here~~
An integer or integer vector indicating the number of factors to fit. Default is 1-2.
}
\item{clas}{
Integer between 0--100 giving the percentage of data taken to be known. Simulates a classification scenario. Additional options to be provided in future updates...
}

  \item{init}{
%%     ~~Describe \code{init} here~~
A list of initializing classification of the form that \code{init[[G]]} contains the initializing vector for all G considered (see example below). Alternatively, the user can use a character string indicating initialization method. Currently the user can choose from \code{"kmeans"} (default), `hard" random - \code{"hard"}, `soft" random - \code{"soft"}, and \code{"uniform"} (classification only).
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
Logical indicating whether or not the function should scale the data. Default is \code{TRUE} and is the prescribed method --- tEIGEN models are not scale invariant.
}
  \item{models}{
%%     ~~Describe \code{models} here~~
A character string or character vector giving the models to fit. See details for instructions on choices.
}
  \item{dfstart}{
%%     ~~Describe \code{dfstart} here~~
The initialized value for the degrees of freedom. The default is 50.
}
\item{dfupdate}{
%%     ~~Describe \code{dfstart} here~~
Character string (\code{"approx"} uses a closed form aqpproximation, \code{"numeric"} uses the \code{uniroot} function) or logical \code{FALSE} to skip degrees of freedom updates. 
}
  \item{known}{
%%     ~~Describe \code{known} here~~
A vector of known classifications that can be numeric or character - optional for clustering, necessary for classification. Must be the same length as the number of rows in the data set.
}
  \item{gauss}{
Logical indicating if the algorithm should use the gaussian distribution. Currently equivalent to setting \code{dfstart=200} and \code{dfupdate=FALSE}. Will be improved in later updates...
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
Tolerance value for the convergence criterion for the AECM algorithm.
}
  \item{parallel.cores}{
%%     ~~Describe \code{eps} here~~
Logical or integer specifying number of computing cores to utilize for coarse-grain parallelization of the algorithm. If \code{FALSE} (default), then the algorithm is run in serial. If \code{TRUE}, then \code{detectCores()} from the \code{parallel} package will determine the number of cores and use all of them.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Model specification (via the \code{models} argument) follows nomenclature developed from the factor analyzer decomposition of the covariance matrix. The nomenclature refers to the decomposition and constraints on the covariance matrix: \deqn{\Sigma_g = \Lambda_g \Lambda_g' + \omega_g \delta_g}

The first letter can be a \code{"C"} (constrained across groups) or \code{"U"} (unconstrained) and that refers to setting \eqn{\Lambda_g = \Lambda} or not, respectively. The second letter has the same choices, resulting in \eqn{\omega_g = \omega} or not. The third letter is permitted as a \code{"C"}, \code{"U"}, or \code{"I"} (constrained to be the identity matrix), applying those constraints to \eqn{\delta_g}. The fourth, and final, letter refers to the degrees of freedom, and again is permitted \code{"C"} or \code{"U"}.

As many models as desired can be selected and ran via the vector supplied to \code{models}. The complete list of possible names is: "UUUU", "UUUC", "UCCU","UCCC", "UUIU", "UUIC", "UCIU", "UCIC", "CUUU", "CUUC", "CCCU", "CCCC", "CUIU", "CUIC", "CCIU", "CCIC", "CUCU", "CUCC", "UUCU", "UUCC", "UCUU", "UCUC", "CCUU", "CCUC".

More commonly, subsets can be called by the following character strings: 
	\code{"all"} runs all 24 MMtFA models (default),
	\code{"dfunconstrained"} runs the 12 unconstrained degrees of freedom models, 
	\code{"dfconstrained"} runs the 12 constrained degrees of freedom models,  
	

Also note that for G=1, several models are equivalent (for example, UUUU and CCCC). Thus, for G=1 only one model from each set of equivalent models will be run.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of "comp1"}
%%  \item{comp2 }{Description of "comp2"}
%% ...
\item{x}{Data used for clustering/classification.}

\item{classification}{Vector of group classifications as determined by the BIC.}

\item{bic}{BIC of the best fitted model.}

\item{modelname}{Name of the best model according to the BIC.}

\item{allbic}{Matrix of BIC values according to model and G. A value of -Inf is returned when the model did not converge. }

\item{bestmodel}{Character string giving best model (BIC) details.}

\item{G}{Value corresponding to the number of components chosen by the BIC. }

\item{tab}{Classification table for BIC-selected model (only available when \code{known} is given). When classification is used the 'known' observations are left out of the table. }

\item{fuzzy}{The fuzzy clustering matrix for the model selected by the BIC. }

\item{logl}{The log-likelihood corresponding to the model with the best BIC. }

\item{iter}{The number of iterations until convergence for the model selected by the BIC.}

\item{parameters}{List containing the fitted parameters: \code{mean} - matrix of means where the rows correspond to the component and the columns are the variables; \code{sigma} - array of covariance matrices (multivariate). Upcoming updates will increase the number of saved parameters...}

\item{iclresults}{List containing all the previous outputs, except \code{x} and \code{index}, pertaining to the model chosen by the best ICL (all under the same name except \code{allicl} and \code{icl} are the equivalent of \code{allbic} and \code{bic}, respectively).}

}
\references{
Andrews JL and McNicholas PD (2011a), 'Extending mixtures of multivariate t-factor analyzers'.
\emph{Statistics and Computing} 21(3), 361--373.

Andrews JL and McNicholas PD (2011b), 'Mixtures of modified t-factor analyzers for model-based
clustering, classification, and discriminant analysis'. \emph{Journal of Statistical Planning and Inference} 141(4), 1479--1486.
}
\author{
%%  ~~who you are~~
Jeffrey L. Andrews, Paul D. McNicholas, and Mathieu Chalifour
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See package manual \code{\link{MMtFA}}
}
\examples{
###Note that only one model is run for each example
###in order to reduce computation time

#Clustering old faithful data with hard random start
tfaith <- mmtfa(faithful, models="UUUU", Gs=1:3, Qs=1, init="hard")

#Clustering old faithful with hierarchical starting values
initial_list <- list()
clustree <- hclust(dist(faithful))
for(i in 1:3){
	initial_list[[i]] <- cutree(clustree,i)
}
tfaith <- mmtfa(faithful, models="CUCU", Gs=1:3, Qs=1, init=initial_list)

#Classification with the iris data set via percentage of data taken to have known membership
tiris <- mmtfa(iris[,-5], Qs=1:2, models="CUIU", init="uniform",clas=50, known=iris[,5])
tiris$tab

}


