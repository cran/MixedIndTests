#'@title Computes the Cramer-von Mises statistic Sn for the test of randomness
#'
#'@description This function computes the Cramer-von Mises statistic for a tests of randomness for p consectives values X(1), ..., X(p). Usefull for traditionall bootstrapping.
#'
#'@param x        Vector
#
#'
#'
#'@return \item{values}{Unique (sorted) values}
#'@return \item{m}{Unique (sorted) values}
#'@return \item{Fn}{Empirical cdf of the values}
#'@return \item{fn}{Empirical pdf of the values}
#'
#'@references B.R Nasri (2021). Test of serial dependence for arbitrary distributions
#'@references C. Genest, J.G. Neslehova, B.N. Remillard and O. Murphy (2019). Testing for independence in arbitrary distributions.
#'

#'#'@examples
#'X <- SimAR1Poisson(c(5,0.2),100)
#'out <- Sn_serial(X,5,3)

#'@keywords internal
#'
#'@export
#'
#'

preparedata =function(x){

  values = sort(unique(x))
  n = length(x)
  m = length(values)


  out0 = .C("prepare_data",
            as.double(x),
            as.integer(n),
            as.double(values),
            as.integer(m),
            Fn = double(m),
            fn = double(m)

  )
  out0$values=values
  out0$m = m
  return(out0)
}






