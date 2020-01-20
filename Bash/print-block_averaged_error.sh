###############################################################################
#
#  given a series of data (second column of input file), the script computes
#  the standard deviation of the mean based on block averages, as described
#  in Allen & Tildesley "Computer simulation of liquids" (page 192), and in
#  Zhu & Hummer, J.Comput.Chem. 33, 453 (2012) (equation 35).
#
#  in practice, average and variance over the whole sample are estimated,
#  then the sample is divided in nb blocks of length lb and the variance
#  of the block average from the global average is computed (sigma2b).
#  The number of correlated samples is given by the sigma2b*lb/sigma2 value
#  at plateau for large enough block length (under such conditions,
#  the data behaves as an uncorrelated set of Gaussian-distributed values), 
#  and the standard error of the mean is given by sqrt(sigma2b/nb) (last column).
#
#  One should plot column 1 vs column 6 and, if a plateau is detected,
#  use the corresponding error.
#
###############################################################################

awk '{

  n++
  x[n]=$2
  ave += $2
  sigma2 += $2*$2

}END{

  ave /= n
  sigma2 = sigma2/n - ave*ave

  print "# length of block,  num. of blocks, sigma2 global, sigma2 over blocks, sigma2b*lb/sigma2, error=sqrt(sigma2b/nb)"

  for (nb=2; nb<=100; nb++) { 
    lb=int(n/nb)
    i=1
    sigma2b = 0.  
    for (ib=0; ib<nb; ib++){ 
      aveb=0.
      for (j=1; j<=lb; j++){
        aveb += x[i]
        i++
      }
      aveb /= lb
      sigma2b += ( aveb - ave )*( aveb - ave )
    } 
    sigma2b /= nb
    print lb,nb,sigma2,sigma2b,sigma2b*lb/sigma2,sqrt(sigma2b/nb),ave
  }

}' $1
