#include <stdlib.h>
#include <cstdio>
#include <stdio.h>
#include <random>
#include <math.h>
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>


// #define ARMA_64BIT_WORD 1
// //// https://stackoverflow.com/questions/40592054/large-matrices-in-rcpparmadillo-via-the-arma-64bit-word-define

//// when armadillo is loaded, //// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]



////////////////////////////////////////////////////////////////////////////////
/// Constants that may be helpful
////////////////////////////////////////////////////////////////////////////////

// #define MIN_NUM = std::numeric_limits<float>::min(); // problem is that it does not know the type
// const double MIN_NUM = std::numeric_limits<float>::min();
// #define PIraw = 3.14159265358979323846;
// const double PI= 3.14159265358979323846;
// # define PI 3.14159265358979323846  /* pi */
//const double MAX_NUM = std::numeric_limits<float>::max();
// double LLMIN = 1.0e+99;

////////////////////////////////////////////////////////////////////////////////
/// PLINK
////////////////////////////////////////////////////////////////////////////////


#define PACK_DENSITY 4
#define PLINK_NA 3

#define PLINK_PHENO_MISSING -9

// The BED file magic numbers
#define PLINK_OFFSET 3

#define COVAR_ACTION_TRAIN_TEST 0
#define COVAR_ACTION_TRAIN_ONLY 1


/* 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions */
#define MASK0 3	  // 3 << 2 * 0 //
#define MASK1 12  // 3 << 2 * 1 //
#define MASK2 48  // 3 << 2 * 2 //
#define MASK3 192 // 3 << 2 * 3 //

#define BUFSIZE 100

/*
 *                   plink BED           sparsnp
 * minor homozyous:  00 => numeric 0     10 => numeric 2
 * heterozygous:     10 => numeric 2     01 => numeric 1
 * major homozygous: 11 => numeric 3     00 => numeric 0
 * missing:          01 => numeric 1     11 => numeric 3
 *
 *
 * http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml says,
 * The bytes in plink are read backwards HGFEDCBA, not GHEFCDAB, but we read
 * them forwards as a character (a proper byte)
 *
 * By default, plink usage dosage of the *major* allele, since allele A1 is
 * usually the minor allele and the code "1" refers to the second allele A2,
 * so that "11" is A2/A2 or major/major.
 *
 * We always use minor allele dosage, to be consistent with the output from
 * plink --recodeA which used minor allele dosage by default.
 *
 * out: array of genotypes
 * in: array of packed genotypes (bytes)
 * n: number of bytes in input
 *
 */
void decode_plink(unsigned char *out,
                  const unsigned char *in, const unsigned int n)
{
  unsigned int i, k;
  unsigned char tmp, geno;
  unsigned int a1, a2;

  for(i = 0 ; i < n ; ++i){
    tmp = in[i];
    k = PACK_DENSITY * i;

    /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
     * allele 2. The final genotype is the sum of the alleles, except for 01
     * which denotes missing.
     */
    geno = (tmp & MASK0);
    a1 = !(geno & 1);
    a2 = !(geno >> 1);
    out[k] = (geno == 1) ? 3 : a1 + a2;
    k++;

    geno = (tmp & MASK1) >> 2;
    a1 = !(geno & 1);
    a2 = !(geno >> 1);
    out[k] = (geno == 1) ? 3 : a1 + a2;
    k++;

    geno = (tmp & MASK2) >> 4;
    a1 = !(geno & 1);
    a2 = !(geno >> 1);
    out[k] = (geno == 1) ? 3 : a1 + a2;
    k++;

    geno = (tmp & MASK3) >> 6;
    a1 = !(geno & 1);
    a2 = !(geno >> 1);
    out[k] = (geno == 1) ? 3 : a1 + a2;
  }
}


// [[Rcpp::export]]
arma::Mat<double> readbed(std::string bedfile,
                          int N, int p,
                          const arma::uvec & myrows,
                          const arma::uvec & mycols,
                          bool verbose=true){
  ////////////////////////////////////////////////////
  // Opening backing file
  arma::Mat<double> X(N,p);
  ////////////////////////////////////////////////////
  // Initializing
  unsigned long len;
  unsigned int np, nsnps;
  ////////////////////////////////////////////////////
  // Opening binary file
  std::ifstream in(bedfile, std::ios::in | std::ios::binary);
  if(!in){
    std::cerr << "[read_bed] Error reading file " << bedfile << std::endl;
    throw std::runtime_error("io error");
  }
  in.seekg(0, std::ifstream::end);
  // file size in bytes, ignoring first 3 bytes (2byte magic number + 1byte mode)
  len = (unsigned int)in.tellg() - 3;
  // Rcpp::Rcout << "The size in bytes of the file is: " << len<< std::endl;

  ////////////////////////////////////////////////////
  // Iterating
  // size of packed data, in bytes, per SNP
  np = (unsigned int) std::ceil((double)N / (double)PACK_DENSITY);
  // Rcpp::Rcout << "Size in bytes of SNP packs is: " << np << std::endl;
  nsnps = len / np;
  in.seekg(3, std::ifstream::beg);

  unsigned char* tmp = new unsigned char[np];

  // Allocate more than the sample size since data must take up whole bytes
  unsigned char* tmp2 = new unsigned char[np * PACK_DENSITY];
  double val=0;
  for(unsigned int j = 0 ; j < nsnps ; j++){
    arma::vec tmp3(N,  arma::fill::zeros);
    in.read((char*)tmp, sizeof(char) * np);
    decode_plink(tmp2, tmp, np);
    for(unsigned int i = 0 ; i < N ; i++){
      val=(double)tmp2[i];
      if(val != PLINK_NA) tmp3[i] = val; // default zeroes
    }
    X.col(j)=tmp3;
  }
  in.close();

  // Subset matrix
  if(myrows.n_elem == X.n_rows){
    X=X.cols(mycols-1);
  }else if(mycols.n_elem == X.n_rows){
    X=X.rows(myrows-1);
  }else{
    X=X.submat(myrows-1,mycols-1);
  }

  return X;
}

// [[Rcpp::export]]
arma::vec sfsloopC(const arma::Mat<int> & X) // genome matrix
  {
 arma::vec n((X).n_cols);
 arma::colvec ones((X).n_rows);
 ones.fill(1);
 n = X.t() * ones;
 return(n);
}
// [[Rcpp::export]]
arma::vec sfsmatC(const arma::Mat<int> & X) // genome matrix
  {
 arma::vec n((X).n_cols);
 arma::colvec ones((X).n_rows);
 ones.fill(1);
 n = X.t() * ones;
 return(n);
}
// [[Rcpp::export]]
arma::vec sfsC(
    const arma::Mat<int> & X, // genome matrix
    const arma::uvec & myrows, // accessions of interests
    const arma::uvec & mycols) // SNPs of interest
  {
  // subset matrix only if necessary
  auto X2 = &X; // create a new name in case is not necessary to subset
  arma::Mat<int> subX; // create new matrix without filling values
  if(myrows.n_elem != X.n_rows & mycols.n_elem != X.n_cols){
    subX=X.submat(myrows-1,mycols-1);
    X2 = &subX;
  }
  // create vector of observations
  arma::vec n((*X2).n_cols);
  n.fill(0); // need to fill otherwise unstable floats
  // iterate and sum number of alleles
  int m,p;
  m= (*X2).n_rows;
  p= (*X2).n_rows;
  for (int i = 0; i < p; i ++) {
    for(int j=0; j < m; j++){
      n(i) += (*X2)(j,i);
    }
  }
  return(n);
}

//Harmonic number now implemented in R
