/**
 * copyright  (C) 2004
 * the icecube collaboration
 * @version $Id$
 * @file Direction.h
 * @date $Date$
 */

//***********************************************************
//-- Created: Dusan Turcan, UMD, Sep 2, 2004
//***********************************************************
#include <units.h>
#include <vector>
#include <iostream>
#include <assert.h>
#include <boost/shared_ptr.hpp>

// $Id$

#ifndef HISTOGRAM_H_INCLUDED
#define HISTOGRAM_H_INCLUDED

template <typename T> 
class Histogram
{
 public:


  /**
   * Additional constructor: 2 arguments mean construct dir. with zen,azi
   */
 
  //template <typename T> 
//  Histogram(unsigned int nbins, T emin, T emax);
    Histogram(unsigned int nbins, T emin, T emax):
    nbins(nbins), emin(emin), emax(emax), entries(0)
  {
    assert(nbins > 0);
    assert(emax > emin);

    T delta = (emax-emin)/nbins;

    for (unsigned ibin = 0; ibin < nbins+1;ibin++) 
    {
        bins.push_back(std::pair<T,double>(emin+delta*ibin,0));
    }
  }




  inline unsigned size() { return bins.size(); }

  void insert(T value, double weight =1.0);
  
  std::vector<double> get_bins();

  std::vector<T> get_edges();

  std::pair<T,double> operator [] (int idx) const;

 private:
    unsigned nbins;
    unsigned entries;
    T emin;
    T emax;
    std::vector<std::pair<T,double> > bins;

};


template <typename T> 
void Histogram<T>::insert(T value, double weight)
{
    unsigned index(0);
    while(bins[index].first < value && index < size() ) index++;

    bins[index].second += weight;
    entries++;
}

template <typename T> 
std::pair<T,double> Histogram<T>::operator [] (int idx) const
{
     return bins[idx];
}

template <typename T> 
std::vector<double> Histogram<T>::get_bins()
{
    std::vector<double> ret;
    //vector<int>::iterator ptr;

    for (unsigned ibin = 0; ibin < bins.size(); ibin++)
    {
        std::cout << bins[ibin].first << ": " ;
        std::cout << bins[ibin].second << ", " <<std::endl;
        ret.push_back(bins[ibin].second);
    }
    return ret;
}

template <typename T> 
std::vector<T> Histogram<T>::get_edges()
{
    std::vector<T> ret;
    //vector<int>::iterator ptr;

    for (unsigned ibin = 0; ibin < bins.size(); ibin++)
        ret.push_back(bins[ibin].first);
    return ret;
}

typedef Histogram<double> DoubleHistogram;
typedef boost::shared_ptr<DoubleHistogram> DoubleHistogramPtr; // Histogram shared pointer

typedef Histogram<float> FloatHistogram;
typedef boost::shared_ptr<FloatHistogram> FloatHistogramPtr; // Histogram shared pointer

#endif //HISTOGRAM_H_INCLUDED

