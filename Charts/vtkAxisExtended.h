#ifndef __vtkAxisExtended_h
#define __vtkAxisExtended_h
#endif

#include "vtkObject.h"
#include "vtkVector.h"
//
#ifndef DBL_EPSILON
#  define VTK_DBL_EPSILON    2.2204460492503131e-16
#else  // DBL_EPSILON
#  define VTK_DBL_EPSILON    DBL_EPSILON
#endif  // DBL_EPSILON



class VTK_CHARTS_EXPORT vtkAxisExtended : public vtkObject
{ 
  public:
    static vtkAxisExtended *New();

    // This method return a value to make step sizes corresponding to low q and j values more preferable
    static double Simplicity(int qIndex, int QLength, int j, double lmin, double lmax, double lstep);

    // This method returns the maximum possible value of simplicity value given q and j
    static double SimplicityMax(int qIndex, int QLength, int j);

    // This method makes the data range approximately same as the labeling range more preferable
    static double Coverage(double dmin, double dmax, double lmin, double lmax);

    //This gives the maximum possible value of coverage given the step size
    static double CoverageMax(double dmin, double dmax, double span);

    // This method return a value to make the density of the labels close to the user given value
    static double Density(int k, double m, double dmin, double dmax, double lmin, double lmax);

    // Derives the maximum values for density given k (number of ticks) and m (user given)
    static double DensityMax(int k, double m);

    // This method implements the algorithm given in the paper
    // The method return the minimum tick position, maximum tick postion and the tick spacing
    static double* GenerateExtendedTickLabels(double dmin, double dmax, double m);

  


  protected:
    vtkAxisExtended() {};
    ~vtkAxisExtended() {};


  private:
    vtkAxisExtended(const vtkAxisExtended&);  // Not implemented.
    void operator=(const vtkAxisExtended&);  // Not implemented.

};