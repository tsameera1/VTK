

#include "math.h"
#include <algorithm>
#include "vtkAxisExtended.h"
#include "vtkStdString.h"
#include "vtksys/ios/sstream"
#include "vtkObjectFactory.h"



// This implements the optimization based tick position calculating algorithm in the paper "An Extension of Wilkinson's Algorithm
// for Positioning Tick Labels on Axes" by Junstin Talbot, Sharon Lin and Pat Hanrahan

vtkStandardNewMacro(vtkAxisExtended);
// This method return a value to make step sizes corresponding to low q and j values more preferable
double vtkAxisExtended::Simplicity(int qIndex, int QLength, int j, double lmin, double lmax, double lstep)
{
  double eps = VTK_DBL_EPSILON * 100;
  int v = 1 ;
  qIndex++;


  double rem = fmod(lmin,lstep);
  if((rem < eps || (lstep - rem ) < eps ) &&  lmin <=0 && lmax >=0 )
    {
    v =0;
    }
  else
    {
    v=1;  // v is 1 is lebelling includes zero
    }

  return 1 - (qIndex-1)/(QLength-1) - j + v; 

}

// This method returns the maximum possible value of simplicity value given q and j
double vtkAxisExtended::SimplicityMax(int qIndex, int QLength, int j)
{
  int v =1;
  qIndex++;
  return 1 - (qIndex-1)/(QLength-1) - j + v; 

}

// This method makes the data range approximately same as the labeling range more preferable
double vtkAxisExtended::Coverage(double dmin, double dmax, double lmin, double lmax)
{
  double coverage = 1- 0.5 * (pow(dmax- lmax, 2) + pow(dmin- lmin, 2) / pow(0.1*(dmax-dmin),2));
  return coverage;

}


//This gives the maximum possible value of coverage given the step size
double vtkAxisExtended::CoverageMax(double dmin, double dmax, double span)
{
  double range = dmax - dmin;
  if (span > range)
    {
    double half = (span - range)/2;
    return 1- 0.5 * (pow(half, 2) + pow(half, 2) / pow(0.1*(range),2));
    }
  else
    {
    return 1.0;
    }

}

// This method return a value to make the density of the labels close to the user given value
double vtkAxisExtended::Density(int k, double m, double dmin, double dmax, double lmin, double lmax)
{
  double r = (k-1)/(lmax-lmin);
  double rt = (m-1) / (std::max(lmax,dmax) - std::min(dmin,lmin));

  return 2 - std::max(r/rt , rt/r);

}

// Derives the maximum values for density given k (number of ticks) and m (user given)
double vtkAxisExtended::DensityMax(int k, double m)
{
  if(k >= m)
    {
    return 2 - (k-1)/(m-1);
    }
  else
    {
    return 1;
    }
}


// This methods gives a weighing factor for each label depending on the range
// The coding for the different formats are
//   1 - Decimal e.g. 5000
//   2 - Factored Decimals e.g. 5 (thousands)
//   3 - K e.g. 5K
//   4 - Factored K e.g. 5(K)
//   5 - M e.g. 5M
//   6 - Factored M e.g. 5(M)
//   7 - Scientific 5 * 10^6
//   8 - Factored Scientific 5 (10^6)


double vtkAxisExtended::FormatLegibilityScore(double n, int format)
{
  switch(format)
  {
    case 1:
      if( abs(n)>0.0001 &&  abs(n) < 1000000)
        {
        return 1.0;
        }
      else
        {
        return 0.0;
        }
      break;
    case 2:
      return 0.5;
      break; 
    case 3:
      if( abs(n)>1000 &&  abs(n) < 1000000)
        {
        return 0.75;
        }
      else
        {
        return 0.0;
        }
      break;
    case 4:
      if( abs(n)>1000 &&  abs(n) < 1000000)
        {
        return 0.4;
        }
      else
        {
        return 0.0;
        }
      break; 
    case 5:
      if( abs(n)>1000000 &&  abs(n) < 1000000000)
        {
        return 0.75;
        }
      else
        {
        return 0.0;
        }
      break;
    case 6:
      if( abs(n)>1000000 &&  abs(n) < 1000000000)
        {
        return 0.4;
        }
      else
        {
        return 0.0;
        }
      break;
    case 7:
      return 0.25;
      break;
    case 8:
      return 0.3;
      break;
    default:
      return 0.0;
      break;
  }
}


// This method returns the length of the label given the format
int vtkAxisExtended::FormatStringLength(int format, double n, int precision)
{
   vtksys_ios::ostringstream ostr;
   ostr.imbue(vtkstd::locale::classic());
      
   ostr.precision(precision);
   ostr.setf(ios::fixed, ios::floatfield);
   ostr << n;

   int numSize = (int) ostr.str().length()-1;  // Gets the length of the string with the current format without the end character
        

  switch(format)
  {
    case 1:
      return numSize;
      break;
    case 2:
      return numSize-3;  // Three 0's get reduced
      break;
    case 3:
      return numSize-2; // minus three zeros + K
      break;
    case 4:
      return numSize-3; // minus three zeros
      break;
    case 5:
      return numSize-5; // minus six zeros 
      break;
    case 6:
      return numSize-6; // minus six zeros + M
      break;
    case 7:
      ostr.setf(vtksys_ios::ios::scientific, vtksys_ios::ios::floatfield);
      ostr<<n;
      numSize = (int) ostr.str().length();
      return numSize;
      break;
    case 8:
      ostr.setf(vtksys_ios::ios::scientific, vtksys_ios::ios::floatfield);
      ostr<<n;
      numSize = (int) ostr.str().length();
      return numSize-3;
      break;
    default:
      return 0;
      break;
  }

}


// This methods determinies the optimum notation, font size and orientation of labels from an exhaustive search
double vtkAxisExtended::Legibility(double lmin, double lmax, double lstep, double scaling)
{
  int numTicks = floor((lmax-lmin)/lstep);
  double* tickPositions = new double[numTicks];
  for(int i = 0; i< numTicks; i++)
    {
    tickPositions[i] = lmax + i*lstep;
    }
  
  this->LabelLegibilityChanged = true;

  int bestFormat = 1;
  int bestOrientation = 0;
  int bestFontSize = this->FontSize;

  int bestLegScore = 0.0;

  for(int iFormat = 1; iFormat<9 ; iFormat++ )
  {
    double formatLegSum = 0.0;
    for(int i = 0; i<numTicks; i++)
      {
      formatLegSum += FormatLegibilityScore(tickPositions[i], iFormat);
      }

    formatLegSum = formatLegSum / numTicks;  // Average of label legibility scores

    double eps = VTK_DBL_EPSILON * 100;
    int v = 1 ;
    double rem = fmod(lmin,lstep);
    if((rem < eps || (lstep - rem ) < eps ) &&  lmin <=0 && lmax >=0 )
      {
      v = 0;
      }
    else
      {
      v = 1;  // v is 1 is lebelling includes zero
      }

    formatLegSum = 0.9 * formatLegSum + 0.1 * v; 
    
    double fontLegSum = 0.0;
    int desiredFontSize = 12;   /// Desired font size is set to 12
    
    for (int iFont = this->FontSize; iFont >= this->MinFontSize ; iFont--)
      {
      if(iFont == desiredFontSize)
        {
        fontLegSum = 1.0;
        }
      else if ( iFont<desiredFontSize && iFont >= this->MinFontSize)
        {
        fontLegSum = 0.2 * (iFont - this->MinFontSize + 1) / (desiredFontSize - this->MinFontSize);
        }
      else
        {
        fontLegSum = -100.0 ;
        }
        
      for(int iOrientation = 0 ; iOrientation <2 ; iOrientation++)
        {
        double orientLegSum = (iOrientation == 0) ? 1 : -0.5; 

        // Here the gap between two consecutive labels is calculated as:
        // 2*Actual distance (in pixels) among two ticks - string lengths of the two largest labels
        double tickDistance = lstep * scaling;
        double fontExtent = (FormatStringLength(iFormat,tickPositions[numTicks-1],this->Precision) + FormatStringLength(iFormat,tickPositions[numTicks-2],this->Precision))*this->FontSize;
        double labelingGap= 2*(tickDistance) - fontExtent;
        double overlapLegSum = 0.0;

        if(labelingGap > 3*this->FontSize)
          {
          overlapLegSum = 1.0;
          }
        else if(labelingGap < 3*this->FontSize && labelingGap > 0)
          {
          overlapLegSum = 2 - 2* fontExtent / tickDistance ;
          }
        else 
          {
          overlapLegSum = -100;
          }

        double legScore = (formatLegSum + fontLegSum + orientLegSum + overlapLegSum)/4;

        if ( legScore > bestLegScore)
          {
          bestFormat = iFormat;
          bestOrientation = iOrientation;
          bestFontSize = iFont;
          }
        }
      }

      
  }
  this->FontSize = bestFontSize;
  this->Orientation = bestOrientation;
  this->LabelFormat = bestFormat;
  delete [] tickPositions;
  return bestLegScore;

}

// This method implements the algorithm given in the paper
double* vtkAxisExtended::GenerateExtendedTickLabels(double dmin, double dmax, double m, double scaling)
{
  double Q[] = {1, 5, 2, 2.5, 4, 3};
  double w[] = {0.25, 0.2, 0.5, 0.05};
  double eps = VTK_DBL_EPSILON * 100;
  double ans[3];

  this->LabelLegibilityChanged = false;
  //vtkVector3d ans
  if(dmin > dmax)
    {
    double temp = dmin;
    dmin = dmax;
    dmax = temp;
    }

  if( dmax - dmin < eps)
    {
    ans[0] = dmin; ans[1]= dmax; ans[2]= m;
    //return Sequence(dmin,dmax, m);
    return ans;
    }

  int QLength = 6;//Q.Length(); // Hard Coded

  //list<double> best; 
  double bestScore = -2;
  double bestLmin(0),bestLmax(0), bestLstep(0);

  const int INF = 100; //INT_MAX;  Again 1000 is hard coded

  int j =1;
  while(j < INF)
    {
    for(int qIndex = 0; qIndex < QLength; qIndex++)
      {
      double sm = SimplicityMax(qIndex, QLength, j);
      if((w[0]*sm + w[1] + w[2] + w[3]) < bestScore)
        {
        j = INF;
        break;
        }

      int k = 2;
      while(k < INF)
        {
        double dm = DensityMax(k,m);
        if((w[0]*sm + w[1] + w[2]*dm + w[3]) < bestScore)
          {
          break;
          }
        double delta = (dmax- dmin)/((k+1)*j*Q[qIndex]) ;
        double z = ceil(log10(delta));
        while( z < INF)
          {
          double step = j*Q[qIndex]*pow(10,z);
          double cm = CoverageMax(dmin, dmax, step*(k-1));
          if((w[0]*sm + w[1] + w[2]*dm + w[3]) < bestScore)
            {
            break;
            }

          double minStart = floor(dmax/step)*j - (k-1)*j;
          double maxStart = ceil(dmin/step)*j;

          if(minStart > maxStart)
            {
            z++;
            continue;
            }

          for(int start = minStart; start <=maxStart; start++)
            {
            double lmin = start * (step/j);
            double lmax = lmin + step*(k-1);
            double lstep = step;

            double s = Simplicity(qIndex, QLength, j, lmin, lmax, lstep);
            double c = Coverage(dmin, dmax, lmin, lmax);
            double g = Density(k,m,dmin, dmax, lmin, lmax);
            

            double score = w[0]*s + w[1]*c + w[2]*g + w[3];

            if(score < bestScore)
               continue;

            double l = Legibility(lmin,lmax,lstep,scaling);  

            score = w[0]*s + w[1]*c + w[2]*g + w[3]*l;

            
            if(score > bestScore)
              {
              bestScore = score;
              bestLmin = lmin;
              bestLmax = lmax;
              bestLstep = lstep;
              }
            }
          z++;
          }
        k++;


        }
      }
    j++;

    }
  ans[0] = bestLmin;
  ans[1] = bestLmax;
  ans[2] = bestLstep;
  //vtkVector3d answers(bestLmin, bestLmax, bestLstep);
  // return Sequence(bestLmin, bestLmax, bestLstep);
  return ans;
}

void vtkAxisExtended::SetFontSize(int fontSize)
{
  this->FontSize = fontSize;

}

void vtkAxisExtended::SetMinFontSize(int minFontSize)
{
  this->MinFontSize = minFontSize;
}

void vtkAxisExtended::SetPrecision(int precision)
{
  this->Precision = precision;
}

void vtkAxisExtended::SetFormat(int format)
{
  this->LabelFormat = format;
}

void vtkAxisExtended::SetOrientation(int orientation)
{
  this->Orientation = orientation;
}
    
int vtkAxisExtended::GetFontSize()
{
  return this->FontSize;
}

int vtkAxisExtended::GetLabelFormat()
{
  return this->LabelFormat;
}

int vtkAxisExtended::GetOrientation()
{
  return this->Orientation;
}

bool vtkAxisExtended::GetLabelLegibilityChanged()
{
  return this->LabelLegibilityChanged;
}