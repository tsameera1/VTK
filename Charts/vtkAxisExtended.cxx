

#include "math.h"
#include <algorithm>
#include "vtkAxisExtended.h"


// This implements the optimization based tick position calculating algorithm in the paper "An Extension of Wilkinson's Algorithm
// for Positioning Tick Labels on Axes" by Junstin Talbot, Sharon Lin and Pat Hanrahan


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
        return 2 - (k-1)/(m-1);
    else
        return 1;
}

//double* Sequence(double lmin, double lmax, double lstep)
//{
//	int numticks = floor((lmax-lmin)/lstep);
//	double* tickPositions = new double[numticks];
//	for(int i =0;i< numticks; i++)
//		tickPositions[i] = lmax + i*lstep;
//
//	return tickPositions;
//
//}


// This method implements the algorithm given in the paper
double* vtkAxisExtended::GenerateExtendedTickLabels(double dmin, double dmax, double m)
{
    double Q[] = {1, 5, 2, 2.5, 4, 3};
    double w[] = {0.25, 0.2, 0.5, 0.05};
    double eps = VTK_DBL_EPSILON * 100;
    double ans[3];
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

    const int INF = 1000; //INT_MAX;  Again 1000 is hard coded

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
                        double l = 1.0;   /// LEGIBILITY function should be implemented and called here

                        double score = w[0]*s + w[1]*c + w[2]*g + w[3]*l;
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
    // return Sequence(bestLmin, bestLmax, bestLstep);
    return ans;
}