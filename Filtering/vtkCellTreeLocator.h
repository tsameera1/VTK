


#ifndef __vtkCellTreeLocator_h
#define __vtkCellTreeLocator_h

#include "vtkAbstractCellLocator.h"

class CellTree;
class CellTreeNode;
class PointTraversal;
class vtkIdTypeArray;
class vtkCellArray;

class VTK_FILTERING_EXPORT vtkCellTreeLocator : public vtkAbstractCellLocator
{
  public:

    vtkTypeMacro(vtkCellTreeLocator,vtkAbstractCellLocator);

    static vtkCellTreeLocator *New();

   // vtkCellTreeLocator( vtkDataSet* ds );
    //~vtkCellTreeLocator();

    vtkIdType FindCell(double pos[3], double vtkNotUsed, vtkGenericCell *cell,  double pcoords[3], 
                                       double* weights );

    void BuildLocator();

    int IntersectWithLine(double a0[3], double a1[3], double tol,
                                      double& t, double x[3], double pcoords[3],
                                      int &subId, vtkIdType &cellId,
                                      vtkGenericCell *cell);

    bool RayMinMaxT(const double origin[3],
                         const double dir[3],
                         double &rTmin,
                         double &rTmax);

    bool RayMinMaxT(const double bounds[6],
                         const double origin[3],
                         const double dir[3],
                         double &rTmin,
                         double &rTmax);

    int getDominantAxis(const double dir[3]);

    void Classify(const double origin[3],
                       const double dir[3],
                       double &rDist,
                       CellTreeNode *&Near, CellTreeNode *&Mid,
                       CellTreeNode *&Far, int &mustCheck);

  int IntersectCellInternal( vtkIdType cell_ID,  const double p1[3],
    const double p2[3],
    const double tol,
    double &t,
    double ipt[3],
    double pcoords[3],
    int &subId);

    void vtkCellTreeLocator::FreeSearchStructure(void);
    void vtkCellTreeLocator::GenerateRepresentation(int level, vtkPolyData *pd);
protected:
  
    vtkCellTreeLocator();
    ~vtkCellTreeLocator();

    void BuildLocatorInternal();
    void vtkCellTreeLocator::BuildLocatorIfNeeded();
    void vtkCellTreeLocator::ForceBuildLocator();
    void Free();

    int MaxCellsPerLeaf;
    int NumberOfBuckets;

    CellTree* Tree;

    vtkIdType* CellArray;
    vtkIdType* Locations;

private:
  vtkCellTreeLocator(const vtkCellTreeLocator&);  // Not implemented.
  void operator=(const vtkCellTreeLocator&);      // Not implemented.
};

#endif
