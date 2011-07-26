


#ifndef __vtkCellTreeLocator_h
#define __vtkCellTreeLocator_h

#include "vtkAbstractCellLocator.h"
#include <vector>

class PointTraversal;
class vtkIdTypeArray;
class vtkCellArray;

class VTK_FILTERING_EXPORT vtkCellTreeLocator : public vtkAbstractCellLocator
{
  public:
    class CellTree;
    class CellTreeNode;

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

    void FreeSearchStructure();
    void GenerateRepresentation(int level, vtkPolyData *pd);
  
    //
    // Internal classes made public to allow subclasses to create
    // customized some traversal algorithms
    //
    class VTK_FILTERING_EXPORT CellTree
    {
      public:
        std::vector<CellTreeNode>  Nodes;
        std::vector<unsigned int> Leaves;
        friend class PointTraversal;
        friend class CellTreeNode;
        friend class CellTreeBuilder;
        friend class vtkCellTreeLocator;

      public:
        float DataBBox[6]; // This store the bounding values of the dataset   
    };

    class VTK_FILTERING_EXPORT CellTreeNode
    {
      public:
       
      protected:
        unsigned int Index;
        float LeftMax;  // left max value
        float RightMin;  // right min value

        unsigned int Sz; // size
        unsigned int St; // start 

        friend class CellTree;
        friend class PointTraversal;
        friend class CellTreeBuilder;
        
      public:
        void MakeNode( unsigned int left, unsigned int d, float b[2] );
        void SetChildren( unsigned int left );
        bool IsNode() const;
        unsigned int GetLeftChildIndex() const;
        unsigned int GetRightChildIndex() const;
        unsigned int GetDimension() const;
        const float& GetLeftMaxValue() const;
        const float& GetRightMinValue() const;
        void MakeLeaf( unsigned int start, unsigned int size );
        bool IsLeaf() const;
        unsigned int Start() const;
        unsigned int Size() const;
    };

protected:
     vtkCellTreeLocator();
    ~vtkCellTreeLocator();

    void BuildLocatorInternal();
    void BuildLocatorIfNeeded();
    void ForceBuildLocator();
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
