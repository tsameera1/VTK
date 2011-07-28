/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCellLocator.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkCellTreeLocator - This class implements the data structures, construction
// algorithms for fast cell location presented in "Fast, Memory-Efficient Cell 
// location in Unstructured Grids for Visualization" by Christop Garth and Kenneth
// I. Joy in VisWeek, 2011.

// .SECTION Description
// Cell Tree is a bounding interval hierarchy based data structure, where child boxes 
// do not form an exact split of the parent boxes along a dimension.  Therefore two axis-
// aligned bounding planes (left max and right min) are stored for each node along a 
// dimension. This class implements the data structure (Cell Tree Node) and its build 
// and traversal algorithms described in the paper.  

// .SECTION Caveats
// 

// .SECTION See Also
// vtkLocator vtkCellLocator vtkModifiedBSPTree


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
    void PrintSelf(ostream& os, vtkIndent indent);

    // Description:
    // Constructor sets the maximum number of cells in a leaf to 8
    // and number of buckets to 5.  Buckets are used in building the cell tree as described in the paper
    static vtkCellTreeLocator *New();

     // Description:
    // Test a point to find if it is inside a cell. Returns the cellId if inside
    // or -1 if not.
    vtkIdType FindCell(double pos[3], double vtkNotUsed, vtkGenericCell *cell,  double pcoords[3], 
                                       double* weights );


        
    // Description:
    // Return intersection point (if any) AND the cell which was intersected by
    // the finite line. The cell is returned as a cell id and as a generic cell.
    // This function is a modification from the vtkModifiedBSPTree class using the
    // data structures in the paper to find intersections.
    int IntersectWithLine(double a0[3], double a1[3], double tol,
                                      double& t, double x[3], double pcoords[3],
                                      int &subId, vtkIdType &cellId,
                                      vtkGenericCell *cell);

    // Description:
    // Return a list of unique cell ids inside of a given bounding box. The
    // user must provide the vtkIdList to populate. This method returns data
    // only after the locator has been built.
    virtual void FindCellsWithinBounds(double *bbox, vtkIdList *cells);


    // Description:
    // Satisfy vtkLocator abstract interface.
    virtual void FreeSearchStructure();
    virtual void GenerateRepresentation(int level, vtkPolyData *pd);
    virtual void BuildLocatorInternal();
    virtual void BuildLocatorIfNeeded();
    virtual void ForceBuildLocator();
    virtual void BuildLocator();
    
  
//BTX
    // Description:
    // Internal classes made public to allow subclasses to create
    // customized some traversal algorithms  
    class VTK_FILTERING_EXPORT CellTree
    {
      public:
        std::vector<CellTreeNode>  Nodes;
        std::vector<unsigned int> Leaves;
        friend class PointTraversal;
        friend class CellTreeNode;
        friend class CellTreeBuilder;
        //friend class vtkCellTreeLocator;

      public:
        float DataBBox[6]; // This store the bounding values of the dataset   
    };

    // Description:
    // This class is the basic building block of the cell tree.  There is a node per dimension. i.e. There are 3 CellTreeNodes
    // in x,y,z directions.  In contrast, vtkModifiedBSPTree class stores the bounding planes for all 3 dimensions in a single node.
    // LeftMax and RightMin defines the bounding planes.
    // start is the location in the cell tree. e.g. for root node start is zero.
    // size is the number of the nodes under the tree 
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
//ETX

protected:
     vtkCellTreeLocator();
    ~vtkCellTreeLocator();

   // Test ray against node BBox : clip t values to extremes
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

  // Order nodes as near/far relative to ray
  void Classify(const double origin[3],
    const double dir[3],
    double &rDist,
    CellTreeNode *&Near, CellTreeNode *&Mid,
    CellTreeNode *&Far, int &mustCheck);

  // From vtkModifiedBSPTRee
  // We provide a function which does the cell/ray test so that
  // it can be overriden by subclasses to perform special treatment
  // (Example : Particles stored in tree, have no dimension, so we must
  // override the cell test to return a value based on some particle size
  virtual int IntersectCellInternal( vtkIdType cell_ID,  const double p1[3],
    const double p2[3],
    const double tol,
    double &t,
    double ipt[3],
    double pcoords[3],
    int &subId);


    int MaxCellsPerLeaf;
    int NumberOfBuckets;

    CellTree* Tree;

    friend class PointTraversal;
    friend class CellTreeNode;
    friend class CellTreeBuilder;
    //friend class vtkCellTreeLocator;

private:
  vtkCellTreeLocator(const vtkCellTreeLocator&);  // Not implemented.
  void operator=(const vtkCellTreeLocator&);      // Not implemented.
};

#endif
