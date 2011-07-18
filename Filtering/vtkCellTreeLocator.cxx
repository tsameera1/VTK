
#include <vtkCellTreeLocator.h>

#include <vtkDataSet.h>
//#include <DebugStream.h>

#include "vtkPolyData.h"

#include <cstdlib>
#include <cassert>
#include <cmath>
#include <stack>
#include <vector>
#include <limits>
#include <algorithm>
#include "vtkObjectFactory.h"
#include "vtkGenericCell.h"
#include "vtkSmartPointer.h"

vtkStandardNewMacro(vtkCellTreeLocator);

const double Epsilon_=1E-8;

enum { POS_X, NEG_X, POS_Y, NEG_Y, POS_Z, NEG_Z };

// -------------------------------------------------------------------------
// This class is the basic building block of the cell tree.  There is a node per dimension. i.e. There are 3 CellTreeNode
// in x,y,z directions.  In contrast, vtkModifiedBSPTree class stores the bounding planes for all 3 dimensions in a single node.
// LeftMax and rm defines the bounding planes.
// start is the location in the cell tree. e.g. for root node start is zero.
// size is the number of the nodes under the tree  


class CellTreeNode
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
    void MakeNode( unsigned int left, unsigned int d, float b[2] ) // b is an array containing left max and right min values
      {
      Index = (d & 3) | (left << 2);
      LeftMax = b[0];
      RightMin = b[1];	
      }

    void SetChildren( unsigned int left )
      {
      Index = GetDimension() | (left << 2); // In index 2 LSBs (Least Significant Bits) store the dimension. MSBs store the position
      }

    bool IsNode() const
      {
      return (Index & 3) != 3;    // For a leaf 2 LSBs in index is 3
      }

    unsigned int GetLeftChildIndex() const
      {
      return (Index >> 2);
      }

    unsigned int GetRightChildIndex() const
      {
      return (Index >> 2) + 1;  // Right child node is adjacent to the Left child node in the data structure
      }

    unsigned int GetDimension() const
      {
      return Index & 3;
      }

    const float& GetLeftMaxValue() const
      {
      return LeftMax;
      }

    const float& GetRightMinValue() const
      {
      return RightMin;
      }

      // ---

    void MakeLeaf( unsigned int start, unsigned int size )
      {
      Index = 3;
      Sz = size;
      St = start;
      }

    bool IsLeaf() const
      {
      return Index == 3;
      }

    unsigned int Start() const
      {
      return St;
      }

    unsigned int Size() const
      {
      return Sz;
      }
};

class CellTree
{
  protected:
    std::vector<CellTreeNode>  Nodes;
    std::vector<unsigned int> Leaves;
    friend class PointTraversal;
    friend class CellTreeNode;
    friend class CellTreeBuilder;
    friend class vtkCellTreeLocator;

  public:
    float DataBBox[6]; // This store the bounding values of the dataset   
};

class PointTraversal
{
  private:
    const CellTree& m_ct;
    unsigned int    m_stack[32]; 
    unsigned int*   m_sp; // stack pointer
    const float*    m_pos; //3-D coordinates of the points

  protected:
   friend class CellTree;
   friend class CellTreeNode;
   friend class CellTreeBuilder;

  public:
    PointTraversal( const CellTree& ct, const float* pos ) :
        m_ct(ct), m_pos(pos)
    {
           
      m_stack[0] = 0; // first element is set to zero

      m_sp = m_stack + 1; //this points to the second element of the stack
    }

    const CellTreeNode* Next()  // this returns n (the location in the CellTree) if it is a leaf or 0 if the point doesnt contain in the data domain
    {
      while( true )
        {
        if( m_sp == m_stack ) //This means the point is not within the domain
          {
          return 0;
          }

        const CellTreeNode* n = &m_ct.Nodes.front() + *(--m_sp);

        if( n->IsLeaf() )
          {
          return n;
          }

        const float p = m_pos[n->GetDimension()];
        const unsigned int left = n->GetLeftChildIndex();

        bool l = p <= n->GetLeftMaxValue(); // Check if the points is within the left sub tree
        bool r = p >= n->GetRightMinValue(); // Check if the point is within the right sub tree

        if( l && r ) //  This means if there is a overlap region both left and right sub trees should be traversed
          {
          if( n->GetLeftMaxValue()-p < p-n->GetRightMinValue() )
            {
            *(m_sp++) = left;
            *(m_sp++) = left+1;
            }
          else
            {
            *(m_sp++) = left+1;
            *(m_sp++) = left;
            }
          }
        else if( l )
          {
          *(m_sp++) = left;
          }
        else if( r )
          {
          *(m_sp++) = left+1;
          }
        }
    }
};



// -------------------------------------------------------------------------

// This class builds the CellTree according to the algorithm given in the paper.

class CellTreeBuilder
{
  private:
      
    struct Bucket
    {
      float        Min;
      float        Max;
      unsigned int Cnt;
    
      Bucket()
      {
        Cnt = 0;
        Min =  std::numeric_limits<float>::max();
        Max = -std::numeric_limits<float>::max();
      }
    
      void Add( const float _min, const float _max )
      {
        ++Cnt;
        
        if( _min < Min )
          {
          Min = _min;
          }
            
        if( _max > Max )
          {
          Max = _max;
          }
      }
    };

    struct PerCell 
    {
      float        min[3];
      float        max[3];
      unsigned int ind;
    };

    struct CenterOrder
    {
      unsigned int d;
    
      CenterOrder( unsigned int _d ) : 
          d(_d)
      {
      }

      bool operator()( const PerCell& pc0, const PerCell& pc1 )
      {
          return (pc0.min[d] + pc0.max[d]) < (pc1.min[d] + pc1.max[d]);
      }
    };

    struct LeftPredicate
    {
      unsigned int       d;
      float              p;
    
      LeftPredicate( unsigned int _d, float _p ) : 
          d(_d), p(2.0f*_p)
      {
      }
   
      bool operator()( const PerCell& pc )
      {
          return (pc.min[d] + pc.max[d]) < p;
      }
    };


    // -------------------------------------------------------------------------

    void FindMinMax( const PerCell* begin, const PerCell* end,  
                       float* min, float* max )
    {
      if( begin == end )
        {
        return;
        }
            
      for( unsigned int d=0; d<3; ++d )
        {
        min[d] = begin->min[d];
        max[d] = begin->max[d];
        }
        
      while( ++begin != end )
        {
        for( unsigned int d=0; d<3; ++d )
          {
          if( begin->min[d] < min[d] )    min[d] = begin->min[d];
          if( begin->max[d] > max[d] )    max[d] = begin->max[d];
          }
        }
    }

    // -------------------------------------------------------------------------
    
    void FindMinD( const PerCell* begin, const PerCell* end,  
                     unsigned int d, float& min )
    {
      min = begin->min[d];
        
      while( ++begin != end )
        {
        if( begin->min[d] < min )
          {    
          min = begin->min[d];
          }
        }
    }

    void FindMaxD( const PerCell* begin, const PerCell* end,  
                     unsigned int d, float& max )
    {
      max = begin->max[d];
        
      while( ++begin != end )
        {
        if( begin->max[d] > max )    
          {
          max = begin->max[d];
          }
        }
    }

    // -------------------------------------------------------------------------

    void Split( unsigned int index, float min[3], float max[3] )
    {
      unsigned int start = m_nodes[index].Start();
      unsigned int size  = m_nodes[index].Size();
        
      if( size < m_leafsize )
        {
        return;
        }

      PerCell* begin = m_pc + start;
      PerCell* end   = m_pc + start + size;
      PerCell* mid = begin;

      const int nbuckets = 6;

      const float ext[3] = { max[0]-min[0], max[1]-min[1], max[2]-min[2] };
      const float iext[3] = { nbuckets/ext[0], nbuckets/ext[1], nbuckets/ext[2] };

      Bucket b[3][nbuckets];
            
      for( const PerCell* pc=begin; pc!=end; ++pc )
        {
        for( unsigned int d=0; d<3; ++d )
          {
          float cen = (pc->min[d] + pc->max[d])/2.0f;
          int   ind = (int)( (cen-min[d])*iext[d] );

          if( ind<0 )
            {
            ind = 0;
            }

          if( ind>=nbuckets )
            {
            ind = nbuckets-1;
            }

          b[d][ind].Add( pc->min[d], pc->max[d] );
          }
        }
        
      float cost = std::numeric_limits<float>::max();
      float plane;
      unsigned int dim;

      for( unsigned int d=0; d<3; ++d )    
        {
        unsigned int sum = 0;
            
        for( unsigned int n=0; n<nbuckets-1; ++n )
          {
          float lmax = -std::numeric_limits<float>::max();
          float rmin =  std::numeric_limits<float>::max();

          for( unsigned int m=0; m<=n; ++m )
            {
            if( b[d][m].Max > lmax )
              {
              lmax = b[d][m].Max;
              }
            }
                
          for( unsigned int m=n+1; m<nbuckets; ++m )
            {
            if( b[d][m].Min < rmin )
              {
              rmin = b[d][m].Min;
              }
            }
                
          sum += b[d][n].Cnt;
                
          float lvol = (lmax-min[d])/ext[d];
          float rvol = (max[d]-rmin)/ext[d];
                
          float c = lvol*sum + rvol*(size-sum);
                
          if( sum > 0 && sum < size && c < cost )
            {
            cost    = c;
            dim     = d;
            plane   = min[d] + (n+1)/iext[d];
            }
          }
        }

      if( cost != std::numeric_limits<float>::max() )
        {
        mid = std::partition( begin, end, LeftPredicate( dim, plane ) );
        }

      // fallback
      if( mid == begin || mid == end )
        {
        dim = std::max_element( ext, ext+3 ) - ext;

        mid = begin + (end-begin)/2;
        std::nth_element( begin, mid, end, CenterOrder( dim ) );
        }

      float lmin[3], lmax[3], rmin[3], rmax[3];

      FindMinMax( begin, mid, lmin, lmax );
      FindMinMax( mid,   end, rmin, rmax );

      float clip[4] = { lmax[dim], rmin[dim], lmin[dim], rmax[dim] }; // lmin and rmax were also added to get the bounding bnox

      CellTreeNode child[2];
      child[0].MakeLeaf( begin - m_pc, mid-begin );
      child[1].MakeLeaf( mid   - m_pc, end-mid );
        
      m_nodes[index].MakeNode( m_nodes.size(), dim, clip );
      m_nodes.insert( m_nodes.end(), child, child+2 );

      Split( m_nodes[index].GetLeftChildIndex(), lmin, lmax );
      Split( m_nodes[index].GetRightChildIndex(), rmin, rmax );
    }
     
  public:

    CellTreeBuilder()
    {
      m_buckets =  5;
      m_leafsize = 8;
    }
    
    void Build( CellTree& ct, vtkDataSet* ds )
    {
      const vtkIdType size = ds->GetNumberOfCells();
      assert( size <= std::numeric_limits<unsigned int>::max() );

      m_pc = new PerCell[size];

      float min[3] = { 
          std::numeric_limits<float>::max(),
          std::numeric_limits<float>::max(),
          std::numeric_limits<float>::max()
      };

      float max[3] = { 
          -std::numeric_limits<float>::max(),
          -std::numeric_limits<float>::max(),
          -std::numeric_limits<float>::max(),
      };
        
        

      for( unsigned int i=0; i<size; ++i )
        {
        m_pc[i].ind = i;

        double bounds[6];
        ds->GetCellBounds( i, bounds );
            
        for( int d=0; d<3; ++d )
          {
          m_pc[i].min[d] = bounds[2*d+0];
          m_pc[i].max[d] = bounds[2*d+1];

          if( m_pc[i].min[d] < min[d] )
            {
            min[d] = m_pc[i].min[d];
            }

          if( m_pc[i].max[d] > max[d] )  /// This can be m_pc[i].max[d] instead of min
            {
            max[d] = m_pc[i].max[d];
            }
          }
        }

      ct.DataBBox[0] = min[0];
      ct.DataBBox[1] = max[0];
      ct.DataBBox[2] = min[1];
      ct.DataBBox[3] = max[1];
      ct.DataBBox[4] = min[2];
      ct.DataBBox[5] = max[2];
                
      CellTreeNode root;
      root.MakeLeaf( 0, size );
      m_nodes.push_back( root );

      Split( 0, min, max );
        
      ct.Nodes.resize( m_nodes.size() );
      ct.Nodes[0] = m_nodes[0];

      std::vector<CellTreeNode>::iterator ni = ct.Nodes.begin();
      std::vector<CellTreeNode>::iterator nn = ct.Nodes.begin()+1;

      for( ; ni!=ct.Nodes.end(); ++ni )
        {
        if( ni->IsLeaf() )
          {
          continue;
          }
            
        *(nn++) = m_nodes[ni->GetLeftChildIndex()];
        *(nn++) = m_nodes[ni->GetRightChildIndex()];

        ni->SetChildren( nn-ct.Nodes.begin()-2 );
        }

      // --- final 
        
      ct.Leaves.resize( size );

      for( int i=0; i<size; ++i )
        {
        ct.Leaves[i] = m_pc[i].ind;
        }
        
      delete[] m_pc;
    }

  public:

    unsigned int                 m_buckets;
    unsigned int                 m_leafsize;
    PerCell*                     m_pc;
    std::vector<CellTreeNode>    m_nodes;
};

// ---------------------------------------------------------------------------

typedef std::stack<CellTreeNode*, std::vector<CellTreeNode*> > nodestack;

vtkCellTreeLocator::vtkCellTreeLocator( ) 
{
    
  this->MaxCellsPerLeaf = 8;
  this->NumberOfBuckets = 5;

  this->CellArray = NULL;
  this->Locations = NULL;

  this->Tree = NULL;

}

// ---------------------------------------------------------------------------

vtkCellTreeLocator::~vtkCellTreeLocator()
{
  Free();
}



void vtkCellTreeLocator::Build()
{
  Free();

  vtkIdType numCells;

  if( !this->DataSet || (numCells = this->DataSet->GetNumberOfCells() < 1) )
    {
    std::cerr << " No Cells in the data set\n";
    return;
    }

  this->Tree = new CellTree;

  CellTreeBuilder builder;

  builder.m_leafsize = MaxCellsPerLeaf;
  builder.m_buckets  = NumberOfBuckets;
  builder.Build( *(Tree), this->DataSet );
}

void vtkCellTreeLocator::BuildLocator()
{
  Build();
}

// ---------------------------------------------------------------------------

void vtkCellTreeLocator::Free()
{
  if( Tree )
    {
    delete Tree;
    Tree = NULL;
    }
}

// ---------------------------------------------------------------------------

vtkIdType vtkCellTreeLocator::FindCell( double pos[3], double vtkNotUsed, vtkGenericCell *cell,  double pcoords[3], 
                                       double* weights )
{
  if( Tree == 0 )
    {
    return -1;
    }

  double closestPoint[3], dist2;
  int subId;

  const float _pos[3] = { pos[0], pos[1], pos[2] };
  PointTraversal pt( *(this->Tree), _pos );

  bool found = false;

  while( const CellTreeNode* n = pt.Next() )
    {
    const unsigned int* begin = &(this->Tree->Leaves[n->Start()]);
    const unsigned int* end   = begin + n->Size();
            
    for( ; begin!=end; ++begin )
      {
      this->DataSet->GetCell(*begin, cell);
      if( cell->EvaluatePosition(pos, closestPoint, subId, pcoords, dist2, weights)==1 )
        {
        return *begin;
        }
      }
    }

  return -1;
}

// ---------------------------------------------------------------------------

double _getMinDistPOS_X(const double origin[3], const double dir[3], const double B[6])
{
  return ((B[0] - origin[0]) / dir[0]);
}
double _getMinDistNEG_X(const double origin[3], const double dir[3], const double B[6])
{
  return ((B[1] - origin[0]) / dir[0]);
}
double _getMinDistPOS_Y(const double origin[3], const double dir[3], const double B[6])
{
  return ((B[2] - origin[1]) / dir[1]);
}
double _getMinDistNEG_Y(const double origin[3], const double dir[3], const double B[6])
{
  return ((B[3] - origin[1]) / dir[1]);
}
double _getMinDistPOS_Z(const double origin[3], const double dir[3], const double B[6])
{
  return ((B[4] - origin[2]) / dir[2]);
}
double _getMinDistNEG_Z(const double origin[3], const double dir[3], const double B[6])
{
  return ((B[5] - origin[2]) / dir[2]);
}

typedef std::pair<double, int> Intersection;

int vtkCellTreeLocator::IntersectWithLine(double p1[3], double p2[3], double tol,
                                      double& t, double x[3], double pcoords[3],
                                      int &subId, vtkIdType &cellIds,
                                      vtkGenericCell *cell)
{
  //
  CellTreeNode  *node, *Near, *Far;
  double    ctmin, ctmax, tmin, tmax, _tmin, _tmax, tDist;
  double    ray_vec[3] = { p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2] };

  double *boundsPtr;
  double cellBounds[6]; 
  
  //
 // this->BuildLocator();
  //
  // Does ray pass through root BBox
  tmin = 0; tmax = 1;
  
  if (!RayMinMaxT(p1, ray_vec, tmin, tmax)) // 0 for root node
    {
    return false;
    }
  // Ok, setup a stack and various params
  nodestack  ns;
  double    closest_intersection = VTK_LARGE_FLOAT;
  bool     HIT = false;
  // setup our axis optimized ray box edge stuff
  int axis = getDominantAxis(ray_vec);
  double (*_getMinDist)(const double origin[3], const double dir[3], const double B[6]);
  switch (axis) 
    {
    case POS_X: _getMinDist = _getMinDistPOS_X; break;
    case NEG_X: _getMinDist = _getMinDistNEG_X; break;
    case POS_Y: _getMinDist = _getMinDistPOS_Y; break;
    case NEG_Y: _getMinDist = _getMinDistNEG_Y; break;
    case POS_Z: _getMinDist = _getMinDistPOS_Z; break;
    default:    _getMinDist = _getMinDistNEG_Z; break;
    }

  //
  // we will sort intersections by t, so keep track using these lists
  //
  std::vector<Intersection> t_list;
  vtkSmartPointer<vtkPoints> tempPoints;
  vtkSmartPointer<vtkIdList>    tempIds;
 /* if (points)
    {
    tempPoints = vtkSmartPointer<vtkPoints>::New();
    }*/
  if (cellIds)
    {
    tempIds = vtkSmartPointer<vtkIdList>::New();
    }
  int icount = 0;
  //
  // OK, lets walk the tree and find intersections
  //
  CellTreeNode* n = &this->Tree->Nodes.front();
  CellTreeNode* n1 = &this->Tree->Nodes.at(1);
  CellTreeNode* n2 = &this->Tree->Nodes.at(2);
  ns.push(n2);   // Has to push 3 elements since there is three dimensions and a node stores only a single dimension
  ns.push(n1);
  ns.push(n);
  while (!ns.empty())
    {
    node = ns.top();
    ns.pop();
    // We do as few tests on the way down as possible, because our BBoxes
    // can be quite tight and we want to reject as many boxes as possible without
    // testing them at all - mainly because we quickly get to a leaf node and
    // test candidates, once we've found a hit, we note the intersection t val,
    // as soon as we pull a BBox of the stack that has a closest point further
    // than the t val, we know we can stop.

    int mustCheck = 0;  // to check if both left and right sub trees need to be checked
  
    //
    while (!node->IsLeaf())
      { // this must be a parent node
      // Which child node is closest to ray origin - given direction
      Classify(p1, ray_vec, tDist, Near, node, Far, mustCheck);
      // if the distance to the far edge of the near box is > tmax, no need to test far box
      // (we still need to test Mid because it may overlap slightly)
      if(mustCheck)
        {
        ns.push(Far);
        node = Near;
        }
      else if ((tDist > tmax) || (tDist <= 0) )
        { //<=0 for ray on edge
        node = Near;
        }
      // if the distance to the far edge of the near box is < tmin, no need to test near box
      else if (tDist < tmin)
        { 
        ns.push(Near);  
        node = Far;
        }
      // All the child nodes may be candidates, keep near, push far then mid
      else
        {
        ns.push(Far);
        node = Near;
        }
      }
    double t_hit, ipt[3];
    // Ok, so we're a leaf node, first check the BBox against the ray
    // then test the candidates in our sorted ray direction order
    _tmin = tmin; _tmax = tmax;
//    if (node->RayMinMaxT(p1, ray_vec, _tmin, _tmax)) {
    // Was the closest point on the box > intersection point
//      if (_tmax>closest_intersection) break;
    //
    for (int i=0; i<node->Size(); i++)
      {
      vtkIdType cell_ID = this->Tree->Leaves[node->Start()+i];
      //
      
      boundsPtr = cellBounds;
      if (this->CellBounds)
      {
      boundsPtr = this->CellBounds[cell_ID];
      }
    else
      {
      this->DataSet->GetCellBounds(cell_ID, cellBounds);
      }
      if (_getMinDist(p1, ray_vec, cellBounds) > closest_intersection)
        {
        break;
        }
      //
      ctmin = _tmin; ctmax = _tmax;
      if (RayMinMaxT(cellBounds, p1, ray_vec, ctmin, ctmax))
        {
        if (this->IntersectCellInternal(cell_ID, p1, p2, tol, t_hit, ipt, pcoords, subId))
          {
          /*if (points)
            {
            tempPoints->InsertNextPoint(ipt);
            }*/
          if (cellIds)
            {
            tempIds->InsertNextId(cell_ID);
            }
          t_list.push_back(Intersection(t_hit, icount++));
          HIT = true;
          if(HIT)
          {
          return 1;
          }
          }
        }
      }
//    }
    }
  if (HIT)
    {
  //  std::sort(t_list.begin(), t_list.end(), Isort());
    int N = static_cast<int>(t_list.size());
    /*if (points)
      {
      points->SetNumberOfPoints(N);
      }*/
   /* if (cellIds)
      {
      cellIds->SetNumberOfIds(N);
      }*/
    //for (int n=0; n<N; n++)
    //  {
    //  Intersection &i = t_list[n];
    //  /*if (points)
    //    {
    //    points->SetPoint(n, tempPoints->GetPoint(i.second));
    //    }*/
    //  if (cellIds)
    //    {
    //    cellIds->SetId(n, tempIds->GetId(i.second));
    //    }
    //  }
    }
  //
  return HIT;

}

bool vtkCellTreeLocator::RayMinMaxT(const double origin[3],
                         const double dir[3],
                         double &rTmin,
                         double &rTmax) 
{
  double tT;
  // X-Axis
  float bounds[6]; 
  
  bounds[0] = this->Tree->DataBBox[0];
  bounds[1] = this->Tree->DataBBox[1];
  bounds[2] = this->Tree->DataBBox[2];
  bounds[3] = this->Tree->DataBBox[3];
  bounds[4] = this->Tree->DataBBox[4];
  bounds[5] = this->Tree->DataBBox[5];

  if (dir[0] < -Epsilon_)
    {                // ray travelling in -x direction
    tT = (bounds[0] - origin[0]) / dir[0];
    if (tT < rTmin)
      {
      return (false);        // ray already left of box. Can't hit
      }
    if (tT <= rTmax)
      {
      rTmax = tT;           // update new tmax
      }
    tT = (bounds[1] - origin[0]) / dir[0]; // distance to right edge
    if (tT >= rTmin)
      {                     // can't see this ever happening
      if (tT > rTmax)
        {
        return false;        // clip start of ray to right edge
        }
      rTmin = tT;
      }
    }
  else if (dir[0] > Epsilon_)
    {
    tT = (bounds[1] - origin[0]) / dir[0];
    if (tT < rTmin)
      {
      return (false);
      }
    if (tT <= rTmax)
      {
      rTmax = tT;
      }
    tT = (bounds[0] - origin[0]) / dir[0];
    if (tT >= rTmin)
      {
      if (tT > rTmax)
        {
        return (false);
        }
      rTmin = tT;
      }
    }
  else if (origin[0] < bounds[0] || origin[0] > bounds[1])
    {
    return (false);
    }
  // Y-Axis
  if (dir[1] < -Epsilon_)
    {
    tT = (bounds[2] - origin[1]) / dir[1];
    if (tT < rTmin)
      {
      return (false);
      }
    if (tT <= rTmax)
      {
      rTmax = tT;
      }
    tT = (bounds[3] - origin[1]) / dir[1];
    if (tT >= rTmin)
      {
      if (tT > rTmax)
        {
        return (false);
        }
      rTmin = tT;
      }
    }
  else if (dir[1] > Epsilon_)
    {
    tT = (bounds[3] - origin[1]) / dir[1];
    if (tT < rTmin)
      {
      return (false);
      }
    if (tT <= rTmax)
      {
      rTmax = tT;
      }
    tT = (bounds[2] - origin[1]) / dir[1];
    if (tT >= rTmin)
      {
      if (tT > rTmax)
        {
        return (false);
        }
      rTmin = tT;
      }
    }
  else if (origin[1] < bounds[2] || origin[1] > bounds[3])
    {
    return (false);
    }
  // Z-Axis
  if (dir[2] < -Epsilon_)
    {
    tT = (bounds[4] - origin[2]) / dir[2];
    if (tT < rTmin)
      {
      return (false);
      }
    if (tT <= rTmax)
      {
      rTmax = tT;
      }
    tT = (bounds[5] - origin[2]) / dir[2];
    if (tT >= rTmin)
      {
      if (tT > rTmax)
        {
        return (false);
        }
      rTmin = tT;
      }
    }
  else if (dir[2] > Epsilon_)
    {
    tT = (bounds[5] - origin[2]) / dir[2];
    if (tT < rTmin)
      {
      return (false);
      }
    if (tT <= rTmax)
      {
      rTmax = tT;
      }
    tT = (bounds[4] - origin[2]) / dir[2];
    if (tT >= rTmin)
      {
      if (tT > rTmax)
        {
        return (false);
        }
      rTmin = tT;
      }
    }
  else if (origin[2] < bounds[4] || origin[2] > bounds[5])
    {
    return (false);
    }
  return (true);
}

bool vtkCellTreeLocator::RayMinMaxT(const double bounds[6],
                         const double origin[3],
                         const double dir[3],
                         double &rTmin,
                         double &rTmax)
{
  double tT;
  // X-Axis
  if (dir[0] < -Epsilon_)
    {                // ray travelling in -x direction
    tT = (bounds[0] - origin[0]) / dir[0]; // Ipoint less than minT - ray outside box!
    if (tT < rTmin)
      {
      return (false);
      }
    if (tT <= rTmax)
      {
      rTmax = tT;           // update new tmax
      }
    tT = (bounds[1] - origin[0]) / dir[0]; // distance to right edge
    if (tT >= rTmin)
      {                     // can't see this ever happening
      if (tT > rTmax)
        {
        return false;        // clip start of ray to right edge
        }
      rTmin = tT;
      }
    }
  else if (dir[0] > Epsilon_)
    {
    tT = (bounds[1] - origin[0]) / dir[0];
    if (tT < rTmin)
      {
      return (false);
      }
    if (tT <= rTmax)
      {
      rTmax = tT;
      }
    tT = (bounds[0] - origin[0]) / dir[0];
    if (tT >= rTmin)
      {
      if (tT > rTmax)
        {
        return (false);
        }
      rTmin = tT;
      }
    }
  else if (origin[0] < bounds[0] || origin[0] > bounds[1])
    {
    return (false);
    }
  // Y-Axis
  if (dir[1] < -Epsilon_)
    {
    tT = (bounds[2] - origin[1]) / dir[1];
    if (tT < rTmin)
      {
      return (false);
      }
    if (tT <= rTmax) rTmax = tT;
    tT = (bounds[3] - origin[1]) / dir[1];
    if (tT >= rTmin)
      {
      if (tT > rTmax)
        {
        return (false);
        }
      rTmin = tT;
      }
    }
  else if (dir[1] > Epsilon_)
    {
    tT = (bounds[3] - origin[1]) / dir[1];
    if (tT < rTmin)
      {
      return (false);
      }
    if (tT <= rTmax)
      {
      rTmax = tT;
      }
    tT = (bounds[2] - origin[1]) / dir[1];
    if (tT >= rTmin)
      {
      if (tT > rTmax)
        {
        return (false);
        }
      rTmin = tT;
      }
    }
  else if (origin[1] < bounds[2] || origin[1] > bounds[3])
    {
    return (false);
    }
  // Z-Axis
  if (dir[2] < -Epsilon_)
    {
    tT = (bounds[4] - origin[2]) / dir[2];
    if (tT < rTmin)
      {
      return (false);
      }
    if (tT <= rTmax)
      {
      rTmax = tT;
      }
    tT = (bounds[5] - origin[2]) / dir[2];
    if (tT >= rTmin)
      {
      if (tT > rTmax)
        {
        return (false);
        }
      rTmin = tT;
      }
    }
  else if (dir[2] > Epsilon_)
    {
    tT = (bounds[5] - origin[2]) / dir[2];
    if (tT < rTmin)
      {
      return (false);
      }
    if (tT <= rTmax)
      {
      rTmax = tT;
      }
    tT = (bounds[4] - origin[2]) / dir[2];
    if (tT >= rTmin)
      {
      if (tT > rTmax)
        {
        return (false);
        }
      rTmin = tT;
      }
    }
  else if (origin[2] < bounds[4] || origin[2] > bounds[5])
    {
    return (false);
    }
  return (true);
}

int vtkCellTreeLocator::getDominantAxis(const double dir[3])
{
  double tX = (dir[0]>0) ? dir[0] : -dir[0];
  double tY = (dir[1]>0) ? dir[1] : -dir[1];
  double tZ = (dir[2]>0) ? dir[2] : -dir[2];
  if (tX > tY && tX > tZ)
    {
    return ((dir[0] > 0) ? POS_X : NEG_X);
    }
  else if ( tY > tZ )
    {
    return ((dir[1] > 0) ? POS_Y : NEG_Y);
    }
  else
    {
    return ((dir[2] > 0) ? POS_Z : NEG_Z);
    }
}

void vtkCellTreeLocator::Classify(const double origin[3],
                       const double dir[3],
                       double &rDist,
                       CellTreeNode *&Near, CellTreeNode *&Parent,
                       CellTreeNode *&Far, int& mustCheck)
{
  double tOriginToDivPlane = Parent->GetLeftMaxValue() - origin[Parent->GetDimension()];
  double tOriginToDivPlane2 = Parent->GetRightMinValue() - origin[Parent->GetDimension()];
  double tDivDirection   = dir[Parent->GetDimension()];
  if ( tOriginToDivPlane2 > 0 )  // origin is right of the rmin
    {
    Near = &this->Tree->Nodes.at(Parent->GetLeftChildIndex());
    Far  = &this->Tree->Nodes.at(Parent->GetLeftChildIndex()+1);
    rDist = (tDivDirection) ? tOriginToDivPlane2 / tDivDirection : VTK_LARGE_FLOAT;
    }
  else if (tOriginToDivPlane < 0)  // origin is left of the lm
    {
    Far  = &this->Tree->Nodes.at(Parent->GetLeftChildIndex());
    Near = &this->Tree->Nodes.at(Parent->GetLeftChildIndex()+1);
    rDist = (tDivDirection) ? tOriginToDivPlane / tDivDirection : VTK_LARGE_FLOAT;
    }
  // Ray was exactly on edge of box, check direction

  else 
    {
    if(tOriginToDivPlane > 0 && tOriginToDivPlane2 < 0)
     {
     mustCheck = 1;
     }

    if ( tDivDirection < 0)
      {
      Near = &this->Tree->Nodes.at(Parent->GetLeftChildIndex());
      Far  = &this->Tree->Nodes.at(Parent->GetLeftChildIndex()+1);
      if(!(tOriginToDivPlane > 0 || tOriginToDivPlane < 0))
        {
        mustCheck=1;
        }
      rDist = (tDivDirection) ? 0 / tDivDirection : VTK_LARGE_FLOAT;
      }
    else
      {
      Far  = &this->Tree->Nodes.at(Parent->GetLeftChildIndex());
      Near = &this->Tree->Nodes.at(Parent->GetLeftChildIndex()+1);
      if(!(tOriginToDivPlane2 > 0 || tOriginToDivPlane2 < 0))
        {
        mustCheck=1;
        }
      rDist = (tDivDirection) ? 0 / tDivDirection : VTK_LARGE_FLOAT;
      }
    }
  
}

int vtkCellTreeLocator::IntersectCellInternal(
  vtkIdType cell_ID,
  const double p1[3],
  const double p2[3],
  const double tol,
  double &t,
  double ipt[3],
  double pcoords[3],
  int &subId)
{
  this->DataSet->GetCell(cell_ID, this->GenericCell);
  return this->GenericCell->IntersectWithLine(const_cast<double*>(p1), const_cast<double*>(p2), tol, t, ipt, pcoords, subId);
}

void vtkCellTreeLocator::FreeSearchStructure(void)
{
 
}

void vtkCellTreeLocator::GenerateRepresentation(int level, vtkPolyData *pd)
{

}