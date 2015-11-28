
#ifndef SORTINGPREDICATES_H
#define SORTINGPREDICATES_H

#include "iostream"
#include "cmath"

namespace tr{
  
  
typedef std::pair< int , float > CameraScorePair;

inline bool cameraScorePairPredicate( CameraScorePair pair1 , CameraScorePair pair2 )
{  
  return pair1.second > pair2.second;  
}


typedef std::pair< int , double > VertexDistancePair;

inline bool vertexDistancePairIncreasingPredicate( VertexDistancePair pair1 , VertexDistancePair pair2 )
{
  return pair1.second < pair2.second;
}

inline bool edgeDistancePredicate( std::pair< int , float > pair1 , std::pair< int , float > pair2 )
{

  return pair1.second < pair2.second;
}

inline bool pointSlopePredicate( std::pair< int , float > pair1 , std::pair< int , float > pair2 )
{
  return pair1.second < pair2.second;
}


}


#endif
