#ifndef TR_TYPES_H
#define TR_TYPES_H

#define _WIN32_WINNT 0x0501
#include "Eigen/Dense"

namespace tr
{
  
  typedef Eigen::Vector2d Vector2;
  typedef Eigen::Vector2d Point2;
  typedef Eigen::Vector2d Point2d;
  typedef Eigen::Vector2f Vector2f;
  typedef Eigen::Vector2d Vector2d;
  typedef Eigen::Matrix< float , 3 , 4 > ProjectionMatrix;
  typedef Eigen::Matrix< float , 3 , 4 > ProjectionMatrixf;
  typedef Eigen::Matrix< double , 3 , 4 > ProjectionMatrixd;
  typedef Eigen::Matrix< double , 3 , 3 > Matrix33d;
  typedef Eigen::Vector3d Vector3d;
  typedef Eigen::Vector3f Vector3f;
  typedef Eigen::Vector2f Point2f;
  typedef Eigen::Vector3d Point3d;
  typedef Eigen::Vector3f Point3f;
  typedef Eigen::Vector4d Point4d;
  typedef Eigen::Vector4d Vector4d;
  typedef Eigen::Vector4f Vector4f;
  typedef Eigen::Vector4f Point4f;  
  
  
  enum{ SP_PR_BG = 0 , SP_PR_FG = 1 , SP_BG = 3 , SP_FG = 4  } ;
  enum{ PR_BG = 0 , PR_FG = 1 , BG = 3 , FG = 4  } ;
  
}

#endif