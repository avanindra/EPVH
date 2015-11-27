/*
 * Copyright 2015 <copyright holder> <email>
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#ifndef CAMERAINFO_H
#define CAMERAINFO_H

#include "vector"
#include "eigenincludes.h"

class CameraInfo
{
  
  std::vector< Eigen::Matrix< float , 3 , 4 > > mProjectionMatrices;
  std::vector< Eigen::Matrix3f > mInvProjectionMatrices;
  std::vector< Eigen::Vector3f > mCameraCenters;

  std::vector< cv::Size > mImageDims;
  

  public:

    CameraInfo();

    int get_image_width( int camId );
    int get_image_height( int camId );

    void projectPoint();

    void getRay( int camId , const Eigen::Vector2d& imagePos , Eigen::Vector3d& rayDir  );

    
    ~CameraInfo();
};

#endif // CAMERAINFO_H
