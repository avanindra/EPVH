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
#include "opencvincludes.h"
#include "silhouette.h"

class CameraInfo
{
  

public:

  std::vector< Eigen::Matrix< double , 3 , 4 > > mProjectionMatrices;
  std::vector< Eigen::Matrix3d > mInvProjectionMatrices;
  std::vector< Eigen::Vector3d > mCameraCenters;

  std::vector< cv::Size > mImageDims;

  std::vector< tr::Silhouette > &mObjectSilhouettes;
  std::vector< std::string > mObjectSilhouettePaths;
  

  public:

	  CameraInfo( std::vector< tr::Silhouette >& silhouettes  );

    int get_image_width( int camId );
    int get_image_height( int camId );

    void projectPoint( Eigen::Vector3d objectPoint , int camId , Eigen::Vector2d& projection );

	void projectHomogeneous(Eigen::Vector4d objectPoint, int camId, Eigen::Vector3d& projection);

    void getRay( int camId , const Eigen::Vector2d& imagePos , Eigen::Vector3d& rayDir  );

	void getRay(int camId, const cv::Point2f& imagePos, Eigen::Vector3d& rayDir);

	int get_cam_count();

	tr::Silhouette& silhouette(int camId);

	void get_camera_center(int camId, Eigen::Vector3d& camCenter);

    
    ~CameraInfo();
};

#endif // CAMERAINFO_H
