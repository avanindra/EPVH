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

#include "camerainfo.h"


CameraInfo::CameraInfo() : mObjectSilhouettes(mDummySilhouettes)
{

}

CameraInfo::CameraInfo(std::vector< tr::Silhouette >& silhouettes) : mObjectSilhouettes( silhouettes )
{

}



int CameraInfo::get_image_width( int camId )
{
  return mImageDims[ camId ].width;
}


int CameraInfo::get_image_height( int camId )
{

  return mImageDims[ camId ].height;
}

void CameraInfo::projectPoint(Eigen::Vector3d objectPoint, int camId, Eigen::Vector2d &projection)
{
  Eigen::Vector4d pos4;

  pos4.block( 0 , 0 , 3 , 1 ) = objectPoint;

  pos4( 3 ) = 1;

  Eigen::Vector3d proj3 = mProjectionMatrices[ camId ] * pos4;

  projection( 0 ) = proj3( 0 ) / proj3( 2 );
  projection( 1 ) = proj3( 1 ) / proj3( 2 );

}

void CameraInfo::projectHomogeneous(Eigen::Vector4d objectPoint, int camId, Eigen::Vector3d& projection)
{

	projection = mProjectionMatrices[camId] * objectPoint;
}

void CameraInfo::getRay(int camId, const Eigen::Vector2d &imagePos, Eigen::Vector3d &rayDir)
{
  Eigen::Vector3d pos3;

  pos3( 0 ) = imagePos( 0 );
  pos3( 1 ) = imagePos( 1 );
  pos3( 2 ) = 1;

  rayDir = mInvProjectionMatrices[ camId ] * pos3;

  rayDir.normalize();
}

void CameraInfo::getRay(int camId, const cv::Point2f& imagePos, Eigen::Vector3d& rayDir)
{
	Eigen::Vector3d pos3;

	pos3(0) = imagePos.x;
	pos3(1) = imagePos.y;
	pos3(2) = 1;

	rayDir = mInvProjectionMatrices[camId] * pos3;

	rayDir.normalize();
}


int CameraInfo::get_cam_count()
{
	return mProjectionMatrices.size();
}

void CameraInfo::get_camera_center(int camId, Eigen::Vector3d& camCenter)
{
	camCenter = mCameraCenters[camId];
}


tr::Silhouette& CameraInfo::silhouette(int camId)
{
	return mObjectSilhouettes[camId];
}

CameraInfo::~CameraInfo()
{

}
