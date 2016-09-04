/*
 * Copyright (c) 2012, avanindra <avanindra.singh@gmail.com>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *     names of its contributors may be used to endorse or promote products
 *     derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY avanindra <email> ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL avanindra <email> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */


#include "contourio.h"
#include "fstream"
#include "iostream"
#include "epvh.h"
#include "camerainfo.h"
#include "opencvincludes.h"
#include "display3droutines.h"
#include "vtkTriangleFilter.h"

#ifdef WIN32
#include "conio.h"
#endif

void readContours( std::string& contourFilePath , std::vector< std::vector< Eigen::Vector2d > >& contours );

void readCameras(std::string& cameraFilePath, std::vector< Eigen::Matrix< double, 3, 4 > >& projectionMatrices);

int main( int argc , char **argv )
{
  
  tr::ContourIO contourReader;

  std::string contourFilePath = DATASET_DIR + std::string("/alien-contour.txt");
  std::string cameraFilepPath = DATASET_DIR + std::string("/alien-camera.txt");

  std::vector< std::vector< Eigen::Vector2d > > contours;
  std::vector< Eigen::Matrix< double, 3, 4 > > projectionMatrices;

  readContours( contourFilePath , contours );
  readCameras(cameraFilepPath, projectionMatrices);

  std::cout << contours.size() << " " << projectionMatrices.size() << std::endl;

  CameraInfo cameraInfo;
  tr::EPVH epvh;

  epvh.loadCameraInfo(&cameraInfo);

  cameraInfo.mProjectionMatrices = projectionMatrices;
  int numCams = projectionMatrices.size();

  cameraInfo.mCameraCenters.resize(numCams);
  cameraInfo.mInvProjectionMatrices.resize(numCams);
  cameraInfo.mImageDims.resize(numCams);

  std::vector< int > silhouetteCameras;

  

  int skip = 0;

  for (int cc = 0; cc < numCams; cc++)
  {
	  Eigen::Vector3d camCenter = -projectionMatrices[cc].block(0, 0, 3, 3).inverse() * projectionMatrices[cc].block(0, 3, 3, 1);

	  cameraInfo.mCameraCenters[cc] = camCenter;

	  cameraInfo.mInvProjectionMatrices[cc] = projectionMatrices[cc].block(0, 0, 3, 3).inverse();

	  std::string imagePath;

	  if ( cc < 10 )
	      imagePath = std::string( DATASET_DIR ) + "/alien-images/000" + std::to_string(cc) + ".jpg";
	  else
		  imagePath = std::string(DATASET_DIR) + "/alien-images/00" + std::to_string(cc) + ".jpg";

	  cv::Mat im = cv::imread(imagePath);

	  cameraInfo.mImageDims[cc] = im.size();

	  //std::cout << im.size() << std::endl;


	  //if ( skip == 0 )
	  silhouetteCameras.push_back( cc );

	  //skip++;

	  //if (skip == 3)
	  //{
		 // skip = 0;
	  //}

  } 

  std::vector< Eigen::Vector3d > colors(numCams, Eigen::Vector3d(0, 1, 0));

  //tr::Display3DRoutines::displayPointSet(cameraInfo.mCameraCenters, colors);

  std::cout << " number of silhouette cams : " << silhouetteCameras.size() << std::endl;

  cameraInfo.mContours = contours;

  epvh.setSilhouetteCameras(silhouetteCameras);

  epvh.compute();

  vtkSmartPointer< vtkPolyData > polygons = vtkSmartPointer< vtkPolyData >::New();

  epvh.generatePolygons(polygons);

  std::cout<< polygons->GetNumberOfPoints() <<"  "<<polygons->GetNumberOfCells() <<std::endl;

  vtkSmartPointer< vtkTriangleFilter > triangulator = vtkSmartPointer< vtkTriangleFilter >::New();

  triangulator->SetInputData(polygons);

  triangulator->Update();

  polygons->DeepCopy( triangulator->GetOutput() );

  tr::Display3DRoutines::displayPolyData(polygons);


#ifdef WIN32
  getch();
#endif
  
  return 0;
}


void readContours( std::string& contourFilePath , std::vector< std::vector< Eigen::Vector2d > >& contours )
{
	std::ifstream fs( contourFilePath );

	int numLines = 0;

	if (fs.is_open())
	{
		std::string str;

		

		while (!fs.eof())
		{
			std::string line1, line2;

			if (std::getline(fs, line1) && std::getline(fs, line2) && std::getline(fs, str))
			{

				int camId, numPoints;

				std::stringstream sstr1(line1);

				sstr1 >> camId >> numPoints;

				std::stringstream sstr2(line2);

				std::vector< Eigen::Vector2d > contour(numPoints);

				for ( int pp = 0; pp < numPoints; pp++ )
				{
					sstr2 >> contour[pp](0) >> contour[pp](1);
				}


				contours.push_back( contour );
			}
			else
			{
				break;
			}

		}
	}
	else
	{
		std::cout << " error opening file "<< contourFilePath << std::endl;
	}


	
	std::cout << " reading completed " <<numLines<< std::endl;

}


void readCameras( std::string& cameraFilePath , std::vector< Eigen::Matrix< double , 3 , 4 > >& projectionMatrices )
{
	std::ifstream fs (cameraFilePath);

	int numLines = 0;

	if (fs.is_open())
	{
		std::string str;

		while (!fs.eof())
		{
			std::string line1, line2 , line3 , line4;

			if (std::getline(fs, line1) && std::getline(fs, line2) && std::getline(fs, line3) && std::getline(fs, line4))
			{
				Eigen::Matrix< double, 3, 4 > projMat;

				std::stringstream sstr1(line1), sstr2(line2), sstr3(line3);

				sstr1 >> projMat(0, 0) >> projMat(0, 1) >> projMat(0, 2) >> projMat(0, 3);
				sstr2>> projMat(1, 0) >> projMat(1, 1) >> projMat(1, 2) >> projMat(1, 3);
				sstr3 >> projMat(2, 0) >> projMat(2, 1) >> projMat(2, 2) >> projMat(2, 3);

				projectionMatrices.push_back(projMat);
			}
			else
			{
				break;
			}
		}
	}

}
