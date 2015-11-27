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

#ifndef EPVH_H
#define EPVH_H

#include "basevh.h"
#include "segmentclipper.h"
#include "vtkincludes.h"


namespace tr{
  
  
struct Edge
{
  Point3d point1 , point2;
  
  double depth1 , depth2;
  
  int leftFirst , rightFirst , leftSecond , rightSecond; 
  
  int camLeft , contourLeft , stripLeft;
  int camRight , contourRight , stripRight;
  int camFirst , contourFirst , stripFirst;
  int camSecond , contourSecond , stripSecond;
  
  tr::Vertex *mVertex1 , *mVertex2;
  
};  


class EPVH : public BaseVH
{  
  
  
  
protected:  
  
  
  std::vector< std::vector< Edge > > mEdges ;

  std::vector< std::vector< std::vector< Edge > > > mCameraViewingEdges;
  
  std::vector< Vertex* > mVertices;
  
  std::vector< SegmentClipper > mClippers;

  int mVertexChainCounter;
  
  bool mDisplayDebug; 
  
  void computeEpipolarRayBins( int camId1 , int camId2  , std::vector< EdgeList > &strips  , std::vector< float > &slopeValues );
  void initiateViewingEdges( int camId1 , int camId2 , std::vector< std::vector< Edge > >  &viewingEdges);
  void initiateViewingEdges2( int camId1 , int camId2 , std::vector< std::vector< Edge > >  &viewingEdges );
  void updateViewingEdges(  int camId1 , int camId2 , std::vector< std::vector< Edge > >  &viewingEdges );
  
  void computeGeneratorNormals();
  void computeGeneratorNormals( int camId );
  void associateViewingEdgesToGenerators();
  void associateViewingEdgesToGenerators( int camId );
  
  std::vector< tr::Edge > computeViewingEdgeIntervalIntersection( std::vector< Edge > &edgeSet1 , std::vector< Edge > &edgeSet2 );
  
  void buildClippers();

  void buildClipper( int camId );
 
  void initiateWithViewingEdges();  
  
  void initializeVertices();

  void clearVertices();
  
  void completeVertex( Vertex *vertex );
  
  
  
  void generateEdgePolygons( vtkSmartPointer< vtkPolyData > &polygons );
  
  tr::Point3d estimateMaximalPoint( tr::Point3d &startPoint ,  tr::Vector3d &direction , 
			                        Generator* gen1 , Generator *gen2 ,  
									uchar &intersectionIndex , int ignoreIndex = -1 );
  
  bool clipEdge( int camId1 , int camId2 , tr::Point3d &point1 , tr::Point3d &point2 , 
		         int &clipCamId , int &clipContourId , int &clipStripId  );
  
  void clipEdge( tr::Edge edge , int camId , std::vector< tr::Edge > &clippedEdges );
  
  
  
  void generatePolyData( std::vector< std::vector< cv::Point2f > >& contours , 
			     std::vector< std::pair< tr::Point2f , tr::Point2f > > &edges, 
			     vtkSmartPointer< vtkPolyData >& edgePolyData);
  
  void generatePolyData(
                         std::vector< std::vector< cv::Point2f > >& contours,
		                 std::vector< std::pair< tr::Point2f, tr::Point2f > >& edges, 
		                 vtkSmartPointer< vtkPolyData >& edgePolyData, 
		                 std::vector< std::pair< int, int > >& selectedPolygonSegments ,
		                 std::vector< cv::Vec3f >& selectedSegmentColors
		               );
  
  void display( std::vector< std::vector< cv::Point2f > >& contours, std::vector< std::pair< tr::Point2f, tr::Point2f > >& edges );
  
  void displayWithColor( int camId ,std::vector< std::vector< cv::Point2f > >& contours, std::vector< std::pair< tr::Point2f, tr::Point2f > >& edges   );
  
  void buildVisualHull( std::vector< int > &silhouetteIds , vtkSmartPointer< vtkPolyData > &visualHull );
  
  double clipWithGenerator( tr::Point3d &p1 , tr::Point3d &p2 , int camId , int contourId , int stripId );
 
  double clipRayWithGenerator( tr::Point3d &o , tr::Vector3d &r , int camId , int contourId , int stripId  );
  
public:
    
  EPVH();
  
  virtual bool compute();

  bool compute( std::vector< int > &silhouetteCameras , std::vector< int > &imageBoundaryCameras );

  bool addSilhouetteCameraAndUpdate( int camId );

  bool generatePolygons( vtkSmartPointer< vtkPolyData > &polygons );

  void copyTo( tr::EPVH &dst );

  void clear();
    
  ~EPVH();

};

}

#endif // EPVH_H
