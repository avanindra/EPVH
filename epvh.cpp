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

#include "epvh.h"
#include "sortingpredicates.h"
#include <math2d.h>
#include "display3droutines.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"
#include "vtkUnsignedCharArray.h"
#include "vtkCellData.h"
#include "vtkSmartPointer.h"
#include "limits.h"
#include "math3d.h"
// #include "opencv/cv".h
namespace tr{

	EPVH::EPVH()
	{

	}

	void EPVH::computeEpipolarRayBins( int camId1 , int camId2 , std::vector< EdgeList >& strips, std::vector< double >& slopeValues)
	{
		std::vector< PriorityEdgeQueue > tempStrips;

		strips.clear();

        Eigen::Vector2d epipole;

		Point3d cameraCenter1;

		mCalibration->get_camera_center( camId1 , cameraCenter1 );

        mCalibration->projectPoint( cameraCenter1 , camId2 , epipole );
		//add check for pole on infinity.

		int numContours = mObjectContours[ camId2 ].size();

		int numTotalPoints = 0 , numSlopePoints = 0;

		for( int contour = 0; contour < numContours ; contour++ )
		{
			numTotalPoints += mObjectContours[ camId2 ][ contour ].size();
			numSlopePoints += mObjectContours[ camId2 ][ contour ].size();
		}

		strips.resize( numSlopePoints );
		tempStrips.resize( numSlopePoints );

		std::vector< std::pair< int , double > > slopes( numSlopePoints ) , backupSlopes;
		std::vector< int > reverseSlopeIndex( numTotalPoints , -1 );

		std::vector< double > edgeDistances( numTotalPoints , 0 );

		int totalStripId = 0 , slopeId = 0;

		double w = mCalibration->get_image_width( camId2 );
		double h = mCalibration->get_image_height( camId2 );

        Eigen::Vector2d refPoint( w / 2 , h/ 2 );

		Vector2 vec1 = refPoint - epipole;

		std::vector< int > offSet( numContours , 0 );

		slopeValues.clear();

		slopeValues.resize( numSlopePoints );

		for( int contour = 0; contour < numContours ; contour++ )
		{
			int numPoints = mObjectContours[ camId2 ][ contour ].size();

			if( contour > 0 )
			{
				offSet[ contour ] = offSet[ contour - 1 ] + mObjectContours[ camId2 ][ contour - 1 ].size();
			}


			for( int pp = 0; pp < numPoints; pp++ , totalStripId++ , slopeId++ )
			{
                Eigen::Vector2d silhouettePoint;

				silhouettePoint.x() = ( mObjectContours[ camId2 ][ contour ][ pp ].x + mOffset[ camId2 ].x ) * mInvScale[ camId2 ];
				silhouettePoint.y() = ( mObjectContours[ camId2 ][ contour ][ pp ].y + mOffset[ camId2 ].y ) * mInvScale[ camId2 ];

				Vector2 vec2 = silhouettePoint - epipole;

				double val1 = vec1.dot( vec2 );
				double val2 = ( vec1.x() * vec2.y() - vec2.x() * vec1.y() );

				slopes[ slopeId ].second = val2 / val1;

				slopes[ slopeId ].first = totalStripId;

				edgeDistances[ totalStripId ] = vec2.norm();

			}
		}

		backupSlopes = slopes;

		std::sort( slopes.begin() , slopes.end() , pointSlopePredicate );

		int numSlopes = slopes.size();

		for( int ss = 0; ss < numSlopes; ss++ )
		{
			int totalStripId1 = slopes[ ss ].first;

			reverseSlopeIndex[ totalStripId1 ] = ss;

			slopeValues[ ss ] = slopes[ ss ].second;

		}

		totalStripId = 0;

		for( int ss = 0; ss < numSlopes; ss++ )
		{
			int stripId = mStripContourMap[ camId2 ][ slopes[ ss ].first ].contourStripId;
			int contourId = mStripContourMap[ camId2 ][ slopes[ ss ].first ].contourId;

			int numStrips = mObjectContours[ camId2 ][ contourId ].size();

			int currentId = slopes[ ss ].first;

			int nextId = offSet[ contourId ] + ( stripId + 1 ) % numStrips;

			int minId , maxId;

			if( ss > reverseSlopeIndex[ nextId ] )
			{
				maxId = ss;
				minId = reverseSlopeIndex[ nextId ];
			}
			else
			{
				maxId = reverseSlopeIndex[ nextId ];
				minId = ss;
			}

			for( int id = minId; id < maxId; id++ )
			{
				EdgeDistancePair edgeDistPair;

				edgeDistPair.first = currentId;

                Eigen::Vector2d end1 , end2;

				int cid1 ,stId1 , cid2 , stId2;

				int tspid1 = slopes[ id ].first;
				int tspid2 = slopes[ id + 1 ].first;

				stId1 = mStripContourMap[ camId2 ][ tspid1 ].contourStripId;
				cid1 = mStripContourMap[ camId2 ][ tspid1 ].contourId;

				stId2 = mStripContourMap[ camId2 ][ tspid2 ].contourStripId;
				cid2 = mStripContourMap[ camId2 ][ tspid2 ].contourId;

				end1.x() = ( mObjectContours[ camId2 ][ cid1 ][ stId1 ].x + mOffset[ camId2 ].x ) * mInvScale[ camId2 ];
				end1.y() = ( mObjectContours[ camId2 ][ cid1 ][ stId1 ].y + mOffset[ camId2 ].y ) * mInvScale[ camId2 ];

				end2.x() = ( mObjectContours[ camId2 ][ cid2 ][ stId2 ].x + mOffset[ camId2 ].x ) * mInvScale[ camId2 ];
				end2.y() = ( mObjectContours[ camId2 ][ cid2 ][ stId2 ].y + mOffset[ camId2 ].y ) * mInvScale[ camId2 ];

                Eigen::Vector2d secondEnd = ( end1 + end2 ) / 2;

				edgeDistPair.second = distanceToTheStrip( epipole , secondEnd , camId2 , contourId , stripId );

				tempStrips[ id ].push( edgeDistPair );
			}

		}

		for( int ss = 0; ss < numSlopes; ss++ )
		{
			int numEdges = tempStrips[ ss ].size();

			for( int ee = 0; ee < numEdges; ee++ )
			{
				strips[ ss ].push_back( tempStrips[ ss ].top() );

				tempStrips[ ss ].pop();
			}
		}

	}


	void EPVH::buildClippers( )
	{
		int numSilhouetteCams = mSilhouetteCameras.size();

		mClippers.resize( mCalibration->get_cam_count() );

		for( int sc = 0; sc < numSilhouetteCams; sc++ )
		{
			int camId = mSilhouetteCameras[ sc ];

			mClippers[ camId ].addContours( mObjectContours[ camId ] );

			mClippers[ camId ].setContourHierarchy( mContourHierarchies[ camId ] );

			mClippers[ camId ].build();
		}
	}


	void EPVH::buildClipper( int camId )
	{
		mClippers[ camId ].addContours( mObjectContours[ camId ] );

		mClippers[ camId ].setContourHierarchy( mContourHierarchies[ camId ] );

		mClippers[ camId ].build();
	}



	void EPVH::initiateWithViewingEdges()
	{

		int numCams = mSilhouetteCameras.size();

		std::vector< std::vector< tr::Edge > > viewingEdges;

		

		for( int ss = 0; ss < numCams; ss++ )
		{
			std::cout << " camera : " << ss << std::endl;

			int sCamId = mSilhouetteCameras[ ss ];

			//displayContour3D(mObjectContours[ss]);

			std::vector< std::vector< tr::Edge > > edges1 ;

			int orthogonalCam = mMostOrthogonalCamera[ sCamId ];

			initiateViewingEdges( sCamId , orthogonalCam , edges1 );

			for( int ss2 = 0; ss2 < numCams; ss2++ )
			{ 

				if( ( ss2 == ss ) || ( mSilhouetteCameras[ ss2 ] == mMostOrthogonalCamera[ sCamId ]  ) )
					continue;

				int sCamId2 = mSilhouetteCameras[ ss2 ];

				updateViewingEdges( sCamId , sCamId2 , edges1 );

			}

			mCameraViewingEdges[ sCamId ] = edges1;

			int numEdgeSets = edges1.size();

			int count = 0;
			for( int ee = 0;  ee < numEdgeSets; ee++ )
			{
				if( edges1[ ee ].size() > 0 )
				{
					mEdges.push_back( edges1[ ee ] );

					int numEdges = edges1[ ee ].size();

					count++;
				}
			}



//            vtkSmartPointer< vtkPolyData > edgeData = vtkSmartPointer< vtkPolyData >::New();

//            generatePolyData(mEdges, edgeData);

//            tr::Display3DRoutines::displayPolyData(edgeData);

		}



	}



	void EPVH::initializeVertices()
	{


		int numCams = mCameraViewingEdges.size();

		for( int cc = 0; cc < numCams; cc++ )
		{

			int numEdgeSets = mCameraViewingEdges[ cc ].size();

			for( int es = 0; es < numEdgeSets; es++ )
			{

				int numEdges = mCameraViewingEdges[ cc ][ es ].size();

				for( int ee = 0; ee < numEdges; ee++ )
				{    
					tr::Edge &edge = mCameraViewingEdges[ cc ][ es ][ ee ];

					Vertex *v1 = new Vertex();

					v1->mCoords = edge.point1;

					v1->mLeftGen = mGenerators[ edge.camRight ][ edge.contourRight ][ edge.stripRight ];
					v1->mRightGen = mGenerators[ edge.camLeft ][ edge.contourLeft ][ edge.stripLeft ];
					v1->mGen = mGenerators[ edge.camFirst ][ edge.contourFirst ][ edge.stripFirst ];

					v1->mIgnoreIndex = 0;

					v1->mLeft = 0;
					v1->mRight = 0; 


					Vertex *v2 = new Vertex();

					v2->mCoords = edge.point2;
					v2->mLeftGen = mGenerators[ edge.camLeft ][ edge.contourLeft ][ edge.stripLeft ];
					v2->mRightGen = mGenerators[ edge.camRight ][ edge.contourRight ][ edge.stripRight ];

					// std::cout<<cc<<" "<<edge.camSecond<<" "<<edge.contourSecond<<" "<<edge.stripSecond<<std::endl; 
					v2->mGen = mGenerators[ edge.camSecond ][ edge.contourSecond ][ edge.stripSecond ];

					v2->mLeft = 0;
					v2->mRight = 0;

					v1->mId = mVertices.size();
					v2->mId = mVertices.size() + 1;

					mVertices.push_back( v1 );
					mVertices.push_back( v2 ); 

					v1->mLeftGen->mVertices.push_back( v1 );
					v2->mRightGen->mVertices.push_back( v2 );
					v1->mLeftGen->mAllVertices.push_back( v1 );
					v2->mRightGen->mAllVertices.push_back( v2 );
					v1->mGen->mAllVertices.push_back( v1 );
					v2->mGen->mAllVertices.push_back( v2 );

                    v1->mLE = ( edge.point1 - edge.point2 );
                    v2->mLE = ( edge.point2 - edge.point1 );

					v1->mLE.normalize();
					v2->mLE.normalize();  

					v2->mIgnoreIndex = 1;

					edge.mVertex1 = v1;
					edge.mVertex2 = v2;

					v1->mOpp = v2;
					v2->mOpp = v1;

					v1->mIsGenVertex = 1;
					v2->mIsGenVertex = 1;

				}

			}

		}

	}


	void EPVH::clearVertices()
	{
		int size1 = mGenerators.size();

		for( int ss = 0; ss < size1 ; ss++ )
		{
			int size2 = mGenerators[ ss ].size();

			for( int ss2 = 0; ss2 < size2; ss2++ )
			{
				int size3 = mGenerators[ ss ][ ss2 ].size();

				for( int ss3 = 0; ss3 < size3 ; ss3++ )
				{
					if( mGenerators[ ss ][ ss2 ][ ss3 ] )
					{
						mGenerators[ ss ][ ss2 ][ ss3 ]->mVertices.clear();
						mGenerators[ ss ][ ss2 ][ ss3 ]->mAllVertices.clear();
						mGenerators[ ss ][ ss2 ][ ss3 ]->mTriplePointVertices.clear();

					}
				}
			}
		}

		int numVertices = mVertices.size();

		for( int vv = 0; vv < numVertices; vv++ )
		{
			if( mVertices[ vv ] )
				delete mVertices[ vv ];

			mVertices[ vv ] = 0;
		}

		mVertices.clear();
	}


	void EPVH::completeVertex( Vertex* vertex )
	{
		bool leftTriplePoint = false;
		bool rightTriplePoint = false;

		mVertexChainCounter++;

		if( mVertexChainCounter > 5000 )
			return;


		if( !vertex->mLeft )
		{
			tr::Vector3d nLeft , nGen , lE , lLeft; 

			nLeft = vertex->mLeftGen->normal;
			nGen = vertex->mGen->normal;

			lLeft = nLeft.cross( nGen );

			int k = 1;

			if( vertex->mLE.cross( lLeft ).dot( nLeft ) < 0 )
			{
				k = -1;
			}

			lLeft *= -1;    

			uchar intersectionIndex = 4;

			int ignoreIndex = -1;

			if( vertex->mIsGenVertex )
			{
				if( vertex->mIgnoreIndex == 1 )
				{
					ignoreIndex = 1;
				}
				else
				{
					ignoreIndex = 0;
				}
			}

            Eigen::Vector3d maximalPoint = estimateMaximalPoint( vertex->mCoords , lLeft , vertex->mLeftGen ,
				vertex->mGen  , intersectionIndex , ignoreIndex );


			Generator *maximalPointGen = 0 , *gen2 = 0;

			bool leftViewingEdge = true;


			if( intersectionIndex == 0 || intersectionIndex == 1  )
			{
				maximalPointGen = vertex->mLeftGen;

				if( intersectionIndex == 1 )
				{ 
					gen2 = maximalPointGen->mRightGen;

					leftViewingEdge = false;
				}
				else
				{
					gen2 = maximalPointGen;
				}
			}
			else if( intersectionIndex == 2 || intersectionIndex == 3 )
			{
				maximalPointGen = vertex->mGen;

				if( intersectionIndex == 3 )
				{
					gen2 = maximalPointGen->mRightGen;

					leftViewingEdge = false;

				}
				else
				{
					gen2 = maximalPointGen; 
				}
			}


			int clipCamId = -1;
			int clipContourId = -1;
			int clipStripId = -1;


			if( clipEdge( vertex->mLeftGen->camId , vertex->mGen->camId , 
				vertex->mCoords , maximalPoint , clipCamId , 
				clipContourId , clipStripId )  )
			{


				Generator *gen = mGenerators[ clipCamId ][ clipContourId ][ clipStripId ];

				int numTriplePointVertices = gen->mTriplePointVertices.size();

				Vertex *newVertex = 0;

				for( int vv = 0; vv < numTriplePointVertices; vv++ )
				{
					Vertex *tpv = gen->mTriplePointVertices[ vv ];

					Generator *ltg = tpv->mLeftGen;
					Generator *rtg = tpv->mRightGen;
					Generator *tg = tpv->mGen;

					Generator *lg = vertex->mLeftGen;
					Generator *g = vertex->mGen;

					if( ( ltg == lg || ltg == gen || ltg == g ) &&
						( rtg == lg || rtg == gen || rtg == g ) &&
						( tg == lg || tg == gen || tg == g )
						)
					{
						vertex->mLeft = tpv;

						if( ( ltg == lg  || ltg == g ) )
						{
							tpv->mLeft = vertex;
						}
						else
						{
							tpv->mRight = vertex;
						}

						newVertex = tpv;

						newVertex->mGen->mAllVertices.push_back( newVertex );
						newVertex->mLeftGen->mAllVertices.push_back( newVertex );
						newVertex->mRightGen->mAllVertices.push_back( newVertex );
					}

				}


				if( !newVertex )
				{
					Vertex *newVertex = new Vertex();

					newVertex->mLeft = 0;
					newVertex->mRight = 0;

					newVertex->mIsGenVertex = 0;

					newVertex->mOpp = vertex;

					newVertex->mCoords = maximalPoint;

					newVertex->mGen = mGenerators[ clipCamId ][ clipContourId ][ clipStripId ];

					newVertex->mLeftGen = vertex->mLeftGen;
					newVertex->mRightGen = vertex->mGen;	

					newVertex->mId = mVertices.size();

					vertex->mLeft = newVertex;

					newVertex->mGen->mAllVertices.push_back( newVertex );
					newVertex->mLeftGen->mAllVertices.push_back( newVertex );
					newVertex->mRightGen->mAllVertices.push_back( newVertex );


					newVertex->mGen->mTriplePointVertices.push_back( newVertex );
					newVertex->mLeftGen->mTriplePointVertices.push_back( newVertex );
					newVertex->mRightGen->mTriplePointVertices.push_back( newVertex );

					newVertex->mIgnoreIndex = -1;

					mVertices.push_back( newVertex );

					leftTriplePoint = true;
				}


			}
			else
			{

				vertex->mLeftCoords = maximalPoint;
				//connect to a generator vertex

				int numGenVertices; 

				if( leftViewingEdge )	
				{
					numGenVertices = maximalPointGen->mVertices.size();
				}	  
				else
				{
					numGenVertices = maximalPointGen->mRightGen->mVertices.size();
				}

				int numConnections = 0;

				for( int vv = 0; vv < numGenVertices; vv++ )
				{
					Vertex *seconVertex;


					if( leftViewingEdge )
						seconVertex = maximalPointGen->mVertices[ vv ];
					else
						seconVertex = maximalPointGen->mRightGen->mVertices[ vv ];

					if( 
						( ( vertex->mLeftGen == seconVertex->mGen ) || ( vertex->mLeftGen == seconVertex->mLeftGen ) || ( vertex->mLeftGen == seconVertex->mRightGen ) ) &&
						( ( vertex->mGen == seconVertex->mGen ) || ( vertex->mGen == seconVertex->mLeftGen ) || ( vertex->mGen == seconVertex->mRightGen ) )
						)
					{

						vertex->mLeft = seconVertex;   

						if( ( seconVertex->mLeftGen == vertex->mLeftGen ) || 
							( seconVertex->mLeftGen == vertex->mGen )    )
						{ 
							if( !seconVertex->mLeft )
								seconVertex->mLeft = vertex;
						}
						else if( ( seconVertex->mRightGen == vertex->mLeftGen ) || 
							( seconVertex->mRightGen == vertex->mGen )    )
						{
							if( !seconVertex->mRight )
								seconVertex->mRight = vertex;	     

						}


						numConnections++;

					}


				}

			}

		}

		if( leftTriplePoint )
		{
			completeVertex( vertex->mLeft );  
		}

		if( !vertex->mRight )
		{
			tr::Vector3d nRight , nGen , lE , lRight; 

			nRight = vertex->mRightGen->normal;
			nGen = vertex->mGen->normal;

			lRight = nGen.cross( nRight );

			int k = 1;

			if( lRight.cross( vertex->mLE ).dot( nGen ) < 0 )
			{
				k = -1;
			}

			lRight *= -1;    

			uchar intersectionIndex = 4;

			int ignoreIndex = -1;

			if( vertex->mIsGenVertex )
			{
				if( vertex->mIgnoreIndex == 1 )
				{
					ignoreIndex = 0;
				}
				else
				{
					ignoreIndex = 1;
				}
			}

            Eigen::Vector3d maximalPoint = estimateMaximalPoint( vertex->mCoords , lRight , vertex->mRightGen ,
				vertex->mGen  , intersectionIndex , ignoreIndex );


			Generator *maximalPointGen = 0, *gen2 = 0;

			bool leftViewingEdge = true;

			if( intersectionIndex == 0 || intersectionIndex == 1 )
			{
				maximalPointGen = vertex->mRightGen;

				if( intersectionIndex == 1 )
				{
					leftViewingEdge = false;

					gen2 = maximalPointGen->mRightGen;	  
				}
				else
				{
					gen2 = maximalPointGen;
				}
			}
			else if( intersectionIndex == 2 || intersectionIndex == 3 )
			{
				maximalPointGen = vertex->mGen;

				if( intersectionIndex == 3 )
				{
					leftViewingEdge = false; 
					gen2 = maximalPointGen->mRightGen;
				}
				else
				{
					gen2 = maximalPointGen;
				}
			}


			int clipCamId = -1;
			int clipContourId = -1;
			int clipStripId = -1;

			int numGenVertices; 



			if( clipEdge( vertex->mRightGen->camId , vertex->mGen->camId , 
				vertex->mCoords , maximalPoint , clipCamId , 
				clipContourId , clipStripId )  )
			{

				Generator *gen = mGenerators[ clipCamId ][ clipContourId ][ clipStripId ];

				int numTriplePointVertices = gen->mTriplePointVertices.size();

				Vertex *newVertex = 0;

				for( int vv = 0; vv < numTriplePointVertices; vv++ )
				{
					Vertex *tpv = gen->mTriplePointVertices[ vv ];

					Generator *ltg = tpv->mLeftGen;
					Generator *rtg = tpv->mRightGen;
					Generator *tg = tpv->mGen;

					Generator *rg = vertex->mRightGen;
					Generator *g = vertex->mGen;

					if( ( ltg == gen || ltg == rg || ltg == g ) &&
						( rtg == gen || rtg == rg || rtg == g ) &&
						( tg == gen || tg == rg || tg == g )
						)
					{
						vertex->mRight = tpv;

						if( (  ltg == rg || ltg == g ) )
						{
							tpv->mLeft = vertex;
						}
						else
						{
							tpv->mRight = vertex;
						}

						newVertex = tpv;

						newVertex->mGen->mAllVertices.push_back( newVertex );
						newVertex->mLeftGen->mAllVertices.push_back( newVertex );
						newVertex->mRightGen->mAllVertices.push_back( newVertex );
					}

				}

				//std::cout<<" triple point "<<clipContourId<<" "<<clipStripId<<std::endl;

				if( !newVertex )
				{
					newVertex = new Vertex();

					newVertex->mIgnoreIndex = -1;

					newVertex->mLeft = 0;
					newVertex->mRight = 0;

					newVertex->mIsGenVertex = 0;

					newVertex->mOpp = vertex;

					newVertex->mCoords = maximalPoint;

					newVertex->mGen = mGenerators[ clipCamId ][ clipContourId ][ clipStripId ];

					newVertex->mLeftGen = vertex->mGen;
					newVertex->mRightGen = vertex->mRightGen;	

					newVertex->mId = mVertices.size();

					vertex->mRight = newVertex;

					newVertex->mGen->mAllVertices.push_back( newVertex );
					newVertex->mLeftGen->mAllVertices.push_back( newVertex );
					newVertex->mRightGen->mAllVertices.push_back( newVertex );

					mVertices.push_back( newVertex );

					rightTriplePoint = true;

					newVertex->mGen->mTriplePointVertices.push_back( newVertex );
					newVertex->mLeftGen->mTriplePointVertices.push_back( newVertex );
					newVertex->mRightGen->mTriplePointVertices.push_back( newVertex );

				}
			}
			else
			{

				vertex->mRightCoords = maximalPoint;
				//connect to a generator vertex

				int numGenVertices;

				if( leftViewingEdge )	  
					numGenVertices = maximalPointGen->mVertices.size();	  
				else	    
					numGenVertices = maximalPointGen->mRightGen->mVertices.size();

				int numConnections = 0;

				for( int vv = 0; vv < numGenVertices; vv++ )
				{
					Vertex *seconVertex;

					if( leftViewingEdge )
						seconVertex = maximalPointGen->mVertices[ vv ];
					else
						seconVertex = maximalPointGen->mRightGen->mVertices[ vv ];

					if( 
						( ( vertex->mRightGen == seconVertex->mGen ) || ( vertex->mRightGen == seconVertex->mLeftGen ) || ( vertex->mRightGen == seconVertex->mRightGen ) ) &&
						( ( vertex->mGen == seconVertex->mGen ) || ( vertex->mGen == seconVertex->mLeftGen ) || ( vertex->mGen == seconVertex->mRightGen ) )
						)
					{

						vertex->mRight = seconVertex;	

						if( ( seconVertex->mLeftGen == vertex->mRightGen ) || 
							( seconVertex->mLeftGen == vertex->mGen )    )
						{ 
							if( !seconVertex->mLeft )
							{
								seconVertex->mLeft = vertex;		
							}

						}
						else if( ( seconVertex->mRightGen == vertex->mLeftGen ) || 
							( seconVertex->mRightGen == vertex->mGen )   )
						{
							if( !seconVertex->mRight )
							{
								seconVertex->mRight = vertex;		
							}


						}

						numConnections++;

					}


				}

			}
		}


		if( rightTriplePoint )
		{
			completeVertex( vertex->mRight );  
		}
	}


	bool EPVH::generatePolygons( vtkSmartPointer< vtkPolyData >& polygons )
	{

		vtkSmartPointer< vtkPolyData > polyData = vtkSmartPointer< vtkPolyData >::New();

		polyData->Reset();

		int numVertices = mVertices.size();

		polyData->Allocate( numVertices * 10 );

		vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();

		points->Allocate( numVertices );


		std::vector< std::vector< uchar > > isVertexUsedInPolygons( numVertices ) ;

		std::vector< bool > vertexUsedFlag( numVertices , false );

		for( int vv = 0; vv < numVertices; vv++ )
		{
			std::vector< uchar > flags( 3 , 0 );//left , right , generator

			isVertexUsedInPolygons[ vv ] = flags;

			points->InsertNextPoint( mVertices[ vv ]->mCoords.data() );
		}


		vtkSmartPointer< vtkIdList > polygon = vtkSmartPointer< vtkIdList >::New();


		int numSilhouetteCams = mSilhouetteCameras.size();

		for( int sc = 0; sc < numSilhouetteCams; sc++ )
		{
			int sCamId = mSilhouetteCameras[ sc ];

			int numContours = mObjectContours[ sCamId ].size();

			for( int cc = 0; cc < numContours; cc++ )
			{
				int numStrips = mObjectContours[ sCamId ][ cc ].size();

				for( int ss = 0; ss < numStrips; ss++ )
				{

					Generator *gen = mGenerators[ sCamId ][ cc ][ ss ];	

					int numVertices =  gen->mAllVertices.size();	

					std::vector< int > resetIds;	

					for( int vv = 0; vv < numVertices; vv++ )	
					{
						Vertex *vert = gen->mAllVertices[ vv ];

						tr::Vertex *prevVertex = 0; 

						//if( !vertexUsedFlag[ vert->mId ] )
						{
							Vertex *vNext = vert->mRight;

							tr::Vertex *currentVertex = vert;      
							tr::Vertex *prevVertex = 0;

							polygon->Reset();

                            if( !vNext )
                            {
                                polygon->Reset();

								std::cout<<" connection error "<<std::endl;

                                break;

                             
                            }

							if( vNext->mLeftGen == gen || vNext->mRightGen == gen || vNext->mGen == gen )
							{ 
								vertexUsedFlag[ vert->mId ] = true;

								prevVertex = vert;

								currentVertex = vert->mRight;

								vertexUsedFlag[ currentVertex->mId ] = true;

								polygon->InsertNextId( vert->mId );

								resetIds.push_back( vert->mId );
								resetIds.push_back( currentVertex->mId );

								mVertexChainCounter = 0;

								while( 1 )
								{
									mVertexChainCounter++;

									bool genCase = true;

									if( mVertexChainCounter > 5000 )
                                    {
                                        polygon->Reset();

										std::cout << " connection error " << std::endl;

                                        break;

                                     }
									polygon->InsertNextId( currentVertex->mId );

									if( currentVertex->mRight != prevVertex  )
									{ 
										tr::Vertex *v = currentVertex->mRight;
                                        if( !v )
                                        {
                                            polygon->Reset();

											std::cout<<" connection error "<<std::endl;

                                            break;

                                         
                                        }


										vertexUsedFlag[ currentVertex->mId ] = true;	  


										if(  v->mRightGen == gen || v->mLeftGen == gen || v->mGen == gen)
										{	    
											genCase = false;

											prevVertex = currentVertex;

											currentVertex = v;

											resetIds.push_back( currentVertex->mId );
										}

									}

									if( currentVertex->mLeft != prevVertex && genCase  )
									{ 



										tr::Vertex *v = currentVertex->mLeft;	  

                                        if( !v )
                                        {
                                            polygon->Reset();

											std::cout<<" connection error "<<std::endl;

                                            break;

                                        
                                        }

										if(  v->mRightGen == gen || v->mLeftGen == gen || v->mGen == gen)
										{	    
											genCase = false;

											prevVertex = currentVertex;

											currentVertex = v;

											vertexUsedFlag[ currentVertex->mId ] = true;

											resetIds.push_back( currentVertex->mId );
										}

									}

									if( currentVertex->mOpp != prevVertex &&  genCase )
									{		  
										tr::Vertex *v = currentVertex->mOpp;	  

										if(  v->mRightGen == gen || v->mLeftGen == gen || v->mGen == gen)
										{	    
											genCase = false;

											prevVertex = currentVertex;

											currentVertex = v;

											vertexUsedFlag[ currentVertex->mId ] = true;

											resetIds.push_back( currentVertex->mId );
										}
									}

									if( currentVertex == vert || genCase )
									{	  
										break;
									}

								}


                                if( polygon->GetNumberOfIds() >= 3 )
								polyData->InsertNextCell( VTK_POLYGON , polygon );
							}

						}
					}

					for( int vv = 0; vv < resetIds.size(); vv++ )
					{
						vertexUsedFlag[ resetIds[ vv ] ] = false;
					}


				}
			}
		}

		polyData->SetPoints( points );

		polygons->DeepCopy( polyData ); 

		return true;

	}

	void EPVH::generateEdgePolygons( vtkSmartPointer< vtkPolyData >& polygons )
	{
		int numVertices = mVertices.size();

		polygons->Allocate( numVertices * 10 );

		vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();

		points->Allocate( numVertices );

		vtkSmartPointer< vtkIdList > edge = vtkSmartPointer< vtkIdList >::New();

		for( int vv = 0; vv < numVertices; vv++ )
		{
			tr::Vertex *vertex = mVertices[ vv ];

			points->InsertNextPoint( vertex->mCoords.data() );
		}

		for( int vv = 0; vv < numVertices; vv++ )
		{
			tr::Vertex *vertex = mVertices[ vv ];
			tr::Generator *gen = 0; 

			if( vertex->mLeft )
			{
				edge->InsertNextId( vertex->mId  );
				edge->InsertNextId( vertex->mLeft->mId  );

				polygons->InsertNextCell( VTK_LINE , edge );
			}

			edge->Reset();

			if( vertex->mRight )
			{
				edge->InsertNextId( vertex->mId  );
				edge->InsertNextId( vertex->mRight->mId  );

				polygons->InsertNextCell( VTK_LINE , edge );
			}

			edge->Reset();


			edge->InsertNextId( vertex->mId  );
			edge->InsertNextId( vertex->mOpp->mId  );

			polygons->InsertNextCell( VTK_LINE , edge );

			edge->Reset();

		}


		polygons->SetPoints( points );

	}


	void EPVH::copyTo( tr::EPVH &dst )
	{
		/* Not complete yet. Do not use it */

		dst.clear();

		dst.mClippers = mClippers;

		dst.mCalibration = mCalibration;

		dst.mObjectContours = mObjectContours;

		dst.mScale = mScale;

		dst.mInvScale = mInvScale;

		dst.mContourHierarchies = mContourHierarchies;

		dst.mIsCameraUsed = mIsCameraUsed;

		dst.mIsOccludingContour = mIsOccludingContour;

		dst.mMostOrthogonalCamera = mMostOrthogonalCamera;

		dst.mOffset = mOffset;

		dst.mSilhouetteCameras = mSilhouetteCameras;

		dst.mCameraViewingEdges = mCameraViewingEdges;

		dst.mStripContourMap = mStripContourMap;

		int numImages = mCalibration->get_cam_count();

		dst.mGenerators.resize( numImages );

		for( int ii = 0; ii < numImages; ii++ )
		{
			int numContours = mGenerators[ ii ].size();

			dst.mGenerators[ ii ].resize( numContours );

			for( int cc = 0; cc < numContours ; cc++ )
			{
				int numGens = mGenerators[ ii ][ cc ].size();

				dst.mGenerators[ ii ][ cc ].resize( numGens );

				for( int gg = 0; gg < numGens ; gg++ )
				{
					dst.mGenerators[ ii ][ cc ][ gg ] = new Generator();

					*( dst.mGenerators[ ii ][ cc ][ gg ] ) = *( mGenerators[ ii ][ cc ][ gg ] );
				}
			}
		}

		dst.clearVertices();	
	}



	Point3d EPVH::estimateMaximalPoint(  Point3d& startPoint, Vector3d& direction, Generator* gen1,
		Generator* gen2, uchar& intersectionIndex , int ignoreIndex)
	{
		Point3d maximalPoint;

		intersectionIndex = 4;

		Point3d camCenter1 , camCenter2;

		mCalibration->get_camera_center( gen1->camId , camCenter1 );
		mCalibration->get_camera_center( gen2->camId , camCenter2 );

		double d1 = tr::Math3D::lineIntersection( startPoint , direction , camCenter1 , gen1->leftRay );
		double d2 = tr::Math3D::lineIntersection( startPoint , direction , camCenter1 , gen1->rightRay );
		double d3 = tr::Math3D::lineIntersection( startPoint , direction , camCenter2 , gen2->leftRay );
		double d4 = tr::Math3D::lineIntersection( startPoint , direction , camCenter2 , gen2->rightRay );

		double d = std::numeric_limits< double >::max(); 

		if( d1 > 0 &&  d1 < d && ignoreIndex != 0    )
		{ 
			d = d1;

			intersectionIndex = 0;
		}


		if( d2 > 0 &&  d2 < d && ignoreIndex != 1   )
		{ 
			d = d2;

			intersectionIndex = 1;
		}

		if( d3 > 0 && d3 < d  )
		{
			intersectionIndex = 2;

			d = d3;
		}

		if( d4 > 0 && d4 < d )
		{
			d = d4 ;

			intersectionIndex = 3;
		}

		maximalPoint = startPoint + d * direction;

		return maximalPoint;

	}

	bool EPVH::clipEdge( int camId1, int camId2, Point3d& point1, Point3d& point2, int& clipCamId, 
		int& clipContourId, int& clipStripId )
	{
		clipCamId = -1;
		clipContourId = -1;
		clipStripId = -1;

		bool isClipped = false;  
		int numSilhouetteCams = mSilhouetteCameras.size();  

        //   tr::Vector3d dir = Eigen::Vector2d - point1;

		for( int sc = 0; sc < numSilhouetteCams; sc++  )
		{  
            Eigen::Vector2d  proj1 ,  proj2 ;

			int camId = mSilhouetteCameras[ sc ];

			if( camId == camId1 || camId == camId2 )
				continue;

            mCalibration->projectPoint( point1 , camId , proj1 );
            mCalibration->projectPoint( point2 , camId , proj2 );

            Eigen::Vector2d clippedEnd;

			std::pair< int , int > strip;


			mClippers[ camId ].firstClip( proj1 , proj2 , clippedEnd , strip  );


			if( strip.first != -1 )
			{

				//       tr::Vector3d ray;
				//       
                //       Eigen::Vector3d cameraCenter;
				//       
				//       mCalibration->get_camera_center( camId , cameraCenter );
				//        
				//       mCalibration->getRay( camId , clippedEnd  , ray );
				// 
				//       double depth = Math3D::lineIntersection( point1 , dir , cameraCenter , ray );

                double depth = clipWithGenerator( point1 , point2 , camId , strip.first , strip.second );

                point2 = point1 + depth * ( point2 - point1 );

				clipCamId = camId;
				clipContourId = strip.first;
				clipStripId = strip.second;

			} 


		}

		return ( clipCamId != -1 );
	}


	void EPVH::clipEdge( Edge edge, int camId, std::vector< Edge >& clippedEdges )
	{

		std::vector< std::pair< int , int >  > strips;
        std::vector< std::pair< Eigen::Vector2d , Eigen::Vector2d > > clippedEdgesTemp , initEdges;

        Eigen::Vector2d end1 , end2;

        mCalibration->projectPoint( edge.point1 , camId , end1 );
        mCalibration->projectPoint( edge.point2 , camId , end2 );

        std::pair< Eigen::Vector2d , Eigen::Vector2d > inititialEdge;

		inititialEdge.first = end1;
		inititialEdge.second = end2;

		initEdges.push_back( inititialEdge );

		//   if( mClippers[ camId ].mDisplayInfo )
		//   displayWithColor( camId , mObjectContours[ camId  , initEdges] );

		mClippers[ camId ].clipSegment2( end1 , end2 , clippedEdgesTemp , strips  );

		//   if( mClippers[ camId ].mDisplayInfo )
		// 	  std::cout<<" num clipped strips : "<<strips.size() << std::endl;

		int numClippedEdges = clippedEdgesTemp.size();

        Vector3d dir = ( edge.point2 - edge.point1 );

		Point3d cameraCenter;

		mCalibration->get_camera_center( camId , cameraCenter );

		for( int ce = 0; ce < numClippedEdges; ce++ )
		{ 
			tr::Edge currentEdge = edge;   

			//if( mClippers[ camId ].mDisplayInfo )
			//  std::cout<<" display strip intersection : "<<strips[ 2 * ce ].first<<" "<<strips[ 2 * ce + 1 ].first << std::endl;

			if( strips[ 2 * ce ].first != -1 )
			{
				currentEdge.camFirst = camId;
				currentEdge.contourFirst = strips[ 2 * ce ].first;
				currentEdge.stripFirst = strips[ 2 * ce ].second;

				/*      tr::Vector3d ray;

				mCalibration->getRay( camId , clippedEdgesTemp[ ce ].first  , ray );

				double depth = Math3D::lineIntersection( edge.point1 , dir , cameraCenter , ray );

				currentEdge.point1 = edge.point1 + depth * dir;   */ 

                double depth = clipWithGenerator( edge.point1 , edge.point2 , camId , strips[ 2 * ce ].first , strips[ 2 * ce ].second );

				currentEdge.point1 = edge.point1 + depth * dir;

			}

			if( strips[ 2 * ce + 1 ].first != -1 )
			{
				currentEdge.camSecond = camId;
				currentEdge.contourSecond = strips[ 2 * ce + 1 ].first;
				currentEdge.stripSecond = strips[ 2 * ce + 1 ].second;

				//       tr::Vector3d ray;
				//        
				//       mCalibration->getRay( camId , clippedEdgesTemp[ ce ].second , ray );
				// 
				//       double depth = Math3D::lineIntersection( edge.point1 , dir , cameraCenter , ray );
				// 
                //       currentEdge.point2 = edge.point1 + depth * dir;

                double depth = clipWithGenerator( edge.point1 , edge.point2 , camId , strips[ 2 * ce + 1 ].first , strips[ 2 * ce + 1 ].second );

                currentEdge.point2 = edge.point1 + depth * dir;

			}

			clippedEdges.push_back( currentEdge );
		}


	}



	void EPVH::initiateViewingEdges( int camId1, int camId2, std::vector< std::vector< Edge > >& viewingEdges )
	{ 
		std::vector< EdgeList > strips;
		std::vector< double > slopeValues;    

		computeEpipolarRayBins( camId1 , camId2 , strips , slopeValues );

		double w = mCalibration->get_image_width( camId2 );
		double h = mCalibration->get_image_height( camId2 );

		//std::cout << w << " " << h << std::endl;

        Eigen::Vector2d epipole , refPoint( w / 2 , h / 2 );

		Point3d cameraCenter1 , cameraCenter2 , epipoleHmg1 , epipoleHmg2 ;

		mCalibration->get_camera_center( camId1 , cameraCenter1 );
		mCalibration->get_camera_center( camId2 , cameraCenter2 );

        mCalibration->projectPoint( cameraCenter1 , camId2 , epipole );

		Vector2 vec1 = refPoint - epipole;

		int numTotalPoints = 0  , numSlopePoints = 0;

		int numContours = mObjectContours[ camId1 ].size();

		for( int contour = 0; contour < numContours ; contour++ )
		{
			numTotalPoints += mObjectContours[ camId1 ][ contour ].size();

			numSlopePoints += mObjectContours[ camId1 ][ contour ].size();
		}


		std::vector< std::pair< int , double > > slopes( numSlopePoints ) ;

        std::vector< Eigen::Vector2d > infinityEnds( numTotalPoints );

		std::vector< Vector3d > raySet( numTotalPoints ) ;

		std::vector< int > offSet( numContours , 0 ) ;

		int totalStripId = 0 , slopeId = 0;

		viewingEdges.resize( numTotalPoints );

		for( int contour = 0; contour < numContours ; contour++ )
		{
			int numPoints = mObjectContours[ camId1 ][ contour ].size();

			if( contour > 0 )
			{
				offSet[ contour ] = offSet[ contour - 1 ] + mObjectContours[ camId1 ][ contour - 1 ].size();
			}

			for( int pp = 0; pp < numPoints; pp++ , totalStripId++ , slopeId++ )
			{
                Eigen::Vector2d silhouettePoint , infiniteEnd;

				silhouettePoint.x() = ( mObjectContours[ camId1 ][ contour ][ pp ].x + mOffset[ camId1 ].x ) * mInvScale[ camId1 ];
				silhouettePoint.y() = ( mObjectContours[ camId1 ][ contour ][ pp ].y + mOffset[ camId1 ].y ) * mInvScale[ camId1 ];

				Vector3d ray , infiniteProjHmg;

				mCalibration->getRay( camId1 , silhouettePoint , ray );

				Vector4d infinitePoint;

				infinitePoint( 0 ) = ray( 0 );
				infinitePoint( 1 ) = ray( 1 );
				infinitePoint( 2 ) = ray( 2 );
				infinitePoint( 3 ) = 0;

				mCalibration->projectHomogeneous( infinitePoint , camId2 , infiniteProjHmg  );

				infiniteEnd.x() = infiniteProjHmg( 0 ) / infiniteProjHmg( 2 );
				infiniteEnd.y() = infiniteProjHmg( 1 ) / infiniteProjHmg( 2 );

				infinityEnds[ totalStripId ] = infiniteEnd;

				raySet[ totalStripId ] = ray;

				mGenerators[ camId1 ][ contour ][ pp ]->leftRay = ray;

				mGenerators[ camId1 ][ contour ][ pp ]->leftViewEdgeId = -1;
				mGenerators[ camId1 ][ contour ][ pp ]->rightViewEdgeId = -1;

				int idx = ( pp - 1 + numPoints ) % numPoints;

				mGenerators[ camId1 ][ contour ][ idx ]->rightRay = ray;

				mGenerators[ camId1 ][ contour ][ idx ]->camId = camId1;

				Vector2 vec2 = infiniteEnd - epipole;

				double val1 = vec1.dot( vec2 );
				double val2 = ( vec1.x() * vec2.y() - vec2.x() * vec1.y() );

				slopes[ slopeId ].second = val2 / val1;

				slopes[ slopeId ].first = totalStripId;

			}
		}


		std::sort( slopes.begin() , slopes.end() , pointSlopePredicate );


		int numSlopes = slopes.size();
		int numBins = slopeValues.size();

		int currentId = 0;

        Eigen::Vector2d intersectionPoint;

		std::pair< double , double > coefficientPair;

		for( int ss = 0; ss < numSlopes; ss++ )
		{
			std::vector< Edge > viewingEdgeSet;

			double slope = slopes[ ss ].second;

			int totalStripId = slopes[ ss ].first;

			for( int id = currentId ; id < numBins - 1; id++ )
			{

				if( slope >= slopeValues[ id ] && slope < slopeValues[ id + 1 ] )
				{
					//possible intersection
					bool edgeStarted = false;

					for( EdgeList::iterator iter = strips[ id ].begin() ; iter != strips[ id ].end(); iter++ )
					{
						int value = iter->first;

						int stripId = mStripContourMap[ camId2 ][ value ].contourStripId;
						int contourId = mStripContourMap[ camId2 ][ value ].contourId;


						stripEdgeIntersection( camId2 , contourId , stripId , epipole ,
							infinityEnds[ totalStripId ] ,
							intersectionPoint , coefficientPair );


						if( edgeStarted )
						{
							//                         Vector3d ray;
							// 
							//                         mCalibration->getRay( camId2 , intersectionPoint , ray );
							// 
							//                         double depth = tr::Math3D::lineIntersection( cameraCenter1 , raySet[ totalStripId ] , cameraCenter2 , ray );

							double depth = clipRayWithGenerator( cameraCenter1 , raySet[ totalStripId ] , camId2 , contourId , stripId );

							Edge &edge = viewingEdgeSet.back();

                            edge.point2 = cameraCenter1 + depth * raySet[ totalStripId ];

							edge.depth2 = depth;

							edge.camSecond = camId2;
							edge.contourSecond = contourId;
							edge.stripSecond = stripId;

							edgeStarted = false;

						}
						else
						{
							//                         Vector3d ray;
							// 
							//                         mCalibration->getRay( camId2 , intersectionPoint , ray );
							// 
							//                         double depth = tr::Math3D::lineIntersection( cameraCenter1 , raySet[ totalStripId ] , cameraCenter2 , ray );

							double depth = clipRayWithGenerator( cameraCenter1 , raySet[ totalStripId ] , camId2 , contourId , stripId );

							Edge edge ;

							edge.point1 = cameraCenter1 + depth * raySet[ totalStripId ];

							edge.depth1 = depth;

							int numStrips = mObjectContours[ camId1 ][ mStripContourMap[ camId1 ][ totalStripId ].contourId ].size();


							edgeStarted = true;

							edge.camLeft = camId1;
							edge.camRight = camId1;
							edge.contourRight = mStripContourMap[ camId1 ][ totalStripId ].contourId;
							edge.stripRight = mStripContourMap[ camId1 ][ totalStripId ].contourStripId;//
							edge.contourLeft = mStripContourMap[ camId1 ][ totalStripId ].contourId;
							edge.stripLeft = ( edge.stripRight - 1 + numStrips ) % numStrips;
							edge.camFirst = camId2;
							edge.contourFirst = contourId;
							edge.stripFirst = stripId;

							edge.leftFirst = -1;
							edge.leftSecond = -1;
							edge.rightFirst = -1;
							edge.rightSecond = -1;

							viewingEdgeSet.push_back( edge );


						}

					}

					assert( !edgeStarted );

					currentId = id;
				}
			}

			viewingEdges[ totalStripId ] = viewingEdgeSet;
		}

	}


	void EPVH::initiateViewingEdges2( int camId1, int camId2, std::vector< std::vector< Edge > >& viewingEdges)
	{
		double w = mCalibration->get_image_width( camId2 );
		double h = mCalibration->get_image_height( camId2 );

        Eigen::Vector2d epipole , refPoint( w / 2 , h / 2 );

		Point3d cameraCenter1 , cameraCenter2 , epipoleHmg1 , epipoleHmg2 ;

		mCalibration->get_camera_center( camId1 , cameraCenter1 );
		mCalibration->get_camera_center( camId2 , cameraCenter2 );

        mCalibration->projectPoint( cameraCenter1 , camId2 , epipole );

		int numContours = mObjectContours[ camId1 ].size();

		for( int contour = 0; contour < numContours; contour++ )
		{
			std::vector< cv::Point2f > &currentContour = mObjectContours[ camId1 ][ contour ];

			int numPoints = currentContour.size();

			for( int pp = 0; pp < numPoints; pp++ )
			{
				std::vector< tr::Edge > edges;

				Vector3d ray1 , infiniteProjHmg;

				ray1 = mGenerators[ camId1 ][ contour ][ pp ]->leftRay;

				//            mCalibration->getRay( camId1 , currentContour[ pp ] , ray1 );


				Vector4d infinitePoint;

				infinitePoint( 0 ) = ray1( 0 );
				infinitePoint( 1 ) = ray1( 1 );
				infinitePoint( 2 ) = ray1( 2 );
				infinitePoint( 3 ) = 0;

				// 	   mGenerators[ camId1 ][ contour ][ pp ]->leftRay = ray1;

				mGenerators[ camId1 ][ contour ][ pp ]->leftViewEdgeId = -1;
				mGenerators[ camId1 ][ contour ][ pp ]->rightViewEdgeId = -1;

				int idx = ( pp - 1 + numPoints ) % numPoints;

				//            mGenerators[ camId1 ][ contour ][ idx ]->rightRay = ray1;
				//            mGenerators[ camId1 ][ contour ][ idx ]->camId = camId1;

				mCalibration->projectHomogeneous( infinitePoint , camId2 , infiniteProjHmg  );

                Eigen::Vector2d infiniteEnd;

				infiniteEnd.x() = infiniteProjHmg( 0 ) / infiniteProjHmg( 2 );
				infiniteEnd.y() = infiniteProjHmg( 1 ) / infiniteProjHmg( 2 );

				std::vector< std::pair< int , int >  > strips;
                std::vector< std::pair< Eigen::Vector2d , Eigen::Vector2d > > clippedEdgesTemp;

                Eigen::Vector2d epipolef;

				epipolef.x() = epipole.x();
				epipolef.y() = epipole.y();

                std::pair< Eigen::Vector2d , Eigen::Vector2d > initialEdge;

				initialEdge.first = epipolef;
				initialEdge.second = infiniteEnd;

				clippedEdgesTemp.push_back( initialEdge );

				clippedEdgesTemp.clear();

				Eigen::Vector2d ray2D = infiniteEnd - epipolef;

				ray2D.normalize();

				// 	   if( pp > 100 )
				// 	     mClippers[ camId2 ].mDisplayInfo = true;
				// 	   else
				mClippers[ camId2 ].mDisplayInfo = false;

				mClippers[ camId2 ].clipRay( epipolef  , ray2D , clippedEdgesTemp , strips  );

				int numClippedEdges = clippedEdgesTemp.size();

				for( int ce = 0; ce < numClippedEdges; ce++ )
				{ 
					tr::Edge currentEdge;  

					currentEdge.camLeft = camId1;
					currentEdge.camRight = camId1;

					currentEdge.contourLeft = contour;
					currentEdge.contourRight = contour;

					currentEdge.stripLeft = ( pp - 1 + numPoints ) % numPoints;
					currentEdge.stripRight = pp;

					if( strips[ 2 * ce ].first != -1 )
					{
						currentEdge.camFirst = camId2;
						currentEdge.contourFirst = strips[ 2 * ce ].first;
						currentEdge.stripFirst = strips[ 2 * ce ].second;

						//                tr::Vector3d ray2;
						//        
						//                mCalibration->getRay( camId2 , clippedEdgesTemp[ ce ].first  , ray2 );
						// 
						//                double depth = Math3D::lineIntersection( cameraCenter2 , ray2 , cameraCenter1 , ray1 );
						// 
						//                currentEdge.point1 = cameraCenter2 + depth * ray2;  

						double depth = clipRayWithGenerator( cameraCenter1 , ray1 , camId2 , currentEdge.contourFirst , currentEdge.stripFirst );

						currentEdge.point1 = cameraCenter1 + depth * ray1;

					}


					if( strips[ 2 * ce + 1 ].first != -1 )
					{
						currentEdge.camSecond = camId2;
						currentEdge.contourSecond = strips[ 2 * ce + 1 ].first;
						currentEdge.stripSecond = strips[ 2 * ce + 1 ].second;

						double depth = clipRayWithGenerator( cameraCenter1 , ray1 , camId2 , currentEdge.contourSecond , currentEdge.stripSecond );

                        currentEdge.point2 = cameraCenter1 + depth * ray1;
					}


					edges.push_back( currentEdge );
				}

				if( edges.size() > 0 )
					viewingEdges.push_back( edges );



			}
		}
	}



	/**
	This method clips viewing edges of camera corresponding to camId1 agaist the silhouettes of
	the camera camId2.
	**/

	void EPVH::updateViewingEdges( int camId1, int camId2, std::vector< std::vector< tr::Edge > >& viewingEdges )
	{ 
		std::vector< std::vector< tr::Edge > > currentPairViewingEdges;

		std::vector< std::vector< tr::Edge > > intersectionEdges;

		int numEdgeSets = viewingEdges.size();

		for( int ee = 0; ee < numEdgeSets; ee++ )
		{
			int numEdges = viewingEdges[ ee ].size();

			std::vector< tr::Edge > clippedEdges;

			for( int ee2 = 0; ee2 < numEdges; ee2++ )
			{
				mClippers[ camId2 ].mDisplayInfo = false;

				tr::Edge &edge = viewingEdges[ ee ][ ee2 ];

				clipEdge( viewingEdges[ ee ][ ee2 ] , camId2 , clippedEdges );	
			}

			intersectionEdges.push_back( clippedEdges );
		}

		viewingEdges = intersectionEdges;

	}


	void EPVH::computeGeneratorNormals()
	{
		int numSilhouetteCams = mSilhouetteCameras.size();

		for( int sc = 0; sc < numSilhouetteCams; sc++ )
		{
			int sCamId = mSilhouetteCameras[ sc ];

			int numContours = mObjectContours[ sCamId ].size();

			for( int contour = 0; contour < numContours; contour++ )
			{
				int numStrips = mObjectContours[ sCamId ][ contour ].size();

				for( int strip = 0; strip < numStrips; strip++ )
				{
					tr::Vector3d& leftRay = mGenerators[ sCamId ][ contour ][ strip ]->leftRay;
					tr::Vector3d& rightRay = mGenerators[ sCamId ][ contour ][ strip ]->rightRay;

					mGenerators[ sCamId ][ contour ][ strip ]->normal = leftRay.cross( rightRay );

					mGenerators[ sCamId ][ contour ][ strip ]->normal.normalize();

				}
			}
		}
	}


	void EPVH::computeGeneratorNormals( int camId )
	{   
		int numContours = mObjectContours[ camId ].size();

		for( int contour = 0; contour < numContours; contour++ )
		{
			int numStrips = mObjectContours[ camId ][ contour ].size();

			for( int strip = 0; strip < numStrips; strip++ )
			{
				tr::Vector3d& leftRay = mGenerators[ camId ][ contour ][ strip ]->leftRay;
				tr::Vector3d& rightRay = mGenerators[ camId ][ contour ][ strip ]->rightRay;

				mGenerators[ camId ][ contour ][ strip ]->normal = leftRay.cross( rightRay );

				mGenerators[ camId ][ contour ][ strip ]->normal.normalize();

			}
		}
	}


	void EPVH::associateViewingEdgesToGenerators()
	{ 


		int numViewingEdgeSets = mEdges.size();

		for( int ee = 0; ee < numViewingEdgeSets; ee++ )
		{
			int camId = mEdges[ ee ].front().camLeft; 

			mGenerators[ camId ][ mEdges[ ee ].front().contourLeft ][ mEdges[ ee ].front().stripRight ]->leftViewEdgeId = ee;

			mGenerators[ camId ][ mEdges[ ee ].front().contourLeft ][ mEdges[ ee ].front().stripLeft ]->rightViewEdgeId = ee;        
		}

		int numCams = mCalibration->get_cam_count();

		for( int cc = 0; cc < numCams; cc++ )
		{
			associateViewingEdgesToGenerators( cc );
		}

	}

	void EPVH::associateViewingEdgesToGenerators( int camId )
	{
		int numViewingEdgeSets = mCameraViewingEdges[ camId ].size();

		for( int ee = 0; ee < numViewingEdgeSets; ee++ )
		{
			if( mCameraViewingEdges[ camId ][ ee ].size() == 0 )
				continue;

			int sCamId = mCameraViewingEdges[ camId ][ ee ].front().camLeft; 

			mGenerators[ camId ][ mCameraViewingEdges[ camId ][ ee ].front().contourLeft ][ mCameraViewingEdges[ camId ][ ee ].front().stripRight ]->leftViewEdgeId = ee;

			mGenerators[ camId ][ mCameraViewingEdges[ camId ][ ee ].front().contourLeft ][ mCameraViewingEdges[ camId ][ ee ].front().stripLeft ]->rightViewEdgeId = ee;   
		}
	}



	std::vector< tr::Edge > EPVH::computeViewingEdgeIntervalIntersection( std::vector< tr::Edge >& edgeSet1, 
		std::vector< tr::Edge >& edgeSet2  )
	{
		std::vector< tr::Edge > intersectedEdgeSet;

		if( edgeSet1.size() == 0 || edgeSet2.size() == 0 )
		{
			return intersectedEdgeSet;
		}


		int numEdges1 = edgeSet1.size();
		int numEdges2 = edgeSet2.size();

		for( int ee1 = 0; ee1 < numEdges1; ee1++ )
		{
			tr::Edge edge1 = edgeSet1[ ee1 ];

			assert( edge1.depth1 < edge1.depth2 );

			for( int ee2 = 0; ee2 < numEdges2; ee2++ )
			{
				tr::Edge &edge2 = edgeSet2[ ee2 ];

				assert( edge2.depth1 < edge2.depth2 );

				if( edge1.depth1 > edge2.depth2 || edge1.depth2 < edge2.depth1 )
					continue;

				if( edge1.depth1 > edge2.depth1 && edge1.depth1 < edge2.depth2 )
				{
					tr::Edge edge = edge1 ;

					if( edge1.depth2 < edge2.depth2 )
					{
						intersectedEdgeSet.push_back( edge );

						break;
					}
					else
					{
						edge.camSecond = edge2.camSecond;
						edge.contourSecond = edge2.contourSecond;
						edge.stripSecond = edge2.stripSecond;
                        edge.point2 = edge2.point2;
						edge.depth2 = edge2.depth2;

						intersectedEdgeSet.push_back( edge );
					}

				}
				else if( edge1.depth1 < edge2.depth1  )
				{
					tr::Edge edge;

					if( edge1.depth2 > edge2.depth1 && edge1.depth2 < edge2.depth2 )
					{
						edge = edge1;

						edge.camFirst = edge2.camFirst;
						edge.contourSecond = edge2.contourFirst;
						edge.stripFirst = edge2.stripFirst;
						edge.point1 = edge2.point1;
						edge.depth1 = edge2.depth1;
						intersectedEdgeSet.push_back( edge );		  

						break;
					}
					else if( edge1.depth2 > edge2.depth2 )
					{
						edge = edge2;
						intersectedEdgeSet.push_back( edge );		  
					}
				}


			}

		}

		return intersectedEdgeSet;
	}

	void EPVH::generatePolyData(std::vector< std::vector< tr::Edge > >& edges, vtkSmartPointer< vtkPolyData >& edgePolyData)
	{

		vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();

		colors->SetName("colors");
		colors->SetNumberOfComponents(3);

		int numEdges = 0;

		for (int ee1 = 0; ee1 < edges.size(); ee1++)
		{
			numEdges += edges[ee1].size();
		}

		vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();

		int numPoints = 2 * numEdges;

		std::vector< std::vector< cv::Vec3f > > colorVecs;

		colorVecs.resize(numEdges);

		points->Allocate(numPoints);

		int numCells = numEdges;

		edgePolyData->Allocate(numCells);

		vtkSmartPointer< vtkIdList > cell = vtkSmartPointer< vtkIdList >::New();

		for (int ee1 = 0; ee1 < edges.size(); ee1++)
			for (int ee2 = 0; ee2 < edges[ee1].size(); ee2++)
			{
				int id1 = points->InsertNextPoint(edges[ee1][ee2].point1.data());
				int id2 = points->InsertNextPoint(edges[ee1][ee2].point2.data());

				cell->Reset();

				cell->InsertNextId(id1);
				cell->InsertNextId(id2);

				edgePolyData->InsertNextCell(VTK_LINE, cell);

				Edge &edge = edges[ee1][ee2];

				if ((edge.camLeft == 15 && (edge.stripRight == 188)) || (edge.camLeft == 0 && (edge.stripRight == 48 || edge.stripRight == 49)))
				{
					const unsigned char _color[] = { 255, 0, 0 };

					colors->InsertNextTypedTuple(_color);
				}
				else
				{
					const unsigned char _color[] = { 0, 0, 255 };

					colors->InsertNextTypedTuple(_color);
				}
			}

		edgePolyData->SetPoints(points);

		edgePolyData->GetCellData()->SetScalars(colors);
	}




	bool EPVH::compute()
	{
		buildPrimitives();

		std::cout << " primitives built " << std::endl;

		computeGeneratorNormals();

		std::cout << " generator normals computed " << std::endl;

		buildClippers();

		std::cout << " clippers built " << std::endl;

		int numCams = mCalibration->get_cam_count();

		mCameraViewingEdges.resize( numCams );

		initiateWithViewingEdges();

		std::cout << " viewing edges initialized " << std::endl;

		associateViewingEdgesToGenerators();



		std::cout << " viewing edges associated with generators " << std::endl;

		initializeVertices();

		std::cout << " vertices initialized , statrting polygon computation " << std::endl;

		int id = 0;

		bool validConnections = true;

		while( 1 )
		{
			if( id >= mVertices.size() )
			{
				break;
			}

			tr::Vertex *vertex = mVertices[ id ];

			mVertexChainCounter = 0;

			completeVertex( vertex );

			if( !vertex->mLeft || !vertex->mRight )
			{
				validConnections = false;

				std::cout << "found invalid vertex connection " << std::endl;
			}

			if( mVertexChainCounter > 5000 )
			{
				validConnections = false;

				std::cout << "found large vertex chain counter : " << mVertexChainCounter << std::endl;

				break;
			}

			//	std::cout<<" incomplete left vertex "<< ( int )vertex->mIsGenVertex <<" "<<vertex->mGen->camId<<" "<<vertex->mGen->stripId <<std::endl;

			//if( !vertex->mRight )
			//	std::cout<<" incomplete right vertex "<< ( int )vertex->mIsGenVertex <<" "<<vertex->mGen->camId<<" "<<vertex->mGen->stripId<<std::endl;

			id++;
		}

		std::cout << " visual hull computed " << std::endl;

	

		return validConnections;  


	}


	bool EPVH::compute( std::vector< int > &silhouetteCameras , std::vector< int > &imageBoundaryCameras )
	{
		int numSilhouetteCameras = silhouetteCameras.size();

		int numImbCams = imageBoundaryCameras.size();

		//remove duplicate cameras
		if( ( numImbCams + numSilhouetteCameras ) < 2 )
		{
			return false;
		}

		for( int ss = 0; ss < numSilhouetteCameras ; ss++ )
		{
			mSilhouetteCameras.push_back( silhouetteCameras[ ss ] );
		}

		for( int bb = 0; bb < numImbCams ; bb++ )
		{
			mSilhouetteCameras.push_back( imageBoundaryCameras[ bb ] );
		}

		int numSilhouttes = mSilhouetteCameras.size();

		int numCams = mCalibration->get_cam_count();

		mScale.resize( numCams , 1.0f );
		mInvScale.resize( numCams , 1.0f );
		mBoundingBoxImages.resize( numCams );
		mGeneratorImages.resize( numCams );
		mOffset.resize( numCams );

		mIsOccludingContour.resize( numCams );
		mObjectContours.resize( numCams );
		mGenerators.resize( numCams );
		mStripContourMap.resize( numCams );
		mContourHierarchies.resize( numCams );

		buildMostOrthogonalCameras();

		for( int ss = 0; ss < numSilhouetteCameras ; ss++ )
		{
			buildPrimitives( mSilhouetteCameras[ ss ] );
		}

		for( int ss = numSilhouetteCameras; ss < mSilhouetteCameras.size() ; ss++ )
		{
			buildPrimitivesFromImageBoundary( mSilhouetteCameras[ ss ] );
		}

		computeGeneratorNormals();

		buildClippers();

		mCameraViewingEdges.resize( numCams );

		initiateWithViewingEdges();

		associateViewingEdgesToGenerators();   

		initializeVertices();

		int id = 0;

		bool validConnections = true;

		while( 1 )
		{
			if( id >= mVertices.size() )
			{
				break;
			}

			tr::Vertex *vertex = mVertices[ id ];

			mVertexChainCounter = 0;

			completeVertex( vertex );

			if( !vertex->mLeft || !vertex->mRight )
			{
				validConnections = false;
			}

			if( mVertexChainCounter > 5000 )
			{
				validConnections = false;
				break;
			}

			id++;
		}

		return validConnections;
	}


	bool EPVH::addSilhouetteCameraAndUpdate( int camId )
	{


		if( std::find( mSilhouetteCameras.begin() , mSilhouetteCameras.end() , camId ) == mSilhouetteCameras.end() )
		{

			//std::cout<<" adding next silhouette camera "<<std::endl;

			mSilhouetteCameras.push_back( camId );
		}
		else
		{
			//no computation needed
			return false;
		}

		//mCalibration->silhoutte( camId ).display_bounding_box_image();

		buildMostOrthogonalCameras();

		buildPrimitives( camId );

		computeGeneratorNormals( camId );

		buildClipper( camId );

		int numCams = mCalibration->get_cam_count();

		int orthogonalCam = mMostOrthogonalCamera[ camId ];

		std::vector< std::vector< Edge > > edges;

		initiateViewingEdges( camId , orthogonalCam , edges );

		int numSilhouetteCams = mSilhouetteCameras.size();

		for( int ss2 = 0; ss2 < numSilhouetteCams; ss2++ )
		{ 

			int sCamId2 = mSilhouetteCameras[ ss2 ];

			if( ( sCamId2 == camId ) || ( sCamId2 == orthogonalCam  ) )
				continue;     

			// std::cout<<" updating viewing edges : "<<camId<<" "<<

			updateViewingEdges( camId , sCamId2 , edges );

		}

		mCameraViewingEdges[ camId ] = edges;

		for( int ss = 0; ss < numSilhouetteCams; ss++ )
		{ 
			int sCamId = mSilhouetteCameras[ ss ];

			if( ( sCamId == camId )  )
				continue;     

			updateViewingEdges( sCamId , camId , mCameraViewingEdges[ sCamId ] );
		}

		associateViewingEdgesToGenerators( camId );



		clearVertices();

		initializeVertices();

		int id = 0;

		bool validConnections = true;

		while( 1 )
		{
			if( id >= mVertices.size() )
			{
				break;
			}

			tr::Vertex *vertex = mVertices[ id ];

			mVertexChainCounter = 0;

			completeVertex( vertex );

			if( !vertex->mLeft || !vertex->mRight )
			{
				validConnections = false;
			}

			if( mVertexChainCounter > 5000 )
			{
				validConnections = false;
				break;
			}

			//if( !vertex->mLeft  || !vertex->mRight )
			//	std::cout<<" vertex completion error "<<std::endl;

			id++;
		}

		return validConnections;

	}



	void EPVH::generatePolyData( std::vector< std::vector< cv::Point2f > >& contours , 
        std::vector< std::pair< Eigen::Vector2d , Eigen::Vector2d > > &edges,
		vtkSmartPointer< vtkPolyData >& edgePolyData )
	{

		vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();    

		colors->SetName("colors");
		colors->SetNumberOfComponents(3);


		int numEdges = 0 ;

		for( int ee1 = 0; ee1 < contours.size(); ee1++ )
		{
			numEdges += contours[ ee1 ].size();
		}

		vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();

		int numPoints = 2 * numEdges + 2 * edges.size();

		points->Allocate( numPoints );

		int numCells = numEdges + edges.size();

		edgePolyData->Allocate( numCells );

		vtkSmartPointer< vtkIdList > cell = vtkSmartPointer< vtkIdList >::New();

		for( int ee1 = 0; ee1 < contours.size(); ee1++ )
			for( int ee2 = 0; ee2 < contours[ ee1 ].size(); ee2++ ) 
			{

				cv::Point2f &pt1f = contours[ ee1 ][ ee2 ];
				cv::Point2f &pt2f = contours[ ee1 ][ ( ee2 + 1 ) % contours[ ee1 ].size() ];

                Eigen::Vector3d pt1( pt1f.x , pt1f.y , 0 ) , pt2( pt2f.x , pt2f.y , 0 );

				int id1 = points->InsertNextPoint( pt1.data() );
				int id2 = points->InsertNextPoint( pt2.data() );

				cell->Reset();

				cell->InsertNextId( id1 );
				cell->InsertNextId( id2 );   

				edgePolyData->InsertNextCell( VTK_LINE , cell );

				const unsigned char _color[] = { 255 , 255 , 255 };

				colors->InsertNextTypedTuple(_color);
			}


			for( int ee = 0; ee < edges.size(); ee++ )
			{

                Eigen::Vector2d &pt1f = edges[ ee ].first;
                Eigen::Vector2d &pt2f = edges[ ee ].second;

                Eigen::Vector3d pt1( pt1f.x() , pt1f.y() , 0 ) , pt2( pt2f.x() , pt2f.y() , 0 );

				int id1 = points->InsertNextPoint( pt1.data() );
				int id2 = points->InsertNextPoint( pt2.data() );

				cell->Reset();

				cell->InsertNextId( id1 );
				cell->InsertNextId( id2 );

				edgePolyData->InsertNextCell( VTK_LINE , cell );

				const unsigned char _color[] = { 255 , 255 , 255 };

				colors->InsertNextTypedTuple(_color);

			}

			edgePolyData->SetPoints( points );

			edgePolyData->GetCellData()->SetScalars( colors );


	}


	void EPVH::generatePolyData(
        std::vector< std::vector< cv::Point2f > >& contours, std::vector< std::pair< Eigen::Vector2d, Eigen::Vector2d > >& edges, vtkSmartPointer< vtkPolyData >& edgePolyData, std::vector< std::pair< int, int > >& selectedPolygonSegments, std::vector< cv::Vec3f >& selectedSegmentColors
		)
	{
		vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();    

		colors->SetName("colors");
		colors->SetNumberOfComponents(3);

		int numEdges = 0 ;

		for( int ee1 = 0; ee1 < contours.size(); ee1++ )
		{
			numEdges += contours[ ee1 ].size();
		}

		vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();

		int numPoints = 2 * numEdges + 2 * edges.size();

		points->Allocate( numPoints );

		int numCells = numEdges + edges.size();

		edgePolyData->Allocate( numCells );

		vtkSmartPointer< vtkIdList > cell = vtkSmartPointer< vtkIdList >::New();

		std::vector< std::vector< cv::Vec3f > > colorVecs;

		colorVecs.resize( contours.size() );

		for( int ee1 = 0; ee1 < contours.size(); ee1++ )
			for( int ee2 = 0; ee2 < contours[ ee1 ].size(); ee2++ ) 
			{
				colorVecs[ ee1 ].push_back( cv::Vec3f( 255 , 255 , 255 ) );  
			}

			int numSelectedPolygonSegments = selectedPolygonSegments.size();

			for( int ss = 0; ss < numSelectedPolygonSegments; ss++ )
			{
				int contourId = selectedPolygonSegments[ ss ].first;
				int segmentId = selectedPolygonSegments[ ss ].second;

				cv::Point2f &pt1f = contours[ contourId ][ segmentId ];
				cv::Point2f &pt2f = contours[ contourId ][ ( segmentId + 1 ) % contours[ contourId ].size() ];

				colorVecs[ contourId ][ segmentId ] = selectedSegmentColors[ ss ];
			}

			for( int ee1 = 0; ee1 < contours.size(); ee1++ )
				for( int ee2 = 0; ee2 < contours[ ee1 ].size(); ee2++ ) 
				{

					cv::Point2f &pt1f = contours[ ee1 ][ ee2 ];
					cv::Point2f &pt2f = contours[ ee1 ][ ( ee2 + 1 ) % contours[ ee1 ].size() ];

                    Eigen::Vector3d pt1( pt1f.x , pt1f.y , 0 ) , pt2( pt2f.x , pt2f.y , 0 );

					int id1 = points->InsertNextPoint( pt1.data() );
					int id2 = points->InsertNextPoint( pt2.data() );

					cell->Reset();

					cell->InsertNextId( id1 );
					cell->InsertNextId( id2 );

					edgePolyData->InsertNextCell( VTK_LINE , cell );

					const unsigned char _color[] = { (unsigned char)colorVecs[ ee1 ][ ee2 ][ 0 ] ,
                                           				 (unsigned char) colorVecs[ ee1 ][ ee2 ][ 1 ] ,
                                           				 (unsigned char) colorVecs[ ee1 ][ ee2 ][ 2 ]};

					colors->InsertNextTypedTuple(_color);
				}


				for( int ee = 0; ee < edges.size(); ee++ )
				{

                    Eigen::Vector2d &pt1f = edges[ ee ].first;
                    Eigen::Vector2d &pt2f = edges[ ee ].second;

                    Eigen::Vector3d pt1( pt1f.x() , pt1f.y() , 0 ) , pt2( pt2f.x() , pt2f.y() , 0 );

					int id1 = points->InsertNextPoint( pt1.data() );
					int id2 = points->InsertNextPoint( pt2.data() );

					cell->Reset();

					cell->InsertNextId( id1 );
					cell->InsertNextId( id2 );

					edgePolyData->InsertNextCell( VTK_LINE , cell );

					const unsigned char _color[] = { 255 , 255 , 255};

					colors->InsertNextTypedTuple(_color);
				}

				edgePolyData->SetPoints( points );

				edgePolyData->GetCellData()->SetScalars( colors );

	}


    void EPVH::display( std::vector< std::vector< cv::Point2f > >& contours, std::vector< std::pair< Eigen::Vector2d, Eigen::Vector2d > >& edges)
	{
		vtkSmartPointer< vtkPolyData > edgePolyData = vtkSmartPointer< vtkPolyData >::New();

		generatePolyData( contours , edges , edgePolyData );

		tr::Display3DRoutines::displayPolyData( edgePolyData );
	}

    void EPVH::displayWithColor( int camId, std::vector< std::vector< cv::Point2f > >& contours, std::vector< std::pair< Eigen::Vector2d, Eigen::Vector2d > >& edges)
	{
		vtkSmartPointer< vtkPolyData > edgePolyData = vtkSmartPointer< vtkPolyData >::New();

		generatePolyData( contours , edges , edgePolyData , mClippers[ camId ].mSelectedSegments , mClippers[ camId ].mColors );

		tr::Display3DRoutines::displayPolyData( edgePolyData );
	}



    double EPVH::clipWithGenerator(Eigen::Vector3d& p1, Eigen::Vector3d& p2, int camId, int contourId, int stripId)
	{

		tr::Generator *gen =  mGenerators[ camId ][ contourId ][ stripId ];

        Eigen::Vector3d &n = gen->normal;

        Eigen::Vector3d camCenter;

		mCalibration->get_camera_center( camId , camCenter );

		tr::Vector3d u = ( p2 - p1 );
		tr::Vector3d w = ( p1 - camCenter );

		double val1 = -n.dot( w );
		double val2 = n.dot( u );

		if( val2 != 0 )
		{
			return ( val1 / val2 );
		}
		else
		{
			return -1.0;
		}

	}


	double EPVH::clipRayWithGenerator( Point3d& o, Vector3d& r, int camId, int contourId, int stripId)
	{

		tr::Generator *gen =  mGenerators[ camId ][ contourId ][ stripId ];

        Eigen::Vector3d &n = gen->normal;

        Eigen::Vector3d camCenter;

		mCalibration->get_camera_center( camId , camCenter );

		tr::Vector3d w = ( o - camCenter );

		double val1 = -n.dot( w );
		double val2 = n.dot( r );

		if( val2 != 0 )
		{
			return ( val1 / val2 );
		}
		else
		{
			return -1.0;
		}

	}




	void EPVH::clear()
	{
		mSilhouetteCameras.clear();
		mMostOrthogonalCamera.clear();
		mIsCameraUsed.clear();
		mStripContourMap.clear();  
		mObjectContours.clear();
		mContourHierarchies.clear();
		mIsOccludingContour.clear();

		clearGenerators();

		mBoundingBoxImages.clear();
		mGeneratorImages.clear();
		mOffset.clear();
		mScale.clear();
		mInvScale.clear() ;

		mEdges.clear();
		mCameraViewingEdges.clear();

		clearVertices();
	}


	EPVH::~EPVH()
	{
		clear();
	}

}
