#include "suVolume.h"
#include "triboxoverlap.h"
#include "suMathBase.h"
#include "overlap.h"
namespace SU {


	/*!
	* \function  LoadMeshFromFile
	*
	* \brief  load stl/obj/ply model by OpenMesh
	* \param  @pNode  [IN]current node
	*         @return  [OUT] if success(true/false).
	* \author Yuan Yao
	* \date 2016/6/16
	*/
	bool suVolume::LoadMeshFromFile(const char *pFileName)
	{
		if (!pFileName)
			return false;

		if (isLoad_)
		{
			clear();
			isLoad_ = false;
		}

		if (OpenMesh::IO::read_mesh(mesh_, pFileName))
		{
			suMesh::ConstVertexIter  v_it(mesh_.vertices_begin()),
				v_end(mesh_.vertices_end());

			bbMin_ = bbMax_ = mesh_.point(*v_it);
			for (; v_it != v_end; ++v_it)
			{
				bbMin_.minimize(mesh_.point(*v_it));
				bbMax_.maximize(mesh_.point(*v_it));
			}
			isLoad_ = true;
		}
		return isLoad_;
	}
	bool suVolume::LoadMeshFromMesh(suMesh & m)
	{
		if (isLoad_)
		{
			clear();
			isLoad_ = false;
		}

		mesh_ = m;
		{
			suMesh::ConstVertexIter  v_it(mesh_.vertices_begin()),
				v_end(mesh_.vertices_end());

			bbMin_ = bbMax_ = mesh_.point(*v_it);
			for (; v_it != v_end; ++v_it)
			{
				bbMin_.minimize(mesh_.point(*v_it));
				bbMax_.maximize(mesh_.point(*v_it));
			}
			isLoad_ = true;
		}
		return isLoad_;
	}

	void suVolume::clear()
	{
		if (pRoot_ != NULL)delete pRoot_;
		pRoot_ = NULL;
		mesh_.clear();
	}

	int suVolume::PartitionSpace(unsigned int level)
	{
		level_ = level;
		leafBoundaryNodes_.clear();
		leafInternalNodes_.clear();

		pRoot_ = new OctNode();
		curLevel_ = 0;

		//All mesh in bounding box, but we need to convenient recursivePartition.
		suMesh::VertexIter   v_it(mesh_.vertices_begin());
		suMesh::VertexIter   v_end(mesh_.vertices_end());
		suMesh::FaceIter     f_it(mesh_.faces_begin());
		suMesh::FaceIter     f_end(mesh_.faces_end());
		if (pRoot_)
		{
			//1.fill node
			pRoot_->max_.x = bbMax_[0];
			pRoot_->max_.y = bbMax_[1];
			pRoot_->max_.z = bbMax_[2];
			pRoot_->min_.x = bbMin_[0];
			pRoot_->min_.y = bbMin_[1];
			pRoot_->min_.z = bbMin_[2];

			///add face handle to node
			for (; f_it != f_end; ++f_it)
			{
				pRoot_->suNode_.AddElement(f_it.handle());
			}
			for (; v_it != v_end; ++v_it)
			{
				suMesh::Point   point = mesh_.point(v_it);
				pRoot_->suNode_.AddElement(v_it);
			}
			pRoot_->label_ = BOUNDARY_CELL;
			pRoot_->level_ = 0;

			if (level >= 9)
			{
				std::cout << "Don't support out-core voxilization\n";
				return 0;
			}

			//2.start partition		
			recursivePartition(pRoot_, level);

			//3.label internal & outernal
			LabelBoundaryNeighbors();

			//4.subdivide all node

			
			//5.foodfill
			floodFill();

		}

		return 0;
	}
	//recursiely partition for each level
	void suVolume::recursivePartition(OctNode * pParent, int nLevel)
	{
		if (pParent->level_ >= nLevel || pParent->label_ != BOUNDARY_CELL) return;

		//partition each boundary node
		Point delta = pParent->max_ - pParent->min_;
		Point halfLen = delta * 0.5f;
		Point mid = halfLen + pParent->min_;

		///new sub oct node
		for (int i = 0; i < 8; i++)
		{
			pParent->children_[i] = new OctNode();
			///fill node
			pParent->children_[i]->parent_ = pParent;
			pParent->children_[i]->level_ = pParent->level_ + 1;

			///generate code
			int xLocCode = i & 0x1;
			int yLocCode = (i >> 1) & 0x1;
			int zLocCode = (i >> 2) & 0x1;

			///update node size
			float minx = (xLocCode > 0) ? (pParent->min_.x + halfLen.x) : pParent->min_.x;
			float miny = (yLocCode > 0) ? (pParent->min_.y + halfLen.y) : pParent->min_.y;
			float minz = (zLocCode > 0) ? (pParent->min_.z + halfLen.z) : pParent->min_.z;

			float maxx = (xLocCode > 0) ? pParent->max_.x : (pParent->min_.x + halfLen.x);
			float maxy = (yLocCode > 0) ? pParent->max_.y : (pParent->min_.y + halfLen.y);
			float maxz = (zLocCode > 0) ? pParent->max_.z : (pParent->min_.z + halfLen.z);

			
			Point newMin = Point(minx, miny, minz);
			Point newMax = Point(maxx, maxy, maxz);

			pParent->children_[i]->max_ = newMax;
			pParent->children_[i]->min_ = newMin;

			pParent->children_[i]->xLocCode_ = (xLocCode << pParent->level_ + 1) | pParent->xLocCode_;
			pParent->children_[i]->yLocCode_ = (yLocCode << pParent->level_ + 1) | pParent->yLocCode_;
			pParent->children_[i]->zLocCode_ = (zLocCode << pParent->level_ + 1) | pParent->zLocCode_;
		}		

		//Label boundary & add faces to node	
		double triverts[3][3];
		double boxCenter[3];
		double boxHalfSize[3];

		for (int i = 0; i < 8; i++)
		{
			delta = pParent->children_[i]->max_ - pParent->children_[i]->min_;
			halfLen = delta * 0.5f;
			mid = pParent->children_[i]->min_ + halfLen;

			boxCenter[0] = mid.x;
			boxCenter[1] = mid.y;
			boxCenter[2] = mid.z;
			boxHalfSize[0] = halfLen.x;
			boxHalfSize[1] = halfLen.y;
			boxHalfSize[2] = halfLen.z;

			//Labeling boundary by faces
			size_t nFaces = pParent->suNode_.FaceVector.size();
			for (size_t k = 0; k < nFaces; k++)
			{
				suMesh::FaceHandle fHandle = pParent->suNode_.FaceVector[k];
				suMesh::FaceVertexIter fv_it = mesh_.fv_iter(fHandle);
				suMesh::Point    point;

				for (int j = 0; j < 3; j++)
				{
					point = mesh_.point(fv_it++);
					triverts[j][0] = point[0];
					triverts[j][1] = point[1];
					triverts[j][2] = point[2];
				}
				/*if (overlap(boxCenter, boxHalfSize, triverts) !=
					triBoxOverlap(boxCenter, boxHalfSize, triverts))
					std::cout << overlap(boxCenter, boxHalfSize, triverts) <<
					triBoxOverlap(boxCenter, boxHalfSize, triverts) << "|";*/
				if (triBoxOverlap(boxCenter, boxHalfSize, triverts))
				{
					pParent->children_[i]->suNode_.AddElement(fHandle);
				}
			}

			//intesect with faces | It is boundary node.
			if (!(pParent->children_[i]->suNode_.FaceVector.empty()))
			{
				pParent->children_[i]->label_ = BOUNDARY_CELL;

				if (nLevel == pParent->children_[i]->level_)
				{
					leafBoundaryNodes_.push_back(pParent->children_[i]);
					//std::cout << leafBoundaryNodes_.size() << std::endl;
				}

				//continue subdivide
				recursivePartition(pParent->children_[i], nLevel);
			}
		}

		return;
	}

	void suVolume::transverse(OctNode * pNode, Callback * cb)
	{
		
	}

	/*\brief 这个函数用于继续划分当前的节点，例如一些位于内部的节点，使其划分到与其它叶节点同样的尺寸，并进行标记
	 *\param 返回包含新叶节点的指针数组
	 *\example
	 * ```
	 *  std::vector<OctNode*> nodes;
	 * patitionToLevel(pNode, 5, pNode->label, nodes);
	 * ```
	 */
	void suVolume::patitionToLevel(OctNode * pNode, int nLevel, SU::NODE_LABEL label, std::vector<OctNode*> &nodeArr)
	{
		if (pNode->level_ >= nLevel || 
			pNode->children_[0] != 0 )    //has been partitioned, not a leafnode.
			return;

		//partition each boundary node
		Point delta = pNode->max_ - pNode->min_;
		Point halfLen = delta * 0.5f;
		Point mid = halfLen + pNode->min_;		

		///new sub oct node
		for (int i = 0; i < 8; i++)
		{
			pNode->children_[i] = new OctNode();
			///fill node
			pNode->children_[i]->parent_ = pNode;
			pNode->children_[i]->level_ = pNode->level_ + 1;

			///generate code
			int xLocCode = i & 0x1;
			int yLocCode = (i >> 1) & 0x1;
			int zLocCode = (i >> 2) & 0x1;

			///update node size
			float minx = (xLocCode > 0) ? (pNode->min_.x + halfLen.x) : pNode->min_.x;
			float miny = (yLocCode > 0) ? (pNode->min_.y + halfLen.y) : pNode->min_.y;
			float minz = (zLocCode > 0) ? (pNode->min_.z + halfLen.z) : pNode->min_.z;

			float maxx = (xLocCode > 0) ? pNode->max_.x : (pNode->min_.x + halfLen.x);
			float maxy = (yLocCode > 0) ? pNode->max_.y : (pNode->min_.y + halfLen.y);
			float maxz = (zLocCode > 0) ? pNode->max_.z : (pNode->min_.z + halfLen.z);


			Point newMin = Point(minx, miny, minz);
			Point newMax = Point(maxx, maxy, maxz);

			pNode->children_[i]->max_ = newMax;
			pNode->children_[i]->min_ = newMin;

			pNode->children_[i]->xLocCode_ = (xLocCode << pNode->level_ + 1) | pNode->xLocCode_;
			pNode->children_[i]->yLocCode_ = (yLocCode << pNode->level_ + 1) | pNode->yLocCode_;
			pNode->children_[i]->zLocCode_ = (zLocCode << pNode->level_ + 1) | pNode->zLocCode_;
		}
				
		if (pNode->children_[0]->level_ == nLevel)
		{
			//has been partitioned  to nLevel
			for (int i = 0; i < 8; i++)
			{
				pNode->children_[i]->label_ = label;
				nodeArr.push_back(pNode->children_[i]);
			}
			return;
		}
		else {
			//recursive partition
			for (int i = 0; i < 8; i++)
			{
				if (pNode->children_[i]->level_ < nLevel)
				{
					patitionToLevel(pNode->children_[i], nLevel, label, nodeArr);
				}
			}
		}
		
		return;
	}

	int suVolume::LabelBoundaryNeighbors()
	{
		for (int i = 0; i < leafBoundaryNodes_.size(); i++)
		{
			std::vector<OctNode*> neighbors;
			int N = get6NeighborNodes(leafBoundaryNodes_[i], neighbors);
			{
				for (int j = 0; j < N; j++)
				{
					if (neighbors[j]->label_ == UNDEFINE_CELL 
						//&& neighbors[j]->level_ == level_ //only for leaf voxels
						) 
					{
						neighbors[j]->label_ = labelNode(neighbors[j]);
						//neighbors[j]->label_ = INTERIOR_CELL;
						if (neighbors[j]->label_ == INTERIOR_CELL)
							leafInternalNodes_.push_back(neighbors[j]);
					}
				}
			}
		}
		return 0;
	}

	/*!
	* \function  labelNode
	*
	* \brief  search the nearest neighbor
	* \param  @pNode  [IN]current node
	*         @return  [OUT]a label.
	* \author Yuan Yao
	* \date 2016/6/16
	*/

	NODE_LABEL suVolume::labelNode(OctNode * pNode)
	{
		//return INTERIOR_CELL;
		Point delta = pNode->max_ - pNode->min_;
		Point halfLen = delta * 0.5f;
		Point mid = halfLen + pNode->min_;

		//find nearest face in neighbor nodes
		std::vector<OctNode*> neighbors;
		std::vector<suMesh::FaceHandle> faces;

		float fMinDist = 3.402823466e+38F;
		suMesh::FaceHandle faceHandle;
		int N = get6NeighborNodes(pNode, neighbors);
		{
			for (int i = 0; i < N; i++)
			{
				//we don't consider about parent node
				if (neighbors[i]->level_ == pNode->level_)
				{
					faces.insert(faces.end(), neighbors[i]->suNode_.FaceVector.begin(),
						neighbors[i]->suNode_.FaceVector.end());
				}
			}


			//for each face surrounding this node
			//find the nereast face
			int f_idx = 0;
			for (int i = 0; i < faces.size(); i++)
			{
				suMesh::FaceHandle fHandle = faces[i];
				suMesh::FaceVertexIter fv_it = mesh_.fv_iter(fHandle);
				suMesh::Point v;
				std::vector<Point> fv;

				for (int j = 0; j < 3; j++)
				{
					v = mesh_.point(fv_it++);
					fv.push_back(Point(v[0], v[1], v[2]));
				}

				Point closestPoint = ClosestPointOnTriangle(fv, mid);

				Point Delta_ = (closestPoint - mid);
				float Dist = Delta_.lenthSqrt();
				if (Dist < fMinDist)
				{
					fMinDist = Dist;
					faceHandle = fHandle;
				}
			}
		}

		suMesh::FaceVertexIter fv_it = mesh_.fv_iter(faceHandle);
		suMesh::Point v;
		std::vector<Point> fv;

		for (int j = 0; j < 3; j++)
		{
			v = mesh_.point(fv_it++);
			fv.push_back(Point(v[0], v[1], v[2]));
		}
		//compute normal			
		Point p1 = fv[2] - fv[0];
		Point p2 = fv[1] - fv[0];
		Point Normal_ = p1 * p2;
		Point closestPoint = ClosestPointOnTriangle(fv, mid);
		Point Delta_ = (closestPoint - mid);

		float fSign = Delta_.dot(Normal_);

		if (fSign < 0) return INTERIOR_CELL;
		else return EXTERIOR_CELL;
		return UNDEFINE_CELL;
	}
	void suVolume::labelNeighbors(OctNode *node)
	{
		std::vector<OctNode*> neighbors;
		int N = get6NeighborNodes(node, neighbors, false);
		for (int j = 0; j < N; j++)
		{
			if (neighbors[j]->label_ == UNDEFINE_CELL)
			{
				neighbors[j]->label_ = node->label_;
				leafInternalNodes_.push_back(neighbors[j]);
				labelNeighbors(neighbors[j]);
			}
		}

	}
	/*!
	* \function  floodFill
	* \brief  floodfill begin with
	*/
	void suVolume::floodFill()
	{
		std::vector<OctNode*> nodeArr = leafInternalNodes_;
		for (size_t i = 0; i < nodeArr.size(); i++)
		{
			std::vector<OctNode*> neighbors;
			int N = get6NeighborNodes(nodeArr[i], neighbors);
			{
				for (int j = 0; j < N; j++)
				{
					if (neighbors[j]->label_ == UNDEFINE_CELL)
					{
						neighbors[j]->label_ = INTERIOR_CELL;

						leafInternalNodes_.push_back(neighbors[j]);
						labelNeighbors(neighbors[j]);
					}
				}
			}

		}
	}

	/*!
	* \function  get6NeighborNodes
	*
	* \brief  return the 6-neighbors node.
	*         if one neighbor node is not exist, add its parent's (or parent's parent's ...) node to array if it exists
	* \param  @pNode  [IN]current node
	*         @nodes  [OUT]neighbor nodes array.
	*         @considerLevel  [IN] true(default):
	* \author Yuan Yao
	* \date 2016/5/25
	*/

	int suVolume::get6NeighborNodes(OctNode * pNode, std::vector<OctNode*>& nodes, bool considerLevel)
	{
		std::vector<OctNode *> neighbors;
		std::vector<Point> trans;
		trans.push_back(Point(1, 0, 0));
		trans.push_back(Point(-1, 0, 0));
		trans.push_back(Point(0, 1, 0));
		trans.push_back(Point(0, -1, 0));
		trans.push_back(Point(0, 0, 1));
		trans.push_back(Point(0, 0, -1));

		OctNode *pNNode_ = NULL;
		for (int i = 0; i < trans.size(); i++)
		{
			switch (i)
			{
				//if on boundary
			case 0: if (isOnTopBoundary(pNode->xLocCode_, pNode->level_)) continue; break;
			case 1: if (isOnBottomBoundary(pNode->xLocCode_, pNode->level_)) continue; break;
			case 2: if (isOnTopBoundary(pNode->yLocCode_, pNode->level_)) continue; break;
			case 3: if (isOnBottomBoundary(pNode->yLocCode_, pNode->level_)) continue; break;
			case 4: if (isOnTopBoundary(pNode->zLocCode_, pNode->level_)) continue; break;
			case 5: if (isOnBottomBoundary(pNode->zLocCode_, pNode->level_)) continue; break;
			}

			unsigned int xx = pNode->xLocCode_;
			unsigned int yy = pNode->yLocCode_;
			unsigned int zz = pNode->zLocCode_;
			if (trans[i].x != 0)
			{
				unsigned int inv_x = reverse_int(pNode->xLocCode_, pNode->level_ + 1);
				inv_x += trans[i].x;
				xx = reverse_int(inv_x, pNode->level_ + 1);
			}
			if (trans[i].y != 0)
			{
				unsigned int inv_y = reverse_int(pNode->yLocCode_, pNode->level_ + 1);
				inv_y += trans[i].y;
				yy = reverse_int(inv_y, pNode->level_ + 1);
			}
			if (trans[i].z != 0)
			{
				unsigned int inv_z = reverse_int(pNode->zLocCode_, pNode->level_ + 1);
				inv_z += trans[i].z;
				zz = reverse_int(inv_z, pNode->level_ + 1);
			}


			pNNode_ = getOctNode(pNode->level_, xx, yy, zz);
			if (pNNode_)
			{
				if (considerLevel)
				{
					if (pNNode_->level_ == pNode->level_)
					{
						nodes.push_back(pNNode_);
					}
				}
				else {
					nodes.push_back(pNNode_);
				}

			}
		}
		return nodes.size();
	}
	/*!
	* \function  getOctNode
	*
	* \brief  return a oct node.
	*         if this node is not exist,
	*             1)return its parent's (or parent's parent's ...) node
	*             2)it is not belong to octree, return null
	* \param  @level  locCode level
	*         @xLocCode, yLocCode, zLocCode  Location code.
	* \return @true: if a special-boundary nodes is found. false: if no such node is found.
	* \author Yuan Yao
	* \date 2016/5/25
	*/

	OctNode* suVolume::getOctNode(unsigned int level,
		unsigned int xLocCode, unsigned int yLocCode, unsigned int zLocCode)
	{
		if (!pRoot_) return NULL;

		OctNode *pNode = pRoot_;

		for (unsigned int i = 1; i <= level; i++)
		{
			unsigned int locCode = (unsigned int)1 << i;

			unsigned int index = ((xLocCode & locCode) >> i) |
				(((yLocCode & locCode) >> i) << 1) |
				(((zLocCode & locCode) >> i) << 2);

			if (pNode->children_[index] == NULL)
			{
				return pNode;
			}
			pNode = pNode->children_[index];
		}
		return pNode;
	}

	/*!
	* \function  getOctNode(x,y,z)
	*
	* \brief  return a a leaf oct node, where (x,y,z) belong to.*
	* \todo: out core compatibility
	* \author Yuan Yao
	* \date 2016/5/25
	*/
	OctNode* suVolume::getOctNode(float x, float y, float z)
	{
		Point p = Point(x, y, z);

		if (!pRoot_) return NULL;

		unsigned int xLocCode = 1;
		unsigned int yLocCode = 1;
		unsigned int zLocCode = 1;

		// If current point is not in the Boundary Box.
		if (p< pRoot_->min_ || p>pRoot_->max_)
			return NULL;

		OctNode *pCurNode = pRoot_;
		while (pCurNode)
		{

			// In current node, determine where the point belong to.	
			Point delta = pCurNode->max_ - pCurNode->min_;
			Point midPoint = pCurNode->min_ + delta / 2;
			xLocCode = (p.x >= midPoint.x) ? 1 : 0;
			yLocCode = (p.y >= midPoint.y) ? 1 : 0;
			zLocCode = (p.z >= midPoint.z) ? 1 : 0;

			OctNode *pNode = pCurNode->children_[((zLocCode << 2) | (yLocCode << 1) | xLocCode)];

			if (!pNode)
			{
				return pCurNode;
			}

			//update current node information	
			pCurNode = pNode;
		}

		return pCurNode;
	}


	//from http://www.gamedev.net/topic/552906-closest-point-on-triangle/
	Point suVolume::ClosestPointOnTriangle(std::vector<Point> fv, Point pos)
	{
		Point edge0 = fv[1] - fv[0];
		Point edge1 = fv[2] - fv[0];
		Point v0 = fv[0] - pos;

		float a = edge0.dot(edge0);
		float b = edge0.dot(edge1);
		float c = edge1.dot(edge1);
		float d = edge0.dot(v0);
		float e = edge1.dot(v0);

		float det = a*c - b*b;
		float s = b*e - c*d;
		float t = b*d - a*e;

		if (s + t < det)
		{
			if (s < 0.0f)
			{
				if (t < 0.0f)
				{
					if (d < 0.0f)
					{
						s = suFSaturate(-d / a);
						t = 0.0f;
					}
					else
					{
						s = 0.0f;
						t = suFSaturate(-e / c);
					}
				}
				else
				{
					s = 0.0f;
					t = suFSaturate(-e / c);
				}
			}
			else if (t < 0.0f)
			{
				s = suFSaturate(-d / a);
				t = 0.f;
			}
			else
			{
				float invDet = 1.0f / det;
				s *= invDet;
				t *= invDet;
			}
		}
		else
		{
			if (s < 0.0f)
			{
				float tmp0 = b + d;
				float tmp1 = c + e;
				if (tmp1 > tmp0)
				{
					float numer = tmp1 - tmp0;
					float denom = a - 2.0f*b + c;
					s = suFSaturate(numer / denom);
					t = 1.0f - s;
				}
				else
				{
					t = suFSaturate(-e / c);
					s = 0.0f;
				}
			}
			else if (t < 0.0f)
			{
				if (a + d > b + e)
				{
					float numer = c + e - b - d;
					float denom = a - 2.0f*b + c;
					s = suFSaturate(numer / denom);
					t = 1.0f - s;
				}
				else
				{
					s = suFSaturate(-e / c);
					t = 0.0f;
				}
			}
			else
			{
				float numer = c + e - b - d;
				float denom = a - 2.0f*b + c;
				s = suFSaturate(numer / denom);
				t = 1.0f - s;
			}
		}
		return fv[0] + edge0 * s + edge1 * t;
	}

	bool suVolume::saveVTK(const char *pVTKFileName, int level, const char *pVTKHead /*= "UnKnown Name"*/, float dx /*= 0*/, float dy /*= 0*/, float dz /*= 0*/)
	{
		if (NULL == pVTKFileName || NULL == pVTKHead)
			return false;


		std::ofstream vtk(pVTKFileName);

		if (!vtk)
			return false;

		vtk << "# vtk DataFile Version 2.0" << std::endl;
		vtk << pVTKHead << std::endl;
		vtk << "ASCII" << std::endl;
		vtk << "DATASET STRUCTURED_POINTS" << std::endl;

		SU::Point bMax(bbMax_[0], bbMax_[1], bbMax_[2]);
		SU::Point bMin(bbMin_[0], bbMin_[1], bbMin_[2]);
		SU::Point vSize = bMax - bMin;

		if (dx == 0)
		{
			float nSize = pow(2, level);
			dx = vSize.x / nSize;
			dy = vSize.y / nSize;
			dz = vSize.z / nSize;
		}

		unsigned int Xdim = (unsigned int)(vSize.x / dx);
		unsigned int Ydim = (unsigned int)(vSize.y / dy);
		unsigned int Zdim = (unsigned int)(vSize.z / dz);

		vtk << "DIMENSIONS " << Xdim << " " << Ydim << " " << Zdim << std::endl;

		// 数据区的其他信息
		vtk << "ASPECT_RATIO 1 1 1" << std::endl;
		vtk << "ORIGIN 0 0 0" << std::endl;
		vtk << "POINT_DATA " << Xdim * Ydim * Zdim << std::endl;
		vtk << "SCALARS volume_scalars char 1" << std::endl;
		vtk << "LOOKUP_TABLE default" << std::endl << std::endl;

		for (unsigned int IndexZ = 0; IndexZ < Zdim; IndexZ++)
		{
			float CellCenter_Z = bMin.z + IndexZ * dz + dz / 2.0f;

			for (unsigned int IndexY = 0; IndexY < Ydim; IndexY++)
			{
				float CellCenter_Y = bMin.y + IndexY * dy + dy / 2.0f;

				for (unsigned int IndexX = 0; IndexX < Xdim; IndexX++)
				{
					float CellCenter_X = bMin.x + IndexX * dx + dx / 2.0f;

					OctNode *pNode = getOctNode(CellCenter_X, CellCenter_Y, CellCenter_Z);
					if (pNode)
					{
#ifndef EXPORT_NODE_VALUE
						// 内部节点输出 0 
						if (INTERIOR_CELL == pNode->label_)
							vtk << 2 << " ";

						// 边界节点输出 1
						else if (BOUNDARY_CELL == pNode->label_ || BOUNDARY_CELL_SPECIAL == pNode->label_)
							vtk << 1 << " ";

						// 外部节点输出 2
						else if (EXTERIOR_CELL == pNode->label_)
							vtk << 3 << " ";

						// 其他情况输出 3 
						else
							vtk << 0 << " "; 
#else
						vtk << pNode->location << " ";
#endif
					}
					else
						vtk << 4 << " ";
				}
				vtk << std::endl;
			}
			vtk << std::endl;
		}

		return 0;
	}





	void SU::voxel_output::output_point1()
	{
		double q = pow(2, level) + 1;
		point_code[0] = (x_code - 1)*q + y_code + z_code *q*q;
		point_code[1] = (x_code - 1)*q + y_code + 1 + z_code *q*q;
		point_code[2] = x_code *q + y_code + 1 + z_code *q*q;
		point_code[3] = x_code *q + y_code + z_code *q*q;

		point_code[4] = (x_code - 1)*q + y_code + (z_code - 1)*q*q;
		point_code[5] = (x_code - 1)*q + y_code + 1 + (z_code - 1)*q*q;
		point_code[6] = x_code *q + y_code + 1 + (z_code - 1)*q*q;
		point_code[7] = x_code *q + y_code + (z_code - 1)*q*q;

		std::fstream outfile;
		outfile.open(address_, std::ios::app);
		for (int i = 0; i <= 7; i++)
		{
			outfile << point_code[i] << "  ";
		}
		//outfile <<"    "<< x_code << "  " << y_code << "  " << z_code;
		outfile << "    	mat 1 crossSect 1	nlgeo 1" << std::endl;

		outfile.close();
	}

	void SU::voxel_output::output_point2(int number_,std::vector<int> forcedPoint_)
	{
		double q = pow(2, level) + 1;
		point_code[0] = (x_code - 1)*q + y_code + z_code *q*q;
		point_code[1] = (x_code - 1)*q + y_code + 1 + z_code *q*q;
		point_code[3] = x_code *q + y_code + 1 + z_code *q*q;
		point_code[2] = x_code *q + y_code + z_code *q*q;

		point_code[4] = (x_code - 1)*q + y_code + (z_code - 1)*q*q;
		point_code[5] = (x_code - 1)*q + y_code + 1 + (z_code - 1)*q*q;
		point_code[7] = x_code *q + y_code + 1 + (z_code - 1)*q*q;
		point_code[6] = x_code *q + y_code + (z_code - 1)*q*q;

		std::fstream outfile;
		outfile.open(address_, std::ios::app);
		/*for (int i = 0; i <= 7; i++)
		{
		outfile << point_code[i] << "  ";
		}
		//outfile <<"    "<< x_code << "  " << y_code << "  " << z_code;
		outfile << "    	mat 1 crossSect 1	nlgeo 1" << endl;*/

		std::vector<int>::iterator forcedPointIt;
		outfile << "LTRSpace  " << 5 * number_ + 1 << "  " << "nodes  4	  " << point_code[0] << "  " << point_code[3] << "  "
			<< point_code[2] << "  " << point_code[6] << "  mat 1 crossSect 1	NIP  1";
		for (forcedPointIt = forcedPoint_.begin(); forcedPointIt != forcedPoint_.end(); forcedPointIt++)
		{
			if ((point_code[0] == *forcedPointIt) || (point_code[3] == *forcedPointIt) || (point_code[2] == *forcedPointIt) || (point_code[6] == *forcedPointIt))
			{
				outfile << " BoundaryLoads 2 2 3";
				break;
			}
		}
		outfile << std::endl;
		outfile << "LTRSpace  " << 5 * number_ + 2 << "  " << "nodes  4	  " << point_code[0] << "  " << point_code[3] << "  "
			<< point_code[6] << "  " << point_code[5] << "  mat 1 crossSect 1	NIP  1";
		for (forcedPointIt = forcedPoint_.begin(); forcedPointIt != forcedPoint_.end(); forcedPointIt++)
		{
			if ((point_code[0] == *forcedPointIt) || (point_code[3] == *forcedPointIt) || (point_code[5] == *forcedPointIt) || (point_code[6] == *forcedPointIt))
			{
				outfile << " BoundaryLoads 2 2 3";
				break;
			}
		}
		outfile << std::endl;
		outfile << "LTRSpace  " << 5 * number_ + 3 << "  " << "nodes  4	  " << point_code[3] << "  " << point_code[6] << "  "
			<< point_code[5] << "  " << point_code[7] << "  mat 1 crossSect 1	NIP  1";
		for (forcedPointIt = forcedPoint_.begin(); forcedPointIt != forcedPoint_.end(); forcedPointIt++)
		{
			if ((point_code[3] == *forcedPointIt) || (point_code[6] == *forcedPointIt) || (point_code[5] == *forcedPointIt) || (point_code[7] == *forcedPointIt))
			{
				outfile << " BoundaryLoads 2 2 3";
				break;
			}
		}
		outfile << std::endl;
		outfile << "LTRSpace  " << 5 * number_ + 4 << "  " << "nodes  4	  " << point_code[0] << "  " << point_code[1] << "  "
			<< point_code[3] << "  " << point_code[5] << "  mat 1 crossSect 1	NIP  1";
		for (forcedPointIt = forcedPoint_.begin(); forcedPointIt != forcedPoint_.end(); forcedPointIt++)
		{
			if ((point_code[0] == *forcedPointIt) || (point_code[1] == *forcedPointIt) || (point_code[3] == *forcedPointIt) || (point_code[5] == *forcedPointIt))
			{
				outfile << " BoundaryLoads 2 2 3";
				break;
			}
		}
		outfile << std::endl;
		outfile << "LTRSpace  " << 5 * number_ + 5 << "  " << "nodes  4	  " << point_code[4] << "  " << point_code[6] << "  "
			<< point_code[5] << "  " << point_code[0] << "  mat 1 crossSect 1	NIP  1";
		for (forcedPointIt = forcedPoint_.begin(); forcedPointIt != forcedPoint_.end(); forcedPointIt++)
		{
			if ((point_code[4] == *forcedPointIt) || (point_code[5] == *forcedPointIt) || (point_code[6] == *forcedPointIt) || (point_code[0] == *forcedPointIt))
			{
				outfile << " BoundaryLoads 2 2 3";
				break;
			}
		}
		outfile << std::endl;
		outfile.close();
	}



}//end namespace SU
