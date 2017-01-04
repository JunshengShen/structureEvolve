#include <windows.h>
#include <GL/gl.h>
#include "suWriteSTL.h"
#include <vector>
#include <algorithm>
#include <fstream>

namespace SU {

	class CmpVec
	{
	public:
		//TODO: double is needed.
		CmpVec(float _eps = 0.00000000001/*OpenMesh::flt_eps__*/) : eps_(_eps) {}

		bool operator()(const OpenMesh::Vec3f& _v0, const OpenMesh::Vec3f& _v1) const
		{
			if (fabs(_v0[0] - _v1[0]) <= eps_)
			{
				if (fabs(_v0[1] - _v1[1]) <= eps_)
				{
					return (_v0[2] < _v1[2] - eps_);
				}
				else return (_v0[1] < _v1[1] - eps_);
			}
			else return (_v0[0] < _v1[0] - eps_);
		}

	private:
		float eps_;
	};



	/*\brief function Write write_metaball_to_stl
	*       Write implicit metaball surface to a stl file
	*\param
	*        @threshold: threshold for metaball surface
	*        @nCubesPerDim: num of cubes per dimension.
	*        @bbox: a bounding box represented by a 2*3 matrix, matrix.row(0)=minPoint, matrix.row(1)=maxPoint
	*/

	bool SU::write_metaball_to_stl(const char * pFilename, std::vector<METABALL> metaballs, float threshold, int nCubesPerDim, std::vector<Eigen::Vector3f> bbox, suMesh originalMesh)
	{
		if (!metaballs.size()) return false;
		suGrid grid;
		grid.bbox.minP = bbox[0];
		grid.bbox.maxP = bbox[1];
		grid.Init(grid.bbox, nCubesPerDim);		//init with grids number per dimension
												//grid.GenField(metaballs);
		grid.GenFieldANN(metaballs);
		grid.GenSurface(threshold);
		grid.SaveStlFile(pFilename,originalMesh);
		return true;
	}

	bool suGrid::openMeshFile(const char * filename)
	{
		suMesh mesh;
		std::string strFileName = std::string(filename);
		OpenMesh::IO::read_mesh(mesh, strFileName);

		//convert obj from openmen to eigen
		suMesh::ConstVertexIter  v_it(mesh.vertices_begin()), v_end(mesh.vertices_end());

		suMesh::Point p;

		Eigen::DenseIndex rows = mesh.n_vertices();
		V.resize(rows, 3);
		rows = mesh.n_faces();
		F.resize(rows, 3);

		//add vertices
		for (; v_it != v_end; ++v_it)
		{
			p = mesh.point(v_it);
			V.row(v_it->idx()) << p[0], p[1], p[2];
		}
		//add face vertice 	
		suMesh::ConstFaceIter f_it(mesh.faces_begin()), f_end(mesh.faces_end());
		int idxFace = 0;
		for (; f_it != f_end; f_it++)
		{
			suMesh::FaceVertexIter fv_it = mesh.fv_begin(f_it.handle());
			suMesh::FaceVertexIter fv_end = mesh.fv_end(f_it.handle());

			int idxFv = 0;  //vertex index on a face
			for (; fv_it != fv_end; ++fv_it)
			{
				F(idxFace, idxFv++) = fv_it.handle().idx();

			}
			idxFace++;
		}


		//get bounding box
		bbox.minP = V.colwise().minCoeff();
		bbox.maxP = V.colwise().maxCoeff();

		isLoad_ = true;
		return isLoad_;
	}

	bool suGrid::Init(suBBox bbox, int num_cubes_per_dim)
	{
		numVertices = (num_cubes_per_dim + 1)*(num_cubes_per_dim + 1)*(num_cubes_per_dim + 1);
		numCubes = num_cubes_per_dim*num_cubes_per_dim*num_cubes_per_dim;

		///create vertices and cubes	
		vertices.resize(numVertices, 3);
		verticesV.resize(numVertices, 1);
		verticesN.resize(numVertices, 3);
		verticesV.setZero();
		verticesN.setZero();

		cubes.resize(numCubes, 8);

		int currentVertex = 0;

		float cube_size = (bbox.maxP - bbox.minP).maxCoeff() / num_cubes_per_dim;

		//init vertices
		for (int i = 0; i<num_cubes_per_dim + 1; i++)
		{
			for (int j = 0; j<num_cubes_per_dim + 1; j++)
			{
				for (int k = 0; k<num_cubes_per_dim + 1; k++)
				{
					Eigen::Vector3f  v1;
					v1 << i * cube_size, j * cube_size, k * cube_size;

					vertices.row(currentVertex) = bbox.minP + v1;

					currentVertex++;
				}
			}
		}

		int currentCube = 0;

		//set 8 indices of vertex for each cube
		for (int i = 0; i<num_cubes_per_dim; i++)
		{
			for (int j = 0; j<num_cubes_per_dim; j++)
			{
				for (int k = 0; k<num_cubes_per_dim; k++)
				{
					cubes.row(currentCube) <<
						(i*(num_cubes_per_dim + 1) + j)*(num_cubes_per_dim + 1) + k,
						(i*(num_cubes_per_dim + 1) + j)*(num_cubes_per_dim + 1) + k + 1,
						(i*(num_cubes_per_dim + 1) + (j + 1))*(num_cubes_per_dim + 1) + k + 1,
						(i*(num_cubes_per_dim + 1) + (j + 1))*(num_cubes_per_dim + 1) + k,
						((i + 1)*(num_cubes_per_dim + 1) + j)*(num_cubes_per_dim + 1) + k,
						((i + 1)*(num_cubes_per_dim + 1) + j)*(num_cubes_per_dim + 1) + k + 1,
						((i + 1)*(num_cubes_per_dim + 1) + (j + 1))*(num_cubes_per_dim + 1) + k + 1,
						((i + 1)*(num_cubes_per_dim + 1) + (j + 1))*(num_cubes_per_dim + 1) + k;

					currentCube++;
				}
			}
		}


		return false;
	}

	bool suGrid::GenField(std::vector<METABALL> metaballs)
	{
		std::vector<Eigen::Vector3f> points;
		for (int i = 0; i < metaballs.size(); i++)
		{
			points.push_back(metaballs[i].position);
		}
#ifndef USE_ANN	
		int numMetaballs = metaballs.size();
		for (int i = 0; i<numMetaballs; i++)
		{
			float squaredRadius = metaballs[i].squaredRadius;
			for (int j = 0; j<numVertices; j++)
			{
				Eigen::Vector3f positionV = vertices.row(j);
				Eigen::Vector3f ballToPoint = positionV - metaballs[i].position;

				//get squared distance from ball to point
				float squaredDistance = ballToPoint.dot(ballToPoint);

				if (squaredDistance == 0.0f)
					squaredDistance = 0.0001f;

				//value = r^2/d^2			
				verticesV.row(j)(0) += squaredRadius / squaredDistance;
				//normal = (r^2 * v)/d^4
				float normalScale = squaredRadius / (squaredDistance*squaredDistance);
				verticesN.row(j) += ballToPoint*normalScale;
			}
		}
#else	
		kdTree_.BuildTree(points);
		float fSearchDist = 0;
		if (metaballs.size()) {
			fSearchDist = sqrt(metaballs[0].squaredRadius) * 4;
		}
		for (int j = 0; j<numVertices; j++)
		{

			Eigen::Vector3f positionV = vertices.row(j);
			std::vector<unsigned int> resultArr_;

			kdTree_.WithinDistance(positionV, fSearchDist, resultArr_);

			for (int i = 0; i < resultArr_.size(); i++)
			{
				int idx_ = resultArr_[i];
				Eigen::Vector3f ballToPoint = positionV - metaballs[idx_].position;
				float squaredRadius = metaballs[idx_].squaredRadius;
				//get squared distance from ball to point
				float squaredDistance = ballToPoint.dot(ballToPoint);

				if (squaredDistance == 0.0f)
					squaredDistance = 0.0001f;

				//value = r^2/d^2			
				verticesV.row(j)(0) += squaredRadius / squaredDistance;
				//normal = (r^2 * v)/d^4
				float normalScale = squaredRadius / (squaredDistance*squaredDistance);
				verticesN.row(j) += ballToPoint*normalScale;
			}
		}
		kdTree_.FreeMemory();
#endif

		std::ofstream file("r:/debug.txt");
		file << verticesV << std::endl;
		return false;
	}

	bool suGrid::GenFieldANN(std::vector<METABALL> metaballs)
	{ 
		//build K-D tree for grid vertices
		std::vector<Eigen::Vector3f> points;		
		for (int i = 0; i < numVertices; i++)
		{
			Eigen::Vector3f _p = vertices.row(i);
			points.push_back(_p);
		}
		kdTree_.BuildTree(points);

		//update potential field near each metaball
		for (int i = 0; i < metaballs.size(); i++)
		{
			//ANN search
			float fSearchDist = sqrt(metaballs[0].squaredRadius) * 5;
			Eigen::Vector3f _p = metaballs[i].position;			
			std::vector<unsigned int> _resultArr;

			kdTree_.WithinDistance(_p, fSearchDist, _resultArr);
			
			//for each vertice within neighbor domain
			float squaredRadius = metaballs[i].squaredRadius;
			for (int j = 0; j < _resultArr.size(); j++)
			{
				int _idx = _resultArr[j];
				Eigen::Vector3f vP= vertices.row(_idx);
				Eigen::Vector3f ballToPoint = vP - _p;
				
				float squaredDistance = ballToPoint.dot(ballToPoint);

				if (squaredDistance == 0.0f)
					squaredDistance = 0.0001f;

				//value = r^2/d^2			
				verticesV.row(_idx)(0) += squaredRadius / squaredDistance;
				//normal = (r^2 * v)/d^4
				float normalScale = squaredRadius / (squaredDistance * squaredDistance);
				verticesN.row(_idx) += ballToPoint * normalScale;
			}
		}
		kdTree_.FreeMemory();

		return true;
	}

	void suGrid::GenSurface(float threshold)
	{
		std::ofstream file("r:/debug-mesh.txt");
		surface.clear();

		numFacesDrawn = 0;

		static SURFACE_VERTEX edgeVertices[12];

		CmpVec comp;
		std::map<OpenMesh::Vec3f, int, CmpVec>            vMap(comp);
		std::map<OpenMesh::Vec3f, int, CmpVec>::iterator  vMapIt;
		unsigned int idxV = 0;   //vertex index


		for (int i = 0; i < numCubes; i++)
		{
			unsigned char cubeIndex = 0;  //for edge table

			if (verticesV.row(cubes.row(i)(0))(0) < threshold)
				cubeIndex |= 1;
			if (verticesV.row(cubes.row(i)(1))(0) < threshold)
				cubeIndex |= 2;
			if (verticesV.row(cubes.row(i)(2))(0) < threshold)
				cubeIndex |= 4;
			if (verticesV.row(cubes.row(i)(3))(0) < threshold)
				cubeIndex |= 8;
			if (verticesV.row(cubes.row(i)(4))(0) < threshold)
				cubeIndex |= 16;
			if (verticesV.row(cubes.row(i)(5))(0) < threshold)
				cubeIndex |= 32;
			if (verticesV.row(cubes.row(i)(6))(0) < threshold)
				cubeIndex |= 64;
			if (verticesV.row(cubes.row(i)(7))(0) < threshold)
				cubeIndex |= 128;

			//look this value up in the edge table to see which edges to interpolate along
			int usedEdges = edgeTable[cubeIndex];
			//if the cube is entirely within/outside surface, no faces			
			if (usedEdges == 0 || usedEdges == 255)
				continue;


			for (int currentEdge = 0; currentEdge<12; currentEdge++)
			{
				if (usedEdges & 1 << currentEdge)
				{
					//if cut edge
					//get indice
					unsigned int idx_v1 = cubes.row(i)(verticesAtEndsOfEdges[currentEdge * 2]);
					unsigned int idx_v2 = cubes.row(i)(verticesAtEndsOfEdges[currentEdge * 2 + 1]);
					//interpolation
					float delta = (threshold - verticesV.row(idx_v1)(0)) / (verticesV.row(idx_v2)(0) - verticesV.row(idx_v1)(0));
					//vertex position
					edgeVertices[currentEdge].position = vertices.row(idx_v1) + delta * (vertices.row(idx_v2) - vertices.row(idx_v1));
					//vertex normal
					edgeVertices[currentEdge].normal = verticesN.row(idx_v1) + delta * (verticesN.row(idx_v2) - verticesN.row(idx_v1));
				}
			}

			//get the triangles
			for (int k = 0; triTable[cubeIndex][k] != -1; k += 3)
			{
				//vertices normal, not used
				edgeVertices[triTable[cubeIndex][k + 0]].normal;
				edgeVertices[triTable[cubeIndex][k + 2]].normal;
				edgeVertices[triTable[cubeIndex][k + 1]].normal;

				//triange vertices
				Eigen::Vector3f vtex[3];
				vtex[0] = edgeVertices[triTable[cubeIndex][k + 0]].position;
				vtex[1] = edgeVertices[triTable[cubeIndex][k + 2]].position;
				vtex[2] = edgeVertices[triTable[cubeIndex][k + 1]].position;

				//add a triangle
				OpenMesh::Vec3f fv;
				std::vector<OpenMesh::VertexHandle > fv_handles;


				for (int ii = 0; ii < 3; ii++)
				{
					for (int jj = 0; jj < 3; jj++)
					{
						fv[jj] = vtex[ii][jj];
					}

					// has vector been referenced before?
					if ((vMapIt = vMap.find(fv)) == vMap.end())
					{
						surface.add_vertex(fv);
						vMap[fv] = idxV;
						OpenMesh::VertexHandle fv_handle(idxV);
						fv_handles.push_back(fv_handle);
						idxV++;
					}
					else {
						OpenMesh::VertexHandle fv_handle(vMapIt->second);
						fv_handles.push_back(fv_handle);
					}
				}

				//for debug
				file << fv_handles[0].idx() << ", " << fv_handles[1].idx() << "," << fv_handles[2].idx() << "\n";

				surface.add_face(fv_handles);


			}
		}


	}

	/*void suGrid::DrawBBox()
	{
		glDisable(GL_LIGHTING);

		glBegin(GL_POINTS);
		{
			glColor3f(1.0f, 0.0f, 0.0f);
			for (int j = 0; j < numVertices; j++)
			{
				glVertex3f(vertices.row(j)(0), vertices.row(j)(1), vertices.row(j)(2));
			}
		}
		glEnd();

		glEnable(GL_LIGHTING);
	}

	void suGrid::DrawMesh()
	{
		glDisable(GL_LIGHTING);

		glBegin(GL_POINTS);
		{
			glColor3f(1.0f, 0.0f, 0.0f);
			for (int j = 0; j < V.rows(); j++)
			{
				glVertex3f(V.row(j)(0), V.row(j)(1), V.row(j)(2));
			}
		}
		glEnd();

		glEnable(GL_LIGHTING);
	}

	void suGrid::DrawSurface()
	{
		std::vector<unsigned int>  indices_;
		suMesh::ConstFaceIter   f_it(surface.faces_sbegin()), f_end(surface.faces_end());
		suMesh::ConstFaceVertexIter  fv_it;

		indices_.clear();
		indices_.reserve(surface.n_faces() * 3);

		for (; f_it != f_end; ++f_it)
		{
			indices_.push_back((fv_it = surface.cfv_iter(f_it)).handle().idx());
			indices_.push_back((++fv_it).handle().idx());
			indices_.push_back((++fv_it).handle().idx());
		}

		glDisable(GL_LIGHTING);
		glColor3f(1.0, 1.0, 1.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

		glEnableClientState(GL_VERTEX_ARRAY);

		int nFaces = indices_.size();
		if (nFaces > 0)
		{
			glVertexPointer(3, GL_FLOAT, 0, surface.points());
			glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);
		}

		glDisableClientState(GL_VERTEX_ARRAY);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
*/
	void suGrid::SaveStlFile(const char * filename,suMesh inputMesh)//add inside meshs(surface) to original meshs(inputMesh)
	{
		/*CmpVec comp;
		std::map<OpenMesh::Vec3f, int, CmpVec>            vMap(comp);
		std::map<OpenMesh::Vec3f, int, CmpVec>::iterator  vMapIt;
		unsigned int idxV = 0;   //vertex index

		suMesh::FaceIter face_it(surface.faces_begin()), face_end(surface.faces_end());
		suMesh test;
		test.clear();
		for (; face_it != face_end; face_it++)
		{
			//triange vertices
			Eigen::Vector3f vtex[3];
			suMesh::FaceVertexIter fv_it = surface.fv_iter(face_it.handle());
			int i = 0;
			for (; fv_it.is_valid(); ++fv_it) {
				//std::cout<<*fv_it;
				vtex[i][0] = surface.point(fv_it).data()[0];
				vtex[i][1] = surface.point(fv_it).data()[1];
				vtex[i][2] = surface.point(fv_it).data()[2];

				i++;

			}			

			//add a triangle
			OpenMesh::Vec3f fv;
			std::vector<OpenMesh::VertexHandle > fv_handles;


			for (int ii = 0; ii < 3; ii++)
			{
				for (int jj = 0; jj < 3; jj++)
				{
					fv[jj] = vtex[ii][jj];
				}

				// has vector been referenced before?
				if ((vMapIt = vMap.find(fv)) == vMap.end())
				{
					surface.add_vertex(fv);
					vMap[fv] = idxV;
					OpenMesh::VertexHandle fv_handle(idxV);
					fv_handles.push_back(fv_handle);
					idxV++;
				}
				else {
					OpenMesh::VertexHandle fv_handle(vMapIt->second);
					fv_handles.push_back(fv_handle);
				}
			}

			//for debug
		
			test.add_face(fv_handles);



		}*/






		//suMesh test;
		suMesh::FaceIter face_it(surface.faces_begin()), face_end(surface.faces_end());
		for (; face_it != face_end; ++face_it) {
			OpenMesh::Vec3f point[3];
			suMesh::FaceVertexIter fv_it = surface.fv_iter(face_it.handle());
			int i = 0;
			for (; fv_it.is_valid(); ++fv_it) {
				//std::cout<<*fv_it;
				point[i][0] = surface.point(fv_it).data()[0];
				point[i][1] = surface.point(fv_it).data()[1];
				point[i][2] = surface.point(fv_it).data()[2];

				//std::cout << point[i].data()[0] << point[i].data()[1] << point[i].data()[2] << std::endl;
				if (fabs(point[0].data()[0])>100 || fabs(point[i].data()[1])>100 || fabs(point[i].data()[2])>100) {
					std::cout << point[i].data()[0] << point[i].data()[1] << point[i].data()[2] << std::endl;
					system("pause");
				}
				i++;
			}
			
			//std::cout << point[0] << point[1] << point[2] << std::endl;
			

			//add vertex
			suMesh::VertexHandle vhandle[3];
			vhandle[0] = inputMesh.add_vertex(point[0]);
			vhandle[1] = inputMesh.add_vertex(point[2]);
			vhandle[2] = inputMesh.add_vertex(point[1]);

			//add face
			std::vector<suMesh::VertexHandle> face_vhandles;
			face_vhandles.clear();
			face_vhandles.push_back(vhandle[0]);
			face_vhandles.push_back(vhandle[1]);
			face_vhandles.push_back(vhandle[2]);
			inputMesh.add_face(face_vhandles);
		}
		OpenMesh::IO::write_mesh(inputMesh, std::string(filename));

	}

}//end of namespace SU