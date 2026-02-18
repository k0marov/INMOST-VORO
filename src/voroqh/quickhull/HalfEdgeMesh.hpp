#include "Structs/Mesh.hpp"

#ifndef HalfEdgeMesh_h
#define HalfEdgeMesh_h

namespace quickhull {
	
	template<typename FloatType, typename IndexType>
	class HalfEdgeMesh {
	public:
		HalfEdgeMesh() = default;
		
		struct HalfEdge {
			IndexType m_endVertex;
			IndexType m_opp;
			IndexType m_face;
			IndexType m_next;
		};
		
		struct Face {
			IndexType m_halfEdgeIndex; // Index of one of the half edges of this face
		};
		
		std::vector<Vector3<FloatType>> m_vertices;
		std::vector<Face> m_faces;
		std::vector<HalfEdge> m_halfEdges;
		
		HalfEdgeMesh(const MeshBuilder<FloatType>& builderObject, const VertexDataSource<FloatType>& vertexData,
                     std::vector<IndexType>& faceMapping, std::vector<IndexType>& halfEdgeMapping, std::vector<IndexType>& vertexMapping)
		{
            if (faceMapping.size() < builderObject.m_faces.size()) faceMapping.resize(builderObject.m_faces.size());
            std::fill(faceMapping.begin(), faceMapping.begin() + builderObject.m_faces.size(), std::numeric_limits<IndexType>::max());

            if (halfEdgeMapping.size() < builderObject.m_halfEdges.size()) halfEdgeMapping.resize(builderObject.m_halfEdges.size());
            std::fill(halfEdgeMapping.begin(), halfEdgeMapping.begin() + builderObject.m_halfEdges.size(), std::numeric_limits<IndexType>::max());

            if (vertexMapping.size() < vertexData.size()) vertexMapping.resize(vertexData.size());
            std::fill(vertexMapping.begin(), vertexMapping.begin() + vertexData.size(), std::numeric_limits<IndexType>::max());
			
			size_t i=0;
			for (const auto& face : builderObject.m_faces) {
				if (!face.isDisabled()) {
					m_faces.push_back({static_cast<IndexType>(face.m_he)});
					faceMapping[i] = m_faces.size()-1;
					
					const auto heIndices = builderObject.getHalfEdgeIndicesOfFace(face);
					for (const auto heIndex : heIndices) {
						const IndexType vertexIndex = builderObject.m_halfEdges[heIndex].m_endVertex;
						if (vertexMapping[vertexIndex] == std::numeric_limits<IndexType>::max()) {
							m_vertices.push_back(vertexData[vertexIndex]);
							vertexMapping[vertexIndex] = m_vertices.size()-1;
						}
					}
				}
				i++;
			}
			
			i=0;
			for (const auto& halfEdge : builderObject.m_halfEdges) {
				if (!halfEdge.isDisabled()) {
					m_halfEdges.push_back({static_cast<IndexType>(halfEdge.m_endVertex),static_cast<IndexType>(halfEdge.m_opp),static_cast<IndexType>(halfEdge.m_face),static_cast<IndexType>(halfEdge.m_next)});
					halfEdgeMapping[i] = m_halfEdges.size()-1;
				}
				i++;
			}
			
			for (auto& face : m_faces) {
				assert(halfEdgeMapping[face.m_halfEdgeIndex] != std::numeric_limits<IndexType>::max());
				face.m_halfEdgeIndex = halfEdgeMapping[face.m_halfEdgeIndex];
			}
			
			for (auto& he : m_halfEdges) {
				he.m_face = faceMapping[he.m_face];
				he.m_opp = halfEdgeMapping[he.m_opp];
				he.m_next = halfEdgeMapping[he.m_next];
				he.m_endVertex = vertexMapping[he.m_endVertex];
			}
		}
		
	};
}


#endif /* HalfEdgeMesh_h */
