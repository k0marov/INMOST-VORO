#include "planar_cutter.h"
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>

PlanarCutter::PlanarCutter() {}

Cell PlanarCutter::Cut(Cell& cell_to_cut, double a, double b, double c, double d, bool cut_positive) {
    if (!cell_to_cut.isValid()) return cell_to_cut;

    Mesh* m = cell_to_cut.GetMeshLink();
    const double epsilon = 1.0e-9;

    auto level_func = [&](const Storage::real p[3]) {
        return a * p[0] + b * p[1] + c * p[2] - d;
    };

    ElementArray<Face> original_faces = cell_to_cut.getFaces();
    ElementArray<Edge> original_edges = cell_to_cut.getEdges();
    
    std::map<HandleType, Node> edge_to_new_node;
    std::set<HandleType> on_plane_nodes;

    for (size_t i = 0; i < original_edges.size(); ++i) {
        Edge edge = original_edges[i];
        if (!edge.isValid()) continue;
        Node n0 = edge.getBeg(), n1 = edge.getEnd();
        double val0 = level_func(n0.Coords().data()), val1 = level_func(n1.Coords().data());
        if (std::fabs(val0) < epsilon) on_plane_nodes.insert(n0.GetHandle());
        if (std::fabs(val1) < epsilon) on_plane_nodes.insert(n1.GetHandle());
        if (val0 * val1 < -epsilon * epsilon) {
            double t = val0 / (val0 - val1);
            Storage::real p[3];
            const Storage::real* c0 = n0.Coords().data();
            const Storage::real* c1 = n1.Coords().data();
            p[0] = c0[0] + t * (c1[0] - c0[0]);
            p[1] = c0[1] + t * (c1[1] - c0[1]);
            p[2] = c0[2] + t * (c1[2] - c0[2]);
            Node new_node = m->CreateNode(p);
            edge_to_new_node[edge.GetHandle()] = new_node;
            Edge::SplitEdge(edge, ElementArray<Node>(m, 1, new_node.GetHandle()), 0);
        }
    }

    ElementArray<Edge> cutting_edges(m);
    for (size_t i = 0; i < original_faces.size(); ++i) {
        Face face = original_faces[i];
        if (!face.isValid()) continue;
        std::vector<HandleType> nodes_on_plane;
        ElementArray<Edge> face_edges = face.getEdges();
        for(size_t j = 0; j < face_edges.size(); ++j) {
            if (edge_to_new_node.count(face_edges[j].GetHandle())) {
                nodes_on_plane.push_back(edge_to_new_node[face_edges[j].GetHandle()].GetHandle());
            }
        }
        ElementArray<Node> face_nodes = face.getNodes();
        for(size_t j = 0; j < face_nodes.size(); ++j) {
            if (on_plane_nodes.count(face_nodes[j].GetHandle())) {
                nodes_on_plane.push_back(face_nodes[j].GetHandle());
            }
        }
        if (nodes_on_plane.size() == 2) {
            Edge new_edge = m->CreateEdge(ElementArray<Node>(m, nodes_on_plane.begin(), nodes_on_plane.end())).first;
            if (new_edge.isValid()) {
                cutting_edges.push_back(new_edge);
                Face::SplitFace(face, ElementArray<Edge>(m, 1, new_edge.GetHandle()), 0);
            }
        }
    }

    if (cutting_edges.size() < 3) return cell_to_cut;

    // 3. Order the cutting_edges to form a loop
    ElementArray<Edge> ordered_edges(m);
    if (!cutting_edges.empty()) {
        std::map<HandleType, std::vector<HandleType>> adj;
        for(size_t i = 0; i < cutting_edges.size(); ++i) {
            adj[cutting_edges[i].getBeg().GetHandle()].push_back(cutting_edges[i].GetHandle());
            adj[cutting_edges[i].getEnd().GetHandle()].push_back(cutting_edges[i].GetHandle());
        }

        std::set<HandleType> visited_edges;
        HandleType first_edge = cutting_edges[0].GetHandle();
        HandleType current_edge = first_edge;
        HandleType current_node = cutting_edges[0].getEnd().GetHandle();
        
        ordered_edges.push_back(Edge(m, current_edge));
        visited_edges.insert(current_edge);

        while(ordered_edges.size() < cutting_edges.size()) {
            bool found = false;
            for(HandleType edge_handle : adj[current_node]) {
                if(visited_edges.find(edge_handle) == visited_edges.end()) {
                    Edge e(m, edge_handle);
                    ordered_edges.push_back(e);
                    visited_edges.insert(edge_handle);
                    current_node = (e.getBeg().GetHandle() == current_node) ? e.getEnd().GetHandle() : e.getBeg().GetHandle();
                    found = true;
                    break;
                }
            }
            if (!found) break; // Should not happen in a simple loop
        }
    }
    
    if (ordered_edges.size() < 3) return cell_to_cut;

    Face cutting_face = m->CreateFace(ordered_edges).first;
    if (!cutting_face.isValid()) return cell_to_cut;

    ElementArray<Cell> new_cells = Cell::SplitCell(cell_to_cut, ElementArray<Face>(m, 1, cutting_face.GetHandle()), 0);

    if (new_cells.empty()) return InvalidCell();
    if (new_cells.size() == 1) return new_cells[0];

    Cell kept_cell = InvalidCell();
    for (size_t i = 0; i < new_cells.size(); ++i) {
        Storage::real centroid[3];
        new_cells[i].Centroid(centroid);
        double val = level_func(centroid);
        bool should_be_kept = cut_positive ? (val <= epsilon) : (val >= -epsilon);
        if (should_be_kept && !kept_cell.isValid()) {
            kept_cell = new_cells[i];
        } else {
            new_cells[i].Delete();
        }
    }
    
//    m->RemoveOrphaned(); this method does not exist
    return kept_cell;
}