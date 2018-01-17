// Utility functions for mesh processing

int SplitEdges(GSVMesh *mesh, RNLength max_edge_length = 0.2, int selected_segment = -1);
int CollapseEdges(GSVMesh *mesh, RNLength min_edge_length = 0.1, int selected_segment = -1);
int SmoothVertices(GSVMesh *mesh, RNScalar sigma = 0.1, int max_iterations = 1, int selected_segment = -1);
int FillHoles(GSVMesh *mesh, RNLength max_edge_length = 20);
