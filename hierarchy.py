from utils import *
import numpy as np
from delaunay2D import Delaunay2D

def polygen_plot(P: Polygon, ocol="k-", show=True):
    for (p1, p2) in P.E:
        plt.plot([p1.x, p2.x], [p1.y, p2.y], ocol)
    if show:
        plt.show()

def polygens_plot(T: set[Polygon], show=False):
    for poly in T:
        polygen_plot(poly, show=False)

class Kirkpatrick:
    def __init__(self) -> None:
        self.previous = Polygon((-200, -200), (-300, -301), (-301, -300))
        self.parent = dict()
        self.adjacent = dict()
        self.node = set()
        self.polygons = list()
        self.user_map = dict()
        self.user_input = dict()
        self.delaunay = set()
        self.super_triangle = Polygon((150, 0),(-150, 150),(-150, -150))
    
    def add_region(self, poly: Polygon, name: str):
        self.user_input[poly] = name
        self.polygons.append(poly)
    
    def add_node(self, poly: Polygon):
        self.node.add(poly) 
        self.adjacent[poly] = []

    def retriangulate(self, poly: Polygon):
        triangulations = set()
        vertices = [v for v in poly.V]
        triangles_found = -1
        while (triangles_found != 0):
            triangles_found = 0
            for index, vertex in enumerate(vertices):
                next = vertices[(index+1)%len(vertices)]
                prev = vertices[index-1]

                v1 = next-vertex
                v2 = prev-vertex

                if (v1.x*v2.y - v1.y*v2.x > 0):
                    com_vertices = getcommon([vertex, prev, next])
                    T = Polygon(*com_vertices)
                    inside = [point_in_triangle(v, T) for v in vertices]
                    if inside.count(True) > 3:
                        continue
                    else:
                        vertices.pop(index)
                        triangulations.add(T)
                        triangles_found+=1
        
        return triangulations

    

    def helper(self, polys: list[Polygon]):
        super_t = self.super_triangle  # Assuming super_triangle is a Polygon object
        dt = Delaunay2D()

        # Extract all points from the input polygons
        points = []
        for poly in polys:
            points.extend([[vertex.x, vertex.y] for vertex in poly.V])  # Evaluate the generator and extend

        # Add points from the super_triangle
        points.extend([[vertex.x, vertex.y] for vertex in super_t.V])  # Evaluate and extend

        # Convert points to a NumPy array
        points = np.array(points)

        # Insert all points into the Delaunay2D object
        for point in points:
            dt.addPoint(point)

        # Get the Delaunay triangles
        delaunay_triangles = dt.exportTriangles()

        # Convert each triangle into a Polygon object
        triangulations = set()
        for triangle in delaunay_triangles:
            vertices = []
            for vertex in triangle:
                # Create Point objects for each vertex
                vertices.append(Point(points[vertex][0], points[vertex][1]))
            # Create a Polygon object for the triangle
            triangulations.add(Polygon(*vertices))

        return triangulations



    def triangulate(self):
        triangulations = self.helper(self.polygons)

        for poly in self.polygons:
            for tri in triangulations:
                vertices = tri.V
                flag = True
                for vertex in vertices:
                    if vertex not in poly.V: 
                        flag = False
                        break
                if flag: self.user_map[tri] = poly
        for tri in triangulations:
            self.node.add(tri)
            self.adjacent[tri] = []
        
        self.delaunay = triangulations
    
    def remove_point(self, p: Point):
        removed = set()
        for tri in self.delaunay:
            if p in tri.V: removed.add(tri)

        v_set = set()
        for tri in removed:
            for v in tri.V:
                v_set.add(v)
            self.delaunay.remove(tri)

        v_set.remove(p)
        v_list = getcommon(list(v_set), p)
        
        print("Independent points : ",v_list)

        new_triangles = self.retriangulate(Polygon(*v_list))

        for tri in new_triangles: 
            self.delaunay.add(tri)
            self.add_node(tri)
        for old in removed:
            self.parent[old] = []

        for new in new_triangles:
            for old in removed:
                self.adjacent[new].append(old)
                self.parent[old].append(new)

    def select_independent_set(self):
        all_vertices = {vertex for tri in self.delaunay for vertex in tri.V}
        ind_set = set()

        if len(all_vertices) >= 72:
            degree = {vertex: set() for vertex in all_vertices}
            for tri in self.delaunay:
                vertices = tri.V
                for i in range(3):
                    v1, v2 = vertices[i], vertices[(i+1) % 3]
                    degree[v1].add(v2)
                    degree[v2].add(v1)

            marked_set = {vertex for vertex in all_vertices if len(degree[vertex]) >= 9 or vertex in self.super_triangle.V}
        else:
            marked_set = {vertex for vertex in self.super_triangle.V}

        for vertex in all_vertices:
            if vertex in marked_set:
                continue
            ind_set.add(vertex)
            marked_set.add(vertex)
            for tri in self.delaunay:
                if vertex in tri.V:
                    for _vertex in tri.V: 
                        marked_set.add(_vertex)
        #print(ind_set)
        return ind_set


    def search_point(self, p: Point,triangulations):
        if not point_in_triangle(p, self.super_triangle):
            return "External Region"

        cur = None
        index = 0

        def plot_util(p:Point):
            nonlocal index
            plt.plot(p.x, p.y, 'ro')
            polygen_plot(self.super_triangle, 'b-', show=False)
            polygen_plot(triangle, 'r-', show=True)
            index+=1
        
        if self.previous != Polygon((-200, -200), (-300, -301), (-301, -300)) and point_in_triangle(p, self.previous):
            for triangle in self.parent[self.previous]:
                for t in self.adjacent[triangle]:
                    if point_in_triangle(p, t):
                        print("Found using parent of previous search!")
                        plot_util(p)

                        cur = triangle
                        return self.user_input.get(self.user_map[self.previous], "External Region")
                        # break
        for triangle in self.delaunay:
            if point_in_triangle(p, triangle):
                cur = triangle

                plot_util(p)

                break
    
        while self.adjacent[cur] != []:
            for triangle in self.adjacent[cur]:
                if point_in_triangle(p, triangle):
                    plot_util(p)

                    cur = triangle
                    
                    break
        
        plt.plot(p.x,p.y,'ro') 
        polygens_plot(triangulations[-1])

        if cur in self.user_map:
            self.previous = cur
            cur = self.user_map[cur]
            polygen_plot(cur, 'g-', show = False)
        plt.show()
        
        #print("ANSWER IS", cur)

        return self.user_input.get(cur, "External Region")
