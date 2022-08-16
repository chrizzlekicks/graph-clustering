import java.util.Arrays;
import java.util.Stack;

// This file (original written by Sedgewick & Wayne, see below) is adapted from
// https://algs4.cs.princeton.edu/code/edu/princeton/cs/algs4/EdgeWeightedGraph.java
// for the use as AlgoDat programming exercise.
// Please also refer to the description on the exercise sheet.

/**
 *  The {@code EdgeWeightedGraph} class represents an edge-weighted
 *  graph of vertices named 0 through <em>V</em> – 1, where each
 *  undirected edge is of type {@link Edge} and has a real-valued weight.
 *  It supports the following two primary operations: add an edge to the graph,
 *  iterate over all of the edges incident to a vertex. It also provides
 *  methods for returning the number of vertices <em>V</em> and the number
 *  of edges <em>E</em>. Parallel edges and self-loops are permitted.
 *  By convention, a self-loop <em>v</em>-<em>v</em> appears in the
 *  adjacency list of <em>v</em> twice and contributes two to the degree
 *  of <em>v</em>.
 *  <p>
 *  This implementation uses an adjacency-lists representation, which 
 *  is a vertex-indexed array of {@link Bag} objects.
 *  All operations take constant time (in the worst case) except
 *  iterating over the edges incident to a given vertex, which takes
 *  time proportional to the number of such edges.
 *  <p>
 *  For additional documentation,
 *  see <a href="https://algs4.cs.princeton.edu/43mst">Section 4.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class EdgeWeightedGraph {
    private static final String NEWLINE = System.getProperty("line.separator");
    private double[][] coord; //saves the coordinates of your nodes/ data points. 
    private final int V;
    private int E;
    private Bag<Edge>[] adj;
    
    /**
     * Initializes an empty edge-weighted graph with {@code V} vertices and 0 edges.
     *
     * @param  V the number of vertices
     * @throws IllegalArgumentException if {@code V < 0}
     */
    public EdgeWeightedGraph(int V) {
        if (V < 0) throw new IllegalArgumentException("Number of vertices must be nonnegative");
        this.V = V;
        this.E = 0;
        adj = (Bag<Edge>[]) new Bag[V];
        for (int v = 0; v < V; v++) {
            adj[v] = new Bag<Edge>();
        }
    }
    
    /**
     * Initializes an empty edge-weighted graph from an input file
     *
     * @param  in input file - for the format see exercise sheet
     * @throws IllegalArgumentException if {@code V < 0}
     */
    public EdgeWeightedGraph(In in) {
        this(in.readInt());
        int E = in.readInt();
        int dim = in.readInt();
        coord = new double[V][dim];
        if (E < 0) throw new IllegalArgumentException("Number of edges must be nonnegative");
        for (int i = 0; i < E; i++) {
            int v = in.readInt();
            int w = in.readInt();
            validateVertex(v);
            validateVertex(w);
	        double weight=0;
	        double[] coordv = new double[dim];
	        double[] coordw = new double[dim];
	
	        for(int j=0; j<dim; j++) {
		        coordv[j]=in.readDouble();
		        coord[v][j]=coordv[j];
            }
	        for(int j=0; j<dim; j++) {
		        coordw[j]=in.readDouble();
		        coord[w][j]=coordw[j];
	        }
	        for(int j=0; j<dim; j++) {
		        weight= weight+Math.pow(coordv[j]-coordw[j],2);
	        }
	        weight=Math.sqrt(weight);
	        Edge e = new Edge(v, w, weight);
	        addEdge(e);
        }
    }
    
    /**
     * @return the coordinates of the nodes in the graph
     */
    public double[][] getCoordinates(){
        /* TODO */
        return this.coord;
    }
    
    public void setCoordinates(double[][] coord) {
        /* TODO */
        this.coord = coord;
    }

    /**
     * Initializes a new edge-weighted graph that is a deep copy of {@code G}.
     *
     * @param  G the edge-weighted graph to copy
     */
    public EdgeWeightedGraph(EdgeWeightedGraph G) {
        this(G.V());
        this.E = G.E();
        for (int v = 0; v < G.V(); v++) {
            // reverse so that adjacency list is in same order as original
            Stack<Edge> reverse = new Stack<Edge>();
            for (Edge e : G.adj[v]) {
                reverse.push(e);
            }
            for (Edge e : reverse) {
                adj[v].add(e);
            }
        }
    }


    /**
     * Returns the number of vertices in this edge-weighted graph.
     *
     * @return the number of vertices in this edge-weighted graph
     */
    public int V() {
        return V;
    }

    /**
     * Returns the number of edges in this edge-weighted graph.
     *
     * @return the number of edges in this edge-weighted graph
     */
    public int E() {
        return E;
    }

    // throw an IllegalArgumentException unless {@code 0 <= v < V}
    private void validateVertex(int v) {
        if (v < 0 || v >= V)
            throw new IllegalArgumentException("vertex " + v + " is not between 0 and " + (V-1));
    }

    /**
     * Adds the undirected edge {@code e} to this edge-weighted graph.
     *
     * @param  e the edge
     * @throws IllegalArgumentException unless both endpoints are between {@code 0} and {@code V-1}
     */
    public void addEdge(Edge e) {
        int v = e.either();
        int w = e.other(v);
        validateVertex(v);
        validateVertex(w);
        adj[v].add(e);
        adj[w].add(e);
        E++;
    }

    /**
     * Returns the edges incident on vertex {@code v}.
     *
     * @param  v the vertex
     * @return the edges incident on vertex {@code v} as an Iterable
     * @throws IllegalArgumentException unless {@code 0 <= v < V}
     */
    public Iterable<Edge> adj(int v) {
        validateVertex(v);
        return adj[v];
    }

    /**
     * Returns the degree of vertex {@code v}.
     *
     * @param  v the vertex
     * @return the degree of vertex {@code v}               
     * @throws IllegalArgumentException unless {@code 0 <= v < V}
     */
    public int degree(int v) {
        validateVertex(v);
        return adj[v].size();
    }

    /**
     * Returns all edges in this edge-weighted graph.
     * To iterate over the edges in this edge-weighted graph, use foreach notation:
     * {@code for (Edge e : G.edges())}.
     *
     * @return all edges in this edge-weighted graph, as an iterable
     */
    public Iterable<Edge> edges() {
        Bag<Edge> list = new Bag<Edge>();
        for (int v = 0; v < V; v++) {
            int selfLoops = 0;
            for (Edge e : adj(v)) {
                if (e.other(v) > v) {
                    list.add(e);
                }
                // add only one copy of each self loop (self loops will be consecutive)
                else if (e.other(v) == v) {
                    if (selfLoops % 2 == 0) list.add(e);
                    selfLoops++;
                }
            }
        }
        return list;
    }

    /**
     * Returns a string representation of the edge-weighted graph.
     * This method takes time proportional to <em>E</em> + <em>V</em>.
     *
     * @return the number of vertices <em>V</em>, followed by the number of edges <em>E</em>,
     *         followed by the <em>V</em> adjacency lists of edges
     */
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(V + " " + E + NEWLINE);
        for (int v = 0; v < V; v++) {
            s.append(v + ": ");
            for (Edge e : adj[v]) {
                s.append(e + "  ");
            }
            s.append(NEWLINE);
        }
        return s.toString();
    }


    public static void main(String[] args) {
		// TODO
        EdgeWeightedGraph ewg = new EdgeWeightedGraph(2);
        ewg.setCoordinates(new double[][]{{5, 6}, {8, 10}});
        System.out.println(Arrays.deepToString(ewg.getCoordinates()));
    }

}

