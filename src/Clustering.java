import java.util.*;
import java.awt.Color;

/**
 * This class solves a clustering problem with the Prim algorithm.
 */
public class Clustering {
	EdgeWeightedGraph G;
	List<List<Integer>> clusters;
	List<List<Integer>> labeled;
	
	/**
	 * Constructor for the Clustering class, for a given EdgeWeightedGraph and no labels.
	 * @param G a given graph representing a clustering problem
	 */
	public Clustering(EdgeWeightedGraph G) {
		this.G = G;
	    this.clusters = new LinkedList<>();
	}
	
    /**
	 * Constructor for the Clustering class, for a given data set with labels
	 * @param in input file for a clustering data set with labels
	 */
	public Clustering(In in) {
            int V = in.readInt();
            int dim= in.readInt();
            this.G = new EdgeWeightedGraph(V);
            this.labeled = new LinkedList<>();
            LinkedList labels = new LinkedList();
            double[][] coord = new double [V][dim];
            for (int v = 0;v<V; v++ ) {
                for(int j=0; j<dim; j++) {
                	coord[v][j]=in.readDouble();
                }
                String label= in.readString();
                    if(labels.contains(label)) {
                    	this.labeled.get(labels.indexOf(label)).add(v);
                    }
                    else {
                    	labels.add(label);
                    	List<Integer> l = new LinkedList<>();
                    	this.labeled.add(l);
                    	this.labeled.get(labels.indexOf(label)).add(v);
                    	System.out.println(label);
                    }                
            }
             
            this.G.setCoordinates(coord);
            for (int w = 0; w < V; w++) {
                for (int v = 0;v < V; v++ ) {
                	if (v != w) {
                	double weight = 0;
                    for(int j = 0; j < dim; j++) {
                    	weight += Math.pow(G.getCoordinates()[v][j]-G.getCoordinates()[w][j],2);
                    }
                    weight = Math.sqrt(weight);
                    Edge e = new Edge(v, w, weight);
                    G.addEdge(e);
                	}
                }
            }
	    this.clusters = new LinkedList<>();
	}
	
    /**
	 * This method finds a specified number of clusters based on a MST.
	 *
	 * It is based in the idea that removing edges from a MST will create a
	 * partition into several connected components, which are the clusters.
	 * @param numberOfClusters number of expected clusters
	 */
	public void findClusters(int numberOfClusters){
		/* TODO */
		// create a MST with given Graph
		PrimMST mst = new PrimMST(this.G);
		// create an empty list and loop over edges returned from MST
		LinkedList<Edge> edges = new LinkedList<>();
		for (Edge e : mst.edges()) {
			edges.add(e);
		}
		// sort edges and remove edges with the highest weight
		Collections.sort(edges);
		for (int i = 0; i < numberOfClusters-1; i++) {
			edges.removeLast();
		}
		// check for connected components and sort clusters
		connectedComponents(edges);
		for (List<Integer> cluster : clusters) {
			Collections.sort(cluster);
		}
	}
	
	/**
	 * This method finds clusters based on a MST and a threshold for the coefficient of variation.
	 *
	 * It is based in the idea that removing edges from a MST will create a
	 * partition into several connected components, which are the clusters.
	 * The edges are removed based on the threshold given. For further explanation see the exercise sheet.
	 *
	 * @param threshold for the coefficient of variation
	 */
	public void findClusters(double threshold){
		/* TODO */
		// create a MST with given Graph
		PrimMST mst = new PrimMST(this.G);
		// create an empty list and loop over edges returned from MST
		LinkedList<Edge> edges = new LinkedList<>();
		for (Edge e : mst.edges()) {
			edges.add(e);
		}
		// sort edges and remove edges where the returned coefficient is higher than the given threshold
		Collections.sort(edges);
		LinkedList<Edge> part = new LinkedList<>();
		for (Edge e : edges) {
			part.add(e);
			double coEff = coefficientOfVariation(part);
			if (coEff > threshold) {
				part.remove(e);
			}
		}
		// check for connected components and sort clusters
		connectedComponents(part);
		for (List<Integer> cluster : clusters) {
			Collections.sort(cluster);
		}
	}
	
	/**
	 * Evaluates the clustering based on a fixed number of clusters.
	 * @return array of the number of the correctly classified data points per cluster
	 */
	public int[] validation() {
		/* TODO */
		int[] clusterNum = new int[clusters.size()];
		for (int c = 0; c < clusters.size(); c++) {
			int count = 0;
			for (int n = 0; n < labeled.get(c).size(); n++) {
				if (clusters.get(c).contains(labeled.get(c).get(n))) {
					count++;
				}
			}
			clusterNum[c] = count;
		}
		return clusterNum;
	}
	
	/**
	 * Calculates the coefficient of variation.
	 * For the formula see exercise sheet.
	 * @param part list of edges
	 * @return coefficient of variation
	 */
	public double coefficientOfVariation(List <Edge> part) {
		/* TODO */
		if (part.isEmpty()) return 0.0;
		return (standardDeviation(part) / average(part));
	}

	/** Source: GeeksForGeeks */
	private double average(List<Edge> part) {
		/* TODO */
		double result = 0;
		for (Edge e : part) {
			result += e.weight();
		}
		return result / part.size();
	}

	/** Source: GeeksForGeeks */
	private double standardDeviation(List<Edge> part) {
		/* TODO */
		double result = 0;
		for (Edge e : part) {
			result += (e.weight() - average(part)) * (e.weight() - average(part));
		}
		return Math.sqrt(result / (part.size()));
	}

	public void connectedComponents(List<Edge> edges) {
		/* TODO */
		// create Union Find
		UF uf = new UF(G.V());
		// loop thru edges and identify representative
		for (Edge e : edges) {
			uf.union(e.either(), e.other(e.either()));
		}
		// create ArrayList to store unique representatives
		ArrayList<Integer> reps = new ArrayList<>();
		for (int i = 0; i < G.V(); i++) {
			int rep = uf.find(i);
			if (!reps.contains(rep)) {
				reps.add(rep);
			}
		}
		// create cluster sub lists based on unique representatives
		for (int j = 0; j < reps.size(); j++) {
			clusters.add(new ArrayList<>());
		}
		// add nodes based on their representatives
		for (int k = 0; k < G.V(); k++) {
			for (int r = 0; r < reps.size(); r++) {
				if (uf.find(k) == reps.get(r)) {
					clusters.get(r).add(k);
				}
			}
		}
	}
	
	/**
	 * Plots clusters in a two-dimensional space.
	 */
	public void plotClusters() {
		int canvas=800;
	    StdDraw.setCanvasSize(canvas, canvas);
	    StdDraw.setXscale(0, 15);
	    StdDraw.setYscale(0, 15);
	    StdDraw.clear(new Color(0,0,0));
		Color[] colors= {new Color(255, 255, 255), new Color(128, 0, 0), new Color(128, 128, 128), 
				new Color(0, 108, 173), new Color(45, 139, 48), new Color(226, 126, 38), new Color(132, 67, 172)};
	    int color=0;
		for(List <Integer> cluster: clusters) {
			if(color>colors.length-1) color=0;
		    StdDraw.setPenColor(colors[color]);
		    StdDraw.setPenRadius(0.02);
		    for(int i: cluster) {
		    	StdDraw.point(G.getCoordinates()[i][0], G.getCoordinates()[i][1]);
		    }
		    color++;
	    }
	    StdDraw.show();
	}

    public static void main(String[] args) {
		// FOR TESTING
//		In inp = new In("graph_small.txt");
		In inp = new In("iris_small.txt");
//		EdgeWeightedGraph G = new EdgeWeightedGraph(inp);
//		Clustering c = new Clustering(G);
		Clustering c = new Clustering(inp);
//		c.findClusters(3);
		c.findClusters(0.4);
		c.validation();
		c.plotClusters();
    }
}

