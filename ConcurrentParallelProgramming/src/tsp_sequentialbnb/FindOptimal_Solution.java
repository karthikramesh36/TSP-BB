/**
 * Traveling salesman problem implemented using the branch and bound method. 
 */
package tsp_sequentialbnb;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;


public class FindOptimal_Solution {
	
	// Adjacency matrix to store the weight/cost of the edges between each pair of vertices
	private int[][] adjacency_Matrix;
	
	private static Random random;
	
	// Flag to check if verbose mode is on
	private boolean isVerbose = true;
	
	// cities entered by a user are stored here
	private List<String> cities;
	
	// Checks if a route containing all the nodes is found or not.
	private boolean isPathTaken = false;
	
	// Stores the cost of the minimum cost route found
	private int minimumCost = Integer.MAX_VALUE;

	// Stores the node that contains the route having the minimum cost
	private Node minCostNode;

	// Nodes part of optimal solution are stored here.
	// Solution with lowest lower bound gets higher priority
	private PriorityQueue<Node> priority_Queue = new PriorityQueue<Node>();
	
	private long startTime;
	private long endTime;
	private int num_Cities;
	private int longestCityName = 0;
	
	/**
	 * Start the main program
	 * @param args
	 */
	public static void main(String[] args) {
		
		FindOptimal_Solution Solution = new FindOptimal_Solution();
				
		// read input graph file.
		Solution.File_Input();
		
		
		Solution.startTime = System.currentTimeMillis();
		System.out.println("Start time is : "+Solution.startTime);
		
		// An Edge_constraints to keep track of the included and excluded edges
		// This matrix is re-calculated for each node. Each node has its own
		// instance of this matrix in order to
		// efficiently track the inclusion/exclusion of edges at each node.
		int[][] Edge_constraints = new int[Solution.num_Cities][Solution.num_Cities];
		
		// Initializing edge matrix to have -1 between those vertices for which no edge exists
		for (int i=0; i<Solution.num_Cities; i++) {
			for (int j=0; j<Solution.num_Cities; j++) {
				if (Solution.adjacency_Matrix[i][j] == -1) Edge_constraints[i][j] = -1;
			}
		}

		// Display the adjacency matrix
		if (Solution.isVerbose) Solution.display_Matrix(Solution.adjacency_Matrix);		
		
		// Finding the lower bound at the root node
		double cost = 0.0;
		for (int v = 0; v < Solution.num_Cities; v++) {
			cost += Solution.findSumOfTwoMinimumCostEdges(v, Edge_constraints);
		}
		
		// Creating the root parent node that considers all the possible routes
		Node_Info value = new Node_Info(Edge_constraints, cost/2);
		Node root = new Node(value, null, null, null, false);

		System.out.println("Estimated lower bound cost at root node is : " + root.value.LowerBound_of);

		// Adding the root node to the queue
		Solution.priority_Queue.add(root);

		// Calling the FindOptimal_Solution function which processes the priority_Queue and finds the optimal route out of all the available routes
		Solution.FindOptimal_Solution(0, 1);		
		Solution.printOptimalTour(Solution.minCostNode);		
		Solution.endTime = System.currentTimeMillis();
		System.out.println("End time is : "+Solution.endTime);		
		System.out.println("Total time taken is : " + (Solution.endTime - Solution.startTime));
	}

	private void File_Input() {
		
		String graph = "tsp_graph_01.txt";
		File file = null;
		Scanner readGraph = null;
		try {
			URL url = getClass().getResource(graph);
			file = new File(url.getPath());
			readGraph = new Scanner(file);
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		String row = readGraph.nextLine();		
		this.num_Cities = row.split(" ").length;

		// Initialize the cities array size
		this.cities = new ArrayList<String>();
		
		// Add cities in the cities array.. starting from A
		for (int b = 0; b < this.num_Cities; b++) {
			this.cities.add((char)(b+65) + "");
		}
		
		this.adjacency_Matrix = new int[this.num_Cities][this.num_Cities];
		try {
			readGraph = new Scanner(file);
			int i=0;
			int j=0;
			while (readGraph.hasNextLine()) {
				row = readGraph.nextLine();
				String[] cols = row.split(" ");
				
				for (String col : cols) {
					int cost = Integer.parseInt(col.trim());
					this.adjacency_Matrix[i][j] = cost;
					j++;
				}
				//i = (i+1) % (this.cities.size()-1);
				i++;
				j = 0;
			}
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}	
		
		for (String name : this.cities) {
			this.longestCityName = Math.max(this.longestCityName, name.length());
		}	
	}

	 //Prints the minimum cost path taken and its cost	 
	private void printOptimalTour(Node node) {
		System.out.println("\n\nThe optimal path will include the below edges of the graph");
		List<String> edgesInvolvedInOptimalTour = new ArrayList<String>();
		int k=0;
		for (int i=0; i < node.value.Edge_constraints[0].length; i++) {
			for (int j=i+1; j < node.value.Edge_constraints[i].length; j++) {
				if (node.value.Edge_constraints[i][j] == 1) {
					edgesInvolvedInOptimalTour.add(this.cities.get(i) + this.cities.get(j));
					System.out.print(" "+this.cities.get(i) + this.cities.get(j));
				}
			}
		}
		System.out.println("\nOptimal Route Cost is : "+node.value.LowerBound_of);
	}

	/**
	 // check the cost of underlying leaf node
		if cost of leaf node is less than its own cost then continue and do
		not expand the tree

	//  check the number of incident edges considered till now
	    represented by i which shows the node index if the number of 
	    incident edges = 2 then break and set isLeafnode = true

	//  if only 2 edges are left and if there are no incident nodes on i
	    and if the calculate cost is less than its leaf node then increment i
	    and proceed to next city and repeat.
	 * @param i
	 * @param j
	 */
	private void FindOptimal_Solution(int i, int j) {		
		
		// Checks if the  priority queue still has some nodes to be processed
		while (!this.priority_Queue.isEmpty()) {
			
			// Remove the top element from the queue
			Node node = this.priority_Queue.remove();
			printNodeInfo(node);
			
			if (this.isPathTaken && node.value.LowerBound_of > minimumCost) {
				System.out.println("\nNode with lower bound "+node.value.LowerBound_of+" is pruned");
				continue;
			}
			
			if (i > this.num_Cities - 1 || j > this.num_Cities - 1) continue;
			
			// if there is no edge between these two vertices, then proceed to
			// next vertex
			if (this.adjacency_Matrix[i][j] == -1) {
				j++;
				continue;
			}
			
			// Keeping track of vertices whose edges are implicitly included/excluded due to the inclusion/exclusion of edge of other vertices
			Set<Integer> Constraints_Reload = new HashSet<Integer>();
			
			// Check if the left child of this node is null. If it is then create the left child
			if (node.left == null) {
				
				// Create the new constraints matrix for the left child
				int[][] Edge_constraints = new int[this.num_Cities][this.num_Cities];
				
				// Copies the constraints matrix from the current node to the left child node
				for (int loop = 0; loop < this.num_Cities; loop++) {
					Edge_constraints[loop] = Arrays.copyOf(node.value.Edge_constraints[loop], this.num_Cities);
				}

				// Considering the edge i,j and updating the state of the constraints matrix
				Edge_constraints[i][j] = Edge_constraints[j][i] = 1;
				updateEdge_constraintsMatrix(Edge_constraints, Constraints_Reload);

				// Updating the edge matrix for the edges present in the Constraints_Reload set
				for (int v : Constraints_Reload) {
					updateEdge_constraintsMatrix(Edge_constraints, v);
				}

				// Clearing the set as all the vertices are now reprocessed
				Constraints_Reload.clear();

				// Check if the route is found at this node
				boolean Path_Taken = isPathTaken(Edge_constraints);

				double PathCost = 0;

				// Update the minimum cost if route is found else find the lower bound and update it in the node
				if (Path_Taken) {
					PathCost = findPathCost(Edge_constraints);
				} else {
					for (int v = 0; v < this.num_Cities; v++) {
						PathCost += this.findSumOfTwoMinimumCostEdges(v, Edge_constraints);
					}
					PathCost = PathCost / 2.0;
				}
				Node_Info value = new Node_Info(Edge_constraints, PathCost);
				node.left = new Node(value, null, null, null, Path_Taken);

				// Prints all the edges considered and not considered for calculation of the lower bound on this node
				if (isVerbose) {
					System.out.println("\n");
					System.out.println("Computed the left node");
					printPathforthisnode(Edge_constraints);
					display_Matrix(Edge_constraints);
					if (Path_Taken) System.out.println("Cost of the route found : "+PathCost);
					else System.out.println("Lower bound of this node : "+PathCost);
				}

				// If a route is found with the minimum cost, update the minimum cost and store this node
				if (Path_Taken && PathCost < minimumCost) {
					this.minimumCost = (int) PathCost;
					minCostNode = node.left;
					this.isPathTaken = Path_Taken;
				} else {					
					// Adds this node in the queue since no route is found till now and so this node is a potential candidate for 
					// future processing.
					this.priority_Queue.add(node.left);
				}
			}
			
			// Check if the right child of this node is null. If yes, then create the right child
			if (node.right == null) {
				// Create the new edge matrix for the left child
				int[][] Edge_constraints = new int[this.num_Cities][this.num_Cities];
				
				// Copies the edge matrix from the current node to the left child node
				for (int loop = 0; loop < this.num_Cities; loop++) {
					Edge_constraints[loop] = Arrays.copyOf(node.value.Edge_constraints[loop], this.num_Cities);
				}

				// Considering the edge i,j and updating the state of the edge matrix
				Edge_constraints[i][j] = Edge_constraints[j][i] = -1;
				updateEdge_constraintsMatrix(Edge_constraints, Constraints_Reload);

				// Updating the edge matrix for the edges present in the Constraints_Reload set
				for (int v : Constraints_Reload) {
					updateEdge_constraintsMatrix(Edge_constraints, v);
				}

				// Clearing the set as all the vertices are now reprocessed
				Constraints_Reload.clear();

				// Check if the route is found at this node
				boolean Path_Taken = isPathTaken(Edge_constraints);

				double PathCost = 0;

				// Update the minimum cost if route is found else find the lower bound and update it in the node
				if (Path_Taken) {
					PathCost = findPathCost(Edge_constraints);
				} else {
					for (int v = 0; v < this.num_Cities; v++) {
						PathCost += this.findSumOfTwoMinimumCostEdges(v, Edge_constraints);
					}
					PathCost = PathCost / 2.0;
				}
				Node_Info value = new Node_Info(Edge_constraints, PathCost);
				node.right = new Node(value, null, null, null, Path_Taken);

				// Prints all the edges considered and not considered for calculation of the lower bound on this node
				if (isVerbose) {
					System.out.println("\n");
					System.out.println("Computed the right node");
					printPathforthisnode(Edge_constraints);
					display_Matrix(Edge_constraints);
					if (Path_Taken) System.out.println("Cost of the route found : "+PathCost);
					else System.out.println("Lower bound of this node : "+PathCost);
				}

				// If a route is found with the minimum cost, update the minimum cost and store this node
				if (Path_Taken && PathCost < minimumCost) {
					this.minimumCost = (int) PathCost;
					minCostNode = node.right;
					this.isPathTaken = Path_Taken;
				} else {					
					// Adds this node in the queue since no route is found till now and so this node is a potential candidate for 
					// future processing.
					this.priority_Queue.add(node.right);
				}
			}
			
			// Checks if all the edges incident with vertex i has been processed
			
			
			
			boolean check = false;
			Node topNodeInPQueue = this.priority_Queue.peek();
			for (int v=0; v < this.num_Cities; v++) {
				if (!isVertexDone(topNodeInPQueue.value.Edge_constraints, v)) {
					i = v;
					for (int col = i+1; col < this.num_Cities; col++) {
						if (topNodeInPQueue.value.Edge_constraints[i][col] == 0) {
							j = col;
							check = true;
							break;
						}
					}
				}
				if (check) break;
			}			
		}
	}

	private void printNodeInfo(Node node) {
		System.out.println("\nBelow node is now popped from the priority_Queue");
		display_Matrix(node.value.Edge_constraints);
		printPathforthisnode(node.value.Edge_constraints);
		System.out.println("Lower bound of this node is : "+node.value.LowerBound_of);
		System.out.println("Does this node has a left child? : "+(node.left != null));
		System.out.println("Does this node has a right child? : "+(node.right != null));
		System.out.println();
	}

	/**
	 * Prints all the edges that are included along with the ones that are not included for 
	 * calculation of lower bound on this node
	 * @param Edge_constraints
	 */
	private void printPathforthisnode(int[][] Edge_constraints) {
		StringBuilder includedEdges = new StringBuilder();
		StringBuilder not_includedEdges = new StringBuilder();
		for (int i = 0; i < this.num_Cities; i++) {
			for (int j = i; j < this.num_Cities; j++) {
				if (Edge_constraints[i][j] == 1) {
					includedEdges.append(" (" + this.cities.get(i) + "-" + this.cities.get(j) + ")");
				} else if (Edge_constraints[i][j] == -1) {
					not_includedEdges.append(" (" + this.cities.get(i) + "-" + this.cities.get(j) + ")");
				}
			}
		}
		System.out.println("Edges included :" + includedEdges.toString());
		System.out.println("Edges Not included :" + not_includedEdges.toString());
	}

	/**
	 * Checks if this vertex has been fully processed
	 * A vertex is fully processed until all the possible paths for 
	 * finding a route passing through that vertex are considered.
	 */
	private boolean isVertexDone(int[][] Edge_constraints, int i) {

		int remainingEdges = 0;
		int includedEdges = 0;
		int not_includedEdges = 0;

		for (int j = 0; j < this.num_Cities; j++) {
			if (i != j && Edge_constraints[i][j] == 0) {
				remainingEdges++;
			} else if (i != j && Edge_constraints[i][j] == -1) {
				not_includedEdges++;
			} else if (i != j && Edge_constraints[i][j] == 1) {
				includedEdges++;
			}
		}
		if (remainingEdges == 0 && not_includedEdges == this.num_Cities - 2 && includedEdges == 2) {
			return true;
		}
		return false;
	}

	// Finding the cost of the path taken.
	private int findPathCost(int[][] Edge_constraints) {
		int cost = 0;
		for (int i = 0; i < this.num_Cities; i++) {
			for (int j = i; j < this.num_Cities; j++) {
				if (Edge_constraints[i][j] == 1) {
					cost += adjacency_Matrix[i][j];
				}
			}
		}
		return cost;
	}

	/**
	 * Checks if a route has been found at this node.
	 * We have found a route if
	 * every vertex has exactly 2 edges incident with it
	 */
	private boolean isPathTaken(int[][] Edge_constraints) {

		for (int i = 0; i < this.num_Cities; i++) {

			// Tracks the total number of edges that should be incident with
			// this vertex i
			int remainingEdges = 0;
			int includedEdges = 0;
			int not_includedEdges = 0;

			for (int j = 0; j < this.num_Cities; j++) {
				if (i != j && Edge_constraints[i][j] == 0) {
					remainingEdges++;
				} else if (i != j && Edge_constraints[i][j] == -1) {
					not_includedEdges++;
				} else if (i != j && Edge_constraints[i][j] == 1) {
					includedEdges++;
				}
			}
			
			if (remainingEdges != 0 || not_includedEdges != this.num_Cities - 3 || includedEdges != 2) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Updates the constraints of all the edges incident with all the vertices in the graph.
	 * Edge_constraints[i,j] == -1 => this edge is not included in the tour
	 * Edge_constraints[i,j] == 1 => this edge must be included in the final tour
	 */
	private boolean updateEdge_constraintsMatrix(int[][] Edge_constraints, Set<Integer> s) {

		for (int i = 0; i < this.num_Cities; i++) {
			// Tracks the total number of edges that should be incident with
			// this vertex i
			int includedEdges = 0;

			// Tracks the remaining edges to be considered
			int remainingEdges = 0;

			for (int j = 0; j < this.num_Cities; j++) {
				if (i != j && Edge_constraints[i][j] == 1) {
					includedEdges++;
				} else if (i != j && Edge_constraints[i][j] == 0) {
					remainingEdges++;
				}
			}
			
			if (includedEdges > 2 || (includedEdges + remainingEdges) < 2) {
				return true;
			}
			
			if (includedEdges < 2) {
				if ((remainingEdges == 2 && includedEdges == 0) || (remainingEdges == 1 && includedEdges == 1)) {
					for (int j = 0; j < this.num_Cities; j++) {
						if (i != j && Edge_constraints[i][j] == 0) {
							Edge_constraints[i][j] = 1;
							Edge_constraints[j][i] = 1;
							s.add(j);
						}
					}
				}
			} else {
				for (int j = 0; j < this.num_Cities; j++) {
					if (i != j && Edge_constraints[i][j] == 0) {
						Edge_constraints[i][j] = -1;
						Edge_constraints[j][i] = -1;
						s.add(j);
					}
				}
			}
		}
		return false;
	}

	
	 // Updates the edge inclusion/exclusion constraints for all the edges incident with this vertex (v)
	 
	private void updateEdge_constraintsMatrix(int[][] Edge_constraints, int v) {

		// Tracks the total number of edges that should be incident with
		// this vertex i
		int includedEdges = 0;

		// Tracks the remaining edges to be considered
		int remainingEdges = 0;

		for (int j = 0; j < this.num_Cities; j++) {
			if (v != j && Edge_constraints[v][j] == 1) {
				includedEdges++;
			} else if (v != j && Edge_constraints[v][j] == 0) {
				remainingEdges++;
			}
		}
		if (includedEdges < 2) {
			if ((remainingEdges == 2 && includedEdges == 0) || (remainingEdges == 1 && includedEdges == 1)) {
				for (int j = 0; j < this.num_Cities; j++) {
					if (v != j && Edge_constraints[v][j] == 0) {
						Edge_constraints[v][j] = 1;
						Edge_constraints[j][v] = 1;
					}
				}
			}
		} else {
			for (int j = 0; j < this.num_Cities; j++) {
				if (v != j && Edge_constraints[v][j] == 0) {
					Edge_constraints[v][j] = -1;
					Edge_constraints[j][v] = -1;
				}
			}
		}
	}
	
	/**
	 * Returns the sum of the cost of two minimum cost edges incident with the
	 * vertex v. The selection of minimum cost edges is based on a certain
	 * heuristic. If we have already selected an edge to be included in the
	 * final tour, then that edge is always considered, regardless of whether it
	 * has a minimum cost or not.

	 */
	private int findSumOfTwoMinimumCostEdges(int v, int[][] Edge_constraints) {

		int firstMin = Integer.MAX_VALUE;
		int secondMin = Integer.MAX_VALUE;

		boolean isFirstMinFound = false;
		boolean isSecondMinFound = false;

		for (int i = 0; i < this.num_Cities; i++) {
			if (i != v && adjacency_Matrix[v][i] != -1 && Edge_constraints[v][i] != -1) {
				if (isFirstMinFound && isSecondMinFound) {
					return firstMin + secondMin;
				}
				if (!isFirstMinFound && Edge_constraints[v][i] == 1) {
					firstMin = adjacency_Matrix[v][i];
					isFirstMinFound = true;
					continue;
				} else if (!isSecondMinFound && Edge_constraints[v][i] == 1) {
					secondMin = adjacency_Matrix[v][i];
					isSecondMinFound = true;
					continue;
				}

				if (!isFirstMinFound && adjacency_Matrix[v][i] < firstMin) {
					secondMin = firstMin;
					firstMin = adjacency_Matrix[v][i];
				} else if (!isSecondMinFound && adjacency_Matrix[v][i] < secondMin) {
					secondMin = adjacency_Matrix[v][i];
				}
			}
		}
		return firstMin + secondMin;
	}
	
	private void generateMatrix() {		
		random = new Random();
		this.adjacency_Matrix = new int[this.num_Cities][this.num_Cities];
		int min = 1;

		for (int i = 0; i < this.num_Cities; i++) {
			for (int j = i + 1; j < this.num_Cities; j++) {
				this.adjacency_Matrix[i][j] = this.adjacency_Matrix[j][i] = random.nextInt(9) + min;
			}
		}
		for (int i = 0; i < this.num_Cities; i++) {
			for (int j = 0; j < this.num_Cities; j++) {
				System.out.print(this.adjacency_Matrix[i][j] + " ");
			}
			System.out.println();
		}
	}

	
	 //Displays the adjacency matrix representing the graph
	private void display_Matrix(int[][] matrix) {
		
		System.out.println("Printing the Reference edge matrix representing the graph.\n");		
		printSpaces(this.longestCityName + 2);
		for (int i = 0; i < this.cities.size(); i++) {
			System.out.print(this.cities.get(i) + " ");
		}

		System.out.println();

		for (int i = 0; i < this.cities.size(); i++) {
			System.out.print(this.cities.get(i));
			printSpaces(this.longestCityName - cities.get(i).length() + 2);
			
			for (int j = 0; j < this.cities.size(); j++) {
				System.out.print(matrix[i][j]);
				printSpaces(this.cities.get(i).length());
			}
			System.out.println();
		}
		System.out.println("\n");

	}

	private void printSpaces(int s) {
		for (int i = 0; i < s; i++) {
			System.out.print(" ");
		}
	}	
	
}

 // Edge_constraints and the lower bound calculated at each of the node stored below
 
class Node_Info {
	int[][] Edge_constraints;
	double LowerBound_of;

	public Node_Info(int[][] Edge_constraints, double cost) {
		super();
		this.Edge_constraints = Edge_constraints;
		this.LowerBound_of = cost;
	}
}

 //Represents the current node and its left and right children
class Node implements Comparable<Node> {
	Node right;
	Node parent;
	Node_Info value;
	Node left;
	boolean isLeafNode;

	public Node() {
	}

	public Node(Node_Info value, Node left, Node right, Node parent, boolean isLeafNode) {
		super();
		this.value = value;
		this.left = left;
		this.right = right;
		this.parent = parent;
		this.isLeafNode = isLeafNode;
	}

	@Override
	public int compareTo(Node node) {
		if (this.value.LowerBound_of < node.value.LowerBound_of) {
			return -1;
		} else if (this.value.LowerBound_of > node.value.LowerBound_of) {
			return 1;
		} else {
			return 0;
		}
	}
}
