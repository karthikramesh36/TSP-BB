
package tsp_threadedbnb;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Scanner;
import java.util.Set;

/**
 * A parallel approach to the traveling salesman problem using threads
 */
public class FindOptimal_Solution {
	
	// To check if the verbose mode is on or off
	private boolean isVerbose = true;
	
	// List to store all the cities entered by a user
	private List<String> cities;
	
	// Adjacency matrix to store the weight/cost of the edges between each pair of vertices
	private int[][] adjacency_Matrix;
	
	// Stores the total number of cities
	private int num_cities;

	private int longestCityName = 0;

	// Checks if a route containing all the nodes is found or not.
	private volatile boolean isPath_taken = false;
	
	// Stores the cost of the minimum cost route found
	private volatile int Current_minimumCost = Integer.MAX_VALUE;

	// Stores the node that contains the route having the minimum cost
	private volatile Node minCostNode;

	// Keeps the nodes that are created during the process of finding an optimal route
	// Priority is given to those nodes having the minimum of all the lower bounds
	private PriorityQueue<Node> priority_queue = new PriorityQueue<Node>();

	// Creating a lock to provide mutually exclusive access the share data structures
	private Object lock1 = new Object();
	private Object lock2 = new Object();
	
	// A counter object that keeps track of the number of waiting threads
	private volatile int count;
	
	// Variable that stores the total number of threads started by the system
	private int Total_Threads = 2;
	
	private long startTime;
	private long endTime;
	
	
	//Start the main program
	public static void main(String[] args) {
		
		final FindOptimal_Solution Solution = new FindOptimal_Solution();
		
		// Reads the input from a file
		Solution.read_inputfile();
		
		Solution.startTime = System.currentTimeMillis();
		System.out.println("Start time is : "+Solution.startTime);
		
		// An constraintsMatrix to keep track of the included and excluded edges
		// Each node has its own instance of this matrix in order to
		// efficiently track the inclusion/exclusion of edges at each node.
		int[][] constraintsMatrix = new int[Solution.num_cities][Solution.num_cities];
		
		// Initializing edge matrix to have -1 between those vertices for which no edge exists
		for (int i=0; i<Solution.num_cities; i++) {
			for (int j=0; j<Solution.num_cities; j++) {
				if (Solution.adjacency_Matrix[i][j] == -1) constraintsMatrix[i][j] = -1;
			}
		}

		// Display the adjacency matrix
		if (Solution.isVerbose) Solution.display_Matrix(Solution.adjacency_Matrix);		
		
		// Finding the lower bound at the root node
		double cost = 0.0;
		for (int v = 0; v < Solution.num_cities; v++) {
			cost += Solution.findSumOfTwoMinimumCostEdges(v, constraintsMatrix);
		}
		
		// Creating the root node that considers all the possible routes
		Node_Info value = new Node_Info(constraintsMatrix, cost/2);
		Node root = new Node(value, null, null, null, false);

		System.out.println("Estimated lower bound cost at root node is : " + root.value.lowerbound_of);

		// Adding the root node to the queue
		Solution.priority_queue.add(root);	
		
		// Create the desired number of threads as entered by the user
		/*Thread th[] = new Thread[Solution.Total_Threads];
		for (int t=0; t < Solution.Total_Threads; t++) {
			th[t] = new Thread(new Runnable() {			
				@Override
				public void run() {
					try {
						Solution.FindOptimal_Solution();
					} catch (InterruptedException e) {
						System.out.println("Thread is interrupted");
					}				
				}
			});
		}*/
		
		// First Thread
		Thread t1 = new Thread(new Runnable() {			
			@Override
			public void run() {
				try {
					Solution.FindOptimal_Solution();
				} catch (InterruptedException e) {
					System.out.println("Thread is interrupted");
					//e.printStackTrace();
				}				
			}
		});
		
		// Second Thread
		Thread t2 = new Thread(new Runnable() {			
			@Override
			public void run() {
				try {
					Solution.FindOptimal_Solution();
				} catch (InterruptedException e) {
					System.out.println("Thread is interrupted");
				}
			}
		});
		
		// Third Thread
		/*Thread t3 = new Thread(new Runnable() {
			@Override
			public void run() {
				try {
					Solution.FindOptimal_Solution();
				} catch (InterruptedException e) {
					System.out.println("Thread is interrupted");
				}
			}
		});*/
		
		try {
			// Start the threads
			/*for (int t=0; t < Solution.Total_Threads; t++) {
				th[t].start();
				th[t].join();
			}*/
			t1.start();
			t2.start();
			//t3.start();
			
			t1.join();
			t2.join();
			//t3.join();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		Solution.PrintOptimalRoute(Solution.minCostNode);
		
		Solution.endTime = System.currentTimeMillis();
		System.out.println("End time is : "+Solution.endTime);
		
		System.out.println("Total time taken is : " + (Solution.endTime - Solution.startTime));
	}
	
	private void FindOptimal_Solution() throws InterruptedException {
		
		/**********************   START of MULTI-THREADED CODE   *************************/
		
		while (count < Total_Threads) {
			
			Node firstNodeinPQ = null;
			
			synchronized (lock1) {
				while (this.priority_queue.isEmpty()) {
					if (count >= Total_Threads) {
						lock1.notifyAll();
						return;
					} else {
						count++;
						System.out.println("Number of threads waiting is : " + count);
						if (count >= Total_Threads) {
							lock1.notifyAll();
							return;
						}
						lock1.wait();
					}
					
					/*System.out.println("Number of threads waiting is : "+count);
					if (count >= Total_Threads) {
						lock1.notifyAll();
						return;
					}
					lock1.wait();*/
				}		
				
				//System.out.println("Decremented Count is : "+count);
				firstNodeinPQ = this.priority_queue.remove();
				
				// Now that it has popped the node from the shared priority_queue, it will notify any waiting threads to access PQ now.
				// But the lock will be released only after the synchronized block ends its execution
				lock1.notifyAll();
			}
			
			if (this.isPath_taken && firstNodeinPQ.value.lowerbound_of > this.Current_minimumCost) {
				System.out.println("\nNode with lower bound " + firstNodeinPQ.value.lowerbound_of + " is pruned");
				continue;
			}
			
			// Gets the edge to be processed based on the edge matrix stored in this node
			int i=0; int j=0;
			boolean check = false;
			for (int v = 0; v < this.num_cities; v++) {
				if (!isVertexDone(firstNodeinPQ.value.constraintsMatrix, v)) {
					i = v;
					for (int col = i + 1; col < this.num_cities; col++) {
						if (firstNodeinPQ.value.constraintsMatrix[i][col] == 0) {
							j = col;
							check = true;
							break;
						}
					}
				}
				if (check)
					break;
			}
			
			// Check if the left child of this node is null. If it is then create the left child
			if (firstNodeinPQ.left == null) {		
				generateLeafNode(firstNodeinPQ, i, j, true);
			}
			
			// Check if the right child of this node is null. If yes, then create the right child
			if (firstNodeinPQ.right == null) {				
				generateLeafNode(firstNodeinPQ, i, j, false);				
			}
		}		
		/**********************   END of MULTI-THREADED CODE   *************************/
		
		
		/**********************   START of SEQUENTIAL CODE   *************************/
		
		// Checks if the queue still has some nodes to be processed
		/*while (!this.priority_queue.isEmpty()) {
			
			// Gets the i,j vertices whose edge is to be processed in the next
			// iteration
			int i=0; int j=0;
			boolean check = false;
			Node firstNodeinPQ = this.priority_queue.peek();
			for (int v = 0; v < this.num_cities; v++) {
				if (!isVertexDone(firstNodeinPQ.value.constraintsMatrix, v)) {
					i = v;
					for (int col = i + 1; col < this.num_cities; col++) {
						if (firstNodeinPQ.value.constraintsMatrix[i][col] == 0) {
							j = col;
							check = true;
							break;
						}
					}
				}
				if (check)
					break;
			}
			
			Node node = null;
			// Remove the top element from the queue
			synchronized (this.priority_queue) {
				node = this.priority_queue.remove();
			}			
			
			// Prints the node info which includes the constraintsMatrix showing edges included and excluded on calculation of the lower bound 
			// of this node
			//printNodeInfo(node);			
			
			if (this.isPath_taken && node.value.lowerbound_of > Current_minimumCost) {
				System.out.println("\nNode with lower bound " + node.value.lowerbound_of + " is pruned");
				continue;
			}
			// Check if the left child of this node is null. If it is then create the left child
			if (node.left == null) {				
				generateLeafNode(node, i, j, true);
			}
			
			// Check if the right child of this node is null. If yes, then create the right child
			if (node.right == null) {				
				generateLeafNode(node, i, j, false);				
			}
			
			System.out.println("PQueue size : "+this.priority_queue.size());
		}*/
		
		/**********************   END of SEQUENTIAL CODE   *************************/
	}
	
	private void generateLeafNode(Node parentNode, int i, int j, boolean isEdgeIncluded) {
		// Create the new edge matrix for the left child
		int[][] constraintsMatrix = new int[this.num_cities][this.num_cities];		

		// Keeping track of vertices whose edges are implicitly included/excluded due to the inclusion/exclusion of edge of other vertices
		Set<Integer> verticesToReprocess = new HashSet<Integer>();
		
		// Copies the edge matrix from the current node to the left child node
		for (int loop = 0; loop < this.num_cities; loop++) {
			constraintsMatrix[loop] = Arrays.copyOf(parentNode.value.constraintsMatrix[loop], this.num_cities);
		}

		// Considering the edge i,j and updating the state of the edge matrix
		if (isEdgeIncluded) constraintsMatrix[i][j] = constraintsMatrix[j][i] = 1;
		else constraintsMatrix[i][j] = constraintsMatrix[j][i] = -1;
		
		updateConstraintsMatrix(constraintsMatrix, verticesToReprocess);

		// Updating the edge matrix for the edges present in the verticesToReprocess set
		for (int v : verticesToReprocess) {
			updateConstraintsMatrix(constraintsMatrix, v);
		}

		// Clearing the set as all the vertices are now reprocessed
		verticesToReprocess.clear();

		// Check if the route is found at this node
		boolean Path_taken = isPath_taken(constraintsMatrix);

		double Path_Cost = 0;

		// Update the minimum cost if route is found else find the lower bound and update it in the node
		if (Path_taken) {
			Path_Cost = findPath_Cost(constraintsMatrix);
		} else {
			for (int v = 0; v < this.num_cities; v++) {
				Path_Cost += this.findSumOfTwoMinimumCostEdges(v, constraintsMatrix);
			}
			Path_Cost = Path_Cost / 2.0;
		}
		Node_Info value = new Node_Info(constraintsMatrix, Path_Cost);
		Node node = new Node(value, null, null, null, Path_taken);
		
		// Appends this node to the left or right of the parent node depending on whether an edge was considered or not
		if (isEdgeIncluded) parentNode.left = node;
		else parentNode.right = node;

		// Prints all the edges considered and not considered for calculation of the lower bound on this node
		if (isVerbose) {
			System.out.println("\n");
			System.out.println("Computing the "+ (isEdgeIncluded ? "left" : "right") + " node");
			printPathforthisNode(constraintsMatrix);
			display_Matrix(constraintsMatrix);
			if (Path_taken) System.out.println("Cost of the route found : "+Path_Cost);
			else System.out.println("Lower bound of this node : "+Path_Cost);
		}

		// If a route is found with the minimum cost, update the minimum cost and store this node
		
		if (Path_taken && Path_Cost < Current_minimumCost) {
			synchronized (lock2) {
				this.Current_minimumCost = (int) Path_Cost;
				this.minCostNode = node;
				this.isPath_taken = Path_taken;
			}
		} else if (Path_Cost < Current_minimumCost) {
			// Adds this node in the queue since no route is found till now and
			// so this node is a potential candidate for future processing.
			// However, this node should only be added if it has a lower bound <
			// cost of the minimum_cost route found till now
			synchronized (lock1) {
				this.priority_queue.add(node);
				/*count--;
				System.out.println("Decremented Count is : "+count);*/
				lock1.notifyAll();				
			}
		}
	}
	
	/**
	 * Checks if this vertex has been fully processed
	 * A vertex is fully processed until all the possible paths for 
	 * finding a route passing through that vertex are considered.
	 */
	private boolean isVertexDone(int[][] constraintsMatrix, int i) {

		int remainingEdges = 0;
		int includedEdges = 0;
		int not_includedEdges = 0;

		for (int j = 0; j < this.num_cities; j++) {
			if (i != j && constraintsMatrix[i][j] == 0) {
				remainingEdges++;
			} else if (i != j && constraintsMatrix[i][j] == -1) {
				not_includedEdges++;
			} else if (i != j && constraintsMatrix[i][j] == 1) {
				includedEdges++;
			}
		}
		if (remainingEdges == 0 && not_includedEdges == this.num_cities - 2 && includedEdges == 2) {
			return true;
		}
		return false;
	}

	// Finding the cost of the path taken.
	private int findPath_Cost(int[][] constraintsMatrix) {
		int cost = 0;
		for (int i = 0; i < this.num_cities; i++) {
			for (int j = i; j < this.num_cities; j++) {
				if (constraintsMatrix[i][j] == 1) {
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
	private boolean isPath_taken(int[][] constraintsMatrix) {

		for (int i = 0; i < this.num_cities; i++) {

			// Tracks the total number of edges that should be incident with
			// this vertex i
			int remainingEdges = 0;
			int includedEdges = 0;
			int not_includedEdges = 0;

			for (int j = 0; j < this.num_cities; j++) {
				if (i != j && constraintsMatrix[i][j] == 0) {
					remainingEdges++;
				} else if (i != j && constraintsMatrix[i][j] == -1) {
					not_includedEdges++;
				} else if (i != j && constraintsMatrix[i][j] == 1) {
					includedEdges++;
				}
			}
			
			if (remainingEdges != 0 || not_includedEdges != this.num_cities - 3 || includedEdges != 2) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Updates the constraints of all the edges incident with all the vertices in the graph.
	 * constraintsMatrix[i,j] == -1 => this edge is not included in the tour
	 * constraintsMatrix[i,j] == 1 => this edge must be included in the final tour
	 * Also takes care of the implicit exclusion inclusion of edges as a result of inclusion and exclusion of edge 
	 * under consideration

	 */
	private void updateConstraintsMatrix(int[][] constraintsMatrix, Set<Integer> s) {

		for (int i = 0; i < this.num_cities; i++) {
			// Tracks the total number of edges that should be incident with
			// this vertex i
			int includedEdges = 0;

			// Tracks the remaining edges to be considered
			int remainingEdges = 0;

			for (int j = 0; j < this.num_cities; j++) {
				if (i != j && constraintsMatrix[i][j] == 1) {
					includedEdges++;
				} else if (i != j && constraintsMatrix[i][j] == 0) {
					remainingEdges++;
				}
			}
			if (includedEdges < 2) {
				if ((remainingEdges == 2 && includedEdges == 0) || (remainingEdges == 1 && includedEdges == 1)) {
					for (int j = 0; j < this.num_cities; j++) {
						if (i != j && constraintsMatrix[i][j] == 0) {
							constraintsMatrix[i][j] = 1;
							constraintsMatrix[j][i] = 1;
							s.add(j);
						}
					}
				}
			} else {
				for (int j = 0; j < this.num_cities; j++) {
					if (i != j && constraintsMatrix[i][j] == 0) {
						constraintsMatrix[i][j] = -1;
						constraintsMatrix[j][i] = -1;
						s.add(j);
					}
				}
			}
		}
	}

// Updates the edge inclusion/exclusion constraints for all the edges incident with this vertex (v)

	private void updateConstraintsMatrix(int[][] constraintsMatrix, int v) {

		// Tracks the total number of edges that should be incident with
		// this vertex i
		int includedEdges = 0;

		// Tracks the remaining edges to be considered
		int remainingEdges = 0;

		for (int j = 0; j < this.num_cities; j++) {
			if (v != j && constraintsMatrix[v][j] == 1) {
				includedEdges++;
			} else if (v != j && constraintsMatrix[v][j] == 0) {
				remainingEdges++;
			}
		}
		if (includedEdges < 2) {
			if ((remainingEdges == 2 && includedEdges == 0) || (remainingEdges == 1 && includedEdges == 1)) {
				for (int j = 0; j < this.num_cities; j++) {
					if (v != j && constraintsMatrix[v][j] == 0) {
						constraintsMatrix[v][j] = 1;
						constraintsMatrix[j][v] = 1;
					}
				}
			}
		} else {
			for (int j = 0; j < this.num_cities; j++) {
				if (v != j && constraintsMatrix[v][j] == 0) {
					constraintsMatrix[v][j] = -1;
					constraintsMatrix[j][v] = -1;
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
	private int findSumOfTwoMinimumCostEdges(int v, int[][] constraintsMatrix) {

		int firstMinimumCost = Integer.MAX_VALUE;
		int secondMinimumCost = Integer.MAX_VALUE;

		boolean isFirstMinFound = false;
		boolean isSecondMinFound = false;

		for (int i = 0; i < this.num_cities; i++) {
			if (i != v && adjacency_Matrix[v][i] != -1 && constraintsMatrix[v][i] != -1) {
				if (isFirstMinFound && isSecondMinFound) {
					return firstMinimumCost + secondMinimumCost;
				}
				if (!isFirstMinFound && constraintsMatrix[v][i] == 1) {
					firstMinimumCost = adjacency_Matrix[v][i];
					isFirstMinFound = true;
					continue;
				} else if (!isSecondMinFound && constraintsMatrix[v][i] == 1) {
					secondMinimumCost = adjacency_Matrix[v][i];
					isSecondMinFound = true;
					continue;
				}

				if (!isFirstMinFound && adjacency_Matrix[v][i] < firstMinimumCost) {
					secondMinimumCost = firstMinimumCost;
					firstMinimumCost = adjacency_Matrix[v][i];
				} else if (!isSecondMinFound && adjacency_Matrix[v][i] < secondMinimumCost) {
					secondMinimumCost = adjacency_Matrix[v][i];
				}
			}
		}
		return firstMinimumCost + secondMinimumCost;
	}
	
// Prints all the edges that are considered along with the ones that are not considered for calculation of lower bound on this node

	private void printPathforthisNode(int[][] constraintsMatrix) {
		StringBuilder includedEdges = new StringBuilder();
		StringBuilder not_includedEdges = new StringBuilder();
		for (int i = 0; i < this.num_cities; i++) {
			for (int j = i; j < this.num_cities; j++) {
				if (constraintsMatrix[i][j] == 1) {
					includedEdges.append(" (" + this.cities.get(i) + "-" + this.cities.get(j) + ")");
				} else if (constraintsMatrix[i][j] == -1) {
					not_includedEdges.append(" (" + this.cities.get(i) + "-" + this.cities.get(j) + ")");
				}
			}
		}
		System.out.println("Edges included :" + includedEdges.toString());
		System.out.println("Edges Not included :" + not_includedEdges.toString());
	}

	 // Displays the adjacency matrix representing the graph
	private void display_Matrix(int[][] matrix) {
		
		System.out.println("Printing the adjacency matrix representing the graph.\n");
		
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
	
	@SuppressWarnings("unused")
	private void printNodeInfo(Node node) {
		System.out.println("\nBelow node is now popped from the priority_queue");
		display_Matrix(node.value.constraintsMatrix);
		printPathforthisNode(node.value.constraintsMatrix);
		System.out.println("Lower bound of this node is : "+node.value.lowerbound_of);
		System.out.println("Does this node has a left child? : "+(node.left != null));
		System.out.println("Does this node has a right child? : "+(node.right != null));
		System.out.println();
	}
	
	// Prints the minimum cost route on the console along with its cost
	
	private void PrintOptimalRoute(Node node) {
		System.out.println("\n\nThe optimal path will include the below edges of the graph");
		List<String> edgesInvolvedInOptimalTour = new ArrayList<String>();
		for (int i=0; i < node.value.constraintsMatrix[0].length; i++) {
			for (int j=i+1; j < node.value.constraintsMatrix[i].length; j++) {
				if (node.value.constraintsMatrix[i][j] == 1) {
					edgesInvolvedInOptimalTour.add(this.cities.get(i) + this.cities.get(j));
					System.out.print(" "+this.cities.get(i) + this.cities.get(j));
				}
			}
		}
		display_Matrix(node.value.constraintsMatrix);
		System.out.println("\nOptimal Route Cost is : "+node.value.lowerbound_of);
	}
	
	@SuppressWarnings("resource")
	private void read_inputfile() {
		Scanner scan = new Scanner(System.in);
		System.out.println("Enter the number of cities");
		this.num_cities = scan.nextInt();

		// Initialize the cities array size
		this.cities = new ArrayList<String>();
		
		// Add cities in the cities array.. starting from A
		for (int b = 0; b < this.num_cities; b++) {
			this.cities.add((char)(b+65) + "");
		}
		
		this.adjacency_Matrix = new int[this.num_cities][this.num_cities];
		
		boolean isValidPath = false;
		String graph = "";
		
		while (!isValidPath) {
			System.out.println("Enter the path to the input text file containing the graph represented as an adjacency matrix");
			graph = scan.next().trim();
			if (graph.contains(".txt")) {
				isValidPath = true;
			} else {
				System.out.println("Invalid path");
			}
		}

		try {
			Scanner readGraph = new Scanner(new File(graph));
			int i=0;
			int j=0;
			while (readGraph.hasNextLine()) {
				String row = readGraph.nextLine();
				String[] cols = row.split(" ");
				
				for (String col : cols) {
					int cost = Integer.parseInt(col.trim());
					this.adjacency_Matrix[i][j] = cost;
					j++;
				}
				i++;
				j = 0;
			}
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}	
		
		for (String name : this.cities) {
			this.longestCityName = Math.max(this.longestCityName, name.length());
		}
		
		System.out.println("Enter the total number of threads that you want ? ");
		this.Total_Threads = Integer.parseInt(scan.next());
		//this.Total_Threads = 2;
		
		System.out.println("Do you want to switch the Verbose mode on (Type Y for yes and N for No) : ");
		isVerbose = scan.next().equalsIgnoreCase("y") ? true : false;
		
	}

	private void printSpaces(int s) {
		for (int i = 0; i < s; i++) {
			System.out.print(" ");
		}
	}
}

 // Stores the constraintsMatrix and the lower bound calculated at each of the node

class Node_Info {
	int[][] constraintsMatrix;
	double lowerbound_of;

	public Node_Info(int[][] constraintsMatrix, double cost) {
		super();
		this.constraintsMatrix = constraintsMatrix;
		this.lowerbound_of = cost;
	}
}

// Represents the current node and its left and right children
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
		if (this.value.lowerbound_of < node.value.lowerbound_of) {
			return -1;
		} else if (this.value.lowerbound_of > node.value.lowerbound_of) {
			return 1;
		} else {
			return 0;
		}
	}
}
