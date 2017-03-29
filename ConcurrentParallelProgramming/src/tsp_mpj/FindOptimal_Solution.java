/**
 * 
 */
package tsp_mpj;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;

import mpi.MPI;

/**
 * @author Manav
 * @author karthik
 *
 */
public class FindOptimal_Solution {
	
	private static Random random;
	
	// Flag to check if verbose mode is on
	private boolean isVerbose = true;
	
	// List to store all the cities entered by a user
	private List<String> cities;
	
	// Adjacency matrix to store the weight/cost of the edges between each pair of vertices
	private int[][] adjacency_Matrix;
	
	// cities entered by a user are stored here
	 private int num_Cities;

	private int longestCity_Name = 0;

	// Checks if a route containing all the nodes is found or not.
	private boolean isPath_Taken = false;
	
	// Stores the cost of the minimum cost route found
	private int minimumCost = Integer.MAX_VALUE;

	// Stores the node that contains the route having the minimum cost
	private Node minCostNode;
	
	// Keeps a count of the number of processors
	private int num_Processors;
	
	// Keeps track of the rank of this processor
	private int My_Rank;

	// Nodes part of optimal solution are stored here.
	// Solution with lowest lower bound gets higher priority
	private PriorityQueue<Node> priority_Queue = new PriorityQueue<Node>();
	
	private long startTime;
	private long endTime;
	
	private Node[] node_Buffer = new Node[1];
	
	
	/**
	 * Start the main program
	 * @param args
	 */
	public static void main(String[] args) {
		
		FindOptimal_Solution Solution = new FindOptimal_Solution();
		
		// Initialize the OpenMP framework
		MPI.Init(args);
		
		Solution.My_Rank = MPI.COMM_WORLD.Rank();
		Solution.num_Processors = MPI.COMM_WORLD.Size();
		
		/*// Initialize the cities array size
		Solution.cities = new ArrayList<String>();
		// Add cities in the cities array.. starting from A
		for (int b = 0; b < Solution.num_Cities; b++) {
			Solution.cities.add((char) (b + 65) + "");
		}*/
		
		// Process 0 is the master process
		// It generates till the number of nodes in the priority queue is equal to the number of processes available to process them
		if (Solution.My_Rank == 0) {
			
			// Generate a random matrix
			// Solution.generate_Matrix();
			
			// Reads the input from a file
			Solution.read_inputFile();
			
			Solution.startTime = System.currentTimeMillis();
			System.out.println("Start time is : "+Solution.startTime);
			
			// An Constraints_Matrix to keep track of the included and excluded edges
			// for each node. 
			int[][] Constraints_Matrix = new int[Solution.num_Cities][Solution.num_Cities];
			
			// Initializing edge matrix to have -1 between those vertices for which no edge exists
			/*for (int i=0; i<Solution.num_Cities; i++) {
				for (int j=0; j<Solution.num_Cities; j++) {
					if (Solution.adjacency_Matrix[i][j] == -1) Constraints_Matrix[i][j] = -1;
				}
			}*/

			// Display the adjacency matrix
			if (Solution.isVerbose) Solution.display_Matrix(Solution.adjacency_Matrix);		
			
			// Finding the lower bound at the root node
			double cost = 0.0;
			for (int v = 0; v < Solution.num_Cities; v++) {
				cost += Solution.findSumOfTwoMinimumCostEdges(v, Constraints_Matrix);
			}
			
			// Creating the root node that considers all the possible routes
			Node_Info value = new Node_Info(Constraints_Matrix, cost/2);
			Node root = new Node(value, null, null, null, false);

			System.out.println("Estimated lower bound cost at root node is : " + root.value.lowerBound_Of);

			// Adding the root node to the queue
			Solution.priority_Queue.add(root);
			
			// Generates the child nodes and puts them in the Priority Queue until the size of the queue is equal to the number of processes
			while (Solution.priority_Queue.size() < Solution.num_Processors - 1) {
				Solution.FindOptimal_Solution();
			}
			
			// Process 0 has generated nodes = number of slave processes available
			// Now send these nodes one by one to the processes
			
			for (int i=1; i < Solution.num_Processors; i++) {
				
				Node n = Solution.priority_Queue.poll();
				
				try {
					// Converting the node to byte array
					ByteArrayOutputStream bos = new ByteArrayOutputStream();
					ObjectOutput out = new ObjectOutputStream(bos);
					out.writeObject(n);
					byte[] node_Bytes = bos.toByteArray();
					
					// Converting the adjacency matrix to byte array
					bos = new ByteArrayOutputStream();
					out = new ObjectOutputStream(bos);
					out.writeObject(Solution.adjacency_Matrix);
					byte[] adjacency_MatrixBytes = bos.toByteArray();
					
					// MPI send the adjacency matrix to each of the slave processors
					MPI.COMM_WORLD.Isend(adjacency_MatrixBytes, 0, adjacency_MatrixBytes.length, MPI.BYTE, i, i);
					
					// MPI send command to send the node n to process i
					MPI.COMM_WORLD.Isend(node_Bytes, 0, node_Bytes.length, MPI.BYTE, i, i);
					
				} catch (IOException ex) {
					ex.printStackTrace();
				}
			}
			
			System.out.println("Master process has send the nodes in the queue to all the processors");
			
			// Master should now wait until all the slaves have finished processing and return their respective results to the master
			for (int i=1; i < Solution.num_Processors; i++) {
				Node minCostNode_Received = null;
				byte minCost_Nodes[] = new byte[5000];
				
				MPI.COMM_WORLD.Recv(minCost_Nodes, 0, 5000, MPI.BYTE, i, i);
				
				ByteArrayInputStream bis = new ByteArrayInputStream(minCost_Nodes);
				ObjectInput in = null;
				try {
					in = new ObjectInputStream(bis);
					Object obj = in.readObject();
					minCostNode_Received = (Node) obj;
					
					if (Solution.minimumCost > minCostNode_Received.value.lowerBound_Of) {
						Solution.minimumCost = (int) minCostNode_Received.value.lowerBound_Of;
						Solution.minCostNode = minCostNode_Received;
					}
					
				} catch (IOException ex) {
					ex.printStackTrace();
				} catch (ClassNotFoundException cnf) {
					cnf.printStackTrace();
				}				
			}
					
			Solution.printRouteTaken(Solution.minCostNode);
			
			Solution.endTime = System.currentTimeMillis();
			System.out.println("End time is : "+Solution.endTime);
			
			System.out.println("Total time taken is : " + (Solution.endTime - Solution.startTime));
		} else {
			
			// Other processors tend to keep listening to the parent processor if it sends any new node
			// Upon receiving the node, they will calculate the lower bound of the node. 
			// Once the lower bound of the node is calculated, they must check that if the lower bound of the new node is greater than the min cost found
			// If yes, this node must not be sent back to the parent process for addition to PQ
			// If not, then send this to parent process so that it can be added to PQ
			
			// PriorityQueue<Node> priority_Queue = new PriorityQueue<Node>();
			
			byte[] node_Bytes = new byte[5000];
			byte[] adjacency_MatrixBytes = new byte[5000];
			Node node_Received = null;
			
			// Each process should receive a node from the parent process
			// After receiving the node, it casts the buffer contents back to a Node class
			// It then adds the node to its own local Priority Queue
			// Then it starts processing the node i.e. finding the subtree rooted at this node.
			for (int i=1; i < Solution.num_Processors; i++) {
				if (Solution.My_Rank == i) {
					
					//System.out.println("Before Receive : "+Solution.My_Rank);
					
					// Receives the adjacency matrix from the parent process
					MPI.COMM_WORLD.Recv(adjacency_MatrixBytes, 0, 5000, MPI.BYTE, 0, i);
					
					// Receives the node from the parent processor
					MPI.COMM_WORLD.Recv(node_Bytes, 0, 5000, MPI.BYTE, 0, i);
					
					//System.out.println("After Receive : "+Solution.My_Rank);
					
					ByteArrayInputStream bis1 = new ByteArrayInputStream(adjacency_MatrixBytes);
					ObjectInput in1 = null;
					try {
						in1 = new ObjectInputStream(bis1);
						Object obj = in1.readObject();
						Solution.adjacency_Matrix = (int[][]) obj;
					} catch (IOException ex) {
						ex.printStackTrace();
					} catch (ClassNotFoundException cnf) {
						cnf.printStackTrace();
					}
					
					ByteArrayInputStream bis2 = new ByteArrayInputStream(node_Bytes);
					ObjectInput in2 = null;
					try {
						in2 = new ObjectInputStream(bis2);
						Object obj = in2.readObject();
						node_Received = (Node) obj;
					} catch (IOException ex) {
						ex.printStackTrace();
					} catch (ClassNotFoundException cnf) {
						cnf.printStackTrace();
					}
					
					Solution.num_Cities = node_Received.value.Constraints_Matrix[0].length;
					// Initialize the cities array size
					Solution.cities = new ArrayList<String>();
					
					// Add cities in the cities array.. starting from A
					for (int b = 0; b < Solution.num_Cities; b++) {
						Solution.cities.add((char)(b+65) + "");
					}
					Solution.priority_Queue.add(node_Received);
				}
			}
			
			
			// By now the node has been received, so this process will start processing the node
			// The process will calculate the subtree rooted at this node and then return the minimum cost of the route found in this branch of
			// the main tree
			
			while (!Solution.priority_Queue.isEmpty()) {
				Solution.FindOptimal_Solution();
			}			
			try {
				ByteArrayOutputStream bos = new ByteArrayOutputStream();
				ObjectOutput out = null;
				out = new ObjectOutputStream(bos);
				out.writeObject(Solution.minCostNode);
				byte[] minCostNode_Bytes = bos.toByteArray();
				
				// MPI send command to send the node n to process i
				MPI.COMM_WORLD.Isend(minCostNode_Bytes, 0, minCostNode_Bytes.length, MPI.BYTE, 0, Solution.My_Rank);
				
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		// MPI.COMM_WORLD.Barrier();
		MPI.Finalize();
	}
	
	private void generate_Matrix() {		
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


	private void FindOptimal_Solution()  {
		
/**********************   START of SEQUENTIAL CODE   *************************/
		
		// Checks if the queue still has some nodes to be processed
		//while (!this.priority_Queue.isEmpty()) {
			
			// Gets the i,j vertices whose edge is to be processed in the next
			// iteration
			int i=0; int j=0;
			boolean check = false;
			Node firstNodeinPQ = this.priority_Queue.peek();
			
			// System.out.println("Process " + this.My_Rank + " is going to print its edge matrix");
			// this.display_Matrix(firstNodeinPQ.value.Constraints_Matrix);
			
			for (int v = 0; v < this.num_Cities; v++) {
				if (!isVertexDone(firstNodeinPQ.value.Constraints_Matrix, v)) {
					i = v;
					for (int col = i + 1; col < this.num_Cities; col++) {
						if (firstNodeinPQ.value.Constraints_Matrix[i][col] == 0) {
							j = col;
							check = true;
							break;
						}
					}
				}
				if (check)
					break;
			}
			
			Node node = this.priority_Queue.remove();		
			
			// Prints the node info which includes the Constraints_Matrix showing edges included and excluded on calculation of the lower bound 
			// of this node
			//printNodeInfo(node);			
			
			if (this.isPath_Taken && node.value.lowerBound_Of > minimumCost) {
				System.out.println("\nNode with lower bound " + node.value.lowerBound_Of + " is pruned");
				return;
			}
			// Check if the left child of this node is null. If it is then create the left child
			if (node.left == null) {				
				generateLeafNode(node, i, j, true);
			}
			
			// Check if the right child of this node is null. If yes, then create the right child
			if (node.right == null) {				
				generateLeafNode(node, i, j, false);				
			}
			
			//System.out.println("priority_Queue size : "+this.priority_Queue.size());
		//}
		
		/**********************   END of SEQUENTIAL CODE   *************************/
		
		
		
	}
	
	private void generateLeafNode(Node parentNode, int i, int j, boolean isEdgeIncluded) {
		// Create the new edge matrix for the left child
		int[][] Constraints_Matrix = new int[this.num_Cities][this.num_Cities];		

		// Keeping track of vertices whose edges are implicitly included/excluded due to the inclusion/exclusion of edge of other vertices
		Set<Integer> constraints_reload = new HashSet<Integer>();
		
		// Copies the edge matrix from the current node to the left child node
		for (int loop = 0; loop < this.num_Cities; loop++) {
			Constraints_Matrix[loop] = Arrays.copyOf(parentNode.value.Constraints_Matrix[loop], this.num_Cities);
		}
		
		System.out.println("i = " + i);
		System.out.println("j = " + j);

		// Considering the edge i,j and updating the state of the edge matrix
		if (isEdgeIncluded) Constraints_Matrix[i][j] = Constraints_Matrix[j][i] = 1;
		else Constraints_Matrix[i][j] = Constraints_Matrix[j][i] = -1;
		
		boolean shouldBePruned =  updateConstraints_Matrix(Constraints_Matrix, constraints_reload);
		
		System.out.println("ShoudBePruned : " + shouldBePruned);
		if (shouldBePruned) return;

		// Updating the edge matrix for the edges present in the constraints_reload set
		for (int v : constraints_reload) {
			updateConstraints_Matrix(Constraints_Matrix, v);
		}

		// Clearing the set as all the vertices are now reprocessed
		constraints_reload.clear();

		// Check if the route is found at this node
		boolean Path_Taken = isPath_Taken(Constraints_Matrix);

		double Path_Cost = 0;

		// Update the minimum cost if route is found else find the lower bound and update it in the node
		if (Path_Taken) {
			Path_Cost = findPath_Cost(Constraints_Matrix);
		} else {
			for (int v = 0; v < this.num_Cities; v++) {
				Path_Cost += this.findSumOfTwoMinimumCostEdges(v, Constraints_Matrix);
			}
			Path_Cost = Path_Cost / 2.0;
		}
		Node_Info value = new Node_Info(Constraints_Matrix, Path_Cost);
		Node node = new Node(value, null, null, null, Path_Taken);
		
		// Appends this node to the left or right of the parent node depending on whether an edge was considered or not
		if (isEdgeIncluded) parentNode.left = node;
		else parentNode.right = node;

		// Prints all the edges considered and not considered for calculation of the lower bound on this node
		if (isVerbose) {
			System.out.println("\n");
			System.out.println("Computing the "+ (isEdgeIncluded ? "left" : "right") + " node");
			printPathforthisNode(Constraints_Matrix);
			display_Matrix(Constraints_Matrix);
			if (Path_Taken) System.out.println("Cost of the route found : "+Path_Cost);
			else System.out.println("Lower bound of this node : "+Path_Cost);
		}

		// If a route is found with the minimum cost, update the minimum cost and store this node
		
		if (Path_Taken && Path_Cost < this.minimumCost) {
			this.minimumCost = (int) Path_Cost;
			this.minCostNode = node;
			this.isPath_Taken = Path_Taken;
			
			// Broadcast the minimum cost found till now to all the processing nodes
			/*for (int rank=1; rank < this.num_Processors; rank++) {
				if (rank == this.My_Rank) {
					continue;
				}
				MPI.COMM_WORLD.Isend(new int[]{this.minimumCost}, 0, 1, MPI.INT, rank, 0);
			}	*/	
		} else if (Path_Cost < minimumCost) {
			// Adds this node in the queue since no route is found till now and
			// so this node is a potential candidate for
			// future processing.
			// However, this node should only be added if it has a lower bound <
			// cost of the mincost route found till now
			this.priority_Queue.add(node);
		}
	}
	
	/**
	 * Checks if this vertex has been fully processed
	 * A vertex is fully processed until all the possible paths for 
	 * finding a route passing through that vertex are considered.
	 */
	private boolean isVertexDone(int[][] Constraints_Matrix, int i) {

		int remainingEdges = 0;
		int includedEdges = 0;
		int not_includedEdges = 0;

		for (int j = 0; j < this.num_Cities; j++) {
			if (i != j && Constraints_Matrix[i][j] == 0) {
				remainingEdges++;
			} else if (i != j && Constraints_Matrix[i][j] == -1) {
				not_includedEdges++;
			} else if (i != j && Constraints_Matrix[i][j] == 1) {
				includedEdges++;
			}
		}
		if (remainingEdges == 0 && not_includedEdges == this.num_Cities - 2 && includedEdges == 2) {
			return true;
		}
		return false;
	}

	//Finds the cost of the route found
	 
	private int findPath_Cost(int[][] Constraints_Matrix) {
		int cost = 0;
		for (int i = 0; i < this.num_Cities; i++) {
			for (int j = i; j < this.num_Cities; j++) {
				if (Constraints_Matrix[i][j] == 1) {
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
	private boolean isPath_Taken(int[][] Constraints_Matrix) {

		for (int i = 0; i < this.num_Cities; i++) {

			// Tracks the total number of edges that should be incident with
			// this vertex i
			int remainingEdges = 0;
			int includedEdges = 0;
			int not_includedEdges = 0;

			for (int j = 0; j < this.num_Cities; j++) {
				if (i != j && Constraints_Matrix[i][j] == 0) {
					remainingEdges++;
				} else if (i != j && Constraints_Matrix[i][j] == -1) {
					not_includedEdges++;
				} else if (i != j && Constraints_Matrix[i][j] == 1) {
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
	 * Updates the state of all the edges incident with all the vertices in the graph.
	 * Constraints_Matrix[i,j] == -1 => this edge is not considered in the tour
	 * Constraints_Matrix[i,j] == 1 => this edge must be considered in the final tour
	 * Also takes care of the implicit exclusion inclusion of edges as a result of inclusion and exclusion of edge 
	 * under consideration
	 * 
	 * @param Constraints_Matrix
	 * @param s
	 */
	private boolean updateConstraints_Matrix(int[][] Constraints_Matrix, Set<Integer> s) {

		for (int i = 0; i < this.num_Cities; i++) {
			// Tracks the total number of edges that should be incident with
			// this vertex i
			int includedEdges = 0;

			// Tracks the remaining edges to be considered
			int remainingEdges = 0;

			for (int j = 0; j < this.num_Cities; j++) {
				if (i != j && Constraints_Matrix[i][j] == 1) {
					includedEdges++;
				} else if (i != j && Constraints_Matrix[i][j] == 0) {
					remainingEdges++;
				}
			}
			
			if (includedEdges > 2 || (includedEdges + remainingEdges) < 2) {
				return true;
			}
			
			if (includedEdges < 2) {
				if ((remainingEdges == 2 && includedEdges == 0) || (remainingEdges == 1 && includedEdges == 1)) {
					for (int j = 0; j < this.num_Cities; j++) {
						if (i != j && Constraints_Matrix[i][j] == 0) {
							Constraints_Matrix[i][j] = 1;
							Constraints_Matrix[j][i] = 1;
							s.add(j);
						}
					}
				}
			} else {
				for (int j = 0; j < this.num_Cities; j++) {
					if (i != j && Constraints_Matrix[i][j] == 0) {
						Constraints_Matrix[i][j] = -1;
						Constraints_Matrix[j][i] = -1;
						s.add(j);
					}
				}
			}
		}
		return false;
	}

	/**
	 * Updates the edge inclusion/exclusion contraints for all the edges incident with this vertex (v)

	 */
	private void updateConstraints_Matrix(int[][] Constraints_Matrix, int v) {

		// Tracks the total number of edges that should be incident with
		// this vertex i
		int includedEdges = 0;

		// Tracks the remaining edges to be considered
		int remainingEdges = 0;

		for (int j = 0; j < this.num_Cities; j++) {
			if (v != j && Constraints_Matrix[v][j] == 1) {
				includedEdges++;
			} else if (v != j && Constraints_Matrix[v][j] == 0) {
				remainingEdges++;
			}
		}
		if (includedEdges < 2) {
			if ((remainingEdges == 2 && includedEdges == 0) || (remainingEdges == 1 && includedEdges == 1)) {
				for (int j = 0; j < this.num_Cities; j++) {
					if (v != j && Constraints_Matrix[v][j] == 0) {
						Constraints_Matrix[v][j] = 1;
						Constraints_Matrix[j][v] = 1;
					}
				}
			}
		} else {
			for (int j = 0; j < this.num_Cities; j++) {
				if (v != j && Constraints_Matrix[v][j] == 0) {
					Constraints_Matrix[v][j] = -1;
					Constraints_Matrix[j][v] = -1;
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
	private int findSumOfTwoMinimumCostEdges(int v, int[][] Constraints_Matrix) {

		int firstMinimumCost = Integer.MAX_VALUE;
		int secondMinimumCost = Integer.MAX_VALUE;

		boolean isFirstMinFound = false;
		boolean isSecondMinFound = false;

		for (int i = 0; i < this.num_Cities; i++) {
			if (i != v && adjacency_Matrix[v][i] != -1 && Constraints_Matrix[v][i] != -1) {
				if (isFirstMinFound && isSecondMinFound) {
					return firstMinimumCost + secondMinimumCost;
				}
				if (!isFirstMinFound && Constraints_Matrix[v][i] == 1) {
					firstMinimumCost = adjacency_Matrix[v][i];
					isFirstMinFound = true;
					continue;
				} else if (!isSecondMinFound && Constraints_Matrix[v][i] == 1) {
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
	

	 // Prints all the edges that are considered/not considered for calculation of lower bound on this node
	 
	private void printPathforthisNode(int[][] Constraints_Matrix) {
		StringBuilder includedEdges = new StringBuilder();
		StringBuilder not_includedEdges = new StringBuilder();
		for (int i = 0; i < this.num_Cities; i++) {
			for (int j = i; j < this.num_Cities; j++) {
				if (Constraints_Matrix[i][j] == 1) {
					includedEdges.append(" (" + this.cities.get(i) + "-" + this.cities.get(j) + ")");
				} else if (Constraints_Matrix[i][j] == -1) {
					not_includedEdges.append(" (" + this.cities.get(i) + "-" + this.cities.get(j) + ")");
				}
			}
		}
		System.out.println("Edges included :" + includedEdges.toString());
		System.out.println("Edges Not included :" + not_includedEdges.toString());
	}

	
	 // Displays the adjacency matrix representing the graph
	 
	private void display_Matrix(int[][] matrix) {
		
		System.out.println("Printing the edge matrix representing the graph.\n");
		
		printSpaces(this.longestCity_Name + 2);
		for (int i = 0; i < this.cities.size(); i++) {
			System.out.print(this.cities.get(i) + " ");
		}

		System.out.println();

		for (int i = 0; i < this.cities.size(); i++) {
			System.out.print(this.cities.get(i));
			printSpaces(this.longestCity_Name - cities.get(i).length() + 2);
			
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
		System.out.println("\nBelow node is now popped from the priority_Queue");
		display_Matrix(node.value.Constraints_Matrix);
		printPathforthisNode(node.value.Constraints_Matrix);
		System.out.println("Lower bound of this node is : "+node.value.lowerBound_Of);
		System.out.println("Does this node has a left child? : "+(node.left != null));
		System.out.println("Does this node has a right child? : "+(node.right != null));
		System.out.println();
	}
	
	/**
	 * Prints the minimum cost route on the console along with its cost
	 */
	private void printRouteTaken(Node node) {
		System.out.println("\n\nThe optimal path will include the below edges of the graph");
		List<String> edgesInvolvedInOptimalTour = new ArrayList<String>();
		for (int i=0; i < node.value.Constraints_Matrix[0].length; i++) {
			for (int j=i+1; j < node.value.Constraints_Matrix[i].length; j++) {
				if (node.value.Constraints_Matrix[i][j] == 1) {
					edgesInvolvedInOptimalTour.add(this.cities.get(i) + this.cities.get(j));
					System.out.print(" "+this.cities.get(i) + this.cities.get(j));
				}
			}
		}
		display_Matrix(node.value.Constraints_Matrix);
		System.out.println("\nOptimal Route Cost is : "+node.value.lowerBound_Of);
	}
	
	@SuppressWarnings("resource")
	private void read_inputFile() {
		/*Scanner scan = new Scanner(System.in);
		System.out.println("Enter the number of cities");
		this.num_Cities = scan.nextInt();*/
		
		String graph = "tsp_graph_02.txt";
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
			this.longestCity_Name = Math.max(this.longestCity_Name, name.length());
		}

		
	}

	private void printSpaces(int s) {
		for (int i = 0; i < s; i++) {
			System.out.print(" ");
		}
	}
}

// Stores the Constraints_Matrix and the lower bound calculated at each of the node
 
class Node_Info implements Serializable {
	int[][] Constraints_Matrix;
	double lowerBound_Of;

	public Node_Info(int[][] Constraints_Matrix, double cost) {
		super();
		this.Constraints_Matrix = Constraints_Matrix;
		this.lowerBound_Of = cost;
	}
}

// Represents the current node and its left and right children
 
class Node implements Comparable<Node>, Serializable {
	Node_Info value;
	Node left;
	Node right;
	Node parent;
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
		if (this.value.lowerBound_Of < node.value.lowerBound_Of) {
			return -1;
		} else if (this.value.lowerBound_Of > node.value.lowerBound_Of) {
			return 1;
		} else {
			return 0;
		}
	}
}
