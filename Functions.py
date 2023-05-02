import networkx as nx
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, odeint
import numpy as np
from matplotlib import animation
import random as rand
import math



class Node:
    def __init__(self, name, conc_u, conc_v, pos, o, lat_pos=0, dir=math.pi/4):
        self.name = name
        self.conc_u = conc_u
        self.conc_v = conc_v
        self.pos = pos
        self.o = o
        self.lat_pos = lat_pos
        self.dir = dir

    def __str__(self):
        """print function"""
        return f"{self.name} - u:{self.conc_u} v:{self.conc_v} {self.pos}"
    

    def update_node(self, new_conc_u, new_conc_v):
        """Update concentrations of nodes"""
        self.conc_u = new_conc_u
        self.conc_v = new_conc_v

    def update_pos(self, new_pos_x, new_pos_y):
        """update positions of nodes"""
        self.pos = (new_pos_x, new_pos_y)

    def update_number(self, new_num):
        """update node ordering"""
        self.name = new_num


class Modular_network:
    
    def initialise(self, m=20, l=20):
        """Initialise network of m modules each containing l nodes"""
        n = m*l
        

        # initialise adjacency matrix
        A = np.zeros((m*l, m*l))

        # prob_in = 1
        # prob_out = 1/64
        # for a in range(m):
        #     for b in range(l):
        #         for c in range(n - ((a+1)*l)):  # connect nodes to other modules
        #             if np.random.choice([0,1], 1, p=[1-prob_out, prob_out])==1 and a != m-1:
        #                     A[a*l + b][(a+1)*l + c] = 1
        #         for d in range(l - b):  # connect nodes in a module to each other
        #             if np.random.choice([0,1], 1, p=[1-prob_in, prob_in]) == 1:
        #                 A[a*l + b][a*l + b + d] = 1

        for a in range(m):
            for b in range(m):
                A[a*l + rand.randint(0,l-1)][b*l + rand.randint(0,l-1)] = 1  # randomly connect one node in each module to another
                for d in range(l):
                    for e in range(l):  
                        A[b*l + d][b*l + e] = 1  # connect nodes in a module to each other
                        


        A = np.tril(A) + np.triu(A.T, 1)  # make symmetrical
        #A = A + A.T - np.diag(np.diag(A))
        np.fill_diagonal(A, 0)  # nodes not connected to themselves
        G = nx.from_numpy_array(A)  # create initial graph G
        pos = nx.spring_layout(G)  # get positions of nodes from initial graph

        np.random.seed(1)
        #r = (np.random.random(N) - 0.5)/10
        r = np.random.normal(0, 0.05, n)
        np.random.seed(2)
        #r2 = (np.random.random(N) - 0.5)/10
        r2 = np.random.normal(0, 0.05, n)

        nodes = {}
        for i in range(m):
            for j in range(l):
                nodes[i*l + j] = Node(i*l + j, 1+r[i*l + j], 1+r2[i*l + j], (pos[i*l + j][0], pos[i*l + j][1]), i)

        self.Adj = A
        self.Nodes = nodes
        self.N = n

        return A, nodes, n
    
    def centroid(self, vertexes):
        """find the center of a polygon with vertexes[]"""
        _x_list = [vertex[0] for vertex in vertexes]
        _y_list = [vertex[1] for vertex in vertexes]
        _len = len(vertexes)
        _x = sum(_x_list) / _len
        _y = sum(_y_list) / _len
        # return center +- a bit
        return(_x + np.random.normal(0.1, 0.1), _y + np.random.normal(0.1, 0.1))
    
    def add_node(self, node):
        """grow the network with a daughter node """
        Adj = self.Adj
        N = self.N
        Nodes = self.Nodes

        prob_in = 8/N
        prob_out = 8/N

        Adj = np.append(Adj, np.zeros([N,1], dtype=int), axis=1)  # add column of zeros
        Adj = np.append(Adj, np.zeros([1,N+1], dtype=int) , axis=0)  # add row of zeros
        Adj[N][node.name] = 1  # connected to parent node
        Adj[node.name][N]
        for n in Nodes:
            if Nodes[n].o == node.o:
                if np.random.choice([0,1], 1, p=[1-prob_in, prob_in])==1:
                    Adj[N][n] = 1
                    Adj[n][N] = 1
                    print("connected in")
            else:
                if np.random.choice([0,1], 1, p=[1-prob_out, prob_out])==1:
                    Adj[N][n] = 1
                    Adj[n][N] = 1
                    print("connected out")
        # connect with probability 4/parent connections)
        # for node in range(len(Adj)): ##############################################################################
        

        neighbours = []
        for i in range(N):
            if Adj[N, i] == 1:
                neighbours.append(Nodes[i].pos)
        print(neighbours)
        c = self.centroid(neighbours)
        node.update_node(node.conc_u/2, node.conc_v/2)
        Nodes[N] = Node(N, node.conc_u, node.conc_v, c, node.o)

        self.N = N + 1
        self.Adj = Adj
        self.Nodes = Nodes

        return Adj
    
    def delete_node(self, node):
        """Remove all connections to and from node but keep entry in dictionary"""
        Adj = self.Adj
        N = self.N
        Nodes = self.Nodes

        # delete row and column from Adjacency matrix
        Adj[0:N][node] = 0
        Adj[node][0:N] = 0
        self.Adj = Adj

        Nodes[node].update_node(0,0)
        self.Nodes = Nodes


class Lattice:

    def initialise(self, N=5, M=5, rigid=True, denser=True):

        G = nx.grid_2d_graph(N,M)
        n = N*M
        A = nx.adjacency_matrix(G).toarray()
        nodes = {}

        if rigid == True:  # strict lattice
            pos = list( (n) for n in G.nodes() )
            for i in range(n):
                nodes[i] = Node(i, np.random.normal(1, 0.01), np.random.normal(1, 0.01), (pos[i][1], pos[i][0]), 1)
                
            self.fixed = []

        else:
            pos_fixed = {} 
            # fix the positions of the nodes on the frame as an axis
            half = int(N / 2)
            for i in range(N):
                pos_fixed[(i, half)] = (i, half)
                pos_fixed[(half, i)] = (half, i)
            # fix corners too
            pos_fixed[(0,0)] = (0,0)
            pos_fixed[(N-1,0)] = (N-1,0)
            pos_fixed[(0,N-1)] = (0,N-1)
            pos_fixed[(N-1,N-1)] = (N-1,N-1)

            fixed = list(pos_fixed.keys())
            self.fixed = fixed
            pos = nx.spring_layout(G, k=1/(N**0.5), fixed=fixed, pos=pos_fixed)  # k term forces optimal distance between nodes

            for i in range(N):
                for j in range(N):
                    nodes[i*N + j] = Node(i*N + j, np.random.normal(1, 0.1), np.random.normal(1, 0.1), (pos[(j,i)][0], pos[(j,i)][1]), 1)

        self.Adj = A
        self.N = n
        self.Nodes = nodes
        self.rigid = rigid
        self.denser = denser
        self.init_dim = N
    
    def update_pos(self):
        G = nx.from_numpy_array(self.Adj)
        pos_fixed = {}
        for n in self.fixed:
            pos_fixed[n[0]*self.init_dim + n[1]] = n
        fixed = list(pos_fixed.keys())
        pos = nx.spring_layout(G, k=1/(self.init_dim**0.5), fixed=fixed, pos=pos_fixed)  # k term forces optimal distance between nodes
        for n in self.Nodes:
            self.Nodes[n].update_pos(pos[n][0], pos[n][1])
        return

    def rand_dir(self):
        x = rand.getrandbits(1)
        if x == 0:
            return -1
        elif x == 1:
            return 1
        
    def check_for_node(self, pos):
        for n in self.Nodes:
            if self.Nodes[n].pos == pos:
                # print("theres a node here called", n)
                return n
        # print("theres no nodes here")
        return "0"
    
    def delete_node(self, node):
        """Remove all connections to and from node but keep entry in dictionary"""
        Adj = self.Adj
        N = self.N
        Nodes = self.Nodes

        # delete row and column from Adjacency matrix
        Adj[0:N][node] = 0
        Adj[node][0:N] = 0

        Nodes[node].update_node(0,0)

        self.Adj = Adj
        return Adj

    def add_node(self, node):
        Adj = self.Adj
        N = self.N
        Nodes = self.Nodes

        #diagonal = True

        neighbour_conc_u = []
        neighbour_conc_v = []

        average_conc_u = 0
        average_conc_v = 0
        for i in Nodes:
            if i != node.name:
                average_conc_u += Nodes[i].conc_u/(N-1)
                average_conc_v += Nodes[i].conc_v/(N-1)

        # Cycle through neighbouring nodes, record average of their concentrations, and dont place nodes on top 
        # of one another, if all eight node spaces taken up then delete the node. Choice whether to grow pattern 
        # inwards as well as outwards (denser option). I think denser and diagonal are equivalent
        new_node_pos = (node.pos[0], node.pos[1]+1)
        n = self.check_for_node(new_node_pos)
        if n != "0":
            neighbour_conc_u.append(Nodes[n].conc_u)
            neighbour_conc_v.append(Nodes[n].conc_v)
            new_node_pos = (node.pos[0]+1, node.pos[1])
            n = self.check_for_node(new_node_pos)
            if n != "0":
                neighbour_conc_u.append(Nodes[n].conc_u)
                neighbour_conc_v.append(Nodes[n].conc_v)
                new_node_pos = (node.pos[0], node.pos[1]-1)
                n = self.check_for_node(new_node_pos)
                if n != "0":
                    neighbour_conc_u.append(Nodes[n].conc_u)
                    neighbour_conc_v.append(Nodes[n].conc_v)
                    new_node_pos = (node.pos[0]-1, node.pos[1])
                    n = self.check_for_node(new_node_pos)
                    if n != "0":
                        neighbour_conc_u.append(Nodes[n].conc_u)
                        neighbour_conc_v.append(Nodes[n].conc_v)

                        if self.denser == False: # maximum reached
                            #self.Adj = self.delete_node(node.name)
                            node.update_node(1, 1)
                            #print("NODE TO DELETE ----------------------------------") # given i cant delete them i can just bring the concentration back down
                            return
                        
                        #if diagonal==True:
                        new_node_pos = (node.pos[0]+0.5, node.pos[1]+0.5)
                        n = self.check_for_node(new_node_pos)
                        if n != "0":
                            neighbour_conc_u.append(Nodes[n].conc_u) 
                            neighbour_conc_v.append(Nodes[n].conc_v)
                            new_node_pos = (node.pos[0]+0.5, node.pos[1]-0.5)
                            n = self.check_for_node(new_node_pos)
                            if n != "0":
                                neighbour_conc_u.append(Nodes[n].conc_u) 
                                neighbour_conc_v.append(Nodes[n].conc_v)
                                new_node_pos = (node.pos[0]-0.5, node.pos[1]-0.5)
                                n = self.check_for_node(new_node_pos)
                                if n != "0":
                                    neighbour_conc_u.append(Nodes[n].conc_u) 
                                    neighbour_conc_v.append(Nodes[n].conc_v)
                                    new_node_pos = (node.pos[0]-0.5, node.pos[1]+0.5)
                                    if self.check_for_node(new_node_pos) != "0":
                                        #self.Adj = self.delete_node(node.name)
                                        node.update_node(1, 1)
                                        #print("NODE TO DELETE ----------------------------------")
                                        return
    
        if len(neighbour_conc_u) > 0:
            new_conc_u = sum(neighbour_conc_u)/len(neighbour_conc_u)
            new_conc_v = sum(neighbour_conc_v)/len(neighbour_conc_v)
        else:
            new_conc_u = 1
            new_conc_v = 1
        

        # Add row and column to adjacency matrix, connect with any neighbours
        Adj = np.append(Adj, np.zeros([N,1], dtype=int), axis=1)  # add column of zeros
        Adj = np.append(Adj, np.zeros([1,N+1], dtype=int) , axis=0)  # add row of zeros

        possible_connection = self.check_for_node((new_node_pos[0]+1, new_node_pos[1]))
        if possible_connection != "0":
            Adj[possible_connection][N] = 1
            Adj[N][possible_connection] = 1
        possible_connection = self.check_for_node((new_node_pos[0]-1, new_node_pos[1]))
        if possible_connection != "0":
            Adj[possible_connection][N] = 1
            Adj[N][possible_connection] = 1
        possible_connection = self.check_for_node((new_node_pos[0], new_node_pos[1]+1))
        if possible_connection != "0":
            Adj[possible_connection][N] = 1
            Adj[N][possible_connection] = 1
        possible_connection = self.check_for_node((new_node_pos[0], new_node_pos[1]-1))
        if possible_connection != "0":
            Adj[possible_connection][N] = 1
            Adj[N][possible_connection] = 1

        if self.denser == True:
            possible_connection = self.check_for_node((new_node_pos[0]+0.5, new_node_pos[1]+0.5))
            if possible_connection != "0":
                Adj[possible_connection][N] = 1
                Adj[N][possible_connection] = 1
            possible_connection = self.check_for_node((new_node_pos[0]+0.5, new_node_pos[1]-0.5))
            if possible_connection != "0":
                Adj[possible_connection][N] = 1
                Adj[N][possible_connection] = 1
            possible_connection = self.check_for_node((new_node_pos[0]-0.5, new_node_pos[1]-0.5))
            if possible_connection != "0":
                Adj[possible_connection][N] = 1
                Adj[N][possible_connection] = 1
            possible_connection = self.check_for_node((new_node_pos[0]-0.5, new_node_pos[1]+0.5))
            if possible_connection != "0":
                Adj[possible_connection][N] = 1
                Adj[N][possible_connection] = 1


        # Average all nodes
        # Nodes[N] = Node(N, np.random.normal(average_conc_u, 0.1), np.random.normal(average_conc_v, 0.1), new_node_pos, 1)
        # node.update_node(np.random.normal(average_conc_u, 0.1), np.random.normal(average_conc_v, 0.1))

        # Average of neighbours
        # Nodes[N] = Node(N, np.random.normal(new_conc_u, 0.1), np.random.normal(new_conc_v, 0.1), new_node_pos, 1)
        # node.update_node(np.random.normal(new_conc_u, 0.1), np.random.normal(new_conc_v, 0.1))

        # Same conc for old, 1 for new
        # Nodes[N] = Node(N, np.random.normal(1, 0.1), np.random.normal(1, 0.1), new_node_pos, 1)

        # 1 for old and new
        Nodes[N] = Node(N, np.random.normal(1, 0.1), np.random.normal(1, 0.1), new_node_pos, 1)
        node.update_node(np.random.normal(1, 0.1), np.random.normal(1, 0.1))


        self.N = N + 1
        self.Adj = Adj
        self.Nodes = Nodes
        
    def add_node_2(self, node):

        Adj = self.Adj
        N = self.N
        Nodes = self.Nodes

        Adj = np.append(Adj, np.zeros([N,1], dtype=int), axis=1)  # add column of zeros
        Adj = np.append(Adj, np.zeros([1,N+1], dtype=int) , axis=0)  # add row of zeros
        Adj[N][node.name] = 1  # connected to parent node
        Adj[node.name][N] = 1 
        
        if node.pos[0] == 0 and node.pos[1] < self.init_dim:
            print(node.pos)

        if node.o == 1 or node.o ==2:  # if daughter of member of original lattice nodes
            x = self.rand_dir()
            x=1
            y = self.rand_dir()
            y=-1
            if 0 <= node.name + x + y*self.init_dim < self.init_dim**2:  # if inside lattice 
                Adj[node.name + x + y*self.init_dim][N] = 1
                Adj[N][node.name + x + y*self.init_dim] = 1
                Nodes[N] = Node(N, np.random.normal(1, 0.01), np.random.normal(1, 0.01), (node.pos[0]+x*0.5, node.pos[1]+y*0.5), 2)
                node.update_node(np.random.normal(1, 0.01), np.random.normal(1, 0.01))
                
              
            elif 0 < node.name + y*self.init_dim < self.init_dim**2:  # if vertical edge % N == 0 or N-1
                Adj

        self.N = N + 1
        self.Adj = Adj
        self.Nodes = Nodes

        if self.rigid == False:
            self.update_pos()


class Tree_model:
    
    def initialise(self, n=3):
        """Initialise line of starting length n"""
        nodes = {}
        for i in range(n):
            nodes[i] = Node(i, np.random.normal(1, 0.01), np.random.normal(1, 0.01), (i,i), 1)

        # initialise adjacency matrix
        A = np.zeros([n,n], dtype=int)  # zero adjacency matrix of size n
        for i in range(n):
            if i != 0:
                A[i,i-1] = 1
            if i != n-1:
                A[i,i+1] = 1

        self.Adj = A
        self.Nodes = nodes
        self.N = n
        return A, nodes, n

    def add_daughter_adj(self, node, Adj, N, backward):
        """ Add column and row in Adjacency matrix to add node to network"""
        #Adj = self.Adj
        #N = self.N

        neighbours = []
        for i in range(N):
            if Adj[node.name, i] == 1:
                neighbours.append(i)
        node.neighbours = len(neighbours)

        if backward == False:
            if node.pos == (0,0): # allow forward growth only
                return Adj
            
        if node.neighbours == 4: # allow no more than 4 neighbours
            return Adj

        Adj = np.append(Adj, np.zeros([N,1], dtype=int), axis=1)  # add column of zeros
        Adj = np.append(Adj, np.zeros([1,N+1], dtype=int) , axis=0)  # add row of zeros
    
        Adj[-1, node.name] = 1  # branch from mother node
        Adj = np.tril(Adj) + np.triu(Adj.T, 1)  # make symmetrical

        self.Adj = Adj

        return Adj

    def add_daughter_node(self, node, Adj, N, node_dict, backward):
        """Split node, assign same concentration as parent node, and position near it"""
        if backward == False:  # allow only forward growth, reset concentration
            if node.pos == (0,0):
                node.update_node(np.random.normal(0.8,0.03), np.random.normal(2,0.03))
                return node_dict, N
            
        # average_conc_u = 0
        # average_conc_v = 0
        # for i in node_dict:
        #     if i != node.name:
        #         average_conc_u += node_dict[i].conc_u/(N-1)
        #         average_conc_v += node_dict[i].conc_v/(N-1)

            
        if node.neighbours < 4:  # allow no more than 4 neighbours

            if node.pos[0] > 0:  # grow in positive direction
                direction = 1
            else:  # grow in negative direction
                direction = -1

            if node.neighbours == 1:  # grow in straight line
                order = node.o
                dir = node.dir
            elif node.neighbours == 2:  # branch right
                order = node.o + 1
                dir = node.dir - math.pi/4
            elif node.neighbours == 3:  # branch left
                order = node.o + 1
                dir = node.dir + math.pi/4
            
            x = node.pos[0] + direction * math.cos(dir)
            y = node.pos[1] + direction * math.sin(dir)

            node.update_node(np.random.normal(0.8,0.03), np.random.normal(2,0.03))
            node_dict[N] = Node(N, np.random.normal(0.8,0.03), np.random.normal(2,0.03), (x,y), order, dir=dir)
            N = N + 1
            self.Nodes = node_dict
            self.N = N
        else:
            node.update_node(np.random.normal(0.8,0.03), np.random.normal(2,0.03))

        return node_dict, N

    def add_node(self, node, backward=False):
        """Add entry to adjacency matrix and node dict"""
        adj = self.Adj
        n = self.N
        nodes_dict = self.Nodes

        Adj = self.add_daughter_adj(node, adj, n, backward)
        Node_dict, N = self.add_daughter_node(node, Adj, n, nodes_dict, backward)
        
        self.Adj = Adj
        self.Nodes = Node_dict
        self.N = N

        return Adj, Node_dict, N


class Ring:
    def initialise(self, N=6):
        np.random.seed(5)
        r = (np.random.random(N) - 0.5)/10
        #r = np.random.normal(0, 0.05, N)
        np.random.seed(2)
        r2 = (np.random.random(N) - 0.5)/10
        #r2 = np.random.normal(0, 0.05, N)

        nodes = {}
        for i in range(N):
            nodes[i] = Node(i, 1+r[i], 1+r2[i], (math.sin(i*2*math.pi/N), math.cos(i*2*math.pi/N)), 1)

        # initialise adjacency matrix
        A = np.zeros([N,N], dtype=int)  # zero adjacency matrix of size n
        for i in range(N):
            A[i,i-1] = 1

            if i == N-1:
                I = -1
            else:
                I = i
            A[i,I+1] = 1

        # add fixed nodes
        pos_fixed = {} 
        for i in range(N):
            pos_fixed[i] = nodes[i].pos
        self.pos_fixed = pos_fixed
        self.fixed = list(pos_fixed.keys())

        self.Adj = A
        self.Nodes = nodes
        self.N = N
        self.OriginalN = N
        return A, nodes, N

    def change_node_num(self, parent_node):
        for node in self.Nodes:
            if node.name > parent_node.name:
                print(node)

    def add_node(self, node):

        Adj = self.Adj
        N = self.N
        Nodes = self.Nodes

        neighbours = []
        for i in range(N):
            if self.Adj[node.name, i] == 1:
                neighbours.append(i)
        
        # new node between itself and highest conc neighbour
        if Nodes[neighbours[0]].conc_u > Nodes[neighbours[1]].conc_u:
            hi = 0
        else:
            hi = 1

        # add entry to adjacency matrix
        Adj = np.append(Adj, np.zeros([N,1], dtype=int), axis=1)  # add column of zeros
        Adj = np.append(Adj, np.zeros([1,N+1], dtype=int) , axis=0)  # add row of zeros
        Adj[node.name][neighbours[hi]] = 0  # break with highest
        Adj[neighbours[hi]][node.name] = 0
        Adj[node.name][N] = 1  # insert new inbetween
        Adj[N][node.name] = 1
        Adj[neighbours[hi]][N] = 1
        Adj[N][neighbours[hi]] = 1

        # add node to node dict
        # new_conc_u = (node.conc_u + Nodes[neighbours[hi]].conc_u)/2
        # new_conc_v = (node.conc_v + Nodes[neighbours[hi]].conc_v)/2
        new_conc_u = node.conc_u/2
        new_conc_v = node.conc_v/2
        node.update_node(node.conc_u/2, node.conc_v/2)

        """Position spring layout"""
        # position
        Nodes[N] = Node(N, new_conc_u, new_conc_v, (0,0), 1) # add random position, it will get rewritten

        if node.name in self.fixed:
            self.fixed.remove(node.name)
            self.pos_fixed.pop(node.name)
        G = nx.from_numpy_array(Adj)
        pos = nx.spring_layout(G, k=0.11, fixed=self.fixed, pos=self.pos_fixed)
        
        for i in range(N+1):
            if i not in self.fixed:
                Nodes[i].update_pos(pos[i][0], pos[i][1])

        """Position half"""
        #Nodes[N] = Node(N, new_conc_u, new_conc_v, ((node.pos[0] + Nodes[neighbours[hi]].pos[0])/2, (node.pos[1] + Nodes[neighbours[hi]].pos[1])/2), 1)
        

        self.Adj = Adj
        self.Nodes = Nodes
        self.N = N+1


class RD:

    def __init__(self, function, a, b, d_u, d_v, network):
        #make du epsilon and dv sigma in inputs
        #outputs epsilon = du, sigma*epsilon = dv, sigma = dv/du only change sigma to change ratio
        #in examples so far epsilon = 0.01, sigma = 20
        self.function = function
        self.step = 0
        self.a = a
        self.b = b
        self.t = 0.01
        self.d_u = d_u
        self.epsilon = d_u
        self.d_v = d_v
        self.sigma = d_v/d_u
        self.network = network

    def update_t(self):
        self.step +=1
    
    def reaction(self, x0, y0):
        a = self.a
        b = self.b
        ti, tf = 0.0, self.t
        sol = odeint(self.function, (x0, y0), (ti,tf), (a,b), tfirst=True)
        return sol[1]

    def diffusion3(self, L, node_dict):
        
        for i in range(len(L)):
            delta_u = [self.epsilon*L[i][j]*(node_dict[j].conc_u-node_dict[i].conc_u) for j in range(len(L))] 
            delta_v = [self.epsilon*self.sigma*L[i][j]*(node_dict[j].conc_v-node_dict[i].conc_v) for j in range(len(L))] 
            node_dict[i].update_node(node_dict[i].conc_u + sum(delta_u), node_dict[i].conc_v + sum(delta_v))

        return node_dict

    def diffusion2(self, L ,node_dict):
        Du=self.d_u
        Dv=self.d_v
        delta_u = []
        delta_v = []
        for i in range(len(L)):
            du=0
            dv=0
            for j in range(len(L)):
                du += Du * L[i][j] * (node_dict[j].conc_u - node_dict[i].conc_u)/-(L[i][i]-1)
                dv += Dv * L[i][j] * (node_dict[j].conc_v - node_dict[i].conc_v)/-(L[i][i]-1)
            delta_u.append(du)
            delta_v.append(dv)
        for i in range(len(delta_u)):
            new_conc_u = node_dict[i].conc_u + delta_u[i]
            new_conc_v = node_dict[i].conc_v + delta_v[i]
            node_dict[i].update_node(new_conc_u, new_conc_v)
        return node_dict

    def diffusion(self, L, node_dict):
        Du=self.d_u, 
        Dv=self.d_v
        delta_u = []
        delta_v = []
        for i in L:
            du = 0
            dv = 0
            for j in range(len(i)):
                du += Du * i[j] * node_dict[j].conc_u
                dv += Dv * i[j] * node_dict[j].conc_v
            delta_u.append(du)
            delta_v.append(dv)
        for i in range(len(delta_u)):
            new_conc_u = node_dict[i].conc_u + delta_u[i]
            new_conc_v = node_dict[i].conc_v + delta_v[i]
            node_dict[i].update_node(new_conc_u, new_conc_v)
        return node_dict
    
    def reaction_diffusion(self, node_dict, L):
        for i in range(len(node_dict)):
            conc_u, conc_v = self.reaction(node_dict[i].conc_u, node_dict[i].conc_v)
            node_dict[i].update_node(conc_u, conc_v)
        node_dict = self.diffusion3(L, node_dict)
        return node_dict


def simple_update(num, n, layout, G, ax):
    ax.clear()

    random_colors = np.random.randint(2, size=n)
    nx.draw(G, pos=layout, node_color=random_colors, ax=ax)
    ax.set_title("Frame {}".format(num))

def simple_animation():

    # Build plot
    fig, ax = plt.subplots(figsize=(6,4))

    # Create a graph and layout
    n = 30 # Number of nodes
    m = 70 # Number of edges
    G = nx.gnm_random_graph(n, m)
    layout = nx.spring_layout(G)

    ani = animation.FuncAnimation(fig, simple_update, frames=10, fargs=(n, layout, G, ax))
    ani.save('animation_1.gif', writer='imagemagick')

    plt.show()

def brusselator(t, X, a, b):
    """Return the derivatives dx/dt and dy/dt."""
    x, y = X
    dxdt = a - (1+b)*x + x**2 * y
    dydt = b*x - x**2 * y
    return dxdt, dydt

def Laplacian(A):
    Laplacian = A.copy()
    for i in range(A.shape[0]):
        Laplacian[i, i] = A[i].sum()*-1
    return Laplacian

def draw_network(node_dict, Adj, filename):

    vmin = 0.5
    vmax = 2

    positions = {}
    attribute = []
    for i in range(len(node_dict)):
        positions[i] = (node_dict[i].pos)
        attribute.append(node_dict[i].conc_u)
    # if max(attribute) < vmax:
    #     vmax = max(attribute)
    # if min(attribute) > vmin:
    #     vmin = min(attribute)
    Adj = np.tril(Adj) + np.triu(Adj.T, 1)  # ensure lower Triangle = upper triangle (not directed graph)
    np.fill_diagonal(Adj, 0)  # ensure nodes not connected to themselves
    G = nx.from_numpy_array(Adj)
    #nx.draw(G, node_size=20, node_color='black', pos=positions)

    #bounds = np.linspace(min(attribute), max(attribute), 100)
    bounds = np.linspace(vmin, vmax, 100)
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=100)

    cmap = mpl.cm.cool
    #colors = plt.get_cmap('cool')(np.linspace(15, 25, 3))
    norm = mpl.colors.Normalize(vmin=min(bounds), vmax=max(bounds))
    nx.draw(G, with_labels=False, pos=positions, node_color=attribute, cmap=cmap, node_size=50, edge_color=(0,0,0,.2), vmin=vmin, vmax=vmax)
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), shrink=0.9)
    plt.savefig(filename)
    plt.show()

    return G

def plot1d(y):
    t = np.arange(0, len(y[0]))
    for i in y:
        plt.plot(t, i)
    plt.show()

def node_map(dict, filename):
    for i in dict:
        print(i)

def plot2d(dict, filename):
    for i in dict:
        plt.plot(dict[i][0], dict[i][1])
    plt.savefig(filename)
    plt.show()

def Eig(G):
    Laplacian = nx.normalized_laplacian_matrix(G)  # Normalised Laplacian
    e = np.linalg.eigvals(Laplacian.A)  # Eigenvalues
    e = sorted(e)
    e = [-i for i in e]
    print("Largest eigenvalue:", max(e))
    print("Smallest eigenvalue:", min(e))
    print(Laplacian)

    plt.scatter(np.linspace(0,len(e), len(e)), e)

    """Histogram"""
    # plt.hist(e, bins=30)  # histogram with 100 bins
    # plt.plot(np.linspace(0,len(e), len(e)), e)
    # plt.xlim(0, 1.5)  # eigenvalues between 0 and 2
    # plt.ylim(0, 5)
    # plt.title("Histogram of Eigenvalues for Balanced Tree Network")
    #plt.savefig(f"{images_dir}/bal_tree_eig.png")


    plt.show()




# G = nx.erdos_renyi_graph(100, 0.1)
# pos = nx.spring_layout(G)
# for node in G.nodes:

# nx.draw(G, node_size=10)
# plt.show()



# Network = Lattice()

# Network.initialise(3)

# Network.add_node_2(Network.Nodes[3])
# Network.add_node_2(Network.Nodes[0])
# Network.add_node_2(Network.Nodes[0])
# Network.add_node_2(Network.Nodes[11])

# G = draw_network(Network.Nodes, Network.Adj)





# L = Laplacian(Network.Adj)
# print(L)
# RD = RD(brusselator)

# to_plot1=[]
# to_plot2=[]

# for i in range(100):
#     if i%20 == 0:
#         draw_network(Network.Nodes, Network.Adj)
#     to_delete = []
#     to_add = []
#     for n in Network.Nodes:
#         print(len(Network.Nodes))
#         if Network.Nodes[n].conc_u < 0:
#             to_delete += [n]
#         elif Network.Nodes[n].conc_u > 5:
#             to_add += [n]
#     print(to_add)
#     for n in to_delete:
#         Network.delete_node(Network.Nodes[n].name)
#         print("deleted")
#     for n in to_add:     # maybe add laplacian here
#         Network.add_node(Network.Nodes[n])
#         print("added")
        
#     RD.reaction_diffusion(node_dict=Network.Nodes, L=L)
#     to_plot1.append(Network.Nodes[0].conc_u)
#     to_plot2.append(Network.Nodes[0].conc_v)
#     print(i)

# G = draw_network(Network.Nodes, Network.Adj)

# plot1d([to_plot1, to_plot2])



# Eig(G)




# def step():
#     nodes = RD.reaction_diffusion(node_dict=nodes, L=L)
#     draw_network(nodes, A, ax)

# anim = animation.FuncAnimation(fig, step, frames=np.arange(10), interval=20)
# anim.save(filename="graph.gif", dpi=60, fps=1, writer='imagemagick')
# plt.close()

# draw_network(nodes, A)

