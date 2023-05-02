from Functions import *

class Params:
    def __init__(self, a, b, d_u, d_v, network):
        # brusselator
        self.a = a
        self.b = b

        # diffusion
        self.d_u = d_u
        self.d_v = d_v

        # concentrations
        self.threshold = 0.5
        
        # run
        self.time_range = 10000
        self.plot_every = int(self.time_range/10)

        # network
        self.network = network
        if network == "Modular":
            self.Network = Modular_network()
        elif network == "Lattice":
            self.Network = Lattice()
        elif network == "Tree":
            self.Network = Tree_model()
        elif network == "Ring":
            self.Network = Ring()
        else:
            self.Network = Modular_network()
            print("please set the network parameter to 'Modular', 'Lattice' or 'Tree'")
            exit()


"""LATTICE"""
def run_lat(grow=True, denser=False):

    Network = params.Network
    
    Network.initialise(10,10, denser=denser)    

    L = Laplacian(Network.Adj)

    rd = RD(brusselator, a=params.a, b=params.b, d_u=params.d_u, d_v=params.d_v, network=params.network)

    to_plot = {}
    to_plotv = {}
    for i in range(10):
        to_plot[i] = [[],[]]
        to_plotv[i]= [[],[]]
    for i in range(params.time_range):

        for n in Network.Nodes:
            if n < 10:
                to_plot[n][0].append(i)
                to_plot[n][1].append(Network.Nodes[n].conc_u)
                to_plotv[n][0].append(i)
                to_plotv[n][1].append(Network.Nodes[n].conc_v)

        L = Laplacian(Network.Adj)

        rd.reaction_diffusion(node_dict=Network.Nodes, L=L)
        
        if grow == True:
            if i%100 == 0:
                to_delete = []
                to_add = []
                for n in Network.Nodes:
                    if Network.Nodes[n].conc_u == 0:
                        to_delete += [n]
                    if Network.Nodes[n].conc_u > 2:
                        to_add += [n]
                for n in to_delete:
                    Network.delete_node(Network.Nodes[n].name)
                for n in to_add: 
                    Network.add_node(Network.Nodes[n])
        if i > 999 and i%params.plot_every == 0:
            draw_network(Network.Nodes, Network.Adj, f"{params.network}/TimeLapse0.8/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}.png")
            plot2d(to_plot, f"{params.network}/TimeLapse0.8/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}-NODES.png")   
            plot2d(to_plotv, f"{params.network}/TimeLapse0.8/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}-VNODES.png")   
        # for n in Network.Nodes:
        #     to_plot[n][0].append(i)
        #     to_plot[n][1].append(Network.Nodes[n].conc_u)
        #     #print(max(to_plot[n][1]))
        print(i)
        max = []
        

    draw_network(Network.Nodes, Network.Adj, f"{params.network}/TimeLapse1/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}.png")

    plot2d(to_plot, f"{params.network}/TimeLapse0.8/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}-NODES.png")
    plot2d(to_plotv, f"{params.network}/TimeLapse0.8/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}-VNODES.png")  
# params = Params(a=1, b=2, d_u=0.004, d_v=0.3, network="Lattice")         
# run_lat()

"""TREE"""
def run_tree(grow=True):
    Network = params.Network

    Network.initialise(30)

    L = Laplacian(Network.Adj)

    rd = RD(brusselator, a=params.a, b=params.b, d_u=params.d_u, d_v=params.d_v, network=params.network)

    to_plot = {}
    to_plot_v = {}
    for i in range(1000):
        to_plot[i] = [[],[]]
        to_plot_v[i] = [[],[]]

    for i in range(params.time_range):
        rd.reaction_diffusion(node_dict=Network.Nodes, L=L)
        L = Laplacian(Network.Adj)
        if grow == True:
            to_delete = []
            to_add = []
            for n in Network.Nodes:
                if Network.Nodes[n].conc_u < 0.1:
                    to_delete += []
                elif Network.Nodes[n].conc_u > 2:
                    to_add += []
            for n in to_delete:
                Network.delete_node(Network.Nodes[n].name)
            for n in to_add: 
                Network.add_node(Network.Nodes[n])
        # if i%params.plot_every == 0 and i>999:
        #     draw_network(Network.Nodes, Network.Adj, f"{params.network}/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}.png")
        #     plot2d(to_plot, f"{params.network}/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}-NODES.png")   
        #     plot2d(to_plot_v, f"{params.network}/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}-VNODES.png") 
        for n in Network.Nodes:
            to_plot[n][0].append(i)
            to_plot[n][1].append(Network.Nodes[n].conc_u)
            to_plot_v[n][0].append(i)
            to_plot_v[n][1].append(Network.Nodes[n].conc_v)
        print(i)
    draw_network(Network.Nodes, Network.Adj, f"{params.network}/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}.png")
    plot2d(to_plot, f"{params.network}/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}-NODES.png")
    plot2d(to_plot_v, f"{params.network}/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}-VNODES.png")
params = Params(a=1, b=1.5, d_u=0.01, d_v=0.2, network="Tree") 
run_tree()

"""MODULAR"""
def run_mod():
    Network = params.Network
    Network.initialise(8,8)
    L = Laplacian(Network.Adj)
    rd = RD(brusselator, a=params.a, b=params.b, d_u=params.d_u, d_v=params.d_v, network=params.network)

    """run"""
    to_plot = {}
    for i in range(800):
        to_plot[i] = [[],[]]
    for i in range(params.time_range):
        for n in Network.Nodes:
            to_plot[n][0].append(i)
            to_plot[n][1].append(Network.Nodes[n].conc_u)
        rd.reaction_diffusion(node_dict=Network.Nodes, L=L)
        L = Laplacian(Network.Adj)
        to_delete = []
        to_add = []
        if i%10 == 0:
            for n in Network.Nodes:
                if Network.Nodes[n].conc_u < 0.2:
                    to_delete += []
                elif Network.Nodes[n].conc_u > 2:
                    to_add += [n]
            for n in to_delete:
                Network.delete_node(Network.Nodes[n].name)
            for n in to_add: 
                Network.add_node(Network.Nodes[n])
                print("added a node", n, to_add)
            if to_add != []:
                draw_network(Network.Nodes, Network.Adj, f"NewDiffusion/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}.pdf")
                #break # add one node at a time - take this away later
        if i%5000 == 0:
            draw_network(Network.Nodes, Network.Adj, f"NewDiffusion/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}.pdf")
            plot2d(to_plot, f"NewDiffusion/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}-NODES.pdf")
        print(i)
    G = draw_network(Network.Nodes, Network.Adj, f"NewDiffusion/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}.pdf")
    plot2d(to_plot, f"NewDiffusion/a={params.a}_b={params.b}_du={params.d_u}_dv={params.d_v}_t={i}-NODES.pdf")
    return G
# params = Params(a=1, b=1.5, d_u=0.004, d_v=0.1, network="Modular") 
# G = run_mod()
#Eig(G)


"""RING"""
def run_ring(number):
    Network = params.Network
    Network.initialise(number)
    #Network.Nodes[0].update_node(1.02,0.98)
    rd = RD(brusselator,a=params.a, b=params.b, d_u=params.d_u, d_v=params.d_v, network=params.Network)

    to_plot = {}
    for i in range(200):
        to_plot[i] = [[],[]]

    # for i in range(1):
    #     Network.add_node(Network.Nodes[0])

    #draw_network(Network.Nodes, Network.Adj, f"rec/ring{i}.png")
    counter = 0
    for i in range(400000):
        for n in Network.Nodes:
            to_plot[n][0].append(i)
            to_plot[n][1].append(Network.Nodes[n].conc_u)
        
        rd.reaction_diffusion(Network.Nodes, L=Laplacian(Network.Adj))
        # if i%10 == 0:
        #     highest = [0,0]
        #     for n in range(Network.N):
        #         if Network.Nodes[n].conc_u > 2 and Network.Nodes[n].conc_u > highest[1]:
        #             highest = [Network.Nodes[n],Network.Nodes[n].conc_u]
        #     if highest[1] != 0:
        #         Network.add_node(highest[0])
        #         #draw_network(Network.Nodes, Network.Adj, f"ring{i}.png")
        #         counter += 1
        #         print(counter)
            
        
        if i %10000 == 0 and i > 19000:
            draw_network(Network.Nodes, Network.Adj, f"rec2/ringstationary20.40.png")
            plot2d(to_plot, "rec2/ringnodes20.40.png")
        print(i)
# params = Params(a=1, b=1.5, d_u=0.01, d_v=0.2, network="Ring") 
# run_ring(15)
