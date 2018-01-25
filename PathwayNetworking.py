'''
Create a random Erdős-Rényi network containing 100 nodes and visualise the network.
What are the lowest, highest and average degree values in this network?
Describe and comment your results in the text box below the code.
'''

import networkx as nx
import boolean2 as b2
import matplotlib.pyplot as plt

G = nx.erdos_renyi_graph(100, 0.2)
nx.draw_networkx(G)
print("Minimum degree:", min(nx.degree(G).values()))  
print("Maximum degree:", max(nx.degree(G).values()))
print("Average degree:", sum(nx.degree(G).values())/nx.number_of_nodes(G))

'''
Create a Barabási-Albert network containing 50 nodes and visualise the network.
Draw a scatterplot for this network where the betweenness centrality of each node is shown along the horizontal axis and the closeness centrality along the vertical axis.
Describe and comment your results in the text box below the code.
'''

ba = nx.barabasi_albert_graph(50, 2)
nx.draw_networkx(ba)
c = nx.closeness_centrality(ba)
b = nx.betweenness_centrality(ba)
lc = plt.subplots()
plt.scatter(list(b.values()),list(c.values()) )

'''
Create a network using the cellcycle_module2.txt file, which contains a list of interacting proteins (each line in the file represents a pair of proteins connected by an edge).
What are the proteins that have the highest degree value in this network?
Explain the role of hub proteins in protein interaction networks.
'''

data = open('cellcycle_module2.txt', 'r')
data.readlines()
G = nx.read_edgelist("cellcycle_module2.txt")
nx.draw_networkx(G)
d = nx.degree(G)
b = nx.betweenness_centrality(G)
print('Thehe proteins that have the highest degree value are:')
for n in G.nodes():
    if (d[n] ==  max(nx.degree(G).values())):
        print (n)
print 'The maximum value of degree in the network is', max(nx.degree(G).values()),'.  '
#The proteins that have the highest degree value are: SMAD3 and ESR1. The maximum value of degree in the network is 10 .  

'''
#Hub nodes tend to be essential nodes. Their selective attack leads to rapid deterioration of information transfer. They are a few nodes have high degree.
CHEK1 is a protein that coordinates the DNA damage response. In the scientific literature we found the following information about proteins involved in the regulation of CHEK1:
ATM activates ATR
ATM activates P53
ATR activates CHEK1
P53 activates PMID
PMID inhibits CHEK1
PMID inhibits P53
Create a Boolean model that contains these five proteins and the appropriate logic rules.
The initial state of ATM is inactive, the state of other proteins is unknown. Run the model and determine if it converges towards a steady state or a limit cycle.
Describe and comment your results in the text box below the code.
'''

model_definition="""
ATM=False
ATR = Random
PMID = Random
P53 = Random
CHEK1 = Random
ATR* = ATM
P53* = ATM
CHEK1* = ATR
PMID* = P53
CHEK1* = not PMID
P53* = not PMID
"""
model = b2.Model(text = model_definition, mode='sync')
model.initialize()
model.iterate(steps = 10)
for state in model.states:
    print state.ATM, state.ATR, state.PMID, state.P53, state.CHEK1
model.report_cycles()
    #The model coverage to a steady state.  
'''
Herceptin is a drug used for the treatment of some types of breast cancer. The action of Herceptin is modulated by dual specificity phosphatases (DUSPs), but their precise function in the cellular response to Herceptin is unknown. We hypothesize that the regulatory mechanism of DUSP1 acts as follows:
Herceptin inhibits HER2
HER2 activates ERK
HER2 activates JNK
HER2 activates P38
ERK activates DUSP1
P38 activates DUSP1
DUSP1 inhibits JNK
DUSP1 inhibits P38
ERK activates Survival
JNK and P38 inhibit Survival
Create a Boolean model that contains these proteins and the appropriate logic rules.
Herceptin is applied to the cells. Run the model and determine if it converges towards a steady state or a limit cycle. Are the cells expected to survive or not after treatment by Herceptin?
New experimental data suggest that in some mutated tumours DUSP1 is no longer able to inhibit P38. Modify the model accordingly and determine if the cells are expected to survive or not after Herceptin treatment.
What are the implications of these results for the use of Herceptin as a drug against breast cancer?
#The network model of the five proteins is appropriate generated synchronously(update the state of all nodes in the same order and the same time) and simulate 10 time steps accroding to the information as requirement( ATM activates ATR; ATM activates P53; ATR activates CHEK1; P53 activates PMID; PMID inhibits CHEK1; PMID inhibits P53; The initial state of ATM is inactive, the state of other proteins is unknown). The model coverage to a steady state at the index(3) of the state of all nodes(step 2). The state is ATM=False, ATR=False, CHEK1=False, P53=False and PMID=False.
'''

from boolean2 import util
model_definition="""
HERCEPTIN = True
HER2 = True
ERK = True
JNK = True
P38 = True
DUSP1 = True
SURVIVAL = True
HER2* = not HERCEPTIN
ERK* = HER2
JNK* = HER2
P38* = HER2
DUSP1* = ERK
DUSP1* = P38
JNK* = not DUSP1
P38* = not DUSP1
SURVIVAL*   =   ERK
SURVIVAL*  =  not JNK or P38
"""
model = b2.Model(text = model_definition, mode='async')
model.initialize()
model.iterate(steps = 30)
for node in model.data:
    print node, model.data[node]
model.report_cycles()
coll = util.Collector()
for i in range(10):
    model = b2.Model(text=model_definition, mode='async')
    model.initialize()
    model.iterate(steps=30)
    coll.collect(states=model.states, nodes=model.nodes)
avgs = coll.get_averages()
image = list()
for node in model.data:
    image.append(avgs[node])
plt.yticks(range(0,5), model.data)
plt.imshow(image, cmap=plt.cm.get_cmap('RdYlGn'), interpolation='none')

print '****************Some mutated tumours DUSP1 is no longer able to inhibit P38. Thus, we generated a new network**:'
model_definition="""
HERCEPTIN = True
HER2 = True
ERK = True
JNK = True
P38 = True
DUSP1A = True
SURVIVAL = True
HER2* = not HERCEPTIN
ERK* = HER2
JNK* = HER2
P38* = HER2
DUSP1A = ERK
DUSP1A = P38
JNK* = not DUSP1A
SURVIVAL*   =   ERK
SURVIVAL*  =  not JNK or P38
"""
model = b2.Model(text = model_definition, mode='sync')
model.initialize()
model.iterate(steps = 15)
for node in model.data:
    print node, model.data[node]
model.report_cycles()
coll = util.Collector()
for i in range(10):
    model = b2.Model(text=model_definition, mode='async')
    model.initialize()
    model.iterate(steps=30)
    coll.collect(states=model.states, nodes=model.nodes)
avgs = coll.get_averages()
image = list()
for node in model.data:
    image.append(avgs[node])
plt.yticks(range(0,5), model.data)
plt.imshow(image, cmap=plt.cm.get_cmap('RdYlGn'), interpolation='none')

'''
#The cells are not going to keep survival due to the application of HERCEPTIN. Here we generate a model of a the scheme of all the vertexs in the network in the synchronous mode in 15 steps. Before theHERCEPTIN is applied. All the  nodes in the model are on and the cells are in survival.After 4 steps of the effect of the HERCEPTIN. The model reach a steady state that all the cells died.  
In the new experimental data suggest that in some mutated tumours DUSP1 is no longer able to inhibit P38. Some cells in this case suvived. In new the model we suggest that we should also find a new medicine to work together with HERCEPTIN to inhibit P38.
'''