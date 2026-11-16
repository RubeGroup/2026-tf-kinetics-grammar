
#########################################################################
#                                                                       #
#                            Shahram Emami                              #
#                                                                       #
# Markov chain (CTMC) model of TF binding to DNA. Only the single-     #
# stranded model (SSM) is implemented here.                             #
#########################################################################

class MarkovChain():
    Q   = None	  #Transition rate matrix
    s2i = None	  #Map from state name to index
    i2s = None	  #Map from index to state name
    L   = None    #Length of the DNA sequence
    Q_trap_ls =  None # Modified Q matrix for the linear systems approach
    Q_trap_int = None # Modified Q matrix for the integration approach

    def __init__(self):
        # print("> Running the constructor of the base class")
        return
