import numpy as np
from .MarkovChain import MarkovChain


#########################################################################
#                                                                       #
#                            Shahram Emami                              #
#                                                                       #
# Markov chain (CTMC) model of TF binding to DNA. Only the single-     #
# stranded model (SSM) is implemented here.                             #
#########################################################################


class SingleStrandChain(MarkovChain):

  #def __init__(self,verbose,L,q01,q10,q11,q12,q21,TF_concentration):
  def __init__(self,L,q01,q10,q11,q12,q21,TF_concentration):
        super().__init__()

        self.L = L #Length of the DNA sequence
        self.q01 = q01
        self.q10 = q10
        self.q11 = q11
        self.q12 = q12
        self.q21 = q21
        self.TF_conc = TF_concentration

        # Construct the Q Matrix for SSM case
        self.__generate_SSM_Q()



        # Create the state to index and index to state dictionaries
        self.__generate_SSM_Dictionaries()

        return
  
#########################################################################
########################################################################

  def __generate_SSM_Q(self):

        # This function generates the Infinitesimal Generator matrix Q for Single Stranded Model.
        # Matrix Q is (2 * L + 1) by (2 * L + 1) matrix, where the
        # labels for the columns and rows are as follows:
        # "Free DNA, NS_1, S_1, NS_2, S_2, ..., NS_L, S_L,"
        # state 0: Free DNA
        # state 1: (Non_Specific Binding) NS
        # state 2: (Specific Binding) SP

        Q = np.zeros((2 * self.L + 1, 2 * self.L + 1))

        # Set the first row of the matrix
        Q[0][0] = self.TF_conc * -1 * self.L * self.q01

        for i in range(1,self.L+1):
          Q[0][i * 2 - 1] = self.TF_conc * self.q01

        # Set the first column of the matrix
        for j in range(1,self.L+1):
          Q[j * 2 - 1][0] = self.q10


        # Now set the Diffusion rates:
        if (self.L > 1):
          # From the first NS position you can only diffuse to the second NS position
          Q[1][3] = self.q11

          # From the Lth NS position you can only diffuse to the (L - 1)th NS position
          Q[2 * self.L - 1][2 * self.L - 1 - 2] = self.q11
        # end of the if statement

        # For the internal NS position (1 < position < L) you could diffuse to
        # either (position - 1) OR (position + 1)
        for i in range(1,self.L - 1):
          Q[2 * i + 1][2 * i + 1 - 2] = self.q11
          Q[2 * i + 1][2 * i + 1 + 2] = self.q11
        # end of for loop

        # Now set the transition rates
        for i in range(0,self.L):
          # The diagonal entry is minus the total rate of leaving NS_i:
          # recognition (q12), unbinding (q10), and sliding (q11) to each
          # neighbouring offset. An offset has 2 neighbours in the interior,
          # 1 at either end of the DNA, and 0 when L == 1.
          n_neighbors = (i > 0) + (i < self.L - 1)
          Q[2*i+1][2*i+1] = - self.q12[i] - self.q10 - n_neighbors * self.q11

          Q[2*i+1][2*i+2] =  self.q12[i]
          # Set the rates from State 2 to State 1
          Q[2*i+2][2*i+1] =  self.q21[i]
          Q[2*i+2][2*i+2] =  -1 * self.q21[i]
        # end for loop

        # Each row of a CTMC generator matrix must sum to zero.
        assert np.allclose(Q.sum(axis=1), 0), "rows of Q must sum to zero"

        self.Q = Q
        return
  

#########################################################################


  def __generate_SSM_Dictionaries(self):
        # This function create the dictionaries for mapping index to states for the
        # indicies in the SSM Q Matrix and the correspoding inverse map.
        # Free = Free State (unbound DNA)
        # NS_A  = NonSpecific (One Aggregate state)
        # SP_x = Specific position x

        self.s2i = {}

        for index in range(1 + 2 * self.L):
          if index == 0:
            self.s2i["Free"] = 0
            continue
          if index % 2 == 0:
            self.s2i["S_" + str(int(index/2))] = index
          else:
            self.s2i["NS_" + str(int((index + 1)/2))] = index

        # Now create the inverse map
        self.i2s = {}
        for key,value in self.s2i.items():
          self.i2s[value] = key

        return
  
#########################################################################

  def get_SSM_Q(self):
    return self.Q
  
#########################################################################

  def s2i_map(self):
    return self.s2i
  
#########################################################################

  def i2s_map(self):
    return self.i2s
  
#########################################################################
