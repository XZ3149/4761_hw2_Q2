class node:
    """
    this is the class for each grid, it contain infor for backtracking and score 
    for diagonal match, insert_1, insert_3, delete_1, delete_3. 
    """
    def __init__(self, position, nucleotide_1, nucleotide_2):
        """
        Parameters
        ----------
        position : list
            In format [i,j], this record the position of this node in a 2D list
        nucletide_1 : str
            nucleotide for read 1 at this position, ATCG or X
        nucletide_2 : TYPE
            nucleotide for read 2 at this position, ATCG or X
        """
        
        self.position = position 
        self.nucleotide_1 = nucleotide_1
        self.nucleotide_2 = nucleotide_2
        self.diagonal_match = [float('-inf'), 3,'diagonal_match', None, None] #[score,preference, type,previous_type, operation]
        # more prefer match and gap of 3 in case of a tie happen.  
        self.insert_1 = [float('-inf'), 1, 'insert_1', None, None] #[score,preference, type,previous_type, operation]
        self.insert_3 = [float('-inf'), 2,'insert_3', None, None] #[score,preference, type,previous_type, operation]
        self.delete_1 = [float('-inf'), 1, 'delete_1', None, None] #[score,preference, type,previous_type, operation]
        self.delete_3 = [float('-inf'), 2, 'delete_3', None, None] #[score,preference, type,previous_type, operation]
        
        
    def __str__(self):
        return f'{self.position},{self.nucleotide_1}, {self.nucleotide_2}'


def dp_aligement(read1, read2):
    """ 
    this is the main function for aligement
    it will take two string as read1 and read2 and using dp to alige them.
    """
    # initiation
    read1 ='-' + read1 # mannually create position 0, the actual read start at 1
    read2 = '-' + read2 # mannually create position 0, the actual read start at 1
    
    Grids = []
    
    #generate a 2D-list that contain node for our computation
    for i in range(len(read1)):
        i_row = []
        for j in range(len(read2)):
            i_row += [node([i,j], read1[i], read2[j])]
        Grids += [i_row]



    #initiate the first row and column
    Grids[0][0].diagonal_match = [0, 3, 'diagonal_match', None, 'diagonal_match']
    Grids[0][0].diagonal_match = [0, 3, 'diagonal_match', None, 'diagonal_match']
    Grids[1][0].insert_1 = [-6, 2, 'insert_1', None, 'start']
    Grids[0][1].delete_1 = [-6, 2, 'delete_1', None, 'start']      
    Grids[3][0].insert_3 = [-11, 2, 'insert_3', None, 'start']
    Grids[0][3].delete_1 = [-11, 2, 'delete_3', None, 'start'] 
    

    for i in range(2,len(read1)):                
        Grids[i][0].insert_1 = [-6 - (i) * 2, 1, 'insert_1', None, 'extent']
        if i%3 ==0 and i != 3:
            Grids[i][0].insert_3 = [-11 - ((i//3) ) * 5, 2, 'insert_3', None, 'extent']
       
    for j in range(2,len(read2)):                
        Grids[0][j].delete_1 = [-6 - (j) * 2, 1, 'delete_1', None, 'extent']
        if j%3 ==0 and i != 3:
            Grids[0][j].delete_3 = [-11 - ((j//3)) * 5, 2, 'delete_3', None, 'extent'] 
            
        
     
    """
    Grids[0][0].insert_1 = [-6, 2, 'insert_1', None, 'start']
    Grids[0][0].delete_1 = [-6, 2, 'delete_1', None, 'start']      
    Grids[0][0].insert_3 = [-11, 2, 'insert_3', None, 'start']
    Grids[0][0].delete_1 = [-11, 2, 'delete_3', None, 'start'] 
    """
       
        
    # computation step, looping over the 2D-list
    for i in range(1,len(read1)):   
        for j in range(1,len(read2)):
            
            if i-1 >= 0 and j-1>=0:
                Grids[i][j].diagonal_match = move_diagonal(Grids[i-1][j-1], Grids[i][j])
            
            if i-1 >= 0:
                Grids[i][j].insert_1 = move_insert_1(Grids[i-1][j], Grids[i][j])
                
            if j-1 >= 0:
                Grids[i][j].delete_1 = move_delete_1(Grids[i][j - 1], Grids[i][j])
                
            if i >= 3:
                Grids[i][j].insert_3 = move_insert_3(Grids[i-3][j], Grids[i][j])
                
            if j >= 3:
                Grids[i][j].delete_3 = move_delete_3(Grids[i][j - 3], Grids[i][j])    
                
            
    
    #backtraking and output step 
    output_1 = ''
    output_2 = ''
    output_3 = ''
    order = ''
    score = 0
    
    current_node = Grids[len(read1) - 1][len(read2) -1] # start at the terminate node
    while current_node.position!= [0,0]:
        

        
        max_path = max(current_node.diagonal_match,current_node.insert_1, \
                       current_node.insert_3, current_node.delete_1, current_node.delete_3)
            
            
            
        if max_path[2] == "diagonal_match":
            output_1 = current_node.nucleotide_1 + output_1
            output_3 = current_node.nucleotide_2 + output_3
            if current_node.nucleotide_1 == current_node.nucleotide_2:
                output_2 = '|' + output_2
            else:
                output_2 = '.' + output_2
                
            current_node = Grids[current_node.position[0] -1][current_node.position[1] -1]
            order = ' diagonal ' + order
            
        elif max_path[2] == "insert_1":
            output_1 = current_node.nucleotide_1 + output_1
            output_3 = '-' + output_3
            output_2 = ' ' + output_2
            current_node = Grids[current_node.position[0] -1][current_node.position[1]]
            order = ' insert_1 ' + order
            
            
        elif max_path[2] == "delete_1":
            output_3 = current_node.nucleotide_2 + output_3
            output_1 = '-' + output_1
            output_2 = ' ' + output_2
            current_node = Grids[current_node.position[0]][current_node.position[1] - 1]
            order = ' delete_1 ' + order
            
        
        elif max_path[2] == "insert_3":
            previous_grid_1 = Grids[current_node.position[0] -1][current_node.position[1]]
            previous_grid_2 = Grids[current_node.position[0] -2][current_node.position[1]]
            output_1 = previous_grid_2.nucleotide_1 + previous_grid_1.nucleotide_1 + \
                current_node.nucleotide_1 + output_1
            output_2 = '   ' + output_2
            output_3 = '---' + output_3
            order = ' insert_3 ' + order
            
            current_node = Grids[current_node.position[0] -3][current_node.position[1]]
    
        elif max_path[2] == "delete_3":
            previous_grid_1 = Grids[current_node.position[0]][current_node.position[1]-1]
            previous_grid_2 = Grids[current_node.position[0]][current_node.position[1]-2]
            output_2 = previous_grid_2.nucleotide_2 + previous_grid_1.nucleotide_2 + \
                current_node.nucleotide_2 + output_2
            
            output_1 = '---' + output_3
            output_2 = '   ' + output_2
            order = ' delete_3 ' + order
            
            current_node = Grids[current_node.position[0]][current_node.position[1] - 3]
    
    
    
    print(output_1)
    print(output_2)
    print(output_3)

    return order

        
def move_diagonal(node_1,node_2):
    """ the helper function for diagonally match
    node_1 = [i-1,j-1], node_2 = [i,j]
    will return a list with [new_score,type,previous_operation, operation]
    """
    
    operation_cost = 0
    if node_2.nucleotide_1 == node_2.nucleotide_2:
        operation_cost = -3
    else:
        operation_cost = 3
    
    
    
    
    pre_infor = max(node_1.diagonal_match,node_1.insert_1, node_1.insert_3, node_1.delete_1, node_1.delete_3)

    score = pre_infor[0] - operation_cost
    pre_type = pre_infor[2]
     
    return [score , 3,'diagonal_match', pre_type, 'diagonal_match']
    
    
def move_insert_1(node_1,node_2): 
    """ the helper function for insert_1 match
    node_1 = [i-1,j], node_2 = [i,j]
    will return a list with [new_score,type,previous_operation, operation]
    """
    open_cost = 6
    extent_cost = 2
    
    
    diagonal_match = node_1.diagonal_match[:]
    diagonal_match[0] -= open_cost # open a new gap
    insert_1 = node_1.insert_1[:]
    
    pre_operation = insert_1[4]
    if pre_operation == None:#initaiton step, dealling with row 0 
        insert_1[0] -= open_cost
    else:
        insert_1[0] -= extent_cost # extend a gap 
        
    delete_1 = node_1.delete_1[:]
    delete_1[0] -= open_cost # open a gap
    
    
    pre_infor = max(diagonal_match,insert_1,delete_1)
    score = pre_infor[0] 
    pre_type = pre_infor[2]
    
    if pre_type == 'insert_1' and pre_operation != None: # if extent a gap work better
        return [score , 1,'insert_1', pre_type, 'extend']
    else:
        return [score , 1,'insert_1', pre_type, 'start']
    

def move_delete_1(node_1,node_2): 
    """ the helper function for insert_1 match
    node_1 = [i,j-1], node_2 = [i,j]
    will return a list with [new_score,type,previous_operation, operation]
    """
    open_cost = 6
    extent_cost = 2
    
    
    diagonal_match = node_1.diagonal_match[:]
    diagonal_match[0] -= open_cost # open a new gap
    insert_1 = node_1.insert_1[:]
    insert_1[0] -= open_cost # open a gap 
    delete_1 = node_1.delete_1[:]  
    pre_operation = delete_1[4]
    if pre_operation == None:#initaiton step, dealling with column 0 
        delete_1[0] -= open_cost
    else:
        delete_1[0] -= extent_cost # extend a gap 
    
    
    pre_infor = max(diagonal_match,insert_1,delete_1)
    score = pre_infor[0] 
    pre_type = pre_infor[2]
    
    if pre_type == 'delete_1'and pre_operation != None: # if extent a gap work better
        return [score , 1,'delete_1', pre_type, 'extend']
    else:
        return [score , 1,'delete_1', pre_type, 'start']    
    
    
    
    
def move_insert_3(node_1,node_2): 
    """ the helper function for insert_1 match
    node_1 = [i-3,j], node_2 = [i,j]
    will return a list with [new_score,type,previous_operation, operation]
    """
    open_cost = 11
    extent_cost = 5
    
    
    diagonal_match = node_1.diagonal_match[:]
    diagonal_match[0] -= open_cost # open a new gap
    insert_3 = node_1.insert_3[:]
    
    pre_operation = insert_3[4]
    if pre_operation == None:#initaiton step, dealling with row 0 
        insert_3[0] -= open_cost
    else:
        insert_3[0] -= extent_cost # extend a gap 
        
    delete_3 = node_1.delete_3[:]
    delete_3[0] -= open_cost # open a gap
    
    
    pre_infor = max(diagonal_match,insert_3,delete_3)
    score = pre_infor[0] 
    pre_type = pre_infor[2]
    
    if pre_type == 'insert_3' and pre_operation != None: # if extent a gap work better
        return [score , 2,'insert_3', pre_type, 'extend']
    else:
        return [score , 2,'insert_3', pre_type, 'start']
      
    
    
def move_delete_3(node_1,node_2): 
    """ the helper function for insert_1 match
    node_1 = [i,j-3], node_2 = [i,j]
    will return a list with [new_score,type,previous_operation, operation]
    """
    open_cost = 11
    extent_cost = 5
    
    
    diagonal_match = node_1.diagonal_match[:]
    diagonal_match[0] -= open_cost # open a new gap
    insert_3 = node_1.insert_3[:]
    insert_3[0] -= open_cost # open a gap 
    delete_3 = node_1.delete_3[:]
    
    pre_operation = delete_3[4]
    if pre_operation == None:#initaiton step, dealling with column 0 
        delete_3[0] -= open_cost
    else:
        delete_3[0] -= extent_cost # extend a gap 
    
    
    pre_infor = max(diagonal_match,insert_3,delete_3)
    score = pre_infor[0] 
    pre_type = pre_infor[2]
    
    if pre_type == 'delete_3'and pre_operation != None: # if extent a gap work better
        return [score , 2,'delete_3', pre_type, 'extend']
    else:
        return [score , 2,'delete_3', pre_type, 'start']    
        
    



