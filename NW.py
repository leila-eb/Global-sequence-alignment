#import inputs

def read_substitution_matrix(fn):
	with open(fn,'r') as reader:
		order = reader.readline().split()
		return { line.split()[0] : { res:int(value) for res,value in zip(order,line.split()[1:]) } for line in reader if line.strip() }


'''PSEUDOCODE
define the function create_matrices. (seq1, seq2, substitution_matrix, gap_penalty)
	initialize a scores matrix (list of lists) with dimension len(seq1)+1 , len(seq2)+1 with all 0
	initialize a bt matrix (list of lists) with dimension len(seq1)+1 , len(seq2)+1 with all ''
	for each cell in the first row
		store gap_penalty multiplied by the column index in the scores matrix
		store an arrow pointing left in the bt matrix
	for each cell in the first column
		store gap_penalty multiplied by the row index in the scores matrix
		store an arrow pointing up in the bt matrix
	for each row starting from the second one
		for each column starting from the second one
			store the cell row-1,column + the gap_penalty
			store the cell row,column-1 + the gap_penalty
			store the cell row-1,column-1 + the substitution score of the corresponding residues in the sequences
			store in the cell of scores the maximum of the three values above and in the cell of bt the corresponding arrow
	return the matrix
'''

def create_matrices(seq1,seq2,substitution_matrix,gap_penalty):
	'''This function creates the scoring matrix and the backtracking matrix'''
	
	scores= [ [0 for col in range(len(seq1)+1)] for row in range(len(seq2)+1) ] #nested list comprehension to initialize the matrix with all 0s
	bt= [ ['' for col in range(len(seq1)+1)] for row in range(len(seq2)+1) ] #backtracing matrix is initialized with empty strings
	
	for j in range(1,len(seq1)+1):
		scores[0][j]= gap_penalty * j #we initialize the first row
		bt    [0][j]= 'l' #we store an l when we come from the left
		
	for i in range(1,len(seq2)+1):
		scores[i][0]= gap_penalty * i #we initialize the first col
		bt    [i][0]='u'
		
	for i in range(1,len(seq2)+1):
		for j in range(1,len(seq1)+1):
			su=scores[i-1][j] + gap_penalty #score from the row above, in the same column
			sl=scores[i][j-1] + gap_penalty #score from same row and column on the left
			sd=scores[i-1][j-1] + substitution_matrix[seq2[i-1]][seq1[j-1]] #score from the diagonal
			l= [su,sl,sd] #list of the 3 values
			scores[i][j]= max(l) #we fix inside the cell, the highest value of the 3
			bt[i][j]= 'uld'[l.index(max(l))]
	
	return scores,bt
	



'''PSEUDOCODE
Define a function nw_backtrack(seq1,seq2,scores,bt)
  initialize aln1 to an empty string
  initialize aln2 to an empty string
  
  position ourself at the last cell of the matrix bt
	while we are not in the first cell of the matrix bt
		if the cell is a d
			(we are matching)
			add to aln1 the corresponding character from the first sequence
			add to aln2 the corresponding character from the second sequence
			move diagonally in the matrix
		otherwise if it is an l
			(we are inserting a gap in the second sequence)
			add to aln1 the corresponding character from the first sequence
			add a gap to aln2
			move to the left in the matrix
		otherwise if it is an u
			(we are inserting a gap in the first sequence)
			add a gap to aln1
			add to aln2 the corresponding character from the second sequence
			move upper in the matrix
  
  reverse aln1
  reverse aln2
  return aln1 and aln2
'''

def nw_backtrack(seq1,seq2,scores,bt):
	aln1,aln2 = '', '' #tuple syntax to initialize the empty strings
	i,j = len(seq2),len(seq1) #indexes to move through the matrix
	while i != 0 and j != 0: #while we're not in the first cell of the matrix
	#while bt[i][j] != '':
		if bt[i][j] == 'd': #we are matching
			aln1 += seq1[j-1]
			aln2 += seq2[i-1]
			i -= 1
			j -= 1 #we move diagonally
		
		elif bt[i][j] == 'l':  #we insert a gap in the second seq
			aln1 += seq1[j-1]
			aln2 += '-' 
			j -= 1 #move to the left
		
		else: #bt[i][j] = 'u':  #we insert a gap in the first seq
			aln1 += '-'
			aln2 += seq2[i-1]
			i -= 1 #move upper
	
	return aln1[::-1],aln2[::-1]


seq1='ACTGG'
seq2='ACCA'
gap_penalty=-2
substitution_matrix= read_substitution_matrix('bases.txt')
scores,bt= create_matrices(seq1,seq2,substitution_matrix,gap_penalty) #I assign the result to 2 variables
for line in scores:
	print(line)
print()
for line in bt:
	print(line)
print()

aln1,aln2= nw_backtrack(seq1,seq2,scores,bt)
print(aln1)
print(aln2)
print('Score:',scores[-1][-1]) #last position of the matrix contains the score


