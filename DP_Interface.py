''' Algorithmic Thinking Part 2: Dynamic Programming - With Rice University on Coursera
Computing Global and Local Pairwise Sequence Alignment Scores of two genes'''

# This Python file uses the following encoding: utf-8

################## 
## INSTRUCTIONS ##
##################
'''
This script takes an input file when run from the command line according to the following scheme:
line 0 = number of Analyses to compute
Line 1 = the alphabet over which the sequences will be analyzed
Line 2 = the 1st gene sequence
Line 3 = the 2nd gene sequence
Line 4 = a boolean value: True for Global Alignments and False for Local Alignments
Ling 5 = integer values to be used in building the scoring matrix corresponding to a letter match, mismatch, dash

Repeat lines 1 - 5 according to the number of analyses that will be run.

Example 
Input:
2
A,C,T,G
ACC
TTTACACGG
True
10 2 -4
A,C,T,G
ACC
TTTACACGG
False
10 2 -4

Output
Computing the Global Alingment of X = ACC and Y = TTTACACGG
The optimal global alignment score is 6
And the resulting aligned sequences are
X:---AC-C--
Y:TTTACACGG


Computing the Local Alingment of ACC and TTTACACGG
Given two sequences, X = ACC and Y = TTTACACGG
The optimal local alignment score is 26
And the resulting aligned sequences are
X:AC-C
Y:ACAC
'''

############### 
## Functions ##
###############

def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
	'''Given a set of letters in 'alphabet', and a scoring methodology defined by 
	values assigned to the respective input scoring parameters, build 
	a scoring matrix that will be used to compute the alignment scores of the two input sequences'''
	_alpha = alphabet | set(["-"])
	# create dictionary to store all individual letter dictionary values to assign when aligning letters
	_scoring_matrix = {}
	for _idx in _alpha:
		# create dictionary to store k,y pairs for each letter
		_letter_d = {}
		for _idy in _alpha:
			# score to assign when two dashes are aligned 
			if _idx == "-" and _idy == "-":
				_score = dash_score
			# score to assign when the letters match
			elif _idx == _idy:
				_score = diag_score
			# score to assign for a deletion (letter,-)
			elif _idx == "-" or _idy == "-":
				_score = dash_score
			# score to assign when the letters mismatch
			else:
				_score = off_diag_score
			# create entry in individual letter dictionary for value to assign when aligning matches
			_letter_d[_idy] = _score
		#print _letter_d
		_scoring_matrix[_idx] = _letter_d
	return _scoring_matrix

def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
	'''Given two gene sequences x,y compute the Dynmamic Programming table S
	which finds the optimal distance for two letters, one in each sequqnce,
	at positions i,j respsectively.  The flag is set to True for a DP table S
	that will be used to compute an optimal Global Alignment.  The flag
	is set to False if the DP table S will be used to compute an optimal 
	local alignment'''

	# insert a dash at the front of sequence x
	_x_prime = "-" + seq_x 
	# insert a dash at the front of sequence y
	_y_prime  = "-" + seq_y
	# populate an empty DP table with a number of rows equal to the length of seq_x + 1
	# the number of columns in each row will be equal to the length of seq_y + 1
	_dp_table = [[]]*len(_x_prime)
	_dp_table[0] = [0]
	# initialize the value at 0,0 in the DP Table to 0
	#_dp_table[0] = [0]
	#print _dp_table

	# Initiliaze values of aligning all letters in X and Y with a dash in DP table to the cumultive negative impact
	# This will create the Global Alignment DP Table to be used to create the Global Alignment
	if global_flag:			
		# initialize i=0 values to 0
		for _idx in range(len(_dp_table)):
			_dp_table[_idx] = [0]
		# initialize values for aligning all of the letters in X with a dash
		for _idx in range(1,len(_x_prime)):
			_xcurrent = _dp_table[_idx-1][0] + scoring_matrix[_x_prime[_idx]]["-"]
			_dp_table[_idx][0] = _xcurrent
		# initialize values for aligning all df the letters in Y with a dash
		for _idy in range(1,len(_y_prime)):
			_ycurrent = _dp_table[0][_idy-1] + scoring_matrix["-"][_y_prime[_idy]]
			_dp_table[0].append(_ycurrent)
		for _idx in range(1,len(_x_prime)):
		# iterate over each of the indices in seq Y, starting at index 1
			for _idy in range(1,len(_y_prime)):
			# determine the value of each of the three possible alignments that can be made: 
			# 1) letter-letter (diag for match or mismatch) 2) letter-dash (insert) 3) dash-letter (delete)
			# value of a letter-letter alignment made by moving diagonally
				_diag = _dp_table[_idx-1][_idy-1]+scoring_matrix[seq_x[_idx-1]][seq_y[_idy-1]]
				# value of a letter-letter alignment made by moving diagonally
				_x_dash = _dp_table[_idx-1][_idy]+scoring_matrix[seq_x[_idx-1]]["-"]
				# value of a letter-letter alignment made by moving diagonally
				_y_dash = _dp_table[_idx][_idy-1]+scoring_matrix["-"][seq_y[_idy-1]]
				# choose the optimal move, which is the maximal value of the three options
				_optimal_move = max(_diag,_x_dash,_y_dash)
				# insert this optimal value into the dp_table at the column index value idy in row idx 
				_dp_table[_idx].append(_optimal_move)
		return _dp_table
	# Initiliaze values of aligning all letters in X and Y with a dash in DP table to 0, essentially ignoring inserts/deletions
	# This will create the Local Alignment DP Table to be used to create the Local Alignment
	else:
		# initialize  an empty table with _idx rows, each with a single value of 0 in each
		for _idx in range(len(_dp_table)):
			_dp_table[_idx] = [0]
		# initialize the first entry in each row to 0
		for _idx in range(1,len(_x_prime)):
			_xcurrent = 0
			_dp_table[_idx][0] = _xcurrent
		# initialize the first entry in each column to 0 in the 0th row
		for _idy in range(1,len(_y_prime)):
			_ycurrent = 0
			_dp_table[0].append(_ycurrent)
	
		# iterate over each of the indices in seq X, starting at index 1 
		for _idx in range(1,len(_x_prime)):
			# iterate over each of the indices in seq Y, starting at index 1
			for _idy in range(1,len(_y_prime)):
				# determine the value of each of the three possible alignments that can be made: 
				# 1) letter-letter (diag for match or mismatch) 2) letter-dash (insert) 3) dash-letter (delete)
				
				# value of a letter-letter alignment made by moving diagonally
				_diag = _dp_table[_idx-1][_idy-1]+scoring_matrix[seq_x[_idx-1]][seq_y[_idy-1]]
				# value of a letter-letter alignment made by moving diagonally
				_x_dash = _dp_table[_idx-1][_idy]+scoring_matrix[seq_x[_idx-1]]["-"]
				# value of a letter-letter alignment made by moving diagonally
				_y_dash = _dp_table[_idx][_idy-1]+scoring_matrix["-"][seq_y[_idy-1]]
				# choose the optimal move, which is the maximal value of the three options
				_optimal_move = max(_diag,_x_dash,_y_dash,0)
				# insert this optimal value into the dp_table at the column index value idy in row idx 
				_dp_table[_idx].append(_optimal_move)
		return _dp_table

def compute_global_alignment(seq_x, seq_y,scoring_matrix,alignment_matrix):
	''' Given two input gene sequences, X and Y, compute the optimal GLOBAL alignement score and
	correpsonding optimally aligned sequences of both X and Y using the scoring matrix, M, and the
	Dynamic Programming Table S, alignment_matrix by iterating backwards through the alignment matrix'''
	
	# get index values of x and y in the dp table (alignment matrix) at the maximum score achieved
	_idx = len(seq_x)
	_idy = len(seq_y)

	# determine what this maximum score is by selecting the value at idx,idy in the table
	_score = alignment_matrix[_idx][_idy]
	
	# Initialize an empty global alignment string for both X and Y
	_align_x = ""
	_align_y = ""
	
	# check to make sure the dp table returned a meaningful value
	if _idx == 0 | _idy == 0:
		_score = 0

	# loop over all of the indices in X and Y until they are both equal to 0
	while _idx != 0 and _idy != 0:
		# if there is a match or a mismatch concatenate the letter at idx/y -1 to the respective alignments for x and y
		if alignment_matrix[_idx][_idy] == alignment_matrix[_idx-1][_idy-1] + scoring_matrix[seq_x[_idx-1]][seq_y[_idy-1]]:
			#_score = max(_score,alignment_matrix[_idx][_idy])
			_align_x = seq_x[_idx-1] + _align_x
			_align_y = seq_y[_idy-1] + _align_y
			# decrement index values for both x and y by one to check the next alignment
			_idx -= 1 
			_idy -= 1

		else:
			# if there is an insertion (letter from x with a dash instead of a ltter from y) 
			# then decrement only the x index value
			if alignment_matrix[_idx][_idy] == alignment_matrix[_idx-1][_idy] + scoring_matrix[seq_x[_idx-1]]["-"]:
			#	_score = max(_score, alignment_matrix[_idx][_idy])
				_align_x = seq_x[_idx-1] + _align_x
				_align_y = "-" + _align_y
				_idx -= 1
			# if there is an deltion (letter from y with a dash instead of a ltter from x) 
			# then decrement only the y index value 
			else:
			#	_score = max(_score, alignment_matrix[_idx][_idy])
				_align_x = "-" + _align_x
				_align_y = seq_y[_idy-1] + _align_y
				_idy -= 1
	# if the y seq letters are exhausted before the x sex letters, insert dashes in y and decrement x
	while _idx != 0:
		#_score = _score + scoring_matrix[seq_x[_idx-1]]["-"]
		_align_x = seq_x[_idx-1] + _align_x
		_align_y = "-" + _align_y
		_idx -= 1
	# if the x seq letters are exhausted before the x sex letters, insert dashes in y and decrement x
	while _idy != 0:
		#_score = _score + scoring_matrix[seq_y[_idy-1]]["-"]
		_align_x = "-" + _align_x
		_align_y = seq_y[_idy-1] + _align_y
		_idy -= 1
	return (_score, _align_x, _align_y)


def compute_local_alignment(seq_x,seq_y,scoring_matrix,alignment_matrix):
	''' Given two input gene sequences, seq_x and seq_y, compute the optimal LOCAL alignment score and
	correpsonding optimally aligned sequences of both X and Y using "scoring matrix" and the
	Dynamic Programming Table, "alignment_matrix"'''
	
	# initialize optimal score and corresponding i,j index values that obtain this maximal score
	_score = 0
	_best_i = 0
	_best_j = 0

	# iterate backwards through the DP table, row by row to find max value in the table
	for _row_idx in range(len((alignment_matrix))-1,0,-1):
		# create a variable to store the current score value by getting max value in row, _row_idx
		_current_score = max(alignment_matrix[_row_idx])
		# determine if existing score of newly observed score are higher 
		if _current_score > _score:
			# set score to max value, and grab  idx i and idx j associated with this value
			_score = _current_score
			_best_i = _row_idx
			_best_j = alignment_matrix[_row_idx].index(_current_score)
	
	_idx = _best_i
	_idy = _best_j
	
	#print "best row score", _score
	#print "i", _idx
	#print "j", _idy
	_align_x = ""
	_align_y = ""
	
	if _idx == 0 | _idy == 0:
		_score = 0
	
	while alignment_matrix[_idx][_idy] != 0:
		# if there is a match or a mismatch of letters associated with the score at the given position
		if alignment_matrix[_idx][_idy] == alignment_matrix[_idx-1][_idy-1] + scoring_matrix[seq_x[_idx-1]][seq_y[_idy-1]]:
			_align_x = seq_x[_idx-1] + _align_x
			_align_y = seq_y[_idy-1] + _align_y
			_idx -= 1 
			_idy -= 1
		else:
			# if there is an insertion
			if alignment_matrix[_idx][_idy] == alignment_matrix[_idx-1][_idy] + scoring_matrix[seq_x[_idx-1]]["-"]:
				_align_x = seq_x[_idx-1] + _align_x
				_align_y = "-" + _align_y
				_idx -= 1
			# if there is a deletion
			else:
				_align_x = "-" + _align_x
				_align_y = seq_y[_idy-1] + _align_y
				_idy -= 1
	return (_score, _align_x, _align_y)

if __name__ == "__main__":
	import sys
	inFile = sys.argv[1]
#	outFile = sys.argv[2]
	lines = ""
	with open(inFile,'r') as i:
		lines = i.readlines()
	# alphabet to use when computing the optimal alignments
	entries = int(lines[0])
	i = 0
	for entry in range(entries):
		_alpha = set(lines[i+1].strip().split(","))
		# X gene sequence to analyze
		_x = lines[i+2].strip().split(" ")[0]
		# Y gene sequence to analyze
		_y = lines[i+3].strip().split(" ")[0]
		# True for Global Alignment
		# False for Local Alignment
		_flag = lines[i+4].strip().split(" ")[0]
		#print _flag
		if _flag:
			if _flag == "True":
				_flag = True
			else:
				_flag = False

		# assign the values for a match, mismatch, insertion and deletion between two sequences	
		_match,_mismatch, _dash = map(int, lines[i+5].strip().split(" "))
		# increment 5 lines to read in the next round of gene sequences to analyze
		i+=5
		#print "x:{}, y:{}".format(_x,_y)
		#print "match:{}, mismatch:{}, insert:{}, delete:{}".format(_match, _mismatch,_dash)#
		#print _alpha

		_M_ = build_scoring_matrix(_alpha, _match, _mismatch, _dash)
		#print "Scoring Matrix"
		#for key, value in _M_.items():
		#	print "Resulting values when pairing with letter {} = {}".format(key,value)
		_S_  = compute_alignment_matrix(_x,_y,_M_,_flag)
		#print type(_flag), _flag
		if _flag:
		# Output the results of the analysis to the console
		#	print type(_flag), _flag
			print ""
			print "Computing the Global Alingment of X = {} and Y = {}".format(_x,_y)
			ans_global = compute_global_alignment(_x,_y,_M_,_S_)
			#print "Given two sequences, X = {} and Y = {}".format(_x,_y)
			print "The optimal global alignment score is {}".format(str(ans_global[0]))
			print "And the resulting aligned sequences are"
			print "X:{}".format(ans_global[1])
			print "Y:{}".format(ans_global[2])
			# Calculate the Edit Distance for word comparisons
			if len(_M_.keys()) == 27:
				print "The Edit Distance between these two words = {}".format(len(_x) + len(_y) - ans_global[0])
			print ""
		else:
			print ""
			print "Computing the Local Alingment of {} and {}".format(_x,_y)
			ans_local = compute_local_alignment(_x,_y,_M_,_S_)
			print "Given two sequences, X = {} and Y = {}".format(_x,_y)
			print "The optimal local alignment score is {}".format(str(ans_local[0]))
			print "And the resulting aligned sequences are"
			print "X:{}".format(ans_local[1])
			print "Y:{}".format(ans_local[2])
			print ""