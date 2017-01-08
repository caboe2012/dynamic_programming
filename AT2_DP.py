''' WTF!!!! '''
import string
# This Python file uses the following encoding: utf-8
#'''seriously?'''
letters = set(["A","C","T","G"])
diag_score1 = 10
off_diag_score1 = 6
dash_score1 = -3
def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
	'''seriously?'''
	_alpha = alphabet | set(["-"])
	_beta = alphabet | set(["-"])
	#print _alpha
	_final_d = {}
	for _idx in _alpha:
		_temp_d = {}
		for _idy in _beta:
			if _idx == "-" and _idy == "-":
				_score = dash_score
			elif _idx == _idy:
				_score = diag_score
			elif _idx == "-" or _idy == "-":
				_score = dash_score
			else:
				_score = off_diag_score
		#	print _idx,_idy, _score
			_temp_d[_idy] = _score
			_final_d[_idx] = _temp_d
	return _final_d
#_M_ = build_scoring_matrix(set(["A","C","T","G"]), 6, 2, -4)
#print result


#X = 'A'
#Y = 'A'
#result = {'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2}, 'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2}, '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4}, 'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2}, 'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}}
#flag = True

#expected [[0, -4], [-4, 6]] but received [[0]]

def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
	'''AGGGHHHHHHHHH'''
#	if len(seq_y) == 0 or len(seq_y) == 0:
#		return [[0]]
	_ecks = "-" + seq_x 
	_why  = "-" + seq_y
	_es = [[]]*len(_ecks)
	_es[0] = [0]
	if global_flag:	
		for _idx in range(len(_es)):
			_es[_idx] = [0]
		for _idx in range(1,len(_ecks)):
			_xcurrent = _es[_idx-1][0] + scoring_matrix[_ecks[_idx]]["-"]
			_es[_idx][0] = _xcurrent
		for _idy in range(1,len(_why)):
			_ycurrent = _es[0][_idy-1] + scoring_matrix["-"][_why[_idy]]
			_es[0].append(_ycurrent)
		for _idx in range(1,len(_ecks)):
			for _idy in range(1,len(_why)):
				_diag = _es[_idx-1][_idy-1]+scoring_matrix[seq_x[_idx-1]][seq_y[_idy-1]]
				_x_dash = _es[_idx-1][_idy]+scoring_matrix[seq_x[_idx-1]]["-"]
				_y_dash = _es[_idx][_idy-1]+scoring_matrix["-"][seq_y[_idy-1]]
				_winner = max(_diag,_x_dash,_y_dash)
				_es[_idx].append(_winner)
		return _es
	else:
		for _idx in range(len(_es)):
			_es[_idx] = [0]
		for _idx in range(1,len(_ecks)):
			_xcurrent = 0
			_es[_idx][0] = _xcurrent
		for _idy in range(1,len(_why)):
			_ycurrent = 0
			_es[0].append(_ycurrent)
		for _idx in range(1,len(_ecks)):
			for _idy in range(1,len(_why)):
				_diag = _es[_idx-1][_idy-1]+scoring_matrix[seq_x[_idx-1]][seq_y[_idy-1]]
				_x_dash = _es[_idx-1][_idy]+scoring_matrix[seq_x[_idx-1]]["-"]
				_y_dash = _es[_idx][_idy-1]+scoring_matrix["-"][seq_y[_idy-1]]
				_winner = max(_diag,_x_dash,_y_dash,0)
				_es[_idx].append(_winner)
		return _es


def compute_global_alignment(seq_x, seq_y,scoring_matrix,alignment_matrix):
	''' give me my points'''
	_idx = len(seq_x)
	_idy = len(seq_y)
	_align_x = ""
	_align_y = ""
	_score = alignment_matrix[_idx][_idy]
	if _idx == 0 | _idy == 0:
		_score = 0
	while _idx != 0 and _idy != 0:
#		print "i",_idx
#		print "j", _idy 
		if alignment_matrix[_idx][_idy] == alignment_matrix[_idx-1][_idy-1] + scoring_matrix[seq_x[_idx-1]][seq_y[_idy-1]]:
			#_score = max(_score,alignment_matrix[_idx][_idy])
			_align_x = seq_x[_idx-1] + _align_x
			_align_y = seq_y[_idy-1] + _align_y
			_idx -= 1 
			_idy -= 1

		else:
			if alignment_matrix[_idx][_idy] == alignment_matrix[_idx-1][_idy] + scoring_matrix[seq_x[_idx-1]]["-"]:
			#	_score = max(_score, alignment_matrix[_idx][_idy])
				_align_x = seq_x[_idx-1] + _align_x
				_align_y = "-" + _align_y
				_idx -= 1
			else:
			#	_score = max(_score, alignment_matrix[_idx][_idy])
				_align_x = "-" + _align_x
				_align_y = seq_y[_idy-1] + _align_y
				_idy -= 1
	while _idx != 0:
		#_score = _score + scoring_matrix[seq_x[_idx-1]]["-"]
		_align_x = seq_x[_idx-1] + _align_x
		_align_y = "-" + _align_y
		_idx -= 1
	while _idy != 0:
		#_score = _score + scoring_matrix[seq_y[_idy-1]]["-"]
		_align_x = "-" + _align_x
		_align_y = seq_y[_idy-1] + _align_y
		_idy -= 1
	return (_score, _align_x, _align_y)

def compute_local_alignment(seq_x,seq_y,scoring_matrix,alignment_matrix):
	''' give me my points'''
	_score = 0
	_best_i = 0
	_best_j = 0
	for _row_idx in range(len((alignment_matrix))-1,0,-1):
		_current_score = max(alignment_matrix[_row_idx])
		if _current_score > _score:
			_score = _current_score
			_best_i = _row_idx
			_best_j = alignment_matrix[_row_idx].index(_current_score)
	_idx = _best_i
	_idy = _best_j
#	print "i", _idx
#	print "j", _idy
	_align_x = ""
	_align_y = ""
	if _idx == 0 | _idy == 0:
		_score = 0
	while alignment_matrix[_idx][_idy] != 0:
		if alignment_matrix[_idx][_idy] == alignment_matrix[_idx-1][_idy-1] + scoring_matrix[seq_x[_idx-1]][seq_y[_idy-1]]:
			_align_x = seq_x[_idx-1] + _align_x
			_align_y = seq_y[_idy-1] + _align_y
			_idx -= 1 
			_idy -= 1
		else:
			if alignment_matrix[_idx][_idy] == alignment_matrix[_idx-1][_idy] + scoring_matrix[seq_x[_idx-1]]["-"]:
				_align_x = seq_x[_idx-1] + _align_x
				_align_y = "-" + _align_y
				_idx -= 1
			else:
				_align_x = "-" + _align_x
				_align_y = seq_y[_idy-1] + _align_y
				_idy -= 1
	return (_score, _align_x, _align_y)






