#!/usr/bin/python3
# -*- coding: utf-8 -*

# Libraries import
import random
import math 			
import matplotlib.pyplot as plt

# Functions
def alphabet(type):
	"""
	Return the alphabet (as a string) corresponding to the type passed in argument
	"""
	if type == "nucleic":
		return "ACGT"
	elif type == "proteic":
		return "ACDEFGHIKLMNPQRSTVWY"
	elif type == "iupac_nucleic":
		return "ACGTUWSMKRYBDHVN"
	elif type == "iupac_proteic":
		return "ARNDCQEGHILKMFPSTWYVX"
	else:
		print("Error : False type.\n")
		return


def randseq(num, alpha):
	"""
	Generate a random sequence of 'num' length
	"""
	# Convert the string alphabet to a list
	alpha_list = list(alpha)	
	seq = str()
	for i in range(num):
		# Get a random value between 0 and the length of the alphabet
		value = random.randint(0, len(alpha_list)-1)
		# Add to the sequence the nucleotide of the alphabet present at the index 'value'
		seq += alpha_list[value]
	return seq


def hamming(seq1, seq2):
	"""
	Calculate the hamming distance between two sequences passed in arguments
	"""
	hamming_distance = 0
	# If sequences have not the same size, the calcul of hammind distance is impossible
	if len(seq1) != len(seq2):
		print("Error : Length of the two sequences are not the same.\n")
	# Convertion of the two sequence string to lists
	seq1_list = list(seq1)
	seq2_list = list(seq2)
	# For each nucleotides, verification of the equality between the two sequences
	for i in range(len(seq1)):
		# If no equality, put +1 to the hamming distance
		if seq1_list[i] != seq2_list[i]:
			hamming_distance += 1
	return hamming_distance


def mutate(seq, num_subs):		
	"""
	Create a mutated sequence from a sequence passed in argument, with 'num_subs' substitutions
	Calculate and return the hamming distance between the sequence passed in argument 
	and the mutated sequence created.
	"""
	mutated_seq = str()
	seq_list = list(seq)
	# Do 'num_subs' substitutions
	for i in range(num_subs):
		# Get a random value
		value = random.randint(0, len(seq_list)-1)
		# Get a random nucleotide
		single_nucleotide = randseq(1, alphabet("nucleic"))
		# Convert the initial nucleotide at the place 'value' to the random nucleotide
		seq_list[value] = single_nucleotide
	# Convert mutatated sequence list to a string
	for nucleotid in seq_list:
		mutated_seq += nucleotid
	return hamming(seq, mutated_seq) 	


def experiments(le,su,nb):
	"""
	Make nb times the experiment and return each hamming distances in a list
	"""
	# Repeat the experiments `R` times
	v = []
	for i in range(nb):
		seq = randseq(le, alphabet('nucleic'))
		v.append(mutate(seq,su))
	return v


def variance(data):
	"""
	Calculate and return the variance of a list with the Welford's method
	"""
	M = 0
	S = 0
	N = len(data)
	for i in range(N):
		#Get the value i
		x = data[i]
		# Add the mean to the old mean (i-1)
		oldM = M
		# Add to the global mean, the mean calculated for the value i
		M += (x-M)/(i+1)
		# Add to the global variance, the variance calculated for the value i
		S += (x-M)*(x-oldM)
	return S/(N-1)


def mean(data):
	"""
	Calculate and return the mean of data present in a list passed in argument
	"""
	return sum(data)/len(data)

def std(data):					
	"""
	Calculate and return the standard deviation of data present in a list passed in argument
	"""
	# Calculation of the square root of the variance (according to the Welford's method)
	# and return it
	return math.sqrt(variance(data))



def distanceJC69(means, L): 		
	"""
	Calculate and return the JC69 distance with a mean of Hamming distance passed in argument
	"""
	# Calculation of the p-distance
	p = means/L
	# Calcultation of the JC69 distance thanks to the equation of Jukes and Cantor (1969).
	d = (-3/4)*math.log(1.0-(p*(4/3))) 
	return d


def generate(le, nb, xx): 
	""""
	Make 'nb' times the experiment and return each hamming distances and JC69 distance in respective list
	"""
	hamming_distances_list = list()
	JC69_distances = list()
	hamming_distances = list()
	# For each substitution number present in the 'xx' list passed in arguemnt
	for nb_substitution in xx:
		# Calculation of the Hamming distance for 'nb__substitution', 'nb' times
		# Put all values in a list
		hamming_distances_list = experiments(le, nb_substitution, nb)
		# Calculation of the mean of the Hamming distances and put it in a list
		# The list contains Hamming distance mean for each 'nb_substitution'
		hamming_distances.append(mean(hamming_distances_list))		
	# Calculation of the JC9 distances for each Hamming distance means
	for m in hamming_distances:
		# Append each JC69 distances calculate in a list
		JC69_distances.append(distanceJC69(m, le))
	return [hamming_distances, JC69_distances]
