import os

cases = [
	'GAW1',
	'GRK1',
	'n0012',
	'n662415',
	'NLF0215',
	'SC2110'
]

for c in cases:
	os.chdir(c)

	os.system('julia %s.jl' % (c,))

	os.chdir('..')
