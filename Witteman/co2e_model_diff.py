def Dab(species1, species2, constants, P, T):

	part1 = constants['a']*(T/(species1['Tc_a']*species2['Tc_b'])**0.5)**constants['b']
	part2 = (species1['Pc_a']*species2['Pc_b'])**(1/3)
	part3 = (species1['Tc_a']*species2['Tc_b'])**(5/12)
	part4 = (1/species1['Mw_a']+1/species2['Mw_b'])**0.5 
	Diff = part1*part2*part3*part4/P
	

	Diff = Diff/1e4

	return Diff