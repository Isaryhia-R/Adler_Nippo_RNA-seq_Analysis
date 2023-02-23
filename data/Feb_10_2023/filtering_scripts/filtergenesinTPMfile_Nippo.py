import sys

infile = open(sys.argv[1])
outfile = open(sys.argv[2], 'w')

for line in infile:
	if line[0:4] == "gene":
		outfile.write(line)
	else:
		f = line.strip().split('\t')
		
		if float(f[1]) >= 1.0 and float(f[2]) >= 1.0 and float(f[3]) >= 1.0:
			outfile.write(line)
		elif float(f[4]) >= 1.0 and float(f[5]) >= 1.0 and float(f[6]) >= 1.0:
			outfile.write(line)
		elif float(f[7]) >= 1.0 and float(f[8]) >= 1.0 and float(f[9]) >= 1.0:
			outfile.write(line)
		elif float(f[10]) >= 1.0 and float(f[11]) >= 1.0 and float(f[12]) >= 1.0:
			outfile.write(line)
		elif float(f[13]) >= 1.0 and float(f[14]) >= 1.0 and float(f[15]) >= 1.0:
			outfile.write(line)
		elif float(f[16]) >= 1.0 and float(f[17]) >= 1.0 and float(f[18]) >= 1.0:
			outfile.write(line)
		elif float(f[19]) >= 1.0 and float(f[20]) >= 1.0 and float(f[21]) >= 1.0:
			outfile.write(line)
		elif float(f[22]) >= 1.0 and float(f[23]) >= 1.0 and float(f[24]) >= 1.0:
			outfile.write(line)
		elif float(f[25]) >= 1.0 and float(f[26]) >= 1.0 and float(f[27]) >= 1.0:
			outfile.write(line)
		elif float(f[28]) >= 1.0 and float(f[29]) >= 1.0 and float(f[30]) >= 1.0:
			outfile.write(line)
		elif float(f[31]) >= 1.0 and float(f[32]) >= 1.0 and float(f[33]) >= 1.0:
			outfile.write(line)
		elif float(f[34]) >= 1.0 and float(f[35]) >= 1.0 and float(f[36]) >= 1.0:
			outfile.write(line)
		elif float(f[37]) >= 1.0 and float(f[38]) >= 1.0 and float(f[39]) >= 1.0:
			outfile.write(line)
		elif float(f[40]) >= 1.0 and float(f[41]) >= 1.0 and float(f[42]) >= 1.0:
			outfile.write(line)
		elif float(f[43]) >= 1.0 and float(f[44]) >= 1.0 and float(f[45]) >= 1.0:
			outfile.write(line)
		elif float(f[46]) >= 1.0 and float(f[47]) >= 1.0 and float(f[48]) >= 1.0:
			outfile.write(line)
		elif float(f[49]) >= 1.0 and float(f[50]) >= 1.0 and float(f[51]) >= 1.0:
			outfile.write(line)
		elif float(f[52]) >= 1.0 and float(f[53]) >= 1.0 and float(f[54]) >= 1.0:
			outfile.write(line)
		elif float(f[55]) >= 1.0 and float(f[56]) >= 1.0 and float(f[57]) >= 1.0:
			outfile.write(line)
		elif float(f[58]) >= 1.0 and float(f[59]) >= 1.0 and float(f[60]) >= 1.0:
			outfile.write(line)
		elif float(f[61]) >= 1.0 and float(f[62]) >= 1.0 and float(f[63]) >= 1.0:
			outfile.write(line)
		elif float(f[64]) >= 1.0 and float(f[65]) >= 1.0 and float(f[66]) >= 1.0:
			outfile.write(line)
		elif float(f[67]) >= 1.0 and float(f[68]) >= 1.0 and float(f[69]) >= 1.0:
			outfile.write(line)
		elif float(f[70]) >= 1.0 and float(f[71]) >= 1.0 and float(f[72]) >= 1.0:
			outfile.write(line)
		elif float(f[73]) >= 1.0 and float(f[74]) >= 1.0 and float(f[75]) >= 1.0:
			outfile.write(line)
		elif float(f[76]) >= 1.0 and float(f[77]) >= 1.0 and float(f[78]) >= 1.0:
			outfile.write(line)
		elif float(f[79]) >= 1.0 and float(f[80]) >= 1.0 and float(f[81]) >= 1.0:
			outfile.write(line)
		elif float(f[82]) >= 1.0 and float(f[83]) >= 1.0 and float(f[84]) >= 1.0:
			outfile.write(line)
		elif float(f[85]) >= 1.0 and float(f[86]) >= 1.0 and float(f[87]) >= 1.0:
			outfile.write(line)
		elif float(f[88]) >= 1.0 and float(f[89]) >= 1.0 and float(f[90]) >= 1.0:
			outfile.write(line)
		elif float(f[91]) >= 1.0 and float(f[92]) >= 1.0 and float(f[93]) >= 1.0:
			outfile.write(line)
		elif float(f[94]) >= 1.0 and float(f[95]) >= 1.0 and float(f[96]) >= 1.0:
			outfile.write(line)
		elif float(f[97]) >= 1.0 and float(f[98]) >= 1.0 and float(f[99]) >= 1.0:
			outfile.write(line)
		elif float(f[100]) >= 1.0 and float(f[101]) >= 1.0 and float(f[102]) >= 1.0:
			outfile.write(line)
		elif float(f[103]) >= 1.0 and float(f[104]) >= 1.0 and float(f[105]) >= 1.0:
			outfile.write(line)
		elif float(f[106]) >= 1.0 and float(f[107]) >= 1.0 and float(f[108]) >= 1.0:
			outfile.write(line)
		elif float(f[109]) >= 1.0 and float(f[110]) >= 1.0 and float(f[111]) >= 1.0:
			outfile.write(line)
		elif float(f[112]) >= 1.0 and float(f[113]) >= 1.0 and float(f[114]) >= 1.0:
			outfile.write(line)
		elif float(f[115]) >= 1.0 and float(f[116]) >= 1.0 and float(f[117]) >= 1.0:
			outfile.write(line)
		elif float(f[118]) >= 1.0 and float(f[119]) >= 1.0 and float(f[120]) >= 1.0:
			outfile.write(line)
		elif float(f[121]) >= 1.0 and float(f[122]) >= 1.0 and float(f[123]) >= 1.0:
			outfile.write(line)
		elif float(f[124]) >= 1.0 and float(f[125]) >= 1.0 and float(f[126]) >= 1.0:
			outfile.write(line)
		elif float(f[127]) >= 1.0 and float(f[128]) >= 1.0 and float(f[129]) >= 1.0:
			outfile.write(line)
		elif float(f[130]) >= 1.0 and float(f[131]) >= 1.0 and float(f[132]) >= 1.0:
			outfile.write(line)
		elif float(f[133]) >= 1.0 and float(f[134]) >= 1.0 and float(f[135]) >= 1.0:
			outfile.write(line)
		elif float(f[136]) >= 1.0 and float(f[137]) >= 1.0 and float(f[138]) >= 1.0:
			outfile.write(line)
		elif float(f[139]) >= 1.0 and float(f[140]) >= 1.0 and float(f[141]) >= 1.0:
			outfile.write(line)
		elif float(f[142]) >= 1.0 and float(f[143]) >= 1.0 and float(f[144]) >= 1.0:
			outfile.write(line)
		elif float(f[145]) >= 1.0 and float(f[146]) >= 1.0 and float(f[147]) >= 1.0:
			outfile.write(line)
		elif float(f[148]) >= 1.0 and float(f[148]) >= 1.0 and float(f[149]) >= 1.0:
			outfile.write(line)
		elif float(f[150]) >= 1.0 and float(f[151]) >= 1.0 and float(f[152]) >= 1.0:
			outfile.write(line)
		elif float(f[153]) >= 1.0 and float(f[154]) >= 1.0 and float(f[155]) >= 1.0:
			outfile.write(line)
		elif float(f[156]) >= 1.0 and float(f[157]) >= 1.0 and float(f[158]) >= 1.0:
			outfile.write(line)
		elif float(f[159]) >= 1.0 and float(f[160]) >= 1.0 and float(f[161]) >= 1.0:
			outfile.write(line)
		elif float(f[162]) >= 1.0 and float(f[163]) >= 1.0 and float(f[164]) >= 1.0:
			outfile.write(line)
		elif float(f[165]) >= 1.0 and float(f[166]) >= 1.0 and float(f[167]) >= 1.0:
			outfile.write(line)
		elif float(f[168]) >= 1.0 and float(f[169]) >= 1.0 and float(f[170]) >= 1.0:
			outfile.write(line)
		elif float(f[171]) >= 1.0 and float(f[172]) >= 1.0 and float(f[173]) >= 1.0:
			outfile.write(line)
		elif float(f[174]) >= 1.0 and float(f[175]) >= 1.0 and float(f[176]) >= 1.0:
			outfile.write(line)
		elif float(f[177]) >= 1.0 and float(f[178]) >= 1.0 and float(f[179]) >= 1.0:
			outfile.write(line)
		elif float(f[180]) >= 1.0 and float(f[181]) >= 1.0 and float(f[182]) >= 1.0:
			outfile.write(line)
		elif float(f[183]) >= 1.0 and float(f[184]) >= 1.0 and float(f[185]) >= 1.0:
			outfile.write(line)
		elif float(f[186]) >= 1.0 and float(f[187]) >= 1.0 and float(f[188]) >= 1.0:
			outfile.write(line)
		elif float(f[189]) >= 1.0 and float(f[190]) >= 1.0 and float(f[191]) >= 1.0:
			outfile.write(line)
		elif float(f[192]) >= 1.0 and float(f[193]) >= 1.0 and float(f[194]) >= 1.0:
			outfile.write(line)
		elif float(f[195]) >= 1.0 and float(f[196]) >= 1.0 and float(f[197]) >= 1.0:
			outfile.write(line)
		elif float(f[198]) >= 1.0 and float(f[199]) >= 1.0 and float(f[200]) >= 1.0:
			outfile.write(line)
		elif float(f[201]) >= 1.0 and float(f[202]) >= 1.0 and float(f[203]) >= 1.0:
			outfile.write(line)
		elif float(f[204]) >= 1.0 and float(f[205]) >= 1.0 and float(f[206]) >= 1.0:
			outfile.write(line)
		elif float(f[207]) >= 1.0 and float(f[208]) >= 1.0 and float(f[209]) >= 1.0:
			outfile.write(line)
		elif float(f[210]) >= 1.0 and float(f[211]) >= 1.0 and float(f[212]) >= 1.0:
			outfile.write(line)
		elif float(f[213]) >= 1.0:
			outfile.write(line)
infile.close()
outfile.close()
