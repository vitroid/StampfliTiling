%.cage:
	genice2 cockayne[$*] -f cage[python] | awk '/^1/&&(NF==4){c[$$1]++}END{for(s in c){print s,c[s]}}'

#test4
#12 136
#14 16
#15 16
#16 56
pep8:
	autopep8 -r -a -a -i ./
