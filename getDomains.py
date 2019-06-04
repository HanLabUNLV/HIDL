import requests, sys, os, re

args = sys.argv
ids = []

try:
	os.mkdir("/home/dbarth/Work/GeneTrees/geneDomains/" + args[1])
except OSError:
	pass

with open("GeneTrees/geneNewicks/" + args[1] + ".newick") as f:
	string = f.readline()
	try:
		ids = re.findall("ENSP0(.+?):",string)
	except AttributeError:
		ids = []
		sys.exit("NO HUMAN PROTEINS IN: "+args[1])
	for i in ids:
		proteinid = "ENSP0"+i
		server = "http://rest.ensembl.org"
		ext = "/overlap/translation/" + proteinid + "?type=Superfamily"

		r = requests.get(server+ext, headers={"Content-Type" : "application/json"})

		if not r.ok:
			r.raise_for_status()
			sys.exit("not r.ok, ID: " + args[1])

		decoded = r.json()
		if decoded != []:
			with open("/home/dbarth/Work/GeneTrees/geneDomains/" + args[1] + "/" + proteinid + ".json","w") as g:
				g.write(repr(decoded))
				g.close()
		else:
			os.rmdir("/home/dbarth/Work/GeneTrees/geneDomains/" + args[1])
