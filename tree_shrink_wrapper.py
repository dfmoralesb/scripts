"""
Detects and trim abnormally long branches using TreeShrink

TreeShrink must be install and on path

"""

import sys, os, shutil

def trim(DIR,tree_file_ending):
	if DIR[-1] != "/": DIR += "/"
	
	filecount = 0
	
	#runs treeshrink
	for i in os.listdir(DIR):
		if i.endswith(tree_file_ending):
			print i
			filecount += 1
			cmd= ["treeshrink.py","-i", DIR+i ,"-o", i+".tt", "-c","-m per-gene", "-d",DIR+i+".tt_dir"]
			print (" ".join(cmd))
			os.system(" ".join(cmd))
	
	#moves output files to DIR and delete treeshrink individual folders
	for j in os.listdir(DIR):
		if j.endswith(".tt_dir"):
			source = DIR+j
			dest = DIR
			files = os.listdir(source)
			for f in files:
				shutil.move(source+"/"+f, dest)
			shutil.rmtree(DIR+j)
	
	#removes single quotes from tip labels
	for k in os.listdir(DIR):
		if k.endswith(".tt"):
			with open(DIR+k, 'r+') as f:
				content = f.read()
				f.seek(0)
				f.truncate()
				f.write(content.replace("'", ""))
			f.close()
	assert filecount > 0, \
		"No file end with "+tree_file_ending+" found in "+DIR
			
			
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "python tree_shrink_wrapper.py DIR tree_file_ending"
		sys.exit(0)

	DIR,tree_file_ending = sys.argv[1:]
	trim(DIR,tree_file_ending)



