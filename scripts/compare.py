#! usr/bin/python
import sys, os, string, re, time, stat

# Check for windows
isWindows = os.name =='nt'

default_css = \
"""
TABLE { border-collapse: collapse; border-spacing: 0px; }
TD.linenum { color: #909090; 
             text-align: right;
             vertical-align: top;
             font-weight: bold;
             border-right: 1px solid black;
             border-left: 1px solid black; }
TD.added { background-color: #DDDDFF; }
TD.modified { background-color: #BBFFBB; }
TD.removed { background-color: #FFCCCC; }
TD.normal { background-color: #FFFFE1; }
"""

def str2html(s) :
    s1 = string.replace(string.rstrip(s), "&", "&amp;")
    if ( len(s1) == 0 ) : return ( s1 ) ;
    s1 = string.replace(s1, "<", "&lt;")
    s1 = string.replace(s1, ">", "&gt;")
    i = 0
    s2 = ""
    while ( s1[i] == " " ) :
        s2 += "&nbsp;"
        i += 1
    s2 += s1[i:]
    return ( s2 )
    

if len(sys.argv) < 4:
	print 'Usage: compare [List File] [Reference Benchmark Directory] [Output Filename] [Optional: Benchmark Directory]'
	sys.exit(1)
		
listFilename = sys.argv[1]
referenceBenchmarkDir = sys.argv[2]
outputFilename = sys.argv[3]
only_changes = 1
external_css = ""
retCode = 0
benchmarkDir = ""
#resultsDir = sys.argv[2]
if len(sys.argv) > 4:
	benchmarkDir = sys.argv[4]

outputFile = open(outputFilename, 'w')

print >>outputFile, """
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN"
 "http://www.w3.org/TR/REC-html40/loose.dtd">
<html>
<head>
    <title>Benchmark results</title>
    <style>%s</style>
</head>
<body>
<h1>Benchmark Results</h1>
<table width="100%%">
	<tr>
		<th width="50%%">Reference benchmark files from %s</th>
		<th width="50%%">New benchmark files</th>
	</tr>
</table>
""" % (default_css, referenceBenchmarkDir)

notExistingFiles = """
	<h1>Missing files</h1>
	<span style="color:red">"""

identicalFiles = """
	<h1>Identical files</h1>
	<span style="color:green">"""

if (os.path.isfile(listFilename)):
	listFile = open(listFilename, 'r')
else:
	print >>outputFile, "List file not found"
	print "List file not found"
	sys.exit(1)
	
for line in listFile:
	if not line.strip():
		continue
	elif line[0] == '#' and line[1] == '#':
		print >>outputFile, """
<h2 style="color:black">%s</h2>""" % line[3:]
		notExistingFiles = notExistingFiles + """
<h2 style="color:black">%s</h2>""" % line[3:]
		identicalFiles= identicalFiles + """
<h2 style="color:black">%s</h2>""" % line[3:]
	# Ignore comments
	elif line[0] == '#':
		continue
	# Ignore CMake variable stuff
	elif (line[0:3] == 'SET'):
		continue
	elif (line[0] == ')'):
		continue
	else:
		line = line.rstrip() # Removes line endings
		file1 = os.path.normpath(referenceBenchmarkDir + line);
		file2 = './' + line
		file2 = os.path.normpath(benchmarkDir + line)
		
		if not os.path.isfile(file1):
			notExistingFiles = notExistingFiles + """
	Reference File %s<br/>""" % file1
			retCode = 1
			continue

		if not os.path.isfile(file2):
			notExistingFiles = notExistingFiles + """
	%s<br/>""" % file2
			retCode = 1
			continue
		
		#line = line.replace('/', '_')
		#resultPath = os.path.normpath(resultsDir + line)
		#retCode = subprocess.Popen('python diff2html.py --only-changes ' + file1 + ' ' + file2 + ' > ' + resultPath + '.html', shell=True )
		
		# Invokes "diff"
		if isWindows:
			diffCommand = '../diff.exe';
		else:
			diffCommand ='diff'
		diff_stdout = os.popen("%s --ignore-all-space %s %s" % (diffCommand, file1.replace('\\', '/'), file2.replace('\\', '/')), "r")
		diff_output = diff_stdout.readlines()
		diff_stdout.close()
	    # Maps to store the reported differences
		changed = {}
		deleted = {}
		added = {}
	    # Magic regular expression
		diff_re = re.compile(
			r"^(?P<f1_start>\d+)(,(?P<f1_end>\d+))?"+ \
			 "(?P<diff>[acd])"+ \
			 "(?P<f2_start>\d+)(,(?P<f2_end>\d+))?")
	    # Now parse the output from "diff"
		for diff_line in diff_output:
			diffs = diff_re.match(string.strip(diff_line))
	        # If the line doesn't match, it's useless for us
			if not ( diffs  == None ) :
				# Retrieving informations about the differences : 
				# starting and ending lines (may be the same)
				f1_start = int(diffs.group("f1_start"))
				if ( diffs.group("f1_end") == None ) :
					f1_end = f1_start
				else :
					f1_end = int(diffs.group("f1_end"))
				f2_start = int(diffs.group("f2_start"))
				if ( diffs.group("f2_end") == None ) :
					f2_end = f2_start
				else :
					f2_end = int(diffs.group("f2_end"))
				f1_nb = (f1_end - f1_start) + 1
				f2_nb = (f2_end - f2_start) + 1
				# Is it a changed (modified) line ?
				if ( diffs.group("diff") == "c" ) :
					# We have to handle the way "diff" reports lines merged
					# or splitted
					if ( f2_nb < f1_nb ) :
						# Lines merged : missing lines are marqued "deleted"
						for lf1 in range(f1_start, f1_start+f2_nb) :
							changed[lf1] = 0
						for lf1 in range(f1_start+f2_nb, f1_end+1) :
							deleted[lf1] = 0
					elif ( f1_nb < f2_nb ) :
						# Lines splitted : extra lines are marqued "added"
						for lf1 in range(f1_start, f1_end+1) :
							changed[lf1] = 0
						for lf2 in range(f2_start+f1_nb, f2_end+1) :
							added[lf2] = 0
					else :
						# Lines simply modified !
						for lf1 in range(f1_start, f1_end+1) :
							changed[lf1] = 0
	            # Is it an added line ?
				elif ( diffs.group("diff") == "a" ) :
					for lf2 in range(f2_start, f2_end+1):
						added[lf2] = 0
				else :
	            # OK, so it's a deleted line
					for lf1 in range(f1_start, f1_end+1) :
						deleted[lf1] = 0
	
		# Storing the two compared files, to produce the HTML output
		f1 = open(file1, "r")
		f1_lines = f1.readlines()
		f1.close()
		f2 = open(file2, "r")
		f2_lines = f2.readlines()
		f2.close()
		
		# Finding some infos about the files
		f1_stat = os.stat(file1)
		f2_stat = os.stat(file2)
		
		# Printing the HTML header, and various known informations
		
		# Preparing the links to changes
		if ( len(changed) == 0 ) :
			changed_lnks = "None"
		else :
			changed_lnks = ""
			keys = changed.keys()
			keys.sort()
			for key in keys :
				changed_lnks += "<a href=\"#F1_%d\">%d</a>, " % (key, key)
			changed_lnks = changed_lnks[:-2]
	    
		if ( len(added) == 0 ) :
			added_lnks = "None"
		else :
			added_lnks = ""
			keys = added.keys()
			keys.sort()
			for key in keys :
				added_lnks += "<a href=\"#F2_%d\">%d</a>, " % (key, key)
			added_lnks = added_lnks[:-2]
	
		if ( len(deleted) == 0 ) :
			deleted_lnks = "None"
		else :
			deleted_lnks = ""
			keys = deleted.keys()
			keys.sort()
			for key in keys :
				deleted_lnks += "<a href=\"#F1_%d\">%d</a>, " % (key, key)
			deleted_lnks = deleted_lnks[:-2]
	        
		if (changed_lnks == 'None' and added_lnks == 'None' and deleted_lnks == 'None'):
			identicalFiles = identicalFiles + """
	%s<br/>""" % line
			continue
		
		retCode = 1
	        
		print >>outputFile, """
	<h3>%s</h3>""" % line
	
#		print >>outputFile, """	
#	<table>
#	    <tr>
#	    	<td width="16">&nbsp;</td>
#	        <td class="modified">Modified lines:&nbsp;</td>
#	        <td class="modified">%s</td>
#	    </tr>
#	    <tr>
#	    	<td width="16">&nbsp;</td>
#	        <td class="added">Added line:&nbsp;</td>
#	        <td class="added">%s</td>
#	    </tr>
#	    <tr>
#	    	<td width="16">&nbsp;</td>
#	        <td class="removed">Removed line:&nbsp;</td>
#	        <td class="removed">%s</td>
#	    </tr>
#	</table>
#	<br/>""" % (changed_lnks, added_lnks, deleted_lnks)
		
		print >>outputFile, """
	<table width="100%%">
    
    <tr>
        <td width="16">&nbsp;</td>
        <td>
        %d lines<br/>
        %d bytes<br/>
        Last modified : %s<br/>
        <hr/>
        </td>
        <td width="16">&nbsp;</td>
        <td width="16">&nbsp;</td>
        <td>
        %d lines<br/>
        %d bytes<br/>
        Last modified : %s<br/>
        <hr/>
        </td>
    </tr>""" % (
       len(f1_lines), f1_stat[stat.ST_SIZE], 
       time.asctime(time.gmtime(f1_stat[stat.ST_MTIME])),
       len(f2_lines), f2_stat[stat.ST_SIZE], 
       time.asctime(time.gmtime(f2_stat[stat.ST_MTIME])))
		
		# Running through the differences...
		nl1 = nl2 = 0
		while not ( (nl1 >= len(f1_lines)) and (nl2 >= len(f2_lines)) ) :
			if ( added.has_key(nl2+1) ) :
				f2_lines[nl2]
	      # This is an added line
				print >>outputFile, """
	    <tr>
	        <td class="linenum">&nbsp;</td>
	        <td class="added">&nbsp;</td>
	        <td width="16">&nbsp;</td>
	        <td class="linenum"><a name="F2_%d">%d</a></td>
	        <td class="added">%s</td>
	    </tr>
	""" % (nl2+1, nl2+1, str2html(f2_lines[nl2]))
				nl2 += 1
			elif ( deleted.has_key(nl1+1) ) :
	      # This is a deleted line
				print >>outputFile, """
	    <tr>
	        <td class="linenum"><a name="F1_%d">%d</a></td>
	        <td class="removed">%s</td>
	        <td width="16">&nbsp;</td>
	        <td class="linenum">&nbsp;</td>
	        <td class="removed">&nbsp;</td>
	    </tr>
	""" % (nl1+1, nl1+1, str2html(f1_lines[nl1]))
				nl1 += 1
			elif ( changed.has_key(nl1+1) ) :
	      # This is a changed (modified) line
				print >>outputFile, """
	    <tr>
	        <td class="linenum"><a name="F1_%d">%d</a></td>
	        <td class="modified">%s</td>
	        <td width="16">&nbsp;</td>
	        <td class="linenum">%d</td>
	        <td class="modified">%s</td>
	    </tr>
	""" % (nl1+1, nl1+1, str2html(f1_lines[nl1]), 
	       nl2+1, str2html(f2_lines[nl2]))
				nl1 += 1
				nl2 += 1
			else :
	      # These lines have nothing special
				if ( not only_changes ) :
					print >>outputFile, """
	    <tr>
	        <td class="linenum">%d</td>
	        <td class="normal">%s</td>
	        <td width="16">&nbsp;</td>
	        <td class="linenum">%d</td>
	        <td class="normal">%s</td>
	    </tr>
	""" % (nl1+1, str2html(f1_lines[nl1]),
	       nl2+1, str2html(f2_lines[nl2]))
				nl1 += 1
				nl2 += 1
	            
	    # And finally, the end of the HTML
		print >>outputFile, """
	</table>
	"""
print >>outputFile, """
%s
	</span>
%s
	</span>
<hr/>
</body>
</html>""" % (notExistingFiles, identicalFiles)

outputFile.close()

# Remove output html file when there were no errors
if (retCode == 0):
	os.remove(outputFilename)

#if listFilename[0:4] == 'temp':
#	os.remove(listFilename)

sys.exit(retCode)
