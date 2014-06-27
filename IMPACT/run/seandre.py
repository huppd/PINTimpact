#!/usr/bin/python
import os, sys
usage = "usage: %s search_text replace_text [infile [outfile]]" %         os.path.basename(sys.argv[0])


def sar( stext, rtext, infile=0, outfile=0):
  #stext = sys.argv[1]
  #rtext = sys.argv[2]
  input = sys.stdin
  output = sys.stdout
  
  if infile != 0:
    input = open( infile )
  
  if outfile != 0 and outfile != infile:
    output = open( outfile, 'w')
  else:
    outfile = infile
    if(infile!=0):
      output = open( outfile+'tmp', 'w')
  
  i = 0
  l = 0
  print 
  print '___________________________________________________________'
  for s in input:
    l += 1
    if( stext in s):
		if( stext in rtext ):
		  print 'replaced', s, 'with', rtext
		  s = rtext
		  output.write( s+"\n" )
		  i += 1
		  #print(s,output)
		else:
		  s = stext+rtext
		  print 'replaced', s, 'with', s
		  output.write( s+"\n" )
		  i += 1
		  #print(s,output)
    else:
      output.write( s )
  
  if infile != 0:
    input.close()
    if outfile != 0:
      output.close()
      if( infile==outfile ):
        os.rename( outfile+'tmp', outfile )
  print 
  print '___________________________________________________________'
  print 'replaced', i, 'of', l, 'lines'
	

if __name__ == "__main__":
	print "There are %s args " %len(sys.argv)
	if len(sys.argv) <= 3:
		print usage
	if len(sys.argv) == 3:
		sar( sys.argv[1], sys.argv[2] )
	if len(sys.argv) == 4:
		sar( sys.argv[1], sys.argv[2], sys.argv[3] )
	if len(sys.argv) > 4:
		sar( sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] )


