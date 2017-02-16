""" 
Reformat Capture Hi-C files by Mifsud 2015 et al.

Ths script takes a promoter-promoter file from the Misfusd 2015 et al. study and
converts it into a tab-separated fromat in which each gene pair is one line.

"""
epilog="""Jonas Ibn-Salem <j.ibn-salem@uni-mainz.de> 2017-02-15"""



import argparse

def commandline():
    """ returns the args from commandline. """
    command_parser = argparse.ArgumentParser(description=__doc__, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    command_parser.add_argument('-i','--input_file', type=str, required=True, help='input file.')
    command_parser.add_argument('-o','--output_file', type=str, required=True, help='output file.')
    args =  command_parser.parse_args()
    return args




def reformat_file(inFile, outFile):
  """
  Ths function takes a promoter-promoter file from the Misfusd 2015 et al. 
  study and converts it into a tab-separated file fromat in which each 
  gene pair is one line.
  """
  
  with open(outFile, "w") as outHandle:
    
		# write header line
		outLine = ["g1", "g2", "raw_count", "log(obs/exp)"]
		
		outHandle.write("\t".join(outLine) + "\n")


		for i, line in enumerate(open(inFile)):
			
			if not i == 0:
				
				sp = line.strip().split("\t")
				
				# get row interaction counts and normalized obs/exp values
				rawCount = sp[12]
				obsExp = sp[13]
				
				genes1 = sp[4].split("|")
				genes2 = sp[10].split("|")
				
				#~ print(g1, g2, rawCount)
				
				# iterate over all pairs
				for g1 in genes1:

					for g2 in genes2:
						
						outLine = [g1, g2, rawCount, obsExp]
						
						outHandle.write("\t".join(outLine) + "\n")
						
				


if __name__ == "__main__":
    args = commandline()
    
    reformat_file(args.input_file, args.output_file)
    
