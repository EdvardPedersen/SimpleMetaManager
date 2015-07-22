### Script to convert dat files (uniprot) to different formats
import dis
import re
from optparse import OptionParser
import sys
import os
import mmap
import subprocess
import shlex
import gdbm as dbm
import glob
import time
import logging
import sqlite3

timer = time.time()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler("logging.log", mode='w')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)




def time_start():
	global timer
	timer = time.time()

def time_stop(text):
	global timer
	time_elapsed = time.time() - timer
	logger.info("Stage " + text + ": " + str(time_elapsed))

def node_list():
	nodes = list()
	for i in range(9):
		nodes.append("compute-0-" + str(i))
	return nodes

def convert_to_fasta(input):
	output = list()
	inputList = input.splitlines(True)
	for line in inputList:
		if(line[:2] == "ID"):
			output.append(">" + line)
		elif(line[:2] == "  "):
			output.append(line.replace(' ', ''))
	return "".join(output)

class UniprotFile:
	def __init__(self, filename="uniprot_sprot.dat"):
		openfile = open(filename, 'r')
		self.file = mmap.mmap(openfile.fileno(), 0, access=mmap.ACCESS_READ)

	def get_entry(self):
		if(self.file.tell() == self.file.size()):
			raise EOFError()
		return self.get_entry_fast()

	def _get_entry_it(self):
		readFunc = self.file.readline
		inLine = readFunc()
		while(not inLine[:2] == "//" or inLine == ''):
			yield inLine
			inLine = readFunc()
		yield inLine

	def get_entry_fast(self):
		readFunc = self.file.readline
		retString = ""
		inLine = readFunc()
		while(not inLine[:2] == "//" or inLine == ''):
			retString = retString + inLine
			inLine = readFunc()
		return retString

	def reset_fp(self):
		self.file.seek(0)
		
class Output:
	def __init__(self):
		pass

	def export_fasta(self, elements, filename, reader):
		outFile = open(filename, 'w')
		for i in xrange(elements):
			try:
				outFile.write(convert_to_fasta(reader.get_entry()))
			except EOFError:
				break
		outFile.close()
		reader.reset_fp()

	def export_dat(self, elements, filename, reader):
		outFile = open(filename, 'w')
		for i in xrange(elements):
			try:
				outFile.write(reader.get_entry())
			except EOFError:
				break
		outFile.close()
		reader.reset_fp()

	def export_dbm(self, elements, filename, reader):
		outFile = dbm.open(filename, 'nf')
		for i in xrange(elements):
			try:
				entry = reader.get_entry()
				id = entry.splitlines(True)[0]
				outFile[id] = entry
			except EOFError:
				break
		outFile.sync()
		outFile.close()
		reader.reset_fp()

	def export_sqlite(self, elements, filename, reader):
		conn = sqlite3.connect(filename)
		conn.text_factory = str
		c = conn.cursor()
		c.execute("CREATE TABLE uniprot (id text, value text)")
		execlist = []
		for i in xrange(elements):
			try:
				entry = reader.get_entry()
				id = entry.splitlines(True)[0]
				execlist.append((id, entry))
				if(len(execlist) > 2000):
					c.executemany("INSERT INTO uniprot VALUES (?,?)", execlist)
					execlist = []
					conn.commit()
			except EOFError:
				break
		conn.commit()
		conn.close()
		reader.reset_fp()
		

def simple_copy(filename):
	with open(filename, 'r') as f:
		print(f.read())

def simple_copy_mmap(filename, outfile):
	with open(filename, 'r') as f:
		with open(outfile, 'w') as of:
			mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
			while(mm.tell() < mm.size()):
				of.write(mm.read(10000))

def formatdb(infile, outfile):
	args = shlex.split("/opt/bio/ncbi/bin/formatdb -t " + outfile + " -p T -i " + infile)
	subprocess.call(args, shell=False, stderr=subprocess.STDOUT)

def copy_to_nodes(infile, dest_path):
	files = glob.glob(infile + "*")
	command_list = list()
	for element in node_list():
		for file in files:
			print file
			lastName = os.path.basename(file)
			command = shlex.split("rsync -z -W --dirs " + file + " " + element + ":" + dest_path)
			command_list.append(command)

	parallel_execute(command_list)
	
	
def parallel_execute(command_list):
	waitList = list()
	for element in command_list:
		try:
			waitList.append(subprocess.Popen(element, bufsize=-1))
			while(len(waitList) > 20):
				for element in waitList:
					if(element.poll() != None):
						print("RSync done, removing from list")
						waitList.remove(element)
				time.sleep(0.2)
		except Exception as E:
			print("Error calling subprocess" + str(E))
			raise
	for element in waitList:
		retcode = element.wait()
		if not retcode == 0:
			print command_list
	

if __name__ == "__main__":
	options = OptionParser()
	options.add_option("-f", "--fasta", dest="fastafilename", action="store", type="string", help="Fasta output file")
	options.add_option("-d", "--dat", dest="datfilename", action="store", type="string", help="DAT output file")
	options.add_option("-i", "--input", dest="input", action="store", type="string", help="Input DAT file")
	options.add_option("-n", "--num", dest="num", action="store", type="int", help="Number of entries to generate, default = all", default=sys.maxsize)
	options.add_option("-D", "--database", dest="database", action="store", type="string", help="Export a gdbm database for annotation")
	(opts, args) = options.parse_args()
	if(not opts.input):
		options.error("No input file specified")
	upFile = UniprotFile(filename=opts.input)
	outputFasta = Output()
	if(opts.fastafilename):
		logger.info("Creating fasta file")
		time_start()
		outputFasta.export_fasta(opts.num, opts.fastafilename, upFile)
		time_stop("Write fasta")

		logger.info("Running formatdb")
		time_start()
		formatdb(opts.fastafilename, "testdb")
		time_stop("formatdb")

		logger.info("Copying files to nodes")
		time_start()
		copy_to_nodes(opts.fastafilename + ".", "/state/partition1/test/")
		time_stop("copy BLAST db to nodes")
	if(opts.datfilename):
		logger.info("Exporting DAT file")
		time_start()
		outputFasta.export_dat(opts.num, opts.datfilename, upFile)
		time_stop("Create dat")
	if(opts.database):
		#logger.info("Exporting database in GDBM format")
		#time_start()
		#outputFasta.export_dbm(opts.num, opts.database, upFile)
		#time_stop("database creation")

		logger.info("Exporting DB in sqlite3 format")
		time_start()
		outputFasta.export_sqlite(opts.num, opts.database + ".sqlite", upFile)
		time_stop("sqlite3 DB")

		logger.info("Copying files to nodes")
		time_start()
		copy_to_nodes(opts.database, "/state/partition1/test/")
		time_stop("copy database to nodes")
