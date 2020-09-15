------------------
Usage Instructions
------------------
Before Running:
	All code is in a single Python file, "ProteinAlignment.py"
	Place ProteinAlignment.py file in the same directory as the protein sequence .fasta files you wish to use
From a terminal:
	Navigate to the directory with the ProteinAlignment.py and .fasta files
	type "python ProteinAlignment.py"
	On the console, the user will be prompted to type the name of a fasta file: "FASTA file for first protein: "
	Type the name of a fasta file within the same directory ie "1_protein.fasta" then push enter
	(alternatively you can type the entire address of a fasta file in another directory)
	The user will be prompted for a second file
	Enter the name or address of the file for the protein you wish to align with the first
	If at either prompt an invalid file name is given, or the file cannot be found, "file not found" will be printed and the program will restart
	If valid files are given, the program will now run the Needleman-Wusch and the Smith-Waterman algorithms on the two given proteins
Output:
	The program will first output the Needlman Wusch gobal alignment results
	First will be the total alignment score, "Score: "
	Next will be the two protein sequences, 80 characters at a time, with gaps represented by "-"
	In between the two sequences there will be a "*" for a mismatch, "|" for a match, and whitespace otherwise
	
	Next the program will output the Smith-Waterman local alignment results 
	First will be the score for the best local alignment
	Next will be the an alignment for part of the protein sequences in the same format as stated above
Exiting:
	After the output, the program will asks for files again
	You can continue matching sequences until you exit
	To exit, when asked for a file, simple type "exit"
	The program will terminate