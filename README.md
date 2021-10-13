READ ME - snp_opdracht

author: Jurrien de Jong
date: 13-10-2021
files: SNP_Annotation.py, seq.txt, data.msf

Quick Note:
	If you are looking for some test data, you can find it at the usage section!

Description of the code:
	The program SNP_Annotation.py will handle a few input arguments:
    -i : Input file in .msf ( so the MSA alignment )
    -seq : The sequence ( DNA or Protein ) to be aligned with the MSA
    -pos : The position in the above sequence where a SNP needs to be implemented
    -SNP : The SNP as a string; like "A"
  When provided with arguments, the program will read the .msf file and safe the alignment in a variable. The translate function accepts both DNA and protein sequences. The only   difference is that DNA need to pe translated first, which it does. After obtaining a proper protein sequence, the SNP will be implemented at the given position ( if valid ) by   another function. The alignment itself will be done by the function 'get_align_score'. It will compare each column of the MSA with the given sequence, and give a score between   0 and the length of the sequence. So 10 matches out of 50 will give a score of 10 / 50. Finally, the scores will be processed by the function 'write_results'. It will           calculate a percentage score between 0 and 100, and print a corresponding message. Afterward it will end the program.

Usage:
	Please open your terminal and head to the folder containing the files, if you cannot get to it
	please visit the links in the support section.

	SNP_Annotation.py can be run with the following syntax:

	python3 SNP_Annotation.py -i [input file in .msf format ] -seq [.txt file with a sequence] -pos [position of SNP] -snp [SNP as a string]
	
  Here are a few examples:

  python3 SNP_Annotation.py -i data.msf -seq -pos 0 "A"
  
  python3 SNP_Annotation.py -i data.msf -seq -pos 0 "M"
  
  python3 SNP_Annotation.py -i data.msf -seq -pos  "A"
  

Support:
	Please view these pages if you are having trouble with cd, running the command or classes
	https://docs.python.org/3/tutorial/classes.html
	
	https://linuxize.com/post/linux-cd-command/
	
	https://realpython.com/run-python-scripts/
  
  For further questions, please send me an email:
  ju.de.jong@st.hanze.nl
