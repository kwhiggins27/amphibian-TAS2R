Script by Matthew Hill, shared with permission
PREREQUISITES:
I recommend running with the tak default python, not the conda python, so "conda deactivate" if you are in a conda instance
If you use a conda python, you may also need to install matplotlib
DNA Features Viewer #install with "pip install dna_features_viewer"

Call "perl GenePlot.pl" with the following arguments:

-i queryfile.file
-g genome.gff3 #Genome file must be in gff3 format
-o output directory for results

-l {optional} #Include only if you want gene labels on your plot; does not need an argument
-p {optional} #Include only if you want PDF output format (vector graphics); does not need an argument

#Example run without labels
perl GenePlot.pl -i queries.txt -g tagetes_patula.gff3 -o /path/to/output/folder

#Example run with labels and pdf output
perl GenePlot.pl -i queries.txt -g tagetes_patula.gff3 -o /path/to/output/folder -l -p

#If GenePlot.pl is in a different folder
perl /path/to/GenePlot.pl -i queries.txt -g tagetes_patula.gff3 -o /path/to/output/folder -l -p

#Query file format: Two gene IDs separated by a dash per line.
#Genes must be located on the same contig or chromosome or unexpected behavior may occur.
#Unexpected behavior may occur with very large (>>100Kbp) ranges.

#Example query file
FUN_000001-FUN_000010
FUN_123456-FUN_123457

A pair of files will be generated for each query-- .py and graphic.
The .py is called by GenePlot.pl to generate the graphic.
You can change the gene colors by editing the .py file and re-running it, which will regenerate the graphic.

To get the length of DNA in the visualization, check the sequence_length field on line 14 of the generated .py file.
