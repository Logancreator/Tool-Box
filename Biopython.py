##Converting other formats to GFF3

from BCBio import GFF
from Bio import SeqIO

in_file = "./GCA_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gbff"
out_file = "./GCA_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gbff.gff"
in_handle = open(in_file)
out_handle = open(out_file, "w")

GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)

in_handle.close()
out_handle.close()
