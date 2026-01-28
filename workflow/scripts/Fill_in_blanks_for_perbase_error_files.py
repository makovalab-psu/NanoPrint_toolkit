import sys

#### Parse argument ####
Perbase_error_file = sys.argv[1]
chr_size = int(sys.argv[2])
out_file = sys.argv[3]

with open(Perbase_error_file, 'r') as perbase_error, open(out_file, 'w') as output_file:
  NEXT_LINE = True
  WRITE_LINE = False
  for i in range(1, chr_size + 1):
      if NEXT_LINE:
        LINE = perbase_error.readline()
        if not LINE:
            NEXT_LINE = False
            OUT_LINE = str(i)+"\n"
        else:
            LINE_split = LINE.strip().split("\t")
            N = int(LINE_split[1])
      if N == i:
        NEXT_LINE = True
        OUT_LINE = LINE
      else:
        NEXT_LINE = False
        OUT_LINE = str(i)+"\n"
      output_file.write(OUT_LINE)

