import csv, argparse, os

parser = argparse.ArgumentParser()
parser.add_argument("input_file")
parser.add_argument("window_size")
parser.add_argument("output_file")

args = parser.parse_args()

targets = []
window_size = args.window_size

class Bed_file:

    def __init__(self,name,chrom,start,end):
	self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name

    def get_new_start(self,start,end,window_size):
        start = int(start)
        end = int(end)
        window_size = int(window_size)
        # If interval is equal or smaller than two windows move to middle of interval and take off one window length for start coordinate
        if (end - start) <= (window_size * 2):
            self.newstart = (start + (end - start)/2) - window_size
        # Otherwise, put a window over midpoint of interval, and remove window lengths until windows overlaps interval
        else:
            mid = (start + (end - start)/2)- (window_size / 2)
            while mid > (start - (window_size / 2)):
                self.newstart = mid
                mid = mid - window_size

    def get_new_end(self,start,end,window_size):
        start = int(start)
        end = int(end)
        window_size = int(window_size)
        # If interval is equal or smaller than two windows move to middle of interval and add on one window length for end coordinate
        if (end - start) <= (window_size * 2):
            self.newend = (start + (end - start)/2) + window_size
        # Otherwise, put a window over midpoint of interval, and add window lengths until windows overlaps interval
        else:
            mid = (start + (end - start)/2)- (window_size / 2)
            while mid < (end + (window_size / 2)):
                self.newend = mid
                mid = mid + window_size
        
with open(args.input_file) as infile:
    # Open bed file and read line into new Class
    csvreader = csv.reader(infile, delimiter = '\t')
    for row in csvreader:
        row[3] = Bed_file(row[3], row[0], row[1], row[2])
        # Alter intervals to make them divisible by desired window length
        row[3].get_new_start(row[1],row[2], window_size)
        row[3].get_new_end(row[1],row[2],window_size)
        # Make array of intervals
        targets.append(row[3])

with open('bed_out.bed', 'wb') as outfile:
    # For all intervals in targets array write to new bed file
    csvwriter = csv.writer(outfile, delimiter='\t')
    for target in targets:
        out_row = [target.chrom, target.newstart, target.newend, target.name]
        csvwriter.writerow(out_row)

# Use bedtools to divide new intervals into desired size and produce final file
#command = "bedtools makewindows -b bed_out.bed -w %s -i src > panel_CNV_%sbp.bed" %(window_size, window_size)
command = "bedtools makewindows -b bed_out.bed -w %s -i src >  %s" %(args.window_size, args.output_file)

os.system(command)



    



