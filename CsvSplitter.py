import os
import sys
import csv
import re


def split(filehandler, delimiter=',', row_limit=10000,
          output_name_template='output_%s.csv', output_path='.', keep_headers=True):
    """
    Splits a CSV file into multiple pieces.

    A quick bastardization of the Python CSV library.

    Arguments:

        `row_limit`: The number of rows you want in each output file. 10,000 by default.
        `output_name_template`: A %s-style template for the numbered output files.
        `output_path`: Where to stick the output files.
        `keep_headers`: Whether or not to print the headers in each output file.

    Example usage:

        >> from toolbox import csv_splitter;
        >> csv_splitter.split(open('/home/ben/input.csv', 'r'));

    """
    reader = csv.reader(filehandler, delimiter=delimiter)
    current_piece = 1
    current_out_path = os.path.join(
        output_path,
        output_name_template % current_piece
    )
    current_out_writer = csv.writer(open(current_out_path, 'w'), delimiter=delimiter)
    current_limit = row_limit
    if keep_headers:
        headers = reader.next()
        current_out_writer.writerow(headers)
    for i, row in enumerate(reader):
        if i + 1 > current_limit:
            current_piece += 1
            current_limit = row_limit * current_piece
            current_out_path = os.path.join(
                output_path,
                output_name_template % current_piece
            )
            current_out_writer = csv.writer(open(current_out_path, 'w'), delimiter=delimiter)
            if keep_headers:
                current_out_writer.writerow(headers)
        current_out_writer.writerow(row)


def breakDownLargeFiles(rootDataDir, rowLimit,size):

    dataFile = ".{1,20}.tsv$"

    for root, dirs, files in os.walk(rootDataDir):
        if root.endswith("mut"):
            for afile in files :
                filePath = root + "/" + afile
                if re.match(dataFile,afile) and os.stat(filePath).st_size >  int(size) :
                    print("Breaking up CSV ")
                    print(filePath)
                    fileHandler = open(filePath, 'r')
                    name = "dataChunk_%s.tsv"
                    split(fileHandler, delimiter='\t', row_limit=int(rowLimit),output_name_template=name, output_path=root, keep_headers=True )


if len(sys.argv) == 4:
    dataDir = sys.argv[1]
    rowLimit = sys.argv[2]
    size = sys.argv[3]
    breakDownLargeFiles(dataDir, rowLimit, size)
else:
    print("please provider the data directory and row limit")