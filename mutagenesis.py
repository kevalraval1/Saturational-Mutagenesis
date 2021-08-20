
# TODO: Clean up, make more concise

from Bio import SeqIO
from tkinter import filedialog
from tkinter import *

window = Tk()
window.title("Saturational Mutagenesis Tool")

def selectFile():
    global filename
    filename = filedialog.askopenfilename(initialdir='/',title='Select File',
    filetypes=(('FASTA File','*.fa'),('All Files', '*.*')))
    if(filename!=''):
        label = Label(frame, text=filename, bg = 'gray')
        spaceText.window_create("end", window=label)
        spaceText.insert('end', '\n')

def main():
    if (intronsFlank.get()) == 1:
        intronSize = 20
    else:
        intronSize = 0
    with open(filename, "r") as file:
        recordList = SeqIO.parse(file, "fasta")
        exonList = []
        intronList = []
        exonStarts = [] 
        for record in recordList:
            if "exon" in str(record.description):
                exonList.append(str(record.seq))
            elif "intron" in str(record.description):
                intronList.append(str(record.seq))
            elif "chromosome" in str(record.description):
                tempList = str(record.description).split(":")
                geneStart = int(tempList[4])
                geneEnd = int(tempList[5])
                chromosome = int(tempList[3])
                strand = int(tempList[6])
        if strand == -1:
            for count, exon in enumerate(exonList):
                if count == 0:
                    exonStarts.append(geneEnd - len(exon) - intronSize)
                    continue
                elif count == (len(exonList) - 1):
                    addLength = exonStarts[count - 1] - len(intronList[count - 1]) - len(exon) + intronSize
                else:
                    addLength = exonStarts[count - 1] - len(intronList[count - 1]) - len(exon)
                exonStarts.append(addLength)
            if intronSize == 20: 
                for count, exons in enumerate(exonList):
                    if count == 0:
                        tempTuple = (exonList[0], intronList[count][0:20])
                        exonList[0] = ''.join(tempTuple)
                    elif count == (len(exonList) - 1):
                        tempTuple = (intronList[count - 1][len(intronList[count-1]) - 20:], exonList[count])
                        exonList[count] = ''.join(tempTuple)
                    else:
                        tempTuple = (intronList[count - 1][len(intronList[count-1]) - 20:], exonList[count], intronList[count][0:20])
                        exonList[count] = ''.join(tempTuple)
            for count, exon in enumerate(exonList):
                TSVFile = open(f"{file.name}_exon_{count + 1}.tsv","w+")
                TSVFile.write("CHROM\tSTART\tEND\tREF\tALT\n")
                for baseNum, base in enumerate(reversed(exon)):
                    if base == "A":
                        tempList = []
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum), str(exonStarts[count] + baseNum + 1), "T", "A")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum), str(exonStarts[count] + baseNum + 1), "T", "C")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum), str(exonStarts[count] + baseNum + 1), "T", "G")
                        tempList.append("\t".join(tempTuple))
                        TSVFile.write("\n".join(tempList))
                    elif base == "T":
                        tempList = []
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum), str(exonStarts[count] + baseNum + 1), "A", "T")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum), str(exonStarts[count] + baseNum + 1), "A", "C")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum), str(exonStarts[count] + baseNum + 1), "A", "G")
                        tempList.append("\t".join(tempTuple))
                        TSVFile.write("\n".join(tempList))
                    elif base == "G":
                        tempList = []
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum), str(exonStarts[count] + baseNum + 1), "C", "A")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum), str(exonStarts[count] + baseNum + 1), "C", "T")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum), str(exonStarts[count] + baseNum + 1), "C", "G")
                        tempList.append("\t".join(tempTuple))
                        TSVFile.write("\n".join(tempList))
                    elif base == "C":
                        tempList = []
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum), str(exonStarts[count] + baseNum + 1), "G", "A")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum), str(exonStarts[count] + baseNum + 1), "G", "T")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum), str(exonStarts[count] + baseNum + 1), "G", "C")
                        tempList.append("\t".join(tempTuple))
                        TSVFile.write("\n".join(tempList))
                    TSVFile.write("\n")
                TSVFile.close()
        else:
            for count, exon in enumerate(exonList):
                if count == 0:
                    exonStarts.append(geneStart)
                    continue
                elif count == 1: 
                    addLength = exonStarts[count - 1] + len(exonList[count - 1]) + len(intronList[count - 1]) - intronSize
                else:
                    addLength = exonStarts[count - 1] + len(exonList[count - 1]) + len(intronList[count - 1])
                exonStarts.append(addLength)
            if intronSize == 20:
                for count, exons in enumerate(exonList):
                    if count == 0:
                        tempTuple = (exonList[0], intronList[count][0:20])
                        exonList[0] = ''.join(tempTuple)
                    elif count == (len(exonList) - 1):
                        tempTuple = (intronList[count - 1][len(intronList[count-1]) - 20:], exonList[count])
                        exonList[count] = ''.join(tempTuple)
                    else:
                        tempTuple = (intronList[count - 1][len(intronList[count-1]) - 20:], exonList[count], intronList[count][0:20])
                        exonList[count] = ''.join(tempTuple)
            for count, exon in enumerate(exonList):
                TSVFile = open(f"{file.name}_exon_{count + 1}.tsv","w+")
                TSVFile.write("CHROM\tSTART\tEND\tREF\tALT\n")
                for baseNum, base in enumerate(exon):
                    if base == "T":
                        tempList = []
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum - 1), str(exonStarts[count] + baseNum), "T", "A")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum - 1), str(exonStarts[count] + baseNum), "T", "C")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum - 1), str(exonStarts[count] + baseNum), "T", "G")
                        tempList.append("\t".join(tempTuple))
                        TSVFile.write("\n".join(tempList))
                    elif base == "A":
                        tempList = []
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum - 1), str(exonStarts[count] + baseNum), "A", "T")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum - 1), str(exonStarts[count] + baseNum), "A", "C")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum - 1), str(exonStarts[count] + baseNum), "A", "G")
                        tempList.append("\t".join(tempTuple))
                        TSVFile.write("\n".join(tempList))
                    elif base == "C":
                        tempList = []
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum - 1), str(exonStarts[count] + baseNum), "C", "A")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum - 1), str(exonStarts[count] + baseNum), "C", "T")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum - 1), str(exonStarts[count] + baseNum), "C", "G")
                        tempList.append("\t".join(tempTuple))
                        TSVFile.write("\n".join(tempList))
                    elif base == "G":
                        tempList = []
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum - 1), str(exonStarts[count] + baseNum), "G", "A")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum - 1), str(exonStarts[count] + baseNum), "G", "T")
                        tempList.append("\t".join(tempTuple))
                        tempTuple = (str(chromosome), str(exonStarts[count] + baseNum - 1), str(exonStarts[count] + baseNum), "G", "C")
                        tempList.append("\t".join(tempTuple))
                        TSVFile.write("\n".join(tempList))
                    TSVFile.write("\n")
                TSVFile.close()
    label = Label(frame, text = "Done", bg = 'green')
    spaceText.window_create("end", window=label)
    spaceText.insert('end', '\n')

canvas = Canvas(window, height = 100, width = 600)
canvas.pack()

frame = Frame(window,relief = 'groove')
frame.place(relx = 0.1, rely = 0.1, relwidth = 0.8, relheight = 0.8)

welcome = Label(frame, text = "Welcome to the Saturational Mutagenesis Tool", fg = "Black")
welcome.pack(side = "top")

spaceText = Text(frame,width=40,height=20, borderwidth=0)
spaceText.pack(side='right',fill='both' ,expand=True)

intronsFlank = IntVar()
c = Checkbutton(window, text='20 BP Intron Flanks', variable = intronsFlank)
c.pack()

openFile = Button(window, text='Open FASTA File', padx = 10, pady = 5, fg = 'black', bg = 'gray', command = selectFile)
openFile.pack()

enterButton = Button(window, text = "Start", padx = 10, pady = 5, fg = "Black", bg = "gray", command = main)
enterButton.pack()

window.mainloop()



