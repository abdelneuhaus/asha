import os
from PyPDF2 import PdfMerger

def merge_subpdf(sublist_of_pdf, index, exp):
    for i in range(len(sublist_of_pdf)):
        merger = PdfMerger()
        for pdf in sublist_of_pdf[i]:
            merger.append(pdf)
        merger.write('results/'+exp+'/'+index[i]+".pdf")
        merger.close()
        

def delete_subpdf(exp):
    for parent, dirnames, filenames in os.walk('results/'+exp+'/'):
        for fn in filenames:
            if fn.lower().endswith('561.pdf') or fn.lower().endswith('405.pdf'):
                os.remove(os.path.join(parent, fn))