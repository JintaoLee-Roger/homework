from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF, renderPM

drawing = svg2rlg("inversion/hw4/hw4e4d3.svg")
renderPDF.drawToFile(drawing, "inversion/hw4/hw4e4d3.pdf")
