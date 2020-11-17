from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF, renderPM

drawing = svg2rlg("hw3e1e.svg")
renderPDF.drawToFile(drawing, "hw3e1e.pdf")
