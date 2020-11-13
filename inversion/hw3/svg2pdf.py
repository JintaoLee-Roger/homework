from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF, renderPM

drawing = svg2rlg("hw3e52.svg")
renderPDF.drawToFile(drawing, "hw3e52.pdf")