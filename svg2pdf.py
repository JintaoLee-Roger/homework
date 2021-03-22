from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF, renderPM

drawing = svg2rlg("inversion/hw5/9.svg")
renderPDF.drawToFile(drawing, "inversion/hw5/9.pdf")
