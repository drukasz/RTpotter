# RTpotter
Open source software for conversion radiotherapy contour sequence data within DICOM files into STL files. It features a stepwise data conversion algorithm, that allows to preserve all the contour points and relations between them (i.e. connections described within curves on subsequent layers) and to remove possible artifacts (such as additional, small contours appearing on different layers as a result of automatic contouring procedures performed by RT planning systems).

INSTRUCTIONS FOR USE:
- download and install Python 3.x
- install pydicom library: in command line type: "pip install pydicom"
- edit RTpotter source file:

