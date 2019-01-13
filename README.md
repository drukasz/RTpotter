# RTpotter
Open source software for conversion radiotherapy contour sequence data within DICOM files into STL files. It features a stepwise data conversion algorithm, that allows to preserve all the contour points and relations between them (i.e. connections described within curves on subsequent layers) and to remove possible artifacts (such as additional, small contours appearing on different layers as a result of automatic contouring procedures performed by RT planning systems).

INSTRUCTIONS FOR USE:
- download and install Python 3.x
- install pydicom library: in command line type: "pip install pydicom"
- edit RTpotter source file:
  * change the name of input DCM file - line 18
  * change the path and filename of output STL file
  * insert contour sequence name from the DCM file - case sensitive! must exactly match the name from RT planning system

