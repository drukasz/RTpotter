# RTpotter
Open source software for conversion radiotherapy contour sequence data within DICOM files into STL files. It features a stepwise data conversion algorithm, that allows to preserve all the contour points and relations between them (i.e. connections described within curves on subsequent layers) and to remove possible artifacts (such as additional, small contours appearing on different layers as a result of automatic contouring procedures performed by RT planning systems).

INSTRUCTIONS FOR USE:
- download and install Python 3.x
- install pydicom library: in command line type: "pip install pydicom"
- edit RTpotter source file:
  * change the name of input DCM file - line 18
  * change the path and filename of output STL file
  * insert contour sequence name from the DCM file - case sensitive! must exactly match the name from RT planning system

IMPORTANT REMARKS AND LIMITATIONS:
- RTpotter is free and open-source software created for development and research purposes. It should not be used for actual medical applications! No warranty is given on the correctness of its operation, you do it at your own risk! The author takes no responsibility for inproper use of this code.
- the software has been only tested on datasets obtained using Eclipse, Varian Treatment Planning System. It worked well for all the investigated structures. 
- The current version is only capable of handling structures with one curve per single layer. In case if more complex geometries would be processed, they should be divided in parts fulfilling the above condition. 



REFERENCES:
Nowak, L. J., & Pawlowska, E. (2019). technical note: an algorithm and software for conversion of radiotherapy contour‐sequence data to ready‐to‐print 3D structures. Medical physics, 46(4), 1829-1832.
