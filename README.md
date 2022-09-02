# Document Background Removing Tool (dbrt)
## Description
C++ project that automatically recognises and whitens parts of a digital image taken of a text document which do not contain any useful information.  
In addition it also improves readability by minimising the difference in lighting conditions across the image (brightening shadows).  
This is meant to be more of a resource than a ready to use program. The example main.cpp is by all means just an example; although it does work, it is very simple and it lacks error handling.  
Main use case is when you take a picture of a text only document with your phone (possibly by using Microsoft Lens or something similar) and the otherwise uniformly white paper is far from pure white and/or there are different colours of white appearing in the image (shadows, reflections and similar effects can cause this easily). This algorithm aims to make compression more efficient (by creating a large continous uniformly white area) and make images more readable.

## Prerequisites (dependencies)
In order to compile and use the program you need the following libraries:
### DevIL (Developer’s Image Library)
http://openil.sourceforge.net/download.php  
for reading and writing images.  
If you're on a popular Linux distribution you can probably find it in a repo using your package manager.  
Windows users will have to download the library files.

### CDT by artem-ogre
https://github.com/artem-ogre/CDT  
for triangulating points on the image. (used during lighting scaling)  
You need to place the CDT directory somewhere where the compiler can find it (if you're not planning on using it for anything else you can just put the folder next to the source files).  
If you wish to use a different triangulation library you can modify the Document::Triangulation function accordingly.

### Boost (Program options)
https://www.boost.org/  
https://www.boost.org/doc/libs/1_80_0/doc/html/program_options.html  
Program options is one of the few boost libraries that do not work just by including the .hpp files, you need to link a library binary.  (The library file is not distributed with the include files, if you need help with building it here is the guide: https://www.boost.org/doc/libs/1_55_0/more/getting_started/unix-variants.html#prepare-to-use-a-boost-library-binary for Linux and https://www.boost.org/doc/libs/1_55_0/more/getting_started/windows.html#prepare-to-use-a-boost-library-binary for Windows)  
This is not necessary for the background recognition and whitening, only for the example program (main.cpp) so if you only plan to use this project as a resource and not as a stand-alone program you don't have to worry about boost.  

## Usage
It is __highly__ recommended to use some kind of automated cropping tool (like the Microsoft Lens app) while taking the photo, or manually cropping the image before processing. Areas of the image that are not part of the document can skew the results (especially the brightness scaling).  

Syntax:  
`dbrt input_image.jpg output_image.png`  
by default the program outputs only a greyscale image, howewer if you want a colour image you can get one by using the -C option and giving a path for the output file like so:  
`dbrt input_image.bmp -C rgb_output.tiff gsc_output.tga`  
You can find out more about the possible options by running `dbrt --help`.  
The first argument should alway be the path to the input image file, and the last argument is path to the greyscale output file.  
Supported file formats are defined by DevIL and you can read a full list of them in the Developer's Image Library Manual under the section "Appendix B: Supported File Formats" (http://openil.sourceforge.net/docs/DevIL%20Manual.pdf)
