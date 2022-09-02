#pragma once

#include <vector>       /* vector */
#include <string>       /* string */
#include <memory>       /* unique_ptr */

/* Libraries: */
#include<IL/il.h>   /* DevIL library */

/* Options: */
//#define RGBAVERAGE

/* Macros: */
#ifndef MIN_ZONE_SIDE_LENGTH
#define MIN_ZONE_SIDE_LENGTH 33
#endif // !MIN_ZONE_SIDE_LENGTH
#ifndef MAX_ZONE_SIDE_LENGTH
#define MAX_ZONE_SIDE_LENGTH 50
#endif // !MAX_ZONE_SIDE_LENGTH
#if (!defined(MIN_ZONES_PER_SIDE) || MIN_ZONES_PER_SIDE < 3)
#define MIN_ZONES_PER_SIDE 3
#endif // (!defined(MIN_ZONES_PER_SIDE) || MIN_ZONES_PER_SIDE < 3)
#define MIN_IMAGE_WIDTH MIN_ZONE_SIDE_LENGTH * MIN_ZONES_PER_SIDE
#define MIN_IMAGE_HEIGHT MIN_IMAGE_WIDTH

/* Typedefs: */
typedef unsigned char ubyte;
typedef std::vector<unsigned long int> frequencies;

/* Exceptions: */
class image_exception : public std::exception
{
protected:
	std::string message;
public:
	image_exception(const char* msg) : std::exception()
	{
		message = msg;
	}
	image_exception(const std::string msg) : std::exception()
	{
		message = msg;
	}
	std::string what()
	{
		return message;
	}
};
class image_too_small : public image_exception
{
public:
	image_too_small(const char* msg) : image_exception(msg) { }
	image_too_small(const std::string msg) : image_exception(msg) { }
};

/* Structs: */
/// <summary>
/// Constants for absolute white and black greyscale pixels.
/// </summary>
struct GScPixel
{
	const static ubyte WHITE = 0xff;
	const static ubyte BLACK = 0x00;
};

struct RGBA
{
public:
	ubyte red;
	ubyte green;
	ubyte blue;
	ubyte alpha;
	RGBA() { red = green = blue = alpha = 0x00; }
	RGBA(ubyte red, ubyte green, ubyte blue, ubyte alpha = 0xff)
	{
		this->red = red;
		this->green = green;
		this->blue = blue;
		this->alpha = alpha;
	}
	ubyte ToGrey()
	{
		//Average:
#ifdef RGBAVERAGE
		return red / 3 + green / 3 + blue / 3;
#endif // RGBAVERAGE
		//Weighted:
#ifndef RGBAVERAGE
		return char((2126 * red + 7152 * green + 722 * blue) / 10000);
#endif // !RGBAVERAGE
	}
	static RGBA WHITE() { return RGBA(0xff, 0xff, 0xff); }
	static RGBA BLACK() { return RGBA(0x00, 0x00, 0x00); }
};

struct zonestat;

/// <summary>
/// Represents the position of a zone or a pixel on the image.
/// </summary>
struct Point {
	int x;
	int y;

	zonestat* zone;

	Point(const int x = 0, const int y = 0, zonestat* zone = nullptr)
	{
		this->x = x;
		this->y = y;
		this->zone = zone;
	}

	bool equals(const Point& that)
	{
		return x == that.x && y == that.y;
	}
	/// <summary>
	/// Distance squared. <br/> since it's only ever used to compare to distances no need to calculate root
	/// </summary>
	static int DistanceSq(const Point& A, const Point& B)
	{
		return (A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y);
	}
};

/// <summary>
/// Stores information for scaling about a zone.
/// </summary>
/// <param name="width">the width of each zone in pixels</param>
/// <param name="height">the heights of each zone in pixels</param>
struct zonestat {
	std::unique_ptr<Point> point = nullptr;

	static int cols;
	static int rows;

	static double width;
	static double height;

	static double max_entropy;
	static double most_common_relative_entropy;
	static ubyte most_common_entropy_byte;

	ubyte darkest, whitest, peak, last_local_max;
	unsigned int darkestval, whitestval, peakval, maxval;
	double entropy = -1;
	double relative_entropy = -1;
	ubyte entropy_byte = 0;

	double scale_factor = 1;

	zonestat(const frequencies frequency)
	{
		darkest = 0, whitest = 255, peak = 0, last_local_max = 0;
		int idx = 0;
		while (idx < 256 && frequency[idx] == 0) idx++;
		darkest = idx;
		for (; idx < 256; idx++)
			if (frequency[idx] >= frequency[peak]) peak = idx;
		idx = 255;
		while (idx > 0 && frequency[idx] == 0) idx--;
		whitest = idx;
		last_local_max = whitest;
		while (idx > 1 && frequency[idx - 1] >= frequency[idx]) idx--;
		last_local_max = idx;

		darkestval = frequency[darkest];
		whitestval = frequency[whitest];
		peakval = frequency[peak];
		maxval = frequency[last_local_max];

		this->point = nullptr;
	}

	zonestat(const Point point, const frequencies frequency)
	{
		darkest = 0, whitest = 255, peak = 0, last_local_max = 0;
		int idx = 0;
		while (idx < 256 && frequency[idx] == 0) idx++;
		darkest = idx;
		for (; idx < 256; idx++)
			if (frequency[idx] >= frequency[peak]) peak = idx;
		idx = 255;
		while (idx > 0 && frequency[idx] == 0) idx--;
		whitest = idx;
		last_local_max = whitest;
		while (idx > 1 && frequency[idx - 1] >= frequency[idx]) idx--;
		last_local_max = idx;

		darkestval = frequency[darkest];
		whitestval = frequency[whitest];
		peakval = frequency[peak];
		maxval = frequency[last_local_max];

		if (point.x >= 0 && point.y >= 0)
			this->point = std::unique_ptr<Point>(new Point(point.x, point.y, this));
	}
};

struct Plane {
	double a;
	double b;
	double c;

	double d;

	Plane()
	{
		a = b = c = d = 0;
	}

	Plane(const Point A, const Point B, const Point C)
	{
		a = (B.y - A.y) * (C.zone->scale_factor - A.zone->scale_factor) -
			(C.y - A.y) * (B.zone->scale_factor - A.zone->scale_factor);
		b = (C.x - A.x) * (B.zone->scale_factor - A.zone->scale_factor) -
			(B.x - A.x) * (C.zone->scale_factor - A.zone->scale_factor);
		c = (B.x - A.x) * (C.y - A.y) - (C.x - A.x) * (B.y - A.y);

		d = -(a * A.x) - (b * A.y) - (c * A.zone->scale_factor);
	}

	double Interpolate(const int x, const int y)
	{
		return -1 * ((a * x + b * y + d) / c);
	}
};

struct Triangle {
	Point A;
	Point B;
	Point C;

	Plane plane;

	Triangle(Point A, Point B, Point C)
	{
		this->A = A;
		this->B = B;
		this->C = C;

		plane = Plane(A, B, C);
	}
};

class Document
{
private:
	ILuint image_name;
	/// <summary>
	/// Converts image data read by DevIL and stores it in the rgbpixels matrix.<br/>
	/// (Used when reading in image.)
	/// </summary>
	void ILDataToRGB();
	/// <summary>
	/// Converts image data stored in rgbpixels to data expected by DevIL and returns a pointer to converted bytes.<br/>
	/// (Used when writing colour image.)
	/// </summary>
	/// <param name="whiten">Wether to colour all pixels that are white in greyscale, white.</param>
	/// <returns>a pointer to converted bytes</returns>
	ILubyte* RGBToILData(bool whiten = false) const;
	/// <summary>
	/// Converts image data stored in gscpixels to data expected by DevIL and returns a pointer to converted bytes.<br/>
	/// (Used when writing greyscale image.)
	/// </summary>
	/// <returns>a pointer to converted bytes</returns>
	ILubyte* GScToILData() const;

	/// <summary>
	/// Vector matrix of the colour pixels of this Document.
	/// </summary>
	std::vector<std::vector<RGBA>> rgbpixels;
	/// <summary>
	/// Vector matrix of the greyscale pixels of this Document.
	/// </summary>
	std::vector<std::vector<ubyte>> gscpixels;

	int height; /*DevIL gives these as int*/
	int width;  /*DevIL gives these as int*/

	/// <summary>
	/// Vector matrix of zonestat objects.
	/// </summary>
	std::vector<std::vector<std::unique_ptr<zonestat>>> zones;
	/// <summary>
	/// The zone used for calculating the scaling factor of other zones.
	/// </summary>
	zonestat* reference_zone;

	bool is_background_removed = false;

	std::string image_path = "";

	/// <summary>
	/// Initialises rgbpixels vector based on image width and height.
	/// </summary>
	void initPixels();
	/// <summary>
	/// Initialises gscpixels vector based on image width and height.
	/// </summary>
	void initGScPixels();
	void initFrequency(frequencies& frequency);

	static bool is_init;

public:
	Document(int width = 0, int height = 0);
	Document(std::string file);
	Document(ILuint image_name, std::string file);
	Document(const Document& that);

	/// <summary>
	/// Initialises the DevIL library for use.
	/// </summary>
	static void init();

	/// <summary>
	/// Checks wether the given min and max x,y coordinates are inside of the Document's coordinates.
	/// If not it corrects them.
	/// Also if max is 0 it corrects it to the maximum width/height.
	/// </summary>
	void ValidateDimensions(int& minWidth, int& minHeight, int& maxWidth, int& maxHeight) const;

	Document& operator=(const Document& that);

	inline int GetHeight() const { return height; }
	inline int GetWidth() const { return width; }

	/// <summary>
	/// Was a background removal performed on this object?
	/// </summary>
	inline bool IsBackgroundRemoved() { return is_background_removed; }

	std::vector<std::vector<RGBA>>& GetRGBPixels();
	std::vector<std::vector<ubyte>>& GetGScPixels();

	/// <summary>
	/// Creates the greyscale pixel matrix from the RGB(A) image.
	/// </summary>
	void RGBtoGSc();

	/// <summary>
	/// Returns a new document object of the specified rectangle inside of this document.
	/// </summary>
	/// <param name="minWidth"> The width (x) position of the upper left corner of the custom rectangle.</param>
	/// <param name="minHeight">The height (y) position of the upper left corner of the custom rectangle.</param>
	/// <param name="maxWidth"> The width (x) position of the lower right corner of the custom rectangle.</param>
	/// <param name="maxHeight">The height (y) position of the lower right corner of the custom rectangle.</param>
	/// <returns></returns>
	Document Crop(int minWidth, int minHeight, int maxWidth, int maxHeight) const;

	/// <summary>
	/// Returns an array containing the amount of pixels that are a certain colour. The array contains all the possible colours and the indices are the colour value.
	/// In case you don't want to get the colour spectrum of the entire Document you can specify a rectangle inside the Document.
	/// </summary>
	/// <param name="minWidth"> The width (x) position of the upper left corner of the custom rectangle.</param>
	/// <param name="minHeight">The height (y) position of the upper left corner of the custom rectangle.</param>
	/// <param name="maxWidth"> The width (x) position of the lower right corner of the custom rectangle.</param>
	/// <param name="maxHeight">The height (y) position of the lower right corner of the custom rectangle.</param>
	/// <returns>An vector of the amount of pixels that are a certain luminance for every luminance value (0-255).</returns>
	frequencies GetGScFrequency(int minWidth = 0, int minHeight = 0, int maxWidth = 0, int maxHeight = 0);

	/// <summary>
	/// Creates the specified number of zones (if it's possible whitin the size constraints),
	/// calculates their values and populates zones vector.<br/>
	/// Also finds the reference zone used for scaling.
	/// </summary>
	/// <param name="amount">The number of zones for the image to be split into. The more zones there are, the smaller their size will be.<br/>
	/// If zone size is not within the limits it will automatically adjust.<br/>If 0 (or not specified) it will create the smallest zones possible. (MIN_ZONE_SIDE_LENGTH � MIN_ZONE_SIDE_LENGTH)</param>
	/// <returns>The actual number of zones created.</returns>
	int InitZones(int amount = 0);

	/// <summary>
	/// Scales the image to its brightest part to minimise the significance of shadows and bad lighting conditions.
	/// </summary>
	void ScaleLighting(bool scale_rgb = false);

	/// <summary>
	/// Returns a vector of pointers to zones that will get used in the triangulation.
	/// </summary>
	std::vector<zonestat*> ZonesToTriangulate();

	/// <summary>
	/// Returns the triangles created by the triangulation process.
	/// </summary>
	/// <param name="zones">The vector of zone pointers to be used to triangulate.</param>
	std::vector<Triangle> Triangulation(std::vector<zonestat*> zones);

	/// <summary>
	/// Can increase or decrease contrast on image.
	/// </summary>
	/// <param name="adjustment">The amount of change. Positive to increase, negative to decrease.</param>
	void AdjustContrast(const char adjustment = 12);

	/// <summary>
	/// Whitens the background of the greyscale Document.
	/// </summary>
	void DeleteBackground();

	//Read:
	/// <summary>
	/// Reads the data from Document file.
	/// </summary>
	/// <returns>true if successful, false if not.</returns>
	bool Read();

	//Writes:
	//Colour
	/// <summary>
	/// Saves document to a colour image file. Every non white pixel is its original colour.
	/// </summary>
	/// <param name="file_to_create">The path of the file to be created/overriden.</param>
	/// <param name="file_type"> The desired file type. Given by DevIL's internal format identifier codes.<br/>DevIL will use the file type given in the path by default.</param>
	/// <returns>true if file creation is successful, false otherwise.</returns>
	bool WriteRGB(std::string file_to_create, int file_type = -1) const;

	//Greyscale
	/// <summary>
	/// Saves document to a greyscale image file.
	/// </summary>
	/// <param name="file_to_create">The path of the file to be created.</param>
	/// <param name="file_type"> The desired file type. Given by DevIL's internal format identifier codes.<br/>DevIL will use the file type given in the path by default.</param>
	/// <returns>true if file creation is successful, false otherwise.</returns>
	bool WriteGSc(std::string file_to_create, int file_type = -1);
};